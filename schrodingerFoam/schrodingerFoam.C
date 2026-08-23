/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    schrodingerFoam

Description
    Unified solver for the Gross-Pitaevskii (non-linear Schroedinger) equation

        i d(Psi)/dt = -D lap(Psi) + ( Vext + g |Psi|^2 ) Psi

    The complex field is split into real / imaginary parts, Psi = Psire + i Psiim.

    Two modes, selected by "mode" in constant/gpProperties:

      * realTime      : unitary (norm-conserving) time evolution using a
                        Crank-Nicolson scheme solved by Picard iteration.
                        Forward Euler is UNCONDITIONALLY unstable for the
                        Schroedinger equation, which is why the naive
                        laplacianFoam modification diverged; Crank-Nicolson
                        has |amplification| = 1 for the linear part.

      * realTimeVisscher : the SAME real-time evolution, but an explicit
                        staggered leap-frog (Visscher, Computers in Physics 5,
                        596, 1991). Re is stored at integer steps, Im at half-
                        integer steps. No Picard iteration (about 1/nCorrectors
                        the cost of realTime) and it EXACTLY conserves the
                        discrete norm  u^2 + v_{n+1/2} v_{n-1/2}. Stable for
                        dt <= 2/||H|| (same restriction as realTime). Add-on
                        alternative; realTime is left unchanged.

      * imaginaryTime : gradient flow  d(Psi)/dtau = -(H - mu) Psi , solved
                        fully implicitly (unconditionally stable, large dtau).
                        Projects an arbitrary start onto the lowest-energy
                        state -> used to prepare the initial condition.

    ---------------------------------------------------------------------------
    Relation to laplacianFoam (what was changed)
    ---------------------------------------------------------------------------
    laplacianFoam solves a single real field:   ddt(T) = laplacian(DT,T).
    This solver differs as follows:

      1. Fields: ONE real field T  ->  TWO real fields Psire, Psiim
         (the real and imaginary parts of the complex wave function).
         [see createFields.H]

      2. Extra physics: an external potential Vext and the GP non-linearity
         g*|Psi|^2 are added as a "potential" term  W = Vext + g|Psi|^2 .

      3. Time scheme: laplacianFoam uses a plain implicit ddt.  The naive
         Schroedinger port (note.com) kept forward Euler via fvc::laplacian,
         which is UNCONDITIONALLY UNSTABLE for i*dPsi/dt = H*Psi.  Here:
           - realTime      -> Crank-Nicolson (norm conserving), Picard loop
           - imaginaryTime -> fully implicit backward Euler + renormalisation

      4. Coupling: the two equations are coupled through H (the i in front of
         d/dt rotates real <-> imaginary), handled by the outer corrector loop.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // ----- solution controls ------------------------------------------- //
    const word mode(gpProperties.get<word>("mode"));

    const label nCorr =
        gpProperties.getOrDefault<label>("nCorrectors", 4);

    const bool normalize =
        gpProperties.getOrDefault<bool>("normalize", false);

    const scalar targetNorm =
        gpProperties.getOrDefault<scalar>("targetNorm", 1.0);

    const scalar convTol =
        gpProperties.getOrDefault<scalar>("convergenceTol", 0.0);

    // Optional: switch the external potential OFF at t >= releaseTime
    // (e.g. releasing a trapped cloud into free expansion). Default: never.
    const scalar releaseTime =
        gpProperties.getOrDefault<scalar>("releaseTime", GREAT);

    // Optional: subtract a chemical potential in realTime, i.e. integrate
    //     i d(Psi)/dt = (H - mu) Psi
    // instead of i d(Psi)/dt = H Psi. This moves to the frame co-rotating with
    // the background global phase e^{-i mu t}, so the PHASE field stops
    // flickering; the density and all vortex/soliton dynamics are unchanged.
    // How mu is chosen (select the option in constant/gpProperties):
    //   * dynamicMu true          -> recompute mu = <Psi|H|Psi>/<Psi|Psi> each
    //                                step (use when mu is not known, e.g. traps)
    //   * dynamicMu false + muShift-> subtract the fixed constant muShift
    //                                (e.g. g*n0 for a uniform bulk)
    // Defaults (dynamicMu false, muShift 0) reproduce plain i dPsi/dt = H Psi.
    const bool dynamicMu =
        gpProperties.getOrDefault<bool>("dynamicMu", false);
    const dimensionedScalar muShift
    (
        "muShift",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        gpProperties.getOrDefault<scalar>("muShift", 0.0)
    );

    Info<< "Solver mode: " << mode << nl << endl;

    // Working copy holding the previous half-step imaginary part v^{n-1/2}
    // (only used by the realTimeVisscher scheme; NO_WRITE).
    volScalarField PsiimPrev("PsiimPrev", Psiim);

    // One-time initialisation for the Visscher staggered leap-frog:
    // shift the imaginary part back by half a step, v^{1/2} = v^0 - (dt/2) H u^0,
    // so that Re lives at integer times and Im at half-integer times.
    if (mode == "realTimeVisscher")
    {
        const dimensionedScalar dt0 = runTime.deltaT();
        const scalar trapOn0 = (runTime.value() < releaseTime) ? 1.0 : 0.0;
        const volScalarField Vnow0(trapOn0*Vext);

        dimensionedScalar muEff0("muEff0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0);
        if (dynamicMu)
        {
            const volScalarField Wm(Vnow0 + g*(sqr(Psire) + sqr(Psiim)));
            const volScalarField Hr(-D*fvc::laplacian(Psire) + Wm*Psire);
            const volScalarField Hi(-D*fvc::laplacian(Psiim) + Wm*Psiim);
            muEff0 = fvc::domainIntegrate(Psire*Hr + Psiim*Hi)
                   / fvc::domainIntegrate(sqr(Psire) + sqr(Psiim));
        }
        else
        {
            muEff0 = muShift;
        }
        const volScalarField W0k(Vnow0 + g*(sqr(Psire) + sqr(Psiim)) - muEff0);
        const volScalarField Hu0k(-D*fvc::laplacian(Psire) + W0k*Psire);
        Psiim = Psiim - 0.5*dt0*Hu0k;    // v^0 -> v^{1/2}
        Psiim.correctBoundaryConditions();
        Info<< "Visscher: staggered Im by -dt/2 (v^{1/2})\n" << endl;
    }

    // ------------------------------------------------------------------- //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Effective potential this step: the trap is released (set to 0) once
        // the run passes releaseTime; otherwise it equals Vext.
        const scalar trapOn = (runTime.value() < releaseTime) ? 1.0 : 0.0;
        const volScalarField Vnow(trapOn*Vext);

        if (mode == "imaginaryTime")
        {
            // d(Psi)/dtau = D lap(Psi) - (W - mu) Psi   (W = Vext + g|Psi|^2)
            volScalarField W("W", Vnow + g*(sqr(Psire) + sqr(Psiim)));

            // chemical potential  mu = <Psi|H|Psi> / <Psi|Psi>
            volScalarField Hre(-D*fvc::laplacian(Psire) + W*Psire);
            volScalarField Him(-D*fvc::laplacian(Psiim) + W*Psiim);

            dimensionedScalar num =
                fvc::domainIntegrate(Psire*Hre + Psiim*Him);
            dimensionedScalar den =
                fvc::domainIntegrate(sqr(Psire) + sqr(Psiim));
            dimensionedScalar mu = num/den;

            // implicit (backward-Euler) update of each component
            solve
            (
                fvm::ddt(Psire)
              - fvm::laplacian(D, Psire)
              + fvm::SuSp(W - mu, Psire)
            );
            solve
            (
                fvm::ddt(Psiim)
              - fvm::laplacian(D, Psiim)
              + fvm::SuSp(W - mu, Psiim)
            );

            // renormalize (optional -- keeps a fixed particle number)
            if (normalize)
            {
                dimensionedScalar curNorm =
                    fvc::domainIntegrate(sqr(Psire) + sqr(Psiim));
                const scalar factor =
                    Foam::sqrt(targetNorm/max(curNorm.value(), SMALL));
                Psire *= factor;
                Psiim *= factor;
            }

            Info<< "  mu = " << mu.value()
                << ",  norm = "
                << fvc::domainIntegrate(sqr(Psire) + sqr(Psiim)).value()
                << endl;
        }
        else if (mode == "realTime")
        {
            const dimensionedScalar dt = runTime.deltaT();

            // start-of-step values
            const volScalarField Psire0(Psire);
            const volScalarField Psiim0(Psiim);

            // Chemical potential to subtract this step (frame co-rotating with
            // the background phase e^{-i mu t}). Either a fixed constant or the
            // instantaneous mu = <Psi|H|Psi>/<Psi|Psi> at the start of the step.
            dimensionedScalar muEff("muEff", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0);
            if (dynamicMu)
            {
                const volScalarField Wmu(Vnow + g*(sqr(Psire0) + sqr(Psiim0)));
                const volScalarField Hre0(-D*fvc::laplacian(Psire0) + Wmu*Psire0);
                const volScalarField Him0(-D*fvc::laplacian(Psiim0) + Wmu*Psiim0);
                const dimensionedScalar num =
                    fvc::domainIntegrate(Psire0*Hre0 + Psiim0*Him0);
                const dimensionedScalar den =
                    fvc::domainIntegrate(sqr(Psire0) + sqr(Psiim0));
                muEff = num/den;
            }
            else
            {
                muEff = muShift;
            }

            // W already carries the -muEff term, so H = -D lap + W integrates
            // (H - muEff) Psi with no other change to the scheme.
            const volScalarField W0
            (
                Vnow + g*(sqr(Psire0) + sqr(Psiim0)) - muEff
            );
            const volScalarField Hv0(-D*fvc::laplacian(Psiim0) + W0*Psiim0);
            const volScalarField Hu0(-D*fvc::laplacian(Psire0) + W0*Psire0);

            // Crank-Nicolson, Picard (Gauss-Seidel) outer iteration:
            //   du/dt = +H v ,   dv/dt = -H u
            for (int corr = 0; corr < nCorr; ++corr)
            {
                volScalarField W
                (
                    Vnow + g*(sqr(Psire) + sqr(Psiim)) - muEff
                );
                const volScalarField Hv(-D*fvc::laplacian(Psiim) + W*Psiim);

                Psire = Psire0 + 0.5*dt*(Hv + Hv0);
                Psire.correctBoundaryConditions();

                W = Vnow + g*(sqr(Psire) + sqr(Psiim)) - muEff;
                const volScalarField Hu(-D*fvc::laplacian(Psire) + W*Psire);

                Psiim = Psiim0 - 0.5*dt*(Hu + Hu0);
                Psiim.correctBoundaryConditions();
            }

            Info<< "  norm = "
                << fvc::domainIntegrate(sqr(Psire) + sqr(Psiim)).value()
                << ",  muShift = " << muEff.value()
                << endl;
        }
        else if (mode == "realTimeVisscher")
        {
            // Explicit staggered leap-frog (Visscher 1991). Re at integer,
            // Im at half-integer steps:
            //   u^{n+1}   = u^n       + dt * H v^{n+1/2}
            //   v^{n+3/2} = v^{n+1/2} - dt * H u^{n+1}
            // No Picard iteration; the potential W (incl. g|Psi|^2) is taken
            // explicitly from the freshest available fields.
            const dimensionedScalar dt = runTime.deltaT();

            // chemical potential to subtract (frame co-rotating with e^{-i mu t})
            dimensionedScalar muEff("muEff", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0);
            if (dynamicMu)
            {
                const volScalarField Wm(Vnow + g*(sqr(Psire) + sqr(Psiim)));
                const volScalarField Hr(-D*fvc::laplacian(Psire) + Wm*Psire);
                const volScalarField Hi(-D*fvc::laplacian(Psiim) + Wm*Psiim);
                muEff = fvc::domainIntegrate(Psire*Hr + Psiim*Hi)
                      / fvc::domainIntegrate(sqr(Psire) + sqr(Psiim));
            }
            else
            {
                muEff = muShift;
            }

            PsiimPrev = Psiim;              // v^{n+1/2}

            // step 1:  u^{n+1} = u^n + dt * H v^{n+1/2}
            {
                const volScalarField W(Vnow + g*(sqr(Psire) + sqr(Psiim)) - muEff);
                const volScalarField Hv(-D*fvc::laplacian(Psiim) + W*Psiim);
                Psire = Psire + dt*Hv;
                Psire.correctBoundaryConditions();
            }
            // step 2:  v^{n+3/2} = v^{n+1/2} - dt * H u^{n+1}
            {
                const volScalarField W(Vnow + g*(sqr(Psire) + sqr(Psiim)) - muEff);
                const volScalarField Hu(-D*fvc::laplacian(Psire) + W*Psire);
                Psiim = Psiim - dt*Hu;
                Psiim.correctBoundaryConditions();
            }

            Info<< "  norm = "
                << fvc::domainIntegrate(sqr(Psire) + Psiim*PsiimPrev).value()
                << ",  muShift = " << muEff.value()
                << endl;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mode " << mode
                << " (use realTime, realTimeVisscher or imaginaryTime)"
                << exit(FatalError);
        }

        // derived fields for post-processing
        if (mode == "realTimeVisscher")
        {
            // integer-time density/phase from the staggered fields:
            //   |Psi|^2 = u^2 + v_{n+1/2} v_{n-1/2}   (Visscher conserved density)
            //   Im at integer time ~ (v_{n+1/2} + v_{n-1/2})/2
            magSqrPsi = sqr(Psire) + Psiim*PsiimPrev;
            phase = atan2(0.5*(Psiim + PsiimPrev), Psire);
        }
        else
        {
            magSqrPsi = sqr(Psire) + sqr(Psiim);
            phase = atan2(Psiim, Psire);
        }

        // convergence stop (mainly for imaginaryTime relaxation)
        if (convTol > 0)
        {
            const scalar res =
                gMax(mag(Psire.primitiveField() - Psire.oldTime().primitiveField()))
              + gMax(mag(Psiim.primitiveField() - Psiim.oldTime().primitiveField()));

            Info<< "  residual = " << res << endl;

            if (res < convTol)
            {
                Info<< "\nConverged (residual < " << convTol << ")\n" << endl;
                runTime.writeAndEnd();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
