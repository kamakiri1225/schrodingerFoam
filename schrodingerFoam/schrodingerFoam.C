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

    Info<< "Solver mode: " << mode << nl << endl;

    // ------------------------------------------------------------------- //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (mode == "imaginaryTime")
        {
            // d(Psi)/dtau = D lap(Psi) - (W - mu) Psi   (W = Vext + g|Psi|^2)
            volScalarField W("W", Vext + g*(sqr(Psire) + sqr(Psiim)));

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

            const volScalarField W0
            (
                Vext + g*(sqr(Psire0) + sqr(Psiim0))
            );
            const volScalarField Hv0(-D*fvc::laplacian(Psiim0) + W0*Psiim0);
            const volScalarField Hu0(-D*fvc::laplacian(Psire0) + W0*Psire0);

            // Crank-Nicolson, Picard (Gauss-Seidel) outer iteration:
            //   du/dt = +H v ,   dv/dt = -H u
            for (int corr = 0; corr < nCorr; ++corr)
            {
                volScalarField W
                (
                    Vext + g*(sqr(Psire) + sqr(Psiim))
                );
                const volScalarField Hv(-D*fvc::laplacian(Psiim) + W*Psiim);

                Psire = Psire0 + 0.5*dt*(Hv + Hv0);
                Psire.correctBoundaryConditions();

                W = Vext + g*(sqr(Psire) + sqr(Psiim));
                const volScalarField Hu(-D*fvc::laplacian(Psire) + W*Psire);

                Psiim = Psiim0 - 0.5*dt*(Hu + Hu0);
                Psiim.correctBoundaryConditions();
            }

            Info<< "  norm = "
                << fvc::domainIntegrate(sqr(Psire) + sqr(Psiim)).value()
                << endl;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mode " << mode
                << " (use realTime or imaginaryTime)"
                << exit(FatalError);
        }

        // derived fields for post-processing
        magSqrPsi = sqr(Psire) + sqr(Psiim);
        phase = atan2(Psiim, Psire);

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
