/*---------------------------------------------------------------------------*\
    ksdftFoam

    Kohn-Sham density-functional theory on an OpenFOAM mesh
    (real-space, finite-volume, atomic units).

    KS equation (single occupied orbital, M1-M3):

        [ -1/2 lap + Veff ] psi = eps psi ,   Veff = Vext + VH + Vxc

    solved as an IMAGINARY-TIME gradient flow (same idea as schrodingerFoam's
    imaginaryTime mode, see blog #5):

        d(psi)/dtau = +1/2 lap(psi) - Veff psi   + renormalisation each step

    Self-consistency: because Veff depends on the density n = f |psi|^2, the
    potential is refreshed from the CURRENT density every imaginary-time step
    (damped SCF-through-time-stepping):

        n      = f psi^2
        VH     : lap(VH) = -4 pi n        (Poisson; switch "hartree", M2)
        Vxc    = -(3 n / pi)^{1/3}        (LDA/Slater; switch "xc", M3)

    Total energy with double-counting corrections:

        E = f eps - E_H + (E_xc - int Vxc n dV) ,  E_H = 1/2 int VH n dV

    Milestones (docs/research_plan_ksdft_openfoam.md):
      M1 hydrogen: hartree off, xc none  ->  eps -> -0.5 Ha (soft-Coulomb a->0)
      M2 helium  : hartree on
      M3 helium  : hartree on, xc slater
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info<< "\nKohn-Sham DFT: imaginary-time relaxation"
        << "  (hartree=" << hartreeOn << ", xc=" << xcType << ")\n" << endl;

    scalar epsPrev = GREAT;

    while (runTime.loop())
    {
        // --- 1. density from the current orbital
        n = occupation*sqr(psi);

        // --- 2. potentials from the current density
        if (hartreeOn)
        {
            // lap(VH) = -4 pi n   (BCs come from the 0/VH file of the case)
            solve(fvm::laplacian(VH) + fourPi*n);
        }
        if (xcType == "slater")
        {
            // LDA exchange  Vx = -(3 n / pi)^{1/3}
            Vxc = -dimensionedScalar("c", Vxc.dimensions(), 1.0)
                *cbrt(3.0*max(n, dimensionedScalar("n0", n.dimensions(), VSMALL))
                      /constant::mathematical::pi);
        }

        const volScalarField Veff(Vext + VH + Vxc);

        // --- 3. one implicit imaginary-time step of the KS orbital
        solve
        (
            fvm::ddt(psi)
         ==
            D*fvm::laplacian(psi)
          - fvm::Sp(Veff, psi)
        );

        // --- 4. renormalise  int psi^2 dV = 1
        const dimensionedScalar norm(fvc::domainIntegrate(sqr(psi)));
        psi *= sqrt(targetNorm/norm);

        // --- 5. orbital energy  eps = <psi|H|psi>  (psi is normalised)
        const volScalarField Hpsi(-D*fvc::laplacian(psi) + Veff*psi);
        const scalar eps =
            (fvc::domainIntegrate(psi*Hpsi)/targetNorm).value();

        // --- 6. total energy with double-counting corrections
        n = occupation*sqr(psi);
        scalar Etot = occupation*eps;
        if (hartreeOn)
        {
            const scalar EH =
                0.5*fvc::domainIntegrate(VH*n).value();
            Etot -= EH;
        }
        if (xcType == "slater")
        {
            // E_x = -3/4 (3/pi)^{1/3} int n^{4/3};  int Vx n = -(3/pi)^{1/3} int n^{4/3}
            const scalar I43 =
                fvc::domainIntegrate
                (
                    pow(max(n, dimensionedScalar("n0", n.dimensions(), VSMALL)),
                        4.0/3.0)
                ).value();
            const scalar c = Foam::cbrt(3.0/constant::mathematical::pi);
            const scalar Ex    = -0.75*c*I43;
            const scalar IVxcn = -c*I43;
            Etot += (Ex - IVxcn);            // = +1/4 c I43
        }

        Info<< "tau = " << runTime.timeName()
            << "  eps = " << eps
            << "  Etot = " << Etot
            << "  norm = " << norm.value() << endl;

        // --- 7. convergence of the orbital energy
        if (convergenceTol > 0 && mag(eps - epsPrev) < convergenceTol)
        {
            Info<< "\nConverged: |d eps| = " << mag(eps - epsPrev)
                << " < " << convergenceTol << nl
                << "Final  eps  = " << eps << " Ha" << nl
                << "Final  Etot = " << Etot << " Ha" << endl;
            runTime.writeAndEnd();
        }
        epsPrev = eps;

        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
