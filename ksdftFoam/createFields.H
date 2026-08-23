// ------------------------------------------------------------------------- //
//  Fields and constants for ksdftFoam (Kohn-Sham DFT, atomic units)
//
//  KS equation      [ -1/2 lap + Veff ] psi = eps psi
//  Veff             = Vext + VH + Vxc
//  Hartree (M2)     lap(VH) = -4 pi n ,   n = f |psi|^2
//  LDA exchange(M3) Vx = -(3 n / pi)^{1/3}
//
//  Units: Hartree atomic units (length = Bohr, energy = Ha). Dimension-wise we
//  reuse the laplacianFoam convention (D ~ m^2/s, potentials ~ 1/s, psi
//  dimensionless) so the checker stays happy; numerically 1 m = 1 Bohr,
//  1 1/s = 1 Ha.
// ------------------------------------------------------------------------- //

Info<< "Reading dftProperties\n" << endl;

IOdictionary dftProperties
(
    IOobject
    (
        "dftProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Kinetic prefactor 1/2 (hbar^2/2m in atomic units).
const dimensionedScalar D
(
    "D",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    dftProperties.getOrDefault<scalar>("D", 0.5)
);

// Occupation of the (single, M1-M3) orbital: 1 for H, 2 for He (closed shell).
const scalar occupation =
    dftProperties.getOrDefault<scalar>("occupation", 1.0);

// Self-consistency switches (M1: both off)
const bool hartreeOn = dftProperties.getOrDefault<bool>("hartree", false);
const word xcType    = dftProperties.getOrDefault<word>("xc", "none");

// Convergence of the orbital energy between imaginary-time steps. 0 = never.
const scalar convergenceTol =
    dftProperties.getOrDefault<scalar>("convergenceTol", 1e-8);

// Orbital normalisation target:  int |psi|^2 dV = 1  (one orbital).
const dimensionedScalar targetNorm
(
    "targetNorm",
    dimensionSet(0, 3, 0, 0, 0, 0, 0),
    dftProperties.getOrDefault<scalar>("targetNorm", 1.0)
);

Info<< "Reading field psi (KS orbital)\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Vext (nuclear potential)\n" << endl;
volScalarField Vext
(
    IOobject
    (
        "Vext",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Electron density  n = f |psi|^2   (dimensionless, like |psi|^2)
volScalarField n
(
    IOobject
    (
        "n",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    occupation*sqr(psi)
);

// Hartree potential (M2). READ_IF_PRESENT so a case may prescribe its BCs;
// fall back to zeroGradient when absent (M1: unused, stays zero).
volScalarField VH
(
    IOobject
    (
        "VH",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("VH", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
    "zeroGradient"
);

// Exchange-correlation potential (M3).
volScalarField Vxc
(
    IOobject
    (
        "Vxc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Vxc", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
    "zeroGradient"
);

// 4*pi with the dimensions that make  lap(VH) == -fourPi*n  consistent.
const dimensionedScalar fourPi
(
    "fourPi",
    dimensionSet(0, -2, -1, 0, 0, 0, 0),
    4.0*Foam::constant::mathematical::pi
);
