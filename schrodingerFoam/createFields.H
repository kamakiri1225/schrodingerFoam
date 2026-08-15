// ------------------------------------------------------------------------- //
//  Fields and physical constants for schrodingerFoam
//
//  Gross-Pitaevskii / non-linear Schroedinger equation
//
//      i d(Psi)/dt = -D lap(Psi) + ( Vext + g |Psi|^2 ) Psi
//
//  The complex wave function Psi is split into real and imaginary parts:
//
//      Psi = Psire + i Psiim
//
//  Dimensions follow laplacianFoam's convention so that the dimension
//  checker stays happy (D ~ m^2/s, g and Vext ~ 1/s, Psi dimensionless).
// ------------------------------------------------------------------------- //

Info<< "Reading gpProperties\n" << endl;

IOdictionary gpProperties
(
    IOobject
    (
        "gpProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// hbar / (2 m). Plays the role of the diffusion coefficient.
dimensionedScalar D("D", dimensionSet(0, 2, -1, 0, 0, 0, 0), gpProperties);

// Non-linear (interaction) coupling constant g.
dimensionedScalar g("g", dimensionSet(0, 0, -1, 0, 0, 0, 0), gpProperties);

Info<< "Reading field Psire (real part)\n" << endl;
volScalarField Psire
(
    IOobject
    (
        "Psire",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Psiim (imaginary part)\n" << endl;
volScalarField Psiim
(
    IOobject
    (
        "Psiim",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// External potential Vext. Read from the time directory if present,
// otherwise default to a homogeneous (zero) potential.
Info<< "Reading/creating field Vext (external potential)\n" << endl;
volScalarField Vext
(
    IOobject
    (
        "Vext",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), Zero)
);

// |Psi|^2 (density) -- written out for post-processing.
volScalarField magSqrPsi
(
    IOobject
    (
        "magSqrPsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqr(Psire) + sqr(Psiim)
);

// Phase arg(Psi) -- written out for post-processing (vortex detection).
volScalarField phase
(
    IOobject
    (
        "phase",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, Zero)
);
