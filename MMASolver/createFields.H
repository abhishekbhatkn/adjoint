volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimTemperature, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField eta_old1
(
    IOobject
    (
        "eta_old1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField eta_old2
(
    IOobject
    (
        "eta_old2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField low
(
    IOobject
    (
        "low",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField upp
(
    IOobject
    (
        "upp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);


volScalarField sens
(
    IOobject
    (
        "sens",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField dfdeta
(
    IOobject
    (
        "dfdeta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
