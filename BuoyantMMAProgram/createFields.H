Info<< "Reading field p\n" << endl;
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
Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

p = thermo.p();

Info<< "Reading field U\n" << endl;
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

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    simple.dict(),
    pRefCell,
    pRefValue
);

mesh.setFluxRequired(p_rgh.name());

//dimensionedScalar initialMass = fvc::domainIntegrate(rho);
//dimensionedScalar totalVolume = sum(mesh.V());

#include "createMRF.H"
#include "createRadiationModel.H"

//const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, simple.dict());
//const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, simple.dict());

#include "createFvOptions.H"

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
/*
dimensionedScalar DTs
(
     transportProperties.lookup("DTs")
);

dimensionedScalar DTf
(
     transportProperties.lookup("DTf")
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
	"T",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField DT
(
    IOobject
    (
	"DT",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimViscosity, 1.0),
    zeroGradientFvPatchScalarField::typeName
);
*/

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

volScalarField eta_MMA
(
    IOobject
    (
        "eta_MMA",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("1.0", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField eta_old
(
    IOobject
    (
	"eta_old",
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

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless/dimTime, 0.0),
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
