/*
Info<< "Reading field T\n" <<endl;
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
*/

/*
dimensionedScalar DTs
(
	transportProperties.lookup("DTs")
);
    
dimensionedScalar DTf
(
	transportProperties.lookup("DTf")
);

dimensionedScalar qf
(
	transportProperties.lookup("qf")
);

dimensionedScalar qs
(
	transportProperties.lookup("qs")
);
*/

/*
dimensionedScalar optStepwidth = dimensionedScalar("optStepwidth", dimless, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("optStepwidth",0.1));
    if (!mesh.solutionDict().subDict("SIMPLE").found("optStepwidth"))
    {
        Info << "Warning: Keyword optStepwidth not found in fvSolution/SIMPLE. Default set to 0.1" << endl;
    }

    #include "Random.H"
    Foam::Random randObj(time(NULL));
    forAll(designSpaceCells,i){
    	const label j = designSpaceCells[i];
       	eta[j] = randObj.bit();
    }
*/


/*
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
    dimensionedScalar("0", dimViscosity, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
	DT = DTs + (DTf - DTs) * eta * (scalar(1) + penalty) / (eta + penalty);
	Info << zeroGradientFvPatchScalarField::typeName << endl;
	
volScalarField heatq
(
    IOobject
    (
        "heatq",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0,1,-2,1,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);
	heatq = qs + (qf - qs) * eta * (scalar(1) + penalty) / (eta + penalty);
	Info << zeroGradientFvPatchScalarField::typeName << endl;
*/
