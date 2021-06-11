/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible flows with turbulence modelling.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "GCMMASolver.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "CheckController.H"
#include "CheckInterface.H"
#include "CheckDict.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class MMAProgram {

private:
    Foam::Time& runTime;
    Foam::fvMesh& mesh;
    Foam::simpleControl& simple;
    Foam::volScalarField& p;
    Foam::volVectorField& U;
    Foam::singlePhaseTransportModel& laminarTransport;
    autoPtr<incompressible::turbulenceModel>& turbulence;
    Foam::IOMRFZoneList& MRF;
    fv::options& fvOptions;
    volScalarField& porosity;
    volScalarField& eta;
    volScalarField& eta_old1;
    volScalarField& eta_old2;
    volScalarField& low;
    volScalarField& upp;
    volScalarField& sens;
    surfaceScalarField& phi;
  
    label&  pRefCell;
    scalar& pRefValue;
    scalar& cumulativeContErr;
    const Foam::wordList costFunctionPatches;
    CheckInterface check;
    CheckDict checkDict;
    CheckDatabase checkDB;
    
    // member variables
    volScalarField dfdeta;
    scalar penalty, volumeConstraint, designVolume, optEpsilon, porosity_s, porosity_f, factorP, asyminit, asymdec, asyminc;
    label nOptSteps, maxPiggyLoop, designSize, maxMMAiter, MMALoop;
    List<label> designSpaceCells;
    List<label> solidSpaceCells;
    
    Foam::scalar fluidCp, fluidRho;
    


public:
    MMAProgram
    (
        Foam::Time& runTime,
        Foam::fvMesh& mesh,
        Foam::simpleControl& simple,
        Foam::volScalarField& p,
        Foam::volVectorField& U,
        Foam::singlePhaseTransportModel& laminarTransport,
        autoPtr<incompressible::turbulenceModel>& turbulence,
        IOMRFZoneList& MRF,
        fv::options& fvOptions,
        volScalarField& porosity,
        volScalarField& eta,
        volScalarField& eta_old1,
        volScalarField& eta_old2,
        volScalarField& low,
        volScalarField& upp,
        volScalarField& sens,
        surfaceScalarField& phi,
        label&  pRefCell,
        scalar& pRefValue,
        scalar& cumulativeContErr
    ) :
        runTime(runTime),
        mesh(mesh),
        simple(simple),
        p(p),
        U(U),
        laminarTransport(laminarTransport),
        turbulence(turbulence),
        MRF(MRF),
        fvOptions(fvOptions),
        porosity(porosity),
        eta(eta),
        eta_old1(eta_old1),
        eta_old2(eta_old2),
        low(low),
        upp(upp),
        sens(sens),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        costFunctionPatches(mesh.solutionDict().subDict("SIMPLE").lookup("costFunctionPatches")),

        check(runTime),
    	checkDict(runTime),
    	checkDB(runTime, checkDict),
    	dfdeta("dfdeta",eta)
{
	label tapeSizeMB = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("tapeSizeMB",4096);
	Info << "Creating Tape, size: " << tapeSizeMB << endl;
	AD::createGlobalTape(tapeSizeMB/Pstream::nProcs());

        AD::switchTapeToPassive();
}

    void runLoop() {
	Info<< "\nStarting time loop\n" << endl;
	label iter = 0;
	bool check = true;
	while (iter < nOptSteps && check && simple.loop())
	{
	    	iter++;
		Info<< "Time = " << runTime.timeName() << nl << endl;
		
		Foam::volScalarField pold("pold",p);
		Foam::volVectorField Uold("Uold",U);
		// --- Pressure-velocity SIMPLE corrector
		{
		    #include "porosity.H"
		    #include "UEqn.H"
		    #include "pEqn.H"
		}
		
		laminarTransport.correct();
		turbulence->correct();

		//runTime.write();
		scalar maxdiff = 0.0;
		forAll(p,i) {
			maxdiff = max(maxdiff,mag(Uold[i].x()-U[i].x()));
			maxdiff = max(maxdiff,mag(Uold[i].y()-U[i].y()));
			maxdiff = max(maxdiff,mag(Uold[i].z()-U[i].z()));
			maxdiff = max(maxdiff,mag(pold[i]-p[i]));
		}
		Foam::reduce(maxdiff,maxOp<scalar>());
		Info << "MaxDiff: " << maxdiff << endl;
		if (maxdiff < optEpsilon && iter > 100) {
			check = false;
		}
		runTime.printExecutionTime(Info);
	}
	//runTime.writeNow();
	Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = "   << runTime.elapsedClockTime() << " s"
	<< nl << endl;
    }
    
    void start() {
        #include "settings.H"
        runLoop();
        Info << "Initialization Complete ! " << endl;

    }
    
    void initialGuess() {
	    forAll(designSpaceCells,i){
		    const label j = designSpaceCells[i];
		    eta[j] = volumeConstraint; //*designVolume/(mesh.V()[j]*designSpaceCells.size());
	    }
	    List<scalar> value(1,0.0);
	    value = calcFval();
	    Info << "Eta Set ! " << endl;

    }
    
    List<scalar> calcFval() {  
	scalar val = 0;
	List<scalar> value(1,0.0);
	forAll(designSpaceCells,i){
	    const label j = designSpaceCells[i];
	    val += eta[j] * mesh.V()[j];
	}
	Foam::reduce(val,sumOp<scalar>());
	value[0] = 100*((volumeConstraint)-(val/designVolume));
	//value[1] = 100*((val/designVolume)-(volumeConstraint));


	Info << "fval = " << value << endl;
	return value;
    }
    
    scalar calcCost(){
        Info << "ACTION::calcPrCost" << endl;
	scalar J = 0;
	forAll(costFunctionPatches,cI)
	{
	Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI] );
	const Foam::fvPatch& patch = mesh.boundary()[patchI];
	J += gSum
		(
		    - phi.boundaryField()[patchI]*(p.boundaryField()[patchI]
		    + 0.5*magSqr(phi.boundaryField()[patchI]/patch.magSf()))
		);
	}
	J = factorP*J;
	Info<< "cost pressure: " << J << endl;
	return J;
    }
       
    void checkpointerLoop() {
	    Info<< "Starting Piggyback Loop: " << "\n" << endl;
	    
	    turbulence->validate();

	    bool firstStep = true;
	    scalar J = 0, oldSensSum = 0;
	    

	    scalar dSensSum = std::numeric_limits<double>::max();
	    //scalar maxSens = 0, minSens = 0;
	    
	    label optStep = 0;
	    bool check = true;
	    while ( optStep < maxPiggyLoop && check && simple.loop())
	    {
		    AD::switchTapeToActive();
		    
		    AD::registerInputVariable(eta.begin(),eta.end());
		    
		    AD::position_t reset_to = AD::getTapePosition();

	    	    checkDB.registerAdjoints(); // store tape indices

		    Info<< "Time = " << runTime.timeName() << nl << endl;
		    
		    // --- Pressure-velocity SIMPLE corrector
		    {
			#include "porosity.H"
			#include "UEqn.H"
			#include "pEqn.H"
		    }
		    
		    laminarTransport.correct();
		    
		    // turbulence correct
		    bool frozenTurbulence = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<bool>("frozenTurbulence",false);
		    if(frozenTurbulence){
			Info << "\nWARNING: Calculating Adjoints with frozen Turbulence assumption\n" << endl;
			AD::switchTapeToPassive();
		    }

		    turbulence->correct();
		    
		    if(frozenTurbulence){
		       AD::switchTapeToActive();
		    }
		    
		    J = calcCost();
		    checkDB.registerAsOutput();
		    if(Pstream::master()){
			AD::derivative(J) = 1.0;
		    }

		    if(!firstStep){
			checkDB.restoreAdjoints();
		    }else{
			firstStep = false;
		    }

	      	    AD::interpretTape(); //_to(interpret_to);
	      	    scalar norm2 = checkDB.calcNormOfStoredAdjoints();
		    checkDB.storeAdjoints();
		    scalar sensSum  = 0.0;

		    forAll(designSpaceCells,i){
			const label j = designSpaceCells[i];
			sens[j] = AD::derivative(eta[j]) * designSize*mesh.V()[j]/designVolume;
			sensSum += mag(sens[j]);
		    }
		    
		    Foam::reduce(sensSum,sumOp<scalar>());
		    
		    dSensSum = mag(sensSum - oldSensSum)/designSize;
		    optStep++;
		    
		    if ( dSensSum < optEpsilon && optStep > 5) {
		    	check = false;
		    }
		    
		    Info << "piggy: " << optStep << " " << runTime.timeName() << " sensSum " << sensSum << " dSensSum " << dSensSum << " norm2: " << norm2 << " cost " << J << endl;
		    
		    oldSensSum = sensSum;

		    AD::resetTapeTo(reset_to);
		    AD::zeroAdjointVector();
		    AD::switchTapeToPassive();
		    
	    }
	    //runTime.writeNow();
	    Info << "avg/max eta,: " << gAverage(eta) << " " << gMax(eta) << endl;
	    Info<< "End of checkpointerLoop\n" << endl;
    }
    
    bool GlobalSolver() {

	start();
    	Info<< "Starting Global Program\n" << endl;
    	
	GCMMASolver gcmma (
		mesh
		, designSpaceCells
		, 1
		, low
		, upp
		, asyminit
		, asymdec
		, asyminc
		, MMALoop
		);
	
	scalar check = 1.0;
	while ( MMALoop < maxMMAiter && check > optEpsilon) {
	   	checkpointerLoop();
		scalar J = calcCost();
		List<scalar> fval = calcFval();
		if (MMALoop == 0) {
			Info<< " Start: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << "\n" << endl; 
			runTime.writeNow();  
		}
		gcmma.MMAUpdate(eta, eta_old1, eta_old2, J, sens, fval, dfdeta);
	    	runLoop();
		J = calcCost();
		fval = calcFval();
		check = etaCheck();
	    	Info<< "Loop "<< MMALoop << " End: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " asyminit: " << asyminit << " asymdec: " << asymdec << " asyminc: " << asyminc << endl;
	    	++MMALoop;
	    	runTime.writeNow();
    	}
	Info<< "GLobal Solver End\n" << endl;
    	return true;
    }
    
    scalar etaCheck() {
    	scalar check = 0.0;
	forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
		check = max(check,mag(eta[j]-eta_old1[j]));
	}
	Foam::reduce(check,maxOp<scalar>());
	return check;
    }
    
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
    
//------------------------------------------------------------------------------------------------------//
    MMAProgram program
    (
         runTime,
         mesh,
         simple,
         p,
         U,
         laminarTransport,
         turbulence,
         MRF,
         fvOptions,
         porosity,
         eta,
         eta_old1,
         eta_old2,
         low,
         upp,
         sens,
         phi,
         pRefCell,
         pRefValue,
         cumulativeContErr
    );
     
    program.GlobalSolver();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

