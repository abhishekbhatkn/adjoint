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
    volScalarField eta_MMA, eta_old, dfdeta;
    scalar penalty, volumeConstraint, designVolume, optEpsilon, porosity_s, porosity_f;
    label nOptSteps, maxoutit, maxPiggyLoop, designSize;
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
        sens(sens),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        costFunctionPatches(mesh.solutionDict().subDict("SIMPLE").lookup("costFunctionPatches")),

        check(runTime),
    	checkDict(runTime),
    	checkDB(runTime, checkDict),
    	eta_MMA("eta_MMA",eta),
    	eta_old("eta_old",eta),
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
	runTime.writeNow();
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
        const Foam::wordList solidSpaceZones(
            mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<Foam::wordList>("solidSpace",Foam::wordList())
    );
    forAll(solidSpaceZones, i){
            const label cellZoneID = mesh.cellZones().findZoneID(solidSpaceZones[i]);
            solidSpaceCells.append( mesh.cellZones()[cellZoneID] );
        }

    forAll(solidSpaceCells,i){
    	const label j = solidSpaceCells[i];
        eta[j] = 0.0;
        }
    runLoop();
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
	//value[1] = 0.0; //((val/designVolume)-(volumeConstraint));


	Info << "fval1 = " << value[0] << " fval2 = " << value[1] << endl;
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
	J = 10000*J;
	Info<< "cost pressure: " << J << endl;
	return 10000*J;
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
		    
		    Info << "piggy: " << optStep << " " << runTime.timeName() << " sensSum " << sensSum << " dSensSum " << dSensSum << " cost " << J << endl;
		    
		    oldSensSum = sensSum;

		    AD::resetTapeTo(reset_to);
		    AD::zeroAdjointVector();
		    AD::switchTapeToPassive();
		    
	    }
	    runTime.writeNow();
	    Info << "avg/max eta,: " << gAverage(eta) << " " << gMax(eta) << endl;
	    Info<< "End of checkpointerLoop\n" << endl;
    }
    
    bool GlobalSolver() {

	start();
    	Info<< "Starting Global Program\n" << endl;
    	scalar J = 0, oldJ = 0.0;
    	List<scalar> fval, oldfval;

    	scalar check = 1.0;
	eta_old = eta;
	eta_old.write();
	
	GCMMASolver gcmma(mesh,designSpaceCells,1);
   	
   	for (int iter = 0; check>optEpsilon && iter < maxoutit; ++iter) {
   	    checkpointerLoop();
	    J = calcCost();
	    fval = calcFval();
	    if (iter == 0) {
   	    	oldJ = J;
		oldfval = fval;
		Info<< " Start: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << "\n" << endl;
	    }   
	    
	    gcmma.MMAUpdate(eta_MMA, eta, J, sens, fval, dfdeta);
	    eta = eta_MMA;
    	    runLoop();
	    J = calcCost();
	    fval = calcFval();

	    check = 0.0;
	    forAll(designSpaceCells,i){
	    	const label j = designSpaceCells[i];
		check = max(check,mag(eta[j]-eta_old[j]));//max(senscheck,sens[j]);
	    }
	    Foam::reduce(check,maxOp<scalar>());
	    eta_old = eta;
	    eta_old.write();
	    oldJ = J;
	    oldfval = fval;
	    
    	    Info<< "Outer Loop "<< iter << " End: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << "\n" << endl;
    	}
	Info<< "GLobal Solver End\n" << endl;
    	return true;
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

