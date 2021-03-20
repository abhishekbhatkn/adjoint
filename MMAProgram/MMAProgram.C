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

#include "costFunctionLibrary.C"
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
    volScalarField& alpha;
    volScalarField& eta;
    volScalarField& eta_MMA;
    volScalarField& sens;
    volScalarField& dfdeta;
    surfaceScalarField& phi;
    
    label&  pRefCell;
    scalar& pRefValue;
    scalar& cumulativeContErr;
    volScalarField& eta_old;
    
    // member variables

    scalar penalty, volumeConstraint, n, optEpsilon, alpha_s, alpha_f;
    label nOptSteps, maxPiggyLoop, maxMMA, maxoutit, maxInnerLoop;
    bool TempSetting;
    List<label> designSpaceCells;
    List<label> solidSpaceCells;
    
    CheckInterface check;
    CheckDict checkDict;
    CheckDatabase checkDB;

public:
    MMAProgram
    (
        Foam::Time& runTime,
        Foam::fvMesh& mesh,
        Foam::simpleControl& simple,
        Foam::volScalarField& p,
        Foam::volVectorField& U,
        //Foam::volVectorField& T,
        Foam::singlePhaseTransportModel& laminarTransport,
        autoPtr<incompressible::turbulenceModel>& turbulence,
        IOMRFZoneList& MRF,
        fv::options& fvOptions,
        volScalarField& alpha,
        volScalarField& eta,
	volScalarField& eta_MMA,
        volScalarField& sens,
        volScalarField& dfdeta,
        surfaceScalarField& phi,
        label&  pRefCell,
        scalar& pRefValue,
        scalar& cumulativeContErr,
        volScalarField& eta_old
    ) :
        runTime(runTime),
        mesh(mesh),
        simple(simple),
        p(p),
        U(U),
        //T(T),
        laminarTransport(laminarTransport),
        turbulence(turbulence),
        MRF(MRF),
        fvOptions(fvOptions),
        alpha(alpha),
        eta(eta),
        eta_MMA(eta_MMA),
        sens(sens),
        dfdeta(dfdeta),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        eta_old(eta_old),

        check(runTime),
    	checkDict(runTime),
    	checkDB(runTime, checkDict)
{
	label tapeSizeMB = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("tapeSizeMB",4096);
	Info << "Creating Tape, size: " << tapeSizeMB << endl;
	AD::createGlobalTape(tapeSizeMB/Pstream::nProcs());
	
	AD::registerInputVariable(eta.begin(),eta.end());

        AD::switchTapeToPassive();
}

    void runLoop() {
	Info<< "\nStarting time loop\n" << endl;
	label iter = 0;
	while (iter < nOptSteps && simple.loop())
	{
	    	iter++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity SIMPLE corrector
		{
		    #include "AlphaDT.H"
		    #include "UEqn.H"
		    if (TempSetting) {
		   	#include "TEqn.H"
		    }
		    #include "pEqn.H"
		}

		laminarTransport.correct();
		turbulence->correct();

		runTime.write();

		runTime.printExecutionTime(Info);
	}

	Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = "   << runTime.elapsedClockTime() << " s"
	<< nl << endl;
    }
    
    void start() {
        #include "settings.H"

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
    }
    
    scalar calcFval() {  
	scalar val;
	forAll(designSpaceCells,i){
	    const label j = designSpaceCells[i];
	    val += eta[j] * mesh.V()[j];
	}
	Foam::reduce(val,sumOp<scalar>());
	val = (val/n) - volumeConstraint;
	    
	Info << "fval = " << val << endl;
	return val;
    }
    
    scalar calcCost(){
        Info << "ACTION::calcCost" << endl;
        return CostFunction(mesh).eval();
    }
       
    void checkpointerLoop() {
	    Info<< "Starting Piggyback Loop: " << "\n" << endl;
	    
	    turbulence->validate();

	    bool firstStep = true;
	    scalar J = 0, oldSensSum = 0;

	    label NN = designSpaceCells.size();
	    Foam::reduce(NN,sumOp<label>());
	    

	    scalar dSensSum = std::numeric_limits<double>::max();
	    
	    label optStep = 0;
	    while ( (dSensSum > optEpsilon || !runTime.writeTime()) && optStep < maxPiggyLoop && simple.loop())
	    {
		    AD::switchTapeToActive();
		    AD::position_t reset_to = AD::getTapePosition();

	    	    checkDB.registerAdjoints(); // store tape indices

		    Info<< "Time = " << runTime.timeName() << nl << endl;

		    // --- Pressure-velocity SIMPLE corrector
		    {
			#include "AlphaDT.H"
			#include "UEqn.H"
			if (TempSetting) {
			    #include "TEqn.H"
			}
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
		    
		    J = CostFunction(mesh).eval();
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
		    scalar norm2 = checkDB.calcNormOfStoredAdjoints();
		    scalar sensSum  = 0.0;

		    forAll(designSpaceCells,i){
			const label j = designSpaceCells[i];
			sens[j] = AD::derivative(eta[j]) / mesh.V()[j];
			sensSum += mag(sens[j]);
		    }

		    Foam::reduce(sensSum,sumOp<scalar>());
		    Foam::reduce(norm2,sumOp<scalar>());
		    dSensSum = mag(sensSum - oldSensSum)/NN;
		    optStep++;
		    Info << "piggy: " << optStep << " " << runTime.timeName() << " " << sensSum << " " << dSensSum << " " << norm2 << " " << J << endl;
		    
		    //Info << "piggy: " << optStep << " " << runTime.timeName() << " " << norm2 << " " << J << endl;
		    oldSensSum = sensSum;

		    AD::resetTapeTo(reset_to);
		    AD::zeroAdjointVector();
		    runTime.write();
		    AD::switchTapeToPassive();	
	    }
	    
	    Info << "avg/max eta,: " << gAverage(eta) << " " << gMax(eta) << endl;

	    //AD::removeTape();
	    Info<< "End\n" << endl;
    }
    
    bool GlobalSolver() {
        MMASolver();
        
    	Info<< "Starting Global Program\n" << endl;
    	scalar J = 0, oldJ = 0.0, fval = 0.0, oldfval = 0.0;
    	
    	GCMMASolver gcmma(mesh,designSpaceCells);

    	scalar ch = 1.0;

   	for (int iter = 0; ch>0.002 && iter < maxoutit; ++iter) {
    	    checkpointerLoop();
	    J = calcCost();
	    fval = calcFval();
	    
	    gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);  
	    eta_MMA.write();

	    eta_old = eta;
	    eta_old.write();
	      
	    oldJ = J;
	    oldfval = fval;

	    eta = eta_MMA;
    	    runLoop();
	    J = calcCost();
	    fval = calcFval();
	    
	    bool conserv = gcmma.ConCheck(J, fval);
	    Info<< "Outer Loop "<< iter << ": " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << "\n" << endl;
	
	    for (int inneriter = 0; !conserv && inneriter < maxInnerLoop; ++inneriter) {
		// Inner iteration update
		gcmma.InnerUpdate(eta_MMA, J, fval, eta, oldJ, sens, oldfval, dfdeta);
		eta_MMA.write();

		eta = eta_MMA;
		runLoop();
		
		J = calcCost();
		fval = calcFval();
		
		Info<< "Outer Loop "<< iter << ": " <<"Inner Sub Loop "<< inneriter << ": " << runTime.timeName() << " cost: "<< J << " fval: " << fval << "\n" << endl;
		
		conserv = gcmma.ConCheck(J, fval);
		//Info<< "Updated Conservative: " << conserv << "\n" << endl;
		if (!conserv) {
	    		eta = eta_old;
	    		J = oldJ;
	    	  	fval = oldfval;
	    	}
	    }
    	    eta = eta_MMA;
	    eta.write();
	    
	    ch = 0.0;
	    forAll(designSpaceCells,i){
	    	const label j = designSpaceCells[i];
		ch = max(ch,mag(eta[j]-eta_old[j]));
	    }
	    Foam::reduce(ch,maxOp<scalar>());
	    eta_old = eta;
    	}
    	return true;
    	Info<< "End\n" << endl;
    }
    
    bool MMASolver() {
    	Info<< "Starting MMA Program\n" << endl;
    	scalar J = 0, oldJ = 0.0, fval = 0.0, oldfval = 0.0;
    	
	start();
	runLoop();
	checkpointerLoop();
	initialGuess();
	GCMMASolver mma(mesh,designSpaceCells);
	
	scalar ch = 1.0;

   	for (int iter = 0; ch>0.002 && iter < maxMMA; ++iter) {
    	    //checkpointerLoop();
	    runLoop();
	    J = calcCost();
	    fval = calcFval();
	    
	    mma.MMAUpdate(eta_MMA, eta, J, sens, fval, dfdeta);  
	    eta_MMA.write();

	    eta_old = eta;
	    eta_old.write();
	      
	    oldJ = J;
	    oldfval = fval;

	    eta = eta_MMA;
    	    runLoop();
	    J = calcCost();
	    fval = calcFval();

    	    eta = eta_MMA;
	    eta.write();
	    
	    ch = 0.0;
	    forAll(designSpaceCells,i){
	    	const label j = designSpaceCells[i];
		ch = max(ch,mag(eta[j]-eta_old[j]));
	    }
	    Foam::reduce(ch,maxOp<scalar>());
	    eta_old = eta;
	    Info<< "MMA Loop "<< iter << ": " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " ch: " << ch <<"\n" << endl;
    	}
    	return true;
    	Info<< "End\n" << endl;
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
         //T,
         laminarTransport,
         turbulence,
         MRF,
         fvOptions,
         alpha,
         eta,
         eta_MMA,
         sens,
         dfdeta,
         phi,
         pRefCell,
         pRefValue,
         cumulativeContErr,
         eta_old
    );
     
    program.GlobalSolver();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

