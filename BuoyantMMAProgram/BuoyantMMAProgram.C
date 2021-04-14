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
#include "radiationModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

//#include "costFunctionLibrary.C"
#include "GCMMASolver.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CheckInterface.H"
#include "CheckDict.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class MMAProgram {

private:
	Foam::Time& runTime;
	Foam::fvMesh& mesh;
	Foam::simpleControl& simple;
	Foam::volScalarField& T;
	Foam::volScalarField& p_rgh;
	Foam::volVectorField& U;
	Foam::surfaceScalarField& phi;
	Foam::singlePhaseTransportModel& laminarTransport;
	const Foam::dimensionedScalar beta;
	const Foam::dimensionedScalar TRef;
	const Foam::dimensionedScalar Pr;
	const Foam::dimensionedScalar Prt;
	const Foam::dimensionedScalar alphaSolid;
	autoPtr<incompressible::turbulenceModel>& turbulence;
	Foam::volScalarField& rhok;
	Foam::volScalarField& alphat;
	Foam::uniformDimensionedVectorField g;
	Foam::volScalarField& gh;
	Foam::surfaceScalarField& ghf;
	Foam::label&  pRefCell;
	Foam::scalar& pRefValue;
	Foam::volScalarField& p;
	Foam::volScalarField& alpha;
	Foam::volScalarField& sens;
	Foam::volScalarField& eta;
	Foam::volScalarField& eta_MMA;
	Foam::volScalarField& eta_old;
	Foam::volScalarField& porosity;
	Foam::volScalarField& dfdeta;
	Foam::IOMRFZoneList& MRF;
	autoPtr<radiation::radiationModel>& radiation;
	const Foam::dimensionedScalar rhoCpRef;
	fv::options& fvOptions;
	Foam::scalar& cumulativeContErr;
	const Foam::wordList costFunctionPatches;
	CheckInterface check;
	CheckDict checkDict;
	CheckDatabase checkDB;
    
    // member variables
    Foam::scalar penalty, volumeConstraint, n, optEpsilon, porosity_s, porosity_f;
    Foam::label nOptSteps, maxPiggyLoop, maxMMA, maxoutit, maxInnerLoop;
    Foam::List<label> designSpaceCells;
    Foam::List<label> solidSpaceCells;
    
    Foam::scalar fluidCp, fluidRho;
    


public:
    MMAProgram
    (
	Foam::Time& runTime,
	Foam::fvMesh& mesh,
	Foam::simpleControl& simple,
	Foam::volScalarField& T,
	Foam::volScalarField& p_rgh,
	Foam::volVectorField& U,
	Foam::surfaceScalarField& phi,
	Foam::singlePhaseTransportModel& laminarTransport,
	Foam::dimensionedScalar beta,
	Foam::dimensionedScalar TRef,
	Foam::dimensionedScalar Pr,
	Foam::dimensionedScalar Prt,
	Foam::dimensionedScalar alphaSolid,
	autoPtr<incompressible::turbulenceModel>& turbulence,
	Foam::volScalarField& rhok,
	Foam::volScalarField& alphat,
	Foam::uniformDimensionedVectorField g,
	Foam::volScalarField& gh,
	Foam::surfaceScalarField& ghf,
	Foam::label&  pRefCell,
	Foam::scalar& pRefValue,
	Foam::volScalarField& p,
	Foam::volScalarField& alpha,
	Foam::volScalarField& sens,
	Foam::volScalarField& eta,
	Foam::volScalarField& eta_MMA,
	Foam::volScalarField& eta_old,
	Foam::volScalarField& porosity,
	Foam::volScalarField& dfdeta,
	Foam::IOMRFZoneList& MRF,
	autoPtr<radiation::radiationModel>& radiation,
	Foam::dimensionedScalar rhoCpRef,
	fv::options& fvOptions,
	Foam::scalar& cumulativeContErr
    ) :
	 runTime(runTime),
	 mesh(mesh),
	 simple(simple),
	 T(T),
	 p_rgh(p_rgh),
	 U(U),
	 phi(phi),
	 laminarTransport(laminarTransport),
	 beta(beta),
	 TRef(TRef),
	 Pr(Pr),
	 Prt(Prt),
	 alphaSolid(alphaSolid),
	 turbulence(turbulence),
	 rhok(rhok),
	 alphat(alphat),
	 g(g),
	 gh(gh),
	 ghf(ghf),
	 pRefCell(pRefCell),
	 pRefValue(pRefValue),
	 p(p),
	 alpha(alpha),
	 sens(sens),
	 eta(eta),
	 eta_MMA(eta_MMA),
	 eta_old(eta_old),
	 porosity(porosity),
	 dfdeta(dfdeta),
	 MRF(MRF),
	 radiation(radiation),
	 rhoCpRef(rhoCpRef),
	 fvOptions(fvOptions),
	 cumulativeContErr(cumulativeContErr),
	 costFunctionPatches(mesh.solutionDict().subDict("SIMPLE").lookup("costFunctionPatches")),
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
		    #include "porosity.H"
		    #include "UEqn.H"
		    #include "TEqn.H"
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
    Info << "Eta Set ! " << gSum(eta)/n << endl;
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
    
    scalar calcPrCost(){
        Info << "ACTION::calcPrCost" << endl;
	scalar Jp = 0;
	forAll(costFunctionPatches,cI)
	{
	Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI] );
	const Foam::fvPatch& patch = mesh.boundary()[patchI];
	Jp += gSum
		(
		    - phi.boundaryField()[patchI]*(p.boundaryField()[patchI]
		    + 0.5*magSqr(phi.boundaryField()[patchI]/patch.magSf()))
		);
	}
	Info<< "cost pressure: " << Jp << endl;
	return Jp;
    }
    
    scalar calcTempCost(){
        Info << "ACTION::calcTempCost" << endl;
        scalar Jt = 0;
	forAll(costFunctionPatches,cI)
	{
	Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
	Jt += gSum
		(
		    phi.boundaryField()[patchI]*(T.boundaryField()[patchI])
		);
	}
	Jt = rhoCpRef.value()*Jt;
	Info<< "cost temp: " << Jt << endl;
	return Jt;
    }
    
    scalar calcCost(){
        scalar J = 0;
	J = calcTempCost() + calcPrCost();
	return J;
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
			#include "porosity.H"
			#include "UEqn.H"
			#include "TEqn.H"
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
	//start();
   	//runLoop();
   	//checkpointerLoop();
    	Info<< "Starting Global Program\n" << endl;
    	scalar J = 0, oldJ = 0.0, fval = 0.0, oldfval = 0.0;
    	
    	GCMMASolver gcmma(mesh,designSpaceCells);

    	scalar ch = 1.0;
	eta_old = eta;
	eta_old.write();
	
	eta = eta_MMA;
    	runLoop();
	oldJ = calcCost();
	oldfval = calcFval();
   	
   	for (int iter = 0; ch>0.002 && iter < maxoutit; ++iter) {
   	    //MMASolver();
    	    checkpointerLoop();
	    J = calcCost();
	    fval = calcFval();
	    
	    gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);  
	    //gcmma.InnerUpdate(eta_MMA, J, fval, eta, oldJ, sens, oldfval, dfdeta);
	    eta_MMA.write();

	    eta = eta_MMA;
    	    runLoop();
	    J = calcCost();
	    fval = calcFval();

	    bool conserv = gcmma.ConCheck(J, fval);
	    if (!conserv) {
		eta = eta_old;
		J = oldJ;
		fval = oldfval;
	    }
	    Info<< "Outer Loop "<< iter << " Start: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << "\n" << endl;
	
	    for (int inneriter = 0; !conserv && inneriter < maxInnerLoop; ++inneriter) {
		// Inner iteration update
		gcmma.InnerUpdate(eta_MMA, J, fval, eta, oldJ, sens, oldfval, dfdeta);
		//gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);
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

	    ch = 0.0;
	    forAll(designSpaceCells,i){
	    	const label j = designSpaceCells[i];
		ch = max(ch,mag(eta[j]-eta_old[j]));
	    }
	    Foam::reduce(ch,maxOp<scalar>());
	    eta_old = eta;
	    eta_old.write();
	    oldJ = J;
	    oldfval = fval;
	    
	    runLoop();
	    J = calcCost();
	    fval = calcFval();
    	    
    	    Info<< "Outer Loop "<< iter << " End: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << " ch: " << ch << "\n" << endl;
    	}
    	return true;
    	Info<< "End\n" << endl;
    }
    
    bool MMASolver() {
    	Info<< "Starting MMA Program\n" << endl;
    	scalar J = 0, oldJ = 0.0, fval = 0.0, oldfval = 0.0;
    	
	start();
	initialGuess();
	runLoop();
	checkpointerLoop();

	GCMMASolver mma(mesh,designSpaceCells);
	
	scalar ch = 1.0;

   	for (int iter = 0; ch>0.002 && iter < maxMMA; ++iter) {
	    //checkpointerLoop();
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
	    eta_old.write();
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
        "Steady-state solver for buoyant, turbulent flow"
        " of incompressible fluids."
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
         T,
         p_rgh,
         U,
         phi,
         laminarTransport,
         beta,
         TRef,
         Pr,
         Prt,
         alphaSolid,
         turbulence,
         rhok,
         alphat,
         g,
         gh,
         ghf,
         pRefCell,
         pRefValue,
         p,
         alpha,
         sens,
         eta,
         eta_MMA,
         eta_old,
         porosity,
         dfdeta,
         MRF,
         radiation,
         rhoCpRef,
         fvOptions,
         cumulativeContErr
    );
    
    program.MMASolver();
    //program.start();
    //program.initialGuess();
    //program.runLoop();
    
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

