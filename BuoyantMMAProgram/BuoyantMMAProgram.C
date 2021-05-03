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
	Foam::volVectorField U_old;
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
    Foam::scalar penalty, volumeConstraint, n, optEpsilon, porosity_s, porosity_f, factorT, factorP, refRho, Cp;
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
	 U_old(U),
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
        initialGuess();

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
    eta.write();
    
    Info << "Eta Set ! " << endl;

    }
    
    List<scalar> calcFval() {  
	scalar val = 0;
	List<scalar> value(2,0.0);
	forAll(designSpaceCells,i){
	    const label j = designSpaceCells[i];
	    val += eta[j] * mesh.V()[j];
	}
	Foam::reduce(val,sumOp<scalar>());
	value[0] = volumeConstraint-(val/n);
	value[1] = (val/n)-volumeConstraint;

	Info << "fval1 = " << value[0] << " fval2 = " << value[1] << endl;
	return value;
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
        scalar Jt = 0, Dt = 0;
	forAll(costFunctionPatches,cI)
	{
	Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
	Dt += gSum
		(
		    phi.boundaryField()[patchI]*(T.boundaryField()[patchI])
		);
	}
	Jt = refRho*Cp*Dt;
	Info<< "cost temp: " << Jt <<" Delta T: : "<< Dt << endl;
	return Jt;
    }
    
    scalar calcCost(){
        scalar J = 0;
	J = -factorT*(calcTempCost()-30)/(5-30) + factorP*(calcPrCost()-0.0006)/(0.000006-0.0006);
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
	    scalar maxSens = 0, minSens = 0;
	    
	    label optStep = 0;
	    while ( (dSensSum > optEpsilon || !runTime.writeTime()) && optStep < maxPiggyLoop && simple.loop())
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
			maxSens = max(maxSens, sens[j]);
			minSens = min(minSens, sens[j]);
		    }
		    
		    Foam::reduce(sensSum,sumOp<scalar>());
		    Foam::reduce(maxSens,maxOp<scalar>());
		    Foam::reduce(minSens,minOp<scalar>());
		    
/*		    forAll(designSpaceCells,i){
			const label j = designSpaceCells[i];
			sens[j] = 1000*(sens[j])/((maxSens-minSens));
		    }
*/		    Foam::reduce(norm2,sumOp<scalar>());
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

        //MMASolver();
	start();
   	//runLoop();
   	//checkpointerLoop();
    	Info<< "Starting Global Program\n" << endl;
    	scalar J = 0, oldJ = 0.0;
    	List<scalar> fval, oldfval;
    	//GCMMASolver gcmma(mesh,designSpaceCells,2);

    	scalar senscheck = 1.0;
	eta_old = eta;
	U_old = U;
	eta_old.write();
	U_old.write();
	
	GCMMASolver gcmma(mesh,designSpaceCells,2);
   	
   	for (int iter = 0; senscheck>0.002 && iter < maxoutit; ++iter) {
   	    runLoop();
   	    //runLoop();
   	    checkpointerLoop();
	    J = calcCost();
	    fval = calcFval();
	    if (iter == 0) {
   	    	oldJ = J;
		oldfval = fval;
	    }   
	    //bool conserv = gcmma.ConCheck(J, fval);
	    gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);  
	    //gcmma.InnerUpdate(eta_MMA, J, fval, eta, oldJ, sens, oldfval, dfdeta);
	    //eta_MMA = eta;
	    eta_MMA.write();
	    bool conserv = gcmma.ConCheck(J, fval);

/*	    eta = eta_MMA;
    	    runLoop();
	    J = calcCost();
	    fval = calcFval();

	    
	    if (!conserv) {
		eta = eta_old;
		U = U_old;
		runLoop();
		J = calcCost();
		fval = calcFval();
	    }*/
	    Info<< "Outer Loop "<< iter << " Start: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << "\n" << endl;
	
	    for (int inneriter = 0; !conserv && inneriter < maxInnerLoop; ++inneriter) {
		// Inner iteration update
		//gcmma.InnerUpdate(eta_MMA, J, fval, eta, oldJ, sens, oldfval, dfdeta);
		//gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);
		gcmma.MMAUpdate(eta_MMA, eta, J, sens, fval, dfdeta);
		eta_MMA.write();

		eta = eta_MMA;
		runLoop();
		J = calcCost();
		fval = calcFval();
		
		conserv = gcmma.ConCheck(J, fval);
/*		
		//Info<< "Updated Conservative: " << conserv << "\n" << endl;
		if (!conserv) {
			eta = eta_old;
			U = U_old;
			runLoop();
			J = calcCost();
			fval = calcFval();
		}*/
		Info<< "Outer Loop "<< iter << ": " <<"Inner Sub Loop "<< inneriter << ": " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << "\n" << endl;
	    }

	    senscheck = 0.0;
	    forAll(designSpaceCells,i){
	    	const label j = designSpaceCells[i];
		senscheck = max(senscheck,sens[j]);
	    }
	    Foam::reduce(senscheck,maxOp<scalar>());
	    eta_old = eta;
	    U_old = U;
	    eta_old.write();
	    U_old.write();
	    oldJ = J;
	    oldfval = fval;
	    
    	    Info<< "Outer Loop "<< iter << " End: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << " senscheck: " << senscheck << "\n" << endl;
    	}
	//MMASolver();
	Info<< "GLobal Solver End\n" << endl;
	fval = calcCost();
    	return true;
    }
    
    bool MMASolver() {
    	Info<< "Starting MMA Program\n" << endl;
    	scalar J = 0;
    	List<scalar> fval;
	start();
	runLoop();
	//checkpointerLoop()
	eta = eta_old;
	eta_old.write();
	scalar ch = 1.0;

   	for (label iter = 0; ch>0.00002 && iter < maxMMA; ++iter) {
	    checkpointerLoop();
	    GCMMASolver mma(mesh,designSpaceCells,2);
	    J = calcCost();
	    fval = calcFval();
	    bool conserv = false;
	    
	    for (label inneriter = 0; !conserv && inneriter < maxInnerLoop; ++inneriter) {
		    mma.MMAUpdate(eta_MMA, eta, J, sens, fval, dfdeta);  
		    eta_MMA.write();
		Info<< "Outer Loop Start "<<  endl;
		    eta = eta_MMA;
		    eta.write();
	    	    runLoop();
		    J = calcCost();
		    fval = calcFval();
		    
		    conserv = mma.ConCheck(J, fval);
		    Info<< "Outer Loop "<< iter << ": " <<"Inner Sub Loop "<< inneriter << ": " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " conserv: " << conserv << "\n" << endl;
	    }
	    ch = 0.0;
	    forAll(designSpaceCells,i){
	    	const label j = designSpaceCells[i];
		ch = max(ch,mag(eta[j]-eta_old[j]));
	    }
	    Foam::reduce(ch,maxOp<scalar>());
	    eta_old = eta;
	    eta_old.write();
	    Info<< "Outer Loop "<< iter << ": End at: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << " ch: " << ch <<"\n" << endl;
    	}
    	Info<< "MMA Solver End\n" << endl;
    	return true;
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
    
    program.GlobalSolver();
    //program.start();
    //program.runLoop();
    //program.MMASolver();
    
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

