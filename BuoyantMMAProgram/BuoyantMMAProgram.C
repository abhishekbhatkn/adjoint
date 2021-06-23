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
	Foam::volScalarField& sens;
	Foam::volScalarField& dfdeta;
	Foam::volScalarField& eta;
	Foam::volScalarField& eta_old1;
	Foam::volScalarField& eta_old2;
	Foam::volScalarField& low;
	Foam::volScalarField& upp;
	Foam::volScalarField& porosity;
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
    Foam::scalar penalty, volumeConstraint, designVolume, optEpsilon, porosity_s, porosity_f, factorT, factorP, refRho, Cp, minTempCost, maxTempCost, minPrCost, maxPrCost, asyminit, asymdec, asyminc;
    Foam::label nOptSteps, maxPiggyLoop, designSize, maxMMAiter, MMALoop;
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
	Foam::volScalarField& dfdeta,
	Foam::volScalarField& eta,
	Foam::volScalarField& eta_old1,
	Foam::volScalarField& eta_old2,
	Foam::volScalarField& low,
	Foam::volScalarField& upp,
	Foam::volScalarField& porosity,
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
	 sens(sens),
	 dfdeta(dfdeta),
	 eta(eta),
	 eta_old1(eta_old1),
	 eta_old2(eta_old2),
	 low(low),
	 upp(upp),
	 porosity(porosity),
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
	#include "settings.H"
	bool MMARun = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<bool>("MMARun",false);
	if (!MMARun) {
		label tapeSizeMB = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("tapeSizeMB",4096);
		Info << "Creating Tape, size: " << tapeSizeMB << endl;
		AD::createGlobalTape(tapeSizeMB/Pstream::nProcs());

		AD::switchTapeToPassive();
	}
}

    void runLoop() {
	Info<< "\nStarting time loop\n" << endl;
	label iter = 0;
	bool check = true;
	while (iter < nOptSteps && check && simple.loop())
	{
	    	iter++;
		Info<< "Time = " << runTime.timeName() << nl << endl;
		
		Foam::volScalarField pold("pold",p), Told("Told",T);
		Foam::volVectorField Uold("Uold",U);
		// --- Pressure-velocity SIMPLE corrector
		{
		    #include "porosity.H"
		    #include "UEqn.H"
		    #include "TEqn.H"
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
			maxdiff = max(maxdiff,mag(Told[i]-T[i]));
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
	//value[1] = 0.0; //((val/designVolume)-(volumeConstraint));


	Info << "fval= " << value << endl;
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
	Jt = Dt;
	Info<< "cost temp: " << Jt << endl;
	return Jt;
    }
    
    scalar calcCost(){
        scalar J = 0;
	J = -factorT*(calcTempCost()-minTempCost)/(maxTempCost-minTempCost) + factorP*(calcPrCost()-minPrCost)/(maxPrCost-minPrCost);
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

		    for(label i = 0; i < eta.size(); i++) {
			sens[j] = AD::derivative(eta[i]) * designVolume / (designSize*mesh.V()[i]);
			sensSum += mag(sens[i]);
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
    
    bool GradientSolver() {

    	Info<< "Starting Gradient Solver\n" << endl;
    	scalar lam = 1, oldJ = 0.0, J = 0.0;
    	runLoop();
	scalar check = 1.0;
	while ( MMALoop < maxMMAiter && check > optEpsilon) {
	   	checkpointerLoop();
		J = calcCost();
		List<scalar> fval = calcFval();
		if (MMALoop == 0) {
			Info<< " Start: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << "\n" << endl; 
			runTime.writeNow();  
		}
		eta_old1 = eta;
		forAll(designSpaceCells,i){
			const label j = designSpaceCells[i];
			if (mag(sens[j]) > optEpsilon) {
				eta[j] = eta[j] - sens[j]*lam;
				eta[j] = max(0.0, eta[j]);
				eta[j] = min(1.0, eta[j]);
			}
		}
		runTime.writeNow();
		check = etaCheck();
		if (check > optEpsilon) {
			oldJ = J;
		    	runLoop();
			J = calcCost();
			fval = calcFval();
			label loop = 0;
			while (J > oldJ && check > optEpsilon) {
				eta = eta_old1;
				lam = lam/2;
				forAll(designSpaceCells,i){
					const label j = designSpaceCells[i];
					if (mag(sens[j]) > optEpsilon) {
						eta[j] = eta[j] - sens[j]*lam;
						eta[j] = max(0.0, eta[j]);
						eta[j] = min(1.0, eta[j]);
					}
				}
				runTime.writeNow();
				check = etaCheck();
				if (check > optEpsilon) {
					runLoop();
					Info << "Outer Loop "<< MMALoop << " Inner Loop "<< loop << " Time: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " lam: " << lam << endl;
					++loop;
				}
			}
		}
	    	Info<< "Outer Loop "<< MMALoop << " End: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " lam: " << lam << endl;
	    	++MMALoop;
	    	runTime.writeNow();
    	}
	Info<< "Gradient Solver End\n" << endl;
    	return true;
    }
    
    bool GlobalSolver() {

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
/*    
    bool GradientSolver() {

    	Info<< "Starting Gradient Solver\n" << endl;
    	scalar lam = 1, oldJ = 0.0, J = 0.0;
    	runLoop();
	scalar check = 1.0;
	while ( MMALoop < maxMMAiter) {
	   	checkpointerLoop();
		J = calcCost();
		List<scalar> fval = calcFval();
		if (MMALoop == 0) {
			Info<< " Start: " << runTime.timeName() << " cost: "<< J << " fval: " << fval << "\n" << endl; 
			runTime.writeNow();  
		}
		eta_old1 = eta;
		forAll(designSpaceCells,i){
			const label j = designSpaceCells[i];
			if (mag(sens[j]) > optEpsilon) {
				eta[j] = eta[j] - sens[j]*lam;
				eta[j] = max(0.0, eta[j]);
				eta[j] = min(1.0, eta[j]);
			}
		}
		runTime.writeNow();
		check = etaCheck();
		if (check > optEpsilon) {
			oldJ = J;
		    	runLoop();
			J = calcCost();
			fval = calcFval();
			label loop = 0;
			while (J > oldJ && check > optEpsilon) {
				eta = eta_old1;
				factorP = factorP/2;
				checkpointerLoop();
				J = calcCost();
				oldJ = J;
				forAll(designSpaceCells,i){
					const label j = designSpaceCells[i];
					if (mag(sens[j]) > optEpsilon) {
						eta[j] = eta[j] - sens[j]*lam;
						eta[j] = max(0.0, eta[j]);
						eta[j] = min(1.0, eta[j]);
					}
				}
				runTime.writeNow();
				check = etaCheck();
				if (check > optEpsilon) {
					runLoop();
					J = calcCost();
					fval = calcFval();
					Info << "Outer Loop "<< MMALoop << " Inner Loop "<< loop << " Time: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " factorP: " << factorP << endl;
					++loop;
				}
			}
		}
	    	Info<< "Outer Loop "<< MMALoop << " End: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " lam: " << lam << endl;
	    	++MMALoop;
	    	runTime.writeNow();
    	}
	Info<< "Gradient Solver End\n" << endl;
    	return true;
    }
 */   
        
    bool MMASolver() {
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
	scalar J = calcCost();
	List<scalar> fval = calcFval();
	gcmma.MMAUpdate(eta, eta_old1, eta_old2, J, sens, fval, dfdeta);
	scalar check = etaCheck();
	Info<< "Loop "<< MMALoop << " Time: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " asyminit: " << asyminit << " asymdec: " << asymdec << " asyminc: " << asyminc << endl;
	runTime.writeNow();
	if (check < optEpsilon) {
		return false;
	}
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
    
    void setTemp() {
    	forAll(designSpaceCells,i){
		const label j = designSpaceCells[i];
		if (eta[j] < 0.1) {
			T[j] = 350;
		}
	}
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
	 dfdeta,
         eta,
         eta_old1,
         eta_old2,
         low,
         upp,
         porosity,
         MRF,
         radiation,
         rhoCpRef,
         fvOptions,
         cumulativeContErr
    );

    bool MMARun = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<bool>("MMARun",false);
    
    if (!MMARun) {
    	program.setTemp();
    	program.runLoop();
    	program.checkpointerLoop();
    	runTime.writeNow();
    }
    else {
    	program.MMASolver();
    }
    //program.GlobalSolver();
    //program.GradientSolver();
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

