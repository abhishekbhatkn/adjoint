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
/*
#include "profiling.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "CheckInterface.H"
#include "CheckController.H"
#include "costFunctionLibrary.C"

#include "GCMMASolver.C"

std::vector<double> tapeMem;

class simpleCheckpointer : public CheckController{
private:
    // inherited from CheckController:
    // Foam::Time& runTime;
    // autoPtr<CheckInterface> checkInterface;

    // references for variables living in main
    Foam::fvMesh& mesh;
    Foam::simpleControl& simple;
    Foam::volScalarField& p;
    Foam::volVectorField& U;
    Foam::singlePhaseTransportModel& laminarTransport;
    autoPtr<incompressible::turbulenceModel>& turbulence;
    Foam::IOMRFZoneList& MRF;
    fv::options& fvOptions;
    volScalarField& alpha;
    volScalarField& sens;
    volScalarField& eta;
    volScalarField& eta_MMA;
    volScalarField& dfdeta;
    surfaceScalarField& phi;

    label&  pRefCell;
    scalar& pRefValue;
    scalar& cumulativeContErr;

    // member variables
    scalar J, penalty, volumeConstraint, fval;
    dimensionedScalar alpha_f, alpha_s;
    bool frozenTurbulence;
    List<label> designSpaceCells;

public:
    simpleCheckpointer
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
        volScalarField& alpha,
        volScalarField& sens,
        volScalarField& eta,
        volScalarField& eta_MMA,
        volScalarField& dfdeta,
        surfaceScalarField& phi,
        label&  pRefCell,
        scalar& pRefValue,
        scalar& cumulativeContErr
    )
    :
        CheckController(runTime),
        mesh(mesh),
        simple(simple),
        p(p),
        U(U),
        laminarTransport(laminarTransport),
        turbulence(turbulence),
        MRF(MRF),
        fvOptions(fvOptions),
        alpha(alpha),
        sens(sens),
        eta(eta),
        eta_MMA(eta_MMA),
        dfdeta(dfdeta),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        alpha_f("alpha_f", dimless/dimTime, 0.0),
        alpha_s("alpha_s", dimless/dimTime,0.0)
        //designSpaceCells(0)

    {
        Info << "Avg Checkpoint Size: "
             << checkInterface().checkDatabase().getCheckDatabaseSize()
             << " MB" << endl;
        Info << "nCheckpoints: " << checkInterface().checkDatabase().nCheckpoints() << endl;

        //CheckObjectScalar* sObj = new CheckObjectScalar(cumulativeContErr);
        checkInterface().checkDatabase().addScalarCheckpoint(cumulativeContErr);
        checkInterface().checkDatabase().addDictionaryCheckpoint(const_cast<Foam::dictionary&>(mesh.solverPerformanceDict()));

        label tapeSizeMB = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("tapeSizeMB",4096);
        Info << "Creating Tape, size: " << tapeSizeMB << endl;
        AD::createGlobalTape(tapeSizeMB/Pstream::nProcs());

        //AD::registerInputVariable(eta.begin(),eta.end());
	AD::registerInputVariable(eta.begin(),eta.end());

        AD::switchTapeToPassive();
        frozenTurbulence = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<bool>("frozenTurbulence",false);
        if(frozenTurbulence){
            Info << "\nWARNING: Calculating Adjoints with frozen Turbulence assumption\n" << endl;
        }
    }

    bool runStep(){
        static bool needsWrite = true;
        bool finished = !simple.loop();
        Info << "stopAt: " << runTime.stopAt() << " " << runTime.endTime() << endl;

        if(finished){
            needsWrite = false;
            return true;
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            alpha = alpha_s+(alpha_f-alpha_s)*eta*(1.0+penalty)/(eta + penalty);
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        bool wasActive = AD::isTapeActive();
        if(frozenTurbulence && wasActive)
            AD::switchTapeToPassive();

        turbulence->correct();

        if(frozenTurbulence && wasActive)
            AD::switchTapeToActive();

        if (needsWrite)
            runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = "   << runTime.elapsedClockTime() << " s"
             << nl << endl;

        return false;
    }

    scalar calcCost(){
        Info << "ACTION::calcCost" << endl;
        return CostFunction(mesh).eval();
    }

    void postInterpret(){
        //Info << "Writing sensitivity for timestep " << runTime.timeName() << endl;
        forAll(designSpaceCells,i){
            const label j = designSpaceCells[i];
        //sens[j] = AD::derivative(alpha[j]) / mesh.V()[j];    
	sens[j] = AD::derivative(eta[j]) / mesh.V()[j];
        }
        Info << "sens Sum: " << runTime.timeName() << "\t" << gSum(sens) << endl;

        //runTime.stopAt(Foam::Time::stopAtControls::saEndTime);
        if(runTime.writeTime()){
            Info << "Post Interpret: Writing sensitivity for timestep " << runTime.timeName() << endl;
            sens.write();
        }
        //runTime.write();
    }
    void write(bool firstRun){
        if(firstRun)
            runTime.write();
        else{
            if(runTime.writeTime() || runTime.stopAt() == Foam::Time::stopAtControls::saWriteNow){
                Info << "Writing sensitivity for timestep " << runTime.timeName() << endl;
            forAll(designSpaceCells,i){
                const label j = designSpaceCells[i];
                sens[j] = AD::derivative(eta[j]) / mesh.V()[j];
            }
                sens.write();
            }
            //Info << "sens Sum: " << runTime.timeName() << "\t" << gSum(sens) << endl;
        }
    }

    void start(){
	    dictionary simpleDict = mesh.solutionDict().subDict("SIMPLE");

	    autoPtr<wordList> costFunctionPatches;
	    costFunctionPatches.reset(new wordList(0));

	    if (simpleDict.found("costFunctionPatches"))
	    {
		costFunctionPatches.reset
		(
		    new wordList(simpleDict.lookup("costFunctionPatches"))
		);
	    }else{
		Info << "Warning: Keyword costFunctionPachtes not found in fvSolution/SIMPLE" << endl;
		Info << "Example:\nSIMPLE{\n  costFunctionPatches (inlet outlet);\n}" << endl;
	    }

	    const Foam::wordList designSpaceZones(
		    simpleDict.lookupOrDefault<Foam::wordList>("designSpace",Foam::wordList())
	    );
	    if(designSpaceZones.size()>0){
		forAll(designSpaceZones, i){
		    const label cellZoneID = mesh.cellZones().findZoneID(designSpaceZones[i]);
		    designSpaceCells.append( mesh.cellZones()[cellZoneID] );
		}
	    }else{ // add all cells
		forAll(eta, i){
		    designSpaceCells.append( i );
		}
	    }

	    volumeConstraint = simpleDict.lookupOrDefault<scalar>("volumeConstraint",0.4);
	    if (!simpleDict.found("volumeConstraint"))
	    {
		Info << "Warning: Keyword volumeConstraint not found in fvSolution/SIMPLE. Default set to 0.4" << endl;
	    }
	    penalty = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("penalty",0.6);
	    if (!mesh.solutionDict().subDict("SIMPLE").found("penalty"))
	    {
		Info << "Warning: Keyword penalty not found in fvSolution/SIMPLE. Default set to 0.6" << endl;
	    }
	    alpha_f =dimensionedScalar("alpha_f", dimless/dimTime, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_f",0.0));
	    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_f"))
	    {
		Info << "Warning: Keyword alpha_f not found in fvSolution/SIMPLE. Default set to 0" << endl;
	    }
	    alpha_s = dimensionedScalar("alpha_s", dimless/dimTime, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_s",10e5));
	    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_s"))
	    {
		Info << "Warning: Keyword alpha_s not found in fvSolution/SIMPLE. Default set to 10e5" << endl;
	    }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    simpleCheckpointer simpleCheck
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
         alpha,
         eta,
         sens,
         dfdeta,
         phi,
         pRefCell,
         pRefValue,
         cumulativeContErr
    );

    simpleCheck.start();
    simpleCheck.run();

    scalar fval = 0;
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
		fval += eta[j] - 1.0; 
    }
    Foam::reduce(fval,sumOp<scalar>());
    label n = designSpaceCells.size();
    Foam::reduce(n,sumOp<label>());
    fval = (fval/n) - volumeConstraint;

    AD::derivative(fval) = 1.0;
    forAll(designSpaceCells,i){
            const label j = designSpaceCells[i];
      	    dfdeta[j] = AD::derivative(eta[j]) / mesh.V()[j];
    }

    Info << "fval: " << fval << " n " << n << endl;

    //runTime.setTime(0,0);
    sens.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Foam::profiling::print(Info);
    
    //GCMMASolver gcmma(mesh,designSpaceCells);
	
    //gcmma.OuterUpdate(eta_MMA, eta, cost, sens, fval, dfdeta);
    runTime.write();

    return 0;
}

*/
//---------------------------------------------------------------------------//


#include "profiling.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "CheckInterface.H"
#include "CheckDict.H"
#include "cellSet.H"

#include "costFunctionLibrary.H"
#include "GCMMASolver.C"

// ************************************************************************* //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    //#include "createMRF.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #include "adjointSettings.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    CheckInterface check(runTime);
    CheckDict checkDict(runTime);
    CheckDatabase checkDB(runTime, checkDict);

    scalar J = 0.0, fval = 0.0;
    scalar oldSensSum = 0;

    Info<< "\nStarting time loop\n" << endl;

    //#include "Random.H"
    //Foam::Random randObj(time(NULL));
    //forAll(designSpaceCells,i){
    //	const label j = designSpaceCells[i];
    //   	eta[j] = randObj.bit();
    //}

    label tapeSizeMB = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("tapeSizeMB",4096);
    Info << "Creating Tape, size: " << tapeSizeMB << endl;
    AD::createGlobalTape(tapeSizeMB/Pstream::nProcs());

    bool firstStep = true;

    turbulence->validate();

    label NN = designSpaceCells.size();

    Foam::reduce(NN,sumOp<label>());
    
    label nOptSteps = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("nOptSteps",100);
    scalar optEpsilon = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("optTolerance",5e-2);

    for(int optStep = 0; optStep < nOptSteps; optStep++){
        //AD::global_tape->reset();
        AD::switchTapeToActive();

        forAll(eta,i){
            AD::registerInputVariable(eta[i]);
        }

        AD::position_t reset_to = AD::getTapePosition();
        scalar dSensSum = std::numeric_limits<double>::max();
	
        while (dSensSum > optEpsilon && simple.loop())
        {
            checkDB.registerAdjoints(); // store tape indices

            Info<< "Time = " << runTime.timeName() << nl << endl;

            // --- Pressure-velocity SIMPLE corrector
            {
		alpha=alpha_s+(alpha_f-alpha_s)*eta*(1.0+penalty)/(eta + penalty);
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
            Info << "piggy: " << optStep << " " << runTime.timeName() << " " << sensSum << " " << dSensSum << " " << norm2 << " " << J << endl;
            
            //Info << "piggy: " << optStep << " " << runTime.timeName() << " " << norm2 << " " << J << endl;
            oldSensSum = sensSum;

            AD::resetTapeTo(reset_to);
            AD::zeroAdjointVector();
            
//--------------------------------------------------------------------------------------------------------------------------------//            
            
            //label n = designSpaceCells.size();
            forAll(designSpaceCells,i){
		const label j = designSpaceCells[i];
		if (eta[j]<0.5)
			fval += 1.0;
     	    }
	    Foam::reduce(fval,sumOp<scalar>());
	    //Foam::reduce(n,sumOp<label>());
            fval = (fval/NN) - volumeConstraint;
            
            forAll(designSpaceCells,i){
                const label j = designSpaceCells[i];
                dfdeta[j] = 1.0/NN;
            }     
            
	    Info << "fval: " << fval << " n " << NN << endl;
            runTime.write();
	    //sens.write();
	    //eta.write();
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = "   << runTime.elapsedClockTime() << " s"
                << nl << endl;
/*
            if(fval < 0.0 || fval > (1-volumeConstraint)) {
		forAll(designSpaceCells,i){
		    const label j = designSpaceCells[i];
		    eta[j] = sens[j];
		    eta[j] =  max(eta[j],scalar(0.0));
		    eta[j] =  min(eta[j],scalar(1.0));
		}
	     }
*/
        }
        AD::switchTapeToPassive();

	GCMMASolver gcmma(mesh,designSpaceCells);
	
	//sens = sens/gAverage(sens);
	//Info << "new avg/max sens: " << gAverage(sens) << " " << gMax(sens) << endl;

	gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);
	
	eta_MMA.write();
	eta = eta_MMA;
    }
    runTime.write();
    Info << "avg/max eta,: " << gAverage(eta) << " " << gMax(eta) << endl;
    Foam::profiling::print(Info);
    //AD::tape_t::remove(AD::global_tape);
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
