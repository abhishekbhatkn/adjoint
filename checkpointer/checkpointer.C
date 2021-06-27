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

#include "profiling.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "CheckInterface.H"
#include "CheckController.H"

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
    volScalarField& porosity;
    volScalarField& eta;
    volScalarField& dfdeta;
    volScalarField& sens;
    surfaceScalarField& phi;

    label&  pRefCell;
    scalar& pRefValue;
    scalar& cumulativeContErr;
    const Foam::wordList costFunctionPatches;

    // member variables
    bool frozenTurbulence, forwardLoop;
    scalar penalty, volumeConstraint, designVolume, optEpsilon, porosity_s, porosity_f, factorP, asyminit, asymdec, asyminc;
    label nOptSteps, maxPiggyLoop, designSize, maxMMAiter, MMALoop, iter;
    List<label> designSpaceCells;
    List<label> solidSpaceCells;


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
        volScalarField& porosity,
        volScalarField& eta,
        volScalarField& dfdeta,
        volScalarField& sens,
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
        porosity(porosity),
        eta(eta),
        dfdeta(dfdeta),
        sens(sens),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        costFunctionPatches(mesh.solutionDict().subDict("SIMPLE").lookup("costFunctionPatches")),
         forwardLoop(true),
	 iter(0)
    {
    	#include "settings.H"
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
        Foam::volScalarField pold("pold",p);
        Foam::volVectorField Uold("Uold",U);
        {
            #include "porosity.H"
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
        
        if (forwardLoop) {
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
			runTime.writeNow();
			forwardLoop = false;
			return true;
		}
		iter++;
	}
        return false;
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
    
    void postInterpret(){
        //Info << "Writing sensitivity for timestep " << runTime.timeName() << endl;
        for(int i = 0; i < eta.size(); i++)
            sens[i] = AD::derivative(eta[i]) / mesh.V()[i];
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
                for(int i = 0; i < eta.size(); i++){
                    sens[i] = AD::derivative(eta[i]) / mesh.V()[i];
                }
                sens.write();
                //runTime.writeNow();
            }
            //Info << "sens Sum: " << runTime.timeName() << "\t" << gSum(sens) << endl;
        }
    }

    void start(){

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
         porosity,
         eta,
         dfdeta,
         sens,
         phi,
         pRefCell,
         pRefValue,
         cumulativeContErr
    );
    
    if( readLabel(mesh.solutionDict().subDict("SIMPLE").lookup("MMALoop"))== 0) {
    	simpleCheck.runLoop();
    }

    dimensionedScalar startTime = runTime.value();
    label maxPiggyLoop = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("maxPiggyLoop",100);
    runTime.setEndTime(startTime+maxPiggyLoop);
    simpleCheck.run();

    runTime.setTime(startTime,0);
    sens.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Foam::profiling::print(Info);

    return 0;
}


// ************************************************************************* //
