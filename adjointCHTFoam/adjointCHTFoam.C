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
#include "radiationModel.H"
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
        //Foam::Time& runTime;
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
	Foam::volScalarField& eta;
	Foam::volScalarField& porosity;
	Foam::IOMRFZoneList& MRF;
	autoPtr<radiation::radiationModel>& radiation;
	const Foam::dimensionedScalar rhoCpRef;
	fv::options& fvOptions;
	Foam::scalar& cumulativeContErr;
	const Foam::wordList costFunctionPatches;

    // member variables
    bool frozenTurbulence;
    Foam::scalar penalty, volumeConstraint, designVolume, optEpsilon, porosity_s, porosity_f, factorT, factorP, refRho, Cp, minTempCost, maxTempCost, minPrCost, maxPrCost, asyminit, asymdec, asyminc;
    Foam::label nOptSteps, maxPiggyLoop, designSize, maxMMAiter, MMALoop;
    Foam::List<label> designSpaceCells;
    Foam::List<label> solidSpaceCells;
    
    Foam::scalar fluidCp, fluidRho;


public:
    simpleCheckpointer
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
	Foam::volScalarField& porosity,
	Foam::IOMRFZoneList& MRF,
	autoPtr<radiation::radiationModel>& radiation,
	Foam::dimensionedScalar rhoCpRef,
	fv::options& fvOptions,
	Foam::scalar& cumulativeContErr
    )
    :
	 CheckController(runTime),
	 //runTime(runTime),
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
	 eta(eta),
	 porosity(porosity),
	 MRF(MRF),
	 radiation(radiation),
	 rhoCpRef(rhoCpRef),
	 fvOptions(fvOptions),
	 cumulativeContErr(cumulativeContErr),
	 costFunctionPatches(mesh.solutionDict().subDict("SIMPLE").lookup("costFunctionPatches"))
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
            #include "UEqn.H"
            #include "TEqn.H"
            #include "pEqn.H"
            Info<< "Done!" << endl;
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
            }
            //Info << "sens Sum: " << runTime.timeName() << "\t" << gSum(sens) << endl;
        }
    }

    void start(){
    #include "settings.H"
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{    
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    
    turbulence->validate();

    simpleCheckpointer simpleCheck
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
         porosity,
         MRF,
         radiation,
         rhoCpRef,
         fvOptions,
         cumulativeContErr
    );
    
    dimensionedScalar startTime = runTime.value();
    //label time = runTime.value();
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
