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
#include "costFunctionLibrary.C"

#include "simpleMMACheckpointer.H"

simpleMMACheckpointer::simpleMMACheckpointer(
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
        volScalarField& eta,
        volScalarField& sens,
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
        eta(eta),
        sens(sens),
        dfdeta(dfdeta),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        alpha_s("alpha_s", dimless/dimTime, 0.0),
        alpha_f("alpha_f", dimless/dimTime, 0.0)
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

    bool simpleMMACheckpointer::runStep(){
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

    scalar simpleMMACheckpointer::calcCost(){
	Info << "ACTION::calcCost J" << endl;
	return CostFunction(mesh).eval();
    }

    void simpleMMACheckpointer::postInterpret(){
        //Info << "Writing sensitivity for timestep " << runTime.timeName() << endl;
        for(int i = 0; i < eta.size(); i++)
    		sens[i] = AD::derivative(eta[i]); // mesh.V()[i];

        Info << "sens Sum: " << runTime.timeName() << "\t" << gSum(sens) << endl;
        //runTime.stopAt(Foam::Time::stopAtControls::saEndTime);
        if(runTime.writeTime()){
	    Info << "Post Interpret: Writing sensitivity for timestep " << runTime.timeName() << endl;
	    sens.write();
        }
        //runTime.write();
    }
    void simpleMMACheckpointer::write(bool firstRun){
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

    void simpleMMACheckpointer::start(){
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

    penalty = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("penalty",0.6);
    if (!mesh.solutionDict().subDict("SIMPLE").found("penalty"))
	    {
		Info << "Warning: Keyword penalty not found in fvSolution/SIMPLE. Default set to 0.6" << endl;
	    }
    alpha_f =dimensionedScalar("alpha_f", dimless/dimTime, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_f",0.1));
    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_f"))
	    {
		Info << "Warning: Keyword alpha_f not found in fvSolution/SIMPLE. Default set to 0.1" << endl;
	    }
    alpha_s = dimensionedScalar("alpha_s", dimless/dimTime, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_s",10e5));
    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_s"))
	    {
		Info << "Warning: Keyword alpha_s not found in fvSolution/SIMPLE. Default set to 10e5" << endl;
	    }
    }
    
    List<label> simpleMMACheckpointer::getDesignSpace(){
    	return designSpaceCells;
    }
    
#include "Step_1.H"
