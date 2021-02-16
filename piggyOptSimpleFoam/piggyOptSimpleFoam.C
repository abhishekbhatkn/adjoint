/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    #include "CHTSettings.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    CheckInterface check(runTime);
    CheckDict checkDict(runTime);
    CheckDatabase checkDB(runTime, checkDict);

    Info<< "\nStarting time loop\n" << endl;

    scalar oldSensSum = 0;
/*
    #include "Random.H"
    Foam::Random randObj(time(NULL));
    forAll(designSpaceCells,i){
    	const label j = designSpaceCells[i];
       	eta[j] = randObj.bit();
    }
*/

    label tapeSizeMB = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("tapeSizeMB",4096);
    Info << "Creating Tape, size: " << tapeSizeMB << endl;
    AD::create_global_tape(tapeSizeMB/Pstream::nProcs());


    bool firstStep = true;

    turbulence->validate();

    label NN = eta.size();

    Foam::reduce(NN,sumOp<label>());
    label nOptSteps = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("nOptSteps",100);
    scalar optEpsilon = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("optTolerance",5e-2);
    scalar optStepwidth = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("optStepwidth",0.1);

    for(int optStep = 0; optStep < nOptSteps; optStep++){
        ADmode::global_tape->reset();
        ADmode::global_tape->switch_to_active();

	//AD::register_input_variable(eta.begin(),eta.end());

        forAll(eta,i){
            ADmode::global_tape->register_variable(eta[i]);
        }

        ADmode::tape_t::position_t reset_to = ADmode::global_tape->get_position();
        scalar dSensSum = std::numeric_limits<double>::max();

	scalar sensSum = 0, J = 0;

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
                ADmode::global_tape->switch_to_passive();
            }

            turbulence->correct();

            if(frozenTurbulence){
                ADmode::global_tape->switch_to_active();
            }

            J = CostFunction(mesh).eval();
	    checkDB.registerAsOutput();
            if(Pstream::master()){
                dco::derivative(J) = 1.0;
            }

            if(!firstStep){
                checkDB.restoreAdjoints();
            }else{
                firstStep = false;
            }

            ADmode::global_tape->interpret_adjoint(); //_to(interpret_to);
            checkDB.storeAdjoints();
            scalar norm2 = checkDB.calcNormOfStoredAdjoints();

            forAll(designSpaceCells,i){
                const label j = designSpaceCells[i];
                sens[j] = dco::derivative(eta[j]) / mesh.V()[j];
                sensSum += std::abs(dco::passive_value(sens[j]));
            }

            Foam::reduce(sensSum,sumOp<scalar>());
            Foam::reduce(norm2,sumOp<scalar>());
            dSensSum = abs(sensSum - oldSensSum)/NN;
            Info << "piggy: " << optStep << " " << runTime.timeName() << " " << sensSum << " " << dSensSum << " " << norm2 << " " << J << endl;
            oldSensSum = sensSum;

            ADmode::global_tape->reset_to(reset_to);
            ADmode::global_tape->zero_adjoints();

            runTime.write();
	    //sens.write();
	    //eta.write();
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = "   << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
        ADmode::global_tape->switch_to_passive();
	
        forAll(designSpaceCells,i){
            const label j = designSpaceCells[i];
	    eta[j] -= sens[j]*optStepwidth;
            eta[j] =  max(eta[j],scalar(0.0));
            eta[j] =  min(eta[j],scalar(1.0));
        }
	alpha = alpha_s + (alpha_f - alpha_s) * eta * (1.0 + penalty) / (eta + penalty);
    }
    runTime.write();
    Info << "avg/max eta,: " << gAverage(eta) << " " << gMax(eta) << endl;

    ADmode::tape_t::remove(ADmode::global_tape);
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
