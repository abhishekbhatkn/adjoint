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
    
    // member variables

    scalar penalty, volumeConstraint;
    label n, nOptSteps;
    //bool frozenTurbulence;
    List<label> designSpaceCells;
    dimensionedScalar alpha_s, alpha_f;

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
        volScalarField& alpha,
        volScalarField& eta,
	volScalarField& eta_MMA,
        volScalarField& sens,
        volScalarField& dfdeta,
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
        alpha(alpha),
        eta(eta),
        eta_MMA(eta_MMA),
        sens(sens),
        dfdeta(dfdeta),
        phi(phi),
        pRefCell(pRefCell),
        pRefValue(pRefValue),
        cumulativeContErr(cumulativeContErr),
        alpha_s(" ",dimless/dimTime, 0.0),
        alpha_f(" ",dimless/dimTime, 0.0)
{ }

    void runLoop() {
	Info<< "\nStarting time loop\n" << endl;
	label iter = 0;
	while (iter < nOptSteps && simple.loop())
	{
	    	iter++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity SIMPLE corrector
		{
		    alpha=alpha_s+(alpha_f-alpha_s)*eta*(1.0+penalty)/(eta + penalty);
		    #include "UEqn.H"
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
    
    scalar calcFval() {  
	scalar val;
	forAll(designSpaceCells,i){
	    const label j = designSpaceCells[i];
	    val += eta[j];
	}
	Foam::reduce(val,sumOp<scalar>());
	val = (val/n) - 1.0 + volumeConstraint;
	    
	Info << "fval = " << val << endl;
	    
	return val;
    }
    
    scalar calcCost(){
        Info << "ACTION::calcCost" << endl;
        return CostFunction(mesh).eval();
    }
    
    List<label> designSp(){
    	return designSpaceCells;
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
 
//-----------------------------------------------------------------------//   
    scalar J = 0, oldJ = 0.0, fval = 0.0, oldfval = 0.0;
      
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
         alpha,
         eta,
         eta_MMA,
         sens,
         dfdeta,
         phi,
         pRefCell,
         pRefValue,
         cumulativeContErr
    );
    
    program.start();
    List<label> designSpaceCells = program.designSp();

    GCMMASolver gcmma(mesh,designSpaceCells);
    
    volScalarField eta_old
    (
	    IOobject
	    (
		"eta_old",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("1.0", dimless, 1.0),
	    zeroGradientFvPatchScalarField::typeName
    );
    
    fval = program.calcFval();
    J = program.calcCost();
    
    gcmma.OuterUpdate(eta_MMA, eta, J, sens, fval, dfdeta);  
    eta_MMA.write();
    
    Info<< "GCMMA Outer Loop completed\n" << endl;

    eta_old = eta;
    eta_old.write();
       
    oldJ = J;
    oldfval = fval;

    eta = eta_MMA;
    program.runLoop();
    J = program.calcCost();
    fval = program.calcFval();
    bool conserv = gcmma.ConCheck(J,fval);
    Info<< "Updated Conservative: " << conserv << "\n" << endl;
    if (conserv) {
    	eta = eta_old;
    	J = oldJ;
    	fval = oldfval;
    }
    
    for (int inneriter = 0; !conserv && inneriter < 15; ++inneriter) {
	// Inner iteration update
	gcmma.InnerUpdate(eta_MMA, J, fval, eta, oldJ, sens, oldfval, dfdeta);
	eta_MMA.write();
	
	eta_old = eta;
	eta_old.write();
    
	// Check conservativity
	eta = eta_MMA;
	program.runLoop();
	J = program.calcCost();
	fval = program.calcFval();
	conserv = gcmma.ConCheck(J,fval);
	if (conserv) {
    		eta = eta_old;
    		J = oldJ;
    		fval = oldfval;
    	}
    }
    runTime.write();
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
