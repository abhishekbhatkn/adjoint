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

#include "GCMMASolver.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class MMAProgram {

private:
    Foam::Time& runTime;
    Foam::fvMesh& mesh;
    Foam::simpleControl& simple;
    Foam::volScalarField& p;
    Foam::volVectorField& U;
    Foam::volScalarField& T;
    volScalarField& eta;
    volScalarField& eta_old1;
    volScalarField& eta_old2;
    volScalarField& low;
    volScalarField& upp;
    volScalarField& sens;
    volScalarField& dfdeta;
    surfaceScalarField& phi;

    const Foam::wordList costFunctionPatches;

    // member variables  
    Foam::scalar volumeConstraint, designVolume, optEpsilon, factorT, factorP, minTempCost, maxTempCost, minPrCost, maxPrCost, asyminit, asymdec, asyminc;
    label maxPiggyLoop, designSize, maxMMAiter, MMALoop;
    List<label> designSpaceCells;
    List<label> solidSpaceCells;
    
    Foam::scalar fluidCp, fluidRho;
    


public:
    MMAProgram
    (
        Foam::Time& runTime,
        Foam::fvMesh& mesh,
        Foam::simpleControl& simple,
        Foam::volScalarField& p,
        Foam::volVectorField& U,
	Foam::volScalarField& T,
        volScalarField& eta,
        volScalarField& eta_old1,
        volScalarField& eta_old2,
        volScalarField& low,
        volScalarField& upp,
        volScalarField& sens,
        volScalarField& dfdeta,
        surfaceScalarField& phi
    ) :
        runTime(runTime),
        mesh(mesh),
        simple(simple),
        p(p),
        U(U),
	T(T),
        eta(eta),
        eta_old1(eta_old1),
        eta_old2(eta_old2),
        low(low),
        upp(upp),
        sens(sens),
        dfdeta(dfdeta),
        phi(phi),
        costFunctionPatches(mesh.solutionDict().subDict("SIMPLE").lookup("costFunctionPatches")) 	
{
	#include "settings.H"
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
	//value[1] = 100*((val/designVolume)-(volumeConstraint));
	Info << "fval = " << value << endl;
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
    	if (factorT!=0) {
		J = -factorT*(calcTempCost()-minTempCost)/(maxTempCost-minTempCost) + factorP*(calcPrCost()-minPrCost)/(maxPrCost-minPrCost);
		return J;
	}
	else {
		J = factorP*(calcPrCost()-minPrCost)/(maxPrCost-minPrCost);
		return J;
	}
    }
    
    bool MMASolve() {

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
	scalar J = calcCost();
	List<scalar> fval = calcFval();
	gcmma.MMAUpdate(eta, eta_old1, eta_old2, J, sens, fval, dfdeta);
	check = etaCheck();
	Info<< "Loop "<< MMALoop << " Time: " << runTime.timeName() << " cost: "<< J << " fval: " << fval  << " check: " << check << " asyminit: " << asyminit << " asymdec: " << asymdec << " asyminc: " << asyminc << endl;
	runTime.writeNow();
	Info<< "MMA Solver End\n" << endl;
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
    
//------------------------------------------------------------------------------------------------------//
    MMAProgram program
    (
         runTime,
         mesh,
         simple,
         p,
         U,
	  T,
         eta,
         eta_old1,
         eta_old2,
         low,
         upp,
         sens,
         dfdeta,
         phi
    );
     
    program.MMASolve();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

