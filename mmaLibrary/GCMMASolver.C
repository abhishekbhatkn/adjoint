////////////////////////////////////////////////////////////////////////////////
// Copyright © 2018 Jérémie Dumas
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

//#include "postProcess.H"
#include "GCMMASolver.H"
//#include "fvCFD.H"
// #include "forces.H"
//#include "turbulentTransportModel.H"
#include "zeroGradientFvPatchFields.H"
//#include "fvcGrad.H"
//#include <algorithm>
//#include <cmath>
//#include <cstdio>
//#include <iostream>

////////////////////////////////////////////////////////////////////////////////
// PUBLIC
////////////////////////////////////////////////////////////////////////////////

GCMMASolver::GCMMASolver(
	  Foam::fvMesh &mesh
	, Foam::List<label> &designSpaceCells
	, Foam::label mm
	, Foam::volScalarField& low
	, Foam::volScalarField& upp
	, Foam::scalar asyminit
	, Foam::scalar asymdec
	, Foam::scalar asyminc
	, Foam::label outeriter
   ) :
	  mesh(mesh)
	, designSpaceCells(designSpaceCells)
	, m(mm)
	, outeriter(outeriter)
	, dummy(IOobject("dummy",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),	mesh, Foam::dimensionedScalar("0", dimless, 0.0), zeroGradientFvPatchScalarField::typeName)
	, raa0eps(1e-6)
	, raaeps(1e-6)
	, xmamieps(1e-5)
	, epsimin(1e-7)
	, move(0.5)
	, albefa(0.1)
	, asyminit(asyminit)
	, asymdec(asymdec)
	, asyminc(asyminc)
	, raa0(1.0)
	, z(0.0)
	, r0(0.0)
	, f0app(0.0)
	, raa(m,1.0)
	, a(m,0.0)
	, c(m,1000.0)
	, d(m,0.0)
	, y(m,0.0)
	, lam(m,0.0)
	, mu(m,0.0)
	, b(m,0.0)
	, grad(m,0.0)
	, r(m,0.0)
	, fapp(m,0.0)
	, s(2*m,0.0)
	, hess(m*m,0.0)
	
	, low(low)
	, upp(upp)
	, alpha("alpha",dummy)
	, beta("beta",dummy)
	, p0("p0",dummy)
	, q0("q0",dummy)
	, pi_0("pi_0",dummy)
	, pi_1("pi_1",dummy)
	, qi_0("qi_0",dummy)
	, qi_1("qi_1",dummy)
{ }

/*
void GCMMASolver::OuterUpdate(Foam::volScalarField& xmma, const Foam::volScalarField& xval, Foam::scalar& f0x, const Foam::volScalarField& df0dx, const Foam::List<Foam::scalar>& fx, const Foam::volScalarField& dfdx)
{
	// Compute asymptotes
	Asymp(xval, df0dx, dfdx); //, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx); //, xmin, xmax);

	// Update xolds
	xold2 = xold1;
	xold1 = xval;

	// Solve the dual with an Foam::labelerior poFoam::label method
	xmma = xval;
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	//SolveDSA(xmma);

	// Compute approximation values
	ComputeApprox(xmma);
	
	//Info<< "OuterUpdate Complete! " << endl;

}
*/
void GCMMASolver::MMAUpdate(Foam::volScalarField& xval, Foam::volScalarField& xval_old1, Foam::volScalarField& xval_old2, Foam::scalar& f0x, const Foam::volScalarField& df0dx, const Foam::List<Foam::scalar>& fx, const Foam::volScalarField& dfdx)
{

	// Compute asymptotes
	Asymp(xval, xval_old1, xval_old2, df0dx, dfdx); //, xmin, xmax);
	
	raa0 = 0.001;
	for (Foam::label j = 0; j < m; ++j) {
		raa[j] = 0.001;
	}

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx); //, xmin, xmax);

	// Update xolds
	xval_old2 = xval_old1;
	xval_old1 = xval;

	// Solve the dual with an Foam::labelerior poFoam::label method
	//SolveDIP(xval);

	// Solve the dual with a steepest ascent method
	SolveDSA(xval);
	
	//Foam::Info<< "MMAUpdate Complete! " << endl;

}
/*
void GCMMASolver::InnerUpdate(Foam::volScalarField &xmma, Foam::scalar& f0xnew, const Foam::List<Foam::scalar>& fxnew, const Foam::volScalarField& xval, Foam::scalar& f0x, const Foam::volScalarField& df0dx, const Foam::List<Foam::scalar>& fx, const Foam::volScalarField& dfdx)
{
	// Update approximation factors
	RaaUpdate(xmma, xval, f0xnew, fxnew); //, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx); //, xmin, xmax);

	// Solve the dual with an Foam::labelerior poFoam::label method
	xmma = xval;
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	//SolveDSA(xmma);

	// Compute approximation values
	ComputeApprox(xmma);
	
	//Info<< "InnerUpdate Complete! " << endl;
}

bool GCMMASolver::ConCheck(Foam::scalar& f0xnew, Foam::List<Foam::scalar>& fxnew) const {
	//Info<< "ConCheck Start! " << endl;
	if (f0app + epsimin < f0xnew) {
		return false;
	}
	for (Foam::label j = 0; j < m; ++j) {
		if (fapp[j] + epsimin < fxnew[j]) {
			return false;
		}
	}
	return true;
}
*/
////////////////////////////////////////////////////////////////////////////////
// PRIVATE
////////////////////////////////////////////////////////////////////////////////

// Compute [low, upp, raa0, raa]
void GCMMASolver::Asymp(const Foam::volScalarField& xval, const Foam::volScalarField& xval_old1, const Foam::volScalarField& xval_old2, const Foam::volScalarField& df0dx, const Foam::volScalarField& dfdx) //, const double *xmin, const double *xmax)
{
	// Forward the iterator
	++outeriter;

	// Set asymptotes
	if (outeriter < 3) {
		forAll(designSpaceCells,jj) {
			const Foam::label i = designSpaceCells[jj];
			low[i] = xval[i] - asyminit;
			upp[i] = xval[i] + asyminit;
		}
	} else {
		forAll(designSpaceCells,jj) {
			const Foam::label i = designSpaceCells[jj];
			Foam::scalar zzz = (xval[i] - xval_old1[i]) * (xval_old1[i] - xval_old2[i]);
			Foam::scalar gamma = 1.0;
			if (zzz < 0.0) {
				gamma = asymdec;
			} else if (zzz > 0.0) {
				gamma = asyminc;
			} else {
				gamma = 1.0;
			}
			low[i] = xval[i] - gamma * (xval_old1[i] - low[i]);
			upp[i] = xval[i] + gamma * (upp[i] - xval_old1[i]);

			Foam::scalar xmami = Foam::max(xmamieps, Foam::scalar(1.0));
			// double xmami = xmax[i] - xmin[i];
			low[i] = Foam::max(low[i], xval[i] - 10.0 * xmami);
			low[i] = Foam::min(low[i], xval[i] - 0.001 * xmami);
			upp[i] = Foam::max(upp[i], xval[i] + 0.001 * xmami);
			upp[i] = Foam::min(upp[i], xval[i] + 10.0 * xmami);

/*			Foam::scalar xmi = 0.0 - 1.0e-6;
			Foam::scalar xma = 1.0 + 1.0e-6;
			if (xval[i] < xmi) {
				low[i] = xval[i] - (xma - xval[i]) / Foam::scalar(0.9);
				upp[i] = xval[i] + (xma - xval[i]) / Foam::scalar(0.9);
			}
			if (xval[i] > xma) {
				low[i] = xval[i] - (xval[i] - xmi) / Foam::scalar(0.9);
				upp[i] = xval[i] + (xval[i] - xmi) / Foam::scalar(0.9);
			}
*/		}
	}
	// Set raa0, raa
/*	raa0 = 0;
	for (Foam::label j = 0; j < m; ++j) {
		raa[j] = 0;
	}
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		Foam::scalar xmami = Foam::max(xmamieps, Foam::scalar(1.0));
		raa0 += Foam::mag(df0dx[i]) * xmami;
		for (Foam::label j = 0; j < m; ++j) {
			raa[j] += Foam::mag(dfdx[i]) * xmami;
		}
	}
	Foam::reduce(raa0,sumOp<scalar>());
	for (Foam::label j = 0; j < m; ++j) {
		Foam::reduce(raa[j],sumOp<scalar>());
	}	
	Foam::label n = designSpaceCells.size();
	Foam::reduce(n,sumOp<label>());
	
	raa0 = Foam::max(raa0eps, (0.1/n) * raa0); //Foam::max(raa0eps, (0.1/n) * raa0);
	for (Foam::label j = 0; j < m; ++j) {
		raa[j] = Foam::max(raaeps, (0.1/n) * raa[j]); //Foam::max(raaeps, (0.1/n) * raa[j]);
	}
	//Info << "raa0 value: " << raa0 <<" raa value: " << sum(raa) <<"\n" << endl;
	*/
}

void GCMMASolver::SolveDIP(Foam::volScalarField& x) {
	
	for (label j = 0; j < m; ++j) {
		lam[j] = c[j] / 2.0;
		mu[j] = 1.0;
	}

	const Foam::scalar tol = epsimin; // 1.0e-9*Foam::sqrt(m+n);
	Foam::scalar epsi = 1.0;
	Foam::scalar err = 1.0;
	Foam::label loop;

	while (epsi > tol) {

		loop = 0;
		while (err > 0.9 * epsi && loop < 1000) {
			loop++;

			// Set up Newton system
			XYZofLAMBDA(x);
			DualGrad(x);
			for (Foam::label j = 0; j < m; ++j) {
				grad[j] = -1.0 * grad[j] - epsi / lam[j];
			}
			DualHess(x);

			// Solve Newton system
			if (m > 1) {
				Factorize(hess);
				Solve(hess, grad);
				for (Foam::label j = 0; j < m; ++j) {
					s[j] = grad[j];
				}
			} else if (m > 0) {
				s[0] = grad[0] / hess[0];
			  }

			// Get the full search direction
			for (Foam::label i = 0; i < m; ++i) {
				s[m+i] = -mu[i] + epsi / lam[i] - s[i] * mu[i] / lam[i];
			}

			// Perform linesearch and update lam and mu
			DualLineSearch();

			XYZofLAMBDA(x);

			// Compute KKT res
			err = DualResidual(x, epsi);
			//Foam::Info << "epsi: " << epsi << " err: " << err << endl;
		}
		epsi = epsi * 0.1;
	}
	//Info<< "SolveDIP Complete! " << endl;
}


void GCMMASolver::XYZofLAMBDA(Foam::volScalarField& x) {

	Foam::scalar lamai = 0.0;
	for (Foam::label j = 0; j < m; ++j) {
		if (lam[j] < 0.0) {
			lam[j] = 0.0;
		}
		y[j] = Foam::max(Foam::scalar(0.0), lam[j] - c[j]); // Note y=(lam-c)/d - however d is fixed at one !!
		lamai += lam[j]* a[j];
	}

	z = Foam::max(Foam::scalar(0.0), 10.0 * (lamai - 1.0)); // SINCE a0 = 1.0

	Foam::scalar pjlam = 0.0, qjlam = 0.0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		pjlam = p0[i];
		qjlam = q0[i];
		for (Foam::label j = 0; j < m; ++j) {
			if (j==0) {
				pjlam += pi_0[i] * lam[j];
				qjlam += qi_0[i] * lam[j];
			}
			if (j==1) {
				pjlam += pi_1[i] * lam[j];
				qjlam += qi_1[i] * lam[j];
			}
		}
		x[i] = (Foam::sqrt(pjlam) * low[i] + Foam::sqrt(qjlam) * upp[i]) / (Foam::sqrt(pjlam) + Foam::sqrt(qjlam));
		if (x[i] < alpha[i]) {
			x[i] = alpha[i];
		}
		if (x[i] > beta[i]) {
			x[i] = beta[i];
		}
	}
	//Info<< "XYZofLAMBDA Complete! " << endl;
}

void GCMMASolver::DualGrad(Foam::volScalarField& x) {
	for (Foam::label j = 0; j < m; ++j) {
		grad[j] = -b[j] - a[j] * z - y[j];
		forAll(designSpaceCells,jj) {
			const Foam::label i = designSpaceCells[jj];
			if (j==0) {
				grad[j] += pi_0[i] / (upp[i] - x[i]) + qi_0[i] / (x[i] - low[i]);
			}
			if (j==1) {
				grad[j] += pi_1[i] / (upp[i] - x[i]) + qi_1[i] / (x[i] - low[i]);
			}
		}
		Foam::reduce(grad[j],sumOp<scalar>());
	}
	//Info<< "DualGrad Complete! " << endl;
}

void GCMMASolver::DualHess(Foam::volScalarField& x) {
	Foam::volScalarField df2("df2", dummy);
	
	Foam::volScalarField PQ1("PQ1", dummy);
	Foam::volScalarField PQ2("PQ2", dummy);

	Foam::scalar pjlam = 0, qjlam = 0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		pjlam = p0[i];
		qjlam = q0[i];
		for (Foam::label j = 0; j < m; ++j) {
			if (j==0) {
				pjlam += pi_0[i] * lam[j];
				qjlam += qi_0[i] * lam[j];
				PQ1[i] = pi_0[i] / Foam::pow(upp[i] - x[i], 2.0) - qi_0[i] / Foam::pow(x[i] - low[i], 2.0);
			}
			if (j==1) {
				pjlam += pi_1[i] * lam[j];
				qjlam += qi_1[i] * lam[j];
				PQ2[i] = pi_1[i] / Foam::pow(upp[i] - x[i], 2.0) - qi_1[i] / Foam::pow(x[i] - low[i], 2.0);
			}
		}
		df2[i] = -1.0 / (2.0 * pjlam / Foam::pow(upp[i] - x[i], 3.0) + 2.0 * qjlam / Foam::pow(x[i] - low[i], 3.0));
		Foam::scalar xp = (Foam::sqrt(pjlam) * low[i] + Foam::sqrt(qjlam) * upp[i]) / (Foam::sqrt(pjlam) + Foam::sqrt(qjlam));
		if (xp < alpha[i]) {
			df2[i] = 0.0;
		}
		if (xp > beta[i]) {
			df2[i] = 0.0;
		}
	}

	
	//Foam::List<Foam::volScalarField> tmp(m, dummy);
	
/*	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		for (Foam::label = 0; j < m; ++j) {
			tmp[j][i] += PQ[j][i] * df2[i];
		}
	}
	for (Foam::label j = 0; j < m; ++j) {
		for (Foam::label k = 0; k < m; ++k) {
			forAll(designSpaceCells,jj) {
				const Foam::label i = designSpaceCells[jj];
				hess[j*m+k] += PQ[j][i] * df2[i] * PQ[k][i];
			}
			Foam::reduce(hess[j*m+k],sumOp<scalar>());
		}
	}
	// Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ*/
	for (Foam::label j = 0; j < m*m; ++j) {
		hess[j] = 0;
	}
	//Info<< "hess: " << hess << endl;
	for (Foam::label j = 0; j < m*m; ++j) {
		forAll(designSpaceCells,jj) {
			const Foam::label i = designSpaceCells[jj];
			if (j==0) {
				hess[j] += PQ1[i] * df2[i] * PQ1[i];
			}
			if (j==1||j==2) {
				hess[j] += PQ2[i] * df2[i] * PQ1[i];
			}
			if (j==3) {
				hess[j] += PQ2[i] * df2[i] * PQ2[i];
			}
		}
	}
	for (Foam::label j = 0; j < m*m; ++j) {
		Foam::reduce(hess[j],sumOp<scalar>());
	}
	Foam::scalar lamai = 0.0;
	for (Foam::label j = 0; j < m; ++j) {
		if (lam[j] < 0.0) {
			lam[j] = 0.0;
		}
		lamai += lam[j] * a[j];
		if (lam[j] > c[j]) {
			hess[j] += -1.0;
		}
		hess[j*m+j] += -mu[j] / lam[j];
	}

	if (lamai > 0.0) {
		for (Foam::label j = 0; j < m; ++j) {
			for (Foam::label k = 0; k < m; ++k) {
				hess[j*m+k] += -10.0 * a[j] * a[k];
			}
		}
	}

	//Info<< "hess: " << hess << endl;
	// pos def check
	Foam::scalar hessTrace = 0.0;
	for (Foam::label j = 0; j < m; ++j) {
		hessTrace += hess[j*m+j];
	}

	Foam::scalar hessCorr = 1e-4 * hessTrace / m;

	if (-1.0 * hessCorr < 1.0e-7) {
		hessCorr = -1.0e-7;
	}

	for (Foam::label j = 0; j < m; ++j) {
		hess[j*m+j] += hessCorr;
	}
	//Info<< "DualHess Complete! " << endl;
}

void GCMMASolver::SolveDSA(Foam::volScalarField& x) {

	lam = 1.0;

	const Foam::scalar tol = epsimin; // 1.0e-9*Foam::sqrt(m+n);
	Foam::scalar err = 1.0;
	Foam::label loop = 0;

	while (err > tol && loop < 500) {
		loop++;
		XYZofLAMBDA(x);
		DualGrad(x);
		Foam::scalar theta = 1.0;
		err = 0.0;
		for (Foam::label j = 0; j < m; ++j) {
			lam[j] = Foam::max(Foam::scalar(0.0), lam[j] + theta * grad[j]);
			err += grad[j] * grad[j];
		}
		err = Foam::sqrt(err);
	}
	//Info<< "SolveDSA Complete! " << endl;
}

Foam::scalar GCMMASolver::DualResidual(Foam::volScalarField& x, Foam::scalar& epsi) {

	Foam::List<Foam::scalar> res(2*m);
	for (Foam::label j = 0; j < m; ++j) {
		res[j] = -b[j] - a[j] * z - y[j] + mu[j];
		res[j + m] = mu[j] * lam[j] - epsi;
		forAll(designSpaceCells,jj) {
			const Foam::label i = designSpaceCells[jj];
			if (j==0) {
				res[j] += pi_0[i] / (upp[i] - x[i]) + qi_0[i] / (x[i] - low[i]);
			}
			if (j==1) {
				res[j] += pi_1[i] / (upp[i] - x[i]) + qi_1[i] / (x[i] - low[i]);
			}
		}
		Foam::reduce(res[j],sumOp<scalar>());
	}

	Foam::scalar nrI = 0.0;
	for (Foam::label i = 0; i < 2 * m; ++i) {
		if (nrI < Foam::mag(res[i])) {
			nrI = Foam::mag(res[i]);
		}
	}
	//Info<< "DualResidual Complete! " << endl;
	return nrI;
}

void GCMMASolver::DualLineSearch() {

	Foam::scalar theta = 1.005;
	for (Foam::label i = 0; i < m; ++i) {
		if (theta < -1.01 * s[i] / lam[i]) {
			theta = -1.01 * s[i] / lam[i];
		}
		if (theta < -1.01 * s[i + m] / mu[i]) {
			theta = -1.01 * s[i + m] / mu[i];
		}
	}
	theta = 1.0 / theta;

	for (Foam::label i = 0; i < m; ++i) {
		lam[i] = lam[i] + theta * s[i];
		mu[i] = mu[i] + theta * s[i + m];
	}
	//Info<< "DualLineSearch Complete! " << endl;
}


// Update [raa0, raa]
void GCMMASolver::RaaUpdate(const Foam::volScalarField& xmma, const Foam::volScalarField& xval, Foam::scalar& f0xnew, const Foam::List<Foam::scalar>& fxnew) //  const double *xmin, const double *xmax)
{
	const Foam::scalar raacofmin = 1e-12;
	Foam::scalar raacof = 0.0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		Foam::scalar xmami = Foam::max(xmamieps, Foam::scalar(1.0));
		Foam::scalar xxux = (xmma[i] - xval[i]) / (upp[i] - xmma[i]);
		Foam::scalar xxxl = (xmma[i] - xval[i]) / (xmma[i] - low[i]);
		Foam::scalar xxul = xxux * xxxl;
		Foam::scalar ulxx = (upp[i] - low[i]) / xmami;
		raacof += xxul * ulxx;
	}
	Foam::reduce(raacof,sumOp<scalar>());
	raacof = Foam::max(raacofmin, raacof);

	if (f0xnew > f0app + 0.5 * epsimin) {
		Foam::scalar deltaraa0 = (1.0 / raacof) * (f0xnew - f0app);
		raa0 = Foam::min(1.1 * (raa0 + deltaraa0), 10.0 * raa0);
	}

	for (Foam::label j = 0; j < m; ++j) {
		if (fxnew[j] > fapp[j] + 0.5 * epsimin) {
			Foam::scalar deltaraa = (1.0 / raacof) * (fxnew[j] - fapp[j]);
			raa[j] = Foam::min(1.1 * (raa[j] + deltaraa), 10.0 * raa[j]);
		}
		//Info<< "raa[" << j << "]: " << raa[j] << endl;
	}
	//Info<< "RaaUpdate Complete! " << endl;
}

void GCMMASolver::GenSub(const Foam::volScalarField& xval, Foam::scalar& f0x, const Foam::volScalarField& df0dx, const Foam::List<Foam::scalar>& fx, const Foam::volScalarField& dfdx)
	//, const double *xmin, const double *xmax)
{
	// Set bounds and the coefficients for the approximation
	r0 = 0.0;
	for(Foam::label j = 0; j < m; ++j) {
		r[j] = 0;
	}
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		// Compute bounds alpha and beta
		alpha[i] = Foam::max(Foam::scalar(0.0), low[i] + albefa * (xval[i] - low[i]));
		alpha[i] = Foam::max(alpha[i], xval[i] - move);
		// alpha[i] = Foam::min(alpha[i], xmax[i]);
		beta[i] = Foam::min(1.0, upp[i] - albefa * (upp[i] - xval[i]));
		beta[i] = Foam::min(beta[i], xval[i] + move);
		// beta[i]  = Foam::max(beta[i], xmin[i]);

		Foam::scalar uxinv = 1.0 / (upp[i] - xval[i]);
		Foam::scalar xlinv = 1.0 / (xval[i] - low[i]);

		// Objective function
		{
			Foam::scalar df0dxp = Foam::max(Foam::scalar(0.0), df0dx[i]);
			Foam::scalar df0dxm = Foam::max(Foam::scalar(0.0), -1.0 * df0dx[i]);
			Foam::scalar xmamiinv = 1.0 / Foam::max(xmamieps, scalar(1.0));
			Foam::scalar pq = 0.001 * Foam::mag(df0dx[i]) + raa0 * xmamiinv;
			p0[i] = Foam::pow(upp[i] - xval[i], 2.0) * (df0dxp + pq);
			q0[i] = Foam::pow(xval[i] - low[i], 2.0) * (df0dxm + pq);
			r0 += p0[i] * uxinv + q0[i] * xlinv;
		}
		// Constraints
		{
			for (Foam::label j = 0; j < m; ++j) {
				Foam::scalar dfdxp = Foam::max(Foam::scalar(0.0), dfdx[i]);
				Foam::scalar dfdxm = Foam::max(Foam::scalar(0.0), -1.0 * dfdx[i]);
				Foam::scalar xmamiinv = 1.0 / Foam::max(xmamieps, 1.0);
				Foam::scalar pq = 0.001 * Foam::mag(dfdx[i]) + raa[j] * xmamiinv;
				if (j==0) {
					pi_0[i] = Foam::pow(upp[i] - xval[i], 2.0) * (dfdxm + pq);
					qi_0[i] = Foam::pow(xval[i] - low[i], 2.0) * (dfdxp + pq);
					r[j] += pi_0[i] * uxinv + qi_0[i] * xlinv;
				}
				if (j==1) {
					pi_1[i] = Foam::pow(upp[i] - xval[i], 2.0) * (dfdxm + pq);
					qi_1[i] = Foam::pow(xval[i] - low[i], 2.0) * (dfdxp + pq);
					r[j] += pi_1[i] * uxinv + qi_1[i] * xlinv;
				}
			}
		}
	}
	Foam::reduce(r0,sumOp<scalar>());
	for (Foam::label j = 0; j < m; ++j) {
		Foam::reduce(r[j],sumOp<scalar>());
	}
	r0 = f0x - r0;
	for (Foam::label j = 0; j < m; ++j) {
		r[j] = fx[j] - r[j];
		b[j] = -r[j]; // The constant for the constraints
	}
	//Info<< "GenSub Complete! " << endl;
}

void GCMMASolver::ComputeApprox(const Foam::volScalarField& xmma) {
	f0app = 0.0;
	for(Foam::label j = 0; j < m; ++j) {
		fapp[j] = 0;
	}
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		Foam::scalar uxinv = 1.0 / (upp[i] - xmma[i]);
		Foam::scalar xlinv = 1.0 / (xmma[i] - low[i]);
		f0app += p0[i] * uxinv + q0[i] * xlinv;
		for (Foam::label j = 0; j < m; ++j) {
			if (j==0) {
				fapp[j] += pi_0[i] * uxinv + qi_0[i] * xlinv;
			}
			if (j==1) {
				fapp[j] += pi_1[i] * uxinv + qi_1[i] * xlinv;
			}
		}
	}
	Foam::reduce(f0app,sumOp<scalar>());
	f0app += r0;
	for (Foam::label j = 0; j < m; ++j) {
		Foam::reduce(fapp[j],sumOp<scalar>());
		fapp[j] += r[j];
	}
	//Info<< "ComputeApprox Complete! " << endl;	
}


void GCMMASolver::Factorize(Foam::List<scalar>& K) {

	for (Foam::label s = 0; s < m - 1; ++s) {
		for (Foam::label i = s + 1; i < m; ++i) {
			K[i * m + s] = K[i * m + s] / K[s * m + s];
			for (Foam::label j = s + 1; j < m; ++j) {
				K[i * m + j] = K[i * m + j] - K[i * m + s] * K[s * m + j];
			}
		}
	}
	//Info<< "Factorize Complete! " << endl;
}

void GCMMASolver::Solve(Foam::List<scalar>& K, Foam::List<scalar>& x) {

	for (Foam::label i = 1; i < m; ++i) {
		Foam::scalar aa = 0.0;
		for (Foam::label j = 0; j < i; ++j) {
			aa = aa - K[i * m + j] * x[j];
		}
		x[i] = x[i] + aa;
	}

	x[m - 1] = x[m - 1] / K[(m - 1) * m + (m - 1)];
	for (Foam::label i = m - 2; i >= 0; --i) {
		Foam::scalar aa = x[i];
		for (Foam::label j = i + 1; j < m; ++j) {
			aa = aa - K[i * m + j] * x[j];
		}
		x[i] = aa / K[i * m + i];
	}
	//Info<< "Solve Complete! " << endl;
}
