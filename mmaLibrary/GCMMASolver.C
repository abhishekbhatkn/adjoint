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

GCMMASolver::GCMMASolver(Foam::fvMesh &mesh, Foam::List<label> &designSpaceCells)
	: mesh(mesh)
	, designSpaceCells(designSpaceCells)
	//, n(designSpaceCells.size())

	, dummy(IOobject("dummy",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),	mesh, Foam::dimensionedScalar("0", dimless, 0.0), zeroGradientFvPatchScalarField::typeName)
	, outeriter(0)
	, raa0eps(1e-6)
	, raaeps(raa0eps)
	, xmamieps(1e-5)
	, epsimin(1e-7)
	, move(0.5)
	, albefa(0.1)
	, asyminit(0.5)
	, asymdec(0.7)
	, asyminc(1.2)
	, raa0(1.0)
	, raa(1.0)
	, a(0.0)
	, c(1000.0)
	, d(0.0)
	, y(0.0)
	, z(0.0)
	, s({0.0,0.0})
	, lam(0.0)
	, mu(0.0)
	, b(0.0)
	, grad(0.0)
	, hess(0.0)
	, r(0.0)
	, fapp(0.0)
	, r0(0.0)
	, f0app(0.0)
	
	, xold1("xold1",dummy)
	, xold2("xold2",dummy)
	, low("low",dummy)
	, upp("upp", dummy)
	, alpha("alpha",dummy)
	, beta("beta",dummy)
	, p0("p0",dummy)
	, q0("q0",dummy)
	, pij("pij",dummy)
	, qij("qij",dummy)
{ }


void GCMMASolver::SetAsymptotes(Foam::scalar init, Foam::scalar decrease, Foam::scalar increase) {
	// Asymptotes initialization and increase/decrease
	asyminit = init;
	asymdec = decrease;
	asyminc = increase;
}

void GCMMASolver::OuterUpdate(Foam::volScalarField& xmma, const Foam::volScalarField& xval, Foam::scalar& f0x,
				const Foam::volScalarField& df0dx, const Foam::scalar& fx, const Foam::volScalarField& dfdx)
{
	// Compute asymptotes
	Asymp(xval, df0dx, dfdx); //, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx); //, xmin, xmax);

	// Update xolds
	xold2 = xold1;
	xold1 = xval;

	// Solve the dual with an interior point method
	xmma = xval;
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	//SolveDSA(xmma);

	// Compute approximation values
	ComputeApprox(xmma);

}

void GCMMASolver::MMAUpdate(Foam::volScalarField& xmma, const Foam::volScalarField& xval, Foam::scalar& f0x,
				const Foam::volScalarField& df0dx, const Foam::scalar& fx, const Foam::volScalarField& dfdx)
{

	// Compute asymptotes
	Asymp(xval, df0dx, dfdx); //, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx); //, xmin, xmax);

	// Update xolds
	xold2 = xold1;
	xold1 = xval;

	// Solve the dual with an interior point method
	xmma = xval;
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	//SolveDSA(xmma);

}

void GCMMASolver::InnerUpdate(Foam::volScalarField &xmma, Foam::scalar& f0xnew, const Foam::scalar& fxnew, const Foam::volScalarField& xval, Foam::scalar& f0x,
                              const Foam::volScalarField& df0dx, const Foam::scalar& fx, const Foam::volScalarField& dfdx)
{
	// Update approximation factors
	RaaUpdate(xmma, xval, f0xnew, fxnew); //, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx); //, xmin, xmax);

	// Solve the dual with an interior point method
	xmma = xval;
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	//SolveDSA(xmma);

	// Compute approximation values
	ComputeApprox(xmma);
}

bool GCMMASolver::ConCheck(Foam::scalar& f0xnew, Foam::scalar& fxnew) const {
	if (f0app + epsimin < f0xnew) {
		return false;
	}
	if (fapp + epsimin < fxnew) {
		return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// PRIVATE
////////////////////////////////////////////////////////////////////////////////

// Compute [low, upp, raa0, raa]
void GCMMASolver::Asymp(const Foam::volScalarField& xval, const Foam::volScalarField& df0dx, const Foam::volScalarField& dfdx) //, const double *xmin, const double *xmax)
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
			Foam::scalar zzz = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]);
			Foam::scalar gamma = 1.0;
			if (zzz < 0.0) {
				gamma = asymdec;
			} else if (zzz > 0.0) {
				gamma = asyminc;
			} else {
				gamma = 1.0;
			}
			low[i] = xval[i] - gamma * (xold1[i] - low[i]);
			upp[i] = xval[i] + gamma * (upp[i] - xold1[i]);

			Foam::scalar xmami = Foam::max(xmamieps, Foam::scalar(1.0));
			// double xmami = xmax[i] - xmin[i];
			low[i] = Foam::max(low[i], xval[i] - 100.0 * xmami);
			low[i] = Foam::min(low[i], xval[i] - 1e-5 * xmami);
			upp[i] = Foam::max(upp[i], xval[i] + 1e-5 * xmami);
			upp[i] = Foam::min(upp[i], xval[i] + 100.0 * xmami);

			Foam::scalar xmi = 0.0 - 1.0e-6;
			Foam::scalar xma = 1.0 + 1.0e-6;
			if (xval[i] < xmi) {
				low[i] = xval[i] - (xma - xval[i]) / Foam::scalar(0.9);
				upp[i] = xval[i] + (xma - xval[i]) / Foam::scalar(0.9);
			}
			if (xval[i] > xma) {
				low[i] = xval[i] - (xval[i] - xmi) / Foam::scalar(0.9);
				upp[i] = xval[i] + (xval[i] - xmi) / Foam::scalar(0.9);
			}
		}
	}
	// Set raa0, raa
	raa0 = 0;
	raa = 0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		Foam::scalar xmami = Foam::max(xmamieps, Foam::scalar(1.0));
		raa0 += Foam::mag(df0dx[i]) * xmami;
		raa += Foam::mag(dfdx[i]) * xmami;
	}
	Foam::reduce(raa0,sumOp<scalar>());
	Foam::reduce(raa,sumOp<scalar>());
	Foam::label n = designSpaceCells.size();
	Foam::reduce(n,sumOp<label>());
	raa0 = Foam::max(raa0eps, (0.1/n) * raa0);
	raa = Foam::max(raaeps, (0.1/n) * raa);
}

void GCMMASolver::SolveDIP(Foam::volScalarField& x) {

	lam = c / 2.0;
	mu = 1.0;

	const Foam::scalar tol = epsimin; // 1.0e-9*Foam::sqrt(m+n);
	Foam::scalar epsi = 1.0;
	Foam::scalar err = 1.0;
	int loop;

	while (epsi > tol) {

		loop = 0;
		while (err > 0.9 * epsi && loop < 100) {
			loop++;

			// Set up Newton system
			XYZofLAMBDA(x);
			DualGrad(x);
			grad = -1.0 * grad - epsi / lam;
			DualHess(x);

			// Solve Newton system
			//if (m > 1) {
			//	Factorize(hess.data(), m);
			//	Solve(hess.data(), grad.data(), m);
			//	for (int j = 0; j < m; ++j) {
			//		s[j] = grad[j];
			//	}
			//} else if (m > 0) {
			s[0] = grad / hess;

			// Get the full search direction
			s[1] = -mu + epsi / lam - s[0] * mu / lam;

			// Perform linesearch and update lam and mu
			DualLineSearch();

			XYZofLAMBDA(x);

			// Compute KKT res
			err = DualResidual(x, epsi);
		}
		epsi = epsi * 0.1;
	}
}


void GCMMASolver::XYZofLAMBDA(Foam::volScalarField& x) {

	Foam::scalar lamai = 0.0;
	if (lam < 0.0) {
		lam = 0.0;
	}
	y = Foam::max(Foam::scalar(0.0), lam - c); // Note y=(lam-c)/d - however d is fixed at one !!
	lamai += lam * a;

	z = Foam::max(Foam::scalar(0.0), 10.0 * (lamai - 1.0)); // SINCE a0 = 1.0

	Foam::scalar pjlam = 0.0, qjlam = 0.0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		pjlam = p0[i] + pij[i] * lam;
		qjlam = q0[i] + qij[i] * lam;

		x[i] = (Foam::sqrt(pjlam) * low[i] + Foam::sqrt(qjlam) * upp[i]) / (Foam::sqrt(pjlam) + Foam::sqrt(qjlam));
		if (x[i] < alpha[i]) {
			x[i] = alpha[i];
		}
		if (x[i] > beta[i]) {
			x[i] = beta[i];
		}
	}
}

void GCMMASolver::DualGrad(Foam::volScalarField& x) {
	grad = -b - a * z - y;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		grad += pij[i] / (upp[i] - x[i]) + qij[i] / (x[i] - low[i]);
	}
	Foam::reduce(grad,sumOp<scalar>());
}

void GCMMASolver::DualHess(Foam::volScalarField& x) {

	Foam::volScalarField df2
	(
		IOobject
		(
		"df2",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
		mesh,
		Foam::dimensionedScalar("0", dimless, 0.0),
		zeroGradientFvPatchScalarField::typeName
	);
	Foam::volScalarField PQ
	(
		IOobject
		(
		"PQ",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
		mesh,
		Foam::dimensionedScalar("0", dimless, 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Foam::scalar pjlam, qjlam;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		pjlam = p0[i];
		qjlam = q0[i];
		pjlam += pij[i] * lam;
		qjlam += qij[i] * lam;
		PQ[i] = pij[i] / Foam::pow(upp[i] - x[i], 2.0) - qij[i] / Foam::pow(x[i] - low[i], 2.0);
		df2[i] = -1.0 / (2.0 * pjlam / Foam::pow(upp[i] - x[i], 3.0) + 2.0 * qjlam / Foam::pow(x[i] - low[i], 3.0));
		Foam::scalar xp = (Foam::sqrt(pjlam) * low[i] + Foam::sqrt(qjlam) * upp[i]) / (Foam::sqrt(pjlam) + Foam::sqrt(qjlam));
		if (xp < alpha[i]) {
			df2[i] = 0.0;
		}
		if (xp > beta[i]) {
			df2[i] = 0.0;
		}
	}

	// Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
	Foam::volScalarField tmp
	(
		IOobject
		(
		"tmp",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
		mesh,
		Foam::dimensionedScalar("0", dimless, 0.0),
		zeroGradientFvPatchScalarField::typeName
	);
	hess = 0.0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		tmp[i] += PQ[i] * df2[i];
		hess += tmp[i] * PQ[i];
	}
	Foam::reduce(hess,sumOp<scalar>());

	Foam::scalar lamai = 0.0;
	if (lam < 0.0) {
		lam = 0.0;
	}
	lamai += lam * a;
	if (lam > c) {
		hess += -1.0;
	}
	hess += -mu / lam;

	if (lamai > 0.0) {
		hess += -10.0 * a * a;
	}

	// pos def check
	Foam::scalar hessTrace = 0.0;
	hessTrace += hess;

	Foam::scalar hessCorr = 1e-4 * hessTrace;

	if (-1.0 * hessCorr < 1.0e-7) {
		hessCorr = -1.0e-7;
	}

	hess += hessCorr;
}

void GCMMASolver::SolveDSA(Foam::volScalarField& x) {

	lam = 1.0;

	const Foam::scalar tol = epsimin; // 1.0e-9*Foam::sqrt(m+n);
	Foam::scalar err = 1.0;
	int loop = 0;

	while (err > tol && loop < 500) {
		loop++;
		XYZofLAMBDA(x);
		DualGrad(x);
		Foam::scalar theta = 1.0;
		err = 0.0;
		lam = Foam::max(Foam::scalar(0.0), lam + theta * grad);
		err += grad * grad;

		err = Foam::sqrt(err);
	}
}

Foam::scalar GCMMASolver::DualResidual(Foam::volScalarField& x, Foam::scalar& epsi) {

	Foam::List<Foam::scalar> res(2);

	res[0] = -b - a * z - y + mu;
	res[1] = mu * lam - epsi;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		res[0] += pij[i] / (upp[i] - x[i]) + qij[i] / (x[i] - low[i]);
	}
	Foam::reduce(res[0],sumOp<scalar>());
	
	Foam::scalar nrI = 0.0;
	for (int i = 0; i < 2; ++i) {
		if (nrI < Foam::mag(res[i])) {
			nrI = Foam::mag(res[i]);
		}
	}

	return nrI;
}

void GCMMASolver::DualLineSearch() {

	Foam::scalar theta = 1.005;
	if (theta < -1.01 * s[0] / lam) {
		theta = -1.01 * s[0] / lam;
	}
	if (theta < -1.01 * s[1] / mu) {
		theta = -1.01 * s[1] / mu;
	}
	theta = 1.0 / theta;

	lam = lam + theta * s[0];
	mu = mu + theta * s[1];
}

// Update [raa0, raa]
void GCMMASolver::RaaUpdate(const Foam::volScalarField& xmma, const Foam::volScalarField& xval, Foam::scalar& f0xnew, const Foam::scalar& fxnew) //  const double *xmin, const double *xmax)
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
	// cout << "raacof: " << raacof << endl;

	if (f0xnew > f0app + 0.5 * epsimin) {
		Foam::scalar deltaraa0 = (1.0 / raacof) * (f0xnew - f0app);
		raa0 = Foam::min(1.1 * (raa0 + deltaraa0), 10.0 * raa0);
	}
	// cout << "raa0: " << raa0 << endl << "raaj: ";
	if (fxnew > fapp + 0.5 * epsimin) {
		Foam::scalar deltaraa = (1.0 / raacof) * (fxnew - fapp);
		raa = Foam::min(1.1 * (raa + deltaraa), 10.0 * raa);
	}
		// cout << raa[j] << ' ';
	// cout << endl;
}


void GCMMASolver::GenSub(const Foam::volScalarField& xval, Foam::scalar& f0x, const Foam::volScalarField& df0dx, const Foam::scalar& fx, const Foam::volScalarField& dfdx)
	//, const double *xmin, const double *xmax)
{
	// Set bounds and the coefficients for the approximation
	r0 = 0.0;
	r = 0.0;
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
			Foam::scalar dfdxp = Foam::max(Foam::scalar(0.0), dfdx[i]);
			Foam::scalar dfdxm = Foam::max(Foam::scalar(0.0), -1.0 * dfdx[i]);
			Foam::scalar xmamiinv = 1.0 / Foam::max(xmamieps, 1.0);
			Foam::scalar pq = 0.001 * Foam::mag(dfdx[i]) + raa * xmamiinv;
			pij[i] = Foam::pow(upp[i] - xval[i], 2.0) * (dfdxp + pq);
			qij[i] = Foam::pow(xval[i] - low[i], 2.0) * (dfdxm + pq);
			r += pij[i] * uxinv + qij[i] * xlinv;
		}
	}
	Foam::reduce(r0,sumOp<scalar>());
	Foam::reduce(r,sumOp<scalar>());
	r0 = f0x - r0;
	r = fx - r;
	b = -r; // The constant for the constraints
}

void GCMMASolver::ComputeApprox(const Foam::volScalarField& xmma) {
	f0app = 0.0;
	fapp = 0.0;
	forAll(designSpaceCells,jj) {
		const Foam::label i = designSpaceCells[jj];
		Foam::scalar uxinv = 1.0 / (upp[i] - xmma[i]);
		Foam::scalar xlinv = 1.0 / (xmma[i] - low[i]);
		f0app += p0[i] * uxinv + q0[i] * xlinv;
		fapp += pij[i] * uxinv + qij[i] * xlinv;
	}
	Foam::reduce(f0app,sumOp<scalar>());
	Foam::reduce(fapp,sumOp<scalar>());
	f0app += r0;
	fapp += r;
}

/*
void GCMMASolver::Factorize(double *K, int n) {

	for (int s = 0; s < n - 1; ++s) {
		for (int i = s + 1; i < n; ++i) {
			K[i * n + s] = K[i * n + s] / K[s * n + s];
			for (int j = s + 1; j < n; ++j) {
				K[i * n + j] = K[i * n + j] - K[i * n + s] * K[s * n + j];
			}
		}
	}
}

void GCMMASolver::Solve(double *K, double *x, int n) {

	for (int i = 1; i < n; ++i) {
		double a = 0.0;
		for (int j = 0; j < i; ++j) {
			a = a - K[i * n + j] * x[j];
		}
		x[i] = x[i] + a;
	}

	x[n - 1] = x[n - 1] / K[(n - 1) * n + (n - 1)];
	for (int i = n - 2; i >= 0; --i) {
		double a = x[i];
		for (int j = i + 1; j < n; ++j) {
			a = a - K[i * n + j] * x[j];
		}
		x[i] = a / K[i * n + i];
	}
}
*/
