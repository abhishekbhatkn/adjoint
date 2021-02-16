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
//
// BETA VERSION  0.99
//
// GCMMA solver using a dual interior point method
//
// Original MMA code by Niels Aage, February 2013
// Extension to GCMMA by Jérémie Dumas, June 2017
//
// The class solves a general non-linear programming problem
// on standard from, i.e. non-linear objective f, m non-linear
// inequality constraints g and box constraints on the n
// design variables xmin, xmax.
//
//        min_x^n f(x)
//        s.t. g_j(x) < 0,   j = 1,m
//        xmin < x_i < xmax, i = 1,n
//
// Each call to Update() sets up and solve the following
// convex subproblem:
//
//   min_x     sum(p0j./(U-x)+q0j./(x-L)) + a0*z + sum(c.*y + 0.5*d.*y.^2)
//
//   s.t.      sum(pij./(U-x)+qij./(x-L)) - ai*z - yi <= bi, i = 1,m
//             Lj < alphaj <=  xj <= betaj < Uj,  j = 1,n
//             yi >= 0, i = 1,m
//             z >= 0.
//
// NOTE: a0 == 1 in this implementation !!!!
//
////////////////////////////////////////////////////////////////////////////////

//#pragma once
#include "fvCFD.H"
//#include <vector>

class GCMMASolver {

public:
	GCMMASolver(Foam::fvMesh &me, Foam::List<label>& designCells, Foam::scalar cost, Foam::scalar constraint);
	void eval();

public:	
	

private:
	Foam::fvMesh& mesh;
	//const Foam::wordList designSpaceZones;
	Foam::List<label>& designSpaceCells;
	int n, outeriter;
	
	Foam::volScalarField dummy;
	
	Foam::volScalarField xold1, xold2;
	
	Foam::scalar cost, constraint;

	const Foam::scalar raa0eps, raaeps, xmamieps, epsimin, move, albefa;

	Foam::scalar asyminit, asymdec, asyminc;

	Foam::scalar raa0, raa;

	Foam::scalar a, c, d, y, z;

	Foam::List<scalar> s;
	Foam::scalar lam, mu, b, grad, hess, r, r0, fapp, f0app;
	
	Foam::volScalarField low, upp, alpha, beta, p0, q0, pij, qij;
	
private:

	void SetAsymptotes(Foam::scalar init, Foam::scalar decrease, Foam::scalar increase);

	// Compute [L, U, raa0, raa], build and solve GCMMA subproblem, compute [f0app, fapp]
	void OuterUpdate(Foam::volScalarField xmma, const Foam::volScalarField xval, Foam::scalar f0x, const Foam::volScalarField df0dx, const Foam::scalar fx, const Foam::volScalarField dfdx);
	//, const double *xmin, const double *xmax);

	// Update [raa0, raa], build and solve GCMMA subproblem, compute [f0app, fapp]
	void InnerUpdate(Foam::volScalarField  xmma, Foam::scalar f0xnew, const Foam::scalar fxnew, const Foam::volScalarField  xval, Foam::scalar f0x,
                              const Foam::volScalarField df0dx, const Foam::scalar fx, const Foam::volScalarField dfdx); //, const double *xmin, const double *xmax);

	// Check whether the new solution is conservative
	bool ConCheck(Foam::scalar f0xnew, const Foam::scalar fxnew) const;

	void Reset() { outeriter = 0; };
	
	// Compute [low, upp, raa0, raa]
	void Asymp(const Foam::volScalarField xval, const Foam::volScalarField df0dx, const Foam::volScalarField dfdx); //, const double *xmin, const double *xmax);

	// Update [raa0, raa]
	void RaaUpdate(const Foam::volScalarField xmma, const Foam::volScalarField xval, Foam::scalar f0xnew, const Foam::scalar fxnew); //, const double *xmin, const double *xmax);

	// Build CGMMA subproblem
	void GenSub(const Foam::volScalarField xval, Foam::scalar f0x, const Foam::volScalarField df0dx, const Foam::scalar fx, const Foam::volScalarField dfdx); //, const double *xmin, const double *xmax);

	// Compute [f0app, fapp]
	void ComputeApprox(const Foam::volScalarField xmma);

	void SolveDSA(Foam::volScalarField x);
	void SolveDIP(Foam::volScalarField x);

	void XYZofLAMBDA(Foam::volScalarField x);

	void DualGrad(Foam::volScalarField x);
	void DualHess(Foam::volScalarField x);
	void DualLineSearch();
	Foam::scalar DualResidual(Foam::volScalarField x, Foam::scalar epsi);

	//static void Factorize(double *K, int n);
	//static void Solve(double *K, double *x, int n); */
};
