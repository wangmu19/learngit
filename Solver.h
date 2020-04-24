/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "SkylineMatrix.h"
#include <vector>
#include <algorithm>


//!	Base class for a solver
/*	New solver should be derived from this base class, and match the storage scheme
	of the global stiffness matrix employed in Domain class. */
class CSolver
{
protected:

	CSkylineMatrix<double>* K;

public:

	CSolver(CSkylineMatrix<double>* K);
    
};

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver : public CSolver
{
public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K) : CSolver(K) {};

//!	Perform L*D*L(T) factorization of the stiffness matrix
	void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	void BackSubstitution(double* Force); 
};

//! Modal solver: modal solver with Lanczos method
class CModal : public CSolver
{
private:

	//! The tolerance for Jacobi
	double Tol_J = 1.0e-9;

	//! The tolerance for Lanczos
	double Tol_L = 1.0e-9;

	//! Number of required eigenvalues *(Input)
	int Nroot = 2;

	//! Maximum number of the restart, set 5 usually
	int N_ite_max = 5;

	//! Number of itegration vectors used
	//! Usually set to be min(2*Nroot, Nroot+8), but less than the freedom of system
	int N_iv ;

	//! The eigenvalue vector
	double* Eig_value;

	//! The eigenvector matrix
	std::vector<double>* Eig_vector;



public:

	//! Constructor
	CModal(CSkylineMatrix<double>* K) : CSolver(K) {};

	//! To solve the generalized eigenproblem with the Lanczos method
	void Lanczos();

	//! To solve the standard eigenproblem with the Jacobi method
	void Jacobi();

	//! Use the Gram-Schmidt orthogonalization
	void Orth();


};