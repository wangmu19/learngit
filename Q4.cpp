/*****************************************************************************/
/*     To acheive a plane strain element                                     */
/*     WangShuai                                                             */
/*     2020/3/14                                                             */
/*****************************************************************************/
#pragma once

#include"Q4.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int CQ4::NG;

//Constructor
CQ4::CQ4()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//Deconstructor
CQ4::~CQ4()
{
}

//!	Read element data from stream Input
bool CQ4::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// The node number in anticlockwise

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//!	Write element data to stream
void CQ4::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber
		<< setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CQ4::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 3; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//!	Calculate element stiffness matrix
void CQ4::ElementStiffness(double* Matrix)
{
	int n = SizeOfStiffnessMatrix();
	clear(Matrix, n);

	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);
	//Calculate the elastric matrix
	Calculate_D(material_,D);

	double KK[8][8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++)
			KK[i][j] = 0;
	}
	double thic = material_->t;
	
	double SS[3]; //! the column of SS = D * B

	for (int Gx = 0; Gx < NG; Gx++) {
		double XI = XG[Gx][NG - 1];
		for (int Gy = 0; Gy < NG; Gy++) {
			double YI = XG[Gy][NG - 1];
			double WT = WGT[Gx][NG-1] * WGT[Gy][NG-1];
			//! Calculate B matrix and Jaccobi matrix
			STDM(nodes_, B, det_J, XI, YI);

			//! Calculate the element stiffness matrix
			for (int j = 0; j < 8; j++) {
				for (int k = 0; k < 3; k++) {
					SS[k] = 0.0;
					for (int l = 0; l < 3; l++) {
						SS[k] += D[k][l] * B[l][j];
					}
				}
				for (int i = j; i < 8; i++) {
					double STIFF = 0.0;
					for (int l = 0; l < 3; l++) {
						STIFF += B[l][i] * SS[l];
					}
					KK[j][i] += STIFF * WT * thic * det_J;
				}
			}
		}
	}

	//! Extend KK into 12*12 matrix
	double K[12][12];
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			K[i][j] = 0.0;
		}
	}
	int m1, m2;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			int k1 = i / 2;
			if (i % 2 == 0)
				m1 = 0;
			else
				m1 = 1;
			int k2 = j / 2;
			if (j % 2 == 0)
				m2 = 0;
			else
				m2 = 1;
			K[3 * k1 + m1][3 * k2 + m2] = KK[i][j];

		}
	}

	//! Store K into Matrix
	int k = 0;
	for (int i = 0; i < 12; i++) {
		for (int j = i; j >= 0; j--) {
			Matrix[k] = K[j][i];
			k += 1;
		}
	}
}

//!	Calculate element stress
void CQ4::ElementStress(double* stress, double* Displacement)
{
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);
	Calculate_D(material_, D);
	// Make up the whole displacement vector
	double dis[8];
	int k = 0;
	for (int i = 0; i < NEN_; i++) {
		for (unsigned int j = 0; j < 2; j++)
		{
			if (nodes_[i]->bcode[j] == 0)
			{
				dis[k] = 0;
			}
			else
			{
				dis[k] = Displacement[nodes_[i]->bcode[j] - 1];
			}
			k += 1;
		}
	}
	k = 0;
	double stan[3], strs[3]; // The strain and stress
	for (int Gx = 0; Gx < NG; Gx++) {
		double XI = XG[Gx][NG - 1];
		for (int Gy = 0; Gy < NG; Gy++) {
			for (int i = 0; i < 3; i++) {
				stan[i] = 0.0;
				strs[i] = 0.0;
			}
			double YI = XG[Gy][NG - 1];
			//! Calculate B matrix and Jaccobi matrix
			STDM(nodes_, B, det_J, XI, YI);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 8; j++) {
					stan[i] += B[i][j] * dis[j];
				}
			}
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					strs[i] += D[i][j] * stan[j];
				}
			}
			for (int i = 0; i < 3; i++) {
				stress[3 * k + i] = strs[i];
			}
			k += 1;
		}
	}
}

//!	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4 node Q4 element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CQ4::SizeOfStiffnessMatrix()
{
	return 78;
}

//! Calculate the elastic matrix
void CQ4:: Calculate_D(CQ4Material* material, double D[][3])
{
	if (material->ss == 0) {}
	else if (material->ss == 1) {
		material->E = material->E * (1 + 2 * material->miu) / (1 + material->miu) / (1 + material->miu);
		material->miu = material->miu / (1 + material->miu);
	}
	else {
		cerr << "Q4 has a wrong input ";
	}

	double F = material->E / (1.0 + material->miu);
	double G = F / (1.0 - 2.0 * material->miu) * material->miu;
	double H = F + G;
	D[0][0] = H;
	D[0][1] = G;
	D[0][2] = 0;
	D[1][0] = G;
	D[1][1] = H;
	D[1][2] = 0;
	D[2][0] = 0;
	D[2][1] = 0;
	D[2][2] = F / 2.0;
}

//! Calculate B matrix and Jaccobi matrix
void CQ4::STDM(CNode** nodes_, double B[][8], double Jac, double R, double S)
{
	//! Construct the shape function
	double RP = 1.0 + R;
	double SP = 1.0 + S;
	double RM = 1.0 - R;
	double SM = 1.0 - S;
	double N[4];
	N[0] = 0.25 * RP * SP;
	N[1] = 0.25 * RM * SP;
	N[2] = 0.25 * RM * SM;
	N[3] = 0.25 * RP * SM;

	//! Calculate the derivative of the shape function
	//! The first row is by R and second row by S
	double P[2][4];
	P[0][0] = 0.25 * SP;
	P[0][1] = -P[0][0];
	P[0][2] = -0.25 * SM;
	P[0][3] = -P[0][2];
	P[1][0] = 0.25* RP;
	P[1][1] = 0.25* RM;
	P[1][2] = -P[1][1];
	P[1][3] = -P[1][0];

	//! Calculate the determinant of Jaccobi matrix
	double XJ[2][2];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			double dum = 0;
			for (int k = 0; k < 4; k++) {
				dum += P[i][k] * nodes_[k]->XYZ[j];
			}
			XJ[i][j] = dum;
		}
	}
	det_J = XJ[0][0] * XJ[1][1] - XJ[1][0] * XJ[0][1];
	if (det_J < 1e-7) {
		cout << "*** Error *** Some element is singular";
		exit(0);
	}

	//! Calculate the inverse of Jaccobi matrix
	double det_JI = 1.0 / det_J; 
	double XJI[2][2];
	XJI[0][0] = XJ[1][1] * det_JI;
	XJI[0][1] = -XJ[0][1] * det_JI;
	XJI[1][0] = -XJ[1][0] * det_JI;
	XJI[1][1] = XJ[0][0] * det_JI;

	//! Calculate the B matrix
	int k2 = -1;
	for (int k = 0; k < 4; k++) {
		k2 += 2;
		B[0][k2 - 1] = 0.0;
		B[0][k2] = 0.0;
		B[1][k2 - 1] = 0.0;
		B[1][k2] = 0.0;
		for (int i = 0; i < 2; i++) {
			B[0][k2 - 1] += XJI[0][i] * P[i][k];
			B[1][k2] += XJI[1][i] * P[i][k];
		}
		B[2][k2] = B[0][k2 - 1];
		B[2][k2 - 1] = B[1][k2];
	}

}