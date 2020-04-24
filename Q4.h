/*****************************************************************************/
/*     To acheive a plane strain element                                     */
/*     WangShuai                                                             */
/*     2020/3/14                                                             */
/*****************************************************************************/
#pragma once

#include "Element.h"

using namespace std;

//Plane strain element class

class CQ4 : public CElement
{
private:
	
	//! Calculate the elastic matrix
	double D[3][3];
	void Calculate_D(CQ4Material* material_,double D[][3]);

	double B[3][8];

	//! The number of Gauss point in one direction
	static int NG;

	int NGauss;

	//! Calculate B matrix and Jaccobi matrix
	void STDM(CNode** nodes_, double B[][8], double Jac, double XX, double YY);


public:


	//!	Constructor
	CQ4();

	//!	Deconstructor
	~CQ4();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

	//! Generate location matrix: the global equation number that corresponding to each DOF of the element
	//	Caution:  Equation number is numbered from 1 !
	virtual void GenerateLocationMatrix();

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

	//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();

	//! Get the number of Gauss point in one direction
	virtual void GetNG(int N_G) { NG = N_G; };

	//! Get the number of nodes in one direction
	virtual int GetNN() { return NEN_; };

};
