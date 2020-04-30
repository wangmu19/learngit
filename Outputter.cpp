/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;
COutputter* COutputter::tec_instance = nullptr;
COutputter* COutputter::vtk_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::Instance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
	return _instance;
}
//! Return the tecplot instance of the class
COutputter* COutputter::Tec_Instance(string FileName)
{
	if (!tec_instance)
		tec_instance = new COutputter(FileName);
	return tec_instance;
}
COutputter* COutputter::vtk_Instance(string FileName)
{
	if(!vtk_instance)
		vtk_instance = new COutputter(FileName);
	return vtk_instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::Instance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this, np);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::Instance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this, np);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::Instance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, Q4 ELEMENTS" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				PrintBarElementData(EleGrp);
				break;
			case ElementTypes::Q4: // Bar element
				PrintQ4ElementData(EleGrp);
				break;
		}
	}
}
//	Output bar element data
void COutputter::PrintBarElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup[Ele].Write(*this, Ele);

	*this << endl;
}

//	Output bar element data
void COutputter::PrintQ4ElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S     POISSON'S     THICKNESS" << endl
		<< " NUMBER     MODULUS           RATIO                    " << endl
		<< "               E               v           t     " << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		<< endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE       MATERIAL" << endl
		<< " NUMBER-N      I        J         K         L       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup[Ele].Write(*this, Ele);

	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::Instance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this, lcase);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement(unsigned int lcase)
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << " LOAD CASE" << setw(5) << lcase + 1 << endl
		  << endl
		  << endl;

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, np, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();
		unsigned int NG = EleGrp.GetNG();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);

					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;

			case ElementTypes::Q4: // Q4 element
				
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double* stress = new double[NG * NG * 3];
					*this << "Element Number" << setw(5) << Ele + 1 << endl;
					//OutputFile << "Element Number" << setw(5) << Ele + 1 << endl;
					Element.ElementStress(stress, Displacement);
					*this << "Stress on GAUSS Points" << endl;
					//OutputFile << "Stress on GAUSS Points" << endl;
					for (int i = 0; i < NG * NG; i++) {
						*this << "Stress on  Point" << setw(5) << i + 1 << endl;

						*this << setw(20) << stress[3 * i + 0] << setw(20) << stress[3 * i + 1] << setw(20) << stress[3 * i + 2] << endl;

					}
					delete[] stress;
				}




				*this << endl;

				break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

// Output into tecplot
void COutputter::OutputTecplot(int step)
{
	CDomain* FEMData = CDomain::Instance();
	double* dis = FEMData->GetDisplacement();
	CNode* NodeList = FEMData->GetNodeList();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUME = 0;
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		NUME += EleGrp.GetNUME();
	}

	if (step == 0) {
		OutputFile << "TITLE = " << "\"" << FEMData->GetTitle() << "\"" << endl;
		OutputFile << "VARIABLES = ";
		OutputFile << setw(6) << "\"X\"" << setw(6) << "\"Y\"" << setw(6) << "\"Z\"" << setw(6) << "\"DIS\"" << endl;
		OutputFile << endl;
	}
	
	OutputFile << setiosflags(ios::right) << setiosflags(ios::fixed);
	OutputFile << "ZONE T=\"Time = " << setw(12) << setprecision(4) << (double)step << "\"" << " F=FEPOINT "<< "N=" << setw(5) << NUMNP << " E=" << setw(5) << NUME << " ET=QUADRILATERAL C=CYAN" << endl;

	// Calculate the position of nodes after simulation. It will not be done here in the dynamics situation
	double* dis_vector = new double[3 * NUMNP];
	int k = 0;
	for (int i = 0; i < 3 * NUMNP; i++)
		dis_vector[i] = 0.0;
	if (step > 0) {
		for (int i = 0; i < NUMNP; i++) {
			for (int j = 0; j < 3; j++) {
				if (NodeList[i].bcode[j] != 0) {
					dis_vector[k] = dis[NodeList[i].bcode[j] - 1] * 1000;		// A 10000 time of the displacement is plotted
					NodeList[i].XYZ[j] += dis_vector[k] ;
				}
				k += 1;
			}
		}
	}

	for (int i = 0; i < NUMNP; i++) {
		double disp = sqrt(dis_vector[3 * i]* dis_vector[3 * i]+ dis_vector[3 * i + 1] * dis_vector[3 * i + 1]+ dis_vector[3 * i + 2] * dis_vector[3 * i + 2]);
		OutputFile << setw(12) << setprecision(4) << NodeList[i].XYZ[0] << setw(12) << setprecision(4) << NodeList[i].XYZ[1] << setw(12) << setprecision(4) << NodeList[i].XYZ[2] << setw(12) << setprecision(4) << disp << endl;
	}
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++) {
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		int Num_Node = EleGrp[0].GetNN();
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			for (int N_n = 0; N_n < Num_Node; N_n++)
			{
				OutputFile << setw(12) << EleGrp[Ele].GetNodes()[N_n]->NodeNumber;
			}
			OutputFile << endl;
		}
	}
}

void COutputter::OutputVTK()//need a reload (double time,double* dis)
{
	OutputVTKNodalDis();
	OutputVTKElemStress();
}

void COutputter::OutputVTK(double time, double* dis)
{
	OutputVTKNodalDis(time, dis);
	OutputVTKElemStress(time, dis);
}

void COutputter::OutputVTKNodalDis(double time, double* dis)
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	unsigned int NUMNP = FEMData->GetNUMNP();

//	allocat a matrix: nodepoint*3 
	double **Dis = new double *[NUMNP];
	for ( int i = 0; i< NUMNP; i++)
		Dis[i]= new double [3];
	int count = 0 ;

	for ( int i =0; i < NUMNP; i++)
	{
		for ( int j = 0; j < 3; j++)
		{
			if( NodeList[i].bcode[j] == 0 )
				Dis[i][j] = 0;
			else
			{
				Dis[i][j] = dis[count];
				count++;
			}
		}
	}

	OutputFile << "POINT_DATA" << setw(8) << NUMNP << endl;
	OutputFile << setiosflags(ios::left);
	OutputFile << "VECTORS DISPLACEMENT_" << setw(8) << time <<" double" << endl;
	OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
	for (unsigned int np = 0; np < NUMNP; np++)
		OutputFile << setw(8) << Dis[np][0] << "    " << setw(8) << Dis[np][1] << "    " << setw(8) << Dis[np][2] << endl;
	OutputFile << endl;
	//delete Dis
	for (int i = 0; i < NUMNP; i++)
		delete []Dis[i];
	delete []Dis;
}

void COutputter::OutputVTKElemStress(double time, double *dis)
{
	CDomain* FEMData = CDomain::Instance();

	unsigned int NUMEG = FEMData->GetNUMEG();
	// count number of all elements
	int NUMEALL = 0;
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		NUMEALL +=EleGrp.GetNUME();
	}
	//Output Elemstress Head Info
	if( (NUMEG!=1) || !(ElementTypes::Bar-1))
	{
		OutputFile << "CELL_DATA " << setw(8) << NUMEALL << endl;
		//Tensor
		OutputFile << setiosflags(ios::left);
		OutputFile << "TENSORS " << "Element_Stress" << setw(8) << time << " double" << endl;
	}

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();
		unsigned int NG = EleGrp.GetNG();
	
		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				if( NUMEG == 1)
				{
// new stress matrix, number of elments*2,column[0] for stress and column[1] for force
					double** stress = new double*[NUME];
					for (int i = 0; i < NUME; i++)
						stress[i] = new double [2];	

					for (unsigned int Ele = 0; Ele < NUME; Ele++)
					{
						CElement& Element = EleGrp[Ele];
						Element.ElementStress(&stress[Ele][0], dis);
						CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
						stress[Ele][1] = stress[Ele][0] * material.Area;
					}
				
					// Output element stress
					OutputFile << "CELL_DATA " << setw(8) << NUME << endl;
					OutputFile << setiosflags(ios::left);
					OutputFile << "SCALARS " << "Stress" << setw(8) << time << " double " << "1" << endl;
					OutputFile << "LOOKUP_TABLE " << "Element_Stress" << setw(8) << time <<endl;
					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
					for (int i = 0; i < NUME; i++)
						OutputFile << setw(8) << stress[i][0] << endl;
					OutputFile << endl;

					// Output element force
					OutputFile << setiosflags(ios::left);
					OutputFile << "SCALARS " << "Force" << setw(8) << time << " double " << "1" << endl;
					OutputFile << "LOOKUP_TABLE " << "Element_Force" << setw(8) << time << endl;
					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
					for (int i = 0; i < NUME; i++)
						OutputFile << setw(8) << stress[i][1] << endl;
					OutputFile << endl;
				//delete stress
				for (int i = 0; i < NUME; i++)
					delete []stress[i];
				delete []stress;
				}
				else
				{
					double stress;
					for (unsigned int Ele = 0; Ele < NUME; Ele++)
					{
						CElement& Element = EleGrp[Ele];
						Element.ElementStress(&stress, dis);
						OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
						OutputFile << setw(8) << stress << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
						OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
						OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
						OutputFile << endl;
					}
				}
				break;

			case ElementTypes::Q4: // Q4 element

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double* stress = new double[NG * NG * 3];
					Element.ElementStress(stress, dis);
				// averagestress
					double averagestress[3]={0};

					for (int i = 0; i < NG * NG; i++) 
						for (int j = 0; j < 3; j++)
							averagestress[j] += stress[ 3 * i + j];

					delete[] stress;
				//output averagestress
					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
					OutputFile << setw(8) << averagestress[0] << " " << setw(8) << averagestress[2] << " " << setw(8) << 0 << endl;
					OutputFile << setw(8) << averagestress[2] << " " << setw(8) << averagestress[1] << " " << setw(8) << 0 << endl;
					OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
					OutputFile << endl;
				}
				break;

//			case ElementTypes::newelement:// new element
//				for (unsigned int Ele = 0; Ele < NUME; Ele++)
//				{
//					CElement& Element = EleGrp[Ele];
//					double* stress = new double[NG * NG * 3];
//					Element.ElementStress(stress, dis);
				// averagestress
//					double averagestress[3]={0};

//					for (int i = 0; i < NG * NG; i++) 
//						for (int j = 0; j < 3; j++)
//							averagestress[j] += stress[ 3 * i + j];

//					delete[] stress;
				//output averagestress
//					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
//					OutputFile << setw(8) << averagestress[0] << " " << setw(8) << averagestress[2] << " " << setw(8) << 0 << endl;
//					OutputFile << setw(8) << averagestress[2] << " " << setw(8) << averagestress[1] << " " << setw(8) << 0 << endl;
//					OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
//					OutputFile << endl;
//				}
//				break;
		}
	}

}

void COutputter::OutputVTKHead()
{
	CDomain* FEMData = CDomain::Instance();
	OutputFile << "# vtk DataFile Version 3.0" <<endl;
	OutputFile << "\"" << FEMData->GetTitle() << "\"" << endl;
	OutputFile << "ASCII" << endl;
	OutputFile << "DATASET UNSTRUCTURED_GRID" << endl
			   << endl;
}

void COutputter::OutputVTKNodes()
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	unsigned int NUMNP = FEMData->GetNUMNP();

	OutputFile << "POINTS" << setw(8) << NUMNP << " double" << endl;
	OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
	for (unsigned int np = 0; np < NUMNP; np++)
		OutputFile << NodeList[np].XYZ[0] << " " << setw(8) << NodeList[np].XYZ[1] << " " << setw(8) << NodeList[np].XYZ[2] << endl;
}

void COutputter::OutputVTKElements()
{
	CDomain* FEMData = CDomain::Instance();
	unsigned int NUMEG = FEMData->GetNUMEG();
	int NUMEALL = 0;//numbers of all types of element
	int nlist = 0;	//nlist for vtk file
	int **ElemInfo = new int *[NUMEG];//allocate a matrix: elment type* nums of this type
	for ( int i = 0; i< NUMEG; i++)
		ElemInfo[i]= new int [2];
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();
		ElemInfo[EleGrp][1] = NUME;
		NUMEALL += NUME;
		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				ElemInfo[EleGrp][0] = 3;
				nlist += 3*NUME;
				break;
			case ElementTypes::Q4: // Q4 element
				ElemInfo[EleGrp][0] = 9;
				nlist += 5*NUME;
				break;
//			case ElementTypes::NewType // NewType 
//				ElemInfo[EleGrp][0] = 23;
//				nlist += 9*NUME;
//				break;
		}
	}
	OutputFile << "CELLS" << setw(8) << NUMEALL << setw(8) << nlist << endl;
	//write cell info
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++) {
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrpIndex].GetNUME();
		int Num_Node = EleGrp[0].GetNN();
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			OutputFile << Num_Node;
			for (int N_n = 0; N_n < Num_Node; N_n++)
			{
				OutputFile << setw(8) << EleGrp[Ele].GetNodes()[N_n]->NodeNumber-1;
			}
			OutputFile << endl;
		}
	}
	OutputFile << "CELL_TYPES" << setw(8) << NUMEALL << endl;
	//write cell_type
	for (int i = 0; i < NUMEG; i++) {
		for (int j =0 ; j< ElemInfo[i][1]; j++){
			OutputFile << ElemInfo[i][0] << endl;
		}
	}
	//delete ElemInfo
	for (int i = 0; i < NUMEG; i++)
		delete []ElemInfo[i];
	delete []ElemInfo;
}

void COutputter::OutputVTKNodalDis()
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();
	unsigned int NUMNP = FEMData->GetNUMNP();

//	allocat a matrix: nodepoint*3 
	double **Dis = new double *[NUMNP];
	for ( int i = 0; i< NUMNP; i++)
		Dis[i]= new double [3];
	int count = 0 ;

	for ( int i =0; i < NUMNP; i++)
	{
		for ( int j = 0; j < 3; j++)
		{
			if( NodeList[i].bcode[j] == 0 )
				Dis[i][j] = 0;
			else
			{
				Dis[i][j] = Displacement[count];
				count++;
			}
		}
	}

	OutputFile << "POINT_DATA" << setw(8) << NUMNP << endl;
	OutputFile << "VECTORS DISPLACEMENT double" << endl;
	OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
	for (unsigned int np = 0; np < NUMNP; np++)
		OutputFile << setw(8) << Dis[np][0] << "    " << setw(8) << Dis[np][1] << "    " << setw(8) << Dis[np][2] << endl;
	OutputFile << endl;
	//delete Dis
	for (int i = 0; i < NUMNP; i++)
		delete []Dis[i];
	delete []Dis;
}

void COutputter::OutputVTKElemStress()
{
	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();
	// count number of all elements
	int NUMEALL = 0;
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		NUMEALL +=EleGrp.GetNUME();
	}
	//Output Elemstress Head Info
	if( (NUMEG!=1) || !(ElementTypes::Bar-1))
	{
		OutputFile << "CELL_DATA " << setw(8) << NUMEALL << endl;
		//Tensor
		OutputFile << "TENSORS " << "Element_Stress " << "double" << endl;
	}

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();
		unsigned int NG = EleGrp.GetNG();
	
		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				if( NUMEG == 1)
				{
// new stress matrix, number of elments*2,column[0] for stress and column[1] for force
					double** stress = new double*[NUME];
					for (int i = 0; i < NUME; i++)
						stress[i] = new double [2];	

					for (unsigned int Ele = 0; Ele < NUME; Ele++)
					{
						CElement& Element = EleGrp[Ele];
						Element.ElementStress(&stress[Ele][0], Displacement);
						CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
						stress[Ele][1] = stress[Ele][0] * material.Area;
					}
				
					// Output element stress
					OutputFile << "CELL_DATA " << setw(8) << NUME << endl;
					OutputFile << "SCALARS " << "Stress " << "double " << "1" << endl;
					OutputFile << "LOOKUP_TABLE " << "Element_Stress" << endl;
					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
					for (int i = 0; i < NUME; i++)
						OutputFile << setw(8) << stress[i][0] << endl;
					OutputFile << endl;

					// Output element force
					OutputFile << "SCALARS " << "Force " << "double " << "1" << endl;
					OutputFile << "LOOKUP_TABLE " << "Element_Force" << endl;
					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
					for (int i = 0; i < NUME; i++)
						OutputFile << setw(8) << stress[i][1] << endl;
					OutputFile << endl;
				//delete stress
				for (int i = 0; i < NUME; i++)
					delete []stress[i];
				delete []stress;
				}
				else
				{
					double stress;
					for (unsigned int Ele = 0; Ele < NUME; Ele++)
					{
						CElement& Element = EleGrp[Ele];
						Element.ElementStress(&stress, Displacement);
						OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
						OutputFile << setw(8) << stress << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
						OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
						OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
						OutputFile << endl;
					}
				}
				break;

			case ElementTypes::Q4: // Q4 element

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double* stress = new double[NG * NG * 3];
					Element.ElementStress(stress, Displacement);
				// averagestress
					double averagestress[3]={0};

					for (int i = 0; i < NG * NG; i++) 
						for (int j = 0; j < 3; j++)
							averagestress[j] += stress[ 3 * i + j];

					delete[] stress;
				//output averagestress
					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
					OutputFile << setw(8) << averagestress[0] << " " << setw(8) << averagestress[2] << " " << setw(8) << 0 << endl;
					OutputFile << setw(8) << averagestress[2] << " " << setw(8) << averagestress[1] << " " << setw(8) << 0 << endl;
					OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
					OutputFile << endl;
				}
				break;

//			case ElementTypes::newelement:// new element
//				for (unsigned int Ele = 0; Ele < NUME; Ele++)
//				{
//					CElement& Element = EleGrp[Ele];
//					double* stress = new double[NG * NG * 3];
//					Element.ElementStress(stress, Displacement);
				// averagestress
//					double averagestress[3]={0};

//					for (int i = 0; i < NG * NG; i++) 
//						for (int j = 0; j < 3; j++)
//							averagestress[j] += stress[ 3 * i + j];

//					delete[] stress;
				//output averagestress
//					OutputFile << setiosflags(ios::right) << setiosflags(ios::scientific);
//					OutputFile << setw(8) << averagestress[0] << " " << setw(8) << averagestress[2] << " " << setw(8) << 0 << endl;
//					OutputFile << setw(8) << averagestress[2] << " " << setw(8) << averagestress[1] << " " << setw(8) << 0 << endl;
//					OutputFile << setw(8) << 0 << " " << setw(8) << 0 << " " << setw(8) << 0 << endl;
//					OutputFile << endl;
//				}
//				break;
		}
	}

}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int H = DiagonalAddress[J] - DiagonalAddress[J - 1];
			if (J - I - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I, J);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement(unsigned int loadcase)
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << "  Load case = " << loadcase << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
