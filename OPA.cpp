/***********************************************************/
/***********************************************************/
//             OPA program C++ version 1.0				   //
/***********************************************************/
/***********************************************************/
// Compilation: g++ -lgsl -lgslcblas -lfftw3 OPAv6.cpp

#include <cmath>





#include <cstdio>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include "utils.cpp"
#include <iomanip>
#include <fftw3.h>


//#include <stdlib.h>
//#include <unistd.h>
//#include <fftw.h>
//#include <string>
//#include <sstream>
//#include "init.h"
//#include <stdio.h>


using namespace std;
using std::vector;

typedef complex<double> dcomplex;
typedef vector<double>* vektor;
typedef const double cdouble;


double thdeg, nOrdPum;
int cCryst, nt;
int warned = 0;
int chirpType;
double *pLambdaj, *iLambdaj, *sLambdaj, *t;
double pLam1, pLam2, sLam1, sLam2, iLam1, iLam2;
double dtps, tlead, dw;	
double alps, alpi, alpp;
double pLambda_nm, sLambda_nm, iLambda_nm, cl;
dcomplex *timeProfSig, *timeProfIdl, *timeProfPum;
double *absTP, dk, dzcm;
double Xcm2;
fftw_plan p;
dcomplex c1(1,0),ci(0,1);


int main(int argc, char *argv[]) {










	int noParm_def, i, j, mode;
	int cStage, nStage;
	int cryst1, cryst2, cryst3;
	int sProf, pProf, iProf;
	int noStep, noStep1, noStep2, noStep3;
	int ppm, frame;
	string version;

	double *nPumj, *nSigj, *nIdlj, *phiPumj, *phiSigj, *phiIdlj, *temp;
	double xeff, gamsij;
	double thdeg1, thdeg2, thdeg3, nColl_deg, gamma_rad; //tanThetaSq;
	double crysLth, crysLth1, crysLth2, crysLth3;
	double pEJ, pEJ1, pEJ2, pEJ3, sEJ, sEJ1, iEJ, iEJ1;
	double sigTra12, sigTra23;
	double dtPumpL, dtPumpL1, dtPumpL2, dtPumpL3;
	double dtPumpT, dtPumpT1, dtPumpT2, dtPumpT3;
	double dtSigL, dtSigL1, dtIdlL, dtIdlL1; 
	double dtSigT, dtSigT1, dtIdlT, dtIdlT1;
	double Xcm2_1, Xcm2_2, Xcm2_3;
	double nOrdSig, nOrdIdl, nExPum, coePum, nonePum;
	double nOrdSigj, nOrdIdlj, coePumj; //nOrdPumj, nExPumj, nonePumj;
	double kSig, kIdl, kPum, knSig, knIdl, knPum, dk_grid;
	double cohLength, tWin;
	double kvPumj, kvSigj, kvIdlj;
	double phasVelPum, phasVelSig, phasVelIdl;
	double phasdtPum, phasdtSig, phasdtIdl;
	double gVelPum, gVelSig, gVelIdl, gtdPum, gtdSig, gtdIdl;
	double grDelPumSig, grDelPumIdl, grDelSigIdl;
	double pa, pb, pc, pd, pu, pw, puw, pAng;
	double omegaSig0, omegaPum0, omegaIdl0;
	double dphm;
	double fwp, fwi, fws;
	double tcPum, tcSig, tcIdl;
	double chpSig, chpPum, chpIdl, chpSig23;
	double chpSig2, chpPum2, chpIdl2, chpSig223;
	double chpSigL, chpPumL, chpIdlL, chpSigNL, chpPumNL, chpIdlNL;
	double phiP, phiS, phiI;
	
	vektor pars;
	dcomplex cPhMisM;
	dcomplex *cPhiPumj, *cPhiSigj, *cPhiIdlj;
	dcomplex *cSig, *cPum, *cIdl;
	
	extern double FindMax(double*,int);
	extern int spectrum(complex<double>*,double*,const char*,int);
	extern int OPA(complex<double>*,complex<double>*,complex<double>*,int,int,double,double,double,int,int,int,int,complex<double>*,complex<double>*,complex<double>*, int);
	extern int writeHeader(const char*, int, ...);
/*=====================================================================*/
// Giving values to some variables & constants


	noParm_def = 66; // Number of parameters
	
	
/*=====================================================================*/

/*=====================================================================*/
// Memory allocation
	cout << "Allocating memory...";
	t = new double[nt];
	nPumj = new double[nt];
	nSigj = new double[nt];
	nIdlj = new double[nt];
	phiPumj = new double[nt];
	phiSigj = new double[nt];
	phiIdlj = new double[nt];
	temp = new double[nt];
	cPhiPumj = new dcomplex[nt];
	cPhiSigj = new dcomplex[nt];
	cPhiIdlj = new dcomplex[nt];
	pLambdaj = new double[nt];
	sLambdaj = new double[nt];
	iLambdaj = new double[nt];
	cSig = new dcomplex[nt];
	cPum = new dcomplex[nt];
	cIdl = new dcomplex[nt];
	timeProfSig = new dcomplex[nt];
	timeProfPum = new dcomplex[nt];
	timeProfIdl = new dcomplex[nt];
	absTP = new double[nt];
	cout << "Done!" << endl;
/*=====================================================================*/




itt



/*=====================================================================*/

/*=====================================================================*/

/*=====================================================================*/

/*=====================================================================*/




/*=====================================================================*/

/*=====================================================================*/

/*=====================================================================*/
	// OPA
		if (cStage==1)
			writeHeader("output//Energies.dat",4,"step","EP[mJ]","ES[mJ]","ET[mJ]");
		OPA(timeProfPum, timeProfSig, timeProfIdl, noStep, cStage, fwp, fws, fwi, chirpType, sProf, pProf, iProf, cPhiSigj, cPhiIdlj, cPhiPumj, cStage-1);
	
	
	
	

		
				
	// Writing data to file
		//writeToFile("phase_sig.dat", sLambdaj, phiSigj);
		

	} // end of stage loop

// Cleaning up.
	delete pars;
	delete t;
//	delete kvPum;
//	delete kvSig;
//	delete kvIdl;

}
/*
 * noisy pulse
 * background noise
 * KDP
 * read in bck
	//double* ggg;
	//ggg = new double[nt];

	//for (j=0;j<nt;j++) {
	//	ggg[j] = j;
	//}
	//	writeToFile("test2.dat", ggg, (timeProfPum));
		//writeToFile("test2.dat", ggg, (timeProfSig));
		//writeToFile("test2.dat", ggg, (timeProfIdl));
//exit(0);
*/
/*diag stuff
 * 
 * 			for (j=0;j<nt;j++) {
				phiPumj[j] = abs(timeProfPum[j]);
			}
			writeToFile("test2.dat", phiPumj, phiPumj);
			exit(0);
 * 
 * 
 * 
 * 
 * 
 * 
*/

