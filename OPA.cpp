/***********************************************************/
/***********************************************************/
//             OPA program C++ version 1.0				   //
/***********************************************************/
/***********************************************************/

#include <cstdio>
#include <cmath>
#include <complex>
//#include <stdlib.h>
#include <cstring>
//#include <unistd.h>
//#include <fftw.h>
#include <fstream>
#include <iostream>
//#include <string>
#include <vector>
//#include <sstream>
//#include "init.h"
#include "utils.cpp"
//#include <stdio.h>
#include <iomanip>

using namespace std;
using std::vector;

typedef complex<double> dcomplex;
typedef vector<double>* vektor;
typedef const double cdouble;

cdouble tPi = 2*4*atan(1);
double thdeg, nOrdPum;
int cCryst;

int main(int argc, char *argv[]) {

// Initialising

	const char *fname;
	const int ndim = 524288;
	cdouble c = 299792458;
	cdouble eps0 = 8.85418782e-12;
	cdouble hpl = 6.626068e-34;
	cdouble hce0=0.5*c*eps0;
	int noParm_def, i, j, mode;
	int cStage, nStage, nt;
	int cryst1, cryst2, cryst3;
	int sProf, pProf, iProf;
	int noStep, noStep1, noStep2, noStep3;
	int ppm;
	string version;

	double xeff;
	double thdeg1, thdeg2, thdeg3, nColl_deg, gamma_rad; //tanThetaSq;
	double crysLth, crysLth1, crysLth2, crysLth3, dzcm;
	double pLambda_nm, sLambda_nm, iLambda_nm;
	double pLam1, pLam2, sLam1, sLam2, iLam1, iLam2;
	double pEJ, pEJ1, pEJ2, pEJ3, sEJ, sEJ1, iEJ, iEJ1;
	double sigTra12, sigTra23, dtps, tlead;
	double dtPumpL, dtPumpL1, dtPumpL2, dtPumpL3;
	double dtPumpT, dtPumpT1, dtPumpT2, dtPumpT3;
	double dtSigL, dtSigL1, dtIdlL, dtIdlL1; 
	double dtSigT, dtSigT1, dtIdlT, dtIdlT1;
	double Xcm2, Xcm2_1, Xcm2_2, Xcm2_3;
	double omegaSig0, omegaPum0, omegaIdl0;
	double nOrdSig, nOrdIdl, nExPum, coePum, nonePum;
	double nOrdSigj, nOrdIdlj, coePumj; //nOrdPumj, nExPumj, nonePumj;
	double kSig, kIdl, kPum, knSig, knIdl, knPum, dk, dk_grid;
	double cohLength, dw, tWin;
	double kvPumj, kvSigj, kvIdlj;
	double pLambdaj, iLambdaj, sLambdaj;
	double nPumj;
	
	vektor pars;
	vektor t;
	dcomplex c1,ci;
	dcomplex cPhMisM;
	
// Giving values to some variables & constants

	version = "1.0";
	noParm_def = 45; // Number of parameters
	c1 = (1,0); ci = (0,1);
	
	cout << "OPA program version " << version << endl;
/***********************************************************/
// In absence of argument, using default input.txt input file.
	if(argc==1){
		cout << "No input argument, using default 'input.txt'" << endl;
		fname = "input.txt";
	}
	else {
		cout << "Input file given in argument -" << argv[1] << "- is used." << endl;
		fname = argv[1];
	}
/***********************************************************/	
// Reading parameters from input file.
	pars = new vector<double>();
	openfile(fname, pars);
	cout << "Read/Required parameters from file: " << pars->size() << "/" << noParm_def << endl;
	if(pars->size()!=noParm_def) errorhl(2);
	mode = (int) pars->at(0);
	nStage = (int) pars->at(1);
	cryst1 = (int) pars->at(2);
	cryst2 = (int) pars->at(3);
	cryst3 = (int) pars->at(4);
	crysLth1 = pars->at(5);
	crysLth2 = pars->at(6);
	crysLth3 = pars->at(7);
	thdeg1 = pars->at(8);
	thdeg2 = pars->at(9);
	thdeg3 = pars->at(10);
	nColl_deg = pars->at(11);
	ppm = pars->at(12);
	pLambda_nm = pars->at(13);
	sLambda_nm = pars->at(14);
	pEJ1 = pars->at(15);
	pEJ2 = pars->at(16);
	pEJ3 = pars->at(17);
	sEJ1 = pars->at(18);
	sigTra12 = pars->at(19);
	sigTra23 = pars->at(20);
	iEJ1 = pars->at(21);
	sProf = (int) pars->at(22);
	pProf = (int) pars->at(23);
	iProf = (int) pars->at(24);
	dtps = pars->at(25);
	tlead = pars->at(26);
	nt = (int) pars->at(27);
	noStep1 = (int) pars->at(28);
	noStep2 = (int) pars->at(29);
	noStep3 = (int) pars->at(30);
	dtPumpL1 = pars->at(31);
	dtPumpL2 = pars->at(32);
	dtPumpL3 = pars->at(33);
	dtPumpT1 = pars->at(34);
	dtPumpT2 = pars->at(35);
	dtPumpT3 = pars->at(36);
	dtSigL1 = pars->at(37);
	dtIdlL1 = pars->at(38);
	dtSigT1 = pars->at(39);
	dtIdlT1 = pars->at(40);
	Xcm2_1 = pars->at(41);
	Xcm2_2 = pars->at(42);
	Xcm2_3 = pars->at(43);
	tWin = pars->at(44);
/***********************************************************/
// Memory allocation
	cout << "Allocating memory...";
	t = new vector<double>();
//	kvPum = new double[nt];
//	kvSig = new double[nt];
//	kvIdl = new double[nt];
	cout << "Done!" << endl;
/***********************************************************/	
// Construct time array.
	for (i=1;i<(ndim+1);i++) {
		t->push_back(i*dtps-tlead);
	}
/***********************************************************/
// Set up stage parameters.
	for (cStage=1;cStage<=nStage;cStage++) {
		if (cStage==1) {
			cCryst = cryst1;
			thdeg = thdeg1;
			xeff = xeffc(cCryst,cStage);
			noStep = noStep1;
			crysLth = crysLth1;
			dtPumpL = dtPumpL1;
			dtPumpT = dtPumpT1;
			dtSigL = dtSigL1;
			dtSigT = dtSigT1;
			dtIdlL = dtIdlL1;
			dtIdlT = dtIdlT1;
			pEJ = pEJ1;
			sEJ = sEJ1;
			iEJ = iEJ1;
			Xcm2 = Xcm2_1;
		}
		else if(cStage==2) {
			cCryst = cryst2;
			thdeg=thdeg2;
			xeff=xeffc(cCryst,cStage);
			noStep = noStep2;
			crysLth = crysLth2;
			dtPumpL = dtPumpL2;
			dtPumpT = dtPumpT2;			
			pEJ = pEJ2;
			Xcm2 = Xcm2_2;
		}
		else if(cStage==3) {
			cCryst = cryst3;
			thdeg = thdeg3;
			xeff = xeffc(cCryst,cStage);
			noStep = noStep3;
			crysLth = crysLth3;
			dtPumpL = dtPumpL3;
			dtPumpT = dtPumpT3;			
			pEJ = pEJ3;
			Xcm2 = Xcm2_3;
		}
/***********************************************************/		
	// Distance step
		dzcm = crysLth/((double) noStep);
	// Calculate idler wavelength
		omegaSig0 = tPi*c/(1e6*sLambda_nm); // in 1/fs
		omegaPum0 = tPi*c/(1e6*pLambda_nm);
		omegaIdl0 = omegaPum0 - omegaSig0;
		iLambda_nm = tPi*c/(1e6*omegaIdl0);
		cout << "Idler wavelength: " << iLambda_nm << " nm." << endl;
	// Calculate refractive indices
		nOrdSig = calcRefInd(cCryst, sLambda_nm, 3);
		nOrdIdl = calcRefInd(cCryst, iLambda_nm, 3);
		nOrdPum = calcRefInd(cCryst, pLambda_nm, 2);
		//nExPum = calcRefInd(cCryst, pLambda_nm, 1);	
/***********************************************************/
// Calculate phase matching
		cout << "\tn_o Signal\t" << "n_o Idler\t" << "n_o Pump\t" << "n_e Pump" << endl;
		cout << "\t" << nOrdSig << "\t\t" << nOrdIdl << "\t\t" << nOrdPum << "\t\t" << nExPum << endl;
	// Angular wavenumbers
		kSig = nOrdSig * tPi*1e3/sLambda_nm; // in 1/micron
		kIdl = nOrdIdl * tPi*1e3/iLambda_nm;
		kPum = nOrdPum * tPi*1e3/pLambda_nm;
	// Effective angular wavenumbers
		gamma_rad = asin(kSig * sin(deg2rad(nColl_deg)/kIdl));
		//tanThetaSq = tan(deg2rad(thdeg)) * tan(deg2rad(thdeg));
		//nonePum = (nOrdPum/nExPum) * (nOrdPum/nExPum);
		//coePum = sqrt((1 + tanThetaSq)/(1 + nonePum*tanThetaSq));
		coePum = CalcPumCo(pLambda_nm);
		knSig = kSig * cos(deg2rad(nColl_deg)); // in 1/micron
		knIdl = kIdl * cos(gamma_rad);
		knPum = kPum * coePum;
	// Phase mismatch and coherence length
		if (ppm == 0) {
			dk = 0;
			cohLength = 0;
		}
		else {
			dk = kPum-kSig-kIdl;
			cohLength = abs(tPi/(2*dk));  // in micron
		}
		cPhMisM = exp(ci*dk*dzcm/2e4);
		// here could end the scan loop
/***********************************************************/
	// Wavelength scan limits
		dw = tPi*1e-3/tWin; // 1/fs
		pLam1 = tPi*c*1e-6/(omegaPum0+nt/2*dw);  // nm
		pLam2 = tPi*c*1e-6/(omegaPum0-nt/2*dw);
		sLam1 = tPi*c*1e-6/(omegaSig0+nt/2*dw);
		sLam2 = tPi*c*1e-6/(omegaSig0-nt/2*dw);		
		iLam1 = tPi*c*1e-6/(omegaIdl0+nt/2*dw);
		iLam2 = tPi*c*1e-6/(omegaIdl0-nt/2*dw);		
		cout << "Scan wavelengths:" << endl;
		//cout << setprecision(3);
		cout << "\tSignal\t" << "\t\tPump\t" << "\t\tIdler\t" << endl;
		cout << "   " << sLam1 << " " << sLam2 << "\t  " << pLam1 << " " << pLam2 << "\t  " << iLam1 << " " << iLam2 << endl;
/***********************************************************/
	// Calculating angular wavenumbers and phases for all the wavelengths in this range.
		dk_grid = 1e9*dw/c; // in microns
		for (j=0;j<nt;j++) {
			kvPumj = tPi*1e3/pLam2 + j*dw; //microns
			kvSigj = tPi*1e3/sLam2 + j*dw;
			kvIdlj = tPi*1e3/iLam2 + j*dw;
			pLambdaj = tPi*1e3/kvPumj; // nm
			sLambdaj = tPi*1e3/kvSigj;
			iLambdaj = tPi*1e3/kvIdlj;
			nPumj = CalcPumCo(pLambdaj);
			nPumj = nPumj * nOrdPum;

			
		
		}
		
	  

	
	
	} // end of stage loop

// Cleaning up.
	delete pars;
	delete t;
//	delete kvPum;
//	delete kvSig;
//	delete kvIdl;

}
/*
- uj pump vagy regi ujrahasznalasa
- temporal widths
- chirp modszerek, mekkora legyen hol stb
- idoablak, grid, hany db
- signal es idler time delayek
- beam diameters
- lepesszam
- kezdofazisok
- koherenciaidok
- hatter? - kiszedheto nem?
- perf phase matching
- mindenfele outputfileok megnyitasa

*/

