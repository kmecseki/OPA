/***********************************************************/
/***********************************************************/
//             OPA program C++ version 1.0				   //
/***********************************************************/
/***********************************************************/
// Compilation: g++ -lgsl -lgslcblas -lfftw3 OPAv3.cpp

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
#include <fftw3.h>

using namespace std;
using std::vector;

typedef complex<double> dcomplex;
typedef vector<double>* vektor;
typedef const double cdouble;

cdouble tPi = 2*4*atan(1);
cdouble small = 1e-20;
double thdeg, nOrdPum;
int cCryst, nt;
int warned = 0;
double *pLambdaj, *iLambdaj, *sLambdaj, *t;
double pLam1, pLam2, sLam1, sLam2, iLam1, iLam2;
double dtps, tlead;
//fftw_complex *timeProfile;
dcomplex *timeProfile;
double *absTP;

int main(int argc, char *argv[]) {

// Initialising

	const char *fname;
	const int ndim = 524288;
	cdouble c = 299792458;
	cdouble eps0 = 8.85418782e-12;
	cdouble hpl = 6.626068e-34;
	cdouble hce0=0.5*c*eps0;
	int noParm_def, i, j, mode;
	int cStage, nStage;
	int cryst1, cryst2, cryst3;
	int sProf, pProf, iProf;
	int noStep, noStep1, noStep2, noStep3;
	int ppm, frame, chirpType;
	string version;

	double *nPumj, *nSigj, *nIdlj, *phiPumj, *phiSigj, *phiIdlj, *temp;
	double xeff, gamsij;
	double thdeg1, thdeg2, thdeg3, nColl_deg, gamma_rad; //tanThetaSq;
	double crysLth, crysLth1, crysLth2, crysLth3, dzcm;
	double pLambda_nm, sLambda_nm, iLambda_nm;
	double pEJ, pEJ1, pEJ2, pEJ3, sEJ, sEJ1, iEJ, iEJ1;
	double sigTra12, sigTra23;
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
	double phasVelPum, phasVelSig, phasVelIdl;
	double phasdtPum, phasdtSig, phasdtIdl;
	double gVelPum, gVelSig, gVelIdl, gtdPum, gtdSig, gtdIdl;
	double grDelPumSig, grDelPumIdl, grDelSigIdl;
	double pa, pb, pc, pd, pu, pw, puw, pAng;
	double dphm, chp1, chp2, chp3;
	double fwp, fwi, fws;
	
	vektor pars;
	dcomplex c1,ci;
	dcomplex cPhMisM;
	dcomplex *cPhiPumj, *cPhiSigj, *cPhiIdlj;
	dcomplex *cSig, *cPum, *cIdl;

// Giving values to some variables & constants

	version = "1.0";
	noParm_def = 49; // Number of parameters
	c1 = (1,0); ci = (0,1);
	
	cout << "OPA program version " << version << endl;
/*=====================================================================*/
// In absence of argument, using default input.txt input file.
	if(argc==1){
		cout << "No input argument, using default 'input.txt'" << endl;
		fname = "input.txt";
	}
	else {
		cout << "Input file given in argument -" << argv[1] << "- is used." << endl;
		fname = argv[1];
	}
/*=====================================================================*/
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
	tlead = pars->at(25);
	dtps = pars->at(26);
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
	frame = (int) pars->at(44);
	chirpType = (int) pars->at(45);
	chp1 = pars->at(46);
	chp2 = pars->at(47);
	chp3 = pars->at(48);
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
	//timeProfile = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
	timeProfile = new dcomplex[nt];
	absTP = new double[nt];
	cout << "Done!" << endl;
/*=====================================================================*/
// Cleaning up old files
	ifstream ifile("phase_sig.dat");
	if (ifile) {
		if (remove("phase_sig.dat")!=0)
			errorhl(8);
	}
/*=====================================================================*/
// Construct time array.
	//for (i=1;i<(ndim+1);i++) {
	for (i=0;i<nt;i++) {
		t[i]=i*dtps-tlead;
	}
/*=====================================================================*/
// Set up stage parameters.
	for (cStage=1;cStage<=nStage;cStage++) {
		warned = 0;
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
/*=====================================================================*/
	// Time window
		tWin = dtps*nt;
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
		nExPum = calcRefInd(cCryst, pLambda_nm, 1);
/*=====================================================================*/
// Calculate phase matching
		cout << "\tn_o Signal\t" << "n_o Idler\t" << "n_o Pump\t" << "n_e Pump" << endl;
		cout << "\t" << nOrdSig << "\t\t" << nOrdIdl << "\t\t" << nOrdPum << "\t\t" << nExPum << endl;
	// Angular wavenumbers
		kSig = nOrdSig * tPi*1e3/sLambda_nm; // in 1/micron
		kIdl = nOrdIdl * tPi*1e3/iLambda_nm;
		kPum = nOrdPum * tPi*1e3/pLambda_nm;
	// Effective angular wavenumbers
		gamma_rad = asin(kSig * sin(deg2rad(nColl_deg))/kIdl);
		coePum = CalcPumCo(pLambda_nm);
		knSig = kSig * cos(deg2rad(nColl_deg)); // in 1/micron
		knIdl = kIdl * cos(gamma_rad);
		knPum = kPum * coePum;
	// Calculate phase matching angle
		pa = cos(deg2rad(nColl_deg))*nOrdSig/sLambda_nm; // nm
		pb = cos(gamma_rad)*nOrdIdl/iLambda_nm;
		pc = nOrdPum/pLambda_nm;
		pd = nExPum/pLambda_nm;
		pu = pow(((pa+pb)/pc),2);
		pw = pow(((pa+pb)/pd),2);
		puw = (1-pu)/(pw-1);
		if (puw<0) errorhl(4);
		pAng = rad2deg(atan(sqrt(puw)));
	// Phase mismatch and coherence length
		if (ppm == 1) {
			dk = 0;
			cohLength = 0;
		}
		else {
			dk = kPum-kSig-kIdl;
			cohLength = abs(tPi/(2*dk));  // in micron
		}
		cPhMisM = exp(ci*dk*dzcm/2e4);
		// here could end the scan loop
/*=====================================================================*/
	// Wavelength scan limits
		dw = tPi*1e-3/tWin; // 1/fs
		pLam1 = tPi*c*1e-6/(omegaPum0+nt/2*dw);  // nm
		pLam2 = tPi*c*1e-6/(omegaPum0-nt/2*dw);
		sLam1 = tPi*c*1e-6/(omegaSig0+nt/2*dw);
		sLam2 = tPi*c*1e-6/(omegaSig0-nt/2*dw);		
		iLam1 = tPi*c*1e-6/(omegaIdl0+nt/2*dw);
		iLam2 = tPi*c*1e-6/(omegaIdl0-nt/2*dw);		
		//cout << setprecision(3);
		cout << "\tSignal\t" << "\t\tPump\t" << "\t\tIdler\t" << endl;
		cout << "Scan wavelengths:" << endl;
		cout << "   " << sLam1 << " " << sLam2 << "\t  " << pLam1 << " " << pLam2 << "\t  " << iLam1 << " " << iLam2 << "\t nm " << endl;
/*=====================================================================*/
	// Calculating angular wavenumbers and phases for all the wavelengths in this range.
		dk_grid = 1e9*dw/c; // in microns
		for (j=0;j<nt;j++) {
			kvPumj = tPi*1e3/pLam2 + j*dk_grid; // microns
			kvSigj = tPi*1e3/sLam2 + j*dk_grid;
			kvIdlj = tPi*1e3/iLam1 - j*dk_grid; // FIX: hopefully it fixes idler issue
			if (kvIdlj<0) kvIdlj=0.1;
			pLambdaj[j] = tPi*1e3/kvPumj; // nm
			sLambdaj[j] = tPi*1e3/kvSigj;
			iLambdaj[j] = tPi*1e3/kvIdlj; // going backwards
			nPumj[j] = CalcPumCo(pLambdaj[j]);
			nPumj[j] = nPumj[j] * calcRefInd(cCryst, pLambdaj[j], 2);
			nSigj[j] = calcRefInd(cCryst, sLambdaj[j], 3);
			nIdlj[j] = calcRefInd(cCryst, iLambdaj[j], 3); // going backwards
			phiPumj[j] = -kvPumj*dzcm*nPumj[j];
			phiSigj[j] = -kvSigj*dzcm*nSigj[j]*cos(deg2rad(nColl_deg));
			gamsij = kvSigj*nSigj[j]*sin(deg2rad(nColl_deg))/(kvIdlj*nIdlj[j]);
			if(abs(gamsij)<=1.0)
				phiIdlj[j] = -kvIdlj*dzcm*nIdlj[j]*cos(asin(gamsij));
			else phiIdlj[j] = 0;
		}
/*=====================================================================*/
/*	// FIX: Reverse again idler indeces??? Ebbe belenezni jobban
		for (j=0;j<nt;j++) {
			temp[j] = phiIdlj[nt-1-j];
		}
		for (j=0;j<nt;j++) {
			phiIdlj[j] = temp[j];
		}
		for (j=0;j<nt;j++) {
			temp[j] = nIdlj[nt-1-j];
		}
		for (j=0;j<nt;j++) {
			nIdlj[j] = temp[j];
		}*/
/*=====================================================================*/
	// Phase velocities
		phasVelPum = c*1e2/nPumj[nt/2]; // cm/sec
		phasVelSig = c*1e2/nSigj[nt/2];
		phasVelIdl = c*1e2/nIdlj[nt/2];
	// Phase delays
		phasdtPum = crysLth*1e12/phasVelPum; // ps
		phasdtSig = crysLth*1e12/phasVelSig;
		phasdtIdl = crysLth*1e12/phasVelIdl;
	// Group velocities
		gVelPum = c*1e2/(nPumj[nt/2]+kPum/nOrdPum*c*(nPumj[nt/2+1]-nPumj[nt/2-1])/(2.0*dw*1e9)); // cm/s
		gVelSig = c*1e2/(nSigj[nt/2]+kSig/nOrdSig*c*(nSigj[nt/2+1]-nSigj[nt/2-1])/(2.0*dw*1e9));
		gVelIdl = c*1e2/(nIdlj[nt/2]+kIdl/nOrdIdl*c*(nIdlj[nt/2-1]-nIdlj[nt/2+1])/(2.0*dw*1e9));
		// FIX: az eggyel fentebb levo sort atirni, ha idler indexet cserelek
	// Group time delays
		gtdPum = 1e12*crysLth/gVelPum; // ps
		gtdSig = 1e12*crysLth/gVelSig;
		gtdIdl = 1e12*crysLth/gVelIdl;
		
		cout << "Phase velocities:" << endl;
		cout << "\t" << phasVelSig << "\t     " << phasVelPum << "\t     " << phasVelIdl << "\t" << "cm/s" << endl;
		cout << "Phase delays:" << endl;
		cout << "\t" << phasdtSig << "\t\t\t" << phasdtPum << "\t\t\t" << phasdtIdl << "\t\t" << "ps" << endl;
		cout << "Group velocities:" << endl;
		cout << "\t" << gVelSig << "\t     " << gVelPum << "\t     " << gVelIdl << "\t" << "cm/s" << endl;
		cout << "Group delays:" << endl;
		cout << "\t" << gtdSig << "\t\t\t" << gtdPum << "\t\t\t" << gtdIdl << "\t\t" << "ps" << endl;
	
	// Relative group-time delay
		grDelPumSig = (gtdPum - gtdSig)* noStep; // ps
		grDelPumIdl = (gtdPum - gtdIdl)* noStep;
		grDelSigIdl = (gtdSig - gtdIdl)* noStep;
		
		cout << "Relative Group time delays" << endl;
		cout << " Signal - Pump \t\t" << -grDelPumSig  << " ps \t Positive for slow signal";
		if (-grDelPumSig<0)
			cout << " -> OK! " << endl;
		cout << " Idler - Pump \t\t" << -grDelPumIdl  << " ps \t Positive for slow idler";
		if (-grDelPumIdl<0) 
			cout << " -> OK! " << endl;
		cout << " Idler - Signal \t" << -grDelSigIdl  << " ps \t Positive for fast signal";
		if (-grDelSigIdl<0) 
			cout << " -> OK! " << endl;
	// Print phase matching info
		cout << "Phase matching summary:" << endl;
		cout << " Phase matching angle: " << pAng << endl;
		cout << " Angle of propagation: " << thdeg << endl;
		cout << " Coherence length: " << cohLength << " um" << endl;
/*=====================================================================*/
	// Calculate at which angle k vector triangle closes properly	
		if (nColl_deg!=0) {
			cout << "Non collinear geometry" << endl;
			double z, alp;
			z = (pow(knPum,2)-pow(kIdl,2)+pow(kSig,2))/(2*knPum*kSig);
			if (z<=1) {
				alp = rad2deg(acos(z));
				cout << " k vector triangle closes @ alpha =" << alp << " degs" << endl;
			}
			else {
				alp = 0;
				errorhl(5);
			}
		// Calculated beta from Ross(15)		
			if (gVelSig/gVelIdl>=1) 
				cout << "No bandwidth optimised alpha exists" << endl;
			else {
				double beta2, alpha2, alpha3, nPumOpt, theOpt;
				beta2 = acos(gVelSig/gVelIdl);
			// From Ross(16) - good approx to the optimal noncoll angle, it will
			// coincide with alpha at optimal phase matchin
				alpha2 = rad2deg(asin(kIdl/knPum)*sin(beta2));
			// Calculated from Geoff's formula, preferable to Ross(16), as it
			// does not depend on knPum, which is weakly dependent on theta
				alpha3 = atan(sin(beta2)/(kSig/kIdl+cos(beta2)));
			// Using Ross(16) backwards to get optimum pump refractive index
				nPumOpt = kIdl*sin(beta2)/(kPum/nOrdPum*sin(alpha3));
				theOpt = rad2deg(asin(sqrt((pow(nOrdPum/nPumOpt,2)-1)/(pow(nOrdPum/nExPum,2)-1))));
				cout << "Optimal alpha: " << rad2deg(alpha3) << "\t Optimal theta: " << theOpt << endl;
				cout << "Actual alpha: " << nColl_deg << "\t Actual theta: " << thdeg << endl;
			}
		}
		else cout << "Collinear geometry" << endl;
		if (mode!=1) exit(0);
/*=====================================================================*/
	// Local time frame selection
	// negative sign because of: gtd = -d(phi)/domega to ensure spatial phase of form exp(-ikz)
		if (frame==1) 
			dphm = -dw*gtdSig;
		else if (frame==2) 
			dphm = -dw*(gtdSig+gtdPum)/2;
		else if (frame==3) 
			dphm = -dw*gtdPum;
		else errorhl(6);
	// Remove central phases and phase gradients
		for (j=0;j<nt;j++) {
			phiPumj[j] = phiPumj[j] - j*dphm - phiPumj[nt/2];
			phiSigj[j] = phiSigj[j] - j*dphm - phiSigj[nt/2];
			phiIdlj[j] = phiIdlj[j] - j*dphm - phiIdlj[nt/2];
			cPhiPumj[j] = exp(ci*phiPumj[j]);
			cPhiSigj[j] = exp(ci*phiSigj[j]);
			cPhiIdlj[j] = exp(ci*phiIdlj[j]);
		}
/*=====================================================================*/
		cPum = 0;
		if (cStage==1) {
			cSig = 0;
			cIdl = 0;
		}
/*=====================================================================*/
	// Create pulses
	cout << "Creating pulses..." << endl;
	// PUMP
		if (pEJ!=0) {
			cout << "Pump:" << endl;
			if (pProf<0) errorhl(10);
			fwp = 0.5*c*eps0*(nOrdPum*coePum)*Xcm2; // Check units if sg is wrong
			GenProfile(pProf, dtPumpL, dtPumpT, Xcm2, pEJ, fwp);
		}
		if (cStage==1) {
			cout << "Signal:" << endl;
			fws = 0.5*c*eps0*(nOrdSig)*Xcm2; // Check units if sg is wrong
			GenProfile(sProf, dtSigL, dtSigT, Xcm2, sEJ, fws);
		}
		else // ide jon majd a signal ujrahasznositasa
		{}
		
		
		
		
	// Writing data to file
		writeToFile("phase_sig.dat", sLambdaj, phiSigj);
		
	
	
	
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
- mindenfele outputfileok megnyitasa

*/

