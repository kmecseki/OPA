/***********************************************************/
/***********************************************************/
//             OPA program C++ version 1.0				   //
/***********************************************************/
/***********************************************************/
// Compilation: g++ -lgsl -lgslcblas -lfftw3 OPAv6.cpp

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
cdouble c = 2.997925000000e8;//299792458;

int main(int argc, char *argv[]) {

// Initialising

	const char *fname;
	const int ndim = 524288;
	cdouble eps0 = 8.854e-12;//8.85418782e-12;
	cdouble hpl = 6.6256e-34;//6.626068e-34;
	cdouble hce0=0.5*c*eps0;
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

	version = "1.0";
	noParm_def = 66; // Number of parameters
	
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
	int o=0;
	cout << "Read/Required parameters from file: " << pars->size() << "/" << noParm_def << endl;
	if(pars->size()!=noParm_def) errorhl(2);
	mode = (int) pars->at(o++);
	nStage = (int) pars->at(o++);
	cryst1 = (int) pars->at(o++);
	cryst2 = (int) pars->at(o++);
	cryst3 = (int) pars->at(o++);
	crysLth1 = pars->at(o++);
	crysLth2 = pars->at(o++);
	crysLth3 = pars->at(o++);
	thdeg1 = pars->at(o++);
	thdeg2 = pars->at(o++);
	thdeg3 = pars->at(o++);
	nColl_deg = pars->at(o++);
	ppm = pars->at(o++);
	pLambda_nm = pars->at(o++);
	sLambda_nm = pars->at(o++);
	pEJ1 = pars->at(o++);
	pEJ2 = pars->at(o++);
	pEJ3 = pars->at(o++);
	sEJ1 = pars->at(o++);
	sigTra12 = pars->at(o++);
	sigTra23 = pars->at(o++);
	iEJ1 = pars->at(o++);
	sProf = (int) pars->at(o++);
	pProf = (int) pars->at(o++);
	iProf = (int) pars->at(o++);
	tlead = pars->at(o++);
	dtps = pars->at(o++);
	nt = (int) pars->at(o++);
	noStep1 = (int) pars->at(o++);
	noStep2 = (int) pars->at(o++);
	noStep3 = (int) pars->at(o++);
	dtPumpL1 = pars->at(o++);
	dtPumpL2 = pars->at(o++);
	dtPumpL3 = pars->at(o++);
	dtPumpT1 = pars->at(o++);
	dtPumpT2 = pars->at(o++);
	dtPumpT3 = pars->at(o++);
	dtSigL1 = pars->at(o++);
	dtIdlL1 = pars->at(o++);
	dtSigT1 = pars->at(o++);
	dtIdlT1 = pars->at(o++);
	Xcm2_1 = pars->at(o++);
	Xcm2_2 = pars->at(o++);
	Xcm2_3 = pars->at(o++);
	frame = (int) pars->at(o++);
	chirpType = (int) pars->at(o++);
	chpSig = pars->at(o++);
	chpPum = pars->at(o++);
	chpIdl = pars->at(o++);
	chpSig23 = pars->at(o++);
	chpSig2 = pars->at(o++);
	chpPum2 = pars->at(o++);
	chpIdl2 = pars->at(o++);
	chpSig223 = pars->at(o++);
	chpSigL = pars->at(o++);
	chpPumL = pars->at(o++);
	chpIdlL = pars->at(o++);
	chpSigNL = pars->at(o++);
	chpPumNL = pars->at(o++);
	chpIdlNL = pars->at(o++);
	tcPum = pars->at(o++);
	tcSig = pars->at(o++);
	tcIdl = pars->at(o++);
	phiP = pars->at(o++);
	phiS = pars->at(o++);
	phiI = pars->at(o++);
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
		/*double freqp, freqs, freqd; // needed for same result
		freqp = c*1e6/(1e-3*pLambda_nm);
		freqs = c*1e6/(1e-3*sLambda_nm);
		freqd = freqp-freqs;
		omegaSig0 = tPi*freqs/1e12; // 1/ps
		omegaPum0 = tPi*freqp/1e12; // 1/ps
		omegaIdl0 = tPi*freqd/1e12; // 1/ps
		omegaSig0 = omegaSig0*1e3; // 1/fs
		omegaPum0 = omegaPum0*1e3; // 1/fs
		omegaIdl0 = omegaIdl0*1e3; // 1/fs
		//had to do it like that do prevent num errors.. :/
		iLambda_nm = (c*1e6)/(freqd)*1e3;*/
		/*
		float a=77617.0;
		float b=33096.0;
		float f=333.75*pow(b,6)+pow(a,2)*(11*pow(a,2)*pow(b,2)-pow(b,6)-121*pow(b,4)-2)+5.5*pow(b,8)+a/(2*b);
		cout << f << endl;
		exit(0);
		*/
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
			//dk = kPum-kSig-kIdl;
			dk = (knPum-knSig-knIdl)*1e4; // in 1/cm
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
/**/// FIX: Reverse again idler indeces??? Ebbe belenezni jobban
	// Uncomment this section for original idler handling
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
		}/**/
/*=====================================================================*/
	// Phase velocities
		phasVelPum = c*1e2/nPumj[nt/2]; // cm/sec
		phasVelSig = c*1e2/nSigj[nt/2];
		phasVelIdl = c*1e2/nIdlj[nt/2-1]; // FIX: Uncomment for ori handling
		// phasVelIdl = c*1e2/nIdlj[nt/2]; // FIX: When it is the new way
	// Phase delays
		phasdtPum = crysLth*1e12/phasVelPum; // ps
		phasdtSig = crysLth*1e12/phasVelSig;
		phasdtIdl = crysLth*1e12/phasVelIdl;
	// Group velocities
		gVelPum = c*1e2/(nPumj[nt/2]+kPum/nOrdPum*c*(nPumj[nt/2+1]-nPumj[nt/2-1])/(2.0*dw*1e9)); // cm/s
		gVelSig = c*1e2/(nSigj[nt/2]+kSig/nOrdSig*c*(nSigj[nt/2+1]-nSigj[nt/2-1])/(2.0*dw*1e9));
		gVelIdl = c*1e2/(nIdlj[nt/2]+kIdl/nOrdIdl*c*(nIdlj[nt/2+1]-nIdlj[nt/2-1])/(2.0*dw*1e9)); // Original idler handling is in use
		//gVelIdl = c*1e2/(nIdlj[nt/2]+kIdl/nOrdIdl*c*(nIdlj[nt/2-1]-nIdlj[nt/2+1])/(2.0*dw*1e9)); // If new idler handling is in use
		// FIX: az eggyel fentebb levo sort atirni, ha idler indexet cserelek
	// Group time delays
		gtdPum = 1e12*dzcm/gVelPum; // ps
		gtdSig = 1e12*dzcm/gVelSig;
		gtdIdl = 1e12*dzcm/gVelIdl;

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
			// coincide with alpha at optimal phase matching
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
		double phiPumjh = phiPumj[nt/2];
		double phiSigjh = phiSigj[nt/2];
		double phiIdljh = phiIdlj[nt/2];
	// Remove central phases and phase gradients
		for (j=0;j<nt;j++) {
			phiPumj[j] = phiPumj[j]*1e4 - (j-nt/2)*dphm*1e3 - phiPumjh*1e4;
			phiSigj[j] = phiSigj[j]*1e4 - (j-nt/2)*dphm*1e3 - phiSigjh*1e4;
			phiIdlj[j] = phiIdlj[j]*1e4 - (j-nt/2)*dphm*1e3 - phiIdljh*1e4;
			cPhiPumj[j] = polar(1.0,phiPumj[j]);
			cPhiSigj[j] = polar(1.0,phiSigj[j]);
			cPhiIdlj[j] = exp(ci*phiIdlj[j]); //does the same, as it should tested
			}
/*=====================================================================*/
		cPum = 0;
		if (cStage==1) {
			cSig = 0;
			cIdl = 0;
		}
/*=====================================================================*/
	// Create pulses
		cout << endl << endl << "*********************************************" << endl;
		cout << "\tCreating pulses" << endl;
		cout << "*********************************************" << endl;
	// PUMP
		if (pEJ!=0) {
			cout << endl << "Pump:" << endl;
			if (pProf<0) errorhl(10);
			fwp = 0.5*c*eps0*(nOrdPum*coePum)*Xcm2; // Check units if sg is wrong
			GenProfile(timeProfPum, pProf, dtPumpL, dtPumpT, Xcm2, pEJ, fwp, tcPum, chpPum, chpPum2, chpPumL, chpPumNL);
		}
		else {
			for (j=0;j<nt;j++) {
				timeProfPum[j] = 0;
			}
		}
		//cout << pow(abs(timeProfPum[nt/2]),2) << " " << fwp << " " << dtps << endl;
		//exit(0);
	// SIGNAL
		if (cStage==1) {
			cout << endl << "Signal:" << endl;
			fws = 0.5*c*eps0*(nOrdSig)*Xcm2; // Check units if sg is wrong
			GenProfile(timeProfSig, sProf, dtSigL, dtSigT, Xcm2, sEJ, fws, tcSig, chpSig, chpSig2, chpSigL, chpSigNL);
		}
//		cout << timeProfPum[nt/2];
//		double* ji = new double[nt];
//		for (j=0;j<nt;j++) {
//			phiPumj[j] = real(timeProfPum[j]);
//			ji[j] = j;
//		}
//		writeToFile("test2.dat", ji, phiPumj);
//		exit(0);
		else if (cStage==2) {
			fws = 0.5*c*eps0*(nOrdSig)*Xcm2;
			for (j=0;j<nt;j++) {
				timeProfSig[j] = timeProfSig[j]*sqrt(Xcm2_1/Xcm2_2);
			}
		}
/*		else if (cStage==3) {
			xarea=xcm23s/xcm22s;
			tts=tts23;
			wsin=tts*wsout;
			ejs=tts*ejsout;
			phots=ejs/hnus;
			ejscm2=ejs/xcm2s;
			cs=cs*sqrt(tts/xarea);
			zjj=-dm2s*(domega**2)/(two*e6);
			chirp(cs,es,nn,imax,zjj,2);
			if(tsdy23.ne.zero)
				tshift(cs,es,nn,tsdy23,irot);
		}
*/
	// IDLER
		if (cStage==1 && iEJ!=0) {
			cout << endl << "Idler:" << endl;
			fwi = 0.5*c*eps0*(nOrdIdl)*Xcm2; // Check units if sg is wrong
			GenProfile(timeProfIdl, iProf, dtIdlL, dtIdlT, Xcm2, iEJ, fwi, tcIdl, chpIdl, chpIdl2, chpIdlL, chpIdlNL);
		}
		else {
			for (j=0;j<nt;j++) {
				timeProfIdl[j] = 0;
				fwi = 0.5*c*eps0*(nOrdIdl)*Xcm2; // Check units if sg is wrong
			}
		}
/*=====================================================================*/
	// Gain predictions - no pump depletion - I.N. Ross 1997
	// Estimation using the analyticial solution of slowly varying envelope approximation
	// Also assuming flat top spatial and temporal profiles + no pump depletion
		if (pEJ!=0) {
			cout << "Gain predictions" << endl;
			double eFieldMaxPum;
			double q1, q2, p1, p2;
			double kvSig, kvIdl, kvPum;
			double dkdzh, sina, cosa, bdz, bl;
			double fdkb, phpk, sxpk, *temp;
			complex<double> csx;
			temp = new double[nt];
			kvSig = tPi*1e7/sLambda_nm; // 1/cm
			kvPum = tPi*1e7/pLambda_nm; // 1/cm
			kvIdl = tPi*1e7/iLambda_nm; // 1/cm
			for (j=0;j<nt;j++) {
				temp[j] = abs(timeProfPum[j]);
			}
			eFieldMaxPum = FindMax(temp, nt);
			alps = -0.5*kvSig*xeff*dzcm/(nOrdSig*1e10);
			alpi = -0.5*kvIdl*xeff*dzcm/(nOrdIdl*1e10);
			alpp = -0.5*kvPum*xeff*dzcm/(nOrdPum*CalcPumCo(pLambda_nm)*1e10);
			dkdzh = dk*dzcm*0.5;
			q1 = pow((sqrt(alps*alpi)*eFieldMaxPum),2);
			q2 = dkdzh*dkdzh;
			sina = sin(dk*crysLth*0.5);
			cosa = cos(dk*crysLth*0.5);
			bdz = sqrt(abs(q1-q2));
			bl = bdz*noStep;
			fdkb = dkdzh/bdz;
			if(q1<q2) {
				cout << "Hyperbolic case" << endl;
				double chbl, shbl;
				chbl = cosh(bl);
				shbl = sinh(bl);
				p1 = (bl*sina*chbl-(dk*crysLth*0.5)*cosa*shbl);
				p2 = (bl*cosa*chbl+(dk*crysLth*0.5)*sina*shbl);
				phpk = atan(p1/p2);
				real(csx) = chbl;
				imag(csx) = fdkb*shbl;
				sxpk = pow(abs(csx),2);
			}
			else {
				cout << "Trigonometric case" << endl;
				double cosbl, sinbl;
				cosbl = cos(bl);
				sinbl = sin(bl);
				p1 = (bl*sina*cosbl-(dk*crysLth*0.5)*cosa*sinbl);
				p2 = (bl*cosa*cosbl+(dk*crysLth*0.5)*sina*sinbl);
				phpk = atan(p1/p2);
				real(csx) = cosbl;
				imag(csx) = fdkb*sinbl;
				sxpk = pow(abs(csx),2);
			}
			// sxpk lehet a gain es phpk a fazis
		}
		else
			double sxpk = 0;
/*=====================================================================*/
	// Applying initial phase
		if (cStage==1) {
			for (j=0;j<nt;j++) {
				timeProfPum[j] = polar(abs(timeProfPum[j]),arg(timeProfPum[j])+tPi/2*phiP);
				timeProfSig[j] = polar(abs(timeProfSig[j]),arg(timeProfSig[j])+tPi/2*phiS);
				timeProfIdl[j] = polar(abs(timeProfIdl[j]),arg(timeProfIdl[j])+tPi/2*phiI);
			}
		}
	// Calculate spectrums
		cout << "Calculating Pump spectrum\t";
		spectrum(timeProfPum, pLambdaj, "output//Spec_pum.dat", pProf);
		cout << "Calculating Signal spectrum\t";
		spectrum(timeProfSig, sLambdaj, "output//Spec_sig.dat", sProf);
		cout << "Calculating Idler spectrum\t";
		spectrum(timeProfIdl, iLambdaj, "output//Spec_idl.dat", iProf);
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

