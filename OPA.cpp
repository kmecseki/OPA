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

