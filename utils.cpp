#include <numbers>
#include <vector>
#include <complex>
#include "utils.h"

double deg2rad(double deg) {
// Returns input deg value in radians.
	return deg / 360 * 2 * std::numbers::pi;
}

double rad2deg(double rad) {
// Returns input deg value in radians.
	return rad * 360 / (2 * std::numbers::pi) ;
}

void cvector_to_fftw(int nt, std::vector<std::complex<double>> in, fftw_complex* out) {
	// Convert complex vector to fftw vector
	for (int j=0; j<nt; j++) {
		out[j][0] = in[j].real();
		out[j][1] = in[j].imag();
	}
}

void fftw_to_cvector(int nt, fftw_complex* in, std::vector<std::complex<double>> out) {
	// Convert fftw vector to complex vector
	for (int j=0; j<nt; j++) {
		out[j] = std::complex(in[j][0],in[j][1]);
	}
}


/*













#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include <cstring>
//#include <unistd.h>
//#include <fftw.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <algorithm>
#include <complex>
#include <fftw3.h>

//#include "init.h"
using std::vector;
using namespace std;
int j;
extern int cCryst;
extern int warned;
extern int nt;
extern int chirpType;
extern const double tPi, small, c;
extern double thdeg, Xcm2;
extern double *pLambdaj, *iLambdaj, *sLambdaj, *t;
extern double *absTP;
extern double pLam1, pLam2, sLam1, sLam2, iLam1, iLam2;
extern double pLambda_nm, sLambda_nm, iLambda_nm;
extern double dtps, tlead, dw;
extern fftw_plan p;
extern complex<double> ci, c0, c1;
extern double alps, alpi, alpp, dk, dzcm;
ifstream ifile("fftplan.dat");

// Function declarations
int errorhl(int);
int warninghl(int);
int openfile(const char,vector<double>*);
int writeToFile(const char*,double*,double*);
int writeToFile(const char*,double,double);
double xeffc(int,int);
double calcRefInd(int,double,int);
double deg2rad(double);
double rad2deg(double);
double CalcPumCo(double);
int GenProfile(complex<double>*,int, double,double,double,double,double,double,double,double,double,double);
double FindMax(double*,int);
double FindMax(complex<double>*,int); 
double rInt(double*);
double cInt(complex<double>*,double);
complex<double>* noise(complex<double>*,int,int,double);
int chirper_norm(complex<double>*,int,int,double,double,fftw_plan);
int fftshift(complex<double>*,int);
int chirper_direct(complex<double>*,int,double,double);
int spectrum(complex<double>*,double*,const char*,int);
double get_FWHM(complex<double>*,double*);
int OPA(complex<double>*,complex<double>*,complex<double>*,int,int,double,double,double,int,int,int,int,complex<double>*,complex<double>*,complex<double>*);
int disperse(complex<double>*,double,int,int);
// not all is here, update list when have time

int errorhl(int h) {
// This handles error codes.
	switch(h) {
		case 1: cout << "Error opening input file!" << endl; break;
		case 2: cout << "Error: wrong number of input parameters!" << endl; break;
		case 3: cout << "Error: KDP not fully supported yet!" << endl; break;
		case 4: cout << "Error: No phase matching angle!" << endl; break;
		case 5: cout << "Error: k vector triangle will not close for this theta!" << endl; break;
		case 6: cout << "Error: iframe value can be only 1,2 or 3" << endl; break;
		case 7: cout << "Error opening output file!" << endl; break;
		case 8: cout << "Error deleting old files!" << endl; break;
		case 9: cout << "Error opening signal profile data file" << endl; break;
		case 10: cout << "Error: Profile reading from file can only be used for signal pulse" << endl; break; 
		case 11: cout << "Error: Noise not tested!" << endl; break;
	}
	exit(0);
}




int writeToFile(const char *ofname, double *data1, double *data2) {
// Opens a file with name ofname, and outputs real data
	ofstream filestr;
	filestr.open(ofname, ios::out | ios::app);
	filestr.precision(15);
	if(filestr.is_open()) {
		for (j=0;j<nt;j++) {
			filestr << data1[j] << "\t" << data2[j] << endl;
		}
	}
	else {
		errorhl(7);
	}
	filestr.close();
}
/*=====================================================================*/
/*
int writeToFile(const char *ofname, double data1, double data2) {
// Opens a file with name ofname, and outputs single real data
	ofstream filestr;
	filestr.open(ofname, ios::out | ios::app);
	if(filestr.is_open()) {
		filestr << data1 << "\t" << data2 << endl;
	}
	else {
		errorhl(7);
	}
	filestr.close();
}

int writeToFile(const char *ofname, double data1, complex<double> data2) {
// Opens a file with name ofname, and outputs single complex data
	ofstream filestr;
	filestr.open(ofname, ios::out | ios::app);
	if(filestr.is_open()) {
		filestr << data1 << "\t" << real(data2) << "\t" << imag(data2) << endl;
	}
	else {
		errorhl(7);
	}
	filestr.close();
}

int writeToFile(const char *ofname, double *data1, complex<double> *data2) {
// Opens a file with name ofname, and outputs complex data
	ofstream filestr;
	filestr.open(ofname, ios::out | ios::app);
	filestr.precision(15);
	if(filestr.is_open()) {
		for (j=0;j<nt;j++) {
			filestr << data1[j] << "\t" << real(data2[j]) << "\t" << imag(data2[j]) << endl;
		}
	}
	else {
		errorhl(7);
	}
	filestr.close();
}










/*=====================================================================

/*=====================================================================
double FindMax(double *vek, int s) {
// Finds maximum of vector 'vek' of size 's'.
	double maxi=0;
	for (j=0;j<s;j++) {
		if(vek[j]>maxi)
			maxi=vek[j];
	}
	return maxi;
}
/*=====================================================================
double FindMax(complex<double> *cfield, int s) {
// Finds maximum of vector 'vek' of size 's'.
	double maxi=0;
	for (j=0;j<s;j++) {
		if(abs(cfield[j])>maxi)
			maxi=abs(cfield[j]);
	}
	return maxi;
}
/*=====================================================================

/*=====================================================================*
double cInt(complex<double> *cfield, double fww) {
// Calculates energy within an intensity profile - given complex field
	double area=0;
	for (j=0;j<nt;j++) {
		area=area+pow(abs(cfield[j]),2);
	}
	area=area*dtps*1e-12*fww;
//	cout << "area" << area << endl;
	return area;
}
/*=====================================================================*
complex<double>* noise (complex<double> *cva, int profType, int s, double fx) {
	complex<double> phas;
	double aj, sumsq=0, factor;
	for (j=0;j<s;j++) {
		real(phas)=0;
		imag(phas)=tPi*(rand()%10)/10;
		if (profType)
		// If Gaussian
			aj=exp(-fx*j*j);
		else
			aj=1/sqrt(1+fx*j*j*8*log(2));
		sumsq=sumsq+aj*aj;
		cva[j]=phas*aj;
	}
	factor=sqrt(1/sumsq);
	 // Missing bit:    call fftnr(eva,n2,factor,+1)
	return cva;
}
/*=====================================================================*
int chirper_norm(complex<double> *timeProfile, int profType, int size, double phi2, double phi3, fftw_plan p) {
// Chirping function - normal chirp
	complex<double>* spec;
	spec = new complex<double>[nt];
// FFT to get to spectral domain:
	ifstream ifile2("fftplan2.dat");
	if(!ifile2) {
		cout << "Planning fastest FFT algorithm for this CPU... This will take a while, but needs to run only once on this machine." << endl;
		p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(timeProfile),reinterpret_cast<fftw_complex*>(spec),FFTW_BACKWARD,FFTW_PATIENT);
		cout << "Plan created, saved as 'fftplan2.dat'" << endl;
		if(!fftw_export_wisdom_to_filename("fftplan2.dat"))
			cout << "Error writing plan" << endl;
	}
	else {
		cout << "Reading FFT plan from file." << endl; 
		fftw_import_wisdom_from_filename("fftplan2.dat");
		p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(timeProfile),reinterpret_cast<fftw_complex*>(spec),FFTW_BACKWARD,FFTW_PATIENT);
	}
	fftw_execute(p);
	if (profType!=0) fftshift(spec,nt);
// Here we have the spectrum, let's chirp.
	for (j=0;j<nt;j++) {
		spec[j] = polar(abs(spec[j]), phi2/2*pow((-nt/2+j)*dw,2)+phi3/6*pow((-nt/2+j)*dw,3));
	}
// Converting back to time
	p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(spec),reinterpret_cast<fftw_complex*>(timeProfile),FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for (j=0;j<nt;j++) {
	// sqrt(nt)*sqrt(nt) division required again after fft
		real(timeProfile[j]) = real(timeProfile[j])/nt; 
		imag(timeProfile[j]) = imag(timeProfile[j])/nt; 
		//timeProfile[j] = polar(abs(timeProfile[j])/nt, arg(timeProfile[j])); 
	}
	fftshift(timeProfile,nt);
}
/*=====================================================================*
int chirper_direct(complex<double> *timeProfile, int size, double chp1, double chp2) {
// This function applies linear (and nonlinear) chirp directly to the complex field
// in the time domain
// Maybe fftshift before and after is needed, check..
	cout << "This type of chirping has not been tested yet.. use with caution" << endl;
	for (j=0;j<nt;j++) {
		double par1 = tPi*pow(dtps,2)*chp1/2000;
		double par2 = tPi*pow(dtps,3)*chp2/6000;
		timeProfile[j] = polar(abs(timeProfile[j]), j*j*(par1+par2*j*j));
	}
}
/*=====================================================================*
int fftshift(complex<double> *vekt, int size) {
// FFT-shift
	complex<double> temp, temp2;
	for(j=0;j<size/4;j++){
		temp = vekt[j];
		vekt[j] = vekt[size/2-(1+j)];
		vekt[size/2-(1+j)] = temp;
		temp2 = vekt[j+size/2];
		vekt[j+size/2] = vekt[size-(1+j)];
		vekt[size-(1+j)] = temp2;
	}

}
/*=====================================================================*
int spectrum(complex<double> *timeProf, double* wl, const char *ofname, int profType) {
// This function calculates the pulse spectrum, FWHM, and writes into file
	complex<double> *spek;
	spek = new complex<double>[nt];
	double *tempspek;
	tempspek = new double[nt];
	double FWHM;
	fftw_plan p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(timeProf),reinterpret_cast<fftw_complex*>(spek),FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	if (profType!=0)
		fftshift(spek,nt);
	for (j=0;j<nt;j++) {
		spek[j] = polar(abs(spek[j])/sqrt(nt), arg(spek[j]));
		tempspek[j] = abs(spek[j]);
    }
    FWHM = get_FWHM(spek, wl); // This can be written out if needed. - not tested
    cout << "FWHM : " << FWHM << " nm" << endl;
    writeToFile(ofname, wl, tempspek);    
}
/*=====================================================================*
double get_FWHM(complex<double> *profile, double* wl) {
// Finds FWHM in indeces, need to scale
	double max = FindMax(profile, nt);
	int k1, k2;
	for (j=0;j<nt;j++) {
		if (abs(profile[j])<max/2 && abs(profile[j+1])>=max/2) {
			k1 = j;
			break;
		}
		else k1 = 0;
	}
	for (j=nt;j>0;j--) {
		if (abs(profile[j])<max/2 && abs(profile[j-1])>=max/2) {
			k2 = j;
			break;
		}
		else k2 = 0;
	}
	return abs(wl[k1]-wl[k2]);

/*=====================================================================*
int disperse(complex<double> *profile, complex<double> *phase, int nt, int proftype) {
// Applies dispersive phase
//fix: proftype not needed - remove
	complex<double>* spek;
	spek = new complex<double>[nt];
	fftshift(profile,nt);
	fftw_plan p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(profile),reinterpret_cast<fftw_complex*>(spek),FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for (j=0;j<nt;j++) {
		real(spek[j]) = real(spek[j])/sqrt(nt); 
		imag(spek[j]) = imag(spek[j])/sqrt(nt); 
	}
	fftshift(phase,nt);
	for (j=0;j<nt;j++) {
		spek[j] = spek[j]*phase[j]; 
	}
	fftshift(phase,nt);
	fftw_plan p2 = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(spek),reinterpret_cast<fftw_complex*>(profile),FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p2);
	for (j=0;j<nt;j++) {
		real(profile[j]) = real(profile[j])/sqrt(nt); 
		imag(profile[j]) = imag(profile[j])/sqrt(nt); 
	}
	fftshift(profile,nt);
}
/*=====================================================================*
double nlindx(int krnum) {
// Returns non linear index of media in cm2/W
	double n2;
	switch(krnum) {
		case 1: // BBO 
			n2 = 2.9e-16;
			break;
		case 2: // LBO
			n2 = 1.1e-16;
			break;
		case 3: // KDP
			n2 = 0;
			errorhl(3);
			break;
	}
	return n2;
}
/*=====================================================================*
int nlshift(int krnum, complex<double>* profile, double fw, double n2) {
// Nonlinear phase shift
	double phaze;
	double dphmax = 0, wcm2mx = 0;
	complex<double> pfj;
	double wcm2, dphas;
	phaze = tPi*n2*dzcm*1.e7/pLambda_nm;
	for (j=0;j<nt;j++) {
		pfj = profile[j];
		wcm2 = fw/Xcm2*pow((abs(pfj)),2);
		dphas = phaze*wcm2;
		if(dphas>dphmax)
			dphmax = dphas;
		if(wcm2>wcm2mx)
			wcm2mx = wcm2;
		profile[j] = pfj*polar(1.0,-dphas);
	}
}
/*=====================================================================*
// OPA step
int OPA(complex<double>* Pum, complex<double>* Sig, complex<double>* Idl, int noStep, int cStage, double fwp, double fws, double fwi, int chirpType, int sProf, int pProf, int iProf, complex<double>* sigPhase, complex<double>* idlPhase, complex<double>* pumPhase) {
// Numerical integration using Runga Kutta 4th sequence
	cout << "entering OPA" << endl;
	double act_EJ6, act_EJ7, act_EJ8;
	complex<double> cds1, cdp1, cdd1, cdk;
	complex<double> cds2, cdp2, cdd2;
	complex<double> cds3, cdp3, cdd3;
	complex<double> cds4, cdp4, cdd4;
	complex<double> csm1, cpm1, cdm1;
	complex<double> csm2, cpm2, cdm2;
	complex<double> csf3, cpf3, cdf3;
	complex<double> clps, clpp, clpi;
	imag(clps) = alps; 
	imag(clpp) = alpp;
	imag(clpi) = alpi;
	complex<double>* ez1 = new complex<double>[nt];
	complex<double>* ez2 = new complex<double>[nt];

	cdk = c1;
	complex<double> cdkl;
	for (int m=0;m<noStep;m++){
		cdkl = cdk;
		for (j=0;j<nt;j++) {
			cdk = cdkl;
			cds1 = clps*Pum[j]*conj(Idl[j])/cdk;
			cdp1 = clpp*Sig[j]*Idl[j]*cdk;
			cdd1 = clpi*Pum[j]*conj(Sig[j])/cdk;

			csm1 = Sig[j]+cds1/2.0;
			cpm1 = Pum[j]+cdp1/2.0;
			cdm1 = Idl[j]+cdd1/2.0;
		
			cdk = cdk*polar(1.0, dk*dzcm*0.5);
			
			cds2 = clps*cpm1*conj(cdm1)/cdk;
			cdp2 = clpp*csm1*cdm1*cdk;
			cdd2 = clpi*cpm1*conj(csm1)/cdk;
			
			csm2 = Sig[j]+cds2/2.0;
			cpm2 = Pum[j]+cdp2/2.0;
			cdm2 = Idl[j]+cdd2/2.0;
			
			cds3 = clps*cpm2*conj(cdm2)/cdk;
			cdp3 = clpp*csm2*cdm2*cdk;
			cdd3 = clpi*cpm2*conj(csm2)/cdk;
			
			csf3 = Sig[j]+cds3;
			cpf3 = Pum[j]+cdp3;
			cdf3 = Idl[j]+cdd3;
	
			cdk = cdk*polar(1.0, dk*dzcm*0.5);

			cds4 = clps*cpf3*conj(cdf3)/cdk;
			cdp4 = clpp*csf3*cdf3*cdk;
			cdd4 = clpi*cpf3*conj(csf3)/cdk;

			Pum[j] = Pum[j]+(cdp1+cdp4+2.0*(cdp2+cdp3))/6.0;
			Sig[j] = Sig[j]+(cds1+cds4+2.0*(cds2+cds3))/6.0;
			Idl[j] = Idl[j]+(cdd1+cdd4+2.0*(cdd2+cdd3))/6.0;

			if(j==nt/2) {
				if(cStage==1) {
					writeToFile("asig1.dat", m, abs(Sig[j])); 
					writeToFile("apum1.dat", m, abs(Pum[j])); 
					writeToFile("aidl1.dat", m, abs(Idl[j])); 
				}
                else if(cStage==2) {
					writeToFile("asig2.dat", m, abs(Sig[j]));
					writeToFile("apum2.dat", m, abs(Pum[j])); 
					writeToFile("aidl2.dat", m, abs(Idl[j])); 
				}
				else {
					writeToFile("asig3.dat", m, abs(Sig[j]));
					writeToFile("apum3.dat", m, abs(Pum[j])); 
					writeToFile("aidl3.dat", m, abs(Idl[j]));
				}
			}
		}
		cout.precision(15);
		act_EJ6 = cInt(Pum, fwp);
		cout << "Energy of Pump:" << act_EJ6*1e3 << "mJ";
		act_EJ7 = cInt(Sig, fws);
		cout << "\t\t Signal:" << act_EJ7*1e3 << "mJ";
		act_EJ8 = cInt(Idl, fwi);
		cout << "\t\tEnergy of Idler:" << act_EJ8*1e3 << "mJ";
		cout << "\t\tTotal:" <<  act_EJ6*1e3+act_EJ7*1e3+act_EJ8*1e3 << " mJ" << endl;	

		if(chirpType==1) {
			disperse(Pum,pumPhase,nt,pProf);
			disperse(Sig,sigPhase,nt,sProf);
			disperse(Idl,idlPhase,nt,iProf);
		}
		// Calculating non-linear phase shift - function for this at one point.
		double n2 = nlindx(cCryst);
		nlshift(cCryst, Pum, fwp, n2);
		nlshift(cCryst, Sig, fws, n2);
		nlshift(cCryst, Idl, fwi, n2);	
	}
}



/*=====================================================================*/
//double phaseWrapper(double angInRad) {
//// Wraps phase between -Pi, Pi
//	return angInRad+tPi/2-tPi*floor((angInRad+tPi/2)/tPi)-tPi/2;
//}
	/*
				double *valami;
				valami = new double[nt];
				double *lanyda;
				lanyda = new double[nt];
				for (j=0;j<nt;j++) {
					//valami[j] = abs(spec[j]);
					valami[j] = pow(abs(timeProfile[j]),2);
					//valami[j] = abs(timeProfile[j]);
					//lanyda[j] = abs(spec[j]);
					lanyda[j] = j;
				}
				writeToFile("gau2.dat", lanyda, valami);


	//	p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(timeProfile),reinterpret_cast<fftw_complex*>(spec),FFTW_BACKWARD,FFTW_ESTIMATE);
	//	fftw_execute(p);
*/
	//for (j=0;j<nt;j++) {
	//	temp[j] = j;
	//}
	//writeToFile("test2.dat", temp, profile);
	//exit(0);

	
	/*
	for (j=0;j<5;j++) {
		cout << spek[j] << endl;
	}
	for (j=0;j<nt;j++) {
		//temp[j] = real(profile[j]);
		//temp2[j] = real(spek[j]);
		temp[j] = j;
		temp2[j] = real(profile[j]);
    
		}
	//writeToFile("prospek.dat", temp, temp2);
	//if (proftype!=0) 
	//fftshift(spek,nt); // If artificial needs to be rotated - check for spek read!	
	for (j=0;j<nt;j++) {
		temp[j] = j;
		temp2[j] = real(spek[j]);
    }
	//writeToFile("prospek2.dat", temp, temp2);
	for (j=0;j<nt;j++) {
		//cout << spek[j] << " " << phase[j] << endl;
		//if(j==50) exit(0);
		spek[j] = spek[j]*phase[j];

	}
	
	for (j=0;j<nt;j++) {
		temp[j] = real(phase[j]);
		temp2[j] = imag(phase[j]);
    }
	//writeToFile("prospek3.dat", temp, temp2);
	p = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(spek),reinterpret_cast<fftw_complex*>(profile),FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for (j=0;j<nt;j++) {
		temp[j] = j;
		temp2[j] = real(profile[j]);
    }
	//writeToFile("prospek4.dat", temp, temp2); // Check with read in spec as well!
	fftshift(profile,nt);
	for (j=0;j<nt;j++) {
		temp[j] = j;
		temp2[j] = real(profile[j]);
    }
	//writeToFile("prospek5.dat", temp, temp2);
	for (j=0;j<nt;j++) {
		real(profile[j]) = real(profile[j])/sqrt(nt); 
		imag(profile[j]) = imag(profile[j])/sqrt(nt); 
		//cout << profile[j] << endl;
	}
	for (j=0;j<nt;j++) {
		temp[j] = j;
		temp2[j] = real(profile[j]);
    }
	//writeToFile("prospek6.dat", temp, temp2);
	//FIX: ... here.
	cout << "Exit at disperse" << endl;
	exit(0);*/



/*	double *valami = new double[nt];
	double *valami2 = new double[nt];
	for (j=0;j<nt;j++) {
		real(ez1[j]) = real(Sig[j]);
		imag(ez1[j]) = imag(Sig[j]);
	}
	
	for (j=0;j<nt;j++) {
		valami[j] = real(ez1[j]);
		valami2[j] = imag(ez1[j]);
	}
	
	writeToFile("test2.dat", valami, valami2);
	exit(0);
*/
	//cout.precision(15);
	
	
				//if (j==nt/2-1) {
				//cout << cdd1 << "\t " << cds1 << "\t " << cdp1 << endl;
				//cout << cdd1 << " " << Pum[j] << " " << Sig[j] << " " << cdk << endl;
			//}
						//writeToFile("test2.dat", (double)j, real(Idl[j]));
			//cout.precision(18);
			//if (m == 0) {
			//cout << dk << " " << dzcm << " " << polar(1.0, dk*dzcm*0.5) << endl;
			//cout << cds1 << " " << cdp1 << " " << cdd1 << " " << cdk << " " << j << endl;
			//cout << csm1 << " " << cpm1 << " " << cdm1 << " " << cdk << " " << j << endl;
			//exit(0);}
