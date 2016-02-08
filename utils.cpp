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
extern const double tPi, small;
extern double thdeg;
extern double *pLambdaj, *iLambdaj, *sLambdaj, *t;
extern double *absTP;
extern double pLam1, pLam2, sLam1, sLam2, iLam1, iLam2;
extern double dtps, tlead;
//extern fftw_complex *timeProfile;
extern complex<double> *timeProfile;

// Function declarations
int errorhl(int);
int warninghl(int);
int openfile(const char, vector<double>*);
int writeToFile(const char, double*, double*);
double xeffc(int, int);
double calcRefInd(int, double, int);
double deg2rad(double);
double rad2deg(double);
double CalcPumCo(double);
int GenProfile(int, double, double, double, double);
double FindMax(double*, int);
double FindMax(complex<double> *, int); 
double rInt(double*);
double cInt(complex<double>*, double);

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
	}
	exit(0);
}
/*=====================================================================*/
int warninghl(int w) {
	switch(w) {
		case 1: cout << "Warning: Wavelength outside transmission window" << endl; break;		
	}
}
/*=====================================================================*/
int openfile(const char *fname, vector<double>* vekto) {
// This opens a file, and reads contents into 
// vector, with number of elements rows.
	double number;
	char delim = '/';
	ifstream filestr;
	string line;
	filestr.open (fname, fstream::in);
	if(filestr.is_open()) {
		while(!filestr.eof()) {
			getline(filestr,line,delim);
			stringstream input(line);
			while(input >> number) {
				vekto->push_back(number);
			}
		}
	}
	else {
		errorhl(1);
	}
	filestr.close();
}
/*=====================================================================*/
int writeToFile(const char *ofname, double *data1, double *data2) {
// Opens a file with name ofname, and outputs data
	ofstream filestr;
	filestr.open(ofname, ios::out | ios::app);
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
double xeffc(int crystnum, int stage) {
// Returns xeff.
	double xeff;
	switch(crystnum) {
		case 1: xeff = 4.04; 
			cout << endl << endl << "Stage " << stage << ": BBO" << endl; break;
		case 2: xeff = 1.676; 
			cout << "Stage " << stage << ": LBO" << endl; break;
		case 3: xeff = 0; errorhl(3); break;
	}
	return xeff;
}
/*=====================================================================*/
double calcRefInd(int crystnum, double lambda, int oe) {
// Returns refractive index for wavelength lambda.
// Input wavelength should be in nm.
// oe = 1 for extraordinary, oe = anything else for ordinary wave.
// oe = 1 for x, oe = 2 for y, oe = 3 for z planes in case of biaxial crystal
	double refInd;
	double min, max;
	double p1, p2, p3, p4, p5;
	double lambdax;
	lambda = lambda/1e3;
	lambdax = lambda*lambda;
	switch(crystnum) {
		case 1: // BBO - transparency 0.189 to 3.5	
				min = 0.189;
				max = 3.5;
				if (lambda < min || lambda > max) {
					if (warned==0) {
						warninghl(1); warned = 1;
					}
					if (lambda < min) lambda = min;
					else lambda = max;
					lambdax = lambda*lambda;
				}
				if (oe != 1) { // From Dmitriev
					p1=2.7359;
					p2=0.01878;
					p3=0.01822;
					p4=-0.01354;
				}
				else {
					p1=2.3753;
					p2=0.01224;
					p3=0.01667;
					p4=-0.01516;
				}
				refInd = sqrt(p1+p2/(lambdax-p3)+p4*lambdax);
				break;
		case 2: // LBO - transparency 0.155 to 3.2	
				// Using XY plane
				min = 0.155;
				max = 3.2;
				if (lambda < min || lambda > max) {
					if (warned==0) {
						warninghl(1); warned = 1;
					}
					if (lambda < min) lambda = min;
					else lambda = max;
					lambdax = lambda*lambda;
				}
				if (oe == 1) {
					p1=2.4542;
					p2=0.01125;
					p3=0.01135;
					p4=-0.01388;
				}
				else if (oe == 2) {
					p1=2.5390;
					p2=0.01277;
					p3=0.01189;
					p4=-0.01848;
				}
				else if (oe == 3) {
					p1=2.5865;
					p2=0.01310;
					p3=0.01223;
					p4=-0.01861;
				}
				refInd = sqrt(p1+p2/(lambdax-p3)+p4*lambdax);
				break;
		case 3: // KDP - not fully supported yet
				// Transparency 0.174 to 1.57
			/*  min = 0.174;
				max = 1.57;
				if (lambda < min || lambda > max) {
					if (warned==0) {
						warninghl(1); warned = 1;
					}
					if (lambda < min) lambda = min;
					else lambda = max;
					lambdax = lambda*lambda;
				}
				if (oe != 1) {
					p1=2.259276;
					p2=0.01008956;
					p3=0.012942625;
					p4=13.00522;
					p5=400;
				}
				else {
					p1=2.132668;
					p2=0.008637494;
					p3=0.012281043;
					p4=3.2279924;
					p5=400;
				}
				refInd = sqrt(p1+p2/(lambdax-p3)+p4*lambdax/(lambdax-p5)); */
				refInd = 0; errorhl(3); break;
	}
	return refInd;
}
/*=====================================================================*/
double deg2rad(double deg) {
// Returns input deg value in radians.
	return deg/360*tPi;
}
/*=====================================================================*/
double rad2deg(double rad) {
// Returns input deg value in radians.
	return rad*360/tPi;
}
/*=====================================================================*/
double CalcPumCo(double lambda) {
// pump wavelength should be in nm	
	double tanThetaSq, nonePum, coePum;	
	double nExPum, nOrPum;
	tanThetaSq = pow(tan(deg2rad(thdeg)),2);
	nOrPum = calcRefInd(cCryst, lambda, 2);
	nExPum = calcRefInd(cCryst, lambda, 1);
	nonePum = pow((nOrPum/nExPum),2);
	coePum = sqrt((1 + tanThetaSq)/(1 + nonePum*tanThetaSq));
	return coePum;
}
/*=====================================================================*/
int GenProfile(int profType, double tl, double tt, double xcm2, double EJ, double fw ) {
// This function generates a skewed gaussian or sech2 profile
// Adds background, noise and chirp if needed.
	switch (profType) {
		case 0: {
		// Reading profile from file (only for signal)
		// Requires a file with wavelengths (nm) and spectrum data
			const char* file = "rspec.txt";
			double maxx = 0;
			double number;
			double *intpSpec, *spec_read, *lamb_read;
			vector<double>* spectrum;
			vector<double>* lambdas;
			//double *lamb_read = &(lambdas[0]);
			//double* spec_read = &spectrum[0];
			spectrum = new vector<double>();
			lambdas = new vector<double>();
			intpSpec = new double[nt];
			int sizeRead;
			ifstream filestr;
			string line;
			filestr.open (file, fstream::in);
			if(filestr.is_open()) {
				while(!filestr.eof()) {
					getline(filestr,line);
					stringstream input(line);
					input >> number;
					lambdas->push_back(number);
					//cout << "elso" << number << endl;
					input >> number;
					// 700 is background now FIXME
					spectrum->push_back(number-700);
					//cout << "masodik" << number << endl;
				}
				cout << lambdas->size() << " elements were read from spectrum file." << endl;
				//reverse(lambdas->begin(),lambdas->end()); 
				//reverse(spectrum->begin(),spectrum->end());
				sizeRead = lambdas->size();
				spec_read = new double[sizeRead];
				lamb_read = new double[sizeRead];
				double *lamb_read = &(lambdas->at(0));
				double *spec_read = &(spectrum->at(0));
			// Little correction for interpolation outside limits
				if (lamb_read[0]>sLam1)
					lamb_read[0] = sLam1;
				if (lamb_read[sizeRead-1]<sLam2)
					lamb_read[sizeRead-1] = sLam2;
			// Interpolating read data to make it equidistant in freq.
				gsl_interp_accel *acc = gsl_interp_accel_alloc();
				gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,sizeRead);
				gsl_spline_init (spline, lamb_read, spec_read, sizeRead);
				for (j=0;j<nt;j++) {
					intpSpec[j] = gsl_spline_eval(spline,sLambdaj[nt-1-j],acc);
				}
				gsl_spline_free(spline);
				gsl_interp_accel_free(acc);
			// Fourier transformation of data
				fftw_plan p;
				ifstream ifile("fftplan.dat");
				if(!ifile) {
					cout << "Planning fastest FFT algorithm for this CPU... This will take a while, but needs to run only once on this machine." << endl;
					p = fftw_plan_dft_r2c_1d(nt,intpSpec,reinterpret_cast<fftw_complex*>(timeProfile), FFTW_PATIENT);
					cout << "Plan created, saved as 'fftplan.dat'" << endl;
					if(!fftw_export_wisdom_to_filename("fftplan.dat"))
						cout << "Error writing plan" << endl;
				}
				else {
					cout << "Reading FFT plan from file." << endl; 
					fftw_import_wisdom_from_filename("fftplan.dat");
					p = fftw_plan_dft_r2c_1d(nt,intpSpec,reinterpret_cast<fftw_complex*>(timeProfile), FFTW_PATIENT);
				}
				for(j=0;j<nt;j++) {
					absTP[j] = abs(pow(timeProfile[j],2));
				}
				maxx = FindMax(absTP, nt);
				for(j=0;j<nt;j++) {
					absTP[j] = absTP[j]/maxx;
				}
				fftw_execute(p);
				fftw_destroy_plan(p);
				//fftw_free(timeProfile);
			}
			else
				errorhl(9);
			filestr.close();
			break;
		}
			
		case 1: 
		// Gaussian temporal intensity profile		
			for(j=0;j<nt/2;j++) {
				absTP[j] = small+exp(-(log(2)*pow(((dtps*j-tlead)/tl),2)));
				//timeProfile[j] = sqrt(absTP[j]);
			//	if (timeProfile[j][0] > 50)
			//		timeProfile[j][0] = small;
			}
			for(j=nt/2;j<nt;j++) {
				absTP[j] = small+exp(-(log(2)*pow(((dtps*j-tlead)/tt),2)));
				//timeProfile[j] = sqrt(absTP[j]);
			//	if (timeProfile[j][0] > 50)
			//		timeProfile[j][0] = small;
			}
			/*double *valami;
			valami = new double[nt];
			for (j=0;j<nt;j++) {
				valami[j] = real(timeProfile[j]);
			}
			writeToFile("gau.dat", t, valami);
			break;*/
		case 2:
		// Sech squared temporal intensity profile
			for(j=0;j<nt/2;j++) {
				absTP[j] = small+pow(2/(exp(-(dtps*j-tlead)/tl))+1/exp(-(dtps*j-tlead)/tl),2);
				//timeProfile[j] = sqrt(absTP[j]);
				//if (timeProfile[j] > 50)
					//timeProfile[j] = small;
			}
			for(j=nt/2;j<nt;j++) {
				absTP[j] = small+pow(2/(exp(-(dtps*j-tlead)/tt))+1/exp(-(dtps*j-tlead)/tt),2);
				//timeProfile[j] = sqrt(absTP[j]);
				//if (timeProfile[j] > 50)
					//timeProfile[j] = small;
			}
			break;
	}
// Basic intensity profile formed.
// Next: Power scaling
	double max=0;
	double act_EJ;
	act_EJ = rInt(absTP);
	for (j=0;j<nt;j++) {
		absTP[j] = (EJ*absTP[j]/act_EJ);
	}
	max = FindMax(absTP, nt);
	act_EJ = rInt(absTP);
	cout << "First energy check: Needed:" << EJ*1e3 << "mJ -> Actual:" << act_EJ*1e3 << "mJ" << endl;
// Complex field profile has to be scaled as well
	double act_EJ2 = cInt(timeProfile, fw);
	if (profType >=0) {
	// For simulated profiles
		for (j=0;j<nt;j++) {
			timeProfile[j] = sqrt(absTP[j]/fw);
		}
	}
	else {
	// In case spectrum is read - use only for signal
		for (j=0;j<nt;j++) {
			timeProfile[j] = timeProfile[j]*sqrt(EJ/act_EJ2);
		}
	}
	double act_EJ3 = cInt(timeProfile, fw);
	cout << "Second energy check: Needed:" << EJ*1e3 << "mJ -> Actual:" << act_EJ3*1e3 << "mJ" << endl;

// Background noise




}
/*=====================================================================*/
double FindMax(double *vek, int s) {
// Finds maximum of vector 'vek' of size 's'.
	double maxi=0;
	for (j=0;j<s;j++) {
		if(vek[j]>maxi)
			maxi=vek[j];
	}
	return maxi;
}
/*=====================================================================*/
double FindMax(complex<double> *cfield, int s) {
// Finds maximum of vector 'vek' of size 's'.
	double maxi=0;
	for (j=0;j<s;j++) {
		if(abs(cfield[j])>maxi)
			maxi=abs(cfield[j]);
	}
	return maxi;
}
/*=====================================================================*/
double rInt(double *absTp) {
// Calculates energy within an intensity profile
	double area=0;
	for (j=0;j<nt;j++) {
		area=+absTp[j];
	}
	area=area*dtps*1e-12;
	return area;
}
/*=====================================================================*/
double cInt(complex<double> *cfield, double fww) {
// Calculates energy within an intensity profile - given complex field
	double area=0;
	for (j=0;j<nt;j++) {
		area=+pow(abs(cfield[j]),2);
	}
	area=area*dtps*1e-12*fww;
	return area;
}
