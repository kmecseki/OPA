#include <numbers>
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
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

void fftw_to_cvector(int nt, fftw_complex* in, std::vector<std::complex<double>> &out) {
	// Convert fftw vector to complex vector
	for (int j=0; j<nt; j++) {
		out[j] = std::complex(in[j][0],in[j][1]);
	}
}

void fftshift(std::vector<std::complex<double>> &vekt, int size) {
// FFT-shift
	std::complex<double> temp, temp2;
	for(int j=0; j<size/4; j++){
		temp = vekt[j];
		vekt[j] = vekt[size/2-(1+j)];
		vekt[size/2-(1+j)] = temp;
		temp2 = vekt[j+size/2];
		vekt[j+size/2] = vekt[size-(1+j)];
		vekt[size-(1+j)] = temp2;
	}
}

void fftshift(fftw_complex* vekt, int size) {
// FFT-shift
	fftw_complex temp;
	fftw_complex temp2;
	for(int j=0; j<size/4; j++){
		temp[0] = vekt[j][0];
		temp[1] = vekt[j][1];
		vekt[j][0] = vekt[size/2-(1+j)][0];
		vekt[j][1] = vekt[size/2-(1+j)][1];
		vekt[size/2-(1+j)][0] = temp[0];
		vekt[size/2-(1+j)][1] = temp[1];
		temp2[0] = vekt[j+size/2][0];
		temp2[1] = vekt[j+size/2][1];
		vekt[j+size/2][0] = vekt[size-(1+j)][0];
		vekt[j+size/2][1] = vekt[size-(1+j)][1];
		vekt[size-(1+j)][0] = temp2[0];
		vekt[size-(1+j)][1] = temp2[1];
	}
}

double FindMax(std::vector<double> &vek) {
// Finds maximum of vector 'vek'.
	int size = vek.size();
	double maxi=0;
	for (int j=0;j <size; j++) {
		if (vek[j]>maxi)
			maxi = vek[j];
	}
	return maxi;
}

double FindMax(const std::vector<std::complex<double>> &vek) {
// Finds maximum of the complex vector 'vek'.
	int size = vek.size();
	double maxi=0;
	for (int j=0;j <size; j++) {
		if (vek[j].real()>maxi)
			maxi = vek[j].real();
	}
	return maxi;
}

double get_FWHM(const std::vector<std::complex<double>> &profile, std::vector<double> wl) {
// Finds FWHM in indeces, need to scale
	double max = FindMax(profile);
	int k1, k2;
	int nt = profile.size();
	for (int j=0; j<nt; j++) {
		if (std::abs(profile[j]) < max / 2 && std::abs(profile[j+1]) >= max / 2) {
			k1 = j;
			break;
		}
		else k1 = 0;
	}
	for (int j=nt; j>0; j--) {
		if (std::abs(profile[j]) < max / 2 && std::abs(profile[j-1]) >= max / 2) {
			k2 = j;
			break;
		}
		else k2 = 0;
	}
	return std::abs(wl[k1]-wl[k2]);
}

void writeToFile(const char *ofname, std::vector<double> &data1, std::vector<std::complex<double>> &data2) {
// Opens a file with name ofname, and outputs complex data
	std::ofstream filestr;
	filestr.open(ofname, std::ios::out | std::ios::app);
	filestr.precision(15);
	int nt = data1.size();
	if (filestr.is_open()) {
		for (int j=0; j<nt; j++) {
			filestr << data1[j] << "\t" << real(data2[j]) << "\t" << imag(data2[j]) << "\n";
		}
	}
	else {
		std::cout << "Error opening output file!" << std::endl;
	}
	filestr.close();
}

void writeToFile(const char *ofname, std::vector<double> &data1, std::vector<double> &data2) {
// Opens a file with name ofname, and outputs complex data
	std::ofstream filestr;
	filestr.open(ofname, std::ios::out | std::ios::app);
	filestr.precision(15);
	int nt = data1.size();
	if (filestr.is_open()) {
		for (int j=0; j<nt; j++) {
			filestr << data1[j] << "\t" << data2[j] << "\n";
		}
	}
	else {
		std::cout << "Error opening output file!" << std::endl;
	}
	filestr.close();
}


/*

#include <vector>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

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


double FindMax(complex<double> *cfield, int s) {
// Finds maximum of vector 'vek' of size 's'.
	double maxi=0;
	for (j=0;j<s;j++) {
		if(abs(cfield[j])>maxi)
			maxi=abs(cfield[j]);
	}
	return maxi;
}

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
