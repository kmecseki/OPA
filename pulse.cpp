#include <numbers>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
#include <complex>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <fftw3.h>
#include "pulse.h"
#include "phys_constants.h"
#include "utils.h"

Crystal::Crystal(int ctype, int stage, int noStep, double thdeg, double cryslth, double xcm2, bool first) : 
    m_firstStage(first),
	m_cType(ctype),
	m_stage(stage),
	m_noStep(noStep),
    m_thdeg(thdeg),
	m_cLength(cryslth),
	m_xcm2(xcm2)
    {	
        m_xeff = calc_xeff();
		m_dzcm = m_cLength/((double) noStep);
    }

double Crystal::calc_xeff() {
	switch(m_cType) {
		case 1: 
			std::cout << "\n\n*********************************************\n";
			std::cout << "Stage " << m_stage << ": BBO" << std::endl;
			std::cout << "*********************************************" << std::endl;
			return 4.04;
		case 2: 
			std::cout << "\n\n*********************************************\n";
			std::cout << "Stage " << m_stage << ": LBO" << std::endl;
			std::cout << "*********************************************" << std::endl;
			return 1.676;
		case 3: 
			std::cout << "Error: KDP not fully supported yet!" << std::endl;
			return 0;
	}
	return 0;
}

double Crystal::calcRefInd(double lambda, int oe) {
	// Returns the refractive index.
	// Input wavelength is in nm from the class.
	// uniaxial crystal: oe = 1 for extraordinary, oe = anything else for ordinary wave.
	// biaxial crystal: oe = 1 for x, oe = 2 for y, oe = 3 for z planes

	double min, max;
	double p1, p2, p3, p4, p5;
	lambda = lambda / 1e3; // convert to microns
	int warned = 0;
	switch(m_cType) {
		case 1: // BBO - transparency 0.189 to 3.5	
			min = 0.189;
			max = 3.5;
			if (lambda < min || lambda > max) {
				if (warned==0) {
					std::cout << "Warning: Wavelength outside transmission window" << std::endl;
					warned = 1;
				}
				if (lambda < min) lambda = min;
				else lambda = max;
			}
			if (oe != 1) { // From Dmitriev
				p1 = 2.7359;
				p2 = 0.01878;
				p3 = 0.01822;
				p4 = -0.01354;
			}
			else {
				p1 = 2.3753;
				p2 = 0.01224;
				p3 = 0.01667;
				p4 = -0.01516;
			}
			return std::sqrt(p1+p2/(lambda*lambda-p3)+p4*lambda*lambda);
		case 2: // LBO - transparency 0.155 to 3.2	
			// Using XY plane
			min = 0.155;
			max = 3.2;
			if (lambda < min || lambda > max) {
				if (warned==0) {
					std::cout << "Warning: Wavelength outside transmission window" << std::endl;
					warned = 1;
				}
				if (lambda < min) lambda = min;
				else lambda = max;
			}
			if (oe == 1) {
				p1 = 2.4542;
				p2 = 0.01125;
				p3 = 0.01135;
				p4 = -0.01388;
			}
			else if (oe == 2) {
				p1 = 2.5390;
				p2 = 0.01277;
				p3 = 0.01189;
				p4 = -0.01848;
			}
			else if (oe == 3) {
				p1 = 2.5865;
				p2 = 0.01310;
				p3 = 0.01223;
				p4 = -0.01861;
			}
			return sqrt(p1+p2/(lambda*lambda-p3)+p4*lambda*lambda);
		case 3: // KDP - not fully supported yet
			// Transparency 0.174 to 1.57
			/*  min = 0.174;
			max = 1.57;
			if (lambda < min || lambda > max) {
				if (warned==0) {
					std::cout << "Warning: Wavelength outside transmission window" << std::endl;
					warned = 1;
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
			refInd = sqrt(p1+p2/(lambda*lambda-p3)+p4*lambda*lambda/(lambda*lambda-p5)); */
			return 0;
	}
}

void Crystal::calc_k(double wlpum, double wlsig, double wlidl) {
	// Calculates and sets the instance's wavenumber k. In 1/micron.
	this->m_kSig = m_nOrdSig * 2 * std::numbers::pi * 1e3 / wlsig;
	this->m_kIdl = m_nOrdIdl * 2 * std::numbers::pi * 1e3 / wlidl;
	this->m_kPum = m_nOrdPum * 2 * std::numbers::pi * 1e3 / wlpum;
}

double Crystal::calcPumCo(double wavelength) {
	// The pump wavelength should be in nm
	// TODO: fix these functions that take lambda, it is now part of the class
	double tanThetaSq, nonePum;	
	double nExPum, nOrPum;
	tanThetaSq = std::pow(tan(deg2rad(m_thdeg)),2);
	nOrPum = calcRefInd(wavelength, 2);
	nExPum = calcRefInd(wavelength, 1);
	nonePum = std::pow((nOrPum/nExPum),2);
	return std::sqrt((1 + tanThetaSq) / (1 + nonePum * tanThetaSq));
}

void Crystal::calc_PM(double nColl_deg, double gamma_rad, bool perfectpm) {
	double pa = std::cos(deg2rad(nColl_deg)) * m_nOrdSig/m_sLam; // nm
	double pb = std::cos(gamma_rad) * m_nOrdIdl / m_iLam;
	double pc = m_nOrdPum / m_pLam;
	double pd = m_xOrdPum / m_pLam;
	double pu = std::pow(((pa + pb) / pc),2);
	double pw = std::pow(((pa + pb) / pd),2);
	double puw = (1 - pu) / (pw - 1);
	if (puw < 0) {
		std::cout << "Error: No phase matching angle!" << std::endl;
		exit(0);
	}
	m_pAng = rad2deg(std::atan(std::sqrt(puw)));

	// Phase mismatch and coherence length
	if (perfectpm) {
		m_dk = 0;
		m_cohL = 0;
	}
	else {
		//dk = kPum-kSig-kIdl;
		m_dk = (m_knPum - m_knSig - m_knIdl) * 1e4; // in 1/cm
		m_cohL = std::abs(2 * std::numbers::pi / (2 * m_dk));  // in micron TODO: is dk used anywhere else?
	}
	using namespace std::complex_literals;
	m_pmism = std::exp(1i * m_dk * m_dzcm / 2e4);
	// We could end the scan loop here.
}

void Crystal::setPhaseVel(double dw, double nColl_deg, std::vector<double> nPumj, std::vector<double> nSigj, std::vector<double> nIdlj, double *gtdPum, double *gtdSig, double *gtdIdl) {
	// Phase velocities
	int nt = nPumj.size();
	double phasVelPum = PhysicalConstants::c0 * 1e2 / nPumj[nt/2]; // cm/sec
	double phasVelSig = PhysicalConstants::c0 * 1e2 / nSigj[nt/2];
	double phasVelIdl = PhysicalConstants::c0 * 1e2 / nIdlj[nt/2-1]; // FIX: Uncomment for ori handling
	// phasVelIdl = c*1e2/nIdlj[nt/2]; // FIX: When it is the new way
	// Phase delays
	double phasdtPum = m_cLength * 1e12 / phasVelPum; // ps
	double phasdtSig = m_cLength * 1e12 / phasVelSig;
	double phasdtIdl = m_cLength * 1e12 / phasVelIdl;
	// Group velocities
	double gVelPum = PhysicalConstants::c0 * 1e2 / (nPumj[nt/2] + m_kPum / m_nOrdPum * PhysicalConstants::c0 *(nPumj[nt/2+1] - nPumj[nt/2-1]) / (2.0 * dw * 1e9)); // cm/s
	double gVelSig = PhysicalConstants::c0 * 1e2 / (nSigj[nt/2] + m_kSig / m_nOrdSig * PhysicalConstants::c0 *(nSigj[nt/2+1] - nSigj[nt/2-1]) / (2.0 * dw * 1e9));
	double gVelIdl = PhysicalConstants::c0 * 1e2 / (nIdlj[nt/2] + m_kIdl / m_nOrdIdl * PhysicalConstants::c0 *(nIdlj[nt/2+1] - nIdlj[nt/2-1]) / (2.0 * dw * 1e9)); // Original idler handling is in use
	//gVelIdl = c*1e2/(nIdlj[nt/2]+kIdl/nOrdIdl*c*(nIdlj[nt/2-1]-nIdlj[nt/2+1])/(2.0*dw*1e9)); // If new idler handling is in use
	// FIX: az eggyel fentebb levo sort atirni, ha idler indexet cserelek
	// Group time delays
	*gtdPum = 1e12 * m_dzcm / gVelPum; // ps
	*gtdSig = 1e12 * m_dzcm / gVelSig;
	*gtdIdl = 1e12 * m_dzcm / gVelIdl;

	std::cout << "Phase velocities:" << std::endl;
	std::cout << "\t" << phasVelSig << "\t     " << phasVelPum << "\t     " << phasVelIdl << "\t" << "cm/s" << std::endl;
	std::cout << "Phase delays:" << std::endl;
	std::cout << "\t" << phasdtSig << "\t\t\t" << phasdtPum << "\t\t\t" << phasdtIdl << "\t\t" << "ps" << std::endl;
	std::cout << "Group velocities:" << std::endl;
	std::cout << "\t" << gVelSig << "\t     " << gVelPum << "\t     " << gVelIdl << "\t" << "cm/s" << std::endl;
	std::cout << "Group delays:" << std::endl;
	std::cout << "\t" << *gtdSig << "\t\t\t" << *gtdPum << "\t\t\t" << *gtdIdl << "\t\t" << "ps" << std::endl;

	// Relative group-time delay
	double grDelPumSig = (gtdPum - gtdSig)* m_noStep; // ps
	double grDelPumIdl = (gtdPum - gtdIdl)* m_noStep;
	double grDelSigIdl = (gtdSig - gtdIdl)* m_noStep;
	std::cout << "Relative Group time delays" << std::endl;
	std::cout << " Signal - Pump \t\t" << -grDelPumSig  << " ps \t Positive for slow signal";
	if (-grDelPumSig<0) {
		std::cout << " -> OK! " << std::endl;
	}
	std::cout << " Idler - Pump \t\t" << -grDelPumIdl  << " ps \t Positive for slow idler";
	if (-grDelPumIdl<0) {
		std::cout << " -> OK! " << std::endl;
	}
	std::cout << " Idler - Signal \t" << -grDelSigIdl  << " ps \t Positive for fast signal";
	if (-grDelSigIdl<0) 
		std::cout << " -> OK! " << std::endl;
	// Print phase matching info
	std::cout << "Phase matching summary:" << std::endl;
	std::cout << " Phase matching angle: " << m_pAng << std::endl;
	std::cout << " Angle of propagation: " << m_thdeg << std::endl;
	std::cout << " Coherence length: " << m_cohL << " um" << std::endl;

	// Calculate at which angle k vector triangle closes properly	
	if (nColl_deg!=0) {
		std::cout << "Non collinear geometry" << std::endl;
		double z, alp;
		z = (std::pow(m_knPum,2) - std::pow(m_kIdl,2) + std::pow(m_kSig,2)) / (2 * m_knPum * m_kSig);
		if (z<=1) {
			alp = rad2deg(std::acos(z));
			std::cout << " k vector triangle closes @ alpha =" << alp << " degs" << std::endl;
		}
		else {
			alp = 0;
			std::cout << "Error: k vector triangle will not close for this theta!" << std::endl;
			exit(0);
		}
		// Calculated beta from Ross(15)		
		if (gVelSig/gVelIdl>=1) 
			std::cout << "No bandwidth optimised alpha exists" << std::endl;
		else {
			double beta2, alpha2, alpha3, nPumOpt, theOpt;
			beta2 = std::acos(gVelSig / gVelIdl);
			// From Ross(16) - good approx to the optimal noncoll angle, it will
			// coincide with alpha at optimal phase matching
			alpha2 = rad2deg(std::asin(m_kIdl/m_knPum) * std::sin(beta2));
			// Calculated from Geoff's formula, preferable to Ross(16), as it
			// does not depend on knPum, which is weakly dependent on theta
			alpha3 = std::atan(std::sin(beta2) / (m_kSig / m_kIdl + std::cos(beta2)));
			// Using Ross(16) backwards to get optimum pump refractive index
			nPumOpt = m_kIdl * std::sin(beta2) / (m_kPum / m_nOrdPum * std::sin(alpha3));
			theOpt = rad2deg(std::asin(std::sqrt((std::pow(m_nOrdPum / nPumOpt,2) - 1) / (std::pow(m_nOrdPum / m_xOrdPum,2) - 1))));
			std::cout << "Optimal alpha: " << rad2deg(alpha3) << "\t Optimal theta: " << theOpt << std::endl;
			std::cout << "Actual alpha: " << nColl_deg << "\t Actual theta: " << m_thdeg << std::endl;
			}
		}
		else std::cout << "Collinear geometry" << std::endl;
		//if (mode!=1) exit(0); Legacy
}

void Crystal::makePhaseRelative(const int frame, const double dw, const double gtdSig, const double gtdPum, const double gtdIdl, Pulse &Pump, Pulse &Signal, Pulse &Idler) {
	// Local time frame selection
	// negative sign because of: gtd = -d(phi)/domega to ensure spatial phase of form exp(-ikz)
	int nt = Pump.m_ctimeProf.size();
	double dphm;
	switch (frame) {
		case 1:
			dphm = -dw * gtdSig;
			break;
		case 2:
			dphm = -dw * (gtdSig + gtdPum) / 2;
			break;
		case 3:
			dphm = -dw * gtdPum;
			break;
		default:
			std::cout << "Error: iframe value can be only 1,2 or 3" << std::endl;
	}
	double phiPumjh = Pump.m_Phij[nt/2];
	double phiSigjh = Signal.m_Phij[nt/2];
	double phiIdljh = Idler.m_Phij[nt/2];
	// Remove central phases and phase gradients
	using namespace std::complex_literals;
	for (int j=0; j<nt; j++) {
		Pump.m_Phij[j] = Pump.m_Phij[j] * 1e4 - (j - nt / 2) * dphm * 1e3 - phiPumjh * 1e4;
		Signal.m_Phij[j] = Signal.m_Phij[j] * 1e4 - (j - nt / 2) * dphm * 1e3 - phiSigjh * 1e4;
		Idler.m_Phij[j] = Idler.m_Phij[j] * 1e4 - (j - nt / 2) * dphm * 1e3 - phiIdljh * 1e4;
		Pump.m_cPhij[j] = std::polar(1.0, Pump.m_Phij[j]); // TODO: add this to pulse class
		Signal.m_cPhij[j] = std::polar(1.0, Signal.m_Phij[j]);
		Idler.m_cPhij[j] = std::exp(1i * Idler.m_Phij[j]); //does the same, as it should tested
	}
	// TODO: attempt to zero pump, but it is not used anywhere. removed for now.
	//m_cPum = 0;
	//if (m_firstStage) {
	//	m_cSig = 0;
	//	m_cIdl = 0;
	//	}
}

// OPA step
int Crystal::OPA(Pulse &Pum, Pulse &Sig, Pulse &Idl, double dtps, int chirpType) {

	//(complex<double>* Pum, complex<double>* Sig, complex<double>* Idl, int noStep, int cStage, double fwp, double fws, double fwi, int chirpType, int sProf, int pProf, int iProf, complex<double>* sigPhase, complex<double>* idlPhase, complex<double>* pumPhase) {
	// Numerical integration using Runga Kutta 4th sequence
	std::cout << "entering OPA" << std::endl;
	//double act_EJ6, act_EJ7, act_EJ8;
	std::complex<double> cds1, cdp1, cdd1, cdk;
	std::complex<double> cds2, cdp2, cdd2;
	std::complex<double> cds3, cdp3, cdd3;
	std::complex<double> cds4, cdp4, cdd4;
	std::complex<double> csm1, cpm1, cdm1;
	std::complex<double> csm2, cpm2, cdm2;
	std::complex<double> csf3, cpf3, cdf3;
	std::complex<double> clps, clpp, clpi;
	clps.imag(Sig.m_alp);
	clpp.imag(Pum.m_alp);
	clpi.imag(Idl.m_alp);
	//complex<double>* ez1 = new complex<double>[nt];
	//complex<double>* ez2 = new complex<double>[nt];
	int nt = Sig.m_ctimeProf.size();
	std::complex<double> cdk(1.0, 0.0);
	std::complex<double> cdkl;
	for (int m=0; m<m_noStep; m++){
		cdkl = cdk;
		for (int j=0; j<nt; j++) {
			cdk = cdkl;
			cds1 = clps * Pum.m_ctimeProf[j] * std::conj(Idl.m_ctimeProf[j]) / cdk;
			cdp1 = clpp * Sig.m_ctimeProf[j] * Idl.m_ctimeProf[j] * cdk;
			cdd1 = clpi * Pum.m_ctimeProf[j] * std::conj(Sig.m_ctimeProf[j]) / cdk;

			csm1 = Sig.m_ctimeProf[j] + cds1 / 2.0;
			cpm1 = Pum.m_ctimeProf[j] + cdp1 / 2.0;
			cdm1 = Idl.m_ctimeProf[j] + cdd1 / 2.0;
		
			cdk = cdk * std::polar(1.0, m_dk * m_dzcm * 0.5);
			
			cds2 = clps * cpm1 * std::conj(cdm1) / cdk;
			cdp2 = clpp * csm1 * cdm1 * cdk;
			cdd2 = clpi * cpm1 * std::conj(csm1) / cdk;
			
			csm2 = Sig.m_ctimeProf[j] + cds2 / 2.0;
			cpm2 = Pum.m_ctimeProf[j] + cdp2 / 2.0;
			cdm2 = Idl.m_ctimeProf[j] + cdd2 / 2.0;
			
			cds3 = clps * cpm2 * std::conj(cdm2) / cdk;
			cdp3 = clpp * csm2 * cdm2 * cdk;
			cdd3 = clpi * cpm2 * std::conj(csm2) / cdk;
			
			csf3 = Sig.m_ctimeProf[j] + cds3;
			cpf3 = Pum.m_ctimeProf[j] + cdp3;
			cdf3 = Idl.m_ctimeProf[j] + cdd3;
	
			cdk = cdk * std::polar(1.0, m_dk * m_dzcm * 0.5);

			cds4 = clps * cpf3 * std::conj(cdf3) / cdk;
			cdp4 = clpp * csf3 * cdf3 * cdk;
			cdd4 = clpi * cpf3 * std::conj(csf3) / cdk;

			Pum.m_ctimeProf[j] = Pum.m_ctimeProf[j] + (cdp1 + cdp4 + 2.0 * (cdp2 + cdp3)) / 6.0;
			Sig.m_ctimeProf[j] = Sig.m_ctimeProf[j] + (cds1 + cds4 + 2.0 * (cds2 + cds3)) / 6.0;
			Idl.m_ctimeProf[j] = Idl.m_ctimeProf[j] + (cdd1 + cdd4 + 2.0 * (cdd2 + cdd3)) / 6.0;

			/*if(j==nt/2) {
				if(cStage==1) {
					writeToFile("asig1.dat", m, abs(Sig.m_ctimeProf[j])); 
					writeToFile("apum1.dat", m, abs(Pum.m_ctimeProf[j])); 
					writeToFile("aidl1.dat", m, abs(Idl.m_ctimeProf[j])); 
				}
                else if(cStage==2) {
					writeToFile("asig2.dat", m, abs(Sig.m_ctimeProf[j]));
					writeToFile("apum2.dat", m, abs(Pum.m_ctimeProf[j])); 
					writeToFile("aidl2.dat", m, abs(Idl.m_ctimeProf[j])); 
				}
				else {
					writeToFile("asig3.dat", m, abs(Sig.m_ctimeProf[j]));
					writeToFile("apum3.dat", m, abs(Pum.m_ctimeProf[j])); 
					writeToFile("aidl3.dat", m, abs(Idl[j]));
				}
			}*/
		}
		//cout.precision(15);
		double act_EJ6 = Pum.cInt(Pum.m_ctimeProf, Pum.m_fw, dtps);
		std::cout << "Energy of Pump:" << act_EJ6 * 1e3 << "mJ";
		double act_EJ7 = Sig.cInt(Sig.m_ctimeProf, Sig.m_fw, dtps); // TODO: do we really need to pass these? Check other uses.
		std::cout << "\t\t Signal:" << act_EJ7 * 1e3 << "mJ";
		double act_EJ8 = Idl.cInt(Idl.m_ctimeProf, Idl.m_fw, dtps);
		std::cout << "\t\tEnergy of Idler:" << act_EJ8 * 1e3 << "mJ";
		std::cout << "\t\tTotal:" <<  act_EJ6 * 1e3 + act_EJ7 * 1e3 + act_EJ8 * 1e3 << " mJ" << std::endl;	

		if(chirpType==1) {
			Pum.disperse();
			Sig.disperse();
			Idl.disperse();
		}
		// Calculating non-linear phase shift - function for this at one point.
		double n2 = nlindx(cCryst);
		nlshift(cCryst, Pum.m_ctimeProf, fwp, n2);
		nlshift(cCryst, Sig, fws, n2);
		nlshift(cCryst, Idl, fwi, n2);	
	}
}

Pulse::Pulse(double dtL, double dtT, double EJ, double Xcm2, double tc, int prof, int nt) :
    m_dtL(dtL),
	m_dtT(dtT),
	m_EJ(EJ),
	m_Xcm2(Xcm2),
	m_tc(tc),
	m_prof(prof),
	m_nt(nt),
	m_absTP(nt),
	m_ctimeProf(nt),
	m_lambdaj(nt)
	{}

double Pulse::calc_omega0() {
	// Calculate angular frequency using the wavelength
	return 2 * std::numbers::pi * PhysicalConstants::c0 / (1e6 * m_wavelength); // 1/fs
}

void Pulse::calc_limits(int nt, double dw) {
	m_lam1 = 2 * std::numbers::pi * PhysicalConstants::c0 * 1e-6 / ( m_omega0 + nt / 2 * dw);
	m_lam2 = 2 * std::numbers::pi * PhysicalConstants::c0 * 1e-6 / ( m_omega0 - nt / 2 * dw);	   
}

double Pulse::rInt(double dtps) {
	// Calculates energy within an intensity profile
	double area = 0;
	int nt = m_absTP.size();
	for (int j=0; j<nt; j++) {
		area = area + m_absTP[j];
	}
	area = area * dtps * 1e-12;
	return area;
}

double Pulse::cInt(std::vector<std::complex<double>> cfield, double fww, double dtps) {
// Calculates energy within an intensity profile - given complex field
	double area = 0;
	int nt = cfield.size();
	std::cout << "cfield size CHECK" << nt << std::endl;
	for (int j=0; j<nt; j++) {
		area = area + std::pow(std::abs(cfield[j]), 2);
	}
	area = area * dtps * 1e-12 * fww;
//	cout << "area" << area << endl;
	return area;
}

void Pulse::GenProfile(Crystal &stage, double dtps, double tlead, int chirpType, double chp1, double chp2, double chd1, double chd2, double dw) {
	//(complex<double> *timeProfile, int profType, double tl, double tt, double xcm2, double EJ, double fw, double tCoh, double cpp1, double cpp2, double cpd1, double cpd2) {
	// This function generates a skewed gaussian or sech2 profile (skewed in time or spec)
	// Adds background, noise and chirp if needed.
	int nt = m_nt;
	switch (m_prof) {
		case 0: {
			// Reading profile from file (only for signal!)
			// Requires a file with wavelengths (nm) and spectrum data
			std::vector<double> spectrum;
			std::vector<double> lambdas;
			std::vector<double> intpSpec;
			const char* file = "rspec.txt";
			std::ifstream filestr;
			double bck = 700; // TODO: Hardcoded for now as bck is always the same on OceanOptics sp I have.
			filestr.open (file, std::fstream::in);
			if(filestr.is_open()) {
				double lambda, intens;
				while(filestr >> lambda >> intens) {
					lambdas.push_back(lambda);
					spectrum.push_back(intens-bck);
				}
			filestr.close();
			}
			else {
				std::cout << "Error opening signal profile data file." << std::endl;
				exit(0);
			}
			std::cout << lambdas.size() << " elements were read from spectrum file." << std::endl;
			//reverse(lambdas->begin(),lambdas->end()); 
			//reverse(spectrum->begin(),spectrum->end());
			int NT = lambdas.size();

			// Correction for interpolation outside limits
			if (lambdas[0]>m_lam1);
				lambdas[0] = m_lam1;
			if (lambdas[NT-1]<m_lam2)
				lambdas[NT-1] = m_lam2;
			// Interpolating read data to make it equidistant in freq.
			gsl_interp_accel *acc = gsl_interp_accel_alloc();
			gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, NT);
			gsl_spline_init(spline, lambdas.data(), spectrum.data(), NT);
			for (int j=0; j<nt; j++) {
				intpSpec.push_back(std::sqrt(gsl_spline_eval(spline,m_lambdaj[nt - 1 - j], acc)));
			}
			gsl_spline_free(spline);
			gsl_interp_accel_free(acc);
			
			// Reverse in place to have the frequency in the right order.
			std::reverse(intpSpec.begin(), intpSpec.end());
			
			// Fourier transformation of data
			fftw_complex *temp = new fftw_complex[nt];
			cvector_to_fftw(nt, m_ctimeProf, temp);
			fftw_plan p0 = fftw_plan_dft_r2c_1d(nt, intpSpec.data(), temp, FFTW_ESTIMATE);
			fftw_execute(p0);
			fftw_to_cvector(nt, temp, m_ctimeProf);
			delete temp;

			// Creating full length of pulse
			double *rea;
			double *ima;
			rea = new double[nt/2+1];
			ima = new double[nt/2+1];
			for(int j=0; j<(nt/2+1); j++) {
				rea[j] = std::real(m_ctimeProf[j]);
				ima[j] = std::imag(m_ctimeProf[j]);
			}
			for(int j=0; j<nt; j++) {
				double R, I;
				R = rea[std::abs(j - nt / 2)];
				if (j<nt/2) { // Check for shift and change to nt/2+1 if necessary..
					I = ima[abs(j - nt / 2)];
				}
				else {
					I = -ima[abs(j - nt / 2)]; // because it only gives back half, and imag part has to be asymmetric
				}
				m_ctimeProf[j] = std::complex(R, I);
			}
			delete rea;
			delete ima;
			double maxx = 0;
			for (int j=0; j<nt; j++) {
				// To calculate the 2nd power of the electric field profile in time
				double number = std::pow(std::abs(m_ctimeProf[j]), 2);
				if (number>maxx) {
					maxx = number;
				}
				m_absTP[j] = number;

			}
			// Normalise
			for (int j=0; j<nt; j++) {
				m_absTP[j] = m_absTP[j]/maxx;
			}
			break;
		}
		
		case 1: 
			// Gaussian temporal intensity profile	
			for(int j=0; j<nt/2; j++) { // Modified original code here to get overlapping profiles (dtps*j instead of j+1)
				m_absTP[j] = PhysicalConstants::small + std::exp(-(std::log(2) * std::pow(((dtps * j - tlead) / stage.m_dtPumpL), 2))); // Make sure tl is dtpumpl
			}
			for(int j=nt/2; j<nt; j++) {
				m_absTP[j] = PhysicalConstants::small + std::exp(-(std::log(2) * std::pow(((dtps * j - tlead) / stage.m_dtPumpT), 2)));
			}
			break;
			
		case 2:
		// Sech squared temporal intensity profile
			for(int j=0; j<nt/2; j++) {
				m_absTP[j] = PhysicalConstants::small + std::pow(2 / (std::exp(-(dtps * (j+1) - tlead) / stage.m_dtPumpL) + 1 / std::exp(-(dtps * (j+1) - tlead) / stage.m_dtPumpL)), 2);
			}
			for(int j=nt/2;j<nt;j++) {
				m_absTP[j] = PhysicalConstants::small + std::pow(2 / (std::exp(-(dtps * (j+1) - tlead) / stage.m_dtPumpT) + 1 / std::exp(-(dtps * (j+1) - tlead) / stage.m_dtPumpT)), 2);
			}
			break;

		default:
			std::cout << "Unsupported profile type." << std::endl;
			exit(0);
	}
	// Basic intensity profile formed in time -> absTP
	// Next: Power scaling
	double act_EJ;
	act_EJ = rInt(dtps);
	for (int j=0; j<nt; j++) {
		m_absTP[j] = (m_EJ * m_absTP[j] / act_EJ);
	}
	auto max = std::max_element(m_absTP.begin(), m_absTP.end());
	//	cout << EJ/act_EJ << " " << EJ << " " << act_EJ <<  endl;
	//	exit(0);
	act_EJ = this->rInt(dtps);
	std::cout << "First energy check: Needed:" << m_EJ * 1e3 << "mJ -> Actual:" << act_EJ * 1e3 << "mJ" << std::endl;
	// Complex field profile has to be scaled as well
	m_fw = 0.5 * PhysicalConstants::c0 * PhysicalConstants::eps0 * (stage.m_nOrdPum * stage.m_coePum) * stage.m_xcm2; // Check units if sg is wrong
	if (m_prof==0) {
		// In this case the spectrum is read from file - (only for signal)
		double act_EJ2 = cInt(m_ctimeProf, m_fw, dtps);
		for (int j=0; j<nt; j++) {
			(m_ctimeProf[j]) = (m_ctimeProf[j]) * sqrt(m_EJ / act_EJ2);
		}
	}
	else {
	// For simulated profiles
		for (int j=0; j<nt; j++) {
			//timeProfile[j] = polar(sqrt((absTP[j]/fw)),0.0);
			m_ctimeProf[j].real(std::sqrt(m_absTP[j] / m_fw));
		}
	}
	//	cout << timeProfile[1] << " " << fw <<  endl;
	double act_EJ3 = cInt(m_ctimeProf, m_fw, dtps);
	//	cout << absTP[nt/2] << " " << timeProfile[nt/2] << " " << fw*dtps*1e-12 << " " << act_EJ3/(fw*dtps*1e-12) << endl;
	//	exit(0);
	std::cout << "Second energy check: Needed:" << m_EJ * 1e3 << "mJ -> Actual:" << act_EJ3 * 1e3 << "mJ" << std::endl;
	
	// Noisy pulse - not background noise - Not tested
	if (m_tc!=0) {
		std::cout << "Error: Noise not implemented!" << std::endl;
		exit(0);
		/*m_tc = 50; // Re-write to a sensible number before test completes.
		double fspec;
		std::vector<complex<double>> cva(nt);
		fspec=dw*dw*tCoh*tCoh/(8*log(2));
		cva=noise(cva,profType,nt,fspec);
		for (j=0;j<nt;j++) {
			imag(m_ctimeProf[j])=imag(m_ctimeProf[j])*imag(cva[j])/(abs(cva[j])+1e-20);
			cout << cva[j] << endl;
		}*/
	}
	// Chirping procedure
	if(chirpType==1) {
		if (chp1!=0 || chp2!=0) {
			// In this chirping procedure the spectral bandwidth is unchanged
			// - pulse duration is increased 
			std::cout << "Normal chirping procedure" << std::endl;
			chirper_norm(m_ctimeProf, m_prof, chp1, chp2, dw);
			double act_EJ4 = cInt(m_ctimeProf, m_fw, dtps);
			std::cout << "Energy check after chirping: Needed:" << m_EJ * 1e3 << "mJ -> Actual:" << act_EJ4 * 1e3 << "mJ" << std::endl;	
			// FIXME: need to check why it changes so much - numerical issue??
			if (m_EJ!=act_EJ4) {
				std::cout << "Inconsistent energies! Correcting..." << std::endl;
				for (int j=0; j<nt; j++) {
					m_ctimeProf[j] = std::polar((std::abs(m_ctimeProf[j]) * std::sqrt(m_EJ / act_EJ4)), std::arg(m_ctimeProf[j]));
				}
				double act_EJ5 = cInt(m_ctimeProf, m_fw, dtps);
				std::cout << "Energy check again: Needed:" << m_EJ * 1e3 << "mJ -> Actual:" << act_EJ5 * 1e3 << "mJ" << std::endl;	
			}
		}
	}
	else {
		if (chd1!=0 || chd2!=0) {
			std::cout << "\tUsing direct chirp" << std::endl;
			chirper_direct(m_ctimeProf, dtps, chd1, chd2);
		}
	}
}

int Pulse::chirper_norm(std::vector<std::complex<double>> &timeProfile, int profType, double phi2, double phi3, double dw) {
	// Chirping function - normal chirp
	// FFT to get to spectral domain:
	std::ifstream ifile2("fftplan2.dat");
	fftw_plan p;
	int nt = timeProfile.size();
	std::unique_ptr<fftw_complex[]> in = std::make_unique<fftw_complex[]>(nt);
	std::unique_ptr<fftw_complex[]> out = std::make_unique<fftw_complex[]>(nt);
	std::vector<std::complex<double>> spec(nt);
	cvector_to_fftw(nt, timeProfile, in.get());
	// TODO: does planning already transform?
	if(!ifile2) {
		std::cout << "Planning fastest FFT algorithm for this CPU... This will take a while, but needs to run only once on this machine." << std::endl;
		p = fftw_plan_dft_1d(nt, in.get(), out.get(), FFTW_BACKWARD, FFTW_PATIENT);
		std::cout << "Plan created, saving it as 'fftplan2.dat'" << std::endl;
		if(!fftw_export_wisdom_to_filename("fftplan2.dat"))
			std::cout << "Error writing plan" << std::endl;
	}
	else {
		std::cout << "Found an FFT plan file. Reading the FFT plan from file." << std::endl; 
		fftw_import_wisdom_from_filename("fftplan2.dat");
		p = fftw_plan_dft_1d(nt, in.get(), out.get(), FFTW_BACKWARD, FFTW_PATIENT);
	}
	fftw_execute(p);
	if (profType!=0) fftshift(out.get(), nt);
	fftw_to_cvector(nt, out.get(), spec); // TODO does this work?
	for (int j=0; j<nt; j++) {
		// this is actually in the spectral domain, just reusing memory.
		// Applying chirp in the spectral domain.
		spec[j] = std::polar(std::abs(spec[j]), phi2 / 2 * std::pow((-nt / 2 + j) * dw, 2) + phi3 / 6 * std::pow((-nt / 2 + j) * dw, 3));
	}
	// Converting back to time
	cvector_to_fftw(nt, spec, in.get());
	p = fftw_plan_dft_1d(nt, in.get(), out.get(), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_to_cvector(nt, out.get(), timeProfile);
	for (int j=0; j<nt; j++) {
	// sqrt(nt)*sqrt(nt) division required again after fft
		timeProfile[j].real(timeProfile[j].real() / nt); 
		timeProfile[j].imag(timeProfile[j].imag() / nt);
		//timeProfile[j] = polar(abs(timeProfile[j])/nt, arg(timeProfile[j])); 
	}
	fftshift(timeProfile, nt);
}

int Pulse::chirper_direct(std::vector<std::complex<double>> &timeProfile, double dtps, double chp1, double chp2) {
	// This function applies linear (and nonlinear) chirp directly to the complex field
	// in the time domain
	// TODO: Maybe fftshift before and after is needed, check..
	int nt = timeProfile.size();
	std::cout << "This type of chirping has not been tested yet." << std::endl;
	for (int j=0; j<nt; j++) {
		double par1 = 2 * std::numbers::pi * std::pow(dtps, 2) * chp1 / 2000;
		double par2 = 2 * std::numbers::pi * std::pow(dtps, 3) * chp2 / 6000; // Fourier expansion
		timeProfile[j] = std::polar(std::abs(timeProfile[j]), j * j * (par1 + par2 * j * j));
	}
}

int Pulse::spectrum(const char *ofname) {
// This function calculates the pulse spectrum, FWHM, and writes into file
	int nt = m_ctimeProf.size();
	std::vector<std::complex<double>> spek(nt);
	std::vector<double> tempspek(nt);
	double FWHM;
	fftw_complex timeProf2[nt];
	cvector_to_fftw(nt, m_ctimeProf, timeProf2);
	fftw_complex spek2[nt];
	cvector_to_fftw(nt, spek, spek2);
	fftw_plan p = fftw_plan_dft_1d(nt, timeProf2, spek2, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	fftw_to_cvector(nt, spek2, spek);
	if (m_prof!=0)
		fftshift(spek,nt);
	for (int j=0; j<nt; j++) {
		spek[j] = std::polar(std::abs(spek[j]) / std::sqrt(nt), arg(spek[j]));
		tempspek[j] = std::abs(spek[j]);
    }
    FWHM = get_FWHM(spek, m_lambdaj); // This can be written out if needed. - not tested
    std::cout << "FWHM : " << FWHM << " nm" << std::endl;
    writeToFile(ofname, m_lambdaj, tempspek);    
}

int Pulse::disperse() {
	// Applies dispersive phase

	int nt = m_ctimeProf.size();
	fftw_complex out[nt];
	fftw_complex in[nt];
	std::vector<std::complex<double>> spek(nt);

	fftshift(m_ctimeProf, nt);
	cvector_to_fftw(nt, m_ctimeProf, in);
	fftw_plan p = fftw_plan_dft_1d(nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_to_cvector(nt, out, spek);
	for (int j=0; j<nt; j++) {
		spek[j].real(spek[j].real() / std::sqrt(nt)); 
		spek[j].imag(spek[j].imag() / std::sqrt(nt)); 
	}

	fftshift(m_cPhij, nt);
	for (int j=0; j<nt; j++) {
		spek[j] = spek[j] * m_cPhij[j]; 
	}
	fftshift(m_cPhij, nt);

	cvector_to_fftw(nt, spek, in);
	fftw_plan p2 = fftw_plan_dft_1d(nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p2);
	fftw_to_cvector(nt, out, m_ctimeProf);
	for (int j=0; j<nt; j++) {
		m_ctimeProf[j].real(m_ctimeProf[j].real() / std::sqrt(nt)); 
		m_ctimeProf[j].real(m_ctimeProf[j].imag() / std::sqrt(nt)); 
	}
	fftshift(m_ctimeProf, nt);
}