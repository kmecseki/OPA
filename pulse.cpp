#include <numbers>
#include <cmath>
#include <iostream>
#include "pulse.h"
#include "phys_constants.h"
#include "utils.h"

Crystal::Crystal(int ctype, int stage, int noStep, double thdeg, double cryslth) : 
    m_cType(ctype),
	m_stage(stage),
	m_noStep(noStep),
    m_thdeg(thdeg),
	m_cLength(cryslth) 
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

double Crystal::calc_k(double wlpum, double wlsig, double wlidl) {
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

void Crystal::setPhaseVel(int nt, double dw, std::vector<double> nPumj, std::vector<double> nSigj, std::vector<double> nIdlj, double *gtdPum, double *gtdSig, double *gtdIdl) {
	// Phase velocities
	double phasVelPum = PhysicalConstants::c0 * 1e2 / nPumj[nt/2]; // cm/sec
	double phasVelSig = PhysicalConstants::c0 * 1e2 / nSigj[nt/2];
	double phasVelIdl = PhysicalConstants::c0 * 1e2 / nIdlj[nt/2-1]; // FIX: Uncomment for ori handling
	// phasVelIdl = c*1e2/nIdlj[nt/2]; // FIX: When it is the new way
	// Phase delays
	double phasdtPum = m_cLength * 1e12 / phasVelPum; // ps
	double phasdtSig = m_cLength * 1e12 / phasVelSig;
	double phasdtIdl = m_cLength * 1e12 / phasVelIdl;
	// Group velocities
	double gVelPum = PhysicalConstants::c0 * 1e2 / (nPumj[nt/2] + m_kPum / m_nOrdPum * PhysicalConstants::c0 *(nPumj[nt/2+1] - nPumj[nt/2-1]) / (2.0  *dw * 1e9)); // cm/s
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
}

Pulse::Pulse(double dtL, double dtT, double EJ, double Xcm2) :
    m_dtL(dtL),
	m_dtT(dtT),
	m_EJ(EJ),
	m_Xcm2(Xcm2)
	{}

double Pulse::calc_omega0() {
	// Calculate angular frequency using the wavelength
	return 2 * std::numbers::pi * PhysicalConstants::c0 / (1e6 * m_wavelength); // 1/fs
}

void Pulse::calc_limits(int nt, double dw) {
	m_lam1 = 2 * std::numbers::pi * PhysicalConstants::c0 * 1e-6 / ( m_omega0 + nt / 2 * dw);
	m_lam2 = 2 * std::numbers::pi * PhysicalConstants::c0 * 1e-6 / ( m_omega0 - nt / 2 * dw);	   
}

