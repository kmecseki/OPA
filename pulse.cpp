#include <numbers>
#include <cmath>
#include <iostream>
#include "pulse.h"
#include "phys_constants.h"

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
	lambda = lambda/1e3; // convert to microns
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

