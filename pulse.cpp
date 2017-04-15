#define _USE_MATH_DEFINES
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

Pulse::Pulse(double dtL, double dtT, double EJ, double Xcm2) :
    m_dtL(dtL),
	m_dtT(dtT),
	m_EJ(EJ),
	m_Xcm2(Xcm2)
	{}

double Pulse::calc_freq(double wavelength_nm) {
	return PhysicalConstants::c0 * 1e6 / (1e-3 * wavelength_nm);
}

double Pulse::calc_omega0() {
	return 2 * M_PI * Pulse::m_freq / 1e9; // 1/fs
}