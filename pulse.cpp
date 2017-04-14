#include <iostream>	
#include "pulse.h"

Crystal::Crystal(int ctype, double thdeg, int stage, int noStep, int cryslth) : 
    m_cType(ctype),
    m_thdeg(thdeg),
	m_stage(stage),
	m_noStep(noStep),
	m_cLength(cryslth) 
    {	
        m_xeff = calc_xeff();
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