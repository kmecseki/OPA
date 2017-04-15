#pragma once
#ifndef PULSE_H
#define PULSE_H

/*=====================================================================*/
// Create crystals
class Crystal {
	public:
		Crystal(int, int, int, double, double);
		int m_cType;
		int m_stage;
		int m_noStep;
		double m_thdeg;
		double m_cLength;
		double m_xeff;
		double m_dtPumpL;
		double m_dtPumpT;
		double m_dzcm;

	private:
		double calc_xeff();
};
    
    /*=====================================================================*/
	// Pulse class

class Pulse {
	public:
		double m_dtL;
		double m_dtT;
		double m_EJ;
		double m_Xcm2;
		double m_freq;
		double m_omega0;
		double m_wavelength;
		Pulse(double, double, double, double);
		double calc_freq(double wavelength_nm);
		double calc_omega0();
};

#endif // PULSE_H