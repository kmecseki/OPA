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
		double m_nOrdSig;
		double m_xOrdSig;
		double m_nOrdIdl;
		double m_xOrdIdl;
		double m_nOrdPum;
		double m_xOrdPum;
		double calcRefInd(double, int);

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
		double m_omega0;
		double m_wavelength;
		double m_nOrd;
		double m_xOrd;
		Pulse(double, double, double, double);
		double calc_omega0();

};

#endif // PULSE_H