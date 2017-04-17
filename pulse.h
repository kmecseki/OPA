#pragma once
#ifndef PULSE_H
#define PULSE_H
#include <complex>

/*=====================================================================*/
// Create crystals
class Crystal {
	public:
		// TODO: check access req.
		Crystal(int, int, int, double, double);
		int m_cType;
		int m_stage;
		int m_noStep;
		double m_thdeg;
		double m_cLength;
		double m_xeff;
		double m_dtPumpL, m_dtPumpT;
		double m_dzcm;
		double m_nOrdSig, m_xOrdSig;
		double m_nOrdIdl, m_xOrdIdl;
		double m_nOrdPum, m_xOrdPum;
		double m_kSig, m_kIdl, m_kPum, m_coePum;
		double m_knSig, m_knPum, m_knIdl;
		double m_pLam, m_sLam, m_iLam;
		double m_pAng; // PM angle
		double m_dk, m_cohL;
		std::complex<double> m_pmism;
		double calcRefInd(double, int);
		double calc_k(double, double, double);
		double calcPumCo();
		void calc_PM(double, double, bool);

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