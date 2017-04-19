#pragma once
#ifndef PULSE_H
#define PULSE_H
#include <complex>
#include <vector>

/*=====================================================================*/
// Create crystals
class Crystal {
	public:
		// TODO: check access req.
		Crystal(int, int, int, double, double, double, bool);
		bool m_firstStage;
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
		double m_xcm2;
		std::complex<double> m_pmism;
		std::vector<std::complex<double>> m_cPum, m_cSig, m_cIdl;
		double calcRefInd(double, int);
		double calc_k(double, double, double);
		double calcPumCo(double);
		void calc_PM(double, double, bool);
		void setPhaseVel(const double, double, std::vector<double>, std::vector<double>, std::vector<double>, double*, double*, double*);
		void makePhaseRelative(const int, const double, const double, const double, const double, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::complex<double>>, std::vector<std::complex<double>>, std::vector<std::complex<double>>);
	private:
		double calc_xeff();
};
    
/*=====================================================================*/
// Pulse class
class Pulse {
	public:
		int m_prof;
		double m_dtL;
		double m_dtT;
		double m_EJ;
		double m_Xcm2;
		double m_omega0;
		double m_wavelength;
		double m_nOrd;
		double m_xOrd;
		double m_lam1, m_lam2;
		std::vector<double> m_lambdaj;
		Pulse(double, double, double, double, int);
		double calc_omega0();
		void calc_limits(int, double);
		void GenProfile(Crystal&);
		

};

#endif // PULSE_H