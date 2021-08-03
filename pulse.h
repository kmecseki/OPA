#pragma once
#ifndef PULSE_H
#define PULSE_H
#include <complex>
#include <vector>

class Pulse;
/*=====================================================================*/
// Group chirp data
struct chirp {
	double chpSig;
	double chpPum;
	double chpIdl;
	double chpSig23;
	double chpSig2;
	double chpPum2;
	double chpIdl2;
	double chpSig223;
	double chpSigL;
	double chpPumL;
	double chpIdlL;
	double chpSigNL;
	double chpPumNL;
	double chpIdlNL;
};

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
		void calc_k(double, double, double);
		double calcPumCo(double);
		void calc_PM(double, double, bool);
		void setPhaseVel(const double, double, std::vector<double>, std::vector<double>, std::vector<double>, double*, double*, double*);
		void makePhaseRelative(const int, const double, const double, const double, const double, Pulse&, Pulse&, Pulse&);
		void OPA(Pulse &Pum, Pulse &Sig, Pulse &Idl, double dtps, int chirpType);
	private:
		double calc_xeff();
		double nlindx(int);
};
    
/*=====================================================================*/
// Pulse class
class Pulse {
	public:
		int m_prof;
		int m_nt;
		double m_dtL;
		double m_dtT;
		double m_EJ;
		double m_Xcm2;
		double m_omega0;
		double m_wavelength;
		double m_nOrd;
		double m_xOrd;
		double m_lam1, m_lam2;
		double m_tc;
		double m_kv;
		double m_alp;
		double m_fw;
		std::vector<std::complex<double>> m_cPhij;
		std::vector<double> m_Phij;
		std::vector<double> m_lambdaj;
		std::vector<double> m_absTP;
		std::vector<std::complex<double>> m_ctimeProf;
		Pulse(double, double, double, double, double, int, int);
		double calc_omega0();
		void calc_limits(int, double);
		void GenProfile(Crystal&, double, double, int, double, double, double, double, double);
		double rInt(double dtps);
		double cInt(std::vector<std::complex<double>>, double, double);
		void spectrum(const char *ofname);
		int disperse();
		void nlshift(int, double, double);	
		// TODO: Move private functions over to private (check)
	private:
		int chirper_norm(std::vector<std::complex<double>>&, int, double, double, double);
		int chirper_direct(std::vector<std::complex<double>>&, double, double, double);
};

#endif // PULSE_H