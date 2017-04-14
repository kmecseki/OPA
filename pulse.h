/*=====================================================================*/
// Create crystals
class Crystal {
	public:
		int m_cType, m_stage, m_noStep;
		double m_thdeg, m_xeff, m_cLength;
		double m_dtPumpL, m_dtPumpT;
		Crystal(int ctype, double thdeg, int stage, int noStep, int cryslth);

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
		Pulse(double dtL, double dtT, double EJ, double Xcm2);
};