
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<memory>
//#include<cstdio>

#include "pulse.h"
#include "phys_constants.h"





int main(int argc, char *argv[]) {
	
    std::string version = "6.0";
    std::cout << "OPA program version " << version << std::endl;
	
	// Helper constants
	const double small = 1e-20;
    const int ndim = 524288; // 2**19
    
    /*=====================================================================*/
    // Checking argument for parameter input file. Fallback to "input.txt".
    const char *fname; 
    if(argc==1){
		std::cout << "No input argument, using default 'input.txt'" << std::endl;
		fname = "input.txt";
	}
	else {
		std::cout << "Input file given in argument -" << argv[1] << "- is used." << std::endl;
		fname = argv[1];
	}

    // Read parameters from input file.
    std::cout << "Reading parameters from input file.." << std::endl;
    std::vector<double> parsed; 
    std::ifstream filestr;
    filestr.open(fname, std::fstream::in);
	if(filestr.is_open()) {
        char delim = '/';
		double read_number;
        std::string line;
		while(!filestr.eof()) {
			std::getline(filestr, line, delim);
			std::stringstream input(line);
			while(input >> read_number) {
				parsed.push_back(read_number);
			}
		}
	}
	else {
		std::cout << "Error opening input file " << fname << "!" << std::endl;
        exit(0);
	}
	filestr.close();
    
    /*=====================================================================*/
	// Read parameters from pars vector.
	int o=0;
	long long unsigned int noParm_def = 66;
	std::cout << "Reading parameters from vector: " << parsed.size() << std::endl;
	if(parsed.size() != noParm_def) {
		std::cout << "Wrong number of parameters in file." << std::endl;
		exit(0);
	}
	int mode = (int) parsed.at(o++);
	int nStage = (int) parsed.at(o++);
	int cryst1 = (int) parsed.at(o++);
	int cryst2 = (int) parsed.at(o++);
	int cryst3 = (int) parsed.at(o++);
	double crysLth1 = parsed.at(o++);
	double crysLth2 = parsed.at(o++);
	double crysLth3 = parsed.at(o++);
	double thdeg1 = parsed.at(o++);
	double thdeg2 = parsed.at(o++);
	double thdeg3 = parsed.at(o++);
	double nColl_deg = parsed.at(o++);
	double ppm = parsed.at(o++);
	double pLambda_nm = parsed.at(o++);
	double sLambda_nm = parsed.at(o++);
	double pEJ1 = parsed.at(o++);
	double pEJ2 = parsed.at(o++);
	double pEJ3 = parsed.at(o++);
	double sEJ1 = parsed.at(o++);
	double sigTra12 = parsed.at(o++);
	double sigTra23 = parsed.at(o++);
	double iEJ1 = parsed.at(o++);
	int sProf = (int) parsed.at(o++);
	int pProf = (int) parsed.at(o++);
	int iProf = (int) parsed.at(o++);
	double tlead = parsed.at(o++);
	double dtps = parsed.at(o++);
	int nt = (int) parsed.at(o++);
	int noStep1 = (int) parsed.at(o++);
	int noStep2 = (int) parsed.at(o++);
	int noStep3 = (int) parsed.at(o++);
	double dtPumpL1 = parsed.at(o++);
	double dtPumpL2 = parsed.at(o++);
	double dtPumpL3 = parsed.at(o++);
	double dtPumpT1 = parsed.at(o++);
	double dtPumpT2 = parsed.at(o++);
	double dtPumpT3 = parsed.at(o++);
	double dtSigL1 = parsed.at(o++);
	double dtIdlL1 = parsed.at(o++);
	double dtSigT1 = parsed.at(o++);
	double dtIdlT1 = parsed.at(o++);
	double Xcm2_1 = parsed.at(o++);
	double Xcm2_2 = parsed.at(o++);
	double Xcm2_3 = parsed.at(o++);
	int frame = (int) parsed.at(o++);
	int chirpType = (int) parsed.at(o++);
	double chpSig = parsed.at(o++);
	double chpPum = parsed.at(o++);
	double chpIdl = parsed.at(o++);
	double chpSig23 = parsed.at(o++);
	double chpSig2 = parsed.at(o++);
	double chpPum2 = parsed.at(o++);
	double chpIdl2 = parsed.at(o++);
	double chpSig223 = parsed.at(o++);
	double chpSigL = parsed.at(o++);
	double chpPumL = parsed.at(o++);
	double chpIdlL = parsed.at(o++);
	double chpSigNL = parsed.at(o++);
	double chpPumNL = parsed.at(o++);
	double chpIdlNL = parsed.at(o++);
	double tcPum = parsed.at(o++);
	double tcSig = parsed.at(o++);
	double tcIdl = parsed.at(o++);
	double phiP = parsed.at(o++);
	double phiS = parsed.at(o++);
	double phiI = parsed.at(o++);

    /*=====================================================================*/
    // Construct time array.
    std::vector<double> t(nt);
	for (int i=0; i<nt ;i++) {
		t[i] = i * dtps - tlead;
	}

	// Time window
		double tWin = dtps * nt;


	/*=====================================================================*/
	// Set up stage parameters.

	auto stage1 = std::make_unique<Crystal>(cryst1, thdeg1, 1, noStep1, crysLth1);
    auto stage2 = std::make_unique<Crystal>(cryst2, thdeg2, 2, noStep2, crysLth2);
    auto stage3 = std::make_unique<Crystal>(cryst3, thdeg3, 3, noStep3, crysLth3);

	/*=====================================================================*/
	// Create pulses
	auto Pump1 = std::make_unique<Pulse>(dtPumpL1, dtPumpT1, pEJ1, Xcm2_1);
	auto Pump2 = std::make_unique<Pulse>(dtPumpL2, dtPumpT2, pEJ2, Xcm2_2);
	auto Pump3 = std::make_unique<Pulse>(dtPumpL3, dtPumpT3, pEJ3, Xcm2_3);
	auto Signal1 = std::make_unique<Pulse>(dtSigL1, dtSigT1, sEJ1, Xcm2_1);
	auto Idler1 = std::make_unique<Pulse>(dtIdlL1, dtIdlT1, iEJ1, Xcm2_1);

	/*=====================================================================*/
	// Calculate idler wavelength and angular frequencies.
	
	Pump1->m_wavelength = pLambda_nm;
	Pump1->m_freq = Pump1->calc_freq(Pump1->m_wavelength);
	Pump2->m_freq = Pump1->m_freq;
	Pump3->m_freq = Pump1->m_freq;
	Signal1->m_wavelength = sLambda_nm;
	Signal1->m_freq = Signal1->calc_freq(Signal1->m_wavelength);
	Idler1->m_freq = Pump1->m_freq - Signal1->m_freq;

	Pump1->m_omega0 = Pump1->calc_omega0();
	Pump2->m_omega0 = Pump1->calc_omega0();
	Pump3->m_omega0 = Pump1->calc_omega0();
	Signal1->m_omega0 = Signal1->calc_omega0();
	Idler1->m_omega0 = Idler1->calc_omega0();
	
	// This is to prevent numerical errors due to large numbers (>1e30+)
	Idler1->m_wavelength = (PhysicalConstants::c0 * 1e6) / (Idler1->m_freq) * 1e3;








}

