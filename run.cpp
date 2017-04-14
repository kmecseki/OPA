
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<memory>
//#include<cstdio>

#include "pulse.h"





int main(int argc, char *argv[]) {
	
    std::string version = "6.0";
    std::cout << "OPA program version " << version << std::endl;

	// Physical constants
	const double eps0 = 8.854e-12; // 8.85418782e-12;
	const double hpl = 6.6256e-34; // 6.626068e-34;
	const double c = 2.997925000000e8;// 299792458;
	const double hce0 = 0.5 * c * eps0;
	
	// Helper constants
	const double tPi = 2 * M_PI; // 2 * 4 * atan(1);
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
	int noParm_def = 66;
	std::cout << "Reading parameters from vector: " << parsed.size() << std::endl;
	if(parsed.size()!=noParm_def) {
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




	/*=====================================================================*/
	// Set up stage parameters.
	Crystal* stage1 = new Crystal(cryst1, thdeg1, 1, noStep1, crysLth1);
	Crystal* stage2 = new Crystal(cryst2, thdeg2, 2, noStep2, crysLth2);
	Crystal* stage3 = new Crystal(cryst3, thdeg3, 3, noStep3, crysLth3);

	/*=====================================================================*/
	// Create pulses
	Pulse* Pump1 = new Pulse(dtPumpL1, dtPumpT1, pEJ1, Xcm2_1);
	Pulse* Signal1 = new Pulse(dtSigL1, dtSigT1, sEJ1, Xcm2_1);
	Pulse* Idler1 = new Pulse(dtIdlL1, dtIdlT1, iEJ1, Xcm2_1);
	
	Pulse* Pump2 = new Pulse(dtPumpL2, dtPumpT2, pEJ2, Xcm2_2);

	Pulse* Pump3 = new Pulse(dtPumpL3, dtPumpT3, pEJ3, Xcm2_3);





	






}

