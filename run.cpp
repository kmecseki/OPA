
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<memory>
//#include<cstdio>





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
    std::unique_ptr<std::vector<double>> parsed = std::make_unique<std::vector<double>>(); 
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
				parsed->push_back(read_number);
			}
		}
	}
	else {
		std::cout << "Error opening input file " << fname << "!" << std::endl;
        exit(0);
	}
	filestr.close();
    
    /*=====================================================================*/
	
    // Construct time array.
    double<vector> t(nt);
	for (int i=0; i<nt ;i++) {
		t[i] = i * dtps - tlead;
	}
    
    
    
	
	
	
	

}

















    openfile(fname, pars);
	int o=0;
	cout << "Read/Required parameters from file: " << pars->size() << "/" << noParm_def << endl;
	if(pars->size()!=noParm_def) errorhl(2);
	mode = (int) pars->at(o++);
	nStage = (int) pars->at(o++);
	cryst1 = (int) pars->at(o++);
	cryst2 = (int) pars->at(o++);
	cryst3 = (int) pars->at(o++);
	crysLth1 = pars->at(o++);
	crysLth2 = pars->at(o++);
	crysLth3 = pars->at(o++);
	thdeg1 = pars->at(o++);
	thdeg2 = pars->at(o++);
	thdeg3 = pars->at(o++);
	nColl_deg = pars->at(o++);
	ppm = pars->at(o++);
	pLambda_nm = pars->at(o++);
	sLambda_nm = pars->at(o++);
	pEJ1 = pars->at(o++);
	pEJ2 = pars->at(o++);
	pEJ3 = pars->at(o++);
	sEJ1 = pars->at(o++);
	sigTra12 = pars->at(o++);
	sigTra23 = pars->at(o++);
	iEJ1 = pars->at(o++);
	sProf = (int) pars->at(o++);
	pProf = (int) pars->at(o++);
	iProf = (int) pars->at(o++);
	tlead = pars->at(o++);
	dtps = pars->at(o++);
	nt = (int) pars->at(o++);
	noStep1 = (int) pars->at(o++);
	noStep2 = (int) pars->at(o++);
	noStep3 = (int) pars->at(o++);
	dtPumpL1 = pars->at(o++);
	dtPumpL2 = pars->at(o++);
	dtPumpL3 = pars->at(o++);
	dtPumpT1 = pars->at(o++);
	dtPumpT2 = pars->at(o++);
	dtPumpT3 = pars->at(o++);
	dtSigL1 = pars->at(o++);
	dtIdlL1 = pars->at(o++);
	dtSigT1 = pars->at(o++);
	dtIdlT1 = pars->at(o++);
	Xcm2_1 = pars->at(o++);
	Xcm2_2 = pars->at(o++);
	Xcm2_3 = pars->at(o++);
	frame = (int) pars->at(o++);
	chirpType = (int) pars->at(o++);
	chpSig = pars->at(o++);
	chpPum = pars->at(o++);
	chpIdl = pars->at(o++);
	chpSig23 = pars->at(o++);
	chpSig2 = pars->at(o++);
	chpPum2 = pars->at(o++);
	chpIdl2 = pars->at(o++);
	chpSig223 = pars->at(o++);
	chpSigL = pars->at(o++);
	chpPumL = pars->at(o++);
	chpIdlL = pars->at(o++);
	chpSigNL = pars->at(o++);
	chpPumNL = pars->at(o++);
	chpIdlNL = pars->at(o++);
	tcPum = pars->at(o++);
	tcSig = pars->at(o++);
	tcIdl = pars->at(o++);
	phiP = pars->at(o++);
	phiS = pars->at(o++);
	phiI = pars->at(o++);
}


