#include <numbers>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <complex>
//#include<cstdio>

#include "pulse.h"
#include "phys_constants.h"
#include "utils.h"

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
	int ppm = parsed.at(o++);
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

	auto stage1 = std::make_unique<Crystal>(cryst1, thdeg1, 1, noStep1, crysLth1, Xcm2_1, true);
    auto stage2 = std::make_unique<Crystal>(cryst2, thdeg2, 2, noStep2, crysLth2, Xcm2_2, false);
    auto stage3 = std::make_unique<Crystal>(cryst3, thdeg3, 3, noStep3, crysLth3, Xcm2_3, false);

	/*=====================================================================*/
	// Create pulses
	auto Pump1 = std::make_unique<Pulse>(dtPumpL1, dtPumpT1, pEJ1, Xcm2_1, pProf);
	auto Pump2 = std::make_unique<Pulse>(dtPumpL2, dtPumpT2, pEJ2, Xcm2_2, pProf);
	auto Pump3 = std::make_unique<Pulse>(dtPumpL3, dtPumpT3, pEJ3, Xcm2_3, pProf);
	auto Signal1 = std::make_unique<Pulse>(dtSigL1, dtSigT1, sEJ1, Xcm2_1, sProf);
	auto Idler1 = std::make_unique<Pulse>(dtIdlL1, dtIdlT1, iEJ1, Xcm2_1, iProf);

	/*=====================================================================*/
	// Calculate idler wavelength and angular frequencies.
	
	Pump1->m_wavelength = pLambda_nm;
	Signal1->m_wavelength = sLambda_nm;

	Pump1->m_omega0 = Pump1->calc_omega0();
	Pump2->m_omega0 = Pump1->calc_omega0();
	Pump3->m_omega0 = Pump1->calc_omega0();
	Signal1->m_omega0 = Signal1->calc_omega0();
	Idler1->m_omega0 = Pump1->m_omega0 - Signal1->m_omega0;
	
	// This is to prevent numerical errors due to large numbers (>1e30+)
	Idler1->m_wavelength = (2 * std::numbers::pi * PhysicalConstants::c0) / (Idler1->m_omega0 * 1e6);
	std::cout << "Idler wavelength: " << Idler1->m_wavelength << " nm." << std::endl;
	stage1->m_sLam = Signal1->m_wavelength;
	stage1->m_iLam = Idler1->m_wavelength;
	stage1->m_pLam = Pump1->m_wavelength;

	/*=====================================================================*/
	// Calculate refractive indices
	// Note: Only the pump travels in the non-ordinary plane for our PM case.
	stage1->m_nOrdSig = stage1->calcRefInd(Signal1->m_wavelength, 3);
	stage1->m_nOrdIdl = stage1->calcRefInd(Idler1->m_wavelength, 3);
	stage1->m_nOrdPum = stage1->calcRefInd(Pump1->m_wavelength, 2);
	stage1->m_xOrdPum = stage1->calcRefInd(Pump1->m_wavelength, 1);
	
	/*=====================================================================*/
	// Calculate phase matching
	std::cout << "\tn_o Signal\t" << "n_o Idler\t" << "n_o Pump\t" << "n_e Pump" << std::endl;
	std::cout << "\t" << stage1->m_nOrdSig << "\t\t" << stage1->m_nOrdIdl << "\t\t" << stage1->m_nOrdPum << "\t\t" << stage1->m_xOrdPum << std::endl;

	// Angular wavenumbers
	stage1->calc_k(Pump1->m_wavelength, Signal1->m_wavelength, Idler1->m_wavelength);
	//stage2->calc_k(Pump2->m_wavelength, Signal2->m_wavelength, Idler2->m_wavelength); 
	//stage3->calc_k(Pump3->m_wavelength, Signal3->m_wavelength, Idler3->m_wavelength); 

	// Effective angular wavenumbers
	// TODO: function for stages
	double gamma_rad1 = std::asin(stage1->m_kSig * std::sin(deg2rad(nColl_deg))/stage1->m_kIdl);
	stage1->m_coePum = stage1->calcPumCo(stage1->m_pLam);
	stage1->m_knSig = stage1->m_kSig * std::cos(deg2rad(nColl_deg)); // in 1/micron
	stage1->m_knIdl = stage1->m_kIdl * cos(gamma_rad1);
	stage1->m_knPum = stage1->m_kPum * stage1->m_coePum;

	double gamma_rad2 = std::asin(stage2->m_kSig * std::sin(deg2rad(nColl_deg))/stage2->m_kIdl);
	stage2->m_coePum = stage2->calcPumCo(stage1->m_pLam);
	stage2->m_knSig = stage2->m_kSig * std::cos(deg2rad(nColl_deg));
	stage2->m_knIdl = stage2->m_kIdl * cos(gamma_rad2);
	stage2->m_knPum = stage2->m_kPum * stage2->m_coePum;	

	double gamma_rad3 = std::asin(stage3->m_kSig * std::sin(deg2rad(nColl_deg))/stage3->m_kIdl);
	stage3->m_coePum = stage3->calcPumCo(stage1->m_pLam);
	stage3->m_knSig = stage3->m_kSig * std::cos(deg2rad(nColl_deg));
	stage3->m_knIdl = stage3->m_kIdl * cos(gamma_rad3);
	stage3->m_knPum = stage3->m_kPum * stage3->m_coePum;

	// Calculate phase matching angle
	// Perfect phase match?
	bool PerfPM = static_cast<bool>(ppm);
	stage1->calc_PM(nColl_deg, gamma_rad1, PerfPM);
	//stage2->calc_PM(nColl_deg, gamma_rad2, PerfPM);
	//stage3->calc_PM(nColl_deg, gamma_rad3, PerfPM);

	// Wavelength scan limits
	double dw = 2 * std::numbers::pi * 1e-3 / tWin; // 1/fs
	Signal1->calc_limits(nt, dw);
	Idler1->calc_limits(nt, dw);
	Pump1->calc_limits(nt, dw);
	std::cout << "\tSignal\t" << "\t\tPump\t" << "\t\tIdler\t" << std::endl;
	std::cout << "Scan wavelengths:" << std::endl;
	std::cout << "   " << Signal1->m_lam1 << " " << Signal1->m_lam2 << "\t  " << Pump1->m_lam1 << " " << Pump1->m_lam2 << "\t  " << Idler1->m_lam1 << " " << Idler1->m_lam2 << "\t nm " << std::endl;

	// Calculating angular wavenumbers and phases for all the wavelengths in this range.

	double dk_grid = 1e9 * dw / PhysicalConstants::c0; // in microns
	double kvPumj, kvSigj, kvIdlj;
	Signal1->m_lambdaj = std::vector<double>(nt);
	Pump1->m_lambdaj = std::vector<double>(nt);
	Idler1->m_lambdaj = std::vector<double>(nt);
	std::vector<double> nPumj(nt);
	std::vector<double> nSigj(nt);
	std::vector<double> nIdlj(nt);
	std::vector<double> phiPumj(nt);
	std::vector<double> phiSigj(nt);
	std::vector<double> phiIdlj(nt);

	for (int j=0; j<nt; j++) {
		kvPumj = 2 * std::numbers::pi * 1e3 / Pump1->m_lam2 + j * dk_grid; // microns
		kvSigj = 2 * std::numbers::pi * 1e3 / Signal1->m_lam2 + j * dk_grid;
		kvIdlj = 2 * std::numbers::pi * 1e3 / Idler1->m_lam1 - j * dk_grid; // FIX: it fixes idler issue(?)
		if (kvIdlj<0) kvIdlj = 0.1;
		Pump1->m_lambdaj[j] = 2 * std::numbers::pi * 1e3 / kvPumj; // nm
		Signal1->m_lambdaj[j] = 2 * std::numbers::pi * 1e3 / kvSigj;
		Idler1->m_lambdaj[j] = 2 * std::numbers::pi * 1e3 / kvIdlj; // going backwards
		nPumj[j] = stage1->calcPumCo(Pump1->m_lambdaj[j]);
		nPumj[j] = nPumj[j] * stage1->calcRefInd(Pump1->m_lambdaj[j], 2);
		nSigj[j] = stage1->calcRefInd(Signal1->m_lambdaj[j], 3);
		nIdlj[j] = stage1->calcRefInd(Idler1->m_lambdaj[j], 3); // going backwards
		phiPumj[j] = -kvPumj * stage1->m_dzcm * nPumj[j];
		phiSigj[j] = -kvSigj * stage1->m_dzcm * nSigj[j] * std::cos(deg2rad(nColl_deg));
		double gamsij = kvSigj * nSigj[j] * std::sin(deg2rad(nColl_deg)) / (kvIdlj * nIdlj[j]);
		if(abs(gamsij)<=1.0)
			phiIdlj[j] = -kvIdlj * stage1->m_dzcm * nIdlj[j] * std::cos(std::asin(gamsij));
		else phiIdlj[j] = 0;
	}

	// TODO: New idler handling, uncomment this section for original idler handling
	std::vector<double> temp(nt);
	for (int j=0; j<nt; j++) {
		temp[j] = phiIdlj[nt-1-j];
	}
	for (int j=0; j<nt; j++) {
		phiIdlj[j] = temp[j];
	}
	for (int j=0; j<nt; j++) {
		temp[j] = nIdlj[nt-1-j];
	}
	for (int j=0; j<nt; j++) {
		nIdlj[j] = temp[j];
	}

	/*=====================================================================*/
	double gtdPum, gtdSig, gtdIdl;
	stage1->setPhaseVel(dw, nColl_deg, nPumj, nSigj, nIdlj, &gtdPum, &gtdSig, &gtdIdl);
	std::vector<std::complex<double>> cPhiPumj, cPhiSigj, cPhiIdlj;
	stage1->makePhaseRelative(frame, dw, gtdSig, gtdPum, gtdIdl, phiPumj, phiSigj, phiIdlj, cPhiPumj, cPhiSigj, cPhiIdlj);

	/*=====================================================================*/
	// Create pulses
	std::cout <<  "\n\n*********************************************";
	std::cout << "\n\tCreating pulses\n";
	std::cout << "*********************************************" << std::endl;
	// PUMP
	if (Pump1->m_EJ!=0) {
		std::cout << "\nPump:" << std::endl;
		if (Pump1->m_prof<1 || Pump1->m_prof>2) {
			std::cout << "Error: Profile reading from file can only be used for signal pulse in the current version" << std::endl;
		}
		double fwp = 0.5 * PhysicalConstants::c0 * PhysicalConstants::eps0 * (stage1->m_nOrdPum * stage1->m_coePum) * stage1->m_xcm2; // Check units if sg is wrong
		std::vector<std::complex<double>> timeProfPum(nt);

		Pump1->GenProfile(timeProfPum);
			
			
			
			
			, dtPumpL, dtPumpT, Xcm2, pEJ, fwp, tcPum, chpPum, chpPum2, chpPumL, chpPumNL);
		
		
		
		
		
		
		}
		else {
			for (j=0;j<nt;j++) {
				timeProfPum[j] = 0;
			}
		}
		//cout << pow(abs(timeProfPum[nt/2]),2) << " " << fwp << " " << dtps << endl;
		//exit(0);

}

