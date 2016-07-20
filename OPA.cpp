/*=====================================================================*/
/*
	OPA program C++ version 0.7

	Calculates blabla using blabla
	Units.
	Based on G New's code

	K. Mecseki <k.mecseki09@imperial.ac.uk>
	June 2012.

	Compilation: g++ -lgsl -lgslcblas -lfftw3 OPAv7.cpp

*/
/*=====================================================================*/

#include "main.h"

using namespace std;
using std::vector;

// Global variables
string version;
int warned=0;
int mode, nStage;
int cryst1, cryst2, cryst3;
int sProf, pProf, iProf;
int nt, noStep1, noStep2, noStep3;
int frame, chirpType;
double crysLth1, crysLth2, crysLth3;
double thdeg1, thdeg2, thdeg3;
double nColl_deg, nColl_deg1, nColl_deg2, nColl_deg3, ppm;
double pEJ1, pEJ2, pEJ3, sEJ1;
double sigTra12, sigTra23, iEJ1;
double tlead, dtps;
double dtPumpL1, dtPumpL2, dtPumpL3, dtPumpT1, dtPumpT2, dtPumpT3;
double dtSigL1, dtIdlL1, dtSigT1, dtIdlT1;
double Xcm2_1, Xcm2_2, Xcm2_3;
double chpSig, chpPum, chpIdl, chpSig23, chpSig2, chpPum2, chpIdl2, chpSig223;
double chpSigL, chpPumL, chpIdlL, chpSigNL, chpPumNL, chpIdlNL;
double tcPum, tcSig, tcIdl, phiP, phiS, phiI;
double tWin, dw;
double *t;
const double tPi = 2*4*atan(1);
const double c = 2.997925000000e8;//299792458;
complex<double> c1(1,0),ci(0,1);
extern Pulse signal, pump, idler;

void init(int, char*[]);
int errorhl(int);
double deg2rad(double);
double rad2deg(double);

/*=====================================================================*/

int main(int argc, char *argv[])
{

// Set up stage parameters
        int cStage, warned;
        int cCryst=0, noStep;
        double thdeg, nColl_deg;
        double xeff;
        double crysLth;
        double dtPumpL, dtPumpT, dtSigL, dtSigT, dtIdlL, dtIdlT;
        double pEJ, sEJ, iEJ, Xcm2;
        double gamma_rad, tanThetaSq, coePum;
        double pa, pb, pc, pd, pu, pw, puw, pAng;
        double cohLength, dk, dzcm;
        double pLam1, pLam2, sLam1, sLam2, iLam1, iLam2;
        complex<double> cPhMisM;

        double xeffc(int,int);
// Initialisation
        init(argc, argv);
// Set up stage parameters.
        for (cStage=1; cStage<=nStage; cStage++) {
                warned = 0;
                if (cStage==1) {
                        cCryst = cryst1;
                        thdeg = thdeg1;
                        nColl_deg = nColl_deg1;
                        xeff = xeffc(cCryst,cStage);
                        noStep = noStep1;
                        crysLth = crysLth1;
                        dtPumpL = dtPumpL1;
                        dtPumpT = dtPumpT1;
                        dtSigL = dtSigL1;
                        dtSigT = dtSigT1;
                        dtIdlL = dtIdlL1;
                        dtIdlT = dtIdlT1;
                        pEJ = pEJ1;
                        sEJ = sEJ1;
                        iEJ = iEJ1;
                        Xcm2 = Xcm2_1;
                } else if(cStage==2) {
                        cCryst = cryst2;
                        thdeg=thdeg2;
                        nColl_deg = nColl_deg2;
                        xeff=xeffc(cCryst,cStage);
                        noStep = noStep2;
                        crysLth = crysLth2;
                        dtPumpL = dtPumpL2;
                        dtPumpT = dtPumpT2;
                        pEJ = pEJ2;
                        Xcm2 = Xcm2_2;
                } else if(cStage==3) {
                        cCryst = cryst3;
                        thdeg = thdeg3;
                        nColl_deg = nColl_deg3;
                        xeff = xeffc(cCryst,cStage);
                        noStep = noStep3;
                        crysLth = crysLth3;
                        dtPumpL = dtPumpL3;
                        dtPumpT = dtPumpT3;
                        pEJ = pEJ3;
                        Xcm2 = Xcm2_3;
                }
                // Calculate refractive indices
                signal.no0 = signal.SetRefind(cCryst,3);
                /*idler.no0 = idler.SetRefind(cCryst,3);
                pump.no0 = pump.SetRefind(cCryst,2);
                pump.ne0 = pump.SetRefind(cCryst,1);
                // Calculate phase matching
                // Angular wavenumbers
                signal.k0 = signal.no0 * tPi*1e3/signal.lambda0; // in 1/micron
                idler.k0 = idler.no0 * tPi*1e3/idler.lambda0;
                pump.k0 = pump.k0 * tPi*1e3/pump.lambda0;
                // Effective angular wavenumbers
                gamma_rad = asin(signal.k0 * sin(deg2rad(nColl_deg))/idler.k0);
                tanThetaSq = pow(tan(deg2rad(thdeg)),2);
                coePum = sqrt((1 + tanThetaSq)/(1 + pow((pump.no0/pump.ne0),2)*tanThetaSq));
                signal.kn0 = signal.k0 * cos(deg2rad(nColl_deg)); // in 1/micron
                idler.kn0 = idler.k0 * cos(gamma_rad);
                pump.kn0 = pump.k0 * coePum;
                // Calculate phase matching angle
                pa = cos(deg2rad(nColl_deg))*signal.no0/signal.lambda0; // nm
                pb = cos(gamma_rad)*idler.no0/idler.lambda0;
                pc = pump.no0/pump.lambda0;
                pd = pump.ne0/pump.lambda0;
                pu = pow(((pa+pb)/pc),2);
                pw = pow(((pa+pb)/pd),2);
                puw = (1-pu)/(pw-1);
                if (puw<0) errorhl(4);
                pAng = rad2deg(atan(sqrt(puw)));
                // Phase mismatch and coherence length
                if (ppm == 1) {
                        dk = 0;
                        cohLength = 0;
                } else {
                        //dk = kPum-kSig-kIdl;
                        dk = (pump.kn0-signal.kn0-idler.kn0)*1e4; // in 1/cm
                        cohLength = abs(tPi/(2*dk));  // in micron
                }
                cPhMisM = exp(ci*dk*dzcm/2e4);
                // Wavelength scan limits
                tWin = dtps*nt;
                dw = tPi*1e-3/tWin; // 1/fs
                pLam1 = tPi*c*1e-6/(pump.omega0+nt/2*dw);  // nm
                pLam2 = tPi*c*1e-6/(pump.omega0-nt/2*dw);
                sLam1 = tPi*c*1e-6/(signal.omega0+nt/2*dw);
                sLam2 = tPi*c*1e-6/(signal.omega0-nt/2*dw);
                iLam1 = tPi*c*1e-6/(idler.omega0+nt/2*dw);
                iLam2 = tPi*c*1e-6/(idler.omega0-nt/2*dw);
                cout << "\tSignal\t" << "\t\tPump\t" << "\t\tIdler\t" << endl;
                cout << "Scan wavelengths:" << endl;
                cout << "   " << sLam1 << " " << sLam2 << "\t  " << pLam1 << " " << pLam2 << "\t  " << iLam1 << " " << iLam2 << "\t nm " << endl;
               */ /*=====================================================================*/
                // Calculating angular wavenumbers and phases for all the wavelengths in this range.
		/*dk_grid = 1e9*dw/c; // in microns
		for (j=0;j<nt;j++) {
			kvPumj = tPi*1e3/pLam2 + j*dk_grid; // microns
			kvSigj = tPi*1e3/sLam2 + j*dk_grid;
			kvIdlj = tPi*1e3/iLam1 - j*dk_grid; // FIX: hopefully it fixes idler issue
			if (kvIdlj<0) kvIdlj=0.1;
			pLambdaj[j] = tPi*1e3/kvPumj; // nm
			sLambdaj[j] = tPi*1e3/kvSigj;
			iLambdaj[j] = tPi*1e3/kvIdlj; // going backwards
			nPumj[j] = CalcPumCo(pLambdaj[j]);
			nPumj[j] = nPumj[j] * calcRefInd(cCryst, pLambdaj[j], 2);
			nSigj[j] = calcRefInd(cCryst, sLambdaj[j], 3);
			nIdlj[j] = calcRefInd(cCryst, iLambdaj[j], 3); // going backwards
			phiPumj[j] = -kvPumj*dzcm*nPumj[j];
			phiSigj[j] = -kvSigj*dzcm*nSigj[j]*cos(deg2rad(nColl_deg));
			gamsij = kvSigj*nSigj[j]*sin(deg2rad(nColl_deg))/(kvIdlj*nIdlj[j]);
			if(abs(gamsij)<=1.0)
				phiIdlj[j] = -kvIdlj*dzcm*nIdlj[j]*cos(asin(gamsij));
			else phiIdlj[j] = 0;
		}*/


        }
}






