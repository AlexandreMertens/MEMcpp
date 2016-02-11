#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "classes/DelphesClasses.h"

#include "src/HelAmps_sm.h"
#include "SubProcesses/P0_Sigma_sm_uux_epvemumvmx/CPPProcess.h"
#include "src/rambo.h"
#include "src/Parameters_sm.h"

#include "TStopwatch.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "cuba.h"

#include "utils.h"
#include "jacobianF.h"

#define M_W 80.419
#define G_W 2.0476

#define SQRT_S 8000 

#define SSTR( x ) dynamic_cast< std::ostringstream & > \
				( std::ostringstream() << std::dec << x ).str()

using namespace LHAPDF;
using namespace std;

int mycount = 0, count_wgt = 0, count_perm=1;

unsigned int setFlags(char verbosity = 0, bool subregion = false, bool retainStateFile = false, unsigned int level = 0, bool smoothing = false, bool takeOnlyGridFromFile = true){
	/* Set option flags for CUBA integrator
	 * 
	 * smoothing only used by Suave
	 * smoothing and takeOnlyGridFromFile only used by Vegas
	 *
	 * verbosity = 0-3
	 * subregion true = only last set of samples is used for final evaluation of integral
	 * smoothing true = smoothe importance function
	 * retainStateFile true => retain state file when integration ends
	 * takeOnlyGridFromFile false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
	 * level determines random-number generator:
	 *	seed = 0 => Sobol quasi-random
	 *	seed != 0 => level is used:
	 *		level = 0 => Mersenne Twister pseudo-random
	 *		level != 0 => Ranlux pseudo-random, level = determines period of generator
	*/

	unsigned int flags = 0;

	unsigned int opt_subregion = 0x04; // bit 2
	unsigned int opt_smoothing = 0x08; // bit 3
	unsigned int opt_retainStateFile = 0x10; // bit 4
	unsigned int opt_takeOnlyGridFromFile = 0x20; // bit 5

	level <<= 8; // bits 8-31
	flags |= level | verbosity; // verbosity: bis 0-1
	if(subregion) flags |= opt_subregion;
	if(!smoothing) flags |= opt_smoothing; // careful true-false inverted
	if(retainStateFile) flags |= opt_retainStateFile;
	if(takeOnlyGridFromFile) flags |= opt_takeOnlyGridFromFile;

	cout << "Integrator flags = ";
	for(int i=31; i>=0; --i){
		bool bit = (flags >> i) & 1;
		cout << bit;
	}
	cout << endl;

	return flags;
}

class MEWeight{
	public:

<<<<<<< HEAD
	MEWeight(const string paramCardPath, const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector met);
=======
  MEWeight(const string paramCardPath, const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met, double q1gen, double q2gen);
	//MEWeight(const string paramCardPath, const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
>>>>>>> 0bcdf71088ef8c00a78a95cc83e15f7c40f138a3
	inline double ComputePdf(const int pid, const double x, const double q2);

	inline TLorentzVector GetP3(void) const { return p3; }
	inline TLorentzVector GetP4(void) const { return p4; }
	inline TLorentzVector GetMet(void) const { return Met; }

	inline void setProcessMomenta(vector<double*> &p){ process.setMomenta(p); }
	inline void computeMatrixElements(){ process.sigmaKin(); }
	inline const double* const getMatrixElements() const { return process.getMatrixElements(); }

  double q1, q2;

	private:

	TLorentzVector p3, p4, Met;

	CPPProcess process;
	PDF *pdf;
};

MEWeight::MEWeight(const string paramCardPath, const TLorentzVector ep, const TLorentzVector mum,  const TLorentzVector met){
	p3 = ep;
	p4 = mum;
	Met = met;

  q1 = q1gen;
  q2 = q2gen;

	process.initProc(paramCardPath);
	pdf = mkPDF("cteq6l1", 0);
}

double MEWeight::ComputePdf(const int pid, const double x, const double q2){
	// return f(pid,x,q2)
	if(x <= 0 || x >= 1)
		return 0.;
	else
		return pdf->xfxQ2(pid, x, q2)/x;
}


int MEFunct(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs){

	//cout << endl << endl << endl << "########## Starting phase-space point ############" << endl << endl;
	//cout << "Inputs = [" << Xarg[0] << "," << Xarg[1] << "," << Xarg[2] << "," << Xarg[3] << endl;

	MEWeight* myWeight = (MEWeight*) inputs;

	TLorentzVector p1 = myWeight->GetP4();
	TLorentzVector p2 = myWeight->GetP6();

	TLorentzVector p3 = myWeight->GetP3();
	TLorentzVector p4 = myWeight->GetP4();
	TLorentzVector Met = myWeight->GetMet();

	*Value = 0.;

	for(int i=0; i<*nDim; ++i){
		if(Xarg[i] == 1. || Xarg[i] == 0.){
			mycount++;
			return 0;
		}
	}

	// We flatten the Breit-Wigners by doing a change of variable for each resonance separately
	// s = M G tan(y) + M^2
	// jac = M G / cos^2(y)
	// ==> BW(s(y))*jac(y) is flat in the variable y, as BW(s) = 1/((s-M^2)^2 - (GM)^2)
	// Where y = -arctan(M/G) + (pi/2+arctan(M/G))*x_foam (x_foam between 0 and 1 => s between 0 and infinity)

	const double range1 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
	const double y1 = - TMath::ATan(M_W/G_W) + range1 * Xarg[0];
	const double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);
         
//        const double range1 = 4000;
//        const double s13 = 4000 + range1 * Xarg[0];
    
        //cout << " s13 : " << s13 << endl;


	//cout << "y1=" << y1 << ", m13=" << TMath::Sqrt(s13) << endl;

	const double range2 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
	const double y2 = - TMath::ATan(M_W/G_W) + range2 * Xarg[1];
	//const double s24 = M_W * G_W * TMath::Tan(y2) + pow(M_W,2.);
	//const double s24 = SQ(77.8949);
	const double s24 = SQ((p2+p4).M());

        //const double range2 = 4000;
        //const double s24 = 4000 + range2 * Xarg[1];


        //cout << " s24 : " << s24 << endl;


        const double range4 = 2000;
        const double range3 = 10;
        


	const double sqrt_shat = range4*Xarg[2];
	//const double rapidity = 0.5*log((shat+171.0)/(shat-171.0)); //-2+2*Xarg[3];
        const double rapidity = -5+range3*Xarg[3];

        //cout << " shat : " << shat << " rapidity : " << rapidity << endl;

        const double q2 = TMath::Exp(-rapidity)*sqrt_shat/SQRT_S;
        const double q1 = TMath::Exp(rapidity)*sqrt_shat/SQRT_S;

        //const double q1 = range3*Xarg[2];
        //const double q2 = range4*Xarg[3];


        //cout <<" q1, q2 " <<  q1 << " " << q2 << endl;

        //const double pz_tot = SQRT_S/2.0 * (q1-q2);
        //cout << " pz   : " << pz_tot << endl;



        //const double Jac_dsdy = 2/(SQRT_S*SQRT_S)*pow(sqrt_shat*sqrt_shat,-3/2);


        const double Jac_dsdy = 2*sqrt(q1*q2)/SQRT_S;

	//cout << "y2=" << y2 << ", M24=" << TMath::Sqrt(s24) << endl;



        // Starting the weight computation:
        // 1) phase-space point kinematic

	// pb = transverse total momentum of the visible particles
	const TLorentzVector pb = p3+p4+p1+p2;
	//const TLorentzVector pb = p3+p4;
	//const TLorentzVector pb = -Met;

	// P2x = a1 E2 + a2 P2y + a3
	// P2z = b1 E2 + b2 P2y + b3

	const double Qm = SQRT_S*(q1-q2)/2.;
	const double Qp = SQRT_S*(q1+q2)/2.;

        const double ka = -4*p4.Px()*p3.Pz()+4*p3.Px()*p4.Pz();
	const double kb = 2*(p3.Pz()*p4.Px()-p3.Px()*p4.Pz());

        const double b1 = -(1./kb) * (2*p4.E()*p3.Px()-2*p3.E()*p4.Px());
        const double b2 = -(1./kb) * (2*p3.Py()*p4.Px()-2*p3.Px()*p4.Py());
        const double b3 = -(1./kb) * (pow(p4.M(),2)*p3.Px()-2*p3.E()*pb.E()*p4.Px()+pow(p3.M(),2)*p4.Px()+2*p3.Px()*p4.Px()*pb.Px()+2*p3.Py()*p4.Px()*pb.Py()+2*p3.Pz()*p4.Px()*pb.Pz()-2*p3.Pz()*p4.Px()*Qm+2*p3.E()*p4.Px()*Qp-p4.Px()*s13-p3.Px()*s24);

        const double a1 = (1./ka)*(-2*p4.Pz()*((-2)*p3.E())-2*p3.Pz()*2*p4.E());
        const double a2 = (1./ka)*(-4*p3.Py()*p4.Pz()-2*p3.Pz()*(-2*p4.Py()));
        const double a3 = -pb.Px()+(1./ka)*(-4*p4.Px()*pb.Px()*p3.Pz()-4*p3.Py()*pb.Py()*p4.Pz()-2*p4.Pz()*(pow(p3.M(),2)+2*p3.Pz()*(pb.Pz()-Qm)-2*p3.E()*(pb.E()-Qp)-s13)-2*p3.Pz()*(pow(p4.M(),2)-s24)) ;

	// 0 = c1 E2 + c2 P2y + c3 
	// E2 = d1 P2y + d2
	
        const double c1 = 2*(Qp-pb.E())+2*pb.Px()*a1-2*(Qm-pb.Pz())*b1;
        const double c2 = 2*pb.Px()*a2+2*pb.Py()-2*(Qm-pb.Pz())*b2;
        const double c3 = -pow(Qp-pb.E(),2)+pow(pb.Px(),2)+2*pb.Px()*a3+pow(pb.Py(),2)+pow(Qm-pb.Pz(),2)-2*(Qm-pb.Pz())*b3;


	//cout << " c : " << c1 << " " << c2 << " " << c3 << endl;
	
	const double d1 = -c2/c1;
	const double d2 = -c3/c1;

	//cout << " d : " << d1 << " " << d2 << endl; 

	// alpha*P2y^2 + beta*P2y + gamma = 0

        const double alpha = pow(a1,2)*pow(d1,2)+pow(a2,2)+2*a1*a2*d1+1+pow(b1,2)*pow(d1,2)+pow(b2,2)+2*b1*b2*d1-pow(d1,2);
        const double beta  = 2*pow(a1,2)*d1*d2+2*a1*a3*d1+2*a1*a2*d2+2*a2*a3+2*pow(b1,2)*d1*d2+2*b1*b3*d1+2*b1*b2*d2+2*b2*b3-2*d1*d2;
        const double gamma = pow(a1,2)*pow(d2,2)+pow(a3,2)+2*a1*a3*d2+pow(b1,2)*pow(d2,2)+pow(b3,2)+2*b1*b3*d2-pow(d2,2);

	
	// Find P2Y
	vector<double> P2Y;

	solveQuadratic(alpha, beta, gamma, P2Y, false);

	// For each solution of P2Y, find the neutrino 4-momenta p1,p2, find the initial quark momenta,
	// evaluate the matrix element and the jacobian
	
	if(P2Y.size() == 0){
		//cout << "No solutions!" << endl;
		mycount++;
		return 0;
	}

	for(unsigned int i=0; i<P2Y.size(); i++){
		const double P2X = a1*d1*P2Y.at(i)+a1*d2+a2*P2Y.at(i)+a3;
		const double P2Z = b1*d1*P2Y.at(i)+b1*d2+b2*P2Y.at(i)+b3;

                const double P1X = -pb.Px()-P2X;
                const double P1Y = -pb.Py()-P2Y.at(i);
                const double P1Z = Qm-pb.Pz()-P2Z;

		const double P2E = sqrt(pow(P2X,2)+pow(P2Y.at(i),2)+pow(P2Z,2));
                const double P1E = sqrt(pow(P1X,2)+pow(P1Y,2)+pow(P1Z,2));

		//cout << "  W24 mass : " << sqrt(s24);
		//cout << "  W13 mass : " << sqrt(s13) << endl;

		TLorentzVector p1,p2;

		p1.SetPxPyPzE( P1X, P1Y, P1Z, P1E);
		p2.SetPxPyPzE( P2X, P2Y.at(i), P2Z, P2E );

                bool debug = 0;

                if (debug){
                  // Test dsig 
                  p1.SetPxPyPzE(5.9056514367095616      ,   16.011232873614162      ,   441.66223458284003      ,   441.99181638773882);
                  p2.SetPxPyPzE(2.2306333952012740      ,  0.73903578464239128      ,   109.40085102684131      ,   109.42608511972350);
                  p3.SetPxPyPzE(18.174801324009138      ,  -14.048879079584797      ,   34.838891037145132      ,   41.730597111209313);
                  p4.SetPxPyPzE(-26.311086155919973      ,  -2.7013895786717557      ,   1.1557566283843559      ,   26.474639445024639);
                }



		const TLorentzVector p13 = p1 + p3;
		const TLorentzVector p24 = p2 + p4;
	
		const TLorentzVector tot = p1 + p2 + p3 + p4;

		const double ETot = tot.E();
		const double PzTot = tot.Pz();

		const double q1Pz = (PzTot + ETot)/2.;
		const double q2Pz = (PzTot - ETot)/2.;
               

                if (debug) cout << "q1Pz : " << q1Pz << endl;
                if (debug) cout << "q2Pz : " << q2Pz << endl;
 

		//cout << "===> Eext=" << ETot << ", Pzext=" << PzTot << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl << endl;
	
	//	if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.){
			//cout << "Fail!" << endl;
	//		mycount++;
			//continue;
	//		break;
	//	}


		// momentum vector definition
		vector<double*> p_1(1, new double[4]);
		p_1[0][0] = TMath::Abs(q1Pz); p_1[0][1] = 0.0; p_1[0][2] = 0.0; p_1[0][3] = q1Pz;
		p_1.push_back(new double[4]);
		p_1[1][0] = TMath::Abs(q2Pz); p_1[1][1] = 0.0; p_1[1][2] = 0.0; p_1[1][3] = q2Pz;
		p_1.push_back(new double[4]);
		p_1[2][0] = p3.E(); p_1[2][1] = p3.Px(); p_1[2][2] = p3.Py(); p_1[2][3] = p3.Pz();
		p_1.push_back(new double[4]);
		p_1[3][0] = p1.E(); p_1[3][1] = p1.Px(); p_1[3][2] = p1.Py(); p_1[3][3] = p1.Pz();
		p_1.push_back(new double[4]);
		p_1[4][0] = p4.E(); p_1[4][1] = p4.Px(); p_1[4][2] = p4.Py(); p_1[4][3] = p4.Pz();
		p_1.push_back(new double[4]);
		p_1[5][0] = p2.E(); p_1[5][1] = p2.Px(); p_1[5][2] = p2.Py(); p_1[5][3] = p2.Pz();


                vector<double*> p_2(1, new double[4]);
                p_2[0][0] = TMath::Abs(q2Pz); p_2[0][1] = 0.0; p_2[0][2] = 0.0; p_2[0][3] = q2Pz;
                p_2.push_back(new double[4]);
                p_2[1][0] = TMath::Abs(q1Pz); p_2[1][1] = 0.0; p_2[1][2] = 0.0; p_2[1][3] = q1Pz;
                p_2.push_back(new double[4]);
                p_2[2][0] = p3.E(); p_2[2][1] = p3.Px(); p_2[2][2] = p3.Py(); p_2[2][3] = p3.Pz();
                p_2.push_back(new double[4]);
                p_2[3][0] = p1.E(); p_2[3][1] = p1.Px(); p_2[3][2] = p1.Py(); p_2[3][3] = p1.Pz();
                p_2.push_back(new double[4]);
                p_2[4][0] = p4.E(); p_2[4][1] = p4.Px(); p_2[4][2] = p4.Py(); p_2[4][3] = p4.Pz();
                p_2.push_back(new double[4]);
                p_2[5][0] = p2.E(); p_2[5][1] = p2.Px(); p_2[5][2] = p2.Py(); p_2[5][3] = p2.Pz();


                //cout << "Phase space point : " << endl;
                //for (int it = 0; it < 6; it++){
                //  cout << p[it][0] << " " << p[it][1] << " " << p[it][2] << " " << p[it][3] << endl;
                //}

                // Set momenta for this event
                myWeight->setProcessMomenta(p_1);
                // Evaluate matrix element
                myWeight->computeMatrixElements();
                const double* const matrix_elements1 = myWeight->getMatrixElements();
                const double ME1 = matrix_elements1[0];
                if (debug) cout << "ME1 : " << ME1 << endl;


                myWeight->setProcessMomenta(p_2);
                // Evaluate matrix element
                myWeight->computeMatrixElements();
                const double* const matrix_elements2 = myWeight->getMatrixElements();
                const double ME2 = matrix_elements2[0];
                if (debug) cout << "ME2 : " << ME2 << endl;
        
 

		// Compute the Pdfs
		const double Q_2 = pow(91.1880,2); //pow(ETot,2)
		const double pdf1_1 = myWeight->ComputePdf(2,TMath::Abs(q1Pz/(SQRT_S/2.)), Q_2);
		const double pdf1_2 = myWeight->ComputePdf(-2,TMath::Abs(q2Pz/(SQRT_S/2.)), Q_2);

                //cout << "pdfs: " << pdf1_1 << " " << pdf1_2 << " " << pdf1_1*pdf1_2 << endl;

                const double pdf2_1 = myWeight->ComputePdf(2,TMath::Abs(q2Pz/(SQRT_S/2.)), Q_2);
                const double pdf2_2 = myWeight->ComputePdf(-2,TMath::Abs(q1Pz/(SQRT_S/2.)), Q_2);

                //cout << "pdfs: " << pdf2_1 << " " << pdf2_2 << " " << pdf2_1*pdf2_2 << endl;
	
		// Compute flux factor 1/(x1*x2*s)
		const double PhaseSpaceIn = 1.0 / ( TMath::Abs(q1Pz/(SQRT_S/2.)) * TMath::Abs(q2Pz/(SQRT_S/2.0)) * pow(SQRT_S,2)); 

                //const double PhaseSpaceIn = 4.0 / (SQRT_S*SQRT_S);

		// Compute finale Phase Space for observed particles (not concerned by the change of variable)
		// dPhi = |P|^2 sin(theta)/(2*E*(2pi)^3)
		const double dPhip3 = pow(p3.P(),2.)*TMath::Sin(p3.Theta())/(2.0*p3.E()*pow(2.*TMath::Pi(),3));
		const double dPhip4 = pow(p4.P(),2.)*TMath::Sin(p4.Theta())/(2.0*p4.E()*pow(2.*TMath::Pi(),3));

		const double PhaseSpaceOut = dPhip3 * dPhip4;


		// Compute jacobian from change of variable:
		vector<TLorentzVector> momenta;
		momenta.push_back(p1);
		momenta.push_back(p2);
		momenta.push_back(p3);
		momenta.push_back(p4);
		const double jac = computeJacobianF(momenta, SQRT_S);
		if(jac <= 0.){
			//cout << "Jac infinite!" << endl;
			mycount++;
			//continue;
			break;
		}


		//cout << "Found PDF1 = " << pdf1_1 << ", PDF2 = " << pdf1_2 << ", PS in = " << PhaseSpaceIn << ", PS out = " << PhaseSpaceOut << ", jac = " << jac << " , dsdy = " <<  endl;
	        //cout << "===> Matrix element = " << matrix_elements1[0] << ", prod = " << PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac << endl << endl ;	
                //cout << "ME values : " << matrix_elements1[0] << " " << matrix_elements1[1] << endl;
		
		*Value += PhaseSpaceIn * (ME1* pdf1_1 * pdf1_2 + ME2 * pdf2_1 * pdf2_2) * PhaseSpaceOut * jac * Jac_dsdy;
                //*Value += PhaseSpaceIn * (matrix_elements1[0] * pdf1_1 * pdf1_2) * PhaseSpaceOut * jac;

                if (debug) {
                  cout << "pdf products : " << pdf1_1 * pdf1_2  << " " << pdf2_1 * pdf2_2 << endl;
                  const double dsig = (ME1 * pdf1_1 * pdf1_2 + ME2 * pdf2_1 * pdf2_2);
                  cout << "dsig : " << (ME1 * pdf1_1 * pdf1_2) << " dsig two inputs " << dsig << endl;
                }

		// free up memory
		for(unsigned int i = 0; i < p_1.size(); ++i){
			delete p_1.at(i); p_1.at(i) = 0;
		}
                for(unsigned int i = 0; i < p_2.size(); ++i){
                        delete p_2.at(i); p_2.at(i) = 0;
                }


	}

	if(*Value <= 0.){
		mycount++;
		//cout << "Zero!" << endl;
		return 0;
	}

	double flatterJac = range1 * range2 * range3 * range4;
	flatterJac *= M_W*G_W * M_W*G_W;
	flatterJac /= pow(TMath::Cos(y1) * TMath::Cos(y2), 2.);

	// ### FOR NWA
	//double flatterJac = pow(TMath::Pi(),4.) * (M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T);

	//cout << "## Phase Space point done. Integrand = " << integrand << ", flatterjac = " << flatterJac << ", prod = " << integrand*flatterJac <<	endl;

        //cout << "## Phase Space point done." << " flatterjac = " << flatterJac <<   endl;

	*Value *= flatterJac;
	
	return 0;
}

double ME(double *error, TLorentzVector ep, TLorentzVector mum, TLorentzVector Met, TLorentzVector nu1, TLorentzVector nu2, double q1, double q2, double *time){
	
	/*TH1D *hst_Wm = new TH1D("test_mu", "test_1D", 150,0,150);
	TH1D *hst_We = new TH1D("test_ep", "test_1D", 150,0,150);
	TH1D *hst_t = new TH1D("test_t", "test_1D", 100,150,250);
	TH1D *hst_tbar = new TH1D("test_tbar", "test_1D", 100,150,250);*/

	cout << "Initializing integration..." << endl;
	
	TStopwatch chrono;
	chrono.Start();


	MEWeight myWeight("/home/fynu/amertens/scratch/MatrixElement/MG5_aMC_v2_2_3/uu_ww_1d_cpp/Cards/param_card.dat", ep, mum, Met);

	int neval, nfail;
	double mcResult=0, mcError=0, prob=0;
	
	char verbosity = 0; // 0-3
	bool subregion = false; // true = only last set of samples is used for final evaluation of integral
	bool smoothing = false;
	bool retainStateFile = false; // false => delete state file when integration ends
	bool takeOnlyGridFromFile = true; // false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
	unsigned int level = 0; 

	unsigned int flags = setFlags(verbosity, subregion, retainStateFile, level, smoothing, takeOnlyGridFromFile);

	cout << "Starting integration..." << endl;
	cubacores(0,0);	
	Vegas(
		4,						// (int) dimensions of the integrated volume
		1,						// (int) dimensions of the integrand
		(integrand_t) MEFunct,	// (integrand_t) integrand (cast to integrand_t)
		//(integrand_t) BWTest,	// (integrand_t) integrand (cast to integrand_t)
		(void*) &myWeight,		// (void*) pointer to additional arguments passed to integrand
		1,						// (int) maximum number of points given the integrand in each invocation (=> SIMD) ==> PS points = vector of sets of points (x[ndim][nvec]), integrand returns vector of vector values (f[ncomp][nvec])
		0.005,					// (double) requested relative accuracy |-> error < max(rel*value,abs)
		0.,						// (double) requested absolute accuracy |
		flags,					// (int) various control flags in binary format, see setFlags function
		8946,						// (int) seed (seed==0 && no control flag => SOBOL; seed!=0 && control flag==0 => Mersenne Twister)
		0,					// (int) minimum number of integrand evaluations
		100000000,					// (int) maximum number of integrand evaluations (approx.!)
		5000000,					// (int) number of integrand evaluations per interations (to start)
		0,						// (int) increase in number of integrand evaluations per interations
		5,					// (int) batch size for sampling
		0,						// (int) grid number, 1-10 => up to 10 grids can be stored, and re-used for other integrands (provided they are not too different)
		"", 					// (char*) name of state file => state can be stored and retrieved for further refinement
		NULL,					// (int*) "spinning cores": -1 || NULL <=> integrator takes care of starting & stopping child processes (other value => keep or retrieve child processes, probably not useful here)
		&neval,					// (int*) actual number of evaluations done
		&nfail,					// 0=desired accuracy was reached; -1=dimensions out of range; >0=accuracy was not reached
		&mcResult,				// (double*) integration result ([ncomp])
		&mcError,				// (double*) integration error ([ncomp])
		&prob					// (double*) Chi-square p-value that error is not reliable (ie should be <0.95) ([ncomp])
	);
	
	cout << "Integration done." << endl;

	/*for(Long_t loop=0; loop<nPoints; loop++){
		int count_old = mycount;
		FoamX->MakeEvent();					// generate MC event
		FoamX->GetMCvect( MCvect );	 // get generated vector (x,y)
	
		double range1 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
		double y1 = - TMath::ATan(M_W/G_W) + range1 * MCvect[0];
		double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);
		//double range1 = 2000;
		//double s13 = pow(M_W,2.) -1000+ range1*MCvect[0];

		double range2 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
		double y2 = - TMath::ATan(M_T/G_T) + range2 * MCvect[1];
		double s134 = M_T * G_T * TMath::Tan(y2) + pow(M_T,2.);
		//double range2 = 2000;
		//double s134 = pow(M_T,2.) -1000+ range2*MCvect[1];

		double range3 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
		double y3 = - TMath::ATan(M_W/G_W) + range3 * MCvect[2];
		double s25 = M_W * G_W * TMath::Tan(y3) + pow(M_W,2.);
		//double range3 = 2000;
		//double s25 = pow(M_W,2.)-1000 + range3*MCvect[2];

		double range4 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
		double y4 = - TMath::ATan(M_T/G_T) + range4 * MCvect[3];
		double s256 = M_T * G_T * TMath::Tan(y4) + pow(M_T,2.);
		//double range4 = 2000;
		//double s256 = pow(M_T,2.)-1000 + range4*MCvect[3];

		if(count_old == mycount){
			hst_We->Fill(TMath::Sqrt(s13));
			hst_Wm->Fill(TMath::Sqrt(s25));
			hst_t->Fill(TMath::Sqrt(s134));
			hst_tbar->Fill(TMath::Sqrt(s256));
		}
	}*/

	*time = chrono.CpuTime();

	cout << "CPU time : " << chrono.CpuTime() << "	Real-time : " << chrono.RealTime() << endl;
	cout << "Status: " << nfail << ", nr. fails: " << mycount << endl;
	mycount = 0;

	cout << " mcResult= " << mcResult << " +- " << mcError << " in " << neval << " evaluations. Chi-square prob. = " << prob << endl;
	
	/*TCanvas *c = new TCanvas("c","Canvas for plotting",600,600);
	c->cd();
	hst_We->Draw();
	//c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_Enu.png");
	delete hst_We; hst_We = 0;
	delete c; c = 0;
	
	c = new TCanvas("c","Canvas for plotting",600,600);
	c->cd();
	hst_Wm->Draw();
	//c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_Munu.png");
	delete hst_Wm; hst_Wm = 0;
	delete c; c = 0;


	c = new TCanvas("c","Canvas for plotting",600,600);
	c->cd();
	hst_t->Draw();
	//c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_t.png");
	delete hst_t; hst_t = 0;
	delete c; c = 0;
	
	c = new TCanvas("c","Canvas for plotting",600,600);
	c->cd();
	hst_tbar->Draw();
	//c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_tbar.png");
	delete hst_tbar; hst_tbar = 0;
	delete c; c = 0;*/

	*error = mcError;
	if(std::isnan(*error))
	*error = 0.;
	if(std::isnan(mcResult))
	mcResult = 0.;
	return mcResult;
}

std::pair<Double_t,Double_t> MwWeight(TString inlhco, int event_no)
{

ifstream infile(inlhco);
Double_t evt;
Double_t test1;
Double_t test2;
Double_t weight, WW_weight;
Double_t error, WW_error;

WW_weight = 0;
WW_error = 0;

std::string line;


//----------------------- load file  --------------------------------------------------------------
while (std::getline(infile, line))
{ 
  std::istringstream iss(line);
  if (iss >> evt >> test1 >> test2 >> weight >> error){
    //std::cout << evt << " " << test1 << " " << test2 << " " << weight << std::endl;
    if (evt == event_no) {
      WW_weight = weight;
      WW_error = error;
      std::cout << " MW weight " << weight << " +- " << error <<  std::endl; 
      continue;
      }
    }
  }
  infile.close();

return std::make_pair(WW_weight,WW_error);
}



int main(int argc, char *argv[])
{
	TString inputFile(argv[1]);
	TString outputFile(argv[2]);
	int start_evt = atoi(argv[3]);
	int end_evt = atoi(argv[4]);

	//gSystem->Load("libDelphes");

	// Create chain of root trees
	TChain chain("Delphes");
	chain.Add(inputFile);

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

	// Get pointers to branches used in this analysis
	TClonesArray *branchGen = treeReader->UseBranch("Particle");

	cout << "Entries:" << chain.GetEntries() << endl;

	TFile* outFile = new TFile(outputFile, "RECREATE");
	TTree* outTree = chain.CloneTree(0);

	double Weight_TT_cpp, Weight_TT_Error_cpp, Weight_MW, Error_MW;
	bool Weighted_TT_cpp;
	double time = 0;
        double combined_error;
	outTree->Branch("Weight_TT_cpp", &Weight_TT_cpp);
	outTree->Branch("Weight_TT_Error_cpp", &Weight_TT_Error_cpp);
	outTree->Branch("Weighted_TT_cpp", &Weighted_TT_cpp);
	outTree->Branch("Weight_TT_cpp_time", &time);

        outTree->Branch("Weight_MW", &Weight_MW);
        outTree->Branch("Error_MW", &Error_MW);
        outTree->Branch("Comb_error", &combined_error);


	//ofstream fout(outputFile);

	// Full weights
	//double madweight1[10] = {1.58495292058e-21, 2.09681384879e-21, 4.34399623629e-22, 1.68163897955e-22, 3.20350498956e-22, 5.22232034307e-22, 6.04738375743e-21, 9.55643564854e-22, 8.12425265344e-22, 5.81210532053e-23};
	//double madweight2[10] = {1.02514966131e-21, 1.45375719248e-21, 1.65080839221e-22, 1.55653414654e-24, 5.60531044001e-25, 1., 9.70526105314e-24, 3.89103636371e-22, 6.38206925825e-23, 9.37189585544e-26};
	/*double madweight1[10] = {1.48990458909e-21,2.00433822978e-21,4.08998078881e-22,1.56339237714e-22,2.98606743727e-22,4.79498317117e-22,5.63645701583e-21,8.99177777775e-22,7.68316733666e-22,5.42606461617e-23};
	double madweight1Err[10] = {8.63813589113e-24,1.08426062115e-23,2.5750146827e-24,7.0506407196e-25,1.10554655068e-24,2.31140842678e-24,2.71677566322e-23,4.8290429288e-24,1.69718762833e-24,2.66346844676e-25};
	double madweight2[10] = {9.62646303545e-22,1.38143123163e-21,1.54526017444e-22,1.45628835295e-24,6.80263123625e-25,0.,1.07797730384e-23,3.61278172744e-22,6.19087950579e-23,7.20276231557e-26};
	double madweight2Err[10] = {2.96180414077e-24,4.8856162625e-24,1.0218999515e-24,1.29754825587e-25,2.72733072519e-25,0.,4.03010515215e-24,4.29592188061e-24,1.67765665953e-24,8.06569780018e-27};*/

	// NWA weights
	//double madweight1[10] = {1.26069381322e-21, 2.85437676736e-21, 4.8136572269e-22, 1., 3.99894656854e-22, 5.7603822256e-22, 6.99323258475e-21, 1.0892124248e-21, 8.28291668972e-22, 1.};
	//double madweight2[10] = {1.46272073513e-21, 1.51733772927e-21, 1.61193875253e-22, 1., 1., 1., 1., 1., 1., 1.};
	// NWA weights, first sol
	//double madweight1[3] = {8.93501179418e-22, 7.42359826601e-22, 1.49577468649e-22};
	//double madweight2[3] = {1.04113131882e-23, 7.04643552065e-22, 4.3214935529e-23};
	// NWA weights, first sol, ME=1
	//double madweight1[3] = {1.9968889994e-17, 1.10734832869e-17, 2.17966664756e-18};
	//double madweight2[3] = {1.2718723458e-19, 2.38734853175e-17, 6.27800021816e-19};

	if(end_evt >= chain.GetEntries())
		end_evt = chain.GetEntries()-1;

	for(int entry = start_evt; entry <= end_evt ; ++entry){
		
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		chain.GetEntry(entry);

		TLorentzVector gen_ep, gen_mum, gen_Met, gen_nm, gen_ne, gen_b, gen_bbar;

		GenParticle *gen;

		//int count_ep=0, count_mum=0;

		for (Int_t i = 0; i < branchGen->GetEntries(); i++){
			gen = (GenParticle *) branchGen->At(i);
			//cout << "Status=" << gen->Status << ", PID=" << gen->PID << ", E=" << gen->P4().E() << endl;
			if (gen->Status == 1){
				if (gen->PID == -11){
					gen_ep = gen->P4();
					//count_ep++;
				}else if (gen->PID == 13){
					gen_mum = gen->P4();
					//count_mum++;
				}
				else if (gen->PID == 12) {gen_Met += gen->P4();/* gen_ne = gen->P4();*/}
				else if (gen->PID == -14) {gen_Met += gen->P4();/* gen_nm = gen->P4();*/}
        if(gen->PID == 5){ gen_b = gen->P4() ; }
        if(gen->PID == -5){ gen_bbar = gen->P4() ; }
        
        if(gen->PID == 24)
          cout << "W+ mass " << gen->P4().M() << endl;
        if(gen->PID == -24)
          cout << "W- mass " << gen->P4().M() << endl;
      }
		}

		cout << "p2x : " << gen_ne.Px() << " , p2y : " << gen_ne.Py() << " , p2z : " << gen_ne.Pz() << " , E2 : " << gen_ne.E() << endl;
		cout << "p2x : " << gen_nm.Px() << " , p2y : " << gen_nm.Py() << " , p2z : " << gen_nm.Pz() << " , E2 : " << gen_nm.E() << endl;

		double Pzext = (gen_Met+gen_ep+gen_mum+gen_b+gen_bbar).Pz();
		double Eext  = (gen_Met+gen_ep+gen_mum+gen_b+gen_bbar).E();

    cout << "Momentum test " << (gen_ep+gen_mum+gen_Met+gen_b+gen_bbar).Px() << endl;

		cout << "pz   : " << (gen_Met+gen_ep+gen_mum+gen_b+gen_bbar).Pz() << endl;
		cout << "etot : " << (gen_Met+gen_ep+gen_mum+gen_b+gen_bbar).E() << endl;
    double q1 = (Pzext+Eext)/(2*4000.);
    double q2 = (-Pzext+Eext)/(2*4000.);
		cout << "q1 : " << q1 << endl;;
		cout << "q2 : " << q2 << endl;;

    gen_Met.SetPtEtaPhiM(gen_Met.Pt(), 0., gen_Met.Phi(), 0.);


		//if(count_ep != 1 || count_mum != 1)
		//	continue;
		//gen_Met.SetPz(0.);
	
		cout << "From MadGraph:" << endl;
		cout << "Electron" << endl;
		cout << gen_ep.Eta() << "," << gen_ep.Phi() << "," << gen_ep.Pt() << "," << gen_ep.M() << endl;
		cout << "Muon" << endl;
		cout << gen_mum.Eta() << "," << gen_mum.Phi() << "," << gen_mum.Pt() << "," << gen_mum.M() << endl;
		cout << "MET" << endl;
		cout << gen_Met.Eta() << "," << gen_Met.Phi() << "," << gen_Met.Pt() << endl;



		const double q1 = (Pzext+Eext)/(2.*4000.);
		const double q2 = -(Pzext-Eext)/(2.*4000.);

                cout << " q1, q2 : " << q1 << " " << q2 << endl;

		const double s13 = pow((gen_ep+gen_ne).M(),2);
		const double s24 = pow((gen_mum+gen_nm).M(),2);

		cout << " s13 : " << s13 << endl;
		cout << " s24 : " << s24 << endl;

		TLorentzVector p3 = gen_ep;
		TLorentzVector p4 = gen_mum;
		TLorentzVector pb = p3+p4;

                */

		Weight_TT_cpp = 0.;
		Weight_TT_Error_cpp = 0.;


		double weight = 0, temp_time;
		double error = 0;

		weight = ME(&error, gen_ep, gen_mum, gen_Met, &temp_time);
		Weight_TT_cpp = weight;
		Weight_TT_Error_cpp = error;
		time = temp_time;

                std::pair<Double_t,Double_t> MW_WE = MwWeight("/home/fynu/amertens/scratch/MatrixElement/MG5_aMC_v2_2_3/uu_ww_2p/Events/MEMpp_test_2P/weights.out", entry);

                Weight_MW = MW_WE.first;
                Error_MW = MW_WE.second;
	
                cout << " ratio : " << Weight_TT_cpp/Weight_MW <<  endl; 

                combined_error = sqrt(Error_MW*Error_MW+Weight_TT_Error_cpp*Weight_TT_Error_cpp);                

                cout << " diff : " << Weight_MW - Weight_TT_cpp << " +- " << combined_error << endl;

		Weight_TT_Error_cpp = TMath::Sqrt(Weight_TT_Error_cpp);
		Weighted_TT_cpp = true;

		outTree->Fill();

		count_wgt++;
		//fout << endl;
		//break;
	}

	outTree->Write();
	outFile->Close();

	delete treeReader; 
}


