#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "cuba.h"

#include "Math/Vector4D.h"
#include "TH1D.h"

#include "MEWeight.h"
#include "MEEvent.h"
#include "jacobianD.h"
#include "transferFunction.h"
#include "utils.h"

#define VEGAS
//#define SUAVE

using namespace std;

MEWeight::MEWeight(CPPProcess &process, const std::string pdfName, const std::string fileTF):
  _process(process),
  _pdf( LHAPDF::mkPDF(pdfName, 0) ),
  _recEvent( new MEEvent() ),
  _TF( new TransferFunction(fileTF) ),
  _isrCorrection(ISRCorrection::transverseISRBoost)
  {

  cout << "Initializing Matrix Element computation with:" << endl;
  cout << "PDF " << pdfName << endl;
  cout << "TF file " << fileTF << endl;
}

MEEvent* MEWeight::GetEvent(){
  return _recEvent;
}

void MEWeight::SetEvent(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met){
  _recEvent->SetVectors(ep, mum, b, bbar, met);
}

void MEWeight::AddTF(const std::string particleName, const std::string histName){
  _TF->DefineComponent(particleName, histName);
}

void MEWeight::AddInitialState(int pid1, int pid2){
  // We must have a quark or a gluon as initial state!
  if( (abs(pid1) > 5 && pid1 != 21) || (abs(pid2) > 5 && pid2 != 21) ){
    cerr << "Warning: initial state (" << pid1 << "," << pid2 << ") is not valid! I'm not including it.\n";
    return;
  }

  // Sort the initial state PIDs in order to check more easily if they are already included
  if(pid1 > pid2)
    swap(pid1, pid2);
  pair<int, int> initialState(pid1, pid2);

  if( find(_initialStates.begin(), _initialStates.end(), initialState) == _initialStates.end() ){
    _initialStates.push_back(initialState);
    // Don't forget to add the crossed possibility if the initial particles are different
    if(pid1 != pid2){
      swap(initialState.first, initialState.second);
      _initialStates.push_back(initialState);
    }
  }else{
    cerr << "Warning: initial state (" << pid1 << "," << pid2 << ") has already been defined. I'm not including it again.\n";
  }
}

void MEWeight::SetISRCorrection(const ISRCorrection newISRCorrection){
  _isrCorrection = newISRCorrection;
}

double MEWeight::ComputeWeight(double &error){
  
  cout << "Initializing integration..." << endl;

  int neval, nfail;
#ifdef SUAVE 
  int nsubregions;
#endif
  double mcResult=0, prob=0;
  
  char verbosity = 3; // 0-3
  bool subregion = false; // true = only last set of samples is used for final evaluation of integral
  bool smoothing = false;
  bool retainStateFile = false; // false => delete state file when integration ends
  bool takeOnlyGridFromFile = true; // false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
  unsigned int level = 0; 

  unsigned int flags = setFlags(verbosity, subregion, retainStateFile, level, smoothing, takeOnlyGridFromFile);

  cout << "Starting integration..." << endl << endl;

  cubacores(0, 0);           // This is mandatory if the integrand wants to *modify* something in the MEWeight object passed as argument
#ifdef VEGAS
  Vegas
#endif
#ifdef SUAVE
  Suave
#endif
  (
    9,                      // (int) dimensions of the integrated volume
    1,                      // (int) dimensions of the integrand
    (integrand_t) CUBAIntegrand,  // (integrand_t) integrand (cast to integrand_t)
    (void*) this,           // (void*) pointer to additional arguments passed to integrand
    1,                      // (int) maximum number of points given the integrand in each invocation (=> SIMD) ==> PS points = vector of sets of points (x[ndim][nvec]), integrand returns vector of vector values (f[ncomp][nvec])
    0.005,                  // (double) requested relative accuracy  /
    0.,                     // (double) requested absolute accuracy /-> error < max(rel*value,abs)
    flags,                  // (int) various control flags in binary format, see setFlags function
    0,                      // (int) seed (seed==0 => SOBOL; seed!=0 && control flag "level"==0 => Mersenne Twister)
    0,                      // (int) minimum number of integrand evaluations
    400000,                 // (int) maximum number of integrand evaluations (approx.!)
#ifdef VEGAS
    25000,                  // (int) number of integrand evaluations per interations (to start)
    0,                      // (int) increase in number of integrand evaluations per interations
    12500,                   // (int) batch size for sampling
    0,                      // (int) grid number, 1-10 => up to 10 grids can be stored, and re-used for other integrands (provided they are not too different)
#endif
#ifdef SUAVE 
    100000,                  // (int) number of new integrand evaluations in each subdivision
    5000,                      // (int) minimum number of samples a previous iteration must contribute to a subregion, to be considered to that subregion's contribution to the integral
    2,                      // (int) exponent in the norm used to compute fluctuations of a sample
#endif
    "",                     // (char*) name of state file => state can be stored and retrieved for further refinement
    NULL,                   // (int*) "spinning cores": -1 || NULL <=> integrator takes care of starting & stopping child processes (other value => keep or retrieve child processes, probably not useful here)
#ifdef SUAVE
    &nsubregions,           // (int*) actual number of subregions needed
#endif
    &neval,                 // (int*) actual number of evaluations done
    &nfail,                 // 0=desired accuracy was reached; -1=dimensions out of range; >0=accuracy was not reached
    &mcResult,              // (double*) integration result ([ncomp])
    &error,                 // (double*) integration error ([ncomp])
    &prob                   // (double*) Chi-square p-value that error is not reliable (ie should be <0.95) ([ncomp])
  );
  
  cout << "Integration done." << endl;

  cout << " mcResult= " << mcResult << " +- " << error << " in " << neval << " evaluations. Chi-square prob. = " << prob << endl << endl;

  if(std::isnan(error))
  error = 0.;
  if(std::isnan(mcResult))
  mcResult = 0.;
  return mcResult;
}

MEWeight::~MEWeight(){
  cout << "Deleting PDF" << endl;
  delete _pdf; _pdf = nullptr;
  cout << "Deleting myEvent" << endl;
  delete _recEvent; _recEvent = nullptr;
  cout << "Deleting myTF" << endl;
  delete _TF; _TF = nullptr;
}

// Wrapper function passed to CUBA, simply calls MEWeight::Integrand (where the MEWeight instance is passed as "input" to the wrapper), passing the PS point and the weight
int CUBAIntegrand(const int *nDim, const double* psPoint, const int *nComp, double *value, void *inputs, const int *nVec, const int *core, const double *weight){
  //cout << endl << endl << endl << "########## Starting phase-space point ############" << endl << endl;

  //cout << "Inputs = [" << Xarg[0] << "," << Xarg[1] << "," << Xarg[2] << "," << Xarg[3] << "," << Xarg[4] << "," << Xarg[5] << "," << Xarg[6] << "," << Xarg[7] << "]" << endl;
  
  *value = static_cast<MEWeight*>(inputs)->Integrand(psPoint, weight);

  return 0;
}
