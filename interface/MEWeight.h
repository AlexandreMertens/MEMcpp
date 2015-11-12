#ifndef _INC_MEWEIGHT
#define _INC_MEWEIGHT

#include <string>
#include <vector>
#include <utility>
#include <memory>

#include "process_base_classes.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "Math/Vector4D.h"

#include "transferFunction.h"
#include "MEEvent.h"

int CUBAIntegrand(const int *nDim, const double* psPoint, const int *nComp, double *value, void *inputs, const int *nVec, const int *core, const double *weight);

enum class ISRCorrection {
  noCorrection = 0,
  transverseISRBoost = 1
};

class MEWeight{
  public:

  double Integrand(const double* psPoint, const double *weight);
  inline double ComputePdf(const int &pid, const double &x, const double &q2);
  inline std::map< std::pair<int, int>, double > getMatrixElements(const std::vector< std::vector<double> > &initialMomenta, const std::vector< std::pair<int, std::vector<double> > > &finalState) const { return _process.sigmaKin(initialMomenta, finalState); }
  double ComputeWeight(double &error);
  MEEvent* GetEvent();
  void SetEvent(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met);
  void AddTF(const std::string particleName, const std::string histName);
  void AddInitialState(int pid1, int pid2);
  void SetISRCorrection(const ISRCorrection newISRCorrection);

  MEWeight(CPPProcess &process, const std::string pdfName, const std::string fileTF);
  ~MEWeight();

  private:

  std::vector< std::pair<int, int> > _initialStates;
  CPPProcess &_process;
  LHAPDF::PDF* _pdf;
  MEEvent* _recEvent;
  TransferFunction* _TF;
  ISRCorrection _isrCorrection;
};

inline double MEWeight::ComputePdf(const int &pid, const double &x, const double &q2){
  // return f(pid,x,q2)
  if(x <= 0 || x >= 1 || q2 <= 0){
    std::cout << "WARNING: PDF x or Q^2 value out of bounds!" << std::endl;
    return 0.;
  }else{
    return _pdf->xfxQ2(pid, x, q2)/x;
  }
}

#endif
