#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES // include M_PI constant
#include <cmath>

#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Boost.h"
#include "TMath.h"

#include "MEWeight.h"
#include "jacobianD.h"
#include "utils.h"

#define M_T 173.
#define G_T 1.4915

#define M_W 80.419
#define G_W 2.0476

#define SQRT_S 13000

double MEWeight::Integrand(const double* psPoint, const double *weight){
  using namespace std;

  double returnValue = 0.;

  for(int i=0; i<4; ++i){
    if(psPoint[i] == 1.)
      return 0;
  }

  // Handle the 2 b-jet permutations using Monte Carlo
  ROOT::Math::PtEtaPhiEVector p3rec;
  ROOT::Math::PtEtaPhiEVector p5rec;
  if(psPoint[8] <= 0.5){
    p3rec = _recEvent->GetP3();
    p5rec = _recEvent->GetP5();
  }else{
    p3rec = _recEvent->GetP5();
    p5rec = _recEvent->GetP3();
  }
  
  const ROOT::Math::PtEtaPhiEVector p4rec( _recEvent->GetP4() );
  const ROOT::Math::PtEtaPhiEVector p6rec( _recEvent->GetP6() );
  const ROOT::Math::PtEtaPhiEVector RecMet( _recEvent->GetMet() );

  // Define ISR vector from the observed particles and MET
  // Do this in PxPyPzE basis as this will be more practical for the next stages (transverse boost)
  ROOT::Math::PxPyPzEVector ISR( -(p3rec + p4rec + p5rec + p6rec + RecMet) );

  ///// Transfer functions

  double TFreturnValue = 1.;

  ROOT::Math::PtEtaPhiEVector p3gen = p3rec;
  const double E3rec = p3rec.E();
  const double p3DeltaRange = _TF->GetDeltaRange("electron", E3rec);
  const double E3gen = E3rec - _TF->GetDeltaMax("electron", E3rec) + p3DeltaRange * psPoint[4];
  const double pt3gen = sqrt( SQ(E3gen) - SQ(p3rec.M()) ) / cosh(p3rec.Eta());
  p3gen.SetCoordinates(pt3gen, p3rec.Eta(), p3rec.Phi(), E3gen);
  if(p3DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("electron", E3rec, E3gen) * p3DeltaRange * dEoverdP(E3gen, p3gen.M());

  ROOT::Math::PtEtaPhiEVector p5gen = p5rec;
  const double E5rec = p5rec.E();
  const double p5DeltaRange = _TF->GetDeltaRange("muon", E5rec);
  const double E5gen = E5rec - _TF->GetDeltaMax("muon", E5rec) + p5DeltaRange * psPoint[6];
  const double pt5gen = sqrt( SQ(E5gen) - SQ(p5rec.M()) ) / cosh(p5rec.Eta());
  p5gen.SetCoordinates(pt5gen, p5rec.Eta(), p5rec.Phi(), E5gen);
  if(p5DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("muon", E5rec, E5gen) * p5DeltaRange * dEoverdP(E5gen, p5gen.M());

  ROOT::Math::PtEtaPhiEVector p4gen = p4rec;
  const double E4rec = p4rec.E();
  const double p4DeltaRange = _TF->GetDeltaRange("jet", E4rec);
  const double E4gen = E4rec - _TF->GetDeltaMax("jet", E4rec) + p4DeltaRange * psPoint[5];
  const double pt4gen = sqrt( SQ(E4gen) - SQ(p4rec.M()) ) / cosh(p4rec.Eta());
  p4gen.SetCoordinates(pt4gen, p4rec.Eta(), p4rec.Phi(), E4gen);
  if(p4DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("jet", E4rec, E4gen) * p4DeltaRange * dEoverdP(E4gen, p4gen.M());

  ROOT::Math::PtEtaPhiEVector p6gen = p6rec;
  const double E6rec = p6rec.E();
  const double p6DeltaRange = _TF->GetDeltaRange("jet", E6rec);
  const double E6gen = p6rec.E() - _TF->GetDeltaMax("jet", E6rec) + p6DeltaRange * psPoint[7];
  const double pt6gen = sqrt( SQ(E6gen) - SQ(p6rec.M()) ) / cosh(p6rec.Eta());
  p6gen.SetCoordinates(pt6gen, p6rec.Eta(), p6rec.Phi(), E6gen);
  if(p6DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("jet", E6rec, E6gen) * p6DeltaRange * dEoverdP(E6gen, p6gen.M());

  // In the following, we want to use PxPyPzE vectors, since the change of variables is done over those variables => we are already in the right basis, no need to recompute quantities every time
  ROOT::Math::PxPyPzEVector p3(p3gen);
  ROOT::Math::PxPyPzEVector p4(p4gen);
  ROOT::Math::PxPyPzEVector p5(p5gen);
  ROOT::Math::PxPyPzEVector p6(p6gen);
  ROOT::Math::PxPyPzEVector Met(RecMet);

  //cout << "Final TF = " << TFreturnValue << endl;

  // We flatten the Breit-Wigners by doing a change of variable for each resonance separately
  // The new integration variables are now the Lorentz invariants of the Breit-Wigners (sXXX)
  // Each transformation also brings about its jacobian factor
  double s13, jac13, s134, jac134, s25, jac25, s256, jac256;
  flattenBW(psPoint[0], M_W, G_W, s13, jac13);
  flattenBW(psPoint[1], M_T, G_T, s134, jac134);
  flattenBW(psPoint[2], M_W, G_W, s25, jac25);
  flattenBW(psPoint[3], M_T, G_T, s256, jac256);
  double flatterJac = jac13 * jac134 * jac25 * jac256;

  if(s13 > s134 || s25 > s256 || s13 < p3.M() || s25 < p5.M() || s134 < p4.M() || s256 < p6.M())
    return 0;

  //cout << "weight = " << *weight << endl;
  
  std::vector<ROOT::Math::PxPyPzEVector> p1vec, p2vec;

  ComputeTransformD(s13, s134, s25, s256,
                    p3, p4, p5, p6, Met, ISR,
                    p1vec, p2vec);

  int countSol = 0;

  for(unsigned short i = 0; i < p1vec.size(); ++i){

    const ROOT::Math::PxPyPzEVector &p1 = p1vec[i];
    const ROOT::Math::PxPyPzEVector &p2 = p2vec[i];

    /*const ROOT::Math::PxPyPzEVector p13 = p1 + p3;
    const ROOT::Math::PxPyPzEVector p134 = p1 + p3 + p4;
    const ROOT::Math::PxPyPzEVector p25 = p2 + p5;
    const ROOT::Math::PxPyPzEVector p256 = p2 + p5 + p6;

    cout << "Solution " << i << ":" << endl;
    cout << "Input: W+ mass=" << sqrt(s13) << ", Top mass=" << sqrt(s134) << ", W- mass=" << sqrt(s25) << ", Anti-top mass=" << sqrt(s256) << endl;
    cout << "Output: W+ mass=" << p13.M() << ", Top mass=" << p134.M() << ", W- mass=" << p25.M() << ", Anti-top mass=" << p256.M() << endl << endl;
    cout << "Accuracies: " << endl;
    cout << "W+ mass=" << (sqrt(s13)-p13.M())/p13.M() << endl;
    cout << "Top mass=" << (sqrt(s134)-p134.M())/p134.M() << endl;
    cout << "W- mass=" << (sqrt(s25)-p25.M())/p25.M() << endl;
    cout << "Anti-top mass=" << (sqrt(s256)-p256.M())/p256.M() << endl;
    cout << "Neutrino mass=" << p1.M() << endl;
    cout << "Anti-neutrino mass=" << p2.M() << endl;
    cout << "Transverse momentum conservation=" << (p1+p2+p3+p4+p5+p6+ISR).Pt() << endl << endl;*/
    
    /*cout << "Electron (E,Px,Py,Pz) = ";
    cout << p3.E() << "," << p3.Px() << "," << p3.Py() << "," << p3.Pz() << endl;
    cout << "Electron neutrino (E,Px,Py,Pz) = ";
    cout << p1.E() << "," << p1.Px() << "," << p1.Py() << "," << p1.Pz() << endl;
    cout << "b quark (E,Px,Py,Pz) = ";
    cout << p4.E() << "," << p4.Px() << "," << p4.Py() << "," << p4.Pz() << endl;
    cout << "Muon (E,Px,Py,Pz) = ";
    cout << p5.E() << "," << p5.Px() << "," << p5.Py() << "," << p5.Pz() << endl;
    cout << "Muon neutrino (E,Px,Py,Pz) = ";
    cout << p2.E() << "," << p2.Px() << "," << p2.Py() << "," << p2.Pz() << endl;
    cout << "Anti b quark (E,Px,Py,Pz) = ";
    cout << p6.E() << "," << p6.Px() << "," << p6.Py() << "," << p6.Pz() << endl << endl;*/
  
    const ROOT::Math::PxPyPzEVector tot = p1 + p2 + p3 + p4 + p5 + p6;

    ROOT::Math::PxPyPzEVector parton1, parton2;
    double q1Pz(0), q2Pz(0), ETot(0), PzTot(0);
    
    /*cout << "**********************" << endl;
    cout << "Total: " << tot  << endl;
    cout << "ISR: " << ISR << endl;*/
   
    if(_isrCorrection == ISRCorrection::transverseISRBoost){
      
      //////////////////////////// ISR CORRECTION ////////////////////////////////
    
      // Define boost that puts the transverse total momentum vector in its CoM frame 
      ROOT::Math::PxPyPzEVector tempTot( tot );
      tempTot.SetPz(0.);
      ROOT::Math::XYZVector isrDeBoostVector( tempTot.BoostToCM() );
      
      // In the "transverse" CoM frame, use total Pz and E to define initial longitudinal quark momenta
      //cout << "Boost: " << isrDeBoostVector << endl;
      const ROOT::Math::Boost isrDeBoost( isrDeBoostVector );
      const ROOT::Math::PxPyPzEVector newTot( isrDeBoost*tot );
      ETot = newTot.E();
      PzTot = newTot.Pz();

      q1Pz = (PzTot + ETot)/2.;
      q2Pz = (PzTot - ETot)/2.;
  
      if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.)
        continue;
      
      parton1.SetCoordinates(0., 0., q1Pz, q1Pz);
      parton2.SetCoordinates(0., 0., q2Pz, abs(q2Pz));

      /*cout << "Before:" << endl;
      cout << " Parton1: " << parton1 << endl;
      cout << " Parton2: " << parton2 << endl;*/
   
      // Boost initial parton momenta by the opposite of the transverse boost needed to put the whole system in its CoM
      const ROOT::Math::Boost isrBoost( -isrDeBoostVector );
      parton1 = isrBoost*parton1;
      parton2 = isrBoost*parton2;

      /*cout << "After:" << endl;
      cout << " Parton1: " << parton1 << endl;
      cout << " Parton2: " << parton2 << endl;*/
   
      //ROOT::Math::PxPyPzEVector testT = parton1 + parton2 + ISR;
      //ROOT::Math::PxPyPzEVector testT = -parton1 - parton2 + p1 + p2 + p3 + p4 + p5 + p6;
      //cout << "Test transverse: " << testT << endl;

    } else if (_isrCorrection == ISRCorrection::noCorrection ){
    
      ///////////////// NO ISR CORRECTION //////////////////////////////////////////
      ETot = tot.E();
      PzTot = tot.Pz();

      q1Pz = (PzTot + ETot)/2.;
      q2Pz = (PzTot - ETot)/2.;
      
      if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.)
        continue;
  
      parton1.SetCoordinates(0., 0., q1Pz, q1Pz);
      parton2.SetCoordinates(0., 0., q2Pz, abs(q2Pz));
      //////////////////////////////////////////////////////////////////////////////   

    }

    q1Pz = parton1.Pz();
    q2Pz = parton2.Pz();
    
    //cout << "===> Eext=" << ETot << ", Pzext=" << PzTot << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl << endl;
    
    // Compute jacobian from change of variable:
    const double jacobian = computeJacobianD(p1, p2, p3, p4, p5, p6, SQRT_S);
    if(jacobian <= 0.){
      cout << "Jac infinite!" << endl;
      continue;
    }

    // Compute the Pdfs
    const double x1 = abs(q1Pz/(SQRT_S/2.));
    const double x2 = abs(q2Pz/(SQRT_S/2.));

    // Compute flux factor 1/(2*x1*x2*s)
    const double phaseSpaceIn = 1.0 / ( 2. * x1 * x2 * SQ(SQRT_S) ); 

    // Compute phase space density for observed particles (not concerned by the change of variable)
    // dPhi = |P|^2 sin(theta)/(2*E*(2pi)^3)
    const double dPhip3 = SQ(p3.P())*sin(p3.Theta())/(2.0*p3.E()*CB(2.*M_PI));
    const double dPhip4 = SQ(p4.P())*sin(p4.Theta())/(2.0*p4.E()*CB(2.*M_PI));
    const double dPhip5 = SQ(p5.P())*sin(p5.Theta())/(2.0*p5.E()*CB(2.*M_PI));
    const double dPhip6 = SQ(p6.P())*sin(p6.Theta())/(2.0*p6.E()*CB(2.*M_PI));
    const double phaseSpaceOut = dPhip5 * dPhip6 * dPhip3 * dPhip4;

    // Define initial momenta to be passed to matrix element
    std::vector< std::vector<double> > initialMomenta = 
    {
      { parton1.E(), parton1.Px(), parton1.Py(), parton1.Pz() },
      { parton2.E(), parton2.Px(), parton2.Py(), parton2.Pz() },
    };
    
    // Define final PID and momenta to be passed to matrix element
    std::vector< std::pair<int, std::vector<double> > > finalState = 
    {
      std::make_pair<int, std::vector<double> >( -11, { p3.E(), p3.Px(), p3.Py(), p3.Pz() } ),
      std::make_pair<int, std::vector<double> >(  12, { p1.E(), p1.Px(), p1.Py(), p1.Pz() } ),
      std::make_pair<int, std::vector<double> >(   5, { p4.E(), p4.Px(), p4.Py(), p4.Pz() } ),
      std::make_pair<int, std::vector<double> >(  13, { p5.E(), p5.Px(), p5.Py(), p5.Pz() } ),
      std::make_pair<int, std::vector<double> >( -14, { p2.E(), p2.Px(), p2.Py(), p2.Pz() } ),
      std::make_pair<int, std::vector<double> >(  -5, { p6.E(), p6.Px(), p6.Py(), p6.Pz() } ),
    };

    // Evaluate matrix element
    std::map< std::pair<int, int>, double > matrixElements = getMatrixElements(initialMomenta, finalState);

    double thisSolResult = phaseSpaceIn * phaseSpaceOut * jacobian * flatterJac * TFreturnValue;

    double pdfMESum = 0.;
    // If no initial states have been defined explicitly, loop over all states returned by the matrix element
    if(!_initialStates.size()){
      for(auto const &me: matrixElements){
        const double pdf1 = ComputePdf(me.first.first, x1, SQ(M_T));
        const double pdf2 = ComputePdf(me.first.second, x2, SQ(M_T));
        pdfMESum += me.second * pdf1 * pdf2;
        //cout << "Initial state (" << me.first.first << ", " << me.first.second << "): " << me.second << endl;
      }
    }else{
      // Otherwise, loop over all states defined by user
      for(auto const &initialState: _initialStates){
        const double pdf1 = ComputePdf(initialState.first, x1, SQ(M_T));
        const double pdf2 = ComputePdf(initialState.second, x2, SQ(M_T));
        pdfMESum += matrixElements[initialState] * pdf1 * pdf2;
        //cout << "Initial state (" << initialState.first << ", " << initialState.second << "): " << matrixElements[initialState] << endl;
      }
    }

    thisSolResult *= pdfMESum;
    returnValue += thisSolResult; 
    
    // Check whether the next solutions for the neutrinos are the same => don't redo all this!
    int countEqualSol = 1;
    for(unsigned int j = i+1; j<p1vec.size(); j++){
      if(p1 == p1vec.at(j) && p2 == p2vec.at(j)){
        returnValue += thisSolResult;
        countEqualSol++;
      }
    }

    // If we have included the next solutions already, skip them!
    i += countEqualSol - 1;
    countSol += countEqualSol;
  }

  //cout << "## Phase Space point done. Integrand = " << returnValue << endl; 

  return returnValue;
}
