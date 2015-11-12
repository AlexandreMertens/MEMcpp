#include <string>
#include <iostream>

#include "classes/DelphesClasses.h"

#include "Math/Vector4D.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "SubProcesses/P1_Sigma_sm_gg_epvebmumvmxbx/cppmem_pp_ttx_epmum.h"
typedef cppmem_pp_ttx_epmum MyTTbarProcess;

#include "MEWeight.h"

using namespace std;

int main(int argc, char *argv[])
{
  std::string inputFile(argv[1]);
  std::string outputFile(argv[2]);
  std::string fileTF(argv[3]);
  int start_evt = atoi(argv[4]);
  int end_evt = atoi(argv[5]);

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile.c_str());

  double MadWeight, MadWeight_Error;
  chain.SetBranchAddress("Weight_TT", &MadWeight);
  chain.SetBranchAddress("Weight_TT_Error", &MadWeight_Error);

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = NULL;
  chain.SetBranchAddress("Particle", &branchGen);

  cout << "Entries:" << chain.GetEntries() << endl;

  TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
  TTree* outTree = chain.CloneTree(0);

  double Weight_TT_cpp, Weight_TT_Error_cpp;
  bool Weighted_TT_cpp;
  double time;
  outTree->Branch("Weight_TT_cpp", &Weight_TT_cpp);
  outTree->Branch("Weight_TT_Error_cpp", &Weight_TT_Error_cpp);
  outTree->Branch("Weighted_TT_cpp", &Weighted_TT_cpp);
  outTree->Branch("Weight_TT_cpp_me", &time);

  if(end_evt >= chain.GetEntries())
    end_evt = chain.GetEntries()-1;

  // Create CPPProcess and MEWeight objects
  Parameters_sm myParams("MatrixElements/cppmem_pp_ttx_epmum/Cards/param_card.dat");
  MyTTbarProcess myProcess(myParams);
  MEWeight* myWeight = new MEWeight(myProcess, "cteq6l1", fileTF);
  
  myWeight->AddTF("electron", "Binned_Egen_DeltaE_Norm_ele");
  myWeight->AddTF("muon", "Binned_Egen_DeltaE_Norm_muon");
  myWeight->AddTF("jet", "Binned_Egen_DeltaE_Norm_jet");

  /*myWeight->AddInitialState(21, 21);
  myWeight->AddInitialState(1, -1);
  myWeight->AddInitialState(2, -2);
  myWeight->AddInitialState(3, -3);
  myWeight->AddInitialState(4, -4);*/

  for(int entry = start_evt; entry <= end_evt ; ++entry){
    // Load selected branches with data from specified event
    chain.GetEntry(entry);

    ROOT::Math::PtEtaPhiEVector gen_ep, gen_mum, gen_b, gen_bbar, gen_nue, gen_num, gen_Met;

    GenParticle *gen;

    for (int i = 0; i < branchGen->GetEntries(); i++){
      gen = (GenParticle*) branchGen->At(i);
      //cout << "Status=" << gen->Status << ", PID=" << gen->PID << ", E=" << gen->P4().E() << endl;
      if (gen->Status == 1){
        if (gen->PID == -11) gen_ep.SetCoordinates(gen->P4().Pt(), gen->P4().Eta(), gen->P4().Phi(), gen->P4().E());
        else if (gen->PID == 13) gen_mum.SetCoordinates(gen->P4().Pt(), gen->P4().Eta(), gen->P4().Phi(), gen->P4().E());
        else if (gen->PID == 12) gen_nue.SetCoordinates(gen->P4().Pt(), gen->P4().Eta(), gen->P4().Phi(), gen->P4().E());
        else if (gen->PID == -14) gen_num.SetCoordinates(gen->P4().Pt(), gen->P4().Eta(), gen->P4().Phi(), gen->P4().E());
        else if (gen->PID == 5) gen_b.SetCoordinates(gen->P4().Pt(), gen->P4().Eta(), gen->P4().Phi(), gen->P4().E());
        else if (gen->PID == -5) gen_bbar.SetCoordinates(gen->P4().Pt(), gen->P4().Eta(), gen->P4().Phi(), gen->P4().E());
      }
    }

    gen_Met = gen_num + gen_nue;

    cout << "From MadGraph:" << endl;
    cout << "Electron" << endl;
    cout << gen_ep.E() << "," << gen_ep.Px() << "," << gen_ep.Py() << "," << gen_ep.Pz() << endl;
    cout << "b quark" << endl;
    cout << gen_b.E() << "," << gen_b.Px() << "," << gen_b.Py() << "," << gen_b.Pz() << endl;
    cout << "Muon" << endl;
    cout << gen_mum.E() << "," << gen_mum.Px() << "," << gen_mum.Py() << "," << gen_mum.Pz() << endl;
    cout << "Anti b quark" << endl;
    cout << gen_bbar.E() << "," << gen_bbar.Px() << "," << gen_bbar.Py() << "," << gen_bbar.Pz() << endl;
    cout << "MET" << endl;
    cout << gen_Met.E() << "," << gen_Met.Px() << "," << gen_Met.Py() << "," << gen_Met.Pz() << endl << endl;

    Weight_TT_cpp = 0.;
    Weight_TT_Error_cpp = 0.;
    time = 0.;
  
    TStopwatch chrono;
    chrono.Start();

    /*for(int permutation = 1; permutation <= 2; permutation++){
      double weight = 0;
      double error = 0;

      if(permutation == 1)
        myWeight->SetEvent(gen_ep, gen_mum, gen_b, gen_bbar, gen_Met);
      if(permutation == 2)
        myWeight->SetEvent(gen_ep, gen_mum, gen_bbar, gen_b, gen_Met);

      weight = myWeight->ComputeWeight(error)/2.; 

      Weight_TT_cpp += weight;
      Weight_TT_Error_cpp += pow(error/2,2.);
    }
    Weight_TT_Error_cpp = TMath::Sqrt(Weight_TT_Error_cpp);*/

    myWeight->SetEvent(gen_ep, gen_mum, gen_b, gen_bbar, gen_Met);
    Weight_TT_cpp = myWeight->ComputeWeight(Weight_TT_Error_cpp);

    time = chrono.CpuTime();
    
    Weighted_TT_cpp = true;

    cout << "====> Event " << entry << ": weight = " << Weight_TT_cpp << " +- " << Weight_TT_Error_cpp << endl;
    cout << "      CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << endl;
    cout << "      MadWeight: " << MadWeight << " +- " << MadWeight_Error << endl << endl;

    outTree->Fill();
  }
  
  outFile->cd();
  outTree->Write();
  
  delete myWeight; myWeight = nullptr;
  delete outFile; outFile = nullptr;
}

