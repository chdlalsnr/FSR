#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/PyMethodBase.h"

int testPyAdaBoostClassification(){

   // Get data file
   std::cout << "Get test data..." << std::endl;
   TString nSignal = "preprocessing/signal.root";
   TFile *fSignal = TFile::Open(nSignal);
   TString nBackground = "preprocessing/background.root";
   TFile *fBackground = TFile::Open(nBackground);

   // Setup PyMVA and factory
   std::cout << "Setup TMVA..." << std::endl;
   TMVA::PyMethodBase::PyInitialize();
   TFile* outputFile = TFile::Open("ResultsTestPyAdaBoostClassification.root", "RECREATE");
   TMVA::Factory *factory = new TMVA::Factory("testPyAdaBoostClassification", outputFile,
      "!V:Silent:Color:!DrawProgressBar:AnalysisType=Classification");

   // Load data
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("datasetTestPyAdaBoostClassification");

   TTree *signal = (TTree*)fSignal->Get("signal;3");
   TTree *background = (TTree*)fBackground->Get("background;3");
   dataloader->AddSignalTree(signal);
   dataloader->AddBackgroundTree(background);

   dataloader->AddVariable("jet1_mass");
   dataloader->AddVariable("jet1_pt");
   dataloader->AddVariable("jet1_eta");
   dataloader->AddVariable("jet1_phi");

   dataloader->AddVariable("jet2_mass");
   dataloader->AddVariable("jet2_pt");
   dataloader->AddVariable("jet2_eta");
   dataloader->AddVariable("jet2_phi");

   dataloader->AddVariable("lep1_mass");
   dataloader->AddVariable("lep1_pt");
   dataloader->AddVariable("lep1_eta");
   dataloader->AddVariable("lep1_phi");

   dataloader->AddVariable("lep2_mass");
   dataloader->AddVariable("lep2_pt");
   dataloader->AddVariable("lep2_eta");
   dataloader->AddVariable("lep2_phi");

   dataloader->AddVariable("dr_ll");
   dataloader->AddVariable("dr_jj");
   dataloader->AddVariable("met_pt");
   dataloader->AddVariable("met_phi");

   dataloader->AddVariable("h1_mass");
   dataloader->AddVariable("h1_pt");
   dataloader->AddVariable("h1_eta");
   dataloader->AddVariable("h1_phi");

   dataloader->AddVariable("h2_mass");
   dataloader->AddVariable("h2_pt");
   dataloader->AddVariable("h2_eta");
   dataloader->AddVariable("h2_phi");

   dataloader->AddVariable("higgsness");
   dataloader->AddVariable("topness");

   dataloader->PrepareTrainingAndTestTree("",
      "SplitMode=Random:NormMode=NumEvents:!V");

   // Book and train method
   factory->BookMethod(dataloader, TMVA::Types::kPyAdaBoost, "PyAdaBoost",
      "H:V:NEstimators=50");
   std::cout << "Train classifier..." << std::endl;
   factory->TrainAllMethods();

   // Clean-up
   delete factory;
   delete dataloader;
   delete outputFile;

   // Setup reader
   UInt_t numEvents = 15000;
   std::cout << "Run reader and classify " << numEvents << " events..." << std::endl;
   TMVA::Reader *reader = new TMVA::Reader("Color:Silent");
   Float_t vars[30];
   reader->AddVariable("jet1_mass", vars+0);
   reader->AddVariable("jet1_pt", vars+1);
   reader->AddVariable("jet1_eta", vars+2);
   reader->AddVariable("jet1_phi", vars+3);
   reader->AddVariable("jet2_mass", vars+4);
   reader->AddVariable("jet2_pt", vars+5);
   reader->AddVariable("jet2_eta", vars+6);
   reader->AddVariable("jet2_phi", vars+7);
   reader->AddVariable("lep1_mass", vars+8);
   reader->AddVariable("lep1_pt", vars+9);
   reader->AddVariable("lep1_eta", vars+10);
   reader->AddVariable("lep1_phi", vars+11);
   reader->AddVariable("lep2_mass", vars+12);
   reader->AddVariable("lep2_pt", vars+13);
   reader->AddVariable("lep2_eta", vars+14);
   reader->AddVariable("lep2_phi", vars+15);
   reader->AddVariable("dr_ll", vars+16);
   reader->AddVariable("dr_jj", vars+17);
   reader->AddVariable("met_pt", vars+18);
   reader->AddVariable("met_phi", vars+19);
   reader->AddVariable("h1_mass", vars+20);
   reader->AddVariable("h1_pt", vars+21);
   reader->AddVariable("h1_eta", vars+22);
   reader->AddVariable("h1_phi", vars+23);
   reader->AddVariable("h2_mass", vars+24);
   reader->AddVariable("h2_pt", vars+25);
   reader->AddVariable("h2_eta", vars+26);
   reader->AddVariable("h2_phi", vars+27);
   reader->AddVariable("higgsness", vars+28);
   reader->AddVariable("topness", vars+29);
   reader->BookMVA("PyAdaBoost", "datasetTestPyAdaBoostClassification/weights/testPyAdaBoostClassification_PyAdaBoost.weights.xml");

   // Get mean response of method on signal and background events
   signal->SetBranchAddress("jet1_mass", vars+0);
   signal->SetBranchAddress("jet1_pt", vars+1);
   signal->SetBranchAddress("jet1_eta", vars+2);
   signal->SetBranchAddress("jet1_phi", vars+3);
   signal->SetBranchAddress("jet2_mass", vars+4);
   signal->SetBranchAddress("jet2_pt", vars+5);
   signal->SetBranchAddress("jet2_eta", vars+6);
   signal->SetBranchAddress("jet2_phi", vars+7);
   signal->SetBranchAddress("lep1_mass", vars+8);
   signal->SetBranchAddress("lep1_pt", vars+9);
   signal->SetBranchAddress("lep1_eta", vars+10);
   signal->SetBranchAddress("lep1_phi", vars+11);
   signal->SetBranchAddress("lep2_mass", vars+12);
   signal->SetBranchAddress("lep2_pt", vars+13);
   signal->SetBranchAddress("lep2_eta", vars+14);
   signal->SetBranchAddress("lep2_phi", vars+15);
   signal->SetBranchAddress("dr_ll", vars+16);
   signal->SetBranchAddress("dr_jj", vars+17);
   signal->SetBranchAddress("met_pt", vars+18);
   signal->SetBranchAddress("met_phi", vars+19);
   signal->SetBranchAddress("h1_mass", vars+20);
   signal->SetBranchAddress("h1_pt", vars+21);
   signal->SetBranchAddress("h1_eta", vars+22);
   signal->SetBranchAddress("h1_phi", vars+23);
   signal->SetBranchAddress("h2_mass", vars+24);
   signal->SetBranchAddress("h2_pt", vars+25);
   signal->SetBranchAddress("h2_eta", vars+26);
   signal->SetBranchAddress("h2_phi", vars+27);
   signal->SetBranchAddress("higgsness", vars+28);
   signal->SetBranchAddress("topness", vars+29);

   background->SetBranchAddress("jet1_mass", vars+0);
   background->SetBranchAddress("jet1_pt", vars+1);
   background->SetBranchAddress("jet1_eta", vars+2);
   background->SetBranchAddress("jet1_phi", vars+3);
   background->SetBranchAddress("jet2_mass", vars+4);
   background->SetBranchAddress("jet2_pt", vars+5);
   background->SetBranchAddress("jet2_eta", vars+6);
   background->SetBranchAddress("jet2_phi", vars+7);
   background->SetBranchAddress("lep1_mass", vars+8);
   background->SetBranchAddress("lep1_pt", vars+9);
   background->SetBranchAddress("lep1_eta", vars+10);
   background->SetBranchAddress("lep1_phi", vars+11);
   background->SetBranchAddress("lep2_mass", vars+12);
   background->SetBranchAddress("lep2_pt", vars+13);
   background->SetBranchAddress("lep2_eta", vars+14);
   background->SetBranchAddress("lep2_phi", vars+15);
   background->SetBranchAddress("dr_ll", vars+16);
   background->SetBranchAddress("dr_jj", vars+17);
   background->SetBranchAddress("met_pt", vars+18);
   background->SetBranchAddress("met_phi", vars+19);
   background->SetBranchAddress("h1_mass", vars+20);
   background->SetBranchAddress("h1_pt", vars+21);
   background->SetBranchAddress("h1_eta", vars+22);
   background->SetBranchAddress("h1_phi", vars+23);
   background->SetBranchAddress("h2_mass", vars+24);
   background->SetBranchAddress("h2_pt", vars+25);
   background->SetBranchAddress("h2_eta", vars+26);
   background->SetBranchAddress("h2_phi", vars+27);
   background->SetBranchAddress("higgsness", vars+28);
   background->SetBranchAddress("topness", vars+29);

   Float_t meanMvaSignal = 0;
   Float_t meanMvaBackground = 0;
   for(UInt_t i=0; i<numEvents; i++){
      signal->GetEntry(i);
      meanMvaSignal += reader->EvaluateMVA("PyAdaBoost");
      background->GetEntry(i);
      meanMvaBackground += reader->EvaluateMVA("PyAdaBoost");
   }
   meanMvaSignal = meanMvaSignal/float(numEvents);
   meanMvaBackground = meanMvaBackground/float(numEvents);

   // Check whether the response is obviously better than guessing

   // NOTE: The scikit-learn AdaBoost classifier returns probability values near
   // to 0.5, thought the separation is quite good. You can check this by looking at
   // the ROC using TMVAGui, which achieves good AUC values.
   // Because of this, the thresholds for passing this test are relaxed to greater or less
   // 0.5 for signal and background respectively.

   std::cout << "Mean MVA response on signal: " << meanMvaSignal << std::endl;
   if(meanMvaSignal < 0.5){
      std::cout << "[ERROR] Mean response on signal is " << meanMvaSignal << " (<0.5)" << std::endl;
      return 1;
   }
   std::cout << "Mean MVA response on background: " << meanMvaBackground << std::endl;
   if(meanMvaBackground > 0.5){
      std::cout << "[ERROR] Mean response on background is " << meanMvaBackground << " (>0.5)" << std::endl;
      return 1;
   }

   return 0;
}

int main(){
   int err = testPyAdaBoostClassification();
   return err;
}
