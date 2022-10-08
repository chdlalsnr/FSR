int tmvaBDTClassifier (TString myMethodList = "") {

  TMVA::Tools::Instance();
  std::map<std::string,int> Use;

  // Use optimization methods
//  Use["Likelihood"] = 0;
//  Use["MLP"] = 0;
  Use["BDT"] = 1;

  // Start point of this module
  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
        std::cout << std::endl;
        return 1;
      }
      Use[regMethod] = 1;
    }
  }

  TFile *input(0);
  // Set input file here
  TString fname = "./tmva_total.root";
  input = TFile::Open( fname );
  std::cout << "--- TMVAClassification : Using input file: " << input->GetName() << std::endl;

  TTree *signalTree = (TTree*)input->Get("signal");
  TTree *background = (TTree*)input->Get("background");

  TString outfileName( "tmvaOutput.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

  dataloader->AddVariable( "HH_mass", 'F' );
  dataloader->AddVariable( "HH_pt",   'F' );
  dataloader->AddVariable( "HH_eta",  'F' );
  dataloader->AddVariable( "HH_phi",  'F' );

//  dataloader->AddVariable( "var4", "Variable 4", "units", 'F' );
//  dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
//  dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  dataloader->AddSignalTree    ( signalTree,     signalWeight );
  dataloader->AddBackgroundTree( background, backgroundWeight );
  // To give different trees for training and testing, do as follows:
  // dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
  // dataloader->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

  //dataloader->SetBackgroundWeightExpression("weight");

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

  // Cut optimisation
//  if (Use["Likelihood"])
//      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
//                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
//  if (Use["MLP"])
//      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
  if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv ) {

  TString methodList;
  for (int i=1; i<argc; i++) {

    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(",");
    methodList += regMethod; }

  return tmvaBDTClassifier(methodList);
}

