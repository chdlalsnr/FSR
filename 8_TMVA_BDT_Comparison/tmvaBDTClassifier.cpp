
int tmvaBDTClassifier (TString myMethodList = "") {

    TMVA::Tools::Instance();
    std::map<std::string,int> Use;
    
    // Use optimization methods
    Use["Likelihood"] = 1;
    Use["MLP"] = 1;
    Use["SVM"] = 1;
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
            } Use[regMethod] = 1;
        }
    }
    
    TFile *f1(0);
    TFile *f2(0);
    // Set input file here
    TString s1 = "signal.root";
    f1 = TFile::Open(s1);
    TString s2 = "background.root";
    f2 = TFile::Open(s2);
    std::cout << "--- TMVAClassification : Using input files: " << f1->GetName() << ", " << f2->GetName() << std::endl;
    
    TTree *signalTree = (TTree*)f1->Get("signal");
    TTree *backgroundTree = (TTree*)f2->Get("background");
    
    TString outfileName( "tmvaOutput.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
    
    dataloader->AddVariable( "jet1_mass", 'F' );
    dataloader->AddVariable( "jet1_pt",   'F' );
    dataloader->AddVariable( "jet1_eta",  'F' );
    dataloader->AddVariable( "jet1_phi",  'F' );
    dataloader->AddVariable( "jet2_mass", 'F' );
    dataloader->AddVariable( "jet2_pt",   'F' );
    dataloader->AddVariable( "jet2_eta",  'F' );
    dataloader->AddVariable( "jet2_phi",  'F' );
//    dataloader->AddVariable( "lep1_mass", 'F' );
    dataloader->AddVariable( "lep1_pt",   'F' );
    dataloader->AddVariable( "lep1_eta",  'F' );
    dataloader->AddVariable( "lep1_phi",  'F' );
//    dataloader->AddVariable( "lep2_mass", 'F' );
    dataloader->AddVariable( "lep2_pt",   'F' );
    dataloader->AddVariable( "lep2_eta",  'F' );
    dataloader->AddVariable( "lep2_phi",  'F' );
    dataloader->AddVariable( "dr_ll",  'F' );
    dataloader->AddVariable( "dr_jj",  'F' );
    dataloader->AddVariable( "met_pt",  'F' );
    dataloader->AddVariable( "met_phi",  'F' );
    dataloader->AddVariable( "h1_mass", 'F' );
    dataloader->AddVariable( "h1_pt",   'F' );
    dataloader->AddVariable( "h1_eta",  'F' );
    dataloader->AddVariable( "h1_phi",  'F' );
    dataloader->AddVariable( "h2_mass", 'F' );
    dataloader->AddVariable( "h2_pt",   'F' );
    dataloader->AddVariable( "h2_eta",  'F' );
    dataloader->AddVariable( "h2_phi",  'F' );
    dataloader->AddVariable( "higgsness",  'F' );
    dataloader->AddVariable( "topness",  'F' );

    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    
    dataloader->AddSignalTree    ( signalTree,     signalWeight );
    dataloader->AddBackgroundTree( backgroundTree, backgroundWeight );
    
    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    
    dataloader->PrepareTrainingAndTestTree( mycuts, mycuts,
                                          "nTrain_Signal=4500:nTrain_Background=4500:SplitMode=Random:NormMode=NumEvents:!V" );

    // Likelihood ("naive Bayes estimator")
    if (Use["Likelihood"])
        factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
    // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if (Use["MLP"])
        factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
    // Support Vector Machine
    if (Use["SVM"])
        factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
    // Adaptive Boost
    if (Use["BDT"])
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:UseRandomisedTrees=True" );
 
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();
    
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    double auc_bdt = factory->GetROCIntegral(dataloader, "BDT");
    double auc_svm = factory->GetROCIntegral(dataloader, "SVM");
    double auc_mlp = factory->GetROCIntegral(dataloader, "MLP");
    double auc_lik = factory->GetROCIntegral(dataloader, "Likelihood");

    std::cout << auc_bdt << "," << auc_svm << "," << auc_mlp << "," << auc_lik << std::endl;

//    std::ofstream ofile("bdt.csv",std::ios::app);
//
//    if (ofile.is_open()) {
//        ofile << auc_bdt << ",";
//        ofile.close();
//    }

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

