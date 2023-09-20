void readSomething() {

  TFile *fIn = new TFile;
  TTree *evt = new TTree;

  TH1F *h_jetpt = new TH1F("jetpt", "jetpt", 0, 100, 100);

  fIn = TFile::Open("F6B8F234-CC43-5D44-8D1E-727E33C38AA7.root");
  evt = (TTree*)fIn->Get("Events");

  TTreeReader evt_r(evt);
  TTreeReaderValue<double> jetpt(evt_r, "Jet_pt");

  std::vector<double> getSomething_;
  while(evt_r.Next()) { getSomething_.push_back(*jetpt); }

  for (int i = 0; i < 10; i++) { std::cout << getSomething_.at(i) << std::endl; }

  return 0;
}
