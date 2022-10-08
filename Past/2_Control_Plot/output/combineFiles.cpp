#include <vector>
#include <string>

void combineFiles()
{

	auto f = new TFile("data/signal.root");
	TDirectory* dir = (TDirectory*) f->Get("signal;1");
	TH1D* ljet_mass = static_cast<TH1D*>(dir->Get("J1_mass"));
	std::vector<double_t> J1_mass;
	//j1_mass->Draw();

	TFile* ff = new TFile("total.root","recreate");
	TTree* tt = new TTree("signal","signal");
	TBranch* br = new TBranch();
	char lfname("J1_mass");
	TLeaf* lf = new TLeaf(br,&lfname,"TH1D");
	//tt->Branch("J1_mass","TH1D",&j1_mass,32000,0);
	
	for (int i=0; i<ljet_mass->GetEntries(); i++) {
		J1_mass.push_back(ljet_mass->GetEntry(i)->GetValue());
	}
	ff->Write();

	f->Close();
	ff->Close();
}
