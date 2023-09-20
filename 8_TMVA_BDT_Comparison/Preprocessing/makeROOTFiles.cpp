//cpp
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
//#include <boost/lexical_cast.hpp>


using namespace std;

void makeROOTFiles() {

	fstream fs;
	string dr = "background";
	fs.open("background.csv",ios::in);
	string line;

	std::vector<Float_t> jet1_mass,jet1_pt,jet1_eta,jet1_phi; 
	std::vector<Float_t> jet2_mass,jet2_pt,jet2_eta,jet2_phi;
	std::vector<Float_t> lep1_mass,lep1_pt,lep1_eta,lep1_phi;
	std::vector<Float_t> lep2_mass,lep2_pt,lep2_eta,lep2_phi;
	std::vector<Float_t> delta_ll,delta_jj,Met_pt,Met_phi;
	std::vector<Float_t> higgs1_mass,higgs1_pt,higgs1_eta,higgs1_phi;
	std::vector<Float_t> higgs2_mass,higgs2_pt,higgs2_eta,higgs2_phi;
	std::vector<Float_t> Higgsness,Topness;

	int cc = 0;
	while(cc < 15001) {

		getline(fs,line);

		istringstream ss(line);
		string buffer;
		std::vector<string> data;

		while(getline(ss,buffer,',')) { data.push_back(buffer); }

		if (cc > 0) {

			Float_t jet1_mass_f,jet1_pt_f,jet1_eta_f,jet1_phi_f;
			Float_t jet2_mass_f,jet2_pt_f,jet2_eta_f,jet2_phi_f;
			Float_t lep1_mass_f,lep1_pt_f,lep1_eta_f,lep1_phi_f;
			Float_t lep2_mass_f,lep2_pt_f,lep2_eta_f,lep2_phi_f;
			Float_t delta_ll_f,delta_jj_f,met_pt_f,met_phi_f;
			Float_t higgs1_mass_f,higgs1_pt_f,higgs1_eta_f,higgs1_phi_f;
			Float_t higgs2_mass_f,higgs2_pt_f,higgs2_eta_f,higgs2_phi_f;
			Float_t higgsness_f,topness_f;

			jet1_mass_f = stof(data.at(0));
			jet1_pt_f   = stof(data.at(1));
			jet1_eta_f  = stof(data.at(2));
			jet1_phi_f  = stof(data.at(3));

			jet2_mass_f = stof(data.at(4));
			jet2_pt_f   = stof(data.at(5));
			jet2_eta_f  = stof(data.at(6));
			jet2_phi_f  = stof(data.at(7));

			lep1_mass_f = stof(data.at(8));
			lep1_pt_f   = stof(data.at(9));
			lep1_eta_f  = stof(data.at(10));
			lep1_phi_f  = stof(data.at(11));

			lep2_mass_f = stof(data.at(12));
			lep2_pt_f   = stof(data.at(13));
			lep2_eta_f  = stof(data.at(14));
			lep2_phi_f  = stof(data.at(15));

			delta_jj_f = stof(data.at(16));
			delta_ll_f = stof(data.at(17));
			met_pt_f   = stof(data.at(18));
			met_phi_f  = stof(data.at(19));

			higgs1_mass_f= stof(data.at(20));
			higgs1_pt_f  = stof(data.at(21));
			higgs1_eta_f = stof(data.at(22));
			higgs1_phi_f = stof(data.at(23));

			higgs2_mass_f= stof(data.at(24));
			higgs2_pt_f  = stof(data.at(25));
			higgs2_eta_f = stof(data.at(26));
			higgs2_phi_f = stof(data.at(27));

			higgsness_f  = stof(data.at(28));
			topness_f    = stof(data.at(29));

			jet1_mass.push_back(jet1_mass_f);
			jet1_pt.push_back(jet1_pt_f);
			jet1_eta.push_back(jet1_eta_f);
			jet1_phi.push_back(jet1_phi_f);

			jet2_mass.push_back(jet2_mass_f);
			jet2_pt.push_back(jet2_pt_f);
			jet2_eta.push_back(jet2_eta_f);
			jet2_phi.push_back(jet2_phi_f);

			lep1_mass.push_back(lep1_mass_f);
			lep1_pt.push_back(lep1_pt_f);
			lep1_eta.push_back(lep1_eta_f);
			lep1_phi.push_back(lep1_phi_f);

			lep2_mass.push_back(lep2_mass_f);
			lep2_pt.push_back(lep2_pt_f);
			lep2_eta.push_back(lep2_eta_f);
			lep2_phi.push_back(lep2_phi_f);
			
			delta_ll.push_back(delta_ll_f);
			delta_jj.push_back(delta_jj_f);
			Met_pt.push_back(met_pt_f);
			Met_phi.push_back(met_phi_f);

			higgs1_mass.push_back(higgs1_mass_f);
			higgs1_pt.push_back(higgs1_pt_f);
			higgs1_eta.push_back(higgs1_eta_f);
			higgs1_phi.push_back(higgs1_phi_f);

			higgs2_mass.push_back(higgs2_mass_f);
			higgs2_pt.push_back(higgs2_pt_f);
			higgs2_eta.push_back(higgs2_eta_f);
			higgs2_phi.push_back(higgs2_phi_f);

			Higgsness.push_back(higgsness_f);
			Topness.push_back(topness_f);
		}
		data.clear();
		cc++;
	}

//	for (vector<Float_t>::iterator it=h1_mass.begin()+1; it!=h1_mass.end(); it++)
//	{ cout << *it << endl; }

	TFile* f = TFile::Open("background.root","recreate");
	TTree* tt = new TTree("background","background");

	Float_t j1_mass,j1_pt,j1_eta,j1_phi;
        Float_t j2_mass,j2_pt,j2_eta,j2_phi;
        Float_t l1_mass,l1_pt,l1_eta,l1_phi;
        Float_t l2_mass,l2_pt,l2_eta,l2_phi;
        Float_t dr_ll,dr_jj,met_pt,met_phi;
        Float_t h1_mass,h1_pt,h1_eta,h1_phi;
        Float_t h2_mass,h2_pt,h2_eta,h2_phi;
        Float_t higgsness,topness;

	tt->Branch("jet1_mass",&j1_mass,"j1_mass/F");
	tt->Branch("jet1_pt",&j1_pt,"j1_pt/F");
	tt->Branch("jet1_eta",&j1_eta,"j1_eta/F");
	tt->Branch("jet1_phi",&j1_phi,"j1_phi/F");

	tt->Branch("jet2_mass",&j2_mass,"j2_mass/F");
	tt->Branch("jet2_pt",&j2_pt,"j2_pt/F");
	tt->Branch("jet2_eta",&j2_eta,"j2_eta/F");
	tt->Branch("jet2_phi",&j2_phi,"j2_phi/F");

	tt->Branch("lep1_mass",&l1_mass,"l1_mass/F");
	tt->Branch("lep1_pt",&l1_pt,"l1_pt/F");
	tt->Branch("lep1_eta",&l1_eta,"l1_eta/F");
	tt->Branch("lep1_phi",&l1_phi,"l1_phi/F");

	tt->Branch("lep2_mass",&l2_mass,"l2_mass/F");
	tt->Branch("lep2_pt",&l2_pt,"l2_pt/F");
	tt->Branch("lep2_eta",&l2_eta,"l2_eta/F");
	tt->Branch("lep2_phi",&l2_phi,"l2_phi/F");

	tt->Branch("dr_ll",&dr_ll,"dr_ll/F");
	tt->Branch("dr_jj",&dr_jj,"dr_jj/F");
	tt->Branch("met_pt",&met_pt,"met_pt/F");
	tt->Branch("met_phi",&met_phi,"met_phi/F");

	tt->Branch("h1_mass",&h1_mass,"h1_mass/F");
	tt->Branch("h1_pt",&h1_pt,"h1_pt/F");
	tt->Branch("h1_eta",&h1_eta,"h1_eta/F");
	tt->Branch("h1_phi",&h1_phi,"h1_phi/F");

	tt->Branch("h2_mass",&h2_mass,"h2_mass/F");
	tt->Branch("h2_pt",&h2_pt,"h2_pt/F");
	tt->Branch("h2_eta",&h2_eta,"h2_eta/F");
	tt->Branch("h2_phi",&h2_phi,"h2_phi/F");

	tt->Branch("higgsness",&higgsness,"higgsness/F");
	tt->Branch("topness",&topness,"hopness/F");

	for (int i=0; i<cc-1; i++) {

		j1_mass = jet1_mass.at(i);
		j1_pt   = jet1_pt.at(i);
		j1_eta  = jet1_eta.at(i);
		j1_phi  = jet1_phi.at(i);

		j2_mass = jet2_mass.at(i);
                j2_pt   = jet2_pt.at(i);
                j2_eta  = jet2_eta.at(i);
                j2_phi  = jet2_phi.at(i);

		l1_mass = lep1_mass.at(i);
                l1_pt   = lep1_pt.at(i);
                l1_eta  = lep1_eta.at(i);
                l1_phi  = lep1_phi.at(i);

                l2_mass = lep2_mass.at(i);
                l2_pt   = lep2_pt.at(i);
                l2_eta  = lep2_eta.at(i);
                l2_phi  = lep2_phi.at(i);

		dr_ll   = delta_ll.at(i);
		dr_jj   = delta_jj.at(i);
		met_pt  = Met_pt.at(i);
		met_phi = Met_phi.at(i);

		h1_mass = higgs1_mass.at(i);
		h1_pt   = higgs1_pt.at(i);
		h1_eta  = higgs1_eta.at(i);
		h1_phi  = higgs1_phi.at(i);

		h2_mass = higgs2_mass.at(i);
                h2_pt   = higgs2_pt.at(i);
                h2_eta  = higgs2_eta.at(i);
                h2_phi  = higgs2_phi.at(i);

		higgsness = Higgsness.at(i);
		topness = Topness.at(i);

		tt->Fill();
	}
	f->cd();
	tt->Write();
	f->Close();

	std::cout << "Successfully generated root file..." << std::endl;
}
