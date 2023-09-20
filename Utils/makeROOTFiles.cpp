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
	fs.open("tt_2l2nu.csv",ios::in);
	string line;

	std::vector<double> jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag; 
	std::vector<double> jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag; 
	std::vector<double> lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge; 
	std::vector<double> lep2_mass,lep2_pt,lep2_eta,lep2_phi,lep2_charge;

	std::vector<double> higgs1_mass,higgs1_pt,higgs1_eta,higgs1_phi;
	std::vector<double> higgs2_mass,higgs2_pt,higgs2_eta,higgs2_phi;
	std::vector<double> hhiggs_mass,hhiggs_pt,hhiggs_eta,hhiggs_phi;
	std::vector<double> deltar_ll,deltar_jj;

	int cc = 0;
	//while(!fs.eof()) {
	//while(cc < 2051) {
	//while(cc < 4001) {
	while(cc < 10000) {

		getline(fs,line);

		istringstream ss(line);
		string buffer;
		std::vector<string> data;

		while(getline(ss,buffer,',')) { data.push_back(buffer); }

		//if (cc > 0) { cout << stod(data.at(0)) << endl; }
		if (cc > 0) {

			double jet1_mass_d,jet1_pt_d,jet1_eta_d,jet1_phi_d,jet1_btag_d;
			double jet2_mass_d,jet2_pt_d,jet2_eta_d,jet2_phi_d,jet2_btag_d;
			double lep1_mass_d,lep1_pt_d,lep1_eta_d,lep1_phi_d;
			double lep2_mass_d,lep2_pt_d,lep2_eta_d,lep2_phi_d;

			double h1_mass_d,h1_pt_d,h1_eta_d,h1_phi_d;
			double h2_mass_d,h2_pt_d,h2_eta_d,h2_phi_d;
			double hh_mass_d,hh_pt_d,hh_eta_d,hh_phi_d;
			double dr_ll_d,dr_jj_d;

			int lep1_charge_i,lep2_charge_i;

			jet1_mass_d = stod(data.at(0));
			jet1_pt_d = stod(data.at(1));
			jet1_eta_d = stod(data.at(2));
			jet1_phi_d = stod(data.at(3));
			jet1_btag_d = stod(data.at(4));

			jet2_mass_d = stod(data.at(5));
			jet2_pt_d = stod(data.at(6));
			jet2_eta_d = stod(data.at(7));
			jet2_phi_d = stod(data.at(8));
			jet2_btag_d = stod(data.at(9));

			lep1_mass_d = stod(data.at(10));
			lep1_pt_d = stod(data.at(11));
			lep1_eta_d = stod(data.at(12));
			lep1_phi_d = stod(data.at(13));
			lep1_charge_i = stoi(data.at(14));

			lep2_mass_d = stod(data.at(15));
			lep2_pt_d = stod(data.at(16));
			lep2_eta_d = stod(data.at(17));
			lep2_phi_d = stod(data.at(18));
			lep2_charge_i = stoi(data.at(19));

			dr_jj_d = stod(data.at(20));
			dr_ll_d = stod(data.at(21));

			h1_mass_d = stod(data.at(22));
			h1_pt_d = stod(data.at(23));
			h1_eta_d = stod(data.at(24));
			h1_phi_d = stod(data.at(25));

			h2_mass_d = stod(data.at(26));
			h2_pt_d = stod(data.at(27));
			h2_eta_d = stod(data.at(28));
			h2_phi_d = stod(data.at(29));

			hh_mass_d = stod(data.at(30));
			hh_pt_d = stod(data.at(31));
			hh_eta_d = stod(data.at(32));
			hh_phi_d = stod(data.at(33));

			jet1_mass.push_back(jet1_mass_d);
			jet1_pt.push_back(jet1_pt_d);
			jet1_eta.push_back(jet1_eta_d);
			jet1_phi.push_back(jet1_phi_d);
			jet1_btag.push_back(jet1_btag_d);

			jet2_mass.push_back(jet2_mass_d);
			jet2_pt.push_back(jet2_pt_d);
			jet2_eta.push_back(jet2_eta_d);
			jet2_phi.push_back(jet2_phi_d);
			jet2_btag.push_back(jet2_btag_d);

			lep1_mass.push_back(lep1_mass_d);
			lep1_pt.push_back(lep1_pt_d);
			lep1_eta.push_back(lep1_eta_d);
			lep1_phi.push_back(lep1_phi_d);
			lep1_charge.push_back(lep1_charge_i);

			lep2_mass.push_back(lep2_mass_d);
			lep2_pt.push_back(lep2_pt_d);
			lep2_eta.push_back(lep2_eta_d);
			lep2_phi.push_back(lep2_phi_d);
			lep2_charge.push_back(lep2_charge_i);

			deltar_jj.push_back(dr_jj_d);
			deltar_ll.push_back(dr_ll_d);

			higgs1_mass.push_back(h1_mass_d);
			higgs1_pt.push_back(h1_pt_d);
			higgs1_eta.push_back(h1_eta_d);
			higgs1_phi.push_back(h1_phi_d);

			higgs2_mass.push_back(h2_mass_d);
			higgs2_pt.push_back(h2_pt_d);
			higgs2_eta.push_back(h2_eta_d);
			higgs2_phi.push_back(h2_phi_d);

			hhiggs_mass.push_back(hh_mass_d);
			hhiggs_pt.push_back(hh_pt_d);
			hhiggs_eta.push_back(hh_eta_d);
			hhiggs_phi.push_back(hh_phi_d);
		}
		data.clear();
		cc++;
	}

//	for (vector<double>::iterator it=h1_mass.begin()+1; it!=h1_mass.end(); it++)
//	{ cout << *it << endl; }

	TFile* f = TFile::Open("tmva_test.root","recreate");
	TTree* test = new TTree("test","test tree");
	for (int i=0; i<cc-1; i++) {

		double j1_mass,j1_pt,j1_eta,j1_phi,j1_btag;
		double j2_mass,j2_pt,j2_eta,j2_phi,j2_btag;
		double l1_mass,l1_pt,l1_eta,l1_phi;
		double l2_mass,l2_pt,l2_eta,l2_phi;
		double dr_ll,dr_jj;
		double h1_mass,h1_pt,h1_eta,h1_phi;
		double h2_mass,h2_pt,h2_eta,h2_phi;
		double hh_mass,hh_pt,hh_eta,hh_phi;
		int l1_charge,l2_charge;

		test->Branch("J1_mass", &j1_mass, "j1_mass/D");
		test->Branch("J1_pt", &j1_pt, "j1_pt/D");
		test->Branch("J1_eta", &j1_eta, "j1_eta/D");
		test->Branch("J1_phi", &j1_phi, "j1_phi/D");
		test->Branch("J1_btag", &j1_btag, "j1_btag/D");

		test->Branch("J2_mass", &j2_mass, "j2_mass/D");
                test->Branch("J2_pt", &j2_pt, "j2_pt/D");
                test->Branch("J2_eta", &j2_eta, "j2_eta/D");
                test->Branch("J2_phi", &j2_phi, "j2_phi/D");
                test->Branch("J2_btag", &j2_btag, "j2_btag/D");

		test->Branch("L1_mass", &l1_mass, "l1_mass/D");
                test->Branch("L1_pt", &l1_pt, "l1_pt/D");
                test->Branch("L1_eta", &l1_eta, "l1_eta/D");
                test->Branch("L1_phi", &l1_phi, "l1_phi/D");
                test->Branch("L1_charge", &l1_charge, "l1_charge/I");

                test->Branch("L2_mass", &l2_mass, "l2_mass/D");
                test->Branch("L2_pt", &l2_pt, "l2_pt/D");
                test->Branch("L2_eta", &l2_eta, "l2_eta/D");
                test->Branch("L2_phi", &l2_phi, "l2_phi/D");
                test->Branch("L2_charge", &l2_charge, "l2_charge/I");

		test->Branch("DeltaR_ll", &dr_ll, "dr_ll/D");
		test->Branch("DeltaR_jj", &dr_jj, "dr_jj/D");

		test->Branch("H1_mass", &h1_mass, "h1_mass/D");
		test->Branch("H1_pt", &h1_pt, "h1_pt/D");
		test->Branch("H1_eta", &h1_eta, "h1_eta/D");
		test->Branch("H1_phi", &h1_phi, "h1_phi/D");

		test->Branch("H2_mass", &h2_mass, "h2_mass/D");
                test->Branch("H2_pt", &h2_pt, "h2_pt/D");
                test->Branch("H2_eta", &h2_eta, "h2_eta/D");
                test->Branch("H2_phi", &h2_phi, "h2_phi/D");

		test->Branch("HH_mass", &hh_mass, "hh_mass/D");
                test->Branch("HH_pt", &hh_pt, "hh_pt/D");
                test->Branch("HH_eta", &hh_eta, "hh_eta/D");
                test->Branch("HH_phi", &hh_phi, "hh_phi/D");

		j1_mass = jet1_mass.at(i);
		j1_pt   = jet1_pt.at(i);
		j1_eta  = jet1_eta.at(i);
		j1_phi  = jet1_phi.at(i);
		j1_btag = jet1_btag.at(i);

		j2_mass = jet2_mass.at(i);
                j2_pt   = jet2_pt.at(i);
                j2_eta  = jet2_eta.at(i);
                j2_phi  = jet2_phi.at(i);
                j2_btag = jet2_btag.at(i);

		l1_mass = lep1_mass.at(i);
                l1_pt   = lep1_pt.at(i);
                l1_eta  = lep1_eta.at(i);
                l1_phi  = lep1_phi.at(i);
                l1_charge = lep1_charge.at(i);

                l2_mass = lep2_mass.at(i);
                l2_pt   = lep2_pt.at(i);
                l2_eta  = lep2_eta.at(i);
                l2_phi  = lep2_phi.at(i);
                l2_charge = lep2_charge.at(i);

		dr_ll   = deltar_ll.at(i);
		dr_jj   = deltar_jj.at(i);

		h1_mass = higgs1_mass.at(i);
		h1_pt   = higgs1_pt.at(i);
		h1_eta  = higgs1_eta.at(i);
		h1_phi  = higgs1_phi.at(i);

		h2_mass = higgs2_mass.at(i);
                h2_pt   = higgs2_pt.at(i);
                h2_eta  = higgs2_eta.at(i);
                h2_phi  = higgs2_phi.at(i);

		hh_mass = hhiggs_mass.at(i);
                hh_pt   = hhiggs_pt.at(i);
                hh_eta  = hhiggs_eta.at(i);
                hh_phi  = hhiggs_phi.at(i);

		test->Fill();
	}
	f->cd();
	test->Write();
	f->Close();

	std::cout << "Successfully generated root file..." << std::endl;
}
