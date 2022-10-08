#include <vector>
#include <iostream>
#include <glob.h>
#include <unistd.h>

using namespace std;

int doubleHiggsAnalyzer_signal() {

	// Signal
	glob_t path;
	int retval;

	path.gl_pathc = 0;
	path.gl_pathv = NULL;
	path.gl_offs = 0;

	retval = glob("/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/64365B41-751A-7B4B-836B-DB4258390231.root", GLOB_NOCHECK | GLOB_NOSORT, NULL, &path);

	int ii(0); int hh(0); int ee(0);
	for (int idx=0; idx<path.gl_pathc; idx++) {
		//cout << path.gl_pathv[idx] << endl;
		
		if (ee > 5000000) break;
		TFile* f = new TFile(path.gl_pathv[idx],"read");
		auto t = f->Get<TTree>("Events");

		bool hlt_mumu_1; bool hlt_mumu_2; bool hlt_mumu_3; bool hlt_mumu_4;
		t->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&hlt_mumu_1);
		t->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&hlt_mumu_2);
		t->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",&hlt_mumu_3);
		t->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&hlt_mumu_4);

		bool flag_good; bool flag_halo; bool flag_hbhen; bool flag_hbiso;
		bool flag_dead; bool flag_badpf; bool flag_ecbad; bool flag_eebad;
		t->SetBranchAddress("Flag_goodVertices",&flag_good);
		t->SetBranchAddress("Flag_Flag_globalSuperTightHalo2016Filter",&flag_halo);
		t->SetBranchAddress("Flag_HBHENoiseFilter",&flag_hbhen);
		t->SetBranchAddress("Flag_HBHENoiseIsoFilter",&flag_hbiso);
		t->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&flag_dead);
		t->SetBranchAddress("Flag_BadPFMuonFilter",&flag_badpf);
		t->SetBranchAddress("Flag_ecalBadCalibFilter",&flag_ecbad);
		t->SetBranchAddress("Flag_eeBadScFilter",&flag_eebad);

		int event_num = 3000;
		//int event_num = t->GetEntries();
		cout << "Processing " << path.gl_pathv[idx] << "..." << endl;
		ii++;
		cout << "Have processed " << ee << " events..." << endl;

		vector<int> hlt_pass;
		for (int i=0; i<event_num; i++) {
			t->GetEntry(i);

			if (hlt_mumu_1 || hlt_mumu_2 || hlt_mumu_3 || hlt_mumu_4) {
				if (flag_good || flag_halo || flag_hbhen || flag_hbiso || flag_dead || flag_badpf || flag_ecbad || flag_eebad) {
					hlt_pass.push_back(i);
				}
			}
		} hh += sizeof(hh)/sizeof(int);
//		vector<int>::iterator it;
//		for (it=hlt_pass.begin(); it!=hlt_pass.end(); it++) {
//			cout << *it << endl;
//		}

		double H1_mass; double H1_pt; double H1_eta; double H1_phi;
		double H2_mass; double H2_pt; double H2_eta; double H2_phi;
		double HH_mass; double HH_pt; double HH_eta; double HH_phi;
		double jet1_mass; double jet1_pt; double jet1_eta; double jet1_phi; double jet1_btag;
		double jet2_mass; double jet2_pt; double jet2_eta; double jet2_phi; double jet2_btag;
		double lep1_mass; double lep1_pt; double lep1_eta; double lep1_phi; double lep1_charge;
		double lep2_mass; double lep2_pt; double lep2_eta; double lep2_phi; double lep2_charge;
		double met_et; double met_pt; double met_phi; double met_mass;
		double dr_jj; double dr_ll;

		for (int i=0; i<event_num; i++) {

			ee++;
			t->GetEntry(i);
			bool flag_hlt = false;
			vector<int>::iterator it;
			for (it=hlt_pass.begin(); it!=hlt_pass.end(); it++) {
				if (*it == i) flag_hlt = true; }
			if (flag_hlt != true) continue;

			TTreeReader reader("Events",f);
			TTreeReaderArray<Float_t> lep_mass(reader,"Muon_mass");
			Float_t temp_mass; double temp_pt; double temp_eta; double temp_phi;
			double temp_dxy; double temp_dxyerr; double temp_dz; double temp_dzerr;
			double temp_charge; double temp_medid; double temp_iso; double temp_pterr;
			double temp_charge_e; double temp_btag;
			int medid_n(0); int btag_n(0);
			int mu_pos(0); int mu_neg(0); int e_pos(0); int e_neg(0);

			bool flag_lep = false; bool flap_jet = false;
			UInt_t nMuon;

			t->SetBranchAddress("nMuon",&nMuon);

				cout << lep_mass.At(j) << endl;	
			}
		}
	}
	return 0;
}

		
