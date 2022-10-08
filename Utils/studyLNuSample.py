from ROOT import *
from multiprocessing import Process, Queue
import numpy as np
import math
import glob
import sys


def studyLNuSample():

        # W->LNu Sample
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/260000/*')

	ii,ee = 0,0
	for fname in flist:

                f = TFile(fname,"read")
                t = f.Get("Events")

                event_num = t.GetEntriesFast()
                print("Processing %s..."%(fname))
                ii = ii+1
                print("Have processed %i events..."%(ee))

                # HLT pass
                hlt_pass = []
                for i in range(event_num):
                        t.GetEntry(i)
                        hlt_mumu_1 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
                        hlt_mumu_2 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
                        hlt_mumu_3 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
                        hlt_mumu_4 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ

                        hlt_doubleMu = False
                        if (hlt_mumu_1 == 1 or hlt_mumu_2 == 1 or hlt_mumu_3 == 1 or hlt_mumu_4 == 1): hlt_doubleMu = True

                        hlt_mu_1 = t.HLT_IsoMu22
                        hlt_mu_2 = t.HLT_IsoTkMu22
                        hlt_mu_3 = t.HLT_IsoMu22_eta2p1
                        hlt_mu_4 = t.HLT_IsoTkMu22_eta2p1
                        hlt_mu_5 = t.HLT_IsoMu24
                        hlt_mu_6 = t.HLT_IsoTkMu24

                        hlt_singleMu = False
                        if (hlt_mu_1 == 1 or hlt_mu_2 == 1 or hlt_mu_3 == 1 or hlt_mu_4 == 1 or hlt_mu_5 == 1 or hlt_mu_6 == 1): hlt_singleMu = True

                        flag_good  = t.Flag_goodVertices
                        flag_halo  = t.Flag_globalSuperTightHalo2016Filter
                        flag_hbhen = t.Flag_HBHENoiseFilter
                        flag_hbiso = t.Flag_HBHENoiseIsoFilter
                        flag_dead  = t.Flag_EcalDeadCellTriggerPrimitiveFilter
                        flag_badpf = t.Flag_BadPFMuonFilter
                        flag_ecbad = t.Flag_ecalBadCalibFilter
                        #flag_eebad = t.Flag_eeBadScFilter

                        hlt_vertex = False
                        if (flag_good == 1 and flag_halo == 1 and flag_hbhen == 1 and flag_hbiso == 1 and flag_dead == 1 and flag_badpf == 1 and flag_ecbad == 1):
                                hlt_vertex = True

                        if (hlt_doubleMu or hlt_singleMu):
                                if hlt_vertex:
                                        hlt_pass.append(i)

		mass,pt,eta,phi = [],[],[],[]
		for i in range(event_num):

                	if (ee > 10000): break
                        ee += 1
                        t.GetEntry(i)
                        if i not in hlt_pass: continue
                        if not (t.nJet >= 3 or (t.nJet >= 1 and t.nFatJet >= 1)): continue

			flag_lep = False
			llep_mass,llep_pt,llep_eta,llep_phi = 0,0,0,0
			temp_mass,temp_pt,temp_eta,temp_phi,temp_medid = [],[],[],[],[]

			for j in range(t.nMuon):

				temp_mass.append(t.Muon_mass.__getitem__(j))
                                temp_pt.append(t.Muon_pt.__getitem__(j))
                                temp_eta.append(t.Muon_eta.__getitem__(j))
                                temp_phi.append(t.Muon_phi.__getitem__(j))
				temp_medid.append(t.Muon_mediumId.__getitem__(j))

			for j in range(len(temp_mass)):
				if (temp_pt[j] > 25 and abs(temp_eta[j]) <= 2.4 and temp_medid[j] == True):

					llep_mass = temp_mass[j]
					llep_pt = temp_pt[j]
					llep_eta = temp_eta[j]
					llep_phi = temp_phi[j]

					flag_lep = True

			lep_p4 = TLorentzVector()
			lep_p4.SetPtEtaPhiM(llep_pt,llep_eta,llep_phi,llep_mass)

			met_et = t.MET_sumEt
			met_pt = t.MET_pt
			met_phi = t.MET_phi

			met_p4 = TLorentzVector()

			mass_cand = np.arange(0,1.0,0.1)
			eta_cand = np.arange(-3.4,3.4,0.4)

			w_mass = 80.379
			w_cand = []
			for j in range(len(eta_cand)):
				temp = []
				for k in range(len(mass_cand)):
		
					met_p4.SetPtEtaPhiM(met_pt,eta_cand[j],met_phi,mass_cand[k])
					w_p4 = lep_p4 + met_p4

					temp.append(abs(w_p4.M()-w_mass))
				w_cand.append(temp)

			w_can = min(w_cand)
			w_min = min(w_can)

			jj,kk = 0,0
			nu_p4 = TLorentzVector()
			for j in range(len(w_cand)):
				for k in range(len(w_cand[j])):
					if (abs(w_cand[j][k]-w_mass) == w_min):

						jj = j
						kk = k

					break
				break

			nu_p4.SetPtEtaPhiE(met_pt,eta_cand[jj],met_phi,mass_cand[kk])

			mass.append(nu_p4.M())
			pt.append(nu_p4.Pt())
			eta.append(nu_p4.Eta())
			phi.append(nu_p4.Phi())


	hist_mass = TH1D("h1","Neutrino Invariant Mass Distribution",40,0,1)
        hist_pt = TH1D("h2","Neutrino Transverse Momentum Distribution",40,0,100)
        hist_eta = TH1D("h3","Neutrino Eta Distribution",40,-2.4,2.4)
        hist_phi = TH1D("h4","Neutrino Phi Distribution",40,-3.14,3.14)

	for i in range(len(mass)):

                hist_mass.Fill(mass[i])
		hist_pt.Fill(pt[i])
		hist_eta.Fill(eta[i])
		hist_phi.Fill(phi[i])


	c = TCanvas("c","Kinematic Variables of RECO Nu", 900,600)
	c.Divide(2,2)
	c.cd(1)

	hist_mass.GetXaxis().SetTitle("Invariant Mass[GeV]")
	hist_mass.GetYaxis().SetTitle("Entries")
	hist_mass.Draw()

	c.cd(2)
	hist_pt.GetXaxis().SetTitle("Transverse Momentum[GeV]")
        hist_pt.GetYaxis().SetTitle("Entries")
        hist_pt.Draw()

	c.cd(3)
	hist_eta.GetXaxis().SetTitle("Eta")
        hist_eta.GetYaxis().SetTitle("Entries")
        hist_eta.Draw()

	c.cd(4)
	hist_phi.GetXaxis().SetTitle("Phi")
        hist_phi.GetYaxis().SetTitle("Entries")
        hist_phi.Draw()

	c.SetGrid()
	c.Draw()
	c.SaveAs("kinematics_nu.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
	studyLNuSample()
