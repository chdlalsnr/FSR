from ROOT import *
from multiprocessing import Process, Queue
import numpy as np
import math
import glob
import sys
import random


def studyLNuSample():

        # ggF H->WW->qqlv Sample
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluHToWWToLNuQQ_M125_13TeV_powheg_JHUGenV628_pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*')

	ii,ee = 0,0

	w1_mass,w1_pt,w1_eta,w1_phi = [],[],[],[]
        w2_mass,w2_pt,w2_eta,w2_phi = [],[],[],[]
        w3_mass,w3_pt,w3_eta,w3_phi = [],[],[],[]

	for fname in flist:

		if (ee > 1000000): break
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


		for i in range(event_num):

                        ee += 1
                        t.GetEntry(i)
                        if i not in hlt_pass: continue
                        if not (t.nJet >= 3): continue

			flag_lep = False
			llep_mass,llep_pt,llep_eta,llep_phi = 0,0,0,0
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
			temp_dxy,temp_dz = [],[]

			for j in range(t.nMuon):

				temp_mass.append(t.Muon_mass.__getitem__(j))
                                temp_pt.append(t.Muon_pt.__getitem__(j))
                                temp_eta.append(t.Muon_eta.__getitem__(j))
                                temp_phi.append(t.Muon_phi.__getitem__(j))
				temp_dxy.append(t.Muon_dxy.__getitem__(j))
                                temp_dz.append(t.Muon_dz.__getitem__(j))


			for j in range(len(temp_mass)):
				if (temp_pt[j] > 25 and abs(temp_eta[j]) <= 2.4 and abs(temp_dxy[j]) <= 0.05 and abs(temp_dz[j]) <= 0.1):

					llep_mass = temp_mass[j]
					llep_pt = temp_pt[j]
					llep_eta = temp_eta[j]
					llep_phi = temp_phi[j]

					flag_lep = True


			temp_mass,temp_pt,temp_eta,temp_phi,temp_jetid = [],[],[],[],[]
                        for j in range(t.nJet): #(t.Jet_mass.__len__())
                                temp_mass.append(t.Jet_mass.__getitem__(j))
                                temp_pt.append(t.Jet_pt.__getitem__(j))
                                temp_eta.append(t.Jet_eta.__getitem__(j))
                                temp_phi.append(t.Jet_phi.__getitem__(j))
				temp_jetid.append(t.Jet_jetId.__getitem__(j))

			ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0
                        sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0
			for k in range(t.nJet):
                                if (temp_pt[k] > 20 and abs(temp_eta[k]) < 2.4 and temp_jetid[k] == 7): # jetid == 1: Loose, jetid == 2: Medium
                                        for l in range(k+1,t.nJet):
                                        	if (temp_pt[l] > 20 and abs(temp_eta[l]) < 2.4 and temp_jetid[k] == 7):

                                                        ljet_mass = temp_mass[k]
                                                        ljet_pt   = temp_pt[k]
                                                        ljet_eta  = temp_eta[k]
                                                        ljet_phi  = temp_phi[k]

                                                        sjet_mass = temp_mass[l]
                                                        sjet_pt   = temp_pt[l]
                                                        sjet_eta  = temp_eta[l]
                                                        sjet_phi  = temp_phi[l]

						break
				break

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()

			jet1_p4.SetPtEtaPhiM(ljet_pt,ljet_eta,ljet_phi,ljet_mass)
			jet2_p4.SetPtEtaPhiM(sjet_pt,sjet_eta,sjet_phi,sjet_mass)

			w_p4 = jet1_p4 + jet2_p4

			if (w_p4.M() == 0): continue

			w1_mass.append(w_p4.M())
			w1_pt.append(w_p4.Pt())
			w1_eta.append(w_p4.Eta())
			w1_phi.append(w_p4.Phi())


#			lep_p4 = TLorentzVector()
#			lep_p4.SetPtEtaPhiM(llep_pt,llep_eta,llep_phi,llep_mass)
#
#			met_pt = t.MET_pt
#			met_phi = t.MET_phi
#
#			met_p1 = reconstructWBoson1(met_pt,met_phi)
#			met_p2 = reconstructWBoson2(met_pt,met_phi)
#			w3_p4 = reconstructWBoson3(met_pt,met_phi,lep_p4)
#
#			w1_p4 = lep_p4 + met_p1
#			w2_p4 = lep_p4 + met_p2
#
#			w1_mass.append(w1_p4.M())
#			w1_pt.append(w1_p4.Pt())
#			w1_eta.append(w1_p4.Eta())
#			w1_phi.append(w1_p4.Phi())
#
#			w2_mass.append(w2_p4.M())
#			w2_pt.append(w2_p4.Pt())
#			w2_eta.append(w2_p4.Eta())
#			w2_phi.append(w2_p4.Phi())
#
#			w3_mass.append(w3_p4.M())
#			w3_pt.append(w3_p4.Pt())
#			w3_eta.append(w3_p4.Eta())
#			w3_phi.append(w3_p4.Phi())


	c = TCanvas("c1","Kinematic Variable Distributions of W_qq", 900,600)
	c.Divide(2,2)

        hist1_mass = TH1D("h1","W_qq Boson Mass Distribution", 40,0,200)
        hist1_pt   = TH1D("h2","W_qq Boson Transeverse Momentum Distribution",40,0,400)
        hist1_eta  = TH1D("h3","W_qq Boson Eta Distribution",40,-2.4,2.4)
        hist1_phi  = TH1D("h4","W_qq Boson Phi Distribution",40,-3.14,3.14)

        for i in range(len(w1_pt)):

                hist1_mass.Fill(w1_mass[i])
                hist1_pt.Fill(w1_pt[i])
                hist1_eta.Fill(w1_eta[i])
                hist1_phi.Fill(w1_phi[i])

	c.cd(1)
	hist1_mass.Draw()

	c.cd(2)
	hist1_pt.Draw()

	c.cd(3)
	hist1_eta.Draw()

	c.cd(4)
	hist1_phi.Draw()

        c.SaveAs("kinematics_wqq.png")

        input("Press Enter to continue...")



def reconstructWBoson1(met_pt,met_phi):

	met_p4 = TLorentzVector()

	met_p4.SetPtEtaPhiM(met_pt,0,met_phi,0)

	return met_p4

def reconstructWBoson2(met_pt,met_phi):

	met_p4 = TLorentzVector()

	eta = random.uniform(-2.4,2.4) 
	met_p4.SetPtEtaPhiM(met_pt,eta,met_phi,0)

	return met_p4

def reconstructWBoson3(met_pt,met_phi,lep_p4):

	met_p4 = TLorentzVector()

        w_cand = []
        w_mass = 80.379
        eta_cand = np.arange(-2.4,2.4,0.4)
        for j in range(len(eta_cand)):

                met_p4.SetPtEtaPhiM(met_pt,eta_cand[j],met_phi,0)
                w_p4 = lep_p4 + met_p4
                w_cand.append(w_p4.M())

        w_res = []
	for j in range(len(w_cand)): w_res.append(abs(w_cand[j]-w_mass))

	w_min = min(w_res)

        jj = 0
        for j in range(len(w_cand)):
                if (abs(w_cand[j]-w_mass) == w_min):
                        jj = j

        met_p4.SetPtEtaPhiM(met_pt,eta_cand[jj],met_phi,0)

	w3_p4 = lep_p4 + met_p4

	return w3_p4


if __name__ == "__main__":
	studyLNuSample()
