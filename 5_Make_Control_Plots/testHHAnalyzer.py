from ROOT import *
from multiprocessing import Process, Queue
import numpy as np
import math
import pandas as pd
import json
import sys
import pickle
import glob


# NOTE: Only considered WW->MuMu SR
def main():

	fname = "/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/64365B41-751A-7B4B-836B-DB4258390231.root"
#	fname = "/xrootd/store/data/Run2016B_ver2/DoubleMuon/NANOAOD/Nano25Oct2019_ver2-v1/100000/4831921A-213B-3541-B295-48486A2E8B9F.root"

	f = TFile.Open(fname,"read")
	t = f.Get("Events")
	r = f.Get("Runs")

	event_num = 1000
#        event_num = t.GetEntriesFast()

        # Only MuMu channel considered at this point
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
		flag_eebad = t.Flag_eeBadScFilter

		met_filter = False
		if (flag_good == 1 and flag_halo == 1 and flag_hbhen == 1 and flag_hbiso == 1 and flag_dead == 1 and flag_badpf == 1 and flag_ecbad == 1):
			met_filter = True

		if (hlt_doubleMu or hlt_singleMu):
			if met_filter: 
				hlt_pass.append(i)

	# =================================================================
	#Remained cut : The lepton isolation, defined as the scalar
	#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
	#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
	# Medium muon discrimination
	for i in range(event_num):

		t.GetEntry(i)
		if i not in hlt_pass: continue
		if not (t.nJet >= 2): continue
		if not (t.nMuon >= 2): continue

		llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
		slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0

		ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
		sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0

		temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
		temp_charge,temp_medid = [],[]
                temp_dxy,temp_dz,temp_btag = [],[],[]
		temp_iso,temp_sip3d,temp_mva = [],[],[]

		flag_lep,flag_jet = False,False
		for j in range(t.nMuon):
			temp_mass.append(t.Muon_mass.__getitem__(j))
                        temp_pt.append(t.Muon_pt.__getitem__(j))
                        temp_eta.append(t.Muon_eta.__getitem__(j))
                        temp_phi.append(t.Muon_phi.__getitem__(j))
                        temp_dxy.append(t.Muon_dxy.__getitem__(j))
                        temp_dz.append(t.Muon_dz.__getitem__(j))
                        temp_charge.append(t.Muon_charge.__getitem__(j))
                        temp_medid.append(t.Muon_mediumId.__getitem__(j))
			temp_iso.append(t.Muon_miniPFRelIso_all.__getitem__(j))
			temp_sip3d.append(t.Muon_sip3d.__getitem__(j))
			temp_mva.append(t.Muon_mvaTTH.__getitem__(j))

		print(t.HLT_LooseIsoPFTau50_Trk30_eta2p1)

		l1,l2 = -1,-1
		for k in range(t.nMuon):
                        if (flag_lep == True): break
			if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_iso[k] < 0.4 and temp_sip3d[k] < 8 and temp_mva[k] > 0.5 and temp_medid[k] == True):
                        	for l in range(k+1,t.nMuon):
					if (flag_lep == True): break
                                	if (temp_pt[l] > 15 and abs(temp_eta[l]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_iso[l] < 0.4 and temp_sip3d[l] < 8 and temp_mva[l] > 0.5 and temp_medid[l] == True):

                                        	if (temp_charge[k] != temp_charge[l]):
							llep_mass = temp_mass[k]
							llep_pt   = temp_pt[k]
							llep_eta  = temp_eta[k]
							llep_phi  = temp_phi[k]
							llep_charge = temp_charge[k]

							slep_mass = temp_mass[l]
							slep_pt   = temp_pt[l]
							slep_eta  = temp_eta[l]
							slep_phi  = temp_phi[l]
							slep_charge = temp_charge[l]
							if (llep_charge == slep_charge): break

							lep1_p4 = TLorentzVector()
							lep2_p4 = TLorentzVector()
							lep1_p4.SetPtEtaPhiM(llep_pt,llep_eta,llep_phi,llep_mass)
							lep2_p4.SetPtEtaPhiM(slep_pt,slep_eta,slep_phi,slep_mass)
							dilep_p4 = lep1_p4 + lep2_p4

							l1 = k
							l2 = l

							if not (dilep_p4.M() > 12): continue
							if (dilep_p4.M() > 12): flag_lep = True


		# =================================================================
		# Jet discrimination
		temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                temp_btag,temp_jetid = [],[]
		temp_rho = t.fixedGridRhoFastjetAll

		temp_jet = []
		for j in range(t.nJet): #(t.Jet_mass.__len__())

			cont = []
			cont.append(t.Jet_mass.__getitem__(j))
			cont.append(t.Jet_pt.__getitem__(j))
			cont.append(t.Jet_eta.__getitem__(j))
			cont.append(t.Jet_phi.__getitem__(j))
			cont.append(t.Jet_btagDeepFlavB.__getitem__(j))
			cont.append(t.Jet_jetId.__getitem__(j))
			cont.append(t.Jet_area.__getitem__(j))

			temp_jet.append(cont)

		# JEC L1Pileup correction
		for j in range(len(temp_jet)):

			temp_jet[j][1] = correctorL1Pileup(temp_jet[j][2],temp_jet[j][1],temp_rho,temp_jet[j][-1])


		jplist = []
		for j in range(len(temp_jet)): jplist.append(temp_jet[j][1])

		jplist.sort(reverse=True)

		j = 0
		sort_jet = []
		while len(jplist) > 0:

			if (jplist[0] != temp_jet[j][1]):
				j += 1
				continue

			if (jplist[0] == temp_jet[j][1]):
				sort_jet.append(temp_jet[j])
				jplist.pop(0)
				temp_jet.pop(j)
				j = 0

		for j in range(t.nJet):
			temp_mass.append(sort_jet[j][0])
			temp_pt.append(sort_jet[j][1])
			temp_eta.append(sort_jet[j][2])
			temp_phi.append(sort_jet[j][3])
			temp_btag.append(sort_jet[j][4])
			temp_jetid.append(sort_jet[j][5])


		# Muon veto
                deltaR2, deltaR = [],[]
                for j in range(t.nJet):
                        for m in range(t.nMuon):
                                if (m is not l1): continue
                                deltaR2 = (temp_eta[j]-t.Muon_eta.__getitem__(m))**2 + (temp_phi[j]-t.Muon_phi.__getitem__(m))**2
                                deltaR  = math.sqrt(deltaR2)
                                if (deltaR < 0.3):
                                        temp_mass.pop(j)
                                        temp_mass.insert(j,0)
                                        temp_pt.pop(j)
                                        temp_pt.insert(j,0)
                                        temp_eta.pop(j)
                                        temp_eta.insert(j,0)
                                        temp_phi.pop(j)
                                        temp_phi.insert(j,0)
                                        temp_btag.pop(j)
                                        temp_btag.insert(j,0)
                                        temp_jetid.pop(j)
                                        temp_jetid.insert(j,0)


#                       # Nearest Higgs mass selection algorithm
#                       mass_cand = []
#                       btag_good = False
#                       for k in range(t.nJet):
#                               if (temp_pt[k] > 20. and abs(temp_eta[k]) < 2.4 and temp_jetid[k] == 7):
#                                       if (temp_btag > 0.3093): btag_good = True
#                                       for l in range(k+1,t.nJet):
#                                               if (temp_pt[l] > 20. and abs(temp_eta[l]) < 2.4 and temp_jetid[l] == 7):
#                                                       if (temp_btag > 0.3093): btag_good = True
#
#                                                       jet1_p4 = TLorentzVector()
#                                                       jet1_p4.SetPtEtaPhiM(temp_pt[k],temp_eta[k],temp_phi[k],temp_mass[k])
#                                                       jet2_p4 = TLorentzVector()
#                                                       jet2_p4.SetPtEtaPhiM(temp_pt[l],temp_eta[l],temp_phi[l],temp_mass[l])
#                                                       jet_p4 = jet1_p4 + jet2_p4
#                                                       mass_cand.append(abs(jet_p4.M()-higgs_mass))
#                                                       if (btag_good == False): mass_cand.pop(-1)
#                       if (mass_cand == []): continue
#                       higgs_select = min(mass_cand)


		j1,j2 = -1,-1
		for k in range(len(temp_mass)):
                        if (flag_jet == True): break
			if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and temp_btag[k] > 0.3093 and temp_jetid[k] >= 7):
                        	for l in range(k+1,len(temp_mass)):
					if (flag_jet == True): break
                                	if (temp_pt[l] > 25 and abs(temp_eta[l]) < 2.4 and temp_btag[l] > 0.3093 and temp_jetid[l] >= 7):

						ljet_mass = temp_mass[k]
						ljet_pt   = temp_pt[k]
						ljet_eta  = temp_eta[k]
						ljet_phi  = temp_phi[k]
						ljet_btag = temp_btag[k]

						sjet_mass = temp_mass[l]
						sjet_pt   = temp_pt[l]
						sjet_eta  = temp_eta[l]
						sjet_phi  = temp_phi[l]
						sjet_btag = temp_btag[l]

						jet1_p4 = TLorentzVector()
                                                jet2_p4 = TLorentzVector()
                                                jet1_p4.SetPtEtaPhiM(ljet_pt,ljet_eta,ljet_phi,ljet_mass)
                                                jet2_p4.SetPtEtaPhiM(sjet_pt,sjet_eta,sjet_phi,sjet_mass)
                                                dijet_p4 = jet1_p4 + jet2_p4

						j1 = k
						j2 = l

						if (dijet_p4.M() == 0): continue
                                                if (dijet_p4.M() > 0): flag_jet = True


	input("Press Enter to continue...")


def correctorL1Pileup(eta,pt,rho,area):

	for i in range(len(jecl1)):

		if (eta > jecl1[i][1]): continue
		corr = 1 - ( jecl1[i][2] + rho*jecl1[i][3]*( 1 + jecl1[i][4]*np.log10(pt).item() ) )*area/pt
		pt_corr = pt*corr

	return pt_corr


if __name__ == "__main__":

	jname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/jec/jetl1corr_2016.txt"
	j = open(jname,"read")

	lines = j.read().split('\n')
	lines.pop(-1)

	j.close()

	jecl1 = []
	for line in lines:
		word = line.split(',')
		for i in range(len(word)):
			word[i] = float(word[i])
		jecl1.append(word)

	main()
