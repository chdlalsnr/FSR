#!/usr/bin/env python
import ROOT
from ROOT import *
from multiprocessing import Process, Queue
import numpy as np
import math
import glob
import sys
from array import array

def main():

	# Signal
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/64365B41-751A-7B4B-836B-DB4258390231.root')

#	# TTBar 
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist.extend(flist_had)
#        flist.extend(flist_semi)
#
#        # DY
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/*/*')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist.extend(flist_m50)
#        flist.extend(flist_0j)
#        flist.extend(flist_1j)
#        flist.extend(flist_2j)
#
#        # WJet
#        flist  = glob.glob('/T2_KR_KISTI/store/mc/RunIISummer16NanoAODv7/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#	flist  = glob.glob('/T2_KR_KISTI/store/mc/RunIISummer16NanoAODv7/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')

	h_J1_mass = TH1D("J1_mass", "Mass Distribution of Leading Jet", 40, 0, 200)
        h_J1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_J1_pt   = TH1D("J1_pt", "Transverse Momentum Distribution of Leading Jet", 40, 0, 400)
        h_J1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_J1_eta  = TH1D("J1_eta", "Eta Distribution of Leading Jet", 40, -2.5, 2.5)
        h_J1_eta.GetXaxis().SetTitle("Eta")
        h_J1_phi  = TH1D("J1_phi", "Phi Distribution of Leading Jet", 40, -3.14, 3.14)
        h_J1_phi.GetXaxis().SetTitle("Phi")	
	h_J1_btag  = TH1D("J1_btag", "btagDeepFlavB Distribution of Leading Jet", 20, 0, 1)
        h_J1_btag.GetXaxis().SetTitle("btagDeepFlavB")

	h_J2_mass = TH1D("J2_mass", "Mass Distribution of Sub-leading Jet", 40, 0, 200)
        h_J2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_J2_pt   = TH1D("J2_pt", "Transverse Momentum Distribution of Sub-leading Jet", 40, 0, 400)
        h_J2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_J2_eta  = TH1D("J2_eta", "Eta Distribution of Sub-leading Jet", 40, -2.5, 2.5)
        h_J2_eta.GetXaxis().SetTitle("Eta")
        h_J2_phi  = TH1D("J2_phi", "Phi Distribution of Sub-leading Jet", 40, -3.14, 3.14)
        h_J2_phi.GetXaxis().SetTitle("Phi")
	h_J2_btag  = TH1D("J2_btag", "btagDeepFlavB Distribution of Sub-leading Jet", 20, 0, 1)
        h_J2_btag.GetXaxis().SetTitle("btagDeepFlavB")

	h_L1_mass = TH1D("L1_mass", "Mass Distribution of Leading Lepton", 40, 0, 200)
        h_L1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_L1_pt   = TH1D("L1_pt", "Transverse Momentum Distribution of Leading Lepton", 40, 0, 400)
        h_L1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_L1_eta  = TH1D("L1_eta", "Eta Distribution of Leading Lepton", 40, -2.5, 2.5)
        h_L1_eta.GetXaxis().SetTitle("Eta")
        h_L1_phi  = TH1D("L1_phi", "Phi Distribution of Leading Lepton", 40, -3.14, 3.14)
        h_L1_phi.GetXaxis().SetTitle("Phi")
	h_L1_charge = TH1D("L1_charge", "Charge Distribution of Leading Lepton", 10, -1, 1)
	h_L1_charge.GetXaxis().SetTitle("Charge")

        h_L2_mass = TH1D("L2_mass", "Mass Distribution of Sub-leading Lepton", 40, 0, 200)
        h_L2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_L2_pt   = TH1D("L2_pt", "Transverse Momentum Distribution of Sub-leading Lepton", 40, 0, 400)
        h_L2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_L2_eta  = TH1D("L2_eta", "Eta Distribution of Sub-leading Lepton", 40, -2.5, 2.5)
        h_L2_eta.GetXaxis().SetTitle("Eta")
        h_L2_phi  = TH1D("L2_phi", "Phi Distribution of Sub-leading Lepton", 40, -3.14, 3.14)
        h_L2_phi.GetXaxis().SetTitle("Phi")
	h_L2_charge = TH1D("L2_charge", "Charge Distribution of Sub-leading Lepton", 10, -1, 1)
        h_L2_charge.GetXaxis().SetTitle("Charge")

        h_H1_mass = TH1D("H1_mass", "Mass Distribution of H1 <- MuMu", 40, 0, 200)
        h_H1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_H1_pt   = TH1D("H1_pt", "Transverse Momentum Distribution of H1 <- MuMu", 40, 0, 400)
        h_H1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_H1_eta  = TH1D("H1_eta", "Eta Distribution of H1 <- MuMu", 40, -2.5, 2.5)
        h_H1_eta.GetXaxis().SetTitle("Eta")
        h_H1_phi  = TH1D("H1_phi", "Phi Distribution of H1 <- MuMu", 40, -3.14, 3.14)
        h_H1_phi.GetXaxis().SetTitle("Phi")

        h_H2_mass = TH1D("H2_mass", "Mass Distribution of H2 <- bb", 40, 0, 200)
        h_H2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_H2_pt   = TH1D("H2_pt", "Transverse Momentum Distribution of H2 <- bb", 40, 0, 400)
        h_H2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_H2_eta  = TH1D("H2_eta", "Eta Distribution of H2 <- bb", 40, -2.5, 2.5)
        h_H2_eta.GetXaxis().SetTitle("Eta")
        h_H2_phi  = TH1D("H2_phi", "Phi Distribution of H2 <- bb", 40, -3.14, 3.14)
        h_H2_phi.GetXaxis().SetTitle("Phi")

	h_HH_mass = TH1D("HH_mass", "Mass Distribution of HH", 60, 0, 600)
        h_HH_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
        h_HH_pt   = TH1D("HH_pt", "Transverse Momentum Distribution of HH", 60, 0, 600)
        h_HH_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
        h_HH_eta  = TH1D("HH_eta", "Eta Distribution of HH", 40, -2.5, 2.5)
        h_HH_eta.GetXaxis().SetTitle("Eta")
        h_HH_phi  = TH1D("HH_phi", "Phi Distribution of HH", 40, -3.14, 3.14)
        h_HH_phi.GetXaxis().SetTitle("Phi")

	h_dr_ll = TH1D("DeltaR_ll", "DeltaR Distribution between Selected Leptons", 14, 0, 7)
	h_dr_ll.GetXaxis().SetTitle("DeltaR")
	h_dr_jj = TH1D("DeltaR_jj", "DeltaR Distribution between Selected Jets", 14, 0, 7)
	h_dr_jj.GetXaxis().SetTitle("DeltaR")

	h_MET_et   = TH1D("MET_et", "Transverse Energy Distribution of MET", 40, 0, 400)
	h_MET_et.GetXaxis().SetTitle("Transverse energy [GeV]")
	h_MET_pt   = TH1D("MET_pt", "Transverse Momentum Distribution of MET", 40, 0, 400)
	h_MET_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
	h_MET_phi  = TH1D("MET_phi", "Phi Distribution of MET", 40, -3.14, 3.14)
	h_MET_phi.GetXaxis().SetTitle("Phi")
	h_MET_mass = TH1D("MET_mass", "Mass_t Distribution of MET", 20, 0, 100)
	h_MET_mass.GetXaxis().SetTitle("Invariant mass [GeV]")

#	h_Higgsness = TH1D("Higgsness", "Distribution of Kinematics:Higgsness",100,0,100)
#	h_Topness = TH1D("Topness", "Distribution of Kinematics:Topness",100,0,100)
#	h_Topness.GetXaxis().SetTitle("Topness")

	ii,hh,ee = 0,0,0
	c_dilep,c_dijet = [],[]
	for fname in flist:

		if (ee > 5000000): break
		f = TFile(fname,"read")
		t = f.Get("Events")

		event_num = 30000
		#event_num = t.GetEntriesFast()
		print("Processing %s..."%(fname))
		ii = ii+1
		print("Have processed %i events..."%(ee))

		# Only MuMu channel considered at this point
		# HLT pass
		hlt_pass = []
		for i in range(event_num):
			t.GetEntry(i)
			hlt_mumu_1 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
			hlt_mumu_2 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
			hlt_mumu_3 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
			hlt_mumu_4 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ

			flag_good  = t.Flag_goodVertices
			flag_halo  = t.Flag_globalSuperTightHalo2016Filter
			flag_hbhen = t.Flag_HBHENoiseFilter
			flag_hbiso = t.Flag_HBHENoiseIsoFilter
			flag_dead  = t.Flag_EcalDeadCellTriggerPrimitiveFilter
			flag_badpf = t.Flag_BadPFMuonFilter
			flag_ecbad = t.Flag_ecalBadCalibFilter
			flag_eebad = t.Flag_eeBadScFilter

			if (hlt_mumu_1 == 1 or hlt_mumu_2 == 1 or hlt_mumu_3 == 1 or hlt_mumu_4 == 1):
				if (flag_good == 1 or flag_halo == 1 or flag_hbhen == 1 or flag_hbiso == 1 or flag_dead == 1 or flag_badpf == 1 or flag_ecbad == 1 or flag_eebad == 1):
					hlt_pass.append(i)
		hh = hh+len(hlt_pass)

		# =================================================================
		#Remained cut : The lepton isolation, defined as the scalar
		#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
		#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
		# Medium muon discrimination
		H1_mass_cand,H1_pt_cand,H1_eta_cand,H1_phi_cand = [],[],[],[]
		H2_mass_cand,H2_pt_cand,H2_eta_cand,H2_phi_cand = [],[],[],[]
		HH_mass_cand,HH_pt_cand,HH_eta_cand,HH_phi_cand = [],[],[],[]
		jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag = [],[],[],[],[]
                jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag = [],[],[],[],[]
                lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge = [],[],[],[],[]
                lep2_mass,lep2_pt,lep2_eta,lep2_phi,lep2_charge = [],[],[],[],[]
		dr_ll,dr_jj = [],[]
		top_cand,w_cand,hl_cand = [],[],[]
		temp_met_et,temp_met_pt,temp_met_eta,temp_met_phi = 0,0,0,0
		higgs_mass = 125.18
		for i in range(event_num):

			ee = ee+1
			t.GetEntry(i)
			if i not in hlt_pass: continue

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0

			ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
			sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0

			temp_pt,temp_eta,temp_phi,temp_mass,temp_pterr = [],[],[],[],[]
			temp_dxy,temp_dxyerr,temp_dz,temp_dzerr = [],[],[],[]
			temp_charge,temp_medid,temp_iso = [],[],[]
			temp_charge_e, temp_btag = [],[]

			flag_lep,flag_jet = False,False
			medid_n,btag_n = 0,0
			mu_pos,mu_neg,e_pos,e_neg = 0,0,0,0
			for j in range(t.nMuon):
				temp_mass.append(t.Muon_mass.__getitem__(j))
				temp_pt.append(t.Muon_pt.__getitem__(j))
				temp_eta.append(t.Muon_eta.__getitem__(j))
				temp_phi.append(t.Muon_phi.__getitem__(j))
				temp_pterr.append(t.Muon_ptErr.__getitem__(j))
				#temp_dxy.append(t.Muon_dxy.__getitem__(j))
				#temp_dxyerr.append(t.Muon_dxyErr.__getitem__(j))
				#temp_dz.append(t.Muon_dz.__getitem__(j))
				#temp_dzerr.append(t.Muon_dzErr.__getitem__(j))
				temp_iso.append(t.Muon_pfRelIso03_all.__getitem__(j))
				temp_charge.append(t.Muon_charge.__getitem__(j))
				if (t.Muon_charge.__getitem__(j) == 1): mu_pos += 1
				elif (t.Muon_charge.__getitem__(j) == -1): mu_neg += 1
				temp_medid.append(t.Muon_mediumId.__getitem__(j))
				if (t.Muon_mediumId.__getitem__(j) == True): medid_n += 1

			for j in range(t.nElectron):
				temp_charge_e.append(t.Electron_charge.__getitem__(j))
				if (t.Electron_charge.__getitem__(j) == 1): e_pos += 1
				elif (t.Electron_charge.__getitem__(j) == -1): e_neg += 1

			for j in range(t.nJet):
				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
				if (t.Jet_btagDeepFlavB.__getitem__(j) > 0.5): btag_n += 1

			for k in range(t.nMuon):
				if (flag_lep == True): break
				if (temp_pt[k] > 20. and abs(temp_eta[k]) < 2.4 and temp_medid[k] == True and temp_iso[k] < 0.15):
					for l in range(k+1,t.nMuon):
						if (temp_pt[l] > 10. and abs(temp_eta[l]) < 2.4 and temp_medid[l] == True and temp_iso[k] < 0.15):
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

								flag_lep = True
								break
				break

			# =================================================================
			# Jet discrimination
			temp_btag,temp_pt,temp_eta,temp_phi,temp_mass = [],[],[],[],[]
			temp_jetid,temp_csv = [],[]
			deltaR2, deltaR = [],[]
			temp_charge,temp_charge_e = [],[]
			temp_pt_l,temp_eta_l,temp_phi_l,temp_mass_l,temp_medid = [],[],[],[],[]
			temp_iso = []

			btag_n = 0
			mu_pos,mu_neg,e_pos,e_neg = 0,0,0,0
			for j in range(t.nJet): #(t.Jet_mass.__len__())
				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
				if (t.Jet_btagDeepFlavB.__getitem__(j) > 0.5): btag_n += 1
				temp_mass.append(t.Jet_mass.__getitem__(j))
				temp_pt.append(t.Jet_pt.__getitem__(j))
				temp_eta.append(t.Jet_eta.__getitem__(j))
				temp_phi.append(t.Jet_phi.__getitem__(j))
				#temp_breg.append(t.Jet_bRegCorr.__getitem__(j))
				temp_jetid.append(t.Jet_jetId.__getitem__(j))
				temp_csv.append(t.Jet_btagCSVV2.__getitem__(j))

			for j in range(t.nMuon):
				temp_charge.append(t.Muon_charge.__getitem__(j))
				if (t.Muon_charge.__getitem__(j) == 1): mu_pos += 1
				elif (t.Muon_charge.__getitem__(j) == -1): mu_neg += 1

				temp_mass_l.append(t.Muon_mass.__getitem__(j))
				temp_pt_l.append(t.Muon_pt.__getitem__(j))
				temp_eta_l.append(t.Muon_eta.__getitem__(j))
				temp_phi_l.append(t.Muon_phi.__getitem__(j))
				temp_medid.append(t.Muon_mediumId.__getitem__(j))
				temp_iso.append(t.Muon_pfRelIso03_all.__getitem__(j))

			for j in range(t.nElectron):
				temp_charge_e.append(t.Electron_charge.__getitem__(j))
				if (t.Electron_charge.__getitem__(j) == 1): e_pos += 1
				elif (t.Electron_charge.__getitem__(j) == -1): e_neg += 1	

			temp_met_et   = t.MET_sumEt
			temp_met_pt   = t.MET_pt
			#temp_met_eta  = t.MET_eta
			temp_met_phi  = t.MET_phi
			temp_met_mass = (temp_met_et**2 - temp_met_pt**2)**0.5

			# Muon veto
			for j in range(t.nJet):
				for m in range(t.nMuon):
					deltaR2 = (temp_eta[j]-t.Muon_eta.__getitem__(m))**2 + (temp_phi[j]-t.Muon_phi.__getitem__(m))**2
					deltaR  = math.sqrt(deltaR2)
					if (deltaR < 0.3):
						temp_btag.pop(j)
						temp_btag.insert(j,0)
						temp_mass.pop(j)
						temp_mass.insert(j,0)
						temp_pt.pop(j)
						temp_pt.insert(j,0)
						temp_eta.pop(j)
						temp_eta.insert(j,0)
						temp_phi.pop(j)
						temp_phi.insert(j,0)

			# Nearest Higgs mass selection algorithm
			mass_cand = []
			for k in range(t.nJet):
				if (temp_pt[k] > 20. and abs(temp_eta[k]) < 2.4 and temp_jetid[k] == 7 and temp_btag > 0.3093):
					for l in range(k+1,t.nJet):
						if (temp_pt[l] > 20. and abs(temp_eta[l]) < 2.4 and temp_jetid[l] == 7 and temp_btag > 0.3093):

							jet1_p4 = TLorentzVector()
							jet1_p4.SetPtEtaPhiM(temp_pt[k],temp_eta[k],temp_phi[k],temp_mass[k])
							jet2_p4 = TLorentzVector()
							jet2_p4.SetPtEtaPhiM(temp_pt[l],temp_eta[l],temp_phi[l],temp_mass[l])
							jet_p4 = jet1_p4 + jet2_p4
							mass_cand.append(abs(jet_p4.M()-higgs_mass))
			if (mass_cand == []): continue
			higgs_select = min(mass_cand)

			for k in range(t.nJet):
				if (flag_jet == True): break
				if (temp_pt[k] > 20. and abs(temp_eta[k]) < 2.4 and temp_jetid[k] == 7 and temp_btag > 0.3093):
					for l in range(k+1,t.nJet):
						if (temp_pt[l] > 20. and abs(temp_eta[l]) < 2.4 and temp_jetid[l] == 7 and temp_btag > 0.3093):

							jet1_p4 = TLorentzVector()
							jet1_p4.SetPtEtaPhiM(temp_pt[k],temp_eta[k],temp_phi[k],temp_mass[k])
							jet2_p4 = TLorentzVector()
							jet2_p4.SetPtEtaPhiM(temp_pt[l],temp_eta[l],temp_phi[l],temp_mass[l])
							jet_p4 = jet1_p4 + jet2_p4

							if (abs(jet_p4.M()-higgs_mass)==higgs_select):

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

								flag_jet = True
						break
				break

			# =================================================================

			dilep,dijet = [],[]
			met_et,met_pt,met_phi,met_mass = [],[],[],[]
			if (flag_lep == True and flag_jet == True):

				jet1_mass.append(ljet_mass)
				jet1_pt.append(ljet_pt)
				jet1_eta.append(ljet_eta)
				jet1_phi.append(ljet_phi)
				jet1_btag.append(ljet_btag)
				jet2_mass.append(sjet_mass)
                                jet2_pt.append(sjet_pt)
                                jet2_eta.append(sjet_eta)
                                jet2_phi.append(sjet_phi)
				jet2_btag.append(sjet_btag)

				lep1_mass.append(llep_mass)
                                lep1_pt.append(llep_pt)
                                lep1_eta.append(llep_eta)
                                lep1_phi.append(llep_phi)
				lep1_charge.append(llep_charge)
                                lep2_mass.append(slep_mass)
                                lep2_pt.append(slep_pt)
                                lep2_eta.append(slep_eta)
                                lep2_phi.append(slep_phi)
				lep2_charge.append(slep_charge)

				dr_ll.append(((llep_eta-slep_eta)**2+(llep_phi-slep_phi)**2)**0.5)
				dr_jj.append(((ljet_eta-sjet_eta)**2+(ljet_phi-sjet_phi)**2)**0.5)

				dijet.append(ljet_mass)
				dijet.append(ljet_pt)
				dijet.append(ljet_eta)
				dijet.append(ljet_phi)
				dijet.append(ljet_btag)
				dijet.append(sjet_mass)
				dijet.append(sjet_pt)
				dijet.append(sjet_eta)
				dijet.append(sjet_phi)
				dijet.append(sjet_btag)

				dilep.append(llep_mass)
				dilep.append(llep_pt)
				dilep.append(llep_eta)
				dilep.append(llep_phi)
				dilep.append(llep_charge)
				dilep.append(slep_mass)
				dilep.append(slep_pt)
				dilep.append(slep_eta)
				dilep.append(slep_phi)
				dilep.append(slep_charge)

				met_et.append(temp_met_et)
				met_pt.append(temp_met_pt)
				met_phi.append(temp_met_phi)
				met_mass.append(temp_met_mass)

				c_dilep.append(dilep)
				c_dijet.append(dijet)

		"""
			# Muons passed jet event selection 
			llep_mass_l,llep_pt_l,llep_eta_l,llep_phi_l = [],[],[],[]
                        slep_mass_l,slep_pt_l,slep_eta_l,slep_phi_l = [],[],[],[]
			for k in range(t.nMuon):
				if (temp_pt_l[k] > 20. and abs(temp_eta_l[k]) < 2.4 and temp_medid[k] == True and temp_iso[k] < 0.15):
					for l in range(k+1,t.nMuon):
						if (temp_pt_l[l] > 10. and abs(temp_eta_l[l]) < 2.4 and temp_medid[l] == True and temp_iso[k] < 0.15):
							if (temp_charge[k] != temp_charge[l]):

								llep_mass_l.append(temp_mass_l[k])
								llep_pt_l.append(temp_pt_l[k])
								llep_eta_l.append(temp_eta_l[k])
								llep_phi_l.append(temp_phi_l[k])

								slep_mass_l.append(temp_mass_l[l])
								slep_pt_l.append(temp_pt_l[l])
								slep_eta_l.append(temp_eta_l[l])
								slep_phi_l.append(temp_phi_l[l])
								lepton_passed.append(i)
						break
				break
		"""

		# =================================================================
		# Object Reconstruction
		top_mass,w_mass,z_mass = 172.26,80.379,91.1876
		higgsness,topness = [],[]
		# Higgs reconstruction
		for i in range(len(c_dilep)):

			lep1_p4 = TLorentzVector()
			lep2_p4 = TLorentzVector()
			lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
			lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])
			#lep1_p4.SetPtEtaPhiM(llep_pt[i],llep_eta[i],llep_phi[i],llep_mass[i])
			#lep2_p4.SetPtEtaPhiM(slep_pt[i],slep_eta[i],slep_phi[i],slep_mass[i])
			lep_p4 = lep1_p4 + lep2_p4

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()
			jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][0])
			jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])
			#jet1_p4.SetPtEtaPhiM(ljet_pt[i],ljet_eta[i],ljet_phi[i],ljet_mass[i])
			#jet2_p4.SetPtEtaPhiM(sjet_pt[i],sjet_eta[i],sjet_phi[i],sjet_mass[i])
			jet_p4 = jet1_p4 + jet2_p4

			HH_p4 = TLorentzVector()
			HH_p4 = lep_p4 + jet_p4

			if (lep_p4.M() > 12. and lep_p4.M() < 76.):

				H1_mass_cand.append(lep_p4.M())
				H1_pt_cand.append(lep_p4.Pt())
				H1_eta_cand.append(lep_p4.Eta())
				H1_phi_cand.append(lep_p4.Phi())

				H2_mass_cand.append(jet_p4.M())
				H2_pt_cand.append(jet_p4.Pt())
				H2_eta_cand.append(jet_p4.Eta())
				H2_phi_cand.append(jet_p4.Phi())

				HH_mass_cand.append(HH_p4.M())
				HH_pt_cand.append(HH_p4.Pt())
				HH_eta_cand.append(HH_p4.Eta())
				HH_phi_cand.append(HH_p4.Phi())

		f.Close()

#			# =================================================================
#			# Define Higgsness and Topness
#
#			met_p4 = TLorentzVector()
#                        met_p4.SetPtEtaPhiM(c_met[i][1],c_met[i][2],c_met[i][3],c_met[i][0]) # no eta for MET
#
#			flag_ht = False
#			top1_p4 = TLorentzVector()
#			top2_p4 = TLorentzVector()
#			wp_p4  = TLorentzVector()
#			wm_p4  = TLorentzVector()
#			# hl_p4 = lep_p4
#
#			if (c_dilep[i][4] > 0):
#				wp_p4 = met_p4 + lep1_p4
#				wm_p4 = met_p4 + lep2_p4
#			else:
#				wp_p4 = met_p4 + lep2_p4
#				wm_p4 = met_p4 + lep1_p4
#
#			top1_p4 = wp_p4 + jet1_p4
#			top2_p4 = wp_p4 + jet2_p4
#			mass_top,mass_w = [],[]
#			mass_top.append(top1_p4.M())
#			mass_top.append(top2_p4.M())
#			mass_w.append(wp_p4.M())
#			mass_w.append(wm_p4.M())
#
#			dev_t = ((mass_top[0]**2+mass_top[1]**2)/2) - ((mass_top[0]+mass_top[1])/2)**2
#			dev_w = ((mass_w[0]**2+mass_w[1]**2)/2) - ((mass_w[0]+mass_w[1])/2)**2
#
#			j1lpmet_p4 = TLorentzVector()
#			j1lmmet_p4 = TLorentzVector()
#			j2lpmet_p4 = TLorentzVector()
#			j2lmmet_p4 = TLorentzVector()
#			lpmet_p4   = TLorentzVector()
#			lmmet_p4   = TLorentzVector()
#
#			if (c_dilep[i][4] > 0):
#				j1lpmet_p4 = jet1_p4 + lep1_p4 + met_p4
#				j1lmmet_p4 = jet1_p4 + lep2_p4 + met_p4
#				j2lpmet_p4 = jet2_p4 + lep1_p4 + met_p4
#				j2lmmet_p4 = jet2_p4 + lep2_p4 + met_p4
#				lpmet_p4   = lep1_p4 + met_p4
#				lmmet_p4   = lep2_p4 + met_p4
#			else:
#				j1lpmet_p4 = jet1_p4 + lep2_p4 + met_p4
#				j1lmmet_p4 = jet1_p4 + lep1_p4 + met_p4
#				j2lpmet_p4 = jet2_p4 + lep2_p4 + met_p4
#				j2lmmet_p4 = jet2_p4 + lep1_p4 + met_p4
#				lpmet_p4   = lep2_p4 + met_p4
#				lmmet_p4   = lep1_p4 + met_p4
#
#			xi_12 = ((j1lpmet_p4.M()**2-top_mass**2)**2/dev_t**4)+((lpmet_p4.M()**2-w_mass**2)**2/dev_w**4)+((j2lmmet_p4.M()**2-top_mass**2)**2/dev_t**4)+((lmmet_p4.M()**2-w_mass**2)**2/dev_w**4) 
#			xi_21 = ((j2lpmet_p4.M()**2-top_mass**2)**2/dev_t**4)+((lpmet_p4.M()**2-w_mass**2)**2/dev_w**4)+((j1lmmet_p4.M()**2-top_mass**2)**2/dev_t**4)+((lmmet_p4.M()**2-w_mass**2)**2/dev_w**4) 
#
#			mass_offw = []
#			if (wp_p4.M() >= 0 and wp_p4.M() < higgs_mass-w_mass):
#				mass_offw.append(wp_p4.M())
#			if (wm_p4.M() >= 0 and wm_p4.M() < higgs_mass-w_mass):
#				mass_offw.append(wm_p4.M())
#
#			temp_offw = []
#			for i in range(len(mass_offw)):
#				temp_offw.append(abs(mass_offw[i]-w_mass))
#			crit_offw = min(temp_offw)
#
#			for i in range(len(mass_offw)):
#				if (abs(mass_offw[i]-w_mass) == crit_offw):
#					offw_mass = mass_off[i]
#					flag_ht = True
#
#			dilepmet_p4 = TLorentzVector()
#			dilepmet_p4 = lep_p4 + met_p4
#
#			if (flag_ht == False): continue
#
#			if (xi_12 > xi_21):
#                                topness.append(xi_21)
#                        else: topness.append(xi_12)
#
	# =================================================================

	dr = "signal"
        ff = TFile("./%s.root"%(dr),"RECREATE")
        tt = TTree(dr,dr)

#	J1_mass,J1_pt,J1_eta,J1_phi,J1_btag = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
#	J2_mass,J2_pt,J2_eta,J2_phi,J2_btag = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
#	L1_mass,L1_pt,L1_eta,L1_phi,L1_charge = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
#        L2_mass,L2_pt,L2_eta,L2_phi,L2_charge = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
#	MET_et,MET_pt,MET_phi,MET_mass = np.array([]),np.array([]),np.array([]),np.array([])
#	DeltaR_jj,DeltaR_ll = np.array([]),np.array([])
#	H1_mass,H1_pt,H1_eta,H1_phi = np.array([]),np.array([]),np.array([]),np.array([])
#	H2_mass,H2_pt,H2_eta,H2_phi = np.array([]),np.array([]),np.array([]),np.array([])
#	HH_mass,HH_pt,HH_eta,HH_phi = np.array([]),np.array([]),np.array([]),np.array([])

	J1_mass,J1_pt,J1_eta,J1_phi,J1_btag = 0,0,0,0,0
	J2_mass,J2_pt,J2_eta,J2_phi,J2_btag = 0,0,0,0,0
	L1_mass,L1_pt,L1_eta,L1_phi,L1_charge = 0,0,0,0,0
        L2_mass,L2_pt,L2_eta,L2_phi,L2_charge = 0,0,0,0,0
	MET_et,MET_pt,MET_phi,MET_mass = 0,0,0,0
	DeltaR_jj,DeltaR_ll = 0,0
	H1_mass,H1_pt,H1_eta,H1_phi = 0,0,0,0
	H2_mass,H2_pt,H2_eta,H2_phi = 0,0,0,0
	HH_mass,HH_pt,HH_eta,HH_phi = 0,0,0,0

	tt.Branch("J1_mass",J1_mass,"J1_mass/F")
	tt.Branch("J1_pt",J1_pt,"J1_pt/F")
	tt.Branch("J1_eta",J1_eta,"J1_eta/F")
	tt.Branch("J1_phi",J1_phi,"J1_phi/F")
	tt.Branch("J1_btag",J1_btag,"J1_btag/F")

	tt.Branch("J2_mass",J2_mass,"J2_mass/F")
        tt.Branch("J2_pt",J2_pt,"J2_pt/F")
        tt.Branch("J2_eta",J2_eta,"J2_eta/F")
        tt.Branch("J2_phi",J2_phi,"J2_phi/F")
        tt.Branch("J2_btag",J2_btag,"J2_btag/F")

	tt.Branch("L1_mass",L1_mass,"L1_mass/F")
        tt.Branch("L1_pt",L1_pt,"L1_pt/F")
        tt.Branch("L1_eta",L1_eta,"L1_eta/F")
        tt.Branch("L1_phi",L1_phi,"L1_phi/F")
        tt.Branch("L1_charge",L1_charge,"L1_charge/I")

        tt.Branch("L2_mass",L2_mass,"L2_mass/F")
        tt.Branch("L2_pt",L2_pt,"L2_pt/F")
        tt.Branch("L2_eta",L2_eta,"L2_eta/F")
        tt.Branch("L2_phi",L2_phi,"L2_phi/F")
        tt.Branch("L2_charge",L2_charge,"L2_charge/I")

	tt.Branch("MET_et",MET_et,"MET_et/F")
	tt.Branch("MET_pt",MET_pt,"MET_pt/F")
	tt.Branch("MET_phi",MET_phi,"MET_phi/F")
	tt.Branch("MET_mass",MET_mass,"MET_mass/F")

	tt.Branch("DeltaR_ll",DeltaR_ll,"DeltaR_ll/F")
	tt.Branch("DeltaR_jj",DeltaR_jj,"DeltaR_jj/F")

	tt.Branch("H1_mass",H1_mass,"H1_mass/F")
	tt.Branch("H1_pt",H1_pt,"H1_pt/F")
	tt.Branch("H1_eta",H1_eta,"H1_eta/F")
	tt.Branch("H1_phi",H1_phi,"H1_phi/F")

	tt.Branch("H2_mass",H2_mass,"H2_mass/F")
        tt.Branch("H2_pt",H2_pt,"H2_pt/F")
        tt.Branch("H2_eta",H2_eta,"H2_eta/F")
        tt.Branch("H2_phi",H2_phi,"H2_phi/F")

	tt.Branch("HH_mass",HH_mass,"HH_mass/F")
        tt.Branch("HH_pt",HH_pt,"HH_pt/F")
        tt.Branch("HH_eta",HH_eta,"HH_eta/F")
        tt.Branch("HH_phi",HH_phi,"HH_phi/F")

#	f_J1_mass = TLeaf(br,"h_J1_mass","TH1D")
#	f_J1_pt   = TLeaf(br,"h_J1_pt","TH1D")
#	f_J1_eta  = TLeaf(br,"h_J1_eta","TH1D")
#	f_J1_phi  = TLeaf(br,"h_J1_phi","TH1D")
#	f_J1_btag = TLeaf(br,"h_J1_btag","TH1D")
#
#	f_J2_mass = TLeaf(br,"h_J2_mass","TH1D")
#        f_J2_pt   = TLeaf(br,"h_J2_pt","TH1D")
#        f_J2_eta  = TLeaf(br,"h_J2_eta","TH1D")
#        f_J2_phi  = TLeaf(br,"h_J2_phi","TH1D")
#        f_J2_btag = TLeaf(br,"h_J2_btag","TH1D")
#
#	f_L1_mass   = TLeaf(br,"h_L1_mass","TH1D")
#        f_L1_pt     = TLeaf(br,"h_L1_pt","TH1D")
#        f_L1_eta    = TLeaf(br,"h_L1_eta","TH1D")
#        f_L1_phi    = TLeaf(br,"h_L1_phi","TH1D")
#        f_L1_charge = TLeaf(br,"h_L1_charge","TH1D")
#
#	f_L2_mass   = TLeaf(br,"h_L2_mass","TH1D")
#        f_L2_pt     = TLeaf(br,"h_L2_pt","TH1D")
#        f_L2_eta    = TLeaf(br,"h_L2_eta","TH1D")
#        f_L2_phi    = TLeaf(br,"h_L2_phi","TH1D")
#        f_L2_charge = TLeaf(br,"h_L2_charge","TH1D")
#
#	f_MET_et   = TLeaf(br,"h_MET_et","TH1D")
#	f_MET_pt   = TLeaf(br,"h_MET_pt","TH1D")
#	f_MET_phi  = TLeaf(br,"h_MET_phi","TH1D")
#	f_MET_mass = TLeaf(br,"h_MET_mass","TH1D")
#
#	f_DeltaR_jj = TLeaf(br,"h_dr_jj","TH1D")
#	f_DeltaR_ll = TLeaf(br,"h_dr_ll","TH1D")
#
#	f_H1_mass = TLeaf(br,"h_H1_mass","TH1D")
#	f_H1_pt   = TLeaf(br,"h_H1_pt","TH1D")
#	f_H1_eta  = TLeaf(br,"h_H1_eta","TH1D")
#	f_H1_phi  = TLeaf(br,"h_H1_phi","TH1D")
#
#	f_H2_mass = TLeaf(br,"h_H2_mass","TH1D")
#        f_H2_pt   = TLeaf(br,"h_H2_pt","TH1D")
#        f_H2_eta  = TLeaf(br,"h_H2_eta","TH1D")
#        f_H2_phi  = TLeaf(br,"h_H2_phi","TH1D")
#
#	f_HH_mass = TLeaf(br,"h_HH_mass","TH1D")
#        f_HH_pt   = TLeaf(br,"h_HH_pt","TH1D")
#        f_HH_eta  = TLeaf(br,"h_HH_eta","TH1D")
#        f_HH_phi  = TLeaf(br,"h_HH_phi","TH1D")

	# Plotting
	for i in range(len(jet1_mass)):
		if (jet1_mass[i] != 0 and jet1_pt[i] != 0):
			J1_mass = jet1_mass[i]
			J1_pt   = jet1_pt[i]
			J1_eta  = jet1_eta[i]
			J1_phi  = jet1_phi[i]
			J1_btag = jet1_btag[i]
			tt.Fill()

	for i in range(len(jet2_mass)):
        	if (jet2_mass[i] != 0 and jet2_pt[i] != 0):
			J2_mass = jet2_mass[i]
			J2_pt   = jet2_pt[i]
			J2_eta  = jet2_eta[i]
			J2_phi  = jet2_phi[i]
			J2_btag = jet2_btag[i]
			tt.Fill()

	for i in range(len(lep1_mass)):
        	if (lep1_mass[i] != 0 and lep1_pt[i] != 0):
			L1_mass = lep1_mass[i]
			L1_pt = lep1_pt[i]
			L1_eta = lep1_eta[i]
			L1_phi = lep1_phi[i]
			L1_charge = lep1_charge[i]
			tt.Fill()

	for i in range(len(lep2_mass)):
		if (lep2_mass[i] != 0 and lep2_pt[i] != 0):
			L2_mass = lep2_mass[i]
			L2_pt = lep2_pt[i]
			L2_eta = lep2_eta[i]
			L2_phi = lep2_phi[i]
			L2_charge = lep2_charge[i]
			tt.Fill()

	for i in range(len(met_et)):
		MET_et = met_et[i]
		MET_pt = met_pt[i]
		MET_phi = met_phi[i]
		MET_mass = met_mass[i]
		tt.Fill()

	for i in range(len(dr_ll)):
		Delta_ll = dr_ll[i]
		Delta_jj = dr_jj[i]
		tt.Fill()

#	H1_mass = np.zeros(len(H1_mass_cand))
#        H1_pt   = np.zeros(len(H1_pt_cand))
#        H1_eta  = np.zeros(len(H1_eta_cand))
#        H1_phi  = np.zeros(len(H1_phi_cand))
	for i in range(len(H1_mass_cand)):
		# considered both bbWW & bbZZ
		if (H1_mass_cand[i] != 0 and H1_pt_cand[i] != 0):
			H1_mass = H1_mass_cand[i]
			H1_pt = H1_pt_cand[i]
			H1_eta = H1_eta_cand[i]
			H1_phi = H1_phi_cand[i]
			tt.Fill()

#	H2_mass = np.zeros(len(H2_mass_cand))
#	H2_pt   = np.zeros(len(H2_pt_cand))
#	H2_eta  = np.zeros(len(H2_eta_cand))
#	H2_phi  = np.zeros(len(H2_phi_cand))
	for i in range(len(H2_mass_cand)):
		if (H2_mass_cand[i] != 0 and H2_pt_cand[i] != 0):
			H2_mass = H2_mass_cand[i]
			H2_pt = H2_pt_cand[i]
			H2_eta = H2_eta_cand[i]
			H2_phi = H2_phi_cand[i]
			tt.Fill()

#	HH_mass = np.zeros(len(HH_mass_cand))
#        HH_pt   = np.zeros(len(HH_pt_cand))
#        HH_eta  = np.zeros(len(HH_eta_cand))
#        HH_phi  = np.zeros(len(HH_phi_cand))
        for i in range(len(HH_mass_cand)):
                if (HH_mass_cand[i] != 0 and HH_pt_cand[i] != 0):
			HH_mass = HH_mass_cand[i]
			HH_pt = HH_pt_cand[i]
			HH_eta = HH_eta_cand[i]
			HH_phi = HH_phi_cand[i]
			tt.Fill()

#	tree.Write("",ROOT.TObject.kOverwrite);
	ff.cd()
	tt.Write()

#	c1 = TCanvas()
#	c1.Divide(2,2)
#	c2 = TCanvas()
#	c2.Divide(2,2)
#
#	c1.cd(1)
#	h_H1_mass.Draw()
#	c1.cd(2)
#	h_H1_pt.Draw()
#	c1.cd(3)
#	h_H1_eta.Draw()
#	c1.cd(4)
#	h_H1_phi.Draw()
#	c1.SaveAs("Higgs1_%s.png"%(dr))
#
#	c2.cd(1)
#	h_H2_mass.Draw()
#	c2.cd(2)
#	h_H2_pt.Draw()
#	c2.cd(3)
#	h_H2_eta.Draw()
#	c2.cd(4)
#	h_H2_phi.Draw()
#	c2.SaveAs("Higgs2_%s.png"%(dr))

        ff.Close()

	print("Have processed %i files..."%(ii))
	print("Have processed total %i events..."%(ee))
	print("Event number passed hlt is %i..."%(hh))
	input("Press Enter to continue...")

if __name__ == "__main__":
	main()

