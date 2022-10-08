from ROOT import *
from multiprocessing import Process, Queue
import numpy as np
import math
import glob
import sys
import pandas as pd

def main():

	# Signal
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/64365B41-751A-7B4B-836B-DB4258390231.root')

#	# TTBar 
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist.extend(flist_had)
#        flist.extend(flist_semi)
#
#        # DY
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/*/*')
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

#	h_Higgsness = TH1D("Higgsness", "Distribution of Kinematics:Higgsness",100,0,100)
#	h_Topness = TH1D("Topness", "Distribution of Kinematics:Topness",100,0,100)
#	h_Topness.GetXaxis().SetTitle("Topness")

	ii,hh,ee = 0,0,0
	data = []
	for fname in flist:

		if (ee > 5000000): break
		f = TFile(fname,"read")
		t = f.Get("Events")

		#event_num = 300000
		event_num = t.GetEntriesFast()
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
		met_et,met_pt,met_eta,met_phi = 0,0,0,0
		c_dilep,c_dijet,c_met = [],[],[]
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

			met_et   = t.MET_sumEt
			met_pt   = t.MET_pt
			#met_eta  = t.MET_eta
			met_phi  = t.MET_phi
			met_mass = (met_et**2 - met_pt**2)**0.5

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

			dilep,dijet,met = [],[],[]
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

				met.append(met_et)
				met.append(met_pt)
				met.append(met_phi)
				met.append(met_mass)

				c_dilep.append(dilep)
				c_dijet.append(dijet)
				c_met.append(met)

		# =================================================================
		# Object Reconstruction
		top_mass,w_mass,z_mass = 172.26,80.379,91.1876
		# Higgs reconstruction
		for i in range(len(c_dilep)):

			lep1_p4 = TLorentzVector()
			lep2_p4 = TLorentzVector()
			lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
			lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])
			lep_p4 = lep1_p4 + lep2_p4

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()
			jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][0])
			jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])
			jet_p4 = jet1_p4 + jet2_p4

			HH_p4 = TLorentzVector()
			HH_p4 = lep_p4 + jet_p4

			#if (lep_p4.M() > 12. and lep_p4.M() < 76.):

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

		# =================================================================
		# Save CSV file
		for i in range(len(c_dilep)):
			if (H1_mass_cand[i] == 0 or H1_pt_cand[i] == 0): continue
			if (H2_mass_cand[i] == 0 or H2_pt_cand[i] == 0): continue
			data.append([jet1_mass[i],jet1_pt[i],jet1_eta[i],jet1_phi[i],jet1_btag[i],
				jet2_mass[i],jet2_pt[i],jet2_eta[i],jet2_phi[i],jet2_btag[i],
				lep1_mass[i],lep1_pt[i],lep1_eta[i],lep1_phi[i],lep1_charge[i],
				lep2_mass[i],lep2_pt[i],lep2_eta[i],lep2_phi[i],lep2_charge[i],
				dr_jj[i],dr_ll[i],
				H1_mass_cand[i],H1_pt_cand[i],H1_eta_cand[i],H1_phi_cand[i],
				H2_mass_cand[i],H2_pt_cand[i],H2_eta_cand[i],H2_phi_cand[i],
				HH_mass_cand[i],HH_pt_cand[i],HH_eta_cand[i],HH_phi_cand[i]])

	df = pd.DataFrame(data, columns=
				['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet1_btag',
				 'jet2_mass','jet2_pt','jet2_eta','jet2_phi','jet2_btag',
				 'lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep1_charge',
				 'lep2_mass','lep2_pt','lep2_eta','lep2_phi','lep2_charge',
				 'dr_jj','dr_ll',
				 'H1_mass','H1_pt','H1_eta','H1_phi',
				 'H2_mass','H2_pt','H2_eta','H2_phi',
				 'HH_mass','HH_pt','HH_eta','HH_phi'])


	df.to_csv("./output/dy_m50.csv", header=True, index=False)

	print("Have processed %i files..."%(ii))
	print("Have processed total %i events..."%(ee))
	print("Event number passed hlt is %i..."%(hh))
	input("Press Enter to continue...")


if __name__ == "__main__":
	main()

