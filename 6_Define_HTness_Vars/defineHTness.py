from ROOT import *
import numpy as np
import math
import glob
import sys


def defineHTness():

	# Signal
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/*')

	# TTBar 
        #flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')

	ii,hh,ee = 0,0,0
        higgsness,topness = [],[]
	for fname in flist:

		if (ee > 3000000): break
		f = TFile(fname,"read")
		t = f.Get("Events")

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
		jet3_mass,jet3_pt,jet3_eta,jet3_phi = [],[],[],[]
		lep3_mass,lep3_pt,lep3_eta,lep3_phi = [],[],[],[]
		dr_ll,dr_jj,genWeight = [],[],[]
                c_dilep,c_dijet,c_met = [],[],[]
                met_et,met_pt,met_phi = 0,0,0
                h_mass = 125.18
		for i in range(event_num):

			ee = ee+1
			t.GetEntry(i)
			if i not in hlt_pass: continue
			if not (t.nJet >= 3 or (t.nJet >= 1 and t.nFatJet >= 1)): continue

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0
			ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
			sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0
			sslep_mass,sslep_pt,sslep_eta,sslep_phi = 0,0,0,0
			ssjet_mass,ssjet_pt,ssjet_eta,ssjet_phi = 0,0,0,0

			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
			temp_dxy,temp_dz,temp_iso = [],[],[]
			temp_charge,temp_medid = [],[]
			temp_charge_e, temp_btag = [],[]

			flag_lep,flag_jet = False,False
                        mu_pos,mu_neg,e_pos,e_neg = 0,0,0,0
                        genweight = t.genWeight
			for j in range(t.nMuon):
				temp_mass.append(t.Muon_mass.__getitem__(j))
				temp_pt.append(t.Muon_pt.__getitem__(j))
				temp_eta.append(t.Muon_eta.__getitem__(j))
				temp_phi.append(t.Muon_phi.__getitem__(j))
				temp_dxy.append(t.Muon_dxy.__getitem__(j))
				temp_dz.append(t.Muon_dz.__getitem__(j))
				temp_iso.append(t.Muon_miniPFRelIso_all.__getitem__(j))
				temp_charge.append(t.Muon_charge.__getitem__(j))
				if (t.Muon_charge.__getitem__(j) == 1): mu_pos += 1
				elif (t.Muon_charge.__getitem__(j) == -1): mu_neg += 1
				temp_medid.append(t.Muon_mediumId.__getitem__(j))

#			for j in range(t.nElectron):
#				temp_charge_e.append(t.Electron_charge.__getitem__(j))
#				if (t.Electron_charge.__getitem__(j) == 1): e_pos += 1
#				elif (t.Electron_charge.__getitem__(j) == -1): e_neg += 1
#
#			for j in range(t.nJet):
#				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
#                                if (t.Jet_btagDeepFlavB.__getitem__(j) > 0.3093): # medium WP
#
#			for j in range(t.nFatJet):
#                                temp_fatpt.append(t.FatJet_pt.__getitem__(j))

			for k in range(t.nMuon):
				if (flag_lep == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) <= 2.4 and abs(temp_dxy[k]) <= 0.05 and abs(temp_dz[k]) <= 0.1 and temp_medid[k] == True):
					for l in range(k+1,t.nMuon):
						if (temp_pt[l] > 15 and abs(temp_eta[l]) < 2.4 and abs(temp_dxy[l]) <= 0.05 and abs(temp_dz[l]) <= 0.1 and temp_medid[l] == True):
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

                                                                lep1_p4 = TLorentzVector()
                                                                lep2_p4 = TLorentzVector()
                                                                lep1_p4.SetPtEtaPhiM(llep_pt,llep_eta,llep_phi,llep_mass)
                                                                lep2_p4.SetPtEtaPhiM(slep_pt,slep_eta,slep_phi,slep_mass)
                                                                dilep_p4 = lep1_p4 + lep2_p4

								if not (dilep_p4.M() > 12): continue
                                                                if (dilep_p4.M() > 12): flag_lep = True

#									sslep_mass = temp_mass[l+1]
#									sslep_pt = temp_pt[l+1]
#									sslep_eta = temp_eta[l+1]
#									sslep_phi = temp_phi[l+1]

                                                                break
                                break


			# =================================================================
			# Jet discrimination
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
			temp_btag,temp_jetid,temp_puid,temp_csv = [],[],[],[]
			temp_charge,temp_charge_e = [],[]

			for j in range(t.nJet): #(t.Jet_mass.__len__())
				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
				temp_mass.append(t.Jet_mass.__getitem__(j))
				temp_pt.append(t.Jet_pt.__getitem__(j))
				temp_eta.append(t.Jet_eta.__getitem__(j))
				temp_phi.append(t.Jet_phi.__getitem__(j))
				temp_jetid.append(t.Jet_jetId.__getitem__(j))
#				temp_csv.append(t.Jet_btagCSVV2.__getitem__(j))

#			temp_pt_l,temp_eta_l,temp_phi_l,temp_mass_l = [],[],[],[],[]
#			for j in range(t.nMuon):
#				temp_charge.append(t.Muon_charge.__getitem__(j))
#				if (t.Muon_charge.__getitem__(j) == 1): mu_pos += 1
#				elif (t.Muon_charge.__getitem__(j) == -1): mu_neg += 1
#
#				temp_mass_l.append(t.Muon_mass.__getitem__(j))
#				temp_pt_l.append(t.Muon_pt.__getitem__(j))
#				temp_eta_l.append(t.Muon_eta.__getitem__(j))
#				temp_phi_l.append(t.Muon_phi.__getitem__(j))
#
#			for j in range(t.nElectron):
#				temp_charge_e.append(t.Electron_charge.__getitem__(j))
#				if (t.Electron_charge.__getitem__(j) == 1): e_pos += 1
#				elif (t.Electron_charge.__getitem__(j) == -1): e_neg += 1	

			met_et   = t.MET_sumEt
			met_pt   = t.MET_pt
			met_phi  = t.MET_phi
			met_mass = (met_et**2 - met_pt**2)**0.5

			# Muon veto
			deltaR2, deltaR = [],[]
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

#			# Nearest Higgs mass selection algorithm
#			mass_cand = []
#			for k in range(t.nJet):
#				if (temp_pt[k] > 20. and abs(temp_eta[k]) < 2.4 and temp_jetid[k] == 7 and temp_btag > 0.3093):
#					for l in range(k+1,t.nJet):
#						if (temp_pt[l] > 20. and abs(temp_eta[l]) < 2.4 and temp_jetid[l] == 7 and temp_btag > 0.3093):
#
#							jet1_p4 = TLorentzVector()
#							jet1_p4.SetPtEtaPhiM(temp_pt[k],temp_eta[k],temp_phi[k],temp_mass[k])
#							jet2_p4 = TLorentzVector()
#							jet2_p4.SetPtEtaPhiM(temp_pt[l],temp_eta[l],temp_phi[l],temp_mass[l])
#							jet_p4 = jet1_p4 + jet2_p4
#							mass_cand.append(abs(jet_p4.M()-higgs_mass))
#			if (mass_cand == []): continue
#			higgs_select = min(mass_cand)

			for k in range(t.nJet):
				if (flag_jet == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4  and temp_jetid[k] == 7 and temp_btag[k] > 0.3093): # jetid == 1: Loose, jetid == 2: Medium
					for l in range(k+1,t.nJet):
						if (temp_pt[l] > 25 and abs(temp_eta[l]) < 2.4  and temp_jetid[l] == 7 and temp_btag[l] > 0.3093):

#							jet1_p4 = TLorentzVector()
#							jet1_p4.SetPtEtaPhiM(temp_pt[k],temp_eta[k],temp_phi[k],temp_mass[k])
#							jet2_p4 = TLorentzVector()
#							jet2_p4.SetPtEtaPhiM(temp_pt[l],temp_eta[l],temp_phi[l],temp_mass[l])
#							jet_p4 = jet1_p4 + jet2_p4
#
#							if (abs(jet_p4.M()-higgs_mass)==higgs_select):

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

							if (dijet_p4.M() == 0): continue
                                                        if (dijet_p4.M() > 0): flag_jet = True
						break
				break

			# =================================================================
			dilep,dijet,met = [],[],[]
			if (flag_lep == True and flag_jet == True):
				hh = hh+1

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

#				jet3_mass.append(ssjet_mass)
#				jet3_pt.append(ssjet_pt)
#				jet3_eta.append(ssjet_eta)
#				jet3_phi.append(ssjet_phi)
#				lep3_mass.append(sslep_mass)
#				lep3_pt.append(sslep_pt)
#				lep3_eta.append(sslep_eta)
#				lep3_phi.append(sslep_phi)

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
#				dijet.append(ssjet_mass)
#				dijet.append(ssjet_pt)
#				dijet.append(ssjet_eta)
#				dijet.append(ssjet_phi)

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
#				dilep.append(sslep_mass)
#				dilep.append(sslep_pt)
#				dilep.append(sslep_eta)
#				dilep.append(sslep_phi)

				met.append(met_mass)
				met.append(met_pt)
				met.append(met_phi)

				c_dilep.append(dilep)
				c_dijet.append(dijet)
				c_met.append(met)

		# =================================================================
		# Object Reconstruction
                t_mass,w_mass,z_mass = 172.26,80.379,91.1876
		# Higgs reconstruction
		for i in range(len(c_dilep)):

			lep1_p4 = TLorentzVector()
			lep2_p4 = TLorentzVector()
			lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
			lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])
			met_p4 = TLorentzVector()
                        met_p4.SetPtEtaPhiM(c_met[i][0],0,c_met[i][1],c_met[i][2])
                        h1_p4 = lep1_p4 + lep2_p4 + met_p4

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()
			jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][0])
			jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])
			h2_p4 = jet1_p4 + jet2_p4

			HH_p4 = TLorentzVector()
			HH_p4 = h1_p4 + h2_p4

			H1_mass_cand.append(h1_p4.M())
			H1_pt_cand.append(h1_p4.Pt())
			H1_eta_cand.append(h1_p4.Eta())
			H1_phi_cand.append(h1_p4.Phi())

			H2_mass_cand.append(h2_p4.M())
			H2_pt_cand.append(h2_p4.Pt())
			H2_eta_cand.append(h2_p4.Eta())
			H2_phi_cand.append(h2_p4.Phi())

			HH_mass_cand.append(HH_p4.M())
			HH_pt_cand.append(HH_p4.Pt())
			HH_eta_cand.append(HH_p4.Eta())
			HH_phi_cand.append(HH_p4.Phi())

		# =================================================================
		# Define Higgsness and Topness
		eta_cand = np.arange(-3.4,3.4,0.4)

		for i in range(len(jet1_mass)):
			if (lep1_charge[i] == lep2_charge[i]): continue

			chi2_cand,H_cand = [],[]
			for j in range(len(eta_cand)):
				for k in range(j,len(eta_cand)):

					lep1_p4 = TLorentzVector()
                                        lep2_p4 = TLorentzVector()

					jet1_p4 = TLorentzVector()
                                        jet2_p4 = TLorentzVector()

                                        nu1_p4 = TLorentzVector()
                                        nu2_p4 = TLorentzVector()

					if (lep1_charge < 0): # w1 is w_minus

						nu1_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[j],c_met[i][2],0)
						nu2_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[k],c_met[i][2],0)

						lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
						lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])

						w1_p4 = lep1_p4 + nu1_p4
						w2_p4 = lep2_p4 + nu2_p4

						jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][4])
						jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])

						t1_p4 = w1_p4 + jet1_p4 # because we don't know charges of jets, I tried to combine w1, w2 with both jet1 and jet2
						t2_p4 = w1_p4 + jet2_p4 # top quark
						t3_p4 = w2_p4 + jet1_p4 # t-bar
						t4_p4 = w2_p4 + jet2_p4 # t-bar

						t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
						w_cand = [w1_p4.M(),w2_p4.M()]

						sig_t = np.std(t_cand)
						sig_w = np.std(w_cand)

						# 4 tops exist, but 2 terms of 4 terms have same top charges
						# there are 2 of W-boson, so when the first top is assigned, its matched W-boson is also assigned
						# by the definition, only a t-bar can come to the first term and only the W+-boson can be assigned to the second term
						todular1 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
						todular2 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
						todular3 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
						todular4 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)

						chi2_cand.append(todular1)
						chi2_cand.append(todular2)
						chi2_cand.append(todular3)
						chi2_cand.append(todular4)

						h1_p4 = TLorentzVector()
						h1_p4 = w1_p4 + w2_p4

						nunu_p4 = TLorentzVector()
						nunu_p4 = nu1_p4 + nu2_p4

						sig_h = 2
                                                sig_nu = 1
						sig_wask = 5
                                                nunu_peak = 30
                                                wask_peak = (1/np.sqrt(3))*np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*h_mass**2*w_mass**2 + w_mass**4) )
                                                # because the condition is lep1_charge is positive, the lep1 goes to only 4th & 5th term
                                                hodular_prev = ((lep2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep1_p4.M()**2-wask_peak**2)**2/sig_wask**4) + ((lep1_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep2_p4.M()**2-wask_peak**2)**2/sig_wask**4)

						hodular1 = ((h1_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev

						H_cand.append(hodular1)


					if (lep1_charge > 0): # w1 is w_minus

						nu1_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[j],c_met[i][2],0)
						nu2_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[k],c_met[i][2],0)

                                                lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
                                                lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])

                                                w1_p4 = lep2_p4 + nu2_p4
						w2_p4 = lep1_p4 + nu1_p4

						jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][4])
                                                jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])

                                                t1_p4 = w1_p4 + jet1_p4 # because we don't know charges of jets, I tried to combine w1, w2 with both jet1 and jet2
                                                t2_p4 = w1_p4 + jet2_p4 # top quark

                                                t3_p4 = w2_p4 + jet1_p4 # t-bar
                                                t4_p4 = w2_p4 + jet2_p4 # t-bar

                                                t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                                                w_cand = [w1_p4.M(),w2_p4.M()]

                                                sig_t = np.std(t_cand)
                                                sig_w = np.std(w_cand)

                                                # 4 tops exist, but 2 terms of 4 terms have same top charges
                                                # there are 2 of W-boson, so when the first top is assigned, its matched W-boson is also assigned
                                                # from the definition, only a t-bar can come to the first term and only the W+-boson can be assigned to the second term
                                                todular1 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                                todular2 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                                todular3 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                                todular4 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)

                                                chi2_cand.append(todular1)
                                                chi2_cand.append(todular2)
                                                chi2_cand.append(todular3)
                                                chi2_cand.append(todular4)

						h1_p4 = TLorentzVector()
                                                h1_p4 = w1_p4 + w2_p4

                                                nunu_p4 = TLorentzVector()
                                                nunu_p4 = nu1_p4 + nu2_p4

						sig_h = 2
                                                sig_nu = 1
						sig_wask = 5
                                                nunu_peak = 30
                                                wask_peak = (1/np.sqrt(3))*np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*h_mass**2*w_mass**2 + w_mass**4) )
                                                # because the condition is lep1_charge is positive, the lep1 goes to only 3rd & 6th term
                                                hodular_prev = ((lep1_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep2_p4.M()**2-wask_peak**2)**2/sig_wask**4) + ((lep2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep1_p4.M()**2-wask_peak**2)**2/sig_wask**4)

                                                hodular1 = ((h1_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev

                                                H_cand.append(hodular1)

			topness.append(np.log(min(chi2_cand)))
			higgsness.append(np.log(min(H_cand)))


	c = TCanvas("c","Higgsness & Topness 2D Histogram_signal_2016", 900,600)
	c.cd()

	# total signal entry = 385
	hist_ht = TH2D("h2","Higgsness & Topness 2D Histogram_signal_2016 (Working in Progress)",50,0,30,50,0,10)
	hist_ht.GetXaxis().SetTitle("Log H")
	hist_ht.GetYaxis().SetTitle("Log T")

	for i in range(len(topness)): hist_ht.Fill(higgsness[i],topness[i])

	hist_ht.Draw("colz")
	c.SaveAs("HTness_2016_sig.png")

	input("Press Enter to continue...")

if __name__ == "__main__":
	defineHTness()
