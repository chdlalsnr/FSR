from ROOT import *
from multiprocessing import Process, Queue
import numpy as np
import math
import glob
import sys


def main():
#	# Signal
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/BC8BEF91-F53A-A641-BB54-EFF05EDD6CB4.root')

        # TTBar 
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist.extend(flist_had)
#        flist.extend(flist_semi)

#        # DY_M10-50
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/*/*')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*.root')
#        flist  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#        flist.extend(flist_m50)
#        flist.extend(flist_0j)
#        flist.extend(flist_1j)
#        flist.extend(flist_2j)

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


	ii,hh,ee = 0,0,0
	c_dilep,c_dijet = [],[]
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
		Higgs_mass = 125.18
		for i in range(event_num):

			ee = ee+1
			t.GetEntry(i)
			if i not in hlt_pass: continue

			llep_mass,llep_pt,llep_eta,llep_phi = 0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi = 0,0,0,0

			ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
			sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0
			llep_mass_l,llep_pt_l,llep_eta_l,llep_phi_l = 0,0,0,0
			slep_mass_l,slep_pt_l,slep_eta_l,slep_phi_l = 0,0,0,0

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
								llep_phi  =temp_phi[k]

								slep_mass = temp_mass[l]
								slep_pt   = temp_pt[l]
								slep_eta  = temp_eta[l]
								slep_phi  = temp_phi[l]

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
							mass_cand.append(abs(jet_p4.M()-Higgs_mass))
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

							if (abs(jet_p4.M()-Higgs_mass)==higgs_select):

								ljet_mass = temp_mass[k]
								ljet_pt   = temp_pt[k]
								ljet_eta  = temp_eta[k]
								ljet_phi  = temp_phi[k]

								sjet_mass = temp_mass[l]
								sjet_pt   = temp_pt[l]
								sjet_eta  = temp_eta[l]
								sjet_phi  = temp_phi[l]

								flag_jet = True
						break
				break

			dilep,dijet = [],[]
			if (flag_lep == True and flag_jet == True):

				dijet.append(ljet_mass)
				dijet.append(ljet_pt)
				dijet.append(ljet_eta)
				dijet.append(ljet_phi)
				dijet.append(sjet_mass)
				dijet.append(sjet_pt)
				dijet.append(sjet_eta)
				dijet.append(sjet_phi)

				dilep.append(llep_mass)
				dilep.append(llep_pt)
				dilep.append(llep_eta)
				dilep.append(llep_phi)
				dilep.append(slep_mass)
				dilep.append(slep_pt)
				dilep.append(slep_eta)
				dilep.append(slep_phi)

				c_dilep.append(dilep)
				c_dijet.append(dijet)

		"""
			# Muons passed jet event selection 
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
		z_mass = 91.19
		# Higgs reconstruction
		for i in range(len(c_dilep)):

			lep1_p4 = TLorentzVector()
			lep2_p4 = TLorentzVector()
			lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
			lep2_p4.SetPtEtaPhiM(c_dilep[i][5],c_dilep[i][6],c_dilep[i][7],c_dilep[i][4])
			#lep1_p4.SetPtEtaPhiM(llep_pt[i],llep_eta[i],llep_phi[i],llep_mass[i])
			#lep2_p4.SetPtEtaPhiM(slep_pt[i],slep_eta[i],slep_phi[i],slep_mass[i])
			lep_p4 = lep1_p4 + lep2_p4

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()
			jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][0])
			jet2_p4.SetPtEtaPhiM(c_dijet[i][5],c_dijet[i][6],c_dijet[i][7],c_dijet[i][4])
			#jet1_p4.SetPtEtaPhiM(ljet_pt[i],ljet_eta[i],ljet_phi[i],ljet_mass[i])
			#jet2_p4.SetPtEtaPhiM(sjet_pt[i],sjet_eta[i],sjet_phi[i],sjet_mass[i])
			jet_p4 = jet1_p4 + jet2_p4

			if (lep_p4.M() > 12. and lep_p4.M() < 76.):

				H1_mass_cand.append(lep_p4.M())
				H1_pt_cand.append(lep_p4.Pt())
				H1_eta_cand.append(lep_p4.Eta())
				H1_phi_cand.append(lep_p4.Phi())

				H2_mass_cand.append(jet_p4.M())
				H2_pt_cand.append(jet_p4.Pt())
				H2_eta_cand.append(jet_p4.Eta())
				H2_phi_cand.append(jet_p4.Phi())

		# =================================================================
		# Plotting
		H1_mass = np.zeros(len(H1_mass_cand))
		H1_pt   = np.zeros(len(H1_pt_cand))
		H1_eta  = np.zeros(len(H1_eta_cand))
		H1_phi  = np.zeros(len(H1_phi_cand))

		for i in range(len(H1_mass)):
			H1_mass[i] = H1_mass_cand[i]
			H1_pt[i]   = H1_pt_cand[i]
			H1_eta[i]  = H1_eta_cand[i]
			H1_phi[i]  = H1_phi_cand[i]

			# considered both bbWW & bbZZ
			if (H1_mass[i] != 0 and H1_pt[i] != 0):
				h_H1_mass.Fill(H1_mass[i])
				h_H1_pt.Fill(H1_pt[i])
				h_H1_eta.Fill(H1_eta[i])
				h_H1_phi.Fill(H1_phi[i])

		H2_mass = np.zeros(len(H2_mass_cand))
		H2_pt   = np.zeros(len(H2_pt_cand))
		H2_eta  = np.zeros(len(H2_eta_cand))
		H2_phi  = np.zeros(len(H2_phi_cand))

		for i in range(len(H2_mass)):
			H2_mass[i] = H2_mass_cand[i]
			H2_pt[i]   = H2_pt_cand[i]
			H2_eta[i]  = H2_eta_cand[i]
			H2_phi[i]  = H2_phi_cand[i]

			if (H2_mass[i] != 0 and H2_pt[i] != 0):
				h_H2_mass.Fill(H2_mass[i])
				h_H2_pt.Fill(H2_pt[i])
				h_H2_eta.Fill(H2_eta[i])
				h_H2_phi.Fill(H2_phi[i])

		f.Close()

#		result.put(H1_mass,H1_pt,H1_eta,H1_phi,H2_mass,H2_pt,H2_eta,H2_phi,ii,ee,hh)
#		return

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

	dr = 'dy_m50'
        ff = TFile("./output/%s.root"%(dr),"recreate")
        ff.cd()
        gDirectory.mkdir(dr)
        ff.cd(dr)

	h_H1_mass.Write()
	h_H1_pt.Write()
	h_H1_eta.Write()
	h_H1_phi.Write()
	
	h_H2_mass.Write()
	h_H2_pt.Write()
	h_H2_eta.Write()
	h_H2_phi.Write()
	
	ff.Close()

	print("Have processed %i files..."%(ii))
	print("Have processed total %i events..."%(ee))
	print("Event number passed hlt is %i..."%(hh))
	input("Press Enter to continue...")

if __name__ == "__main__":
	main()

#	# Signal
#	flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/BC8BEF91-F53A-A641-BB54-EFF05EDD6CB4.root') 
#
#	# TTBar 
#	flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root') 
#	flist_had  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#	flist_semi = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#	flist.extend(flist_had)
#	flist.extend(flist_semi)
#
#	# DY_M10-50
#	flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#	flist_m50 = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/*/*')
#	flist_0j  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*.root')
#	flist_1j  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/*/*.root')
#	flist_2j  = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root')
#	flist.extend(flist_m50)
#	flist.extend(flist_0j)
#	flist.extend(flist_1j)
#	flist.extend(flist_2j)
#
#	h_H1_mass = TH1D("H1_mass", "Mass Distribution of H1 <- MuMu", 40, 0, 200)
#	h_H1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
#	h_H1_pt   = TH1D("H1_pt", "Transverse Momentum Distribution of H1 <- MuMu", 40, 0, 400)
#	h_H1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
#	h_H1_eta  = TH1D("H1_eta", "Eta Distribution of H1 <- MuMu", 40, -2.5, 2.5)
#	h_H1_eta.GetXaxis().SetTitle("Eta")
#	h_H1_phi  = TH1D("H1_phi", "Phi Distribution of H1 <- MuMu", 40, -3.14, 3.14)
#	h_H1_phi.GetXaxis().SetTitle("Phi")
#
#	h_H2_mass = TH1D("H2_mass", "Mass Distribution of H2 <- bb", 40, 0, 200)
#	h_H2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
#	h_H2_pt   = TH1D("H2_pt", "Transverse Momentum Distribution of H2 <- bb", 40, 0, 400)
#	h_H2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
#	h_H2_eta  = TH1D("H2_eta", "Eta Distribution of H2 <- bb", 40, -2.5, 2.5)
#	h_H2_eta.GetXaxis().SetTitle("Eta")
#	h_H2_phi  = TH1D("H2_phi", "Phi Distribution of H2 <- bb", 40, -3.14, 3.14)
#	h_H2_phi.GetXaxis().SetTitle("Phi")
#
#	ii,ee,hh = 0,0,0
#
#	START,END = 0,len(flist)
#	result = Queue()
#	th1 = Process(target=main, args=(1, START, END//2, result))
#	th2 = Process(target=main, args=(2, END//2, END, result))
#
#	th1.start()
#	th2.start()
#	th1.join()
#	th2.join()
#
#	result.put('STOP')
#	while True:
#		tmp = result.get()
#		if (tmp == 'STOP'): break
#		else:
#			for i in range(len(tmp.H1_mass)):
#				h_H1_mass.Fill(tmp.H1_mass[i])
#				h_H1_pt.Fill(tmp.H1_pt[i])
#				h_H1_eta.Fill(tmp.H1_eta[i])
#				h_H1_phi.Fill(tmp.H1_phi[i])
#			for i in range(len(tmp.H2_mass)):
#				h_H2_mass.Fill(tmp.H2_mass[i])
#				h_H2_pt.Fill(tmp.H2_pt[i])
#				h_H2_eta.Fill(tmp.H2_eta[i])
#				h_H2_phi.Fill(tmp.H2_phi[i])
#
#			ii += tmp.ii
#			ee += tmp.ee
#			hh += tmp.hh
#
#	dr = 'ttbar'
#        c1 = TCanvas()
#        c1.Divide(2,2)
#        c2 = TCanvas()
#        c2.Divide(2,2)
#
#        c1.cd(1)
#        h_H1_mass.Draw()
#        c1.cd(2)
#        h_H1_pt.Draw()
#        c1.cd(3)
#        h_H1_eta.Draw()
#        c1.cd(4)
#        h_H1_phi.Draw()
#        c1.SaveAs("Higgs1_%s.png"%(dr))
#
#        c2.cd(1)
#        h_H2_mass.Draw()
#        c2.cd(2)
#        h_H2_pt.Draw()
#        c2.cd(3)
#        h_H2_eta.Draw()
#        c2.cd(4)
#        h_H2_phi.Draw()
#        c2.SaveAs("Higgs2_%s.png"%(dr))
#
#        h_H1_mass.Write()
#        h_H1_pt.Write()
#        h_H1_eta.Write()
#        h_H1_phi.Write()
#        
#        h_H2_mass.Write()
#        h_H2_pt.Write()
#        h_H2_eta.Write()
#        h_H2_phi.Write()
#        
#        ff.Close()
#
#	print("Finished whole processes...")
