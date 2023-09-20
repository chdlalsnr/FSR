from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
from scipy.optimize import minimize,NonlinearConstraint


def defineHTness():

	# Signal
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph*/*/*/*/*')

	# TTBar 
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')

	ii,hh,ee = 0,0,0
	h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876
	hgg,w_on,w_off,dl = [],[],[],[]
	tt,ww = [],[]
	data = []

	for fname in flist:

		if ee > 1000000: break
		f = TFile(fname,"read")
		t = f.Get("Events")

		event_num = t.GetEntriesFast()
		print("Processing %s..."%(fname))
		ii = ii+1
		print("Have processed %i events..."%(ee))

		# Only MuMu channel considered at this point
		# HLT conditions
		for i in range(event_num):

                        t.GetEntry(i)
			hlt_mumu,hlt_elel,hlt_elmu = False,False,False
			hlt_mu,hlt_el = False,False
			good_event = False

			# DoubleMuon
                        hlt_mumu_1 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
                        hlt_mumu_2 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
                        hlt_mumu_3 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
                        hlt_mumu_4 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ

                        if (hlt_mumu_1 == 1 or hlt_mumu_2 == 1 or hlt_mumu_3 == 1 or hlt_mumu_4 == 1): hlt_mumu = True

			# DoubleElectron
			hlt_elel_1 = t.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ

			if (hlt_elel_1 == 1): hlt_elel = True

			# Muon + Electron
                        hlt_elmu_1 = t.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
                        hlt_elmu_2 = t.HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL
                        hlt_elmu_3 = t.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL

                        if (hlt_elmu_1 == 1 or hlt_elmu_2 == 1 or hlt_elmu_3 == 1): hlt_elmu == True

                        # SingleMuon
                        hlt_mu_1 = t.HLT_IsoMu22
                        hlt_mu_2 = t.HLT_IsoMu22_eta2p1
                        hlt_mu_3 = t.HLT_IsoTkMu22
                        hlt_mu_4 = t.HLT_IsoTkMu22_eta2p1
                        hlt_mu_5 = t.HLT_IsoMu24
                        hlt_mu_6 = t.HLT_IsoTkMu24

                        if (hlt_mu_1 == 1 or hlt_mu_2 == 1 or hlt_mu_3 == 1 or hlt_mu_4 == 1 or hlt_mu_5 == 1 or hlt_mu_6 == 1): hlt_mu = True

                        # SingleElectron
                        hlt_el_1 = t.HLT_Ele27_WPTight_Gsf
                        hlt_el_2 = t.HLT_Ele25_eta2p1_WPTight_Gsf
                        hlt_el_3 = t.HLT_Ele27_eta2p1_WPLoose_Gsf

                        if (hlt_el_1 == 1 or hlt_el_2 == 1 or hlt_el_3 == 1): hlt_el = True

			# MET Filter
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

                        if (hlt_mumu or hlt_elel or hlt_elmu or hlt_mu or hlt_el):
                                if met_filter: good_event = True

			# =================================================================
			#Remained cut : The lepton isolation, defined as the scalar
			#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
			#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
			# Medium muon discrimination

			ee = ee+1
			if not good_event: continue
			if (t.nJet < 2): continue
			if (t.nMuon < 2): continue

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0
			ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
			sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0
			dr_ll,dr_jj,met_pt,met_phi = 0,0,0,0

			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                        temp_charge,temp_medid = [],[]
                        temp_dxy,temp_dz,temp_btag = [],[],[]
                        temp_iso,temp_sip3d,temp_mva = [],[],[]

			flag_lep,flag_jet = False,False
                        genweight = t.genWeight
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

			l1,l2 = -1,-1
			for k in range(t.nMuon):
				if (flag_lep == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_medid[k] == True and temp_iso[k] < 0.4 and temp_sip3d[k] < 8 and temp_mva[k] > 0.5):
					for l in range(k+1,t.nMuon):
						if (flag_lep == True): break
						if (temp_pt[l] > 15 and abs(temp_eta[l]) < 2.4 and abs(temp_dxy[l]) < 0.05 and abs(temp_dz[l]) < 0.1 and temp_medid[l] == True and temp_iso[l] < 0.4 and temp_sip3d[l] < 8 and temp_mva[l] > 0.5):
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

								if (dilep_p4.M() > z_mass-5 and dilep_p4.M() < z_mass+5): continue
								if not (dilep_p4.M() > 12): continue
                                                                else: flag_lep = True

								l1 = k
								l2 = l

			# =================================================================
			# Jet discrimination
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                        temp_btag,temp_jetid = [],[]
			met_pt  = t.MET_pt
			met_phi = t.MET_phi

			for j in range(t.nJet): #(t.Jet_mass.__len__())
				temp_mass.append(t.Jet_mass.__getitem__(j))
				temp_pt.append(t.Jet_pt.__getitem__(j))
				temp_eta.append(t.Jet_eta.__getitem__(j))
				temp_phi.append(t.Jet_phi.__getitem__(j))
				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
				temp_jetid.append(t.Jet_jetId.__getitem__(j))

			# Muon veto
			deltaR2, deltaR = [],[]
			for j in range(t.nJet):
				for m in range(t.nMuon):
					if (m != l1): continue
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

			for k in range(t.nJet):
				if (flag_jet == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and temp_btag[k] > 0.3093 and temp_jetid[k] >= 7): # jetid == 1: Loose, 3: Tight, 7: TightLeptonVeto
					for l in range(k+1,t.nJet):
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

							if (dijet_p4.M() == 0): continue
                                                        if (dijet_p4.M() > 0): flag_jet = True

			if not (flag_lep == True and flag_jet == True): continue

			dr_ll = ((llep_eta-slep_eta)**2+(llep_phi-slep_phi)**2)**0.5
                        dr_jj = ((ljet_eta-sjet_eta)**2+(ljet_phi-sjet_phi)**2)**0.5

			# =================================================================
			# Higgs reconstruction

			lep1_p4 = TLorentzVector()
			lep2_p4 = TLorentzVector()

			lep1_p4.SetPtEtaPhiM(llep_pt,llep_eta,llep_phi,llep_mass)
			lep2_p4.SetPtEtaPhiM(slep_pt,slep_eta,slep_phi,slep_mass)

			met_p4 = TLorentzVector()
			met_p4.SetPtEtaPhiM(met_pt,0,met_phi,0)
			h1_p4 = lep1_p4 + lep2_p4 + met_p4

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()

			jet1_p4.SetPtEtaPhiM(ljet_pt,ljet_eta,ljet_phi,ljet_mass)
			jet2_p4.SetPtEtaPhiM(sjet_pt,sjet_eta,sjet_phi,sjet_mass)
			h2_p4 = jet1_p4 + jet2_p4

			hh_p4 = TLorentzVector()
			hh_p4 = h1_p4 + h2_p4

			# =================================================================
			# Define Higgsness and Topness

			nu1_p4 = TLorentzVector()
			nu2_p4 = TLorentzVector()

			def defineHiggsness(x):

				nu1_p4.SetPtEtaPhiM(x[0], x[1], x[2], 0)
                                nu2_p4.SetPtEtaPhiM(x[3], x[4], x[5], 0)

                                w1_p4 = lep1_p4 + nu1_p4
                                w2_p4 = lep2_p4 + nu2_p4

				h1_p4 = w1_p4 + w2_p4
				dilep_p4 = lep1_p4 + lep2_p4

                                sig_h = 30
                                sig_on,sig_off = 30,30
                                sig_lep = 30
                                peak_dilep = 30
                                peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

                                inner = [((w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4), ((w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4)]

                                higgsness = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (dilep_p4.M()**2-peak_dilep**2)**2/sig_lep**4  + min(inner)

                                return higgsness

			ini = np.array([0,0,0, 0,0,0])
			bnd = [(0,100),(-2.4,2.4),(-3.14,3.14), (0,100),(-2.4,2.4),(-3.14,3.14)]

			con_eq = lambda x: x[0] + x[3]
			cons = NonlinearConstraint(con_eq, met_pt, np.inf) 

			opt = minimize(defineHiggsness, ini, constraints=cons, bounds=bnd, method='SLSQP', options={'ftol':20}) #'disp':True})

			if opt.success: hh += 1
                        if not opt.success: continue

                        nu1_p4.SetPtEtaPhiM(opt.x[0], opt.x[1], opt.x[2], 0)
                        nu2_p4.SetPtEtaPhiM(opt.x[3], opt.x[4], opt.x[5], 0)

                        w1_p4 = lep1_p4 + nu1_p4
                        w2_p4 = lep2_p4 + nu2_p4

                        h1_p4 = w1_p4 + w2_p4

                        def findOnshellW():

                                flag_on = False

                                sig_h = 30
                                sig_on,sig_off = 30,30
                                sig_lep = 30
                                peak_dilep = 30
                                peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

                                opt1 = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (dilep_p4.M()**2-peak_dilep**2)**2/sig_lep**4 + (w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4
                                opt2 = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (dilep_p4.M()**2-peak_dilep**2)**2/sig_lep**4 + (w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4

                                if (abs(opt2-opt.fun) >= abs(opt1-opt.fun)): flag_on = True
                                if (abs(opt1-opt.fun) > abs(opt2-opt.fun)) : flag_on = False

                                return flag_on

                        if findOnshellW():
                                w_on.append(w1_p4.M())
                                w_off.append(w2_p4.M())

                        if not findOnshellW():
                                w_on.append(w2_p4.M())
                                w_off.append(w1_p4.M())

			hgg.append(h1_p4.M())
                        dl.append(dilep_p4.M())

        c1 = TCanvas("c1","Higgs Boson Invariant Mass Distribution_DL", 900,600)
        c2 = TCanvas("c2","W_onshell Invariant Mass Distribution_DL", 900,600)
        c3 = TCanvas("c3","W_offshell Invariant Mass Distribution_DL", 900,600)
        c4 = TCanvas("c4","Dilepton Invariant Mass Distribution_DL", 900,600)

        hist_h = TH1D("hist_h","Higgs Boson Invariant Mass Distribution_DL", 20,0,300)
        hist_on = TH1D("hist_on","W_onshell Invariant Mass Distribution_DL", 20,10,130)
        hist_off = TH1D("hist_off","W_offshell Invariant Mass Distribution_DL", 20,10,200)
        hist_dl = TH1D("hist_dl","Dilepton Invariant Mass Distribution_DL", 20,10,200)

        for i in range(len(hgg)):
                hist_h.Fill(hgg[i])
                hist_on.Fill(w_on[i])
                hist_off.Fill(w_off[i])
                hist_dl.Fill(dl[i])

        hist_h.GetXaxis().SetTitle("H Mass [GeV]")
        hist_h.GetYaxis().SetTitle("Entries")

        hist_on.GetXaxis().SetTitle("W_onshell Mass [GeV]")
        hist_on.GetYaxis().SetTitle("Entries")

        hist_off.GetXaxis().SetTitle("W_offshell Mass [GeV]")
        hist_off.GetYaxis().SetTitle("Entries")

        hist_dl.GetXaxis().SetTitle("Dilepton Mass [GeV]")
        hist_dl.GetYaxis().SetTitle("Entries")

        c1.cd()
        hist_h.Draw()
        c1.SaveAs("h_mass.png")

        c2.cd()
        hist_on.Draw()
        c2.SaveAs("w_on_mass.png")

        c3.cd()
        hist_off.Draw()
        c3.SaveAs("w_off_mass.png")

        c4.cd()
        hist_dl.Draw()
        c4.SaveAs("dl_mass.png")

#			def defineTopness(y):
#
#				nu1_p4.SetPtEtaPhiM(y[0], y[1], y[2], 0)
#				nu2_p4.SetPtEtaPhiM(y[3], y[4], y[5], 0)
#
#				w1_p4 = lep1_p4 + nu1_p4
#				w2_p4 = lep2_p4 + nu2_p4
#
#                                chi2 = []
#
#                                # (t1, t2) and (t5, t6) have anti-particle relations
#                                t1_p4 = w1_p4 + jet1_p4
#                                t2_p4 = w1_p4 + jet2_p4
#                                t3_p4 = w2_p4 + jet1_p4
#                                t4_p4 = w2_p4 + jet2_p4
#
#                                t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
#                                w_cand = [w1_p4.M(),w2_p4.M()]
#
#                                sig_t = 30
#                                sig_w = 30
#
#				chi2.append( ((t_cand[0]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[3]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )
#				chi2.append( ((t_cand[1]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[2]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )
#
#                                topness = min(chi2)
#
#                                return topness
#
#			ini = np.array([0,0,0, 0,0,0])
#                        bnd = [(0,100),(-2.4,2.4),(-3.14,3.14), (0,100),(-2.4,2.4),(-3.14,3.14)]
#
#                        con_eq = lambda y: y[0] + y[3]
#                        cons = NonlinearConstraint(con_eq, met_pt, np.inf)
#
#                        opt2 = minimize(defineTopness, ini, constraints=cons, bounds=bnd, method='SLSQP', options={'ftol':5}) #'disp':True})
#
#                        if opt.success and opt2.success: hh += 1
#                        else: continue
#
#			# =================================================================
#			# Save CSV file
#        		data.append([ljet_mass,ljet_pt,ljet_eta,ljet_phi,
#				     sjet_mass,sjet_pt,sjet_eta,sjet_phi,
#				     llep_mass,llep_pt,llep_eta,llep_phi,
#				     slep_mass,slep_pt,slep_eta,slep_phi,
#				     dr_ll,dr_jj,met_pt,met_phi,
#				     opt.fun,opt2.fun,1])
#
#	df = pd.DataFrame(data, columns=
#			['jet1_mass','jet1_pt','jet1_eta','jet1_phi',
#			'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
#			'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
#			'lep2_mass','lep2_pt','lep2_eta','lep2_phi',
#			'dr_ll','dr_jj','met_pt','met_phi',
#			'higgsness','topness','target'])
#
#	df.to_csv("htness_hh.csv", header=True, index=False)
#
#        print("Have processed %i files..."%(ii))
#        print("Have processed total %i events..."%(ee))
#        print("Event number finally passed is %i..."%(hh))

	input("Press Enter to continue...")


if __name__ == "__main__":
	defineHTness()
