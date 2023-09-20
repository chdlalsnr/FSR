from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
from scipy.optimize import minimize


def main():

	# Signal SL
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2WToLNu2J_node_SM_TuneCUETP8M1_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')

	# TTBar 
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')

	ii,hh,ee = 0,0,0
	h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876
	hgg,w_on,w_off = [],[],[]
        tt,ww = [],[]
	data = []

	for fname in flist:

		if ee > 1000000: break
		f = TFile(fname,"read")
		t = f.Get("Events")

		try: event_num = t.GetEntriesFast()
		except AttributeError: continue
		print("Processing %s..."%(fname))
		ii = ii+1
		print("Have processed %i events..."%(ee))

		# Only Mu channel considered at this point
		# HLT conditions
		for i in range(event_num):

                        t.GetEntry(i)
			hlt_mu,hlt_el = False,False
			good_event = False

			# SingleMuon
                        hlt_mu_1 = t.HLT_IsoMu22
                        hlt_mu_2 = t.HLT_IsoTkMu22
                        hlt_mu_3 = t.HLT_IsoMu22_eta2p1
                        hlt_mu_4 = t.HLT_IsoTkMu22_eta2p1
                        hlt_mu_5 = t.HLT_IsoMu24
                        hlt_mu_6 = t.HLT_IsoTkMu24

                        if (hlt_mu_1 == 1 or hlt_mu_2 == 1 or hlt_mu_3 == 1 or hlt_mu_4 == 1 or hlt_mu_5 == 1 or hlt_mu_6 == 1): hlt_mu = True

			# SingleElectron
                        hlt_el_1 = t.HLT_Ele27_WPTight_Gsf
                        hlt_el_2 = t.HLT_Ele25_eta2p1_WPTight_Gsf
                        hlt_el_3 = t.HLT_Ele27_eta2p1_WPLoose_Gsf

                        if (hlt_el_1 == 1 or hlt_el_2 == 1 or hlt_el_3 == 1): hlt_el = True

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

			if ((hlt_mu or hlt_el) and met_filter): good_event = True

			# =================================================================
			#Remained cut : The lepton isolation, defined as the scalar
			#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
			#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
			# Medium muon discrimination

			ee += 1
			if not good_event: continue
			if t.nJet < 4: continue
			pt,phi = 0,0

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0
			ljet_mass,ljet_pt,ljet_eta,ljet_phi,ljet_btag = 0,0,0,0,0
			sjet_mass,sjet_pt,sjet_eta,sjet_phi,sjet_btag = 0,0,0,0,0
			tjet_mass,tjet_pt,tjet_eta,tjet_phi = 0,0,0,0
			fjet_mass,fjet_pt,fjet_eta,fjet_phi = 0,0,0,0

			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                        temp_charge,temp_medid = [],[]
                        temp_dxy,temp_dz,temp_btag = [],[],[]
                        temp_iso,temp_sip3d,temp_mva = [],[],[]

			flag_lep,flag_jet1,flag_jet2 = False,False,False
                        genweight = t.genWeight
			tightMuon = 0
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
				if (t.Muon_tightId.__getitem__(j) == True): tightMuon += 1

			if tightMuon >= 2: continue

			l1 = -1
			for k in range(t.nMuon):
				if (flag_lep == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_medid[k] == True and temp_iso[k] < 0.4 and temp_sip3d[k] < 8 and temp_mva[k] > 0.5):

					llep_mass = temp_mass[k]
					llep_pt   = temp_pt[k]
					llep_eta  = temp_eta[k]
					llep_phi  = temp_phi[k]
					llep_charge = temp_charge[k]

					l1 = k
					flag_lep = True

			# =================================================================
			# Jet discrimination
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                        temp_btag,temp_jetid,temp_charge = [],[],[]
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
					if m is not l1: continue
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

			j1,j2 = -1,-1
			for k in range(t.nJet):
				if (flag_jet1 == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4  and temp_jetid[k] >= 7 and temp_btag[k] > 0.3093): # jetid == 1: Loose, 3: Tight, 7: TightLeptonVeto
                                        for l in range(k+1,t.nJet):
                                                if (flag_jet1 == True): break
                                                if (temp_pt[l] > 25 and abs(temp_eta[l]) < 2.4  and temp_jetid[l] >= 7 and temp_btag[l] > 0.3093):

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
                                                        if (dijet_p4.M() > 0): flag_jet1 = True
							j1,j2 = k,l


			if not flag_jet1: continue
			if flag_jet1:
				temp_mass.pop(j1)
				temp_mass.pop(j2-1)
				temp_pt.pop(j1)
				temp_pt.pop(j2-1)
				temp_eta.pop(j1)
				temp_eta.pop(j2-1)
				temp_phi.pop(j1)
				temp_phi.pop(j2-1)
				temp_btag.pop(j1)
				temp_btag.pop(j2-1)
				temp_jetid.pop(j1)
				temp_jetid.pop(j2-1)

			for k in range(len(temp_mass)-1,0,-1):
				if (flag_jet2 == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and temp_jetid[k] >= 1):
					for l in range(k-1,0,-1):
						if (flag_jet2 == True): break
						if (temp_pt[l] > 25 and abs(temp_eta[l]) < 2.4 and temp_jetid[l] >= 1):

							tjet_mass = temp_mass[k]
							tjet_pt   = temp_pt[k]
							tjet_eta  = temp_eta[k]
							tjet_phi  = temp_phi[k]

							fjet_mass = temp_mass[l]
							fjet_pt   = temp_pt[l]
							fjet_eta  = temp_eta[l]
							fjet_phi  = temp_phi[l]

							jet3_p4 = TLorentzVector()
							jet4_p4 = TLorentzVector()
							jet3_p4.SetPtEtaPhiM(tjet_pt,tjet_eta,tjet_phi,tjet_mass)
							jet4_p4.SetPtEtaPhiM(fjet_pt,fjet_eta,fjet_phi,fjet_mass)
							dijet2_p4 = jet3_p4 + jet4_p4

							if (dijet2_p4.M() == 0): continue
							if (dijet2_p4.M() > 0 and dijet2_p4.M() < 120): flag_jet2 = True

			if not (flag_lep == True and flag_jet1 == True and flag_jet2 == True): continue

			dr_ll = ((llep_eta-slep_eta)**2+(llep_phi-slep_phi)**2)**0.5
                        dr_jj = ((ljet_eta-sjet_eta)**2+(ljet_phi-sjet_phi)**2)**0.5

			# =================================================================
                        # Higgs reconstruction

			lep1_p4 = TLorentzVector()
			lep1_p4.SetPtEtaPhiM(llep_pt,llep_eta,llep_phi,llep_mass)

			jet3_p4 = TLorentzVector()
			jet4_p4 = TLorentzVector()
			jet3_p4.SetPtEtaPhiM(tjet_pt,tjet_eta,tjet_phi,tjet_mass)
			jet4_p4.SetPtEtaPhiM(fjet_pt,fjet_eta,fjet_phi,fjet_mass)

			met_p4 = TLorentzVector()
			met_p4.SetPtEtaPhiM(met_pt,0,met_phi,0)
			h1_p4 = lep1_p4 + met_p4 + jet3_p4 + jet4_p4

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

                        def defineHiggsness(x):

                                nu1_p4.SetPtEtaPhiM(met_pt, x, met_phi, 0)

                                w1_p4 = lep1_p4 + nu1_p4
                                w2_p4 = jet3_p4 + jet4_p4

                                h1_p4 = w1_p4 + w2_p4

				sig_h = 48.2
                                sig_on,sig_off = 20.6,22.4
                                peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

                                inner = [((w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4), ((w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4)]
                                higgsness = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + min(inner)

                                return higgsness

                        ini = np.array([0])
                        bnd = [(-2.4,2.4)]

                        opt = minimize(defineHiggsness, ini, bounds=bnd, method='SLSQP', options={'ftol':20}) #'disp':True})

                        if opt.success: hh += 1

                        nu1_p4.SetPtEtaPhiM(met_pt, opt.x[0], met_phi, 0)

                        w1_p4 = lep1_p4 + nu1_p4
                        w2_p4 = jet3_p4 + jet4_p4

                        h1_p4 = w1_p4 + w2_p4

			def findOnshellW():

				flag_on = False

				sig_h = 48.2
                                sig_on,sig_off = 20.6,22.4
				peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

				opt1 = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4
				opt2 = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4

				if (abs(opt2-opt.fun) >= abs(opt1-opt.fun)): flag_on = True
				if (abs(opt1-opt.fun) > abs(opt2-opt.fun)) : flag_on = False

				return flag_on

                        def defineTopness():

                                chi2 = []

				sig_t = 44.2
                                sig_w = 19.3

                                # because we don't know charges of jets, I tried to combine w1, w2 with both jet1 and jet2
                                t1_p4 = w1_p4 + jet1_p4 # top quark
                                t2_p4 = w1_p4 + jet2_p4 # top quark
                                t3_p4 = w2_p4 + jet1_p4 # t-bar
                                t4_p4 = w2_p4 + jet2_p4 # t-bar

                                t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                                w_cand = [w1_p4.M(),w2_p4.M()]

                                # if taken t1_p4 as the first top quark, the second top quark should be t4_p4
                                # otherwise, if t2_p4 is the first top quark, the second should be t3_p4
                                chi2.append( ((t_cand[0]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[3]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )
                                chi2.append( ((t_cand[1]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[2]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )

                                topness = min(chi2)

				return topness

			def findTopCombination():

				flag_top = False

				sig_t = 44.2
                                sig_w = 19.3

				t1_p4 = w1_p4 + jet1_p4
				t2_p4 = w1_p4 + jet2_p4
				t3_p4 = w2_p4 + jet1_p4
				t4_p4 = w2_p4 + jet2_p4

				t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
				w_cand = [w1_p4.M(),w2_p4.M()]

				opt1 = ((t_cand[0]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[3]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4)
				opt2 = ((t_cand[1]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[2]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4)

				if (abs(opt2-defineTopness()) >= abs(opt1-defineTopness())): flag_top = True
				if (abs(opt1-defineTopness()) > abs(opt2-defineTopness())) : flag_top = False

				return flag_top

			if findOnshellW(): # w1 is onshell W
				ww.append(w1_p4.M())
				if findTopCombination(): # jet1
					t1_p4 = w1_p4 + jet1_p4
					tt.append(t1_p4.M())
				else: # jet2
					t1_p4 = w1_p4 + jet2_p4
					tt.append(t1_p4.M())

			else: # w2 is onshell W
				ww.append(w2_p4.M())
				if findTopCombination():
					t1_p4 = w2_p4 + jet1_p4
					tt.append(t1_p4.M())
				else:
					t1_p4 = w2_p4 + jet2_p4
					tt.append(t1_p4.M())

	c1 = TCanvas("c1","Top Quark Invariant Mass Distribution_SL", 900,600)
	c2 = TCanvas("c2","W Boson Invariant Mass Distribution_SL", 900,600)

	hist_t = TH1D("hist_t","Top Quark Invariant Mass Distribution_SL", 20,0,300)
	hist_w = TH1D("hist_w","W Boson Invariant Mass Distribution_SL", 20,10,130)

	for i in range(len(tt)):
		hist_t.Fill(tt[i])
		hist_w.Fill(ww[i])

	hist_t.GetXaxis().SetTitle("Top Mass [GeV]")
	hist_t.GetYaxis().SetTitle("Entries")

	hist_w.GetXaxis().SetTitle("W Mass [GeV]")
	hist_w.GetYaxis().SetTitle("Entries")

	c1.cd()
	hist_t.Draw()
	c1.SaveAs("t_mass.png")

	c2.cd()
	hist_w.Draw()
	c2.SaveAs("w_mass.png")


#			# =================================================================
#			# Save CSV file
#			data.append([ljet_mass,ljet_pt,ljet_eta,ljet_phi,
#                                     sjet_mass,sjet_pt,sjet_eta,sjet_phi,
#				     tjet_mass,tjet_pt,tjet_eta,tjet_phi,
#				     fjet_mass,fjet_pt,fjet_eta,fjet_phi,
#                                     llep_mass,llep_pt,llep_eta,llep_phi,
#				     dr_ll,dr_jj,met_pt,met_phi,
#				     opt.fun,defineTopness(),1])
#
#        df = pd.DataFrame(data, columns=
#                                ['jet1_mass','jet1_pt','jet1_eta','jet1_phi',
#                                 'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
#                                 'jet3_mass','jet3_pt','jet3_eta','jet3_phi',
#                                 'jet4_mass','jet4_pt','jet4_eta','jet4_phi',
#				 'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
#				 'dr_ll','dr_jj','met_pt','met_phi',
#				 'higgsness','topness','target'])
#
#        df.to_csv("htness_hh.csv", header=True, index=False)
#
#        print("Have processed %i files..."%(ii))
#        print("Have processed total %i events..."%(ee))
#        print("Event number finally passed is %i..."%(hh))

	input("Press Enter to continue...")		


if __name__ == "__main__":
	main()
