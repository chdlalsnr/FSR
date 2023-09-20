from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
from scipy.optimize import minimize


def defineHTness():

	# Signal
#        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/*')

	# TTBar 
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*')

	ii,hh,ee = 0,0,0
	h_mass,t_mass,w_mass,z_mass = 172.26,80.379,91.1876,125.18
	higgsness,topness = [],[]


	H1_mass_cand,H1_pt_cand,H1_eta_cand,H1_phi_cand = [],[],[],[]
        H2_mass_cand,H2_pt_cand,H2_eta_cand,H2_phi_cand = [],[],[],[]
        HH_mass_cand,HH_pt_cand,HH_eta_cand,HH_phi_cand = [],[],[],[]
        jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag = [],[],[],[],[]
        jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag = [],[],[],[],[]
        lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge = [],[],[],[],[]
        lep2_mass,lep2_pt,lep2_eta,lep2_phi,lep2_charge = [],[],[],[],[]
        dr_ll,dr_jj,met_pt,met_phi = [],[],[],[]
	for fname in flist:

		if ee > 1000000 or hh > 500: break
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
                        hlt_mumu,hlt_mu,hlt_el = False,False,False

                        # DoubleMuon
                        hlt_mumu_1 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
                        hlt_mumu_2 = t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
                        hlt_mumu_3 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
                        hlt_mumu_4 = t.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ

                        if (hlt_mumu_1 == 1 or hlt_mumu_2 == 1 or hlt_mumu_3 == 1 or hlt_mumu_4 == 1): hlt_mumu = True

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

			# MET Filter
                        flag_good  = t.Flag_goodVertices
                        flag_halo  = t.Flag_globalSuperTightHalo2016Filter
                        flag_hbhen = t.Flag_HBHENoiseFilter
                        flag_hbiso = t.Flag_HBHENoiseIsoFilter
                        flag_dead  = t.Flag_EcalDeadCellTriggerPrimitiveFilter
                        flag_badpf = t.Flag_BadPFMuonFilter
                        flag_ecbad = t.Flag_ecalBadCalibFilter
                        flag_eebad = t.Flag_eeBadScFilter

                        hlt_vertex = False
                        if (flag_good == 1 and flag_halo == 1 and flag_hbhen == 1 and flag_hbiso == 1 and flag_dead == 1 and flag_badpf == 1 and flag_ecbad == 1):
                                hlt_vertex = True

                        if (hlt_mumu or hlt_mu or hlt_el):
                                if hlt_vertex: hlt_pass.append(i)

		# =================================================================
		#Remained cut : The lepton isolation, defined as the scalar
		#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
		#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
		# Medium muon discrimination
		for i in range(event_num):

			ee = ee+1
                        t.GetEntry(i)
                        if i not in hlt_pass: continue
                        if (t.nJet < 2): continue
                        if (t.nMuon < 2): continue
                        pt,phi = 0,0

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			slep_mass,slep_pt,slep_eta,slep_phi,slep_charge = 0,0,0,0,0
			ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
			sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0

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

								if not (dilep_p4.M() > 12): continue
                                                                if (dilep_p4.M() > 12): flag_lep = True

								l1 = k
								l2 = l


			# =================================================================
			# Jet discrimination
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
                        temp_btag,temp_jetid,temp_charge = [],[],[]

			for j in range(t.nJet): #(t.Jet_mass.__len__())
				temp_mass.append(t.Jet_mass.__getitem__(j))
				temp_pt.append(t.Jet_pt.__getitem__(j))
				temp_eta.append(t.Jet_eta.__getitem__(j))
				temp_phi.append(t.Jet_phi.__getitem__(j))
				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
				temp_jetid.append(t.Jet_jetId.__getitem__(j))


			pt   = t.MET_pt
			phi  = t.MET_phi

			# Muon veto
			deltaR2, deltaR = [],[]
			for j in range(t.nJet):
				for m in range(t.nMuon):
					if (m != l1): continue
					deltaR2 = (temp_eta[j]-t.Muon_eta.__getitem__(m))**2 + (temp_phi[j]-t.Muon_phi.__getitem__(m))**2
					deltaR  = math.sqrt(deltaR2)
					if (deltaR < 0.4):
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
						temp_btag.pop(j)
                                                temp_btag.insert(j,0)
						temp_jetid.pop(j)
						temp_jetid.insert(j,0)


			for k in range(t.nJet):
				if (flag_jet == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4  and temp_jetid[k] >= 7 and temp_btag[k] > 0.3093): # jetid == 1: Loose, jetid == 2: Medium
					for l in range(k+1,t.nJet):
						if (flag_jet == True): break
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
                                                        if (dijet_p4.M() > 0): flag_jet = True

			# =================================================================
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

				met_pt.append(pt)
				met_phi.append(phi)

	# =================================================================
	# Define Higgsness and Topness
	data = []
	for i in range(len(jet1_mass)):

		lep1_p4 = TLorentzVector()
                lep2_p4 = TLorentzVector()

		jet1_p4 = TLorentzVector()
                jet2_p4 = TLorentzVector()

		lep1_p4.SetPtEtaPhiM(lep1_pt[i],lep1_eta[i],lep1_phi[i],lep1_mass[i])
		lep2_p4.SetPtEtaPhiM(lep2_pt[i],lep2_eta[i],lep2_phi[i],lep2_mass[i])

		jet1_p4.SetPtEtaPhiM(jet1_pt[i],jet1_eta[i],jet1_phi[i],jet1_mass[i])
                jet2_p4.SetPtEtaPhiM(jet2_pt[i],jet2_eta[i],jet2_phi[i],jet2_mass[i])

                nu1_p4 = TLorentzVector()
                nu2_p4 = TLorentzVector()

		pt = met_pt[i]
		phi = met_phi[i]

		sig_t = 46.7
                sig_w = 20.6

		chi2 = []
		def defineTopness(x):

                        chi2 = []

                        pt1,eta1,phi1 = x[0],x[1],x[2]
                        pt2,eta2,phi2 = x[3],x[4],x[5]

                        nu1_p4.SetPtEtaPhiM(x[0], x[1], x[2], 0)
                        nu2_p4.SetPtEtaPhiM(x[3], x[4], x[5], 0)

                        # w1 and w2 have different electric charges
                        w1_p4 = lep1_p4 + nu1_p4
                        w2_p4 = lep2_p4 + nu2_p4

                        # (t1, t2) and (t5, t6) have anti-particle relations
                        t1_p4 = w1_p4 + jet1_p4
                        t2_p4 = w1_p4 + jet2_p4
                        t3_p4 = w2_p4 + jet1_p4
                        t4_p4 = w2_p4 + jet2_p4

                        t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                        w_cand = [w1_p4.M(),w2_p4.M()]

                        chi2.append( ((t_cand[0]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[3]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )
                        chi2.append( ((t_cand[1]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[2]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )

                        topness = min(chi2)

                        return topness

		ini = np.array([0,0,0, 0,0,0])
                bnd = [(0,50),(-2.4,2.4),(-3.14,3.14), (0,50),(-2.4,2.4),(-3.14,3.14)]

                cons = ({'type':'eq',
                         'fun': lambda x: x[0]+x[3] - pt},
                        {'type':'eq',
                         'fun': lambda x: x[2]+x[5] - phi})

                opt = minimize(defineTopness, ini, bounds=bnd, method='SLSQP', constraints=cons, options={'ftol':5}) #'disp':True

                if opt.success:
                        hh += 1
                        topness.append(opt.fun)

                if not opt.success: continue

		cc = -1
		for j in range(len(chi2)):
			if topness == chi2[j]: cc = j

		nu1_p4.SetPtEtaPhiM(opt.x[0],opt.x[1],opt.x[2], 0)
                nu2_p4.SetPtEtaPhiM(opt.x[3],opt.x[4],opt.x[5], 0)

		w1_p4 = lep1_p4 + nu1_p4
                w2_p4 = lep2_p4 + nu2_p4

                t1_p4 = w1_p4 + jet1_p4
                t2_p4 = w1_p4 + jet2_p4
                t3_p4 = w2_p4 + jet1_p4
                t4_p4 = w2_p4 + jet2_p4

                t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                w_cand = [w1_p4.M(),w2_p4.M()]

		a,b,c,d = 0,0,0,0
		if cc == 0:
			a = ((t_cand[0]**2-t_mass**2)**2/sig_t**4)
			b = ((t_cand[3]**2-t_mass**2)**2/sig_t**4)
		if cc == 1:
			a = ((t_cand[1]**2-t_mass**2)**2/sig_t**4)
			b = ((t_cand[2]**2-t_mass**2)**2/sig_t**4)

		c = ((w_cand[0]**2-w_mass**2)**2/sig_w**4)
		d = ((w_cand[1]**2-w_mass**2)**2/sig_w**4)

                dilep_p4 = lep1_p4 + lep2_p4

                w1_p4 = lep1_p4 + nu1_p4
                w2_p4 = lep2_p4 + nu2_p4

                h1_p4 = w1_p4 + w2_p4

                sig_h  = 46.9
                sig_on,sig_off = 21.0,25.4
		sig_lep = 20.1
                peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )
		peak_dilep = 30

                inner = [((w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4), ((w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4)]
                higgsness = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (dilep_p4.M()**2-peak_dilep**2)**2/sig_lep**4  + min(inner)

		e,f,g,h = 0,0,0,0
		e = (w1_p4.M()**2-w_mass**2)**2/sig_on**4
		f = (w2_p4.M()**2-peak_off**2)**2/sig_off**4
		g = (w2_p4.M()**2-w_mass**2)**2/sig_on**4
		h = (w1_p4.M()**2-peak_off**2)**2/sig_off**4


		# Save CSV file
                data.append([jet1_p4.M(),jet1_p4.Pt(),jet1_p4.Eta(),jet1_p4.Phi(),
                        jet2_p4.M(),jet2_p4.Pt(),jet2_p4.Eta(),jet2_p4.Phi(),
                        lep1_p4.M(),lep1_p4.Pt(),lep1_p4.Eta(),lep1_p4.Phi(),
                        lep2_p4.M(),lep2_p4.Pt(),lep2_p4.Eta(),lep2_p4.Phi(),
                        a,b,c,d,e,f,g,h,0])

        df = pd.DataFrame(data, columns=
                                ['jet1_mass','jet1_pt','jet1_eta','jet1_phi',
                                 'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
                                 'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
                                 'lep2_mass','lep2_pt','lep2_eta','lep2_phi',
                                 'T_1','T_2','T_3','T_4','H_1','H_2','H_3','H_4','target'])

        df.to_csv("htness_tt.csv", header=False, index=False)

        print("Have processed %i files..."%(ii))
        print("Have processed total %i events..."%(ee))
        print("Event number finally passed is %i..."%(hh))

	input("Press Enter to continue...")



if __name__ == "__main__":
	defineHTness()
