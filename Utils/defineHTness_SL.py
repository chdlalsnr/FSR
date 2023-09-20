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
        higgsness,topness = [],[]

	H1_mass_cand,H1_pt_cand,H1_eta_cand,H1_phi_cand = [],[],[],[]
        H2_mass_cand,H2_pt_cand,H2_eta_cand,H2_phi_cand = [],[],[],[]
        HH_mass_cand,HH_pt_cand,HH_eta_cand,HH_phi_cand = [],[],[],[]
        jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag = [],[],[],[],[]
        jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag = [],[],[],[],[]
        jet3_mass,jet3_pt,jet3_eta,jet3_phi = [],[],[],[]
        jet4_mass,jet4_pt,jet4_eta,jet4_phi = [],[],[],[]
        lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge = [],[],[],[],[]
        dr_ll,dr_jj,genWeight = [],[],[]
        nu_pt, nu_phi = [],[]
	for fname in flist:

		if ee > 500000 or hh > 500: break
		f = TFile(fname,"read")
		t = f.Get("Events")

		try: event_num = t.GetEntriesFast()
		except AttributeError: continue
		print("Processing %s..."%(fname))
		ii = ii+1
		print("Have processed %i events..."%(ee))

		# Only MuMu channel considered at this point
		# HLT pass
		hlt_pass = []
		for i in range(event_num):
                        t.GetEntry(i)
                        hlt_el_1 = t.HLT_Ele25_eta2p1_WPTight_Gsf
                        hlt_el_2 = t.HLT_Ele27_WPTight_Gsf
                        hlt_el_3 = t.HLT_Ele27_eta2p1_WPLoose_Gsf

                        hlt_singleEl = False
                        if (hlt_el_1 == 1 or hlt_el_2 == 1 or hlt_el_3 == 1): hlt_singleEl = True

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

			if hlt_singleMu == True:
				if hlt_vertex == True:
					hlt_pass.append(i)

		# =================================================================
		#Remained cut : The lepton isolation, defined as the scalar
		#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
		#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
		# Medium muon discrimination
		for i in range(event_num):

			ee += 1
			t.GetEntry(i)
			if i not in hlt_pass: continue
			if t.nJet < 4: continue

			llep_mass,llep_pt,llep_eta,llep_phi,llep_charge = 0,0,0,0,0
			lele_mass,lele_pt,lele_eta,lele_phi,lele_charge = 0,0,0,0,0
			ljet_btag,ljet_mass,ljet_pt,ljet_eta,ljet_phi = 0,0,0,0,0
			sjet_btag,sjet_mass,sjet_pt,sjet_eta,sjet_phi = 0,0,0,0,0
			lq_mass,lq_pt,lq_eta,lq_phi = 0,0,0,0
			sq_mass,sq_pt,sq_eta,sq_phi = 0,0,0,0
                	met_pt,met_phi = 0,0

			temp_mass,temp_pt,temp_eta,temp_phi,temp_charge = [],[],[],[],[]
			temp_dxy,temp_dz,temp_medid = [],[],[]

			flag_lep,flag_jet = False,False
                        genweight = t.genWeight
			tight = 0
			for j in range(t.nMuon):
				temp_mass.append(t.Muon_mass.__getitem__(j))
				temp_pt.append(t.Muon_pt.__getitem__(j))
				temp_eta.append(t.Muon_eta.__getitem__(j))
				temp_phi.append(t.Muon_phi.__getitem__(j))
				temp_dxy.append(t.Muon_dxy.__getitem__(j))
				temp_dz.append(t.Muon_dz.__getitem__(j))
				temp_charge.append(t.Muon_charge.__getitem__(j))
				if (t.Muon_tightId.__getitem__(j) == True): tight += 1
				temp_medid.append(t.Muon_mediumId.__getitem__(j))

			if tight >= 2: continue

			ll = 100
			for k in range(t.nMuon):
				if (flag_lep == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) <= 2.4 and abs(temp_dxy[k]) <= 0.05 and abs(temp_dz[k]) <= 0.1 and temp_medid[k] == True):

					llep_mass = temp_mass[k]
					llep_pt   = temp_pt[k]
					llep_eta  = temp_eta[k]
					llep_phi  = temp_phi[k]

					ll = k
					flag_lep = True
				break

			# =================================================================
			# Jet discrimination
			temp_mass,temp_pt,temp_eta,temp_phi = [],[],[],[]
			temp_btag,temp_jetid,temp_puid,temp_csv = [],[],[],[]
			temp_charge,temp_charge_e = [],[]

			for j in range(t.nJet): #(t.Jet_mass.__len__())
				temp_mass.append(t.Jet_mass.__getitem__(j))
				temp_pt.append(t.Jet_pt.__getitem__(j))
				temp_eta.append(t.Jet_eta.__getitem__(j))
				temp_phi.append(t.Jet_phi.__getitem__(j))
				temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
				temp_jetid.append(t.Jet_jetId.__getitem__(j))


			met_pt   = t.MET_pt
			met_phi  = t.MET_phi

			# Muon veto
			deltaR2, deltaR = [],[]
			for j in range(t.nJet):
				for m in range(t.nMuon):
					if m is not ll: continue
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

			flag_jet1,flag_jet2 = False,False
			j1,j2 = 0,0
			for k in range(t.nJet):
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and temp_btag[k] > 0.3093 and temp_jetid[k] == 7):
					for l in range(k+1,t.nJet):
						if (temp_pt[l] > 25 and abs(temp_eta[l]) < 2.4 and temp_btag[l] > 0.3093 and temp_jetid[l] == 7):

							ljet_mass = temp_mass[k]
							ljet_pt   = temp_pt[k]
							ljet_eta  = temp_eta[k]
							ljet_phi  = temp_phi[k]

							sjet_mass = temp_mass[l]
							sjet_pt   = temp_pt[l]
							sjet_eta  = temp_eta[l]
							sjet_phi  = temp_phi[l]

							j1 = k
							j2 = l

							flag_jet1 = True
						break
				break

			temp_mass.pop(j1)
			temp_mass.pop(j2)
			temp_pt.pop(j1)
			temp_pt.pop(j2)
			temp_eta.pop(j1)
			temp_eta.pop(j2)
			temp_phi.pop(j1)
			temp_phi.pop(j2)
			temp_jetid.pop(j1)
			temp_jetid.pop(j2)

			if not flag_jet1: continue
			for k in range(len(temp_mass)-1,0,-1):
				if (temp_pt[k] > 0 and abs(temp_eta[k]) < 2.4 and temp_jetid[k] == 7):
					for l in range(k-1,0,-1):
						if (temp_pt[l] > 0 and abs(temp_eta[l]) < 2.4 and temp_jetid[l] == 7):

							lq_mass = temp_mass[k]
							lq_pt   = temp_pt[k]
							lq_eta  = temp_eta[k]
							lq_phi  = temp_phi[k]

							sq_mass = temp_mass[l]
							sq_pt   = temp_pt[l]
							sq_eta  = temp_eta[l]
							sq_phi  = temp_phi[l]

							flag_jet2 = True
						break
				break

			if (flag_jet1 == True and flag_jet2 == True): flag_jet = True
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

				jet3_mass.append(lq_mass)
				jet3_pt.append(lq_pt)
				jet3_eta.append(lq_eta)
				jet3_phi.append(lq_phi)

				jet4_mass.append(sq_mass)
                                jet4_pt.append(sq_pt)
                                jet4_eta.append(sq_eta)
                                jet4_phi.append(sq_phi)

				lep1_mass.append(llep_mass)
                                lep1_pt.append(llep_pt)
                                lep1_eta.append(llep_eta)
                                lep1_phi.append(llep_phi)
				lep1_charge.append(llep_charge)

				nu_pt.append(met_pt)
				nu_phi.append(met_phi)


	# =================================================================
	data = []
	# Higgs reconstruction
	for i in range(len(jet1_mass)):

		lep1_p4 = TLorentzVector()
		lep1_p4.SetPtEtaPhiM(lep1_pt[i],lep1_eta[i],lep1_phi[i],lep1_mass[i])

		jet1_p4 = TLorentzVector()
		jet2_p4 = TLorentzVector()
		jet1_p4.SetPtEtaPhiM(jet1_pt[i],jet1_eta[i],jet1_phi[i],jet1_mass[i])
		jet2_p4.SetPtEtaPhiM(jet2_pt[i],jet2_eta[i],jet2_phi[i],jet2_mass[i])

		jet3_p4 = TLorentzVector()
		jet4_p4 = TLorentzVector()
		jet3_p4.SetPtEtaPhiM(jet3_pt[i],jet3_eta[i],jet3_phi[i],jet3_mass[i])
		jet4_p4.SetPtEtaPhiM(jet4_pt[i],jet4_eta[i],jet4_phi[i],jet4_mass[i])

		pt  = nu_pt[i]
		phi = nu_phi[i]

		cc = -1
		def defineTopness(x):

			chi2 = []

			eta = x
			nu_p4 = TLorentzVector()
			nu_p4.SetPtEtaPhiM(pt, x, phi, 0)

			w1_p4 = lep1_p4 + nu_p4
			w2_p4 = jet3_p4 + jet4_p4

			# because we don't know charges of jets, I tried to combine w1, w2 with both jet1 and jet2
			t1_p4 = w1_p4 + jet1_p4 # top quark
			t2_p4 = w1_p4 + jet2_p4 # top quark
			t3_p4 = w2_p4 + jet1_p4 # t-bar
			t4_p4 = w2_p4 + jet2_p4 # t-bar

			t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
			w_cand = [w1_p4.M(),w2_p4.M()]

			sig_t = 45.0
			sig_w = 21.0

			# if taken t1_p4 as the first top quark, the second top quark should be t4_p4
			# otherwise, if t2_p4 is the first top quark, the second should be t3_p4
			chi2.append( ((t_cand[0]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[3]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )
			chi2.append( ((t_cand[1]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[2]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )

			topness = min(chi2)

			return topness


		ini = np.array([0])
		bnd = [(-2.4,2.4)]
#		from scipy.optimize import Bounds
#		Bounds(-2.4,2.4)

		opt = minimize(defineTopness, ini, bounds=bnd, method='SLSQP', options={'ftol':5}) #options={'disp':True})

		if opt.success:
			hh += 1
			topness.append(opt.fun)

		if not opt.success: continue

		nu_p4 = TLorentzVector()
                nu_p4.SetPtEtaPhiM(pt, opt.x[0], phi, 0)

                w1_p4 = lep1_p4 + nu_p4
                w2_p4 = jet3_p4 + jet4_p4

		h1_p4 = w1_p4 + w2_p4

		def defineHiggsness():

			sig_h  = 47.6
			sig_on,sig_off = 22.8,30.6
			peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

			inner = [((w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4), ((w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4)]
			higgsness = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + min(inner)

			return higgsness


		# Save CSV file
		data.append([jet1_p4.M(),jet1_p4.Pt(),jet1_p4.Eta(),jet1_p4.Phi(),
			jet2_p4.M(),jet2_p4.Pt(),jet2_p4.Eta(),jet2_p4.Phi(),
			jet3_p4.M(),jet3_p4.Pt(),jet3_p4.Eta(),jet3_p4.Phi(),
			jet4_p4.M(),jet4_p4.Pt(),jet4_p4.Eta(),jet4_p4.Phi(),
			lep1_p4.M(),lep1_p4.Pt(),lep1_p4.Eta(),lep1_p4.Phi(),
			pt,opt.x,phi,np.log(defineHiggsness()),np.log(opt.fun),0])

        df = pd.DataFrame(data, columns=
                                ['jet1_mass','jet1_pt','jet1_eta','jet1_phi',
                                 'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
                                 'jet3_mass','jet3_pt','jet3_eta','jet3_phi',
                                 'jet4_mass','jet4_pt','jet4_eta','jet4_phi',
				 'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
				 'nu_pt','nu_eta','nu_phi','higgsness','topness','target'])


        df.to_csv("htness_tt.csv", header=False, index=False)

        print("Have processed %i files..."%(ii))
        print("Have processed total %i events..."%(ee))
        print("Event number finally passed is %i..."%(hh))	

	input("Press Enter to continue...")		


if __name__ == "__main__":
	main()
