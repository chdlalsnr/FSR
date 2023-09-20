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
def main(dataset): #index 

        with open("dataset_2016/%s.dat"%(dataset),'rb') as js:
                flist = pickle.load(js)

        ii,hh,ee = 0,0,0
	genweightSum = 0
        data = []

	h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876
	H1_mass_cand,H1_pt_cand,H1_eta_cand,H1_phi_cand = [],[],[],[]
        H2_mass_cand,H2_pt_cand,H2_eta_cand,H2_phi_cand = [],[],[],[]
        HH_mass_cand,HH_pt_cand,HH_eta_cand,HH_phi_cand = [],[],[],[]
        jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag = [],[],[],[],[]
        jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag = [],[],[],[],[]
        lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge = [],[],[],[],[]
        lep2_mass,lep2_pt,lep2_eta,lep2_phi,lep2_charge = [],[],[],[],[]
        dr_ll,dr_jj,genWeight = [],[],[]
        c_dilep,c_dijet,c_met = [],[],[]
        for fname in flist:

		f = TFile.Open(fname,"read")
		t = f.Get("Events")
		r = f.Get("Runs")

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

			met_filter = False
			if (flag_good == 1 and flag_halo == 1 and flag_hbhen == 1 and flag_hbiso == 1 and flag_dead == 1 and flag_badpf == 1 and flag_ecbad == 1):
				met_filter = True

			if (hlt_mumu or hlt_mu or hlt_el):
				if met_filter: hlt_pass.append(i)

		# =================================================================
		#Remained cut : The lepton isolation, defined as the scalar
		#p T sum of all particle candidates, excluding the lepton, in a cone around the lepton, divided by
		#the lepton p T , is required to be < 0.04 ( < 0.15) for electrons (muons)
		# Medium muon discrimination
		for i in range(event_num):

			ee += 1
			met_pt,met_phi = 0,0
			t.GetEntry(i)
			genweight = t.genWeight
			genweightSum += genweight
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


			l1,l2 = -1,-1
			for k in range(t.nMuon):
                                if (flag_lep == True): break
				if (temp_pt[k] > 25 and abs(temp_eta[k]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_medid[k] == True and temp_iso[k] < 0.4 and temp_sip3d[k] < 8 and temp_mva[k] > 0.5):
                                        for l in range(k+1,t.nMuon):
						if (flag_lep == True): break
                                                if (temp_pt[l] > 15 and abs(temp_eta[l]) < 2.4 and abs(temp_dxy[k]) < 0.05 and abs(temp_dz[k]) < 0.1 and temp_medid[l] == True and temp_iso[l] < 0.4 and temp_sip3d[l] < 8 and temp_mva[l] > 0.5):

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
                        temp_btag,temp_jetid,temp_charge = [],[],[]
			temp_rho = t.fixedGridRhoFastjetAll
			met_pt  = t.MET_pt
			met_phi = t.MET_phi

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


			# JEC L1L2L3 Correction
			for j in range(len(temp_jet)):

				if temp_rho < 20:

					if (abs(temp_jet[j][2]) < 0.8): temp_jet[j][1] = temp_jet[j][1]
					if (abs(temp_jet[j][2]) > 0.8 and abs(temp_jet[j][2]) < 1.2): temp_jet[j][1] = temp_jet[j][1]*1.05
					if (abs(temp_jet[j][2]) > 1.2 and abs(temp_jet[j][2]) < 1.8): temp_jet[j][1] = temp_jet[j][1]*1.1
					if (abs(temp_jet[j][2]) > 1.8 and abs(temp_jet[j][2]) < 2.15): temp_jet[j][1] = temp_jet[j][1]*1.15
					if (abs(temp_jet[j][2]) > 2.15 and abs(temp_jet[j][2]) < 2.4): temp_jet[j][1] = temp_jet[j][1]*1.1

				if (temp_rho > 20 and temp_rho < 40):

					if temp_jet[j][1] < 30:

						if (abs(temp_jet[j][2]) < 0.8): temp_jet[j][1] = temp_jet[j][1]*0.8
						if (abs(temp_jet[j][2]) > 0.8 and abs(temp_jet[j][2]) < 1.2): temp_jet[j][1] = temp_jet[j][1]*0.85
						if (abs(temp_jet[j][2]) > 1.2 and abs(temp_jet[j][2]) < 1.8): temp_jet[j][1] = temp_jet[j][1]*0.95
						if (abs(temp_jet[j][2]) > 1.8 and abs(temp_jet[j][2]) < 2.15): temp_jet[j][1] = temp_jet[j][1]*1.05
						if (abs(temp_jet[j][2]) > 2.15 and abs(temp_jet[j][2]) < 2.4): temp_jet[j][1] = temp_jet[j][1]*0.95

					if (temp_jet[j][1] > 30 and temp_jet[j][1] < 100):

                                                if (abs(temp_jet[j][2]) < 0.8): temp_jet[j][1] = temp_jet[j][1]*0.95
                                                if (abs(temp_jet[j][2]) > 0.8 and abs(temp_jet[j][2]) < 1.2): temp_jet[j][1] = temp_jet[j][1]
                                                if (abs(temp_jet[j][2]) > 1.2 and abs(temp_jet[j][2]) < 1.8): temp_jet[j][1] = temp_jet[j][1]*1.05
                                                if (abs(temp_jet[j][2]) > 1.8 and abs(temp_jet[j][2]) < 2.15): temp_jet[j][1] = temp_jet[j][1]*1.1
                                                if (abs(temp_jet[j][2]) > 2.15 and abs(temp_jet[j][2]) < 2.4): temp_jet[j][1] = temp_jet[j][1]*1.05

					if (temp_jet[j][1] > 100):

                                                if (abs(temp_jet[j][2]) < 0.8): temp_jet[j][1] = temp_jet[j][1]*1.05
                                                if (abs(temp_jet[j][2]) > 0.8 and abs(temp_jet[j][2]) < 1.2): temp_jet[j][1] = temp_jet[j][1]*1.05
                                                if (abs(temp_jet[j][2]) > 1.2 and abs(temp_jet[j][2]) < 1.8): temp_jet[j][1] = temp_jet[j][1]*1.1
                                                if (abs(temp_jet[j][2]) > 1.8 and abs(temp_jet[j][2]) < 2.15): temp_jet[j][1] = temp_jet[j][1]*1.1
                                                if (abs(temp_jet[j][2]) > 2.15 and abs(temp_jet[j][2]) < 2.4): temp_jet[j][1] = temp_jet[j][1]*1.05

				if temp_rho > 40:

                                        if (abs(temp_jet[j][2]) < 0.8): temp_jet[j][1] = temp_jet[j][1]*0.9
                                        if (abs(temp_jet[j][2]) > 0.8 and abs(temp_jet[j][2]) < 1.2): temp_jet[j][1] = temp_jet[j][1]*0.95
                                        if (abs(temp_jet[j][2]) > 1.2 and abs(temp_jet[j][2]) < 1.8): temp_jet[j][1] = temp_jet[j][1]*1.05
                                        if (abs(temp_jet[j][2]) > 1.8 and abs(temp_jet[j][2]) < 2.15): temp_jet[j][1] = temp_jet[j][1]*1.1
                                        if (abs(temp_jet[j][2]) > 2.15 and abs(temp_jet[j][2]) < 2.4): temp_jet[j][1] = temp_jet[j][1]


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

							if (dijet_p4.M() == 0): continue
                                                        if (dijet_p4.M() > 0): flag_jet = True

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

				met.append(met_pt)
				met.append(met_phi)

				c_dilep.append(dilep)
				c_dijet.append(dijet)
				c_met.append(met)

				genWeight.append(genweight)
		# =================================================================
		# Higgs reconstruction
		for i in range(len(c_dilep)):

			lep1_p4 = TLorentzVector()
			lep2_p4 = TLorentzVector()
			lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
			lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])
			met_p4 = TLorentzVector()
			met_p4.SetPtEtaPhiM(c_met[i][0],0,c_met[i][1],0)
			h1_p4 = lep1_p4 + lep2_p4 + met_p4

			jet1_p4 = TLorentzVector()
			jet2_p4 = TLorentzVector()
			jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][0])
			jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])
			h2_p4 = jet1_p4 + jet2_p4

			hh_p4 = TLorentzVector()
			hh_p4 = h1_p4 + h2_p4

			H1_mass_cand.append(h1_p4.M())
			H1_pt_cand.append(h1_p4.Pt())
			H1_eta_cand.append(h1_p4.Eta())
			H1_phi_cand.append(h1_p4.Phi())

			H2_mass_cand.append(h2_p4.M())
			H2_pt_cand.append(h2_p4.Pt())
			H2_eta_cand.append(h2_p4.Eta())
			H2_phi_cand.append(h2_p4.Phi())

			HH_mass_cand.append(hh_p4.M())
			HH_pt_cand.append(hh_p4.Pt())
			HH_eta_cand.append(hh_p4.Eta())
			HH_phi_cand.append(hh_p4.Phi())

		# =================================================================
		# Save CSV file
		for i in range(len(c_dilep)):
			data.append([jet1_mass[i],jet1_pt[i],jet1_eta[i],jet1_phi[i],jet1_btag[i],
				jet2_mass[i],jet2_pt[i],jet2_eta[i],jet2_phi[i],jet2_btag[i],
				lep1_mass[i],lep1_pt[i],lep1_eta[i],lep1_phi[i],lep1_charge[i],
				lep2_mass[i],lep2_pt[i],lep2_eta[i],lep2_phi[i],lep2_charge[i],
				dr_jj[i],dr_ll[i],
				H1_mass_cand[i],H1_pt_cand[i],H1_eta_cand[i],H1_phi_cand[i],
				H2_mass_cand[i],H2_pt_cand[i],H2_eta_cand[i],H2_phi_cand[i],
				HH_mass_cand[i],HH_pt_cand[i],HH_eta_cand[i],HH_phi_cand[i],genWeight[i]])

	df = pd.DataFrame(data, columns=
				['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet1_btag',
				 'jet2_mass','jet2_pt','jet2_eta','jet2_phi','jet2_btag',
				 'lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep1_charge',
				 'lep2_mass','lep2_pt','lep2_eta','lep2_phi','lep2_charge',
				 'dr_jj','dr_ll',
				 'H1_mass','H1_pt','H1_eta','H1_phi',
				 'H2_mass','H2_pt','H2_eta','H2_phi',
				 'HH_mass','HH_pt','HH_eta','HH_phi','genWeight'])

	df.to_csv("./%s.csv"%(dataset), header=False, index=False)
#	df.to_csv("./%s.csv"%(key), header=True, index=False)
#	xf = open("./%s.txt"%(key),'w')
#	xf.write("%s"%(ee))
#	xf.write("%s"%(hh))
#	xf.close()

	print("Have processed %i files..."%(ii))
	print("Have processed total %i events..."%(ee))
	print("Event number finally passed is %i..."%(hh))
	print("genweightSum is %f..."%(genweightSum))
#	input("Press Enter to continue...")


if __name__ == "__main__":

	main(sys.argv[1])
