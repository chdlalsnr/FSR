from ROOT import *
import numpy as np
#from array import array

f = TFile("F6B8F234-CC43-5D44-8D1E-727E33C38AA7.root", "read") #nanoAOD
#f = TFile("36B3D9CB-E060-E711-B73A-0026B92785F6.root", "read") #miniAOD
t = f.Get("Events")

evt_num = 100
#evt_num = t.GetEntriesFast()
filter_nJet = 0
filter_bScore = 0
filtered = filter_nJet #+ filter_bScore
filtered_idx = []
jet_btag,jet_pt,jet_eta,jet_phi,jet_mass = [],[],[],[],[]

for i in range(evt_num):
	t.GetEntry(i)
	temp_btag,temp_pt,temp_eta,temp_phi,temp_mass = [],[],[],[],[]
	temp_breg = []

	if (t.nJet < 2):
		filter_nJet += 1
		filtered_idx.append(i)
		jet_btag.append([])
                jet_mass.append([])
                jet_pt.append([])
                jet_eta.append([])
                jet_phi.append([])
		continue

	btag_n = 0
	for j in range(t.nJet): #(t.Jet_mass.__len__())
		temp_btag.append(t.Jet_btagDeepFlavB.__getitem__(j))
		if (t.Jet_btagDeepFlavB.__getitem__(j) > 0.8): btag_n += 1
		temp_mass.append(t.Jet_mass.__getitem__(j))
		temp_pt.append(t.Jet_pt.__getitem__(j))
		temp_eta.append(t.Jet_pt.__getitem__(j))
		temp_phi.append(t.Jet_pt.__getitem__(j))
		temp_breg.append(t.Jet_bRegCorr.__getitem__(j))

		temp_mass[j] *= temp_breg[j]
		temp_pt[j] *= temp_breg[j]
	temp_btag.sort(reverse=True)

	jet_btag.append(temp_btag)
        jet_mass.append(temp_mass)
        jet_pt.append(temp_pt)
        jet_eta.append(temp_eta)
        jet_phi.append(temp_phi)

H1_m = []
H1_p = []
H1_e = []
for i in range(evt_num):
        t.GetEntry(i)

	jet1_p4 = TLorentzVector()
	jet2_p4 = TLorentzVector()
	if i not in filtered_idx:
		jet1_p4.SetPtEtaPhiM(jet_pt[i][0],jet_eta[i][0],jet_phi[i][0],jet_mass[i][0])
        	jet2_p4.SetPtEtaPhiM(jet_pt[i][1],jet_eta[i][1],jet_phi[i][1],jet_mass[i][1])
        jet_p4 = jet1_p4 + jet2_p4
        H1_m.append(jet_p4.M())
	H1_p.append(jet_p4.P())
	H1_e.append(jet_p4.E())

#print(H1_m)
#print("==============")
#print(H1_p)
#print("==============")
#print(H1_e)

h_H1_mass = TH1D("H1_mass_test", "Mass Distribution of H1<-bb", 100, 0, 400)
c = TCanvas()
for i in range(len(H1_m)):
        if i not in filtered_idx:
		h_H1_mass.Fill(H1_m[i])

h_H1_mass.Draw()
c.SaveAs("Higgs1_mass_test.png")
input("Press Enter to continue...")

