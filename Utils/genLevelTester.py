from ROOT import *
import numpy as np
#from array import array

f = TFile("F6B8F234-CC43-5D44-8D1E-727E33C38AA7.root", "read") #nanoAOD
#f = TFile("36B3D9CB-E060-E711-B73A-0026B92785F6.root", "read") #miniAOD
t = f.Get("Events")

jet_hf, jet_pf, jet_eta = [],[],[]
jet_mass, jet_phi, jet_pt = [],[],[]
genjet_hf, genjet_pf, genjet_eta = [],[],[]
genjet_mass, genjet_phi, genjet_pt = [],[],[]

#evt_num = t.GetEntriesFast()
evt_num = 1000

for i in range(evt_num):
        t.GetEntry(i)

	temp_breg = []
	temp_hf, temp_pf, temp_eta = [],[],[]
	temp_mass, temp_phi, temp_pt = [],[],[]
	temp_hf_gj, temp_pf_gj, temp_eta_gj = [],[],[]
        temp_mass_gj, temp_phi_gj, temp_pt_gj = [],[],[]

	for j in range(t.nJet): #(t.Jet_mass.__len__())
                temp_hf.append(t.Jet_hadronFlavour.__getitem__(j))
		temp_pf.append(t.Jet_partonFlavour.__getitem__(j))
		temp_eta.append(t.Jet_eta.__getitem__(j))
		temp_mass.append(t.Jet_mass.__getitem__(j))
		temp_phi.append(t.Jet_phi.__getitem__(j))
		temp_pt.append(t.Jet_pt.__getitem__(j))
		temp_breg.append(t.Jet_bRegCorr.__getitem__(j))

		temp_mass[j] *= temp_breg[j]
		temp_pt[j] *= temp_breg[j]

	for j in range(t.GenJet_mass.__len__()):
		temp_hf_gj.append(t.GenJet_hadronFlavour.__getitem__(j))
                temp_pf_gj.append(t.GenJet_partonFlavour.__getitem__(j))
                temp_eta_gj.append(t.GenJet_eta.__getitem__(j))
                temp_mass_gj.append(t.GenJet_mass.__getitem__(j))
                temp_phi_gj.append(t.GenJet_phi.__getitem__(j))
                temp_pt_gj.append(t.GenJet_pt.__getitem__(j))

	jet_hf.append(temp_hf)
	jet_pf.append(temp_pf)
	jet_eta.append(temp_eta)
	jet_mass.append(temp_mass)
	jet_phi.append(temp_phi)
	jet_pt.append(temp_pt)

	genjet_hf.append(temp_hf_gj)
	genjet_pf.append(temp_pf_gj)
	genjet_eta.append(temp_eta_gj)
	genjet_mass.append(temp_mass_gj)
	genjet_phi.append(temp_phi_gj)
	genjet_pt.append(temp_pt_gj)

print(genjet_hf)

h1 = TH1D("h1","jet_hf_corr",10,0,1)
#h2 = TH1D("h2","jet_pf_corr",10,0,1)
h3 = TH1D("h3","jet_eta",80,-4,4)
h4 = TH1D("h4","jet_mass",50,0,50)
h5 = TH1D("h5","jet_phi",63,-3.15,3.15)
h6 = TH1D("h6","jet_pt",40,0,400)

g1 = TH1D("g1","genjet_hf",10,0,1)
#g2 = TH1D("g2","genjet_pf",10,0,1)
g3 = TH1D("g3","genjet_eta",80,-4,4)
g4 = TH1D("g4","genjet_mass",50,0,50)
g5 = TH1D("g5","genjet_phi",63,-3.15,3.15)
g6 = TH1D("g6","genjet_pt",40,0,400)

for i in range(evt_num):
	t.GetEntry(i)
	for j in range(t.nJet):
		h3.Fill(jet_eta[i][j])
		h4.Fill(jet_mass[i][j])
		h5.Fill(jet_phi[i][j])
		h6.Fill(jet_pt[i][j])
		if (jet_hf[i][j] == 5): h1.Fill(1)
		else: h1.Fill(0)
		#if (jet_pf[i][j] == 'b'): h2.Fill(1) ...dunno
		#else: h2.Fill(0)

        for j in range(t.GenJet_mass.__len__()):
                g3.Fill(genjet_eta[i][j])
                g4.Fill(genjet_mass[i][j])
                g5.Fill(genjet_phi[i][j])
                g6.Fill(genjet_pt[i][j])		
		if (genjet_hf[i][j] == 5): g1.Fill(1)
		else: g1.Fill(0)
		#if (genjet_pf[i][j] == 'b'): g2.Fill(1)
		#else: g2.Fill(0)

a = TCanvas()
a.Divide(2,1)
a.cd(1)
h1.Draw()
a.cd(2)
g1.Draw("same")
a.SaveAs("hadflv.png")

#b = TCanvas()
#h2.Draw()
#b.SaveAs("parflv.png")
"""
c = TCanvas()
c.Divide(2,1)
c.cd(1)
h3.Draw()
c.cd(2)
g3.Draw("same")
c.SaveAs("eta.png")

d = TCanvas()
d.Divide(2,1)
d.cd(1)
h4.Draw()
d.cd(2)
g4.Draw("same")
d.SaveAs("mass.png")

e = TCanvas()
e.Divide(2,1)
e.cd(1)
h5.Draw()
e.cd(2)
g5.Draw("same")
e.SaveAs("phi.png")

g = TCanvas()
g.Divide(2,1)
g.cd(1)
h6.Draw()
g.cd(2)
g6.Draw("same")
g.SaveAs("pt.png")
"""
input("Press Enter to continue...")
