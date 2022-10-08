#python
# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np 
import math
import csv 
import sys
import ast

def main():


	flist = glob('*.csv')
	for fname in flist:

		dr = fname[:-4]
		f = open('%s.csv'%(dr),'r')
		readline = csv.reader(f)
		print("Processing %s..."%(dr))

		jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag = [],[],[],[],[]
		jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag = [],[],[],[],[]
		lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge = [],[],[],[],[]
		lep2_mass,lep2_pt,lep2_eta,lep2_phi,lep2_charge = [],[],[],[],[]
		H1_mass,H1_pt,H1_eta,H1_phi = [],[],[],[]
		H2_mass,H2_pt,H2_eta,H2_phi = [],[],[],[]
		HH_mass,HH_pt,HH_eta,HH_phi = [],[],[],[]
		dr_ll,dr_jj = [],[]

		ll = 0
		for line in readline:

			jet1_mass.append(eval(line[0]))
			jet1_pt.append(eval(line[1]))
			jet1_eta.append(eval(line[2]))
			jet1_phi.append(eval(line[3]))
			jet1_btag.append(eval(line[4]))
			jet2_mass.append(eval(line[5]))
			jet2_pt.append(eval(line[6]))
			jet2_eta.append(eval(line[7]))
			jet2_phi.append(eval(line[8]))
			jet2_btag.append(eval(line[9]))
			lep1_mass.append(eval(line[10]))
			lep1_pt.append(eval(line[11]))
			lep1_eta.append(eval(line[12]))
			lep1_phi.append(eval(line[13]))
			lep1_charge.append(eval(line[14]))
			lep2_mass.append(eval(line[15]))
			lep2_pt.append(eval(line[16]))
			lep2_eta.append(eval(line[17]))
			lep2_phi.append(eval(line[18]))
			lep2_charge.append(eval(line[19]))
			dr_jj.append(eval(line[20]))
			dr_ll.append(eval(line[21]))
			H1_mass.append(eval(line[22]))
			H1_pt.append(eval(line[23]))
			H1_eta.append(eval(line[24]))
			H1_phi.append(eval(line[25]))
			H2_mass.append(eval(line[26]))
			H2_pt.append(eval(line[27]))
			H2_eta.append(eval(line[28]))
			H2_phi.append(eval(line[29]))
			HH_mass.append(eval(line[30]))
			HH_pt.append(eval(line[31]))
			HH_eta.append(eval(line[32]))
			HH_phi.append(eval(line[33]))
			ll = ll+1

		f.close()

		h_J1_mass = TH1D("J1_mass", "Mass Distribution of Leading Jet", 40, 0, 200)
		h_J1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_J1_pt   = TH1D("J1_pt", "Transverse Momentum Distribution of Leading Jet", 40, 0, 400)
		h_J1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_J1_eta  = TH1D("J1_eta", "Eta Distribution of Leading Jet", 40, -3.4, 3.4)
		h_J1_eta.GetXaxis().SetTitle("Eta")
		h_J1_phi  = TH1D("J1_phi", "Phi Distribution of Leading Jet", 40, -3.14, 3.14)
		h_J1_phi.GetXaxis().SetTitle("Phi")
		h_J1_btag  = TH1D("J1_btag", "btagDeepFlavB Distribution of Leading Jet", 40, 0, 1.01)
		h_J1_btag.GetXaxis().SetTitle("btagDeepFlavB")

		h_J2_mass = TH1D("J2_mass", "Mass Distribution of Sub-leading Jet", 40, 0, 200)
		h_J2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_J2_pt   = TH1D("J2_pt", "Transverse Momentum Distribution of Sub-leading Jet", 40, 0, 400)
		h_J2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_J2_eta  = TH1D("J2_eta", "Eta Distribution of Sub-leading Jet", 40, -3.4, 3.4)
		h_J2_eta.GetXaxis().SetTitle("Eta")
		h_J2_phi  = TH1D("J2_phi", "Phi Distribution of Sub-leading Jet", 40, -3.14, 3.14)
		h_J2_phi.GetXaxis().SetTitle("Phi")
		h_J2_btag  = TH1D("J2_btag", "btagDeepFlavB Distribution of Sub-leading Jet", 40, 0, 1.01)
		h_J2_btag.GetXaxis().SetTitle("btagDeepFlavB")

		h_L1_mass = TH1D("L1_mass", "Mass Distribution of Leading Lepton", 40, 0, 200)
		h_L1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_L1_pt   = TH1D("L1_pt", "Transverse Momentum Distribution of Leading Lepton", 40, 0, 400)
		h_L1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_L1_eta  = TH1D("L1_eta", "Eta Distribution of Leading Lepton", 40, -3.4, 3.4)
		h_L1_eta.GetXaxis().SetTitle("Eta")
		h_L1_phi  = TH1D("L1_phi", "Phi Distribution of Leading Lepton", 40, -3.14, 3.14)
		h_L1_phi.GetXaxis().SetTitle("Phi")
		h_L1_charge = TH1D("L1_charge", "Charge Distribution of Leading Lepton", 40, -1, 1.01)
		h_L1_charge.GetXaxis().SetTitle("Charge")

		h_L2_mass = TH1D("L2_mass", "Mass Distribution of Sub-leading Lepton", 40, 0, 200)
		h_L2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_L2_pt   = TH1D("L2_pt", "Transverse Momentum Distribution of Sub-leading Lepton", 40, 0, 400)
		h_L2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_L2_eta  = TH1D("L2_eta", "Eta Distribution of Sub-leading Lepton", 40, -3.4, 3.4)
		h_L2_eta.GetXaxis().SetTitle("Eta")
		h_L2_phi  = TH1D("L2_phi", "Phi Distribution of Sub-leading Lepton", 40, -3.14, 3.14)
		h_L2_phi.GetXaxis().SetTitle("Phi")
		h_L2_charge = TH1D("L2_charge", "Charge Distribution of Sub-leading Lepton", 40, -1, 1.01)
		h_L2_charge.GetXaxis().SetTitle("Charge")

		h_H1_mass = TH1D("H1_mass", "Mass Distribution of H1 <- MuMu", 40, 0, 200)
		h_H1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_H1_pt   = TH1D("H1_pt", "Transverse Momentum Distribution of H1 <- MuMu", 40, 0, 400)
		h_H1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_H1_eta  = TH1D("H1_eta", "Eta Distribution of H1 <- MuMu", 40, -3.4, 3.4)
		h_H1_eta.GetXaxis().SetTitle("Eta")
		h_H1_phi  = TH1D("H1_phi", "Phi Distribution of H1 <- MuMu", 40, -3.14, 3.14)
		h_H1_phi.GetXaxis().SetTitle("Phi")

		h_H2_mass = TH1D("H2_mass", "Mass Distribution of H2 <- bb", 40, 0, 200)
		h_H2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_H2_pt   = TH1D("H2_pt", "Transverse Momentum Distribution of H2 <- bb", 40, 0, 400)
		h_H2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_H2_eta  = TH1D("H2_eta", "Eta Distribution of H2 <- bb", 40, -3.4, 3.4)
		h_H2_eta.GetXaxis().SetTitle("Eta")
		h_H2_phi  = TH1D("H2_phi", "Phi Distribution of H2 <- bb", 40, -3.14, 3.14)
		h_H2_phi.GetXaxis().SetTitle("Phi")

		h_HH_mass = TH1D("HH_mass", "Mass Distribution of HH", 40, 0, 600)
		h_HH_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_HH_pt   = TH1D("HH_pt", "Transverse Momentum Distribution of HH", 40, 0, 600)
		h_HH_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_HH_eta  = TH1D("HH_eta", "Eta Distribution of HH", 40, -3.4, 3.4)
		h_HH_eta.GetXaxis().SetTitle("Eta")
		h_HH_phi  = TH1D("HH_phi", "Phi Distribution of HH", 40, -3.14, 3.14)
		h_HH_phi.GetXaxis().SetTitle("Phi")

		h_dr_ll = TH1D("DeltaR_ll", "DeltaR Distribution between Selected Leptons", 40, 0, 10)
		h_dr_ll.GetXaxis().SetTitle("DeltaR")
		h_dr_jj = TH1D("DeltaR_jj", "DeltaR Distribution between Selected Jets", 40, 0, 10)
		h_dr_jj.GetXaxis().SetTitle("DeltaR")

		for i in range(ll):

			if (jet1_mass[i] >= h_J1_mass.GetBinCenter(40)): jet1_mass[i] = h_J1_mass.GetBinCenter(40)
			h_J1_mass.Fill(jet1_mass[i])
			if (jet1_pt[i] >= h_J1_pt.GetBinCenter(40)): jet1_pt[i] = h_J1_pt.GetBinCenter(40)
			h_J1_pt.Fill(jet1_pt[i])
			h_J1_eta.Fill(jet1_eta[i])
			h_J1_phi.Fill(jet1_phi[i])
			h_J1_btag.Fill(jet1_btag[i])

			if (jet2_mass[i] >= h_J2_mass.GetBinCenter(40)): jet2_mass[i] = h_J2_mass.GetBinCenter(40)
                        h_J2_mass.Fill(jet2_mass[i])
                        if (jet2_pt[i] >= h_J2_pt.GetBinCenter(40)): jet2_pt[i] = h_J2_pt.GetBinCenter(40)
                        h_J2_pt.Fill(jet2_pt[i])
                        h_J2_eta.Fill(jet2_eta[i])
                        h_J2_phi.Fill(jet2_phi[i])
                        h_J2_btag.Fill(jet2_btag[i])

			if (lep1_mass[i] >= h_L1_mass.GetBinCenter(40)): lep1_mass[i] = h_L1_mass.GetBinCenter(40)
                        h_L1_mass.Fill(lep1_mass[i])
                        if (lep1_pt[i] >= h_L1_pt.GetBinCenter(40)): lep1_pt[i] = h_L1_pt.GetBinCenter(40)
                        h_L1_pt.Fill(lep1_pt[i])
                        h_L1_eta.Fill(lep1_eta[i])
                        h_L1_phi.Fill(lep1_phi[i])
                        h_L1_charge.Fill(lep1_charge[i])

			if (lep2_mass[i] >= h_L2_mass.GetBinCenter(40)): lep2_mass[i] = h_L2_mass.GetBinCenter(40)
                        h_L2_mass.Fill(lep2_mass[i])
                        if (lep2_pt[i] >= h_L2_pt.GetBinCenter(40)): lep2_pt[i] = h_L2_pt.GetBinCenter(40)
                        h_L2_pt.Fill(lep2_pt[i])
                        h_L2_eta.Fill(lep2_eta[i])
                        h_L2_phi.Fill(lep2_phi[i])
                        h_L2_charge.Fill(lep2_charge[i])

                        if (dr_ll[i] >= h_dr_ll.GetBinCenter(40)): dr_ll[i] = h_dr_ll.GetBinCenter(40)
			h_dr_ll.Fill(dr_ll[i])
                        if (dr_jj[i] >= h_dr_jj.GetBinCenter(40)): dr_jj[i] = h_dr_jj.GetBinCenter(40)
			h_dr_jj.Fill(dr_jj[i])

			if (H1_mass[i] >= h_H1_mass.GetBinCenter(40)): H1_mass[i] = h_H1_mass.GetBinCenter(40)
                        h_H1_mass.Fill(H1_mass[i])
                        if (H1_pt[i] >= h_H1_pt.GetBinCenter(40)): H1_pt[i] = h_H1_pt.GetBinCenter(40)
                        h_H1_pt.Fill(H1_pt[i])
                        h_H1_eta.Fill(H1_eta[i])
                        h_H1_phi.Fill(H1_phi[i])

			if (H2_mass[i] >= h_H2_mass.GetBinCenter(40)): H2_mass[i] = h_H2_mass.GetBinCenter(40)
                        h_H2_mass.Fill(H2_mass[i])
                        if (H2_pt[i] >= h_H2_pt.GetBinCenter(40)): H2_pt[i] = h_H2_pt.GetBinCenter(40)
                        h_H2_pt.Fill(H2_pt[i])
                        h_H2_eta.Fill(H2_eta[i])
                        h_H2_phi.Fill(H2_phi[i])

			if (HH_mass[i] >= h_HH_mass.GetBinCenter(40)): HH_mass[i] = h_HH_mass.GetBinCenter(40)
                        h_HH_mass.Fill(HH_mass[i])
                        if (HH_pt[i] >= h_HH_pt.GetBinCenter(40)): HH_pt[i] = h_HH_pt.GetBinCenter(40)
                        h_HH_pt.Fill(HH_pt[i])
                        h_HH_eta.Fill(HH_eta[i])
                        h_HH_phi.Fill(HH_phi[i])


		ff = TFile("./%s.root"%(dr),"recreate")
		ff.cd()
		gDirectory.mkdir(dr)
		ff.cd(dr)
		
		h_J1_mass.Write()
		h_J1_pt.Write()
		h_J1_eta.Write()
		h_J1_phi.Write()
		h_J1_btag.Write()

		h_J2_mass.Write()
		h_J2_pt.Write()
		h_J2_eta.Write()
		h_J2_phi.Write()
		h_J2_btag.Write()

		h_L1_mass.Write()
		h_L1_pt.Write()
		h_L1_eta.Write()
		h_L1_phi.Write()
		h_L1_charge.Write()

		h_L2_mass.Write()
		h_L2_pt.Write()
		h_L2_eta.Write()
		h_L2_phi.Write()
		h_L2_charge.Write()

		h_dr_ll.Write()
		h_dr_jj.Write()

		h_H1_mass.Write()
		h_H1_pt.Write()
		h_H1_eta.Write()
		h_H1_phi.Write()

		h_H2_mass.Write()
		h_H2_pt.Write()
		h_H2_eta.Write()
		h_H2_phi.Write()

		h_HH_mass.Write()
		h_HH_pt.Write()
		h_HH_eta.Write()
		h_HH_phi.Write()

		ff.Close()
		print("Processed %d events..."%(ll))
		print("Successfully generated root file...")


if __name__ == "__main__":
        main()
