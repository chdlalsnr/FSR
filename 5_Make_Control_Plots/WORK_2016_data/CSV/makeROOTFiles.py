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
		if "run" not in dr: continue
		f = open('%s.csv'%(dr),'r')
		readline = csv.reader(f)

		h_J1_mass = TH1D("J1_mass", "Mass Distribution of Leading Jet", 20, 0, 200)
		h_J1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_J1_pt   = TH1D("J1_pt", "Transverse Momentum Distribution of Leading Jet", 20, 0, 200)
		h_J1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_J1_eta  = TH1D("J1_eta", "Eta Distribution of Leading Jet", 20, -2.4, 2.4)
		h_J1_eta.GetXaxis().SetTitle("Eta")
		h_J1_phi  = TH1D("J1_phi", "Phi Distribution of Leading Jet", 20, -3.14, 3.14)
		h_J1_phi.GetXaxis().SetTitle("Phi")
		h_J1_btag  = TH1D("J1_btag", "btagDeepFlavB Distribution of Leading Jet", 20, 0, 1.01)
		h_J1_btag.GetXaxis().SetTitle("btagDeepFlavB")

		h_J2_mass = TH1D("J2_mass", "Mass Distribution of Sub-leading Jet", 20, 0, 200)
		h_J2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_J2_pt   = TH1D("J2_pt", "Transverse Momentum Distribution of Sub-leading Jet", 20, 0, 200)
		h_J2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_J2_eta  = TH1D("J2_eta", "Eta Distribution of Sub-leading Jet", 20, -2.4, 2.4)
		h_J2_eta.GetXaxis().SetTitle("Eta")
		h_J2_phi  = TH1D("J2_phi", "Phi Distribution of Sub-leading Jet", 20, -3.14, 3.14)
		h_J2_phi.GetXaxis().SetTitle("Phi")
		h_J2_btag  = TH1D("J2_btag", "btagDeepFlavB Distribution of Sub-leading Jet", 20, 0, 1.01)
		h_J2_btag.GetXaxis().SetTitle("btagDeepFlavB")

		h_L1_mass = TH1D("L1_mass", "Mass Distribution of Leading Lepton", 20, 0, 200)
		h_L1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_L1_pt   = TH1D("L1_pt", "Transverse Momentum Distribution of Leading Lepton", 20, 0, 200)
		h_L1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_L1_eta  = TH1D("L1_eta", "Eta Distribution of Leading Lepton", 20, -2.4, 2.4)
		h_L1_eta.GetXaxis().SetTitle("Eta")
		h_L1_phi  = TH1D("L1_phi", "Phi Distribution of Leading Lepton", 20, -3.14, 3.14)
		h_L1_phi.GetXaxis().SetTitle("Phi")
		h_L1_charge = TH1D("L1_charge", "Charge Distribution of Leading Lepton", 20, -1, 1.01)
		h_L1_charge.GetXaxis().SetTitle("Charge")

		h_L2_mass = TH1D("L2_mass", "Mass Distribution of Sub-leading Lepton", 20, 0, 200)
		h_L2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_L2_pt   = TH1D("L2_pt", "Transverse Momentum Distribution of Sub-leading Lepton", 20, 0, 200)
		h_L2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_L2_eta  = TH1D("L2_eta", "Eta Distribution of Sub-leading Lepton", 20, -2.4, 2.4)
		h_L2_eta.GetXaxis().SetTitle("Eta")
		h_L2_phi  = TH1D("L2_phi", "Phi Distribution of Sub-leading Lepton", 20, -3.14, 3.14)
		h_L2_phi.GetXaxis().SetTitle("Phi")
		h_L2_charge = TH1D("L2_charge", "Charge Distribution of Sub-leading Lepton", 20, -1, 1.01)
		h_L2_charge.GetXaxis().SetTitle("Charge")

		h_dr_ll = TH1D("DeltaR_ll", "DeltaR Distribution between Selected Leptons", 20, 0, 10)
		h_dr_ll.GetXaxis().SetTitle("DeltaR")
		h_dr_jj = TH1D("DeltaR_jj", "DeltaR Distribution between Selected Jets", 20, 0, 10)
		h_dr_jj.GetXaxis().SetTitle("DeltaR")

		h_met_pt = TH1D("MET_pt", "Transverse Momentum Distribution of MET", 20, 0, 200)
		h_met_pt.GetXaxis().SetTitle("Transverse Momentum [GeV]")
		h_met_phi = TH1D("MET_phi", "Phi Distribution of MET", 20, -3.14, 3.14)

		h_H1_mass = TH1D("H1_mass", "Mass Distribution of H1 <- MuMu", 20, 0, 400)
		h_H1_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_H1_pt   = TH1D("H1_pt", "Transverse Momentum Distribution of H1 <- MuMu", 20, 0, 400)
		h_H1_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_H1_eta  = TH1D("H1_eta", "Eta Distribution of H1 <- MuMu", 20, -2.4, 2.4)
		h_H1_eta.GetXaxis().SetTitle("Eta")
		h_H1_phi  = TH1D("H1_phi", "Phi Distribution of H1 <- MuMu", 20, -3.14, 3.14)
		h_H1_phi.GetXaxis().SetTitle("Phi")

		h_H2_mass = TH1D("H2_mass", "Mass Distribution of H2 <- bb", 20, 0, 400)
		h_H2_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_H2_pt   = TH1D("H2_pt", "Transverse Momentum Distribution of H2 <- bb", 20, 0, 400)
		h_H2_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_H2_eta  = TH1D("H2_eta", "Eta Distribution of H2 <- bb", 20, -2.4, 2.4)
		h_H2_eta.GetXaxis().SetTitle("Eta")
		h_H2_phi  = TH1D("H2_phi", "Phi Distribution of H2 <- bb", 20, -3.14, 3.14)
		h_H2_phi.GetXaxis().SetTitle("Phi")

		h_HH_mass = TH1D("HH_mass", "Mass Distribution of HH", 20, 0, 800)
		h_HH_mass.GetXaxis().SetTitle("Invariant mass [GeV]")
		h_HH_pt   = TH1D("HH_pt", "Transverse Momentum Distribution of HH", 20, 0, 800)
		h_HH_pt.GetXaxis().SetTitle("Transverse momentum [GeV]")
		h_HH_eta  = TH1D("HH_eta", "Eta Distribution of HH", 20, -2.4, 2.4)
		h_HH_eta.GetXaxis().SetTitle("Eta")
		h_HH_phi  = TH1D("HH_phi", "Phi Distribution of HH", 20, -3.14, 3.14)
		h_HH_phi.GetXaxis().SetTitle("Phi")

		jet1_mass,jet1_pt,jet1_eta,jet1_phi,jet1_btag = 0,0,0,0,0
		jet2_mass,jet2_pt,jet2_eta,jet2_phi,jet2_btag = 0,0,0,0,0
		lep1_mass,lep1_pt,lep1_eta,lep1_phi,lep1_charge = 0,0,0,0,0
		lep2_mass,lep2_pt,lep2_eta,lep2_phi,lep2_charge = 0,0,0,0,0
		dr_ll,dr_jj,met_pt,met_phi  = 0,0,0,0
		H1_mass,H1_pt,H1_eta,H1_phi = 0,0,0,0
		H2_mass,H2_pt,H2_eta,H2_phi = 0,0,0,0
		HH_mass,HH_pt,HH_eta,HH_phi = 0,0,0,0

		ll = 0
		for line in readline:

			ll += 1
			jet1_mass = float(line[0])
			jet1_pt = float(line[1])
			jet1_eta = float(line[2])
			jet1_phi = float(line[3])
			jet1_btag = float(line[4])

			jet2_mass = float(line[5])
			jet2_pt = float(line[6])
			jet2_eta = float(line[7])
			jet2_phi = float(line[8])
			jet2_btag = float(line[9])

			lep1_mass = float(line[10])
			lep1_pt = float(line[11])
			lep1_eta = float(line[12])
			lep1_phi = float(line[13])
			lep1_charge = float(line[14])

			lep2_mass = float(line[15])
			lep2_pt = float(line[16])
			lep2_eta = float(line[17])
			lep2_phi = float(line[18])
			lep2_charge = float(line[19])

			dr_jj = float(line[20])
			dr_ll = float(line[21])
			met_pt = float(line[22])
			met_phi = float(line[23])

			H1_mass = float(line[24])
			H1_pt = float(line[25])
			H1_eta = float(line[26])
			H1_phi = float(line[27])

			H2_mass = float(line[28])
			H2_pt = float(line[29])
			H2_eta = float(line[30])
			H2_phi = float(line[31])

			HH_mass = float(line[32])
			HH_pt = float(line[33])
			HH_eta = float(line[34])
			HH_phi = float(line[35])

			if (jet1_mass >= h_J1_mass.GetBinCenter(20)): jet1_mass = h_J1_mass.GetBinCenter(20)
			h_J1_mass.Fill(jet1_mass)
			if (jet1_pt >= h_J1_pt.GetBinCenter(20)): jet1_pt = h_J1_pt.GetBinCenter(20)
			h_J1_pt.Fill(jet1_pt)
			h_J1_eta.Fill(jet1_eta)
			h_J1_phi.Fill(jet1_phi)
			h_J1_btag.Fill(jet1_btag)

			if (jet2_mass >= h_J2_mass.GetBinCenter(20)): jet2_mass = h_J2_mass.GetBinCenter(20)
                        h_J2_mass.Fill(jet2_mass)
                        if (jet2_pt >= h_J2_pt.GetBinCenter(20)): jet2_pt = h_J2_pt.GetBinCenter(20)
                        h_J2_pt.Fill(jet2_pt)
                        h_J2_eta.Fill(jet2_eta)
                        h_J2_phi.Fill(jet2_phi)
                        h_J2_btag.Fill(jet2_btag)

			if (lep1_mass >= h_L1_mass.GetBinCenter(20)): lep1_mass = h_L1_mass.GetBinCenter(20)
                        h_L1_mass.Fill(lep1_mass)
                        if (lep1_pt >= h_L1_pt.GetBinCenter(20)): lep1_pt = h_L1_pt.GetBinCenter(20)
                        h_L1_pt.Fill(lep1_pt)
                        h_L1_eta.Fill(lep1_eta)
                        h_L1_phi.Fill(lep1_phi)
                        h_L1_charge.Fill(lep1_charge)

			if (lep2_mass >= h_L2_mass.GetBinCenter(20)): lep2_mass = h_L2_mass.GetBinCenter(20)
                        h_L2_mass.Fill(lep2_mass)
                        if (lep2_pt >= h_L2_pt.GetBinCenter(20)): lep2_pt = h_L2_pt.GetBinCenter(20)
                        h_L2_pt.Fill(lep2_pt)
                        h_L2_eta.Fill(lep2_eta)
                        h_L2_phi.Fill(lep2_phi)
                        h_L2_charge.Fill(lep2_charge)

                        if (dr_ll >= h_dr_ll.GetBinCenter(20)): dr_ll = h_dr_ll.GetBinCenter(20)
			h_dr_ll.Fill(dr_ll)
                        if (dr_jj >= h_dr_jj.GetBinCenter(20)): dr_jj = h_dr_jj.GetBinCenter(20)
			h_dr_jj.Fill(dr_jj)
			if (met_pt >= h_met_pt.GetBinCenter(20)): met_pt = h_met_pt.GetBinCenter(20)
			h_met_pt.Fill(met_pt)
			if (met_phi >= h_met_phi.GetBinCenter(20)): met_phi = h_met_phi.GetBinCenter(20)
			h_met_phi.Fill(met_phi)

			if (H1_mass >= h_H1_mass.GetBinCenter(20)): H1_mass = h_H1_mass.GetBinCenter(20)
                        h_H1_mass.Fill(H1_mass)
                        if (H1_pt >= h_H1_pt.GetBinCenter(20)): H1_pt = h_H1_pt.GetBinCenter(20)
                        h_H1_pt.Fill(H1_pt)
                        h_H1_eta.Fill(H1_eta)
                        h_H1_phi.Fill(H1_phi)

			if (H2_mass >= h_H2_mass.GetBinCenter(20)): H2_mass = h_H2_mass.GetBinCenter(20)
                        h_H2_mass.Fill(H2_mass)
                        if (H2_pt >= h_H2_pt.GetBinCenter(20)): H2_pt = h_H2_pt.GetBinCenter(20)
                        h_H2_pt.Fill(H2_pt)
                        h_H2_eta.Fill(H2_eta)
                        h_H2_phi.Fill(H2_phi)

			if (HH_mass >= h_HH_mass.GetBinCenter(20)): HH_mass = h_HH_mass.GetBinCenter(20)
                        h_HH_mass.Fill(HH_mass)
                        if (HH_pt >= h_HH_pt.GetBinCenter(20)): HH_pt = h_HH_pt.GetBinCenter(20)
                        h_HH_pt.Fill(HH_pt)
                        h_HH_eta.Fill(HH_eta)
                        h_HH_phi.Fill(HH_phi)

		f.close()

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
		h_met_pt.Write()
		h_met_phi.Write()

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