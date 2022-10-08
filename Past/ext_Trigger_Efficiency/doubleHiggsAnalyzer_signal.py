#!/usr/bin/env python
import ROOT
from ROOT import *
import numpy as np
import math
import glob
import sys

def main():

	# Signal
        flist = glob.glob('/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/130000/64365B41-751A-7B4B-836B-DB4258390231.root')

	ii,hh,ee = 0,0,0
	jet_r,jet_g = [],[]
	for fname in flist:

		if (ee > 5000000): break
		f = TFile(fname,"read")
		t = f.Get("Events")

		#event_num = 10000
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
		for i in range(event_num):

			ee = ee+1
			t.GetEntry(i)
			if i not in hlt_pass: continue

			for j in range(t.nJet):
				jet_r.append(t.Jet_pt.__getitem__(j))
			for j in range(t.GenJet_pt.__len__()):
				jet_g.append(t.GenJet_pt.__getitem__(j))


	reco_jet = TH1D("reco_jet","Reconstructed Jets",20,0,200)
	for i in range(len(jet_r)): reco_jet.Fill(jet_r[i])
	gen_jet  = TH1D("gen_jet","Generated Jets",20,0,200)
	for i in range(len(jet_g)): gen_jet.Fill(jet_g[i])

	trig_eff = TH1D("trig_eff","Trigger Efficiency Jet",20,0,200)
	trig_eff.Divide(gen_jet,reco_jet)

	cvs = TCanvas("cvs","Trigger Efficiency Jet",900,600)
	cvs.cd()

	trig_eff.SetMarkerColor(kBlue)
        trig_eff.SetMarkerSize(1.5)
        trig_eff.SetMarkerStyle(21)
	trig_eff.GetXaxis().SetRangeUser(20,200)
	trig_eff.GetYaxis().SetRangeUser(0,1.1)

	trig_eff.Draw("P")

	cvs.SaveAs("Trigger_Efficiency_(Jet).png")

	print("Have processed %i files..."%(ii))
	print("Have processed total %i events..."%(ee))
	print("Event number passed hlt is %i..."%(hh))
	input("Press Enter to continue...")

if __name__ == "__main__":
	main()

