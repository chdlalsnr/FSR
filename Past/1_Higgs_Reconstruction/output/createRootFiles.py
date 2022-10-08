#!/usr/bin/python
from ROOT import *
import numpy as np
import os
import sys
import yaml
import json
import glob
import math

def main():

	#flist = glob.glob('signal.root')
        #for fname in flist:
        f_input = TFile("signal.root","read")
	f_input.cd('signal')

	f = TFile("test.root","recreate")	
	t = TTree("signal","signal")
	h1 = gROOT.FindObject("h1_mass")

	h1_mass = TH1D("h1_mass","h1_mass",40,0,200)
	t.Branch("H1_mass","h1_mass")

	for i in range(1,h1_mass.GetNbinsX()+1):
		h1_mass.SetBinContent(i,(h1_mass[i])*1.639*100)

	cvs = TCanvas()
	cvs.cd()
	h1_mass.Draw()

	t.Fill()
	t.Write()
	f.Close()

	input("Press Enter to continue...")

 
if __name__ == '__main__':
        main()
