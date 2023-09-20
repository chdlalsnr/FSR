# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
import sys


def main():

	topness_hh,topness_tt = [],[]
	higgsness_hh,higgsness_tt = [],[]

	hh = open('htness_hh.csv','r')
	hhline = csv.reader(hh)

	tt = open('htness_tt.csv','r')
	ttline = csv.reader(tt)

	for line in hhline:
		if 'jet1_mass' in line: continue
		topness_hh.append(np.log10(float(line[-2])))
		higgsness_hh.append(np.log10(float(line[-3])))

	for line in ttline:
		if 'jet1_mass' in line: continue
		topness_tt.append(np.log10(float(line[-2])))
		higgsness_tt.append(np.log10(float(line[-3])))

	c = TCanvas("c","Histogram in (H,T) Space HH Events_DL", 900,600)

	h2 = TH2D("h2","Histogram in (H,T) Space HH Events_DL", 50,0,10, 50,0,5)

	for i in range(len(topness_hh)):
		h2.Fill(higgsness_hh[i],topness_hh[i])

	c.cd()

	h2.GetXaxis().SetTitle("Log H")
	h2.GetYaxis().SetTitle("Log T")
	h2.Draw("colz")

	c.SaveAs("htness_hh.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()


