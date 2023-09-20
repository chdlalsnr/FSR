# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
import sys

def main():

	nominal = open('nominal.csv','r')
	nnline = csv.reader(nominal)

	up = open('up.csv','r')
	uuline = csv.reader(up)

	down = open('down.csv','r')
	ddline = csv.reader(down)

	nn,uu,dd = [],[],[]
	for line in nnline:
		if 'jet1_mass' in line: continue
		nn.append(float(line[24]))

	for line in uuline:
		if 'jet1_mass' in line: continue
		uu.append(float(line[24]))

	for line in ddline:
		if 'jet1_mass' in line: continue
		dd.append(float(line[24]))

	c = TCanvas("c","Histogram of Invariant Mass of H->bb", 900,600)
	h1 = TH1D("h1","Histogram of Invariant Mass of H->bb", 50,0,200)
	h2 = TH1D("h2","Histogram of Invariant Mass of H->bb", 50,0,200)
	h3 = TH1D("h3","Histogram of Invariant Mass of H->bb", 50,0,200)

	for i in range(len(nn)):
		h1.Fill(nn[i])

	for i in range(len(uu)):
		h2.Fill(uu[i])

	for i in range(len(dd)):
		h3.Fill(dd[i])

	nn = sorted(nn)
	uu = sorted(uu)
	dd = sorted(dd)

	median_nn = nn[int(len(nn)/2)]
	median_uu = uu[int(len(uu)/2)]
	median_dd = dd[int(len(dd)/2)]

        l1 = TLine(median_nn,0, median_nn,1850)
        l2 = TLine(median_uu,0, median_uu,1850)
        l3 = TLine(median_dd,0, median_dd,1850)

	tt = TText(140,1400, "JEC Uncertainty: <{0:.1f}%".format( (abs(median_uu-median_dd)/2)*100/median_nn ))

	h1.SetLineColor(kRed)
	l1.SetLineColor(kRed)
	h2.SetLineColor(kBlue)
	l2.SetLineColor(kBlue)
	h3.SetLineColor(kGreen)
	l3.SetLineColor(kGreen)

	l1.SetLineStyle(9)
	l2.SetLineStyle(9)
	l3.SetLineStyle(9)

	c.cd()

	h3.Draw("hist")
	h1.Draw("same")
	h2.Draw("same")

	l1.Draw("same")
	l2.Draw("same")
	l3.Draw("same")

	tt.Draw("same")

	h3.GetXaxis().SetTitle("Higgs Mass[GeV]")
	h3.GetYaxis().SetTitle("Entries")

	c.SaveAs("jec.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()

