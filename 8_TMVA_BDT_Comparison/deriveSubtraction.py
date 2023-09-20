# -*- coding: utf-8 -*-
from ROOT import *
import numpy as np
import pandas
import math
import csv
from array import array

def main():

        xgb = pandas.read_csv('../7_Full_Data_Analysis/xgb_w.csv', sep=',', header=None)
        bdt = pandas.read_csv('bdt_w.csv', sep=',', header=None)

	xch = np.array(xgb, dtype=float)
	bch = np.array(bdt, dtype=float)

	subt = []
	for i in range(len(xch[0])): subt.append(xch[0][i]-bch[0][i])

	h1 = TH1D("h1","AUC Subtraction Distribution (XGB - BDT)", 30,-0.005,0.01) 

	c1 = TCanvas("c1","AUC Subtraction Distribution XGB - BDT", 900,600)

	for i in range(len(subt)):
		h1.Fill(subt[i])
		if subt[i] < 0: print(subt[i])

	h1.SetMarkerStyle(kFullSquare)
	h1.SetMarkerColor(kBlue)
	h1.SetMarkerSize(1)

	h1.GetXaxis().SetTitle("AUC Subtraction (XGB - BDT)")
	h1.GetYaxis().SetTitle("Entry")

#	gStyle.SetOptFit(111)
#	f1 = h1.Fit("gaus")
#	f2 = TF1("fit",f1)
#	f2.SetFillColor(6)
#	f2.SetFillStyle(3011)

	ran = [-0.001,-0.0005,0]
	ara = array('d',ran)

	h2 = TH1D("h2","AUC Subtraction Distribution (XGB - BDT)", 2,ara)
	h2.SetBinContent(1, 2)
	h2.SetBinContent(2, 3)
	h2.SetLineColor(6)
	h2.SetFillColor(6)
	h2.SetLineStyle(7)
	h2.SetFillStyle(3002)

	t1 = TLine(0, 0, 0, 26.8)
	t1.SetLineColor(6)
	t1.SetLineStyle(7)

	h1.Draw("hist e1")
	h2.Draw("same hist")
	t1.Draw("same")

	c1.Draw()
	c1.SaveAs("subt.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
	main()
