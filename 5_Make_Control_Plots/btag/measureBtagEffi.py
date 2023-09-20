#python
# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
from array import array


def main():

	gStyle.SetOptStat(0)

        f1 = open('data.csv','r')
	f2 = open('mc.csv','r')

        line1 = csv.reader(f1)
	line2 = csv.reader(f2)

	h1_loose  = TH1D("h1_loose", "Mass Distribution of Leading Jet", 14,0,140)
	h1_medium = TH1D("h1_medium","Mass Distribution of Leading Jet", 14,0,140)
	h2_loose  = TH1D("h2_loose", "Mass Distribution of Leading Jet", 14,0,140)
	h2_medium = TH1D("h2_medium","Mass Distribution of Leading Jet", 14,0,140)

	jet1_pt,jet2_pt = 0,0
	jet1_btag,jet2_btag = 0,0
	for line in line1:

		jet1_pt = float(line[1])
		jet2_pt = float(line[6])
		jet1_btag = float(line[4])
		jet2_btag = float(line[9])

		h1_loose.Fill(jet1_pt)
		if (jet1_btag > 0.3093): h1_medium.Fill(jet1_pt)
		h1_loose.Fill(jet2_pt)
		if (jet2_btag > 0.3093): h1_medium.Fill(jet2_pt)

	jet1_pt,jet2_pt = 0,0
        jet1_btag,jet2_btag = 0,0
        for line in line2:

                jet1_pt = float(line[1])
                jet2_pt = float(line[6])
                jet1_btag = float(line[4])
                jet2_btag = float(line[9])

                h2_loose.Fill(jet1_pt)
                if (jet1_btag > 0.3093): h2_medium.Fill(jet1_pt)
                h2_loose.Fill(jet2_pt)
                if (jet2_btag > 0.3093): h2_medium.Fill(jet2_pt)

	xbins = [0,10,20,30,40,50,60,80,100,120,160]
	abins = array('d',xbins)

	f1 = TH1D("f1", "B-tagging Efficiency Measurement", 10,abins)
	f2 = TH1D("f2", "B-tagging Efficiency Measurement", 10,abins)

	l1 = TLine(20,0,20,1)
	l1.SetLineColor(kBlue)
	l1.SetLineStyle(7)

	for i in range(1,h1_medium.GetNbinsX()+1):
		try:
			f1.SetBinContent(i, h1_medium.GetBinContent(i)/h1_loose.GetBinContent(i))
			f1.SetBinError(i, 0.001)
		except ZeroDivisionError: f1.SetBinContent(i,0)
		
		try: 
			f2.SetBinContent(i, h2_medium.GetBinContent(i)/h2_loose.GetBinContent(i))
			f2.SetBinError(i, 0.001)
		except ZeroDivisionError: f2.SetBinContent(i,0)

	c1 = TCanvas("c1","B-tagging Efficiency Measurement", 900,600)

	a,b = [],[]
	for i in range(7,10):
		a.append(f1.GetBinContent(i))
		b.append(f2.GetBinContent(i))

	t1 = TText(80, 0.5, "Efficiency in data = {0:.1f}%".format(np.mean(a)*100))
	t2 = TText(80, 0.6, "Efficiency in MC = {0:.1f}%".format(np.mean(b)*100))

	l2 = TLine(0,np.mean(a),160,np.mean(a))
	l2.SetLineColor(kBlack)
	l2.SetLineStyle(7)

	l3 = TLine(0,np.mean(b),160,np.mean(b))
	l3.SetLineColor(kRed)
	l3.SetLineStyle(7)

	t2.SetTextColor(kRed)
#	t1.SetTextSize(0.1)

	f2.Draw("hist e1")
	f1.Draw("same hist e1")
	t1.Draw("same")
	t2.Draw("same")
	l1.Draw("same")
	l2.Draw("same")
	l3.Draw("same")

	f2.SetLineColor(kRed)
        f2.SetMarkerColor(kRed)

        f1.SetMarkerStyle(kFullSquare)
        f2.SetMarkerStyle(kFullSquare)
        f1.SetMarkerSize(1)
        f2.SetMarkerSize(1)

        f2.GetYaxis().SetRangeUser(0,1)
        f2.GetXaxis().SetTitle("Transverse momentum [GeV]")

	c1.SetGrid()
	c1.Draw()

	input("Press Enter to continue...")


if __name__ == "__main__":
	main()
