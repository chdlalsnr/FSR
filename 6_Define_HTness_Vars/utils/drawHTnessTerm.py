# -*- coding: utf-8 -*-
from ROOT import *
from glob import glob
import numpy as np
import math
import csv
import sys


def main():


	t1_hh,t2_hh,t3_hh,t4_hh = [],[],[],[]
	t1_tt,t2_tt,t3_tt,t4_tt = [],[],[],[]
	h1_hh,h2_hh,h3_hh,h4_hh = [],[],[],[]
	h1_tt,h2_tt,h3_tt,h4_tt = [],[],[],[]

	hh = open('htness_hh.csv','r')
	hhline = csv.reader(hh)

	tt = open('htness_tt.csv','r')
	ttline = csv.reader(tt)

	for line in hhline:
		t1_hh.append(float(line[-9]))
		t2_hh.append(float(line[-8]))
		t3_hh.append(float(line[-7]))
		t4_hh.append(float(line[-6]))
		h1_hh.append(float(line[-5]))
		h2_hh.append(float(line[-4]))
		h3_hh.append(float(line[-3]))
		h4_hh.append(float(line[-2]))

	for line in ttline:
		t1_tt.append(float(line[-9]))
                t2_tt.append(float(line[-8]))
                t3_tt.append(float(line[-7]))
                t4_tt.append(float(line[-6]))
                h1_tt.append(float(line[-5]))
                h2_tt.append(float(line[-4]))
                h3_tt.append(float(line[-3]))
		h4_tt.append(float(line[-2]))

	t1 = t1_hh + t1_tt
	t2 = t2_hh + t2_tt
	t3 = t3_hh + t3_tt
	t4 = t4_hh + t4_tt
	h1 = h1_hh + h1_tt
	h2 = h2_hh + h2_tt
	h3 = h3_hh + h3_tt
	h4 = h4_hh + h4_tt

	t1.sort(reverse=True)
	t2.sort(reverse=True)
	t3.sort(reverse=True)
	t4.sort(reverse=True)
	h1.sort(reverse=True)
	h2.sort(reverse=True)
	h3.sort(reverse=True)
	h4.sort(reverse=True)

	t1_max = t1[len(t1)//10]
	t2_max = t2[len(t2)//10]
	t3_max = t3[len(t3)//10]
	t4_max = t4[len(t4)//10]
	h1_max = h1[len(h1)//10]
	h2_max = h2[len(h2)//10]
	h3_max = h3[len(h3)//10]
	h4_max = h4[len(h4)//10]

	c1 = TCanvas("c1","Distributions on Terms of Topness_DL", 900,600)
	c2 = TCanvas("c2","Distributions on Terms of Higgsness_DL", 900,600)

	hh_t1 = TH1D("hist_t1","Topness First Term Distribution", 40,0,t1_max)
	hh_t2 = TH1D("hist_t2","Topness Second Term Distribution", 40,0,t2_max)
	hh_t3 = TH1D("hist_t3","Topness Third Term Distribution", 40,0,t3_max)
	hh_t4 = TH1D("hist_t4","Topness Fourth Term Distribution", 40,0,t4_max)

	tt_t1 = TH1D("hist_t1","Topness First Term Distribution", 40,0,t1_max)
        tt_t2 = TH1D("hist_t2","Topness Second Term Distribution", 40,0,t2_max)
        tt_t3 = TH1D("hist_t3","Topness Third Term Distribution", 40,0,t3_max)
        tt_t4 = TH1D("hist_t4","Topness Fourth Term Distribution", 40,0,t4_max)

	hh_h1 = TH1D("hist_h1","Higgsness First Term Distribution", 40,0,h1_max)
        hh_h2 = TH1D("hist_h2","Higgsness Second Term Distribution", 40,0,h2_max)
        hh_h3 = TH1D("hist_h3","Higgsness Third Term Distribution", 40,0,h3_max)
	hh_h4 = TH1D("hist_h4","Higgsness Fourth Term Distribution", 40,0,h4_max)

        tt_h1 = TH1D("hist_h1","Higgsness First Term Distribution", 40,0,h1_max)
        tt_h2 = TH1D("hist_h2","Higgsness Second Term Distribution", 40,0,h2_max)
        tt_h3 = TH1D("hist_h3","Higgsness Third Term Distribution", 40,0,h3_max)
	tt_h4 = TH1D("hist_h4","Higgsness Fourth Term Distribution", 40,0,h4_max)

	for i in range(len(h1_hh)):
		if h1_hh[i] > h1_max: h1_hh[i] = h1_max - 0.01
		if h2_hh[i] > h2_max: h2_hh[i] = h2_max - 0.01
		if h3_hh[i] > h3_max: h3_hh[i] = h3_max - 0.01
		if h4_hh[i] > h4_max: h4_hh[i] = h4_max - 0.01

	for i in range(len(h1_tt)):
		if h1_tt[i] > h1_max: h1_tt[i] = h1_max - 0.01
		if h2_tt[i] > h2_max: h2_tt[i] = h2_max - 0.01
		if h3_tt[i] > h3_max: h3_tt[i] = h3_max - 0.01
		if h4_tt[i] > h4_max: h4_tt[i] = h4_max - 0.01

	for i in range(len(t1_hh)):
		if t1_hh[i] > t1_max: t1_hh[i] = t1_max - 0.01
		if t2_hh[i] > t2_max: t2_hh[i] = t2_max - 0.01
		if t3_hh[i] > t3_max: t3_hh[i] = t3_max - 0.01
		if t4_hh[i] > t4_max: t4_hh[i] = t4_max - 0.01

	for i in range(len(t1_tt)):
		if t1_tt[i] > t1_max: t1_tt[i] = t1_max - 0.01
		if t2_tt[i] > t1_max: t1_tt[i] = t1_max - 0.01
		if t3_tt[i] > t1_max: t1_tt[i] = t1_max - 0.01
		if t4_tt[i] > t1_max: t1_tt[i] = t1_max - 0.01


	for i in range(len(t1_hh)):
		hh_t1.Fill(t1_hh[i])
		hh_t2.Fill(t2_hh[i])
		hh_t3.Fill(t3_hh[i])
		hh_t4.Fill(t4_hh[i])
		hh_h1.Fill(h1_hh[i])
		hh_h2.Fill(h2_hh[i])
		hh_h3.Fill(h3_hh[i])
		hh_h4.Fill(h4_hh[i])

	for i in range(len(t1_tt)):
		tt_t1.Fill(t1_tt[i])
                tt_t2.Fill(t2_tt[i])
                tt_t3.Fill(t3_tt[i])
                tt_t4.Fill(t4_tt[i])
                tt_h1.Fill(h1_tt[i])
                tt_h2.Fill(h2_tt[i])
                tt_h3.Fill(h3_tt[i]) 
		tt_h4.Fill(h4_tt[i])

	hh_t1.Scale(1/hh_t1.GetEntries())
	hh_t2.Scale(1/hh_t2.GetEntries())
	hh_t3.Scale(1/hh_t3.GetEntries())
	hh_t4.Scale(1/hh_t4.GetEntries())

	tt_t1.Scale(1/tt_t1.GetEntries())
	tt_t2.Scale(1/tt_t2.GetEntries())
	tt_t3.Scale(1/tt_t3.GetEntries())
	tt_t4.Scale(1/tt_t4.GetEntries())

	hh_h1.Scale(1/hh_h1.GetEntries())
	hh_h2.Scale(1/hh_h2.GetEntries())
	hh_h3.Scale(1/hh_h3.GetEntries())
	hh_h4.Scale(1/hh_h4.GetEntries())

	tt_h1.Scale(1/tt_h1.GetEntries())
	tt_h2.Scale(1/tt_h2.GetEntries())
	tt_h3.Scale(1/tt_h3.GetEntries())
	tt_h4.Scale(1/tt_h4.GetEntries())

	hh_t1_scaled = TH1D("hist_t1","Topness First Term Distribution", 40,0,t1_max)
	hh_t2_scaled = TH1D("hist_t2","Topness Second Term Distribution", 40,0,t2_max)
	hh_t3_scaled = TH1D("hist_t3","Topness Third Term Distribution", 40,0,t3_max)
	hh_t4_scaled = TH1D("hist_t4","Topness Fourth Term Distribution", 40,0,t4_max)

	tt_t1_scaled = TH1D("hist_t1","Topness First Term Distribution", 40,0,t1_max)
        tt_t2_scaled = TH1D("hist_t2","Topness Second Term Distribution", 40,0,t2_max)
        tt_t3_scaled = TH1D("hist_t3","Topness Third Term Distribution", 40,0,t3_max)
        tt_t4_scaled = TH1D("hist_t4","Topness Fourth Term Distribution", 40,0,t4_max)

        hh_h1_scaled = TH1D("hist_h1","Higgsness First Term Distribution", 40,0,h1_max)
        hh_h2_scaled = TH1D("hist_h2","Higgsness Second Term Distribution", 40,0,h2_max)
        hh_h3_scaled = TH1D("hist_h3","Higgsness Third Term Distribution", 40,0,h3_max)
	hh_h4_scaled = TH1D("hist_h4","Higgsness Fourth Term Distribution", 40,0,h4_max)

        tt_h1_scaled = TH1D("hist_h1","Higgsness First Term Distribution", 40,0,h1_max)
        tt_h2_scaled = TH1D("hist_h2","Higgsness Second Term Distribution", 40,0,h2_max)
        tt_h3_scaled = TH1D("hist_h3","Higgsness Third Term Distribution", 40,0,h3_max)
	tt_h4_scaled = TH1D("hist_h4","Higgsness Fourth Term Distribution", 40,0,h4_max)


	for i in range(1,hh_h1.GetNbinsX()+1):
		hh_h1_scaled.SetBinContent(i,hh_h1.GetBinContent(i))
		hh_h2_scaled.SetBinContent(i,hh_h2.GetBinContent(i))
		hh_h3_scaled.SetBinContent(i,hh_h3.GetBinContent(i))
		hh_h4_scaled.SetBinContent(i,hh_h4.GetBinContent(i))

	for i in range(1,tt_h1.GetNbinsX()+1):
		tt_h1_scaled.SetBinContent(i,tt_h1.GetBinContent(i))
		tt_h2_scaled.SetBinContent(i,tt_h2.GetBinContent(i))
		tt_h3_scaled.SetBinContent(i,tt_h3.GetBinContent(i))
		tt_h4_scaled.SetBinContent(i,tt_h4.GetBinContent(i))

	for i in range(1,hh_t1.GetNbinsX()+1):
		hh_t1_scaled.SetBinContent(i,hh_t1.GetBinContent(i))
		hh_t2_scaled.SetBinContent(i,hh_t2.GetBinContent(i))
		hh_t3_scaled.SetBinContent(i,hh_t3.GetBinContent(i))
		hh_t4_scaled.SetBinContent(i,hh_t4.GetBinContent(i))

	for i in range(1,tt_t1.GetNbinsX()+1):
		tt_t1_scaled.SetBinContent(i,tt_t1.GetBinContent(i))
		tt_t2_scaled.SetBinContent(i,tt_t2.GetBinContent(i))
		tt_t3_scaled.SetBinContent(i,tt_t3.GetBinContent(i))
		tt_t4_scaled.SetBinContent(i,tt_t4.GetBinContent(i))


	hh_t1_scaled.GetXaxis().SetTitle("T First [GeV]")
	hh_t1_scaled.GetYaxis().SetTitle("Entries")
	hh_t1_scaled.SetLineColor(2)
	hh_t1_scaled.SetFillColor(2)
	hh_t1_scaled.SetFillStyle(3004)
	tt_t1_scaled.SetFillColor(4)
	tt_t1_scaled.SetFillStyle(3005)

	hh_t2_scaled.GetXaxis().SetTitle("T Second [GeV]")
        hh_t2_scaled.GetYaxis().SetTitle("Entries")
	hh_t2_scaled.SetLineColor(2)
        hh_t2_scaled.SetFillColor(2)
        hh_t2_scaled.SetFillStyle(3004)
        tt_t2_scaled.SetFillColor(4)
        tt_t2_scaled.SetFillStyle(3005)

	hh_t3_scaled.GetXaxis().SetTitle("T Third [GeV]")
        hh_t3_scaled.GetYaxis().SetTitle("Entries")
	hh_t3_scaled.SetLineColor(2)
        hh_t3_scaled.SetFillColor(2)
        hh_t3_scaled.SetFillStyle(3004)
        tt_t3_scaled.SetFillColor(4)
        tt_t3_scaled.SetFillStyle(3005)

	hh_t4_scaled.GetXaxis().SetTitle("T Fourth [GeV]")
        hh_t4_scaled.GetYaxis().SetTitle("Entries")
	hh_t4_scaled.SetLineColor(2)
        hh_t4_scaled.SetFillColor(2)
        hh_t4_scaled.SetFillStyle(3004)
        tt_t4_scaled.SetFillColor(4)
        tt_t4_scaled.SetFillStyle(3005)

	hh_h1_scaled.GetXaxis().SetTitle("H First [GeV]")
	hh_h1_scaled.GetYaxis().SetTitle("Entries")
	hh_h1_scaled.SetLineColor(2)
        hh_h1_scaled.SetFillColor(2)
        hh_h1_scaled.SetFillStyle(3004)
        tt_h1_scaled.SetFillColor(4)
        tt_h1_scaled.SetFillStyle(3005)

	hh_h2_scaled.GetXaxis().SetTitle("H Second [GeV]")
        hh_h2_scaled.GetYaxis().SetTitle("Entries")
	hh_h2_scaled.SetLineColor(2)
        hh_h2_scaled.SetFillColor(2)
        hh_h2_scaled.SetFillStyle(3004)
        tt_h2_scaled.SetFillColor(4)
        tt_h2_scaled.SetFillStyle(3005)

	hh_h3_scaled.GetXaxis().SetTitle("H Third [GeV]")
        hh_h3_scaled.GetYaxis().SetTitle("Entries")
	hh_h3_scaled.SetLineColor(2)
        hh_h3_scaled.SetFillColor(2)
        hh_h3_scaled.SetFillStyle(3004)
        tt_h3_scaled.SetFillColor(4)
        tt_h3_scaled.SetFillStyle(3005)

	hh_h4_scaled.GetXaxis().SetTitle("H Fourth [GeV]")
        hh_h4_scaled.GetYaxis().SetTitle("Entries")
	hh_h4_scaled.SetLineColor(2)
        hh_h4_scaled.SetFillColor(2)
        hh_h4_scaled.SetFillStyle(3004)
        tt_h4_scaled.SetFillColor(4)
        tt_h4_scaled.SetFillStyle(3005)


	c1.Divide(2,2)

	c1.cd(1)
	hh_t1_scaled.Draw()
	tt_t1_scaled.Draw("same")

	c1.cd(2)
	hh_t2_scaled.Draw()
	tt_t2_scaled.Draw("same")

	c1.cd(3)
	hh_t3_scaled.Draw()
	tt_t3_scaled.Draw("same")

	c1.cd(4)
	hh_t4_scaled.Draw()
	tt_t4_scaled.Draw("same")

	c1.SaveAs("topness_terms.png")

	c2.Divide(2,2)

	c2.cd(1)
	hh_h1_scaled.Draw()
	tt_h1_scaled.Draw("same")

	c2.cd(2)
	hh_h2_scaled.Draw()
	tt_h2_scaled.Draw("same")

	c2.cd(3)
	hh_h3_scaled.Draw()
	tt_h3_scaled.Draw("same")

	c2.cd(4)
	hh_h4_scaled.Draw()
	tt_h4_scaled.Draw("same")

	c2.SaveAs("higgsness_terms.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()


