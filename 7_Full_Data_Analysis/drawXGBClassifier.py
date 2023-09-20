# -*- coding: utf-8 -*-
from ROOT import *
import numpy as np
import pandas
import math
import csv


def main():

	gStyle.SetOptStat(0)

        auc1 = pandas.read_csv('xgb_w.csv', sep=',', header=None)
        auc2 = pandas.read_csv('xgb_wo.csv', sep=',', header=None)

	chun1 = np.array(auc1, dtype=float)
	chun2 = np.array(auc2, dtype=float)

	x,y1,y2 = [],[],[]
	for i in range(len(chun1[0])):
		x.append(float(i))
		y1.append(chun1[0][i])
		y2.append(chun2[0][i])

	g1 = TH1D("g1","AUC Values from Multiple Trials_XGB", len(chun1[0]),0,len(chun1[0]))
	g2 = TH1D("g2","AUC Values from Multiple Trials_XGB", len(chun2[0]),0,len(chun2[0]))

	c = TCanvas("c","AUC Values from Multiple Trials_XGB", 900,600)

	for i in range(1,len(chun1[0])+1):
		g1.SetBinContent(i,y1[i-1])
		g1.SetBinError(i,1e-5)
		g2.SetBinContent(i,y2[i-1])
		g2.SetBinError(i,1e-5)

	l1 = TLine(0,np.mean(y1),100,np.mean(y1))
        l1.SetLineColor(kBlack)

	l2 = TLine(0,np.mean(y2),100,np.mean(y2))
        l2.SetLineColor(kRed)

	t1 = TText(50, 0.95, "AUC with H & T = {0:.3f}".format(np.mean(y1)))
	t2 = TText(50, 0.94, "AUC without H & T = {0:.3f}".format(np.mean(y2)))

	g1.SetMarkerStyle(kFullSquare)
	g1.SetMarkerSize(0.5)
	g1.SetMarkerColor(kBlack)
	g1.SetLineColor(kBlack)

	g2.SetMarkerStyle(kFullSquare)
	g2.SetMarkerSize(0.5)
	g2.SetMarkerColor(kRed)
	g2.SetLineColor(kRed)

	g1.GetYaxis().SetRangeUser(0.9,1)
        g1.GetXaxis().SetTitle("Trial")
        g1.GetYaxis().SetTitle("AUC")

	t2.SetTextColor(kRed)
	t1.SetTextSize(0.05)
	t2.SetTextSize(0.05)

	g1.Draw()
	g2.Draw("same")
	t1.Draw("same")
	t2.Draw("same")
	l1.Draw("same")
	l2.Draw("same")

	c.SetGrid()
	c.Draw()
	c.SaveAs("xgb.png")

	input("Press Enter to continue...")


if __name__ == "__main__":
	main()
