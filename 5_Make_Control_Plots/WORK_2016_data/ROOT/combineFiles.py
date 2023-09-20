from ROOT import *
import numpy as np
import math
import glob


def main():

	hreal_j1_mass = TH1F("hreal_j1_mass","Real J1 Mass Distribution", 40,0,200)
        hreal_j1_pt   = TH1F("hreal_j1_pt","Real J1 Pt Distribution", 40,0,200)
        hreal_j1_eta  = TH1F("hreal_j1_eta","Real J1 Eta Distribution", 40,-3.4,3.4)
        hreal_j1_phi  = TH1F("hreal_j1_phi","Real J1 Phi Distribution", 40,-3.14,3.14)
        hreal_j1_btag = TH1F("hreal_j1_btag","Real J1 Btag Distribution", 40,0,1.01)

        hreal_j2_mass = TH1F("hreal_j2_mass","Real J2 Mass Distribution", 40,0,200)
        hreal_j2_pt   = TH1F("hreal_j2_pt","Real J2 Pt Distribution", 40,0,400)
        hreal_j2_eta  = TH1F("hreal_j2_eta","Real J2 Eta Distribution", 40,-3.4,3.4)
        hreal_j2_phi  = TH1F("hreal_j2_phi","Real J2 Phi Distribution", 40,-3.14,3.14)
        hreal_j2_btag = TH1F("hreal_j2_btag","Real J2 Btag Distribution", 40,0,1.01)

	hreal_l1_mass = TH1F("hreal_l1_mass","Real L1 Mass Distribution", 40,0,200)
        hreal_l1_pt   = TH1F("hreal_l1_pt","Real L1 Pt Distribution", 40,0,400)
        hreal_l1_eta  = TH1F("hreal_l1_eta","Real L1 Eta Distribution", 40,-3.4,3.4)
        hreal_l1_phi  = TH1F("hreal_l1_phi","Real L1 Phi Distribution", 40,-3.14,3.14)
        hreal_l1_charge = TH1F("hreal_l1_charge","Real L1 Charge Distribution", 40,0,1.01)

        hreal_l2_mass = TH1F("hreal_l2_mass","Real L2 Mass Distribution", 40,0,200)
        hreal_l2_pt   = TH1F("hreal_l2_pt","Real L2 Pt Distribution", 40,0,400)
        hreal_l2_eta  = TH1F("hreal_l2_eta","Real L2 Eta Distribution", 40,-3.4,3.4)
        hreal_l2_phi  = TH1F("hreal_l2_phi","Real L2 Phi Distribution", 40,-3.14,3.14)
        hreal_l2_charge = TH1F("hreal_l2_charge","Real L2 Charge Distribution", 40,0,1.01)

	hreal_dr_ll   = TH1F("hreal_dr_ll","Real DeltaR Between LL", 40,0,10)
	hreal_dr_jj   = TH1F("hreal_dr_jj","Real DeltaR Between JJ", 40,0,10)

	hreal_h1_mass = TH1F("hreal_h1_mass","Real H1 Mass Distribution", 40,0,200)
	hreal_h1_pt   = TH1F("hreal_h1_pt","Real H1 Pt Distribution", 40,0,400)
	hreal_h1_eta  = TH1F("hreal_h1_eta","Real H1 Eta Distribution", 40,-3.4,3.4)
	hreal_h1_phi  = TH1F("hreal_h1_phi","Real H1 Phi Distribution", 40,-3.14,3.14)

	hreal_h2_mass = TH1F("hreal_h2_mass","Real H2 Mass Distribution", 40,0,200)
	hreal_h2_pt   = TH1F("hreal_h2_pt","Real H2 Pt Distribution", 40,0,400)
	hreal_h2_eta  = TH1F("hreal_h2_eta","Real H2 Eta Distribution", 40,-3.4,3.4)
	hreal_h2_phi  = TH1F("hreal_h2_phi","Real H2 Phi Distribution", 40,-3.14,3.14)

	hreal_hh_mass = TH1F("hreal_hh_mass","Real HH Mass Distribution", 40,0,600)
        hreal_hh_pt   = TH1F("hreal_hh_pt","Real HH Pt Distribution", 40,0,600)
        hreal_hh_eta  = TH1F("hreal_hh_eta","Real HH Eta Distribution", 40,-3.4,3.4)
        hreal_hh_phi  = TH1F("hreal_hh_phi","Real HH Phi Distribution", 40,-3.14,3.14)


	real_j1_mass = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j1_pt   = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j1_eta  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j1_phi  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j1_btag = np.zeros(hreal_j1_btag.GetNbinsX()+1)

        real_j2_mass = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j2_pt   = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j2_eta  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j2_phi  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_j2_btag = np.zeros(hreal_h1_mass.GetNbinsX()+1)

	real_l1_mass = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l1_pt   = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l1_eta  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l1_phi  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l1_charge = np.zeros(hreal_l1_charge.GetNbinsX()+1)

        real_l2_mass = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l2_pt   = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l2_eta  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l2_phi  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_l2_charge = np.zeros(hreal_l1_charge.GetNbinsX()+1)

	real_dr_ll   = np.zeros(hreal_dr_ll.GetNbinsX()+1)
	real_dr_jj   = np.zeros(hreal_dr_ll.GetNbinsX()+1)

	real_h1_mass = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_h1_pt   = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_h1_eta  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_h1_phi  = np.zeros(hreal_h1_mass.GetNbinsX()+1)

        real_h2_mass = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_h2_pt   = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_h2_eta  = np.zeros(hreal_h1_mass.GetNbinsX()+1)
        real_h2_phi  = np.zeros(hreal_h1_mass.GetNbinsX()+1)

	real_hh_mass = np.zeros(hreal_hh_mass.GetNbinsX()+1)
        real_hh_pt   = np.zeros(hreal_hh_mass.GetNbinsX()+1)
        real_hh_eta  = np.zeros(hreal_hh_mass.GetNbinsX()+1)
        real_hh_phi  = np.zeros(hreal_hh_mass.GetNbinsX()+1)


        flist = glob.glob('*.root')
	for fname in flist:

        	f = TFile(fname,"read")
		dr = fname[:-5]
		print("Processing %s file..."%(dr))
        	f.cd(dr)

		j1_mass = gROOT.FindObject("J1_mass")
		j1_pt   = gROOT.FindObject("J1_pt")
		j1_eta  = gROOT.FindObject("J1_eta")
		j1_phi  = gROOT.FindObject("J1_phi")
		j1_btag = gROOT.FindObject("J1_btag")

		j2_mass = gROOT.FindObject("J2_mass")
                j2_pt   = gROOT.FindObject("J2_pt")
                j2_eta  = gROOT.FindObject("J2_eta")
                j2_phi  = gROOT.FindObject("J2_phi")
                j2_btag = gROOT.FindObject("J2_btag")

		l1_mass = gROOT.FindObject("L1_mass")
                l1_pt   = gROOT.FindObject("L1_pt")
                l1_eta  = gROOT.FindObject("L1_eta")
                l1_phi  = gROOT.FindObject("L1_phi")
                l1_charge = gROOT.FindObject("L1_charge")

                l2_mass = gROOT.FindObject("L2_mass")
                l2_pt   = gROOT.FindObject("L2_pt")
                l2_eta  = gROOT.FindObject("L2_eta")
                l2_phi  = gROOT.FindObject("L2_phi")
                l2_charge = gROOT.FindObject("L2_charge")

		deltaR_ll = gROOT.FindObject("DeltaR_ll")
		deltaR_jj = gROOT.FindObject("DeltaR_jj")

		h1_mass = gROOT.FindObject("H1_mass")
                h1_pt   = gROOT.FindObject("H1_pt")
                h1_eta  = gROOT.FindObject("H1_phi")
                h1_phi  = gROOT.FindObject("H1_eta")

                h2_mass = gROOT.FindObject("H2_mass")
                h2_pt   = gROOT.FindObject("H2_pt")
                h2_eta  = gROOT.FindObject("H2_eta")
                h2_phi  = gROOT.FindObject("H2_phi")

		hh_mass = gROOT.FindObject("HH_mass")
		hh_pt   = gROOT.FindObject("HH_pt")
		hh_eta  = gROOT.FindObject("HH_eta")
		hh_phi  = gROOT.FindObject("HH_phi")

		# Real
		for i in range(1,h1_mass.GetNbinsX()+1):

			real_j1_mass[i] = j1_mass.GetBinContent(i)
                        real_j1_pt[i]   = j1_pt.GetBinContent(i)
                        real_j1_eta[i]  = j1_eta.GetBinContent(i)
                        real_j1_phi[i]  = j1_phi.GetBinContent(i)
                        real_j1_btag[i] = j1_btag.GetBinContent(i)

                        real_j2_mass[i] = j2_mass.GetBinContent(i)
                        real_j2_pt[i]   = j2_pt.GetBinContent(i)
                        real_j2_eta[i]  = j2_eta.GetBinContent(i)
                        real_j2_phi[i]  = j2_phi.GetBinContent(i)
                        real_j2_btag[i] = j2_btag.GetBinContent(i)

			real_l1_mass[i] = l1_mass.GetBinContent(i)
                        real_l1_pt[i]   = l1_pt.GetBinContent(i)
                        real_l1_eta[i]  = l1_eta.GetBinContent(i)
                        real_l1_phi[i]  = l1_phi.GetBinContent(i)
                        real_l1_charge[i] = l1_charge.GetBinContent(i)

                        real_l2_mass[i] = l2_mass.GetBinContent(i)
                        real_l2_pt[i]   = l2_pt.GetBinContent(i)
                        real_l2_eta[i]  = l2_eta.GetBinContent(i)
                        real_l2_phi[i]  = l2_phi.GetBinContent(i)
                        real_l2_charge[i] = l2_charge.GetBinContent(i)

			real_dr_ll[i]   = deltaR_ll.GetBinContent(i)
			real_dr_jj[i]   = deltaR_jj.GetBinContent(i)

			real_h1_mass[i] = h1_mass.GetBinContent(i)
                        real_h1_pt[i]   = h1_pt.GetBinContent(i)
                        real_h1_eta[i]  = h1_eta.GetBinContent(i)
                        real_h1_phi[i]  = h1_phi.GetBinContent(i)

			real_h2_mass[i] = h2_mass.GetBinContent(i)
                        real_h2_pt[i]   = h2_pt.GetBinContent(i)
                        real_h2_eta[i]  = h2_eta.GetBinContent(i)
                        real_h2_phi[i]  = h2_phi.GetBinContent(i)

			real_hh_mass[i] = hh_mass.GetBinContent(i)
			real_hh_pt[i]   = hh_pt.GetBinContent(i)
			real_hh_eta[i]  = hh_eta.GetBinContent(i)
			real_hh_phi[i]  = hh_phi.GetBinContent(i)

		for i in range(1,h1_mass.GetNbinsX()+1):

			hreal_j1_mass.SetBinContent(i,real_j1_mass[i])
                        hreal_j1_pt.SetBinContent(i,real_j1_pt[i])
                        hreal_j1_eta.SetBinContent(i,real_j1_eta[i])
                        hreal_j1_phi.SetBinContent(i,real_j1_phi[i])
                        hreal_j1_btag.SetBinContent(i,real_j1_btag[i])

                        hreal_j2_mass.SetBinContent(i,real_j2_mass[i])
                        hreal_j2_pt.SetBinContent(i,real_j2_pt[i])
                        hreal_j2_eta.SetBinContent(i,real_j2_eta[i])
                        hreal_j2_phi.SetBinContent(i,real_j2_phi[i])
                        hreal_j2_btag.SetBinContent(i,real_j2_btag[i])

			hreal_l1_mass.SetBinContent(i,real_l1_mass[i])
                        hreal_l1_pt.SetBinContent(i,real_l1_pt[i])
                        hreal_l1_eta.SetBinContent(i,real_l1_eta[i])
                        hreal_l1_phi.SetBinContent(i,real_l1_phi[i])
                        hreal_l1_charge.SetBinContent(i,real_l1_charge[i])

                        hreal_l2_mass.SetBinContent(i,real_l2_mass[i])
                        hreal_l2_pt.SetBinContent(i,real_l2_pt[i])
                        hreal_l2_eta.SetBinContent(i,real_l2_eta[i])
                        hreal_l2_phi.SetBinContent(i,real_l2_phi[i])
                        hreal_l2_charge.SetBinContent(i,real_l2_charge[i])

			hreal_dr_ll.SetBinContent(i,real_dr_ll[i])
			hreal_dr_jj.SetBinContent(i,real_dr_jj[i])

			hreal_h1_mass.SetBinContent(i,real_h1_mass[i])
			hreal_h1_pt.SetBinContent(i,real_h1_pt[i])
			hreal_h1_eta.SetBinContent(i,real_h1_eta[i])
			hreal_h1_phi.SetBinContent(i,real_h1_phi[i])
	
			hreal_h2_mass.SetBinContent(i,real_h2_mass[i])
			hreal_h2_pt.SetBinContent(i,real_h2_pt[i])
			hreal_h2_eta.SetBinContent(i,real_h2_eta[i])
			hreal_h2_phi.SetBinContent(i,real_h2_phi[i])

			hreal_hh_mass.SetBinContent(i,real_hh_mass[i])
			hreal_hh_pt.SetBinContent(i,real_hh_pt[i])
			hreal_hh_eta.SetBinContent(i,real_hh_eta[i])
			hreal_hh_phi.SetBinContent(i,real_hh_phi[i])

		#NOTE: Error should not be symmetric up and down
		hreal_h1_mass.SetBinError(i,real_h1_mass[i]**0.5)
		hreal_h1_pt.SetBinError(i,  real_h1_pt[i]**0.5)
		hreal_h1_eta.SetBinError(i, real_h1_eta[i]**0.5)
		hreal_h1_phi.SetBinError(i, real_h1_phi[i]**0.5)

		hreal_h2_mass.SetBinError(i,real_h2_mass[i]**0.5)
		hreal_h2_pt.SetBinError(i,  real_h2_pt[i]**0.5)
		hreal_h2_eta.SetBinError(i, real_h2_eta[i]**0.5)
		hreal_h2_phi.SetBinError(i, real_h2_phi[i]**0.5)

		hreal_hh_mass.SetBinError(i,real_hh_mass[i]**0.5)
                hreal_hh_pt.SetBinError(i,  real_hh_pt[i]**0.5)
                hreal_hh_eta.SetBinError(i, real_hh_eta[i]**0.5)
                hreal_hh_phi.SetBinError(i, real_hh_phi[i]**0.5)


	lst_j1_mass = TLegend(0.80,0.70,0.95,0.90)
	lst_j1_pt   = TLegend(0.80,0.70,0.95,0.90)
	lst_j1_eta  = TLegend(0.80,0.70,0.95,0.90)
	lst_j1_phi  = TLegend(0.80,0.70,0.95,0.90)
	lst_j1_btag = TLegend(0.80,0.70,0.95,0.90)

	lst_j2_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_pt   = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_eta  = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_phi  = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_btag = TLegend(0.80,0.70,0.95,0.90)

	lst_l1_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_pt   = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_eta  = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_phi  = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_charge = TLegend(0.80,0.70,0.95,0.90)

        lst_l2_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_pt   = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_eta  = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_phi  = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_charge = TLegend(0.80,0.70,0.95,0.90)

	lst_dr_ll = TLegend(0.80,0.70,0.95,0.90)
	lst_dr_jj = TLegend(0.80,0.70,0.95,0.90)

	lst_h1_mass = TLegend(0.80,0.70,0.95,0.90)
	lst_h1_pt   = TLegend(0.80,0.70,0.95,0.90)
	lst_h1_eta  = TLegend(0.80,0.70,0.95,0.90)
	lst_h1_phi  = TLegend(0.80,0.70,0.95,0.90)

	lst_h2_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_pt   = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_eta  = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_phi  = TLegend(0.80,0.70,0.95,0.90)

	lst_hh_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_hh_pt   = TLegend(0.80,0.70,0.95,0.90)
        lst_hh_eta  = TLegend(0.80,0.70,0.95,0.90)
        lst_hh_phi  = TLegend(0.80,0.70,0.95,0.90)


	lst_j1_mass.AddEntry(hreal_j1_mass,"Real","l")
	lst_j1_pt.AddEntry(hreal_j1_pt,"Real","l")
	lst_j1_eta.AddEntry(hreal_j1_eta,"Real","l")
	lst_j1_phi.AddEntry(hreal_j1_phi,"Real","l")
	lst_j1_btag.AddEntry(hreal_j1_btag,"Real","l")

	lst_j2_mass.AddEntry(hreal_j2_mass,"Real","l")
        lst_j2_pt.AddEntry(hreal_j2_pt,"Real","l")
        lst_j2_eta.AddEntry(hreal_j2_eta,"Real","l")
        lst_j2_phi.AddEntry(hreal_j2_phi,"Real","l")
        lst_j2_btag.AddEntry(hreal_j2_btag,"Real","l")

	lst_l1_mass.AddEntry(hreal_l1_mass,"Real","l")
        lst_l1_pt.AddEntry(hreal_l1_pt,"Real","l")
        lst_l1_eta.AddEntry(hreal_l1_eta,"Real","l")
        lst_l1_phi.AddEntry(hreal_l1_phi,"Real","l")
        lst_l1_charge.AddEntry(hreal_l1_charge,"Real","l")

        lst_l2_mass.AddEntry(hreal_l2_mass,"Real","l")
        lst_l2_pt.AddEntry(hreal_l2_pt,"Real","l")
        lst_l2_eta.AddEntry(hreal_l2_eta,"Real","l")
        lst_l2_phi.AddEntry(hreal_l2_phi,"Real","l")
        lst_l2_charge.AddEntry(hreal_l2_charge,"Real","l")

	lst_dr_ll.AddEntry(hreal_dr_ll,"Real","l")
	lst_dr_jj.AddEntry(hreal_dr_jj,"Real","l")

	lst_h1_mass.AddEntry(hreal_h1_mass,"Real","l")
        lst_h1_pt.AddEntry(hreal_h1_pt,"Real","l")
        lst_h1_eta.AddEntry(hreal_h1_eta,"Real","l")
        lst_h1_phi.AddEntry(hreal_h1_phi,"Real","l")

	lst_h2_mass.AddEntry(hreal_h2_mass,"Real","l")
        lst_h2_pt.AddEntry(hreal_h2_pt,"Real","l")
        lst_h2_eta.AddEntry(hreal_h2_eta,"Real","l")
        lst_h2_phi.AddEntry(hreal_h2_phi,"Real","l")

	lst_hh_mass.AddEntry(hreal_hh_mass,"Real","l")
        lst_hh_pt.AddEntry(hreal_hh_pt,"Real","l")
        lst_hh_eta.AddEntry(hreal_hh_eta,"Real","l")
        lst_hh_phi.AddEntry(hreal_hh_phi,"Real","l")


	lst_j1_mass.SetTextAlign(13)
	lst_j1_pt.SetTextAlign(13)
	lst_j1_eta.SetTextAlign(13)
	lst_j1_phi.SetTextAlign(13)
	lst_j1_btag.SetTextAlign(13)

	lst_j2_mass.SetTextAlign(13)
        lst_j2_pt.SetTextAlign(13)
        lst_j2_eta.SetTextAlign(13)
        lst_j2_phi.SetTextAlign(13)
        lst_j2_btag.SetTextAlign(13)

	lst_l1_mass.SetTextAlign(13)
        lst_l1_pt.SetTextAlign(13)
        lst_l1_eta.SetTextAlign(13)
        lst_l1_phi.SetTextAlign(13)
        lst_l1_charge.SetTextAlign(13)

        lst_l2_mass.SetTextAlign(13)
        lst_l2_pt.SetTextAlign(13)
        lst_l2_eta.SetTextAlign(13)
        lst_l2_phi.SetTextAlign(13)
        lst_l2_charge.SetTextAlign(13)

	lst_dr_ll.SetTextAlign(13)
	lst_dr_jj.SetTextAlign(13)

	lst_h1_mass.SetTextAlign(13)
        lst_h1_pt.SetTextAlign(13)
        lst_h1_eta.SetTextAlign(13)
        lst_h1_phi.SetTextAlign(13)

	lst_h2_mass.SetTextAlign(13)
        lst_h2_pt.SetTextAlign(13)
        lst_h2_eta.SetTextAlign(13)
        lst_h2_phi.SetTextAlign(13)

	lst_hh_mass.SetTextAlign(13)
        lst_hh_pt.SetTextAlign(13)
        lst_hh_eta.SetTextAlign(13)
        lst_hh_phi.SetTextAlign(13)


#	dr = 'combine'
#	ff = TFile("%s.root"%(dr),"recreate")
#        ff.cd()
#        gDirectory.mkdir(dr)
#        ff.cd(dr)
#
#	hst_j1_mass.Write()
#        hst_j1_pt.Write()
#        hst_j1_eta.Write()
#        hst_j1_phi.Write()
#        hst_j1_btag.Write()
#
#	hst_j2_mass.Write()
#        hst_j2_pt.Write()
#        hst_j2_eta.Write()
#        hst_j2_phi.Write()
#        hst_j2_btag.Write()
#
#	hst_l1_mass.Write()
#        hst_l1_pt.Write()
#        hst_l1_eta.Write()
#        hst_l1_phi.Write()
#        hst_l1_charge.Write()
#
#        hst_l2_mass.Write()
#        hst_l2_pt.Write()
#        hst_l2_eta.Write()
#        hst_l2_phi.Write()
#        hst_l2_charge.Write()
#
#	hst_dr_ll.Write()
#	hst_dr_jj.Write()
#
#	hst_h1_mass.Write()
#        hst_h1_pt.Write()
#        hst_h1_eta.Write()
#        hst_h1_phi.Write()
#
#        hst_h2_mass.Write()
#        hst_h2_pt.Write()
#        hst_h2_eta.Write()
#        hst_h2_phi.Write()
#
#        hst_hh_mass.Write()
#        hst_hh_pt.Write()
#        hst_hh_eta.Write()
#        hst_hh_phi.Write()
#
#	ff.Close()


	c1 = TCanvas("c1","Real Histograms of Basic Kinematics hh", 900,600)

	c1.Divide(2,2)
        c1.cd(1)
        hreal_hh_mass.Draw("hist e1")
	hreal_hh_mass.GetXaxis().SetTitle("Invariant Mass[GeV]")
	hreal_hh_mass.GetYaxis().SetTitle("Events/5GeV")
	hreal_hh_mass.SetLineColor(kBlack)

	c1.cd(2)
	hreal_hh_pt.Draw("hist e1")
        hreal_hh_pt.GetXaxis().SetTitle("Transverse Momentum[GeV]")
        hreal_hh_pt.GetYaxis().SetTitle("Events/10GeV")
	hreal_hh_pt.SetLineColor(kBlack)

	c1.cd(3)
	hreal_hh_eta.Draw("hist e1")
        hreal_hh_eta.GetXaxis().SetTitle("Eta")
        hreal_hh_eta.GetYaxis().SetTitle("Events/.17")
	hreal_hh_eta.SetLineColor(kBlack)

	c1.cd(4)
	hreal_hh_phi.Draw("hist e1")
        hreal_hh_phi.GetXaxis().SetTitle("Phi")
        hreal_hh_phi.GetYaxis().SetTitle("Events/.157deg")
	hreal_hh_phi.SetLineColor(kBlack)

        c1.SetGrid()
	c1.Draw()
	c1.SaveAs("hh.png")

	input("Press Enter to continue...")


if __name__ == '__main__':
        main()
