from ROOT import *
import numpy as np
import math
import glob


def main():

	hsig_j1_mass = TH1F("hsig_j1_mass","Signal J1 Mass Distribution", 40,0,200)
        hsig_j1_pt   = TH1F("hsig_j1_pt","Signal J1 Pt Distribution", 40,0,200)
        hsig_j1_eta  = TH1F("hsig_j1_eta","Signal J1 Eta Distribution", 40,-3.4,3.4)
        hsig_j1_phi  = TH1F("hsig_j1_phi","Signal J1 Phi Distribution", 40,-3.14,3.14)
        hsig_j1_btag = TH1F("hsig_j1_btag","Signal J1 Btag Distribution", 40,0,1.01)

        hsig_j2_mass = TH1F("hsig_j2_mass","Signal J2 Mass Distribution", 40,0,200)
        hsig_j2_pt   = TH1F("hsig_j2_pt","Signal J2 Pt Distribution", 40,0,400)
        hsig_j2_eta  = TH1F("hsig_j2_eta","Signal J2 Eta Distribution", 40,-3.4,3.4)
        hsig_j2_phi  = TH1F("hsig_j2_phi","Signal J2 Phi Distribution", 40,-3.14,3.14)
        hsig_j2_btag = TH1F("hsig_j2_btag","Signal J2 Btag Distribution", 40,0,1.01)

	hsig_l1_mass = TH1F("hsig_l1_mass","Signal L1 Mass Distribution", 40,0,200)
        hsig_l1_pt   = TH1F("hsig_l1_pt","Signal L1 Pt Distribution", 40,0,400)
        hsig_l1_eta  = TH1F("hsig_l1_eta","Signal L1 Eta Distribution", 40,-3.4,3.4)
        hsig_l1_phi  = TH1F("hsig_l1_phi","Signal L1 Phi Distribution", 40,-3.14,3.14)
        hsig_l1_charge = TH1F("hsig_l1_charge","Signal L1 Charge Distribution", 40,0,1.01)

        hsig_l2_mass = TH1F("hsig_l2_mass","Signal L2 Mass Distribution", 40,0,200)
        hsig_l2_pt   = TH1F("hsig_l2_pt","Signal L2 Pt Distribution", 40,0,400)
        hsig_l2_eta  = TH1F("hsig_l2_eta","Signal L2 Eta Distribution", 40,-3.4,3.4)
        hsig_l2_phi  = TH1F("hsig_l2_phi","Signal L2 Phi Distribution", 40,-3.14,3.14)
        hsig_l2_charge = TH1F("hsig_l2_charge","Signal L2 Charge Distribution", 40,0,1.01)

	hsig_dr_ll   = TH1F("hsig_dr_ll","Signal DeltaR Between LL", 40,0,10)
	hsig_dr_jj   = TH1F("hsig_dr_jj","Signal DeltaR Between JJ", 40,0,10)

	hsig_h1_mass = TH1F("hsig_h1_mass","Signal H1 Mass Distribution", 40,0,200)
	hsig_h1_pt   = TH1F("hsig_h1_pt","Signal H1 Pt Distribution", 40,0,400)
	hsig_h1_eta  = TH1F("hsig_h1_eta","Signal H1 Eta Distribution", 40,-3.4,3.4)
	hsig_h1_phi  = TH1F("hsig_h1_phi","Signal H1 Phi Distribution", 40,-3.14,3.14)

	hsig_h2_mass = TH1F("hsig_h2_mass","Signal H2 Mass Distribution", 40,0,200)
	hsig_h2_pt   = TH1F("hsig_h2_pt","Signal H2 Pt Distribution", 40,0,400)
	hsig_h2_eta  = TH1F("hsig_h2_eta","Signal H2 Eta Distribution", 40,-3.4,3.4)
	hsig_h2_phi  = TH1F("hsig_h2_phi","Signal H2 Phi Distribution", 40,-3.14,3.14)

	hsig_hh_mass = TH1F("hsig_hh_mass","Signal HH Mass Distribution", 40,0,600)
        hsig_hh_pt   = TH1F("hsig_hh_pt","Signal HH Pt Distribution", 40,0,600)
        hsig_hh_eta  = TH1F("hsig_hh_eta","Signal HH Eta Distribution", 40,-3.4,3.4)
        hsig_hh_phi  = TH1F("hsig_hh_phi","Signal HH Phi Distribution", 40,-3.14,3.14)


	htt_j1_mass = TH1F("htt_j1_mass","TTbar J1 Mass Distribution", 40,0,200)
        htt_j1_pt   = TH1F("htt_j1_pt","TTbar J1 Pt Distribution", 40,0,400)
        htt_j1_eta  = TH1F("htt_j1_eta","TTbar J1 Eta Distribution", 40,-3.4,3.4)
        htt_j1_phi  = TH1F("htt_j1_phi","TTbar J1 Phi Distribution", 40,-3.14,3.14)
        htt_j1_btag = TH1F("htt_j1_btag","TTbar J1 Btag Distribution", 40,0,1.01)

        htt_j2_mass = TH1F("htt_j2_mass","TTbar J2 Mass Distribution", 40,0,200)
        htt_j2_pt   = TH1F("htt_j2_pt","TTbar J2 Pt Distribution", 40,0,400)
        htt_j2_eta  = TH1F("htt_j2_eta","TTbar J2 Eta Distribution", 40,-3.4,3.4)
        htt_j2_phi  = TH1F("htt_j2_phi","TTbar J2 Phi Distribution", 40,-3.14,3.14)
        htt_j2_btag = TH1F("htt_j2_btag","TTbar J2 Btag Distribution", 40,0,1.01)

        htt_l1_mass = TH1F("htt_l1_mass","TTbar L1 Mass Distribution", 40,0,200)
        htt_l1_pt   = TH1F("htt_l1_pt","TTbar L1 Pt Distribution", 40,0,400)
        htt_l1_eta  = TH1F("htt_l1_eta","TTbar L1 Eta Distribution", 40,-3.4,3.4)
        htt_l1_phi  = TH1F("htt_l1_phi","TTbar L1 Phi Distribution", 40,-3.14,3.14)
        htt_l1_charge = TH1F("htt_l1_charge","TTbar L1 Charge Distribution", 40,-1,1.01)

        htt_l2_mass = TH1F("htt_l2_mass","TTbar L2 Mass Distribution", 40,0,200)
        htt_l2_pt   = TH1F("htt_l2_pt","TTbar L2 Pt Distribution", 40,0,400)
        htt_l2_eta  = TH1F("htt_l2_eta","TTbar L2 Eta Distribution", 40,-3.4,3.4)
        htt_l2_phi  = TH1F("htt_l2_phi","TTbar L2 Phi Distribution", 40,-3.14,3.14)
        htt_l2_charge = TH1F("htt_l2_charge","TTbar L2 Charge Distribution", 40,-1,1.01)

        htt_dr_ll   = TH1F("htt_dr_ll","TTbar DeltaR Between LL", 40,0,10)
        htt_dr_jj   = TH1F("htt_dr_jj","TTbar DeltaR Between JJ", 40,0,10)

	htt_h1_mass = TH1F("htt_h1_mass","TTbar H1 Mass Distribution", 40,0,200)
	htt_h1_pt   = TH1F("htt_h1_pt","TTbar H1 Pt Distribution", 40,0,400)
	htt_h1_eta  = TH1F("htt_h1_eta","TTbar H1 Eta Distribution", 40,-3.4,3.4)
	htt_h1_phi  = TH1F("htt_h1_phi","TTbar H1 Phi Distribution", 40,-3.14,3.14)

	htt_h2_mass = TH1F("htt_h2_mass","TTbar H2 Mass Distribution", 40,0,200)
	htt_h2_pt   = TH1F("htt_h2_pt","TTbar H2 Pt Distribution", 40,0,400)
	htt_h2_eta  = TH1F("htt_h2_eta","TTbar H2 Eta Distribution", 40,-3.4,3.4)
	htt_h2_phi  = TH1F("htt_h2_phi","TTbar H2 Phi Distribution", 40,-3.14,3.14)

	htt_hh_mass  = TH1F("htt_hh_mass","TTbar HH Mass Distribution", 40,0,600)
        htt_hh_pt    = TH1F("htt_hh_pt","TTbar HH Pt Distribution", 40,0,600)
        htt_hh_eta   = TH1F("htt_hh_eta","TTbar HH Eta Distribution", 40,-3.4,3.4)
        htt_hh_phi   = TH1F("htt_hh_phi","TTbar HH Phi Distribution", 40,-3.14,3.14)


	hdy_j1_mass = TH1F("hdy_j1_mass","DY J1 Mass Distribution", 40,0,200)
        hdy_j1_pt   = TH1F("hdy_j1_pt","DY J1 Pt Distribution", 40,0,400)
        hdy_j1_eta  = TH1F("hdy_j1_eta","DY J1 Eta Distribution", 40,-3.4,3.4)
        hdy_j1_phi  = TH1F("hdy_j1_phi","DY J1 Phi Distribution", 40,-3.14,3.14)
        hdy_j1_btag = TH1F("hdy_j1_btag","DY J1 Btag Distribution", 40,0,1.01)

        hdy_j2_mass = TH1F("hdy_j2_mass","DY J2 Mass Distribution", 40,0,200)
        hdy_j2_pt   = TH1F("hdy_j2_pt","DY J2 Pt Distribution", 40,0,400)
        hdy_j2_eta  = TH1F("hdy_j2_eta","DY J2 Eta Distribution", 40,-3.4,3.4)
        hdy_j2_phi  = TH1F("hdy_j2_phi","DY J2 Phi Distribution", 40,-3.14,3.14)
        hdy_j2_btag = TH1F("hdy_j2_btag","DY J2 Btag Distribution", 40,0,1.01)

        hdy_l1_mass = TH1F("hdy_l1_mass","DY L1 Mass Distribution", 40,0,200)
        hdy_l1_pt   = TH1F("hdy_l1_pt","DY L1 Pt Distribution", 40,0,400)
        hdy_l1_eta  = TH1F("hdy_l1_eta","DY L1 Eta Distribution", 40,-3.4,3.4)
        hdy_l1_phi  = TH1F("hdy_l1_phi","DY L1 Phi Distribution", 40,-3.14,3.14)
        hdy_l1_charge = TH1F("hdy_l1_charge","DY L1 Charge Distribution", 40,-1,1.01)

        hdy_l2_mass = TH1F("hdy_l2_mass","DY L2 Mass Distribution", 40,0,200)
        hdy_l2_pt   = TH1F("hdy_l2_pt","DY L2 Pt Distribution", 40,0,400)
        hdy_l2_eta  = TH1F("hdy_l2_eta","DY L2 Eta Distribution", 40,-3.4,3.4)
        hdy_l2_phi  = TH1F("hdy_l2_phi","DY L2 Phi Distribution", 40,-3.14,3.14)
        hdy_l2_charge = TH1F("hdy_l2_charge","DY L2 Charge Distribution", 40,-1,1.01)

        hdy_dr_ll   = TH1F("hdy_dr_ll","DY DeltaR Between LL", 40,0,10)
        hdy_dr_jj   = TH1F("hdy_dr_jj","DY DeltaR Between JJ", 40,0,10)

	hdy_h1_mass = TH1F("hdy_h1_mass","DY H1 Mass Distribution", 40,0,200)
	hdy_h1_pt   = TH1F("hdy_h1_pt","DY H1 Pt Distribution", 40,0,400)
	hdy_h1_eta  = TH1F("hdy_h1_eta","DY H1 Eta Distribution", 40,-3.4,3.4)
	hdy_h1_phi  = TH1F("hdy_h1_phi","DY H1 Phi Distribution", 40,-3.14,3.14)

	hdy_h2_mass = TH1F("hdy_h2_mass","DY H2 Mass Distribution", 40,0,200)
	hdy_h2_pt   = TH1F("hdy_h2_pt","DY H2 Pt Distribution", 40,0,400)
	hdy_h2_eta  = TH1F("hdy_h2_eta","DY H2 Eta Distribution", 40,-3.4,3.4)
	hdy_h2_phi  = TH1F("hdy_h2_phi","DY H2 Phi Distribution", 40,-3.14,3.14)

	hdy_hh_mass  = TH1F("hdy_hh_mass","DY HH Mass Distribution", 40,0,600)
	hdy_hh_pt    = TH1F("hdy_hh_pt","DY HH Pt Distribution", 40,0,600)
	hdy_hh_eta   = TH1F("hdy_hh_eta","DY HH Eta Distribution", 40,-3.4,3.4)
	hdy_hh_phi   = TH1F("hdy_hh_phi","DY HH Phi Distribution", 40,-3.14,3.14)


	hot_j1_mass = TH1F("hot_j1_mass","Others J1 Mass Distribution", 40,0,200)
        hot_j1_pt   = TH1F("hot_j1_pt","Others J1 Pt Distribution", 40,0,400)
        hot_j1_eta  = TH1F("hot_j1_eta","Others J1 Eta Distribution", 40,-3.4,3.4)
        hot_j1_phi  = TH1F("hot_j1_phi","Others J1 Phi Distribution", 40,-3.14,3.14)
        hot_j1_btag = TH1F("hot_j1_btag","Others J1 Btag Distribution", 40,0,1.01)

        hot_j2_mass = TH1F("hot_j2_mass","Others J2 Mass Distribution", 40,0,200)
        hot_j2_pt   = TH1F("hot_j2_pt","Others J2 Pt Distribution", 40,0,400)
        hot_j2_eta  = TH1F("hot_j2_eta","Others J2 Eta Distribution", 40,-3.4,3.4)
        hot_j2_phi  = TH1F("hot_j2_phi","Others J2 Phi Distribution", 40,-3.14,3.14)
        hot_j2_btag = TH1F("hot_j2_btag","Others J2 Btag Distribution", 40,0,1.01)

        hot_l1_mass = TH1F("hot_l1_mass","Others L1 Mass Distribution", 40,0,200)
        hot_l1_pt   = TH1F("hot_l1_pt","Others L1 Pt Distribution", 40,0,400)
        hot_l1_eta  = TH1F("hot_l1_eta","Others L1 Eta Distribution", 40,-3.4,3.4)
        hot_l1_phi  = TH1F("hot_l1_phi","Others L1 Phi Distribution", 40,-3.14,3.14)
        hot_l1_charge = TH1F("hot_l1_charge","Others L1 Charge Distribution", 40,-1,1.01)

        hot_l2_mass = TH1F("hot_l2_mass","Others L2 Mass Distribution", 40,0,200)
        hot_l2_pt   = TH1F("hot_l2_pt","Others L2 Pt Distribution", 40,0,400)
        hot_l2_eta  = TH1F("hot_l2_eta","Others L2 Eta Distribution", 40,-3.4,3.4)
        hot_l2_phi  = TH1F("hot_l2_phi","Others L2 Phi Distribution", 40,-3.14,3.14)
        hot_l2_charge = TH1F("hot_l2_charge","Others L2 Charge Distribution", 40,-1,1.01)

        hot_dr_ll   = TH1F("hot_dr_ll","Others DeltaR Between LL", 40,0,10)
        hot_dr_jj   = TH1F("hot_dr_jj","Others DeltaR Between JJ", 40,0,10)

	hot_h1_mass = TH1F("hot_h1_mass","Others H1 Mass Distribution", 40,0,200)
	hot_h1_pt   = TH1F("hot_h1_pt","Others H1 Pt Distribution", 40,0,400)
	hot_h1_eta  = TH1F("hot_h1_eta","Others H1 Eta Distribution", 40,-3.4,3.4)
	hot_h1_phi  = TH1F("hot_h1_phi","Others H1 Phi Distribution", 40,-3.14,3.14)

	hot_h2_mass = TH1F("hot_h2_mass","Others H2 Mass Distribution", 40,0,200)
        hot_h2_pt   = TH1F("hot_h2_pt","Others H2 Pt Distribution", 40,0,400)
        hot_h2_eta  = TH1F("hot_h2_eta","Others H2 Eta Distribution", 40,-3.4,3.4)
        hot_h2_phi  = TH1F("hot_h2_phi","Others H2 Phi Distribution", 40,-3.14,3.14)

	hot_hh_mass = TH1F("hot_hh_mass","Others HH Mass Distribution", 40,0,600)
        hot_hh_pt   = TH1F("hot_hh_pt","Others HH Pt Distribution", 40,0,600)
        hot_hh_eta  = TH1F("hot_hh_eta","Others HH Eta Distribution", 40,-3.4,3.4)
        hot_hh_phi  = TH1F("hot_hh_phi","Others HH Phi Distribution", 40,-3.14,3.14)


	htot_j1_mass = TH1F("htot_j1_mass","Total J1 Mass Distribution", 40,0,200)
        htot_j1_pt   = TH1F("htot_j1_pt","Total J1 Pt Distribution", 40,0,400)
        htot_j1_eta  = TH1F("htot_j1_eta","Total J1 Eta Distribution", 40,-3.4,3.4)
        htot_j1_phi  = TH1F("htot_j1_phi","Total J1 Phi Distribution", 40,-3.14,3.14)
        htot_j1_btag = TH1F("htot_j1_btag","Total J1 Btag Distribution", 40,0,1.01)

        htot_j2_mass = TH1F("htot_j2_mass","Total J2 Mass Distribution", 40,0,200)
        htot_j2_pt   = TH1F("htot_j2_pt","Total J2 Pt Distribution", 40,0,400)
        htot_j2_eta  = TH1F("htot_j2_eta","Total J2 Eta Distribution", 40,-3.4,3.4)
        htot_j2_phi  = TH1F("htot_j2_phi","Total J2 Phi Distribution", 40,-3.14,3.14)
        htot_j2_btag = TH1F("htot_j2_btag","Total J2 Btag Distribution", 40,0,1.01)

        htot_l1_mass = TH1F("htot_l1_mass","Total L1 Mass Distribution", 40,0,200)
        htot_l1_pt   = TH1F("htot_l1_pt","Total L1 Pt Distribution", 40,0,400)
        htot_l1_eta  = TH1F("htot_l1_eta","Total L1 Eta Distribution", 40,-3.4,3.4)
        htot_l1_phi  = TH1F("htot_l1_phi","Total L1 Phi Distribution", 40,-3.14,3.14)
        htot_l1_charge = TH1F("htot_l1_charge","Total L1 Charge Distribution", 40,-1,1.01)

        htot_l2_mass = TH1F("htot_l2_mass","Total L2 Mass Distribution", 40,0,200)
        htot_l2_pt   = TH1F("htot_l2_pt","Total L2 Pt Distribution", 40,0,400)
        htot_l2_eta  = TH1F("htot_l2_eta","Total L2 Eta Distribution", 40,-3.4,3.4)
        htot_l2_phi  = TH1F("htot_l2_phi","Total L2 Phi Distribution", 40,-3.14,3.14)
        htot_l2_charge = TH1F("htot_l2_charge","Total L2 Charge Distribution", 40,-1,1.01)

        htot_dr_ll   = TH1F("htot_dr_ll","Total DeltaR Between LL", 40,0,10)
        htot_dr_jj   = TH1F("htot_dr_jj","Total DeltaR Between JJ", 40,0,10)

        htot_h1_mass = TH1F("htot_h1_mass","Total H1 Mass Distribution", 40,0,200)
        htot_h1_pt   = TH1F("htot_h1_pt","Total H1 Pt Distribution", 40,0,400)
        htot_h1_eta  = TH1F("htot_h1_eta","Total H1 Eta Distribution", 40,-3.4,3.4)
        htot_h1_phi  = TH1F("htot_h1_phi","Total H1 Phi Distribution", 40,-3.14,3.14)

        htot_h2_mass = TH1F("htot_h2_mass","Total H2 Mass Distribution", 40,0,200)
        htot_h2_pt   = TH1F("htot_h2_pt","Total H2 Pt Distribution", 40,0,400)
        htot_h2_eta  = TH1F("htot_h2_eta","Total H2 Eta Distribution", 40,-3.4,3.4)
        htot_h2_phi  = TH1F("htot_h2_phi","Total H2 Phi Distribution", 40,-3.14,3.14)

        htot_hh_mass = TH1F("htot_hh_mass","Total HH Mass Distribution", 40,0,600)
        htot_hh_pt   = TH1F("htot_hh_pt","Total HH Pt Distribution", 40,0,600)
        htot_hh_eta  = TH1F("htot_hh_eta","Total HH Eta Distribution", 40,-3.4,3.4)
        htot_hh_phi  = TH1F("htot_hh_phi","Total HH Phi Distribution", 40,-3.14,3.14)


	sig_j1_mass = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j1_pt   = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j1_eta  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j1_phi  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j1_btag = np.zeros(hsig_j1_btag.GetNbinsX()+1)

        sig_j2_mass = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j2_pt   = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j2_eta  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j2_phi  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_j2_btag = np.zeros(hsig_h1_mass.GetNbinsX()+1)

	sig_l1_mass = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l1_pt   = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l1_eta  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l1_phi  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l1_charge = np.zeros(hsig_l1_charge.GetNbinsX()+1)

        sig_l2_mass = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l2_pt   = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l2_eta  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l2_phi  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_l2_charge = np.zeros(hsig_l1_charge.GetNbinsX()+1)

	sig_dr_ll   = np.zeros(hsig_dr_ll.GetNbinsX()+1)
	sig_dr_jj   = np.zeros(hsig_dr_ll.GetNbinsX()+1)

	sig_h1_mass = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_h1_pt   = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_h1_eta  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_h1_phi  = np.zeros(hsig_h1_mass.GetNbinsX()+1)

        sig_h2_mass = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_h2_pt   = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_h2_eta  = np.zeros(hsig_h1_mass.GetNbinsX()+1)
        sig_h2_phi  = np.zeros(hsig_h1_mass.GetNbinsX()+1)

	sig_hh_mass = np.zeros(hsig_hh_mass.GetNbinsX()+1)
        sig_hh_pt   = np.zeros(hsig_hh_mass.GetNbinsX()+1)
        sig_hh_eta  = np.zeros(hsig_hh_mass.GetNbinsX()+1)
        sig_hh_phi  = np.zeros(hsig_hh_mass.GetNbinsX()+1)


	tt_j1_mass = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j1_pt   = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j1_eta  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j1_phi  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j1_btag = np.zeros(htt_j1_btag.GetNbinsX()+1)

        tt_j2_mass = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j2_pt   = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j2_eta  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j2_phi  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_j2_btag = np.zeros(htt_h1_mass.GetNbinsX()+1)

        tt_l1_mass = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l1_pt   = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l1_eta  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l1_phi  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l1_charge = np.zeros(htt_l1_charge.GetNbinsX()+1)

        tt_l2_mass = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l2_pt   = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l2_eta  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l2_phi  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_l2_charge = np.zeros(htt_l1_charge.GetNbinsX()+1)

        tt_dr_ll   = np.zeros(htt_dr_ll.GetNbinsX()+1)
        tt_dr_jj   = np.zeros(htt_dr_ll.GetNbinsX()+1)

        tt_h1_mass = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_h1_pt   = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_h1_eta  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_h1_phi  = np.zeros(htt_h1_mass.GetNbinsX()+1)

        tt_h2_mass = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_h2_pt   = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_h2_eta  = np.zeros(htt_h1_mass.GetNbinsX()+1)
        tt_h2_phi  = np.zeros(htt_h1_mass.GetNbinsX()+1)

        tt_hh_mass = np.zeros(htt_hh_mass.GetNbinsX()+1)
        tt_hh_pt   = np.zeros(htt_hh_mass.GetNbinsX()+1)
        tt_hh_eta  = np.zeros(htt_hh_mass.GetNbinsX()+1)
        tt_hh_phi  = np.zeros(htt_hh_mass.GetNbinsX()+1)


	dy_j1_mass = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j1_pt   = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j1_eta  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j1_phi  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j1_btag = np.zeros(hdy_j1_btag.GetNbinsX()+1)

        dy_j2_mass = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j2_pt   = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j2_eta  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j2_phi  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_j2_btag = np.zeros(hdy_h1_mass.GetNbinsX()+1)

        dy_l1_mass = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l1_pt   = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l1_eta  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l1_phi  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l1_charge = np.zeros(hdy_l1_charge.GetNbinsX()+1)

        dy_l2_mass = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l2_pt   = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l2_eta  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l2_phi  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_l2_charge = np.zeros(hdy_l1_charge.GetNbinsX()+1)

        dy_dr_ll   = np.zeros(hdy_dr_ll.GetNbinsX()+1)
        dy_dr_jj   = np.zeros(hdy_dr_ll.GetNbinsX()+1)

        dy_h1_mass = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_h1_pt   = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_h1_eta  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_h1_phi  = np.zeros(hdy_h1_mass.GetNbinsX()+1)

        dy_h2_mass = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_h2_pt   = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_h2_eta  = np.zeros(hdy_h1_mass.GetNbinsX()+1)
        dy_h2_phi  = np.zeros(hdy_h1_mass.GetNbinsX()+1)

        dy_hh_mass = np.zeros(hdy_hh_mass.GetNbinsX()+1)
        dy_hh_pt   = np.zeros(hdy_hh_mass.GetNbinsX()+1)
        dy_hh_eta  = np.zeros(hdy_hh_mass.GetNbinsX()+1)
        dy_hh_phi  = np.zeros(hdy_hh_mass.GetNbinsX()+1)


	wj_j1_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j1_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j1_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j1_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j1_btag = np.zeros(hot_j1_btag.GetNbinsX()+1)

        wj_j2_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j2_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j2_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j2_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_j2_btag = np.zeros(hot_h1_mass.GetNbinsX()+1)

        wj_l1_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l1_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l1_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l1_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l1_charge = np.zeros(hot_l1_charge.GetNbinsX()+1)

        wj_l2_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l2_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l2_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l2_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_l2_charge = np.zeros(hot_l1_charge.GetNbinsX()+1)

        wj_dr_ll   = np.zeros(hot_dr_ll.GetNbinsX()+1)
        wj_dr_jj   = np.zeros(hot_dr_ll.GetNbinsX()+1)

        wj_h1_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_h1_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_h1_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_h1_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)

        wj_h2_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_h2_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_h2_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        wj_h2_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)

        wj_hh_mass = np.zeros(hot_hh_mass.GetNbinsX()+1)
        wj_hh_pt   = np.zeros(hot_hh_mass.GetNbinsX()+1)
        wj_hh_eta  = np.zeros(hot_hh_mass.GetNbinsX()+1)
        wj_hh_phi  = np.zeros(hot_hh_mass.GetNbinsX()+1)


	st_j1_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j1_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j1_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j1_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j1_btag = np.zeros(hot_j1_btag.GetNbinsX()+1)

        st_j2_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j2_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j2_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j2_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_j2_btag = np.zeros(hot_h1_mass.GetNbinsX()+1)

        st_l1_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l1_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l1_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l1_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l1_charge = np.zeros(hot_l1_charge.GetNbinsX()+1)

        st_l2_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l2_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l2_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l2_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_l2_charge = np.zeros(hot_l1_charge.GetNbinsX()+1)

        st_dr_ll   = np.zeros(hot_dr_ll.GetNbinsX()+1)
        st_dr_jj   = np.zeros(hot_dr_ll.GetNbinsX()+1)

        st_h1_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_h1_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_h1_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_h1_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)

        st_h2_mass = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_h2_pt   = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_h2_eta  = np.zeros(hot_h1_mass.GetNbinsX()+1)
        st_h2_phi  = np.zeros(hot_h1_mass.GetNbinsX()+1)

        st_hh_mass = np.zeros(hot_hh_mass.GetNbinsX()+1)
        st_hh_pt   = np.zeros(hot_hh_mass.GetNbinsX()+1)
        st_hh_eta  = np.zeros(hot_hh_mass.GetNbinsX()+1)
        st_hh_phi  = np.zeros(hot_hh_mass.GetNbinsX()+1)


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

		# Signal
		sf_signal = (100/299999.0)
		if (dr == 'signal'):
			for i in range(1,h1_mass.GetNbinsX()+1):

				sig_j1_mass[i] = j1_mass.GetBinContent(i)*sf_signal
                                sig_j1_pt[i]   = j1_pt.GetBinContent(i)*sf_signal
                                sig_j1_eta[i]  = j1_eta.GetBinContent(i)*sf_signal
                                sig_j1_phi[i]  = j1_phi.GetBinContent(i)*sf_signal
                                sig_j1_btag[i] = j1_btag.GetBinContent(i)*sf_signal

                                sig_j2_mass[i] = j2_mass.GetBinContent(i)*sf_signal
                                sig_j2_pt[i]   = j2_pt.GetBinContent(i)*sf_signal
                                sig_j2_eta[i]  = j2_eta.GetBinContent(i)*sf_signal
                                sig_j2_phi[i]  = j2_phi.GetBinContent(i)*sf_signal
                                sig_j2_btag[i] = j2_btag.GetBinContent(i)*sf_signal

				sig_l1_mass[i] = l1_mass.GetBinContent(i)*sf_signal
                                sig_l1_pt[i]   = l1_pt.GetBinContent(i)*sf_signal
                                sig_l1_eta[i]  = l1_eta.GetBinContent(i)*sf_signal
                                sig_l1_phi[i]  = l1_phi.GetBinContent(i)*sf_signal
                                sig_l1_charge[i] = l1_charge.GetBinContent(i)*sf_signal

                                sig_l2_mass[i] = l2_mass.GetBinContent(i)*sf_signal
                                sig_l2_pt[i]   = l2_pt.GetBinContent(i)*sf_signal
                                sig_l2_eta[i]  = l2_eta.GetBinContent(i)*sf_signal
                                sig_l2_phi[i]  = l2_phi.GetBinContent(i)*sf_signal
                                sig_l2_charge[i] = l2_charge.GetBinContent(i)*sf_signal

				sig_dr_ll[i]   = deltaR_ll.GetBinContent(i)*sf_signal
				sig_dr_jj[i]   = deltaR_jj.GetBinContent(i)*sf_signal

				sig_h1_mass[i] = h1_mass.GetBinContent(i)*sf_signal
                                sig_h1_pt[i]   = h1_pt.GetBinContent(i)*sf_signal
                                sig_h1_eta[i]  = h1_eta.GetBinContent(i)*sf_signal
                                sig_h1_phi[i]  = h1_phi.GetBinContent(i)*sf_signal

				sig_h2_mass[i] = h2_mass.GetBinContent(i)*sf_signal
                                sig_h2_pt[i]   = h2_pt.GetBinContent(i)*sf_signal
                                sig_h2_eta[i]  = h2_eta.GetBinContent(i)*sf_signal
                                sig_h2_phi[i]  = h2_phi.GetBinContent(i)*sf_signal

				sig_hh_mass[i] = hh_mass.GetBinContent(i)*sf_signal
				sig_hh_pt[i]   = hh_pt.GetBinContent(i)*sf_signal
				sig_hh_eta[i]  = hh_eta.GetBinContent(i)*sf_signal
				sig_hh_phi[i]  = hh_phi.GetBinContent(i)*sf_signal

		for i in range(1,h1_mass.GetNbinsX()+1):

			hsig_j1_mass.SetBinContent(i,sig_j1_mass[i])
                        hsig_j1_pt.SetBinContent(i,sig_j1_pt[i])
                        hsig_j1_eta.SetBinContent(i,sig_j1_eta[i])
                        hsig_j1_phi.SetBinContent(i,sig_j1_phi[i])
                        hsig_j1_btag.SetBinContent(i,sig_j1_btag[i])

                        hsig_j2_mass.SetBinContent(i,sig_j2_mass[i])
                        hsig_j2_pt.SetBinContent(i,sig_j2_pt[i])
                        hsig_j2_eta.SetBinContent(i,sig_j2_eta[i])
                        hsig_j2_phi.SetBinContent(i,sig_j2_phi[i])
                        hsig_j2_btag.SetBinContent(i,sig_j2_btag[i])

			hsig_l1_mass.SetBinContent(i,sig_l1_mass[i])
                        hsig_l1_pt.SetBinContent(i,sig_l1_pt[i])
                        hsig_l1_eta.SetBinContent(i,sig_l1_eta[i])
                        hsig_l1_phi.SetBinContent(i,sig_l1_phi[i])
                        hsig_l1_charge.SetBinContent(i,sig_l1_charge[i])

                        hsig_l2_mass.SetBinContent(i,sig_l2_mass[i])
                        hsig_l2_pt.SetBinContent(i,sig_l2_pt[i])
                        hsig_l2_eta.SetBinContent(i,sig_l2_eta[i])
                        hsig_l2_phi.SetBinContent(i,sig_l2_phi[i])
                        hsig_l2_charge.SetBinContent(i,sig_l2_charge[i])

			hsig_dr_ll.SetBinContent(i,sig_dr_ll[i])
			hsig_dr_jj.SetBinContent(i,sig_dr_jj[i])

			hsig_h1_mass.SetBinContent(i,sig_h1_mass[i])
			hsig_h1_pt.SetBinContent(i,sig_h1_pt[i])
			hsig_h1_eta.SetBinContent(i,sig_h1_eta[i])
			hsig_h1_phi.SetBinContent(i,sig_h1_phi[i])
	
			hsig_h2_mass.SetBinContent(i,sig_h2_mass[i])
			hsig_h2_pt.SetBinContent(i,sig_h2_pt[i])
			hsig_h2_eta.SetBinContent(i,sig_h2_eta[i])
			hsig_h2_phi.SetBinContent(i,sig_h2_phi[i])

			hsig_hh_mass.SetBinContent(i,sig_hh_mass[i])
			hsig_hh_pt.SetBinContent(i,sig_hh_pt[i])
			hsig_hh_eta.SetBinContent(i,sig_hh_eta[i])
			hsig_hh_phi.SetBinContent(i,sig_hh_phi[i])

		#NOTE: Error should not be symmetric up and down
			hsig_h1_mass.SetBinError(i,sig_h1_mass[i]**0.5)
			hsig_h1_pt.SetBinError(i,  sig_h1_pt[i]**0.5)
			hsig_h1_eta.SetBinError(i, sig_h1_eta[i]**0.5)
			hsig_h1_phi.SetBinError(i, sig_h1_phi[i]**0.5)

			hsig_h2_mass.SetBinError(i,sig_h2_mass[i]**0.5)
			hsig_h2_pt.SetBinError(i,  sig_h2_pt[i]**0.5)
			hsig_h2_eta.SetBinError(i, sig_h2_eta[i]**0.5)
			hsig_h2_phi.SetBinError(i, sig_h2_phi[i]**0.5)

			hsig_hh_mass.SetBinError(i,sig_hh_mass[i]**0.5)
                        hsig_hh_pt.SetBinError(i,  sig_hh_pt[i]**0.5)
                        hsig_hh_eta.SetBinError(i, sig_hh_eta[i]**0.5)
                        hsig_hh_phi.SetBinError(i, sig_hh_phi[i]**0.5)


		# TT
		sf_tt_2l2nu = (3211572/4891620336.1)
		sf_tt_semi  = (13296780/32366940321.3)
		sf_tt_had   = (13732740/21432310862.9)
		if (dr == 'tt_2l2nu' or dr == "tt_semi" or dr == "tt_had"):
			if (dr == 'tt_2l2nu'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					tt_j1_mass[i] += j1_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_tt_2l2nu
					tt_j1_btag[i] += j1_btag.GetBinContent(i)*sf_tt_2l2nu

					tt_j2_mass[i] += j2_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_tt_2l2nu
					tt_j2_btag[i] += j2_btag.GetBinContent(i)*sf_tt_2l2nu

					tt_l1_mass[i] += l1_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_tt_2l2nu
					tt_l1_charge[i] += l1_charge.GetBinContent(i)*sf_tt_2l2nu

					tt_l2_mass[i] += l2_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_tt_2l2nu
					tt_l2_charge[i] += l2_charge.GetBinContent(i)*sf_tt_2l2nu

					tt_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_tt_2l2nu
					tt_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_tt_2l2nu

					tt_h1_mass[i] += h1_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_tt_2l2nu
					
					tt_h2_mass[i] += h2_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_tt_2l2nu

					tt_hh_mass[i] += hh_mass.GetBinContent(i)*sf_tt_2l2nu
					tt_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_tt_2l2nu
					tt_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_tt_2l2nu
					tt_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_tt_2l2nu

			if (dr == 'tt_semi'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					tt_j1_mass[i] += j1_mass.GetBinContent(i)*sf_tt_semi
                                        tt_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_tt_semi
                                        tt_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_tt_semi
                                        tt_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_tt_semi
                                        tt_j1_btag[i] += j1_btag.GetBinContent(i)*sf_tt_semi

                                        tt_j2_mass[i] += j2_mass.GetBinContent(i)*sf_tt_semi
                                        tt_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_tt_semi
                                        tt_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_tt_semi
                                        tt_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_tt_semi
                                        tt_j2_btag[i] += j2_btag.GetBinContent(i)*sf_tt_semi

                                        tt_l1_mass[i] += l1_mass.GetBinContent(i)*sf_tt_semi
                                        tt_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_tt_semi
                                        tt_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_tt_semi
                                        tt_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_tt_semi
                                        tt_l1_charge[i] += l1_charge.GetBinContent(i)*sf_tt_semi

                                        tt_l2_mass[i] += l2_mass.GetBinContent(i)*sf_tt_semi
                                        tt_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_tt_semi
                                        tt_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_tt_semi
                                        tt_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_tt_semi
                                        tt_l2_charge[i] += l2_charge.GetBinContent(i)*sf_tt_semi

                                        tt_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_tt_semi
                                        tt_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_tt_semi

					tt_h1_mass[i] += h1_mass.GetBinContent(i)*sf_tt_semi
                                        tt_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_tt_semi
                                        tt_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_tt_semi
                                        tt_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_tt_semi
					
                                        tt_h2_mass[i] += h2_mass.GetBinContent(i)*sf_tt_semi
                                        tt_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_tt_semi
                                        tt_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_tt_semi
                                        tt_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_tt_semi

					tt_hh_mass[i] = hh_mass.GetBinContent(i)*sf_tt_semi
					tt_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_tt_semi
                                        tt_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_tt_semi
                                        tt_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_tt_semi

			if (dr == 'tt_had'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					tt_j1_mass[i] += j1_mass.GetBinContent(i)*sf_tt_had
                                        tt_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_tt_had
                                        tt_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_tt_had
                                        tt_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_tt_had
                                        tt_j1_btag[i] += j1_btag.GetBinContent(i)*sf_tt_had

                                        tt_j2_mass[i] += j2_mass.GetBinContent(i)*sf_tt_had
                                        tt_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_tt_had
                                        tt_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_tt_had
                                        tt_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_tt_had
                                        tt_j2_btag[i] += j2_btag.GetBinContent(i)*sf_tt_had

                                        tt_l1_mass[i] += l1_mass.GetBinContent(i)*sf_tt_had
                                        tt_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_tt_had
                                        tt_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_tt_had
                                        tt_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_tt_had
                                        tt_l1_charge[i] += l1_charge.GetBinContent(i)*sf_tt_had

                                        tt_l2_mass[i] += l2_mass.GetBinContent(i)*sf_tt_had
                                        tt_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_tt_had
                                        tt_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_tt_had
                                        tt_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_tt_had
                                        tt_l2_charge[i] += l2_charge.GetBinContent(i)*sf_tt_had

                                        tt_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_tt_had
                                        tt_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_tt_had

					tt_h1_mass[i] += h1_mass.GetBinContent(i)*sf_tt_had
                                        tt_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_tt_had
                                        tt_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_tt_had
                                        tt_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_tt_had

                                        tt_h2_mass[i] += h2_mass.GetBinContent(i)*sf_tt_had
                                        tt_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_tt_had
                                        tt_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_tt_had
                                        tt_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_tt_had

					tt_hh_mass[i] += hh_mass.GetBinContent(i)*sf_tt_had
                                        tt_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_tt_had
                                        tt_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_tt_had
                                        tt_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_tt_had

			for i in range(1,h1_mass.GetNbinsX()+1):

				htt_j1_mass.SetBinContent(i,(tt_j1_mass[i]))
				htt_j1_pt.SetBinContent(i,(tt_j1_pt[i]))
				htt_j1_eta.SetBinContent(i,(tt_j1_eta[i]))
				htt_j1_phi.SetBinContent(i,(tt_j1_phi[i]))
				htt_j1_btag.SetBinContent(i,(tt_j1_btag[i]))

				htt_j2_mass.SetBinContent(i,(tt_j2_mass[i]))
				htt_j2_pt.SetBinContent(i,(tt_j2_pt[i]))
				htt_j2_eta.SetBinContent(i,(tt_j2_eta[i]))
				htt_j2_phi.SetBinContent(i,(tt_j2_phi[i]))
				htt_j2_btag.SetBinContent(i,(tt_j2_btag[i]))

				htt_l1_mass.SetBinContent(i,(tt_l1_mass[i]))
				htt_l1_pt.SetBinContent(i,(tt_l1_pt[i]))
				htt_l1_eta.SetBinContent(i,(tt_l1_eta[i]))
				htt_l1_phi.SetBinContent(i,(tt_l1_phi[i]))
				htt_l1_charge.SetBinContent(i,(tt_l1_charge[i]))

				htt_l2_mass.SetBinContent(i,(tt_l2_mass[i]))
				htt_l2_pt.SetBinContent(i,(tt_l2_pt[i]))
				htt_l2_eta.SetBinContent(i,(tt_l2_eta[i]))
				htt_l2_phi.SetBinContent(i,(tt_l2_phi[i]))
				htt_l2_charge.SetBinContent(i,(tt_l2_charge[i]))

				htt_dr_ll.SetBinContent(i,(tt_dr_ll[i]))
				htt_dr_jj.SetBinContent(i,(tt_dr_jj[i]))

				htt_h1_mass.SetBinContent(i,(tt_h1_mass[i]))
				htt_h1_pt.SetBinContent(i,(tt_h1_pt[i]))
				htt_h1_eta.SetBinContent(i,(tt_h1_eta[i]))
				htt_h1_phi.SetBinContent(i,(tt_h1_phi[i]))

				htt_h2_mass.SetBinContent(i,(tt_h2_mass[i]))
				htt_h2_pt.SetBinContent(i,(tt_h2_pt[i]))
				htt_h2_eta.SetBinContent(i,(tt_h2_eta[i]))
				htt_h2_phi.SetBinContent(i,(tt_h2_phi[i]))

				htt_hh_mass.SetBinContent(i,(tt_hh_mass[i]))
				htt_hh_pt.SetBinContent(i,(tt_hh_pt[i]))
				htt_hh_eta.SetBinContent(i,(tt_hh_eta[i]))
				htt_hh_phi.SetBinContent(i,(tt_hh_phi[i]))

			#NOTE: Error should not be symmetric up and down
				htt_h1_mass.SetBinError(i,tt_h1_mass[i]**0.5)
				htt_h1_pt.SetBinError(i,  tt_h1_pt[i]**0.5)
				htt_h1_eta.SetBinError(i, tt_h1_eta[i]**0.5)
				htt_h1_phi.SetBinError(i, tt_h1_phi[i]**0.5)

				htt_h2_mass.SetBinError(i,tt_h2_mass[i]**0.5)
				htt_h2_pt.SetBinError(i,  tt_h2_pt[i]**0.5)
				htt_h2_eta.SetBinError(i, tt_h2_eta[i]**0.5)
				htt_h2_phi.SetBinError(i, tt_h2_phi[i]**0.5)

				htt_hh_mass.SetBinError(i,tt_hh_mass[i]**0.5)
				htt_hh_pt.SetBinError(i,  tt_hh_pt[i]**0.5)
				htt_hh_eta.SetBinError(i, tt_hh_eta[i]**0.5)
				htt_hh_phi.SetBinError(i, tt_hh_phi[i]**0.5)


		# DY
		sf_dy_m50 = (220886400/1.89714162185e+12)
		sf_dy_0j  = (175837200/2.83853150541e+11)
		sf_dy_1j  = (32624340/2.24631958677e+11)
		sf_dy_2j  = (12206880/1.02265928169e+11)
		if (dr == 'dy_m50' or dr == 'dy_0j' or dr == 'dy_1j' or dr == 'dy_2j'):

			if (dr == 'dy_m50'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					dy_j1_mass[i] += j1_mass.GetBinContent(i)*sf_dy_m50
                                        dy_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_dy_m50
                                        dy_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_dy_m50
                                        dy_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_dy_m50
                                        dy_j1_btag[i] += j1_btag.GetBinContent(i)*sf_dy_m50

                                        dy_j2_mass[i] += j2_mass.GetBinContent(i)*sf_dy_m50
                                        dy_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_dy_m50
                                        dy_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_dy_m50
                                        dy_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_dy_m50
                                        dy_j2_btag[i] += j2_btag.GetBinContent(i)*sf_dy_m50

                                        dy_l1_mass[i] += l1_mass.GetBinContent(i)*sf_dy_m50
                                        dy_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_dy_m50
                                        dy_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_dy_m50
                                        dy_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_dy_m50
                                        dy_l1_charge[i] += l1_charge.GetBinContent(i)*sf_dy_m50

                                        dy_l2_mass[i] += l2_mass.GetBinContent(i)*sf_dy_m50
                                        dy_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_dy_m50
                                        dy_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_dy_m50
                                        dy_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_dy_m50
                                        dy_l2_charge[i] += l2_charge.GetBinContent(i)*sf_dy_m50

                                        dy_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_dy_m50
                                        dy_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_dy_m50

					dy_h1_mass[i] += h1_mass.GetBinContent(i)*sf_dy_m50
                                        dy_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_dy_m50
                                        dy_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_dy_m50
                                        dy_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_dy_m50
					
                                        dy_h2_mass[i] += h2_mass.GetBinContent(i)*sf_dy_m50
                                        dy_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_dy_m50
                                        dy_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_dy_m50
                                        dy_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_dy_m50

					dy_hh_mass[i] += hh_mass.GetBinContent(i)*sf_dy_m50
                                        dy_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_dy_m50
                                        dy_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_dy_m50
                                        dy_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_dy_m50

			if (dr == 'dy_0j'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					dy_j1_mass[i] += j1_mass.GetBinContent(i)*sf_dy_0j
                                        dy_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_dy_0j
                                        dy_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_dy_0j
                                        dy_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_dy_0j
                                        dy_j1_btag[i] += j1_btag.GetBinContent(i)*sf_dy_0j

                                        dy_j2_mass[i] += j2_mass.GetBinContent(i)*sf_dy_0j
                                        dy_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_dy_0j
                                        dy_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_dy_0j
                                        dy_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_dy_0j
                                        dy_j2_btag[i] += j2_btag.GetBinContent(i)*sf_dy_0j

                                        dy_l1_mass[i] += l1_mass.GetBinContent(i)*sf_dy_0j
                                        dy_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_dy_0j
                                        dy_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_dy_0j
                                        dy_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_dy_0j
                                        dy_l1_charge[i] += l1_charge.GetBinContent(i)*sf_dy_0j

                                        dy_l2_mass[i] += l2_mass.GetBinContent(i)*sf_dy_0j
                                        dy_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_dy_0j
                                        dy_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_dy_0j
                                        dy_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_dy_0j
                                        dy_l2_charge[i] += l2_charge.GetBinContent(i)*sf_dy_0j

                                        dy_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_dy_0j
                                        dy_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_dy_0j

					dy_h1_mass[i] += h1_mass.GetBinContent(i)*sf_dy_0j
                                        dy_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_dy_0j
                                        dy_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_dy_0j
                                        dy_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_dy_0j
					
                                        dy_h2_mass[i] += h2_mass.GetBinContent(i)*sf_dy_0j
                                        dy_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_dy_0j
                                        dy_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_dy_0j
                                        dy_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_dy_0j

					dy_hh_mass[i] += hh_mass.GetBinContent(i)*sf_dy_0j
                                        dy_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_dy_0j
                                        dy_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_dy_0j
                                        dy_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_dy_0j

			if (dr == 'dy_1j'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

					dy_j1_mass[i] += j1_mass.GetBinContent(i)*sf_dy_1j
                                        dy_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_dy_1j
                                        dy_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_dy_1j
                                        dy_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_dy_1j
                                        dy_j1_btag[i] += j1_btag.GetBinContent(i)*sf_dy_1j

                                        dy_j2_mass[i] += j2_mass.GetBinContent(i)*sf_dy_1j
                                        dy_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_dy_1j
                                        dy_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_dy_1j
                                        dy_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_dy_1j
                                        dy_j2_btag[i] += j2_btag.GetBinContent(i)*sf_dy_1j

                                        dy_l1_mass[i] += l1_mass.GetBinContent(i)*sf_dy_1j
                                        dy_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_dy_1j
                                        dy_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_dy_1j
                                        dy_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_dy_1j
                                        dy_l1_charge[i] += l1_charge.GetBinContent(i)*sf_dy_1j

                                        dy_l2_mass[i] += l2_mass.GetBinContent(i)*sf_dy_1j
                                        dy_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_dy_1j
                                        dy_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_dy_1j
                                        dy_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_dy_1j
                                        dy_l2_charge[i] += l2_charge.GetBinContent(i)*sf_dy_1j

                                        dy_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_dy_1j
                                        dy_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_dy_1j

                                        dy_h1_mass[i] += h1_mass.GetBinContent(i)*sf_dy_1j
                                        dy_h1_pt[i] += h1_pt.GetBinContent(i)*sf_dy_1j
                                        dy_h1_eta[i] += h1_eta.GetBinContent(i)*sf_dy_1j
                                        dy_h1_phi[i] += h1_phi.GetBinContent(i)*sf_dy_1j
					
                                        dy_h2_mass[i] += h2_mass.GetBinContent(i)*sf_dy_1j
                                        dy_h2_pt[i] += h2_pt.GetBinContent(i)*sf_dy_1j
                                        dy_h2_eta[i] += h2_eta.GetBinContent(i)*sf_dy_1j
                                        dy_h2_phi[i] += h2_phi.GetBinContent(i)*sf_dy_1j

					dy_hh_mass[i] += hh_mass.GetBinContent(i)*sf_dy_1j
                                        dy_hh_pt[i] += hh_pt.GetBinContent(i)*sf_dy_1j
                                        dy_hh_eta[i] += hh_eta.GetBinContent(i)*sf_dy_1j
                                        dy_hh_phi[i] += hh_phi.GetBinContent(i)*sf_dy_1j

			if (dr == 'dy_2j'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

					dy_j1_mass[i] += j1_mass.GetBinContent(i)*sf_dy_2j
                                        dy_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_dy_2j
                                        dy_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_dy_2j
                                        dy_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_dy_2j
                                        dy_j1_btag[i] += j1_btag.GetBinContent(i)*sf_dy_2j

                                        dy_j2_mass[i] += j2_mass.GetBinContent(i)*sf_dy_2j
                                        dy_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_dy_2j
                                        dy_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_dy_2j
                                        dy_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_dy_2j
                                        dy_j2_btag[i] += j2_btag.GetBinContent(i)*sf_dy_2j

                                        dy_l1_mass[i] += l1_mass.GetBinContent(i)*sf_dy_2j
                                        dy_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_dy_2j
                                        dy_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_dy_2j
                                        dy_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_dy_2j
                                        dy_l1_charge[i] += l1_charge.GetBinContent(i)*sf_dy_2j

                                        dy_l2_mass[i] += l2_mass.GetBinContent(i)*sf_dy_2j
                                        dy_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_dy_2j
                                        dy_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_dy_2j
                                        dy_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_dy_2j
                                        dy_l2_charge[i] += l2_charge.GetBinContent(i)*sf_dy_2j

                                        dy_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_dy_2j
                                        dy_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_dy_2j

                                        dy_h1_mass[i] += h1_mass.GetBinContent(i)*sf_dy_2j
                                        dy_h1_pt[i] += h1_pt.GetBinContent(i)*sf_dy_2j
                                        dy_h1_eta[i] += h1_eta.GetBinContent(i)*sf_dy_2j
                                        dy_h1_phi[i] += h1_phi.GetBinContent(i)*sf_dy_2j
					
                                        dy_h2_mass[i] += h2_mass.GetBinContent(i)*sf_dy_2j
                                        dy_h2_pt[i] += h2_pt.GetBinContent(i)*sf_dy_2j
                                        dy_h2_eta[i] += h2_eta.GetBinContent(i)*sf_dy_2j
                                        dy_h2_phi[i] += h2_phi.GetBinContent(i)*sf_dy_2j

					dy_hh_mass[i] += hh_mass.GetBinContent(i)*sf_dy_2j
                                        dy_hh_pt[i] += hh_pt.GetBinContent(i)*sf_dy_2j
                                        dy_hh_eta[i] += hh_eta.GetBinContent(i)*sf_dy_2j
                                        dy_hh_phi[i] += hh_phi.GetBinContent(i)*sf_dy_2j

			for i in range(1,h1_mass.GetNbinsX()+1):

				hdy_j1_mass.SetBinContent(i,dy_j1_mass[i])
				hdy_j1_pt.SetBinContent(i,dy_j1_pt[i])
				hdy_j1_eta.SetBinContent(i,dy_j1_eta[i])
				hdy_j1_phi.SetBinContent(i,dy_j1_phi[i])
				hdy_j1_btag.SetBinContent(i,dy_j1_btag[i])

				hdy_j2_mass.SetBinContent(i,dy_j2_mass[i])
				hdy_j2_pt.SetBinContent(i,dy_j2_pt[i])
				hdy_j2_eta.SetBinContent(i,dy_j2_eta[i])
				hdy_j2_phi.SetBinContent(i,dy_j2_phi[i])
				hdy_j2_btag.SetBinContent(i,dy_j2_btag[i])

				hdy_l1_mass.SetBinContent(i,dy_l1_mass[i])
				hdy_l1_pt.SetBinContent(i,dy_l1_pt[i])
				hdy_l1_eta.SetBinContent(i,dy_l1_eta[i])
				hdy_l1_phi.SetBinContent(i,dy_l1_phi[i])
				hdy_l1_charge.SetBinContent(i,dy_l1_charge[i])

				hdy_l2_mass.SetBinContent(i,dy_l2_mass[i])
				hdy_l2_pt.SetBinContent(i,dy_l2_pt[i])
				hdy_l2_eta.SetBinContent(i,dy_l2_eta[i])
				hdy_l2_phi.SetBinContent(i,dy_l2_phi[i])
				hdy_l2_charge.SetBinContent(i,dy_l2_charge[i])

				hdy_dr_ll.SetBinContent(i,dy_dr_ll[i])
				hdy_dr_jj.SetBinContent(i,dy_dr_jj[i])

				hdy_h1_mass.SetBinContent(i,dy_h1_mass[i])
				hdy_h1_pt.SetBinContent(i,dy_h1_pt[i])
				hdy_h1_eta.SetBinContent(i,dy_h1_eta[i])
				hdy_h1_phi.SetBinContent(i,dy_h1_phi[i])
				
				hdy_h2_mass.SetBinContent(i,dy_h2_mass[i])
				hdy_h2_pt.SetBinContent(i,dy_h2_pt[i])
				hdy_h2_eta.SetBinContent(i,dy_h2_eta[i])
				hdy_h2_phi.SetBinContent(i,dy_h2_phi[i])

				hdy_hh_mass.SetBinContent(i,dy_hh_mass[i])
                                hdy_hh_pt.SetBinContent(i,dy_hh_pt[i])
                                hdy_hh_eta.SetBinContent(i,dy_hh_eta[i])
                                hdy_hh_phi.SetBinContent(i,dy_hh_phi[i])

			#NOTE: Error should not be symmetric up and down
                                hdy_h1_mass.SetBinError(i,dy_h1_mass[i]**0.5)
                                hdy_h1_pt.SetBinError(i,  dy_h1_pt[i]**0.5)
                                hdy_h1_eta.SetBinError(i, dy_h1_eta[i]**0.5)
                                hdy_h1_phi.SetBinError(i, dy_h1_phi[i]**0.5)

                                hdy_h2_mass.SetBinError(i,dy_h2_mass[i]**0.5)
                                hdy_h2_pt.SetBinError(i,  dy_h2_pt[i]**0.5)
                                hdy_h2_eta.SetBinError(i, dy_h2_eta[i]**0.5)
                                hdy_h2_phi.SetBinError(i, dy_h2_phi[i]**0.5)

                                hdy_hh_mass.SetBinError(i,dy_hh_mass[i]**0.5)
                                hdy_hh_pt.SetBinError(i,  dy_hh_pt[i]**0.5)
                                hdy_hh_eta.SetBinError(i, dy_hh_eta[i]**0.5)
                                hdy_hh_phi.SetBinError(i, dy_hh_phi[i]**0.5)

		# WJ
		sf_wj_jets   = (2234295000/4.01939869276e+13)
		sf_wj_1j     = 342955200/45283121.0
		sf_wj_2j     = 118072500/60438768.0
		sf_wj_3j     = 41779500/59228199.0
		sf_wj_4j     = 23033220/27868119.0
		if (dr == 'wj_jets' or dr == 'wj_1j' or dr == 'wj_2j' or dr == 'wj_3j' or dr == 'wj_4j'):

			if (dr == 'wj_jets'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					wj_j1_mass[i] += j1_mass.GetBinContent(i)*sf_wj_jets
                                        wj_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_wj_jets
                                        wj_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_wj_jets
                                        wj_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_wj_jets
                                        wj_j1_btag[i] += j1_btag.GetBinContent(i)*sf_wj_jets

                                        wj_j2_mass[i] += j2_mass.GetBinContent(i)*sf_wj_jets
                                        wj_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_wj_jets
                                        wj_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_wj_jets
                                        wj_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_wj_jets
                                        wj_j2_btag[i] += j2_btag.GetBinContent(i)*sf_wj_jets

                                        wj_l1_mass[i] += l1_mass.GetBinContent(i)*sf_wj_jets
                                        wj_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_wj_jets
                                        wj_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_wj_jets
                                        wj_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_wj_jets
                                        wj_l1_charge[i] += l1_charge.GetBinContent(i)*sf_wj_jets

                                        wj_l2_mass[i] += l2_mass.GetBinContent(i)*sf_wj_jets
                                        wj_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_wj_jets
                                        wj_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_wj_jets
                                        wj_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_wj_jets
                                        wj_l2_charge[i] += l2_charge.GetBinContent(i)*sf_wj_jets

                                        wj_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_wj_jets
                                        wj_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_wj_jets

					wj_h1_mass[i] += h1_mass.GetBinContent(i)*sf_wj_jets
					wj_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_wj_jets
					wj_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_wj_jets
					wj_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_wj_jets

					wj_h2_mass[i] += h2_mass.GetBinContent(i)*sf_wj_jets
                                        wj_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_wj_jets
                                        wj_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_wj_jets
                                        wj_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_wj_jets

					wj_hh_mass[i] += hh_mass.GetBinContent(i)*sf_wj_jets
                                        wj_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_wj_jets
                                        wj_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_wj_jets
                                        wj_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_wj_jets

			if (dr == 'wj_1j'):
				for i in range(1,h1_mass.GetNbinsX()+1):

					wj_j1_mass[i] += j1_mass.GetBinContent(i)*sf_wj_1j
                                        wj_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_wj_1j
                                        wj_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_wj_1j
                                        wj_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_wj_1j
                                        wj_j1_btag[i] += j1_btag.GetBinContent(i)*sf_wj_1j

                                        wj_j2_mass[i] += j2_mass.GetBinContent(i)*sf_wj_1j
                                        wj_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_wj_1j
                                        wj_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_wj_1j
                                        wj_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_wj_1j
                                        wj_j2_btag[i] += j2_btag.GetBinContent(i)*sf_wj_1j

                                        wj_l1_mass[i] += l1_mass.GetBinContent(i)*sf_wj_1j
                                        wj_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_wj_1j
                                        wj_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_wj_1j
                                        wj_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_wj_1j
                                        wj_l1_charge[i] += l1_charge.GetBinContent(i)*sf_wj_1j

                                        wj_l2_mass[i] += l2_mass.GetBinContent(i)*sf_wj_1j
                                        wj_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_wj_1j
                                        wj_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_wj_1j
                                        wj_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_wj_1j
                                        wj_l2_charge[i] += l2_charge.GetBinContent(i)*sf_wj_1j

                                        wj_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_wj_1j
                                        wj_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_wj_1j

					wj_h1_mass[i] += h1_mass.GetBinContent(i)*sf_wj_1j
                                        wj_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_wj_1j
                                        wj_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_wj_1j
                                        wj_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_wj_1j

                                        wj_h2_mass[i] += h2_mass.GetBinContent(i)*sf_wj_1j
                                        wj_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_wj_1j
                                        wj_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_wj_1j
                                        wj_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_wj_1j

                                        wj_hh_mass[i] += hh_mass.GetBinContent(i)*sf_wj_1j
                                        wj_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_wj_1j
                                        wj_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_wj_1j
                                        wj_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_wj_1j

			if (dr == 'wj_2j'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

					wj_j1_mass[i] += j1_mass.GetBinContent(i)*sf_wj_2j
                                        wj_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_wj_2j
                                        wj_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_wj_2j
                                        wj_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_wj_2j
                                        wj_j1_btag[i] += j1_btag.GetBinContent(i)*sf_wj_2j

                                        wj_j2_mass[i] += j2_mass.GetBinContent(i)*sf_wj_2j
                                        wj_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_wj_2j
                                        wj_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_wj_2j
                                        wj_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_wj_2j
                                        wj_j2_btag[i] += j2_btag.GetBinContent(i)*sf_wj_2j

                                        wj_l1_mass[i] += l1_mass.GetBinContent(i)*sf_wj_2j
                                        wj_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_wj_2j
                                        wj_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_wj_2j
                                        wj_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_wj_2j
                                        wj_l1_charge[i] += l1_charge.GetBinContent(i)*sf_wj_2j

                                        wj_l2_mass[i] += l2_mass.GetBinContent(i)*sf_wj_2j
                                        wj_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_wj_2j
                                        wj_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_wj_2j
                                        wj_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_wj_2j
                                        wj_l2_charge[i] += l2_charge.GetBinContent(i)*sf_wj_2j

                                        wj_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_wj_2j
                                        wj_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_wj_2j

                                        wj_h1_mass[i] += h1_mass.GetBinContent(i)*sf_wj_2j
                                        wj_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_wj_2j
                                        wj_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_wj_2j
                                        wj_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_wj_2j

                                        wj_h2_mass[i] += h2_mass.GetBinContent(i)*sf_wj_2j
                                        wj_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_wj_2j
                                        wj_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_wj_2j
                                        wj_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_wj_2j

                                        wj_hh_mass[i] += hh_mass.GetBinContent(i)*sf_wj_2j
                                        wj_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_wj_2j
                                        wj_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_wj_2j
                                        wj_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_wj_2j

			if (dr == 'wj_3j'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

					wj_j1_mass[i] += j1_mass.GetBinContent(i)*sf_wj_3j
                                        wj_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_wj_3j
                                        wj_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_wj_3j
                                        wj_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_wj_3j
                                        wj_j1_btag[i] += j1_btag.GetBinContent(i)*sf_wj_3j

                                        wj_j2_mass[i] += j2_mass.GetBinContent(i)*sf_wj_3j
                                        wj_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_wj_3j
                                        wj_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_wj_3j
                                        wj_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_wj_3j
                                        wj_j2_btag[i] += j2_btag.GetBinContent(i)*sf_wj_3j

                                        wj_l1_mass[i] += l1_mass.GetBinContent(i)*sf_wj_3j
                                        wj_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_wj_3j
                                        wj_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_wj_3j
                                        wj_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_wj_3j
                                        wj_l1_charge[i] += l1_charge.GetBinContent(i)*sf_wj_3j

                                        wj_l2_mass[i] += l2_mass.GetBinContent(i)*sf_wj_3j
                                        wj_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_wj_3j
                                        wj_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_wj_3j
                                        wj_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_wj_3j
                                        wj_l2_charge[i] += l2_charge.GetBinContent(i)*sf_wj_3j

                                        wj_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_wj_3j
                                        wj_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_wj_3j

                                        wj_h1_mass[i] += h1_mass.GetBinContent(i)*sf_wj_3j
                                        wj_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_wj_3j
                                        wj_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_wj_3j
                                        wj_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_wj_3j

                                        wj_h2_mass[i] += h2_mass.GetBinContent(i)*sf_wj_3j
                                        wj_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_wj_3j
                                        wj_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_wj_3j
                                        wj_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_wj_3j

                                        wj_hh_mass[i] += hh_mass.GetBinContent(i)*sf_wj_3j
                                        wj_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_wj_3j
                                        wj_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_wj_3j
                                        wj_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_wj_3j

			if (dr == 'wj_4j'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

					wj_j1_mass[i] += j1_mass.GetBinContent(i)*sf_wj_4j
                                        wj_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_wj_4j
                                        wj_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_wj_4j
                                        wj_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_wj_4j
                                        wj_j1_btag[i] += j1_btag.GetBinContent(i)*sf_wj_4j

                                        wj_j2_mass[i] += j2_mass.GetBinContent(i)*sf_wj_4j
                                        wj_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_wj_4j
                                        wj_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_wj_4j
                                        wj_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_wj_4j
                                        wj_j2_btag[i] += j2_btag.GetBinContent(i)*sf_wj_4j

                                        wj_l1_mass[i] += l1_mass.GetBinContent(i)*sf_wj_4j
                                        wj_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_wj_4j
                                        wj_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_wj_4j
                                        wj_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_wj_4j
                                        wj_l1_charge[i] += l1_charge.GetBinContent(i)*sf_wj_4j

                                        wj_l2_mass[i] += l2_mass.GetBinContent(i)*sf_wj_4j
                                        wj_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_wj_4j
                                        wj_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_wj_4j
                                        wj_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_wj_4j
                                        wj_l2_charge[i] += l2_charge.GetBinContent(i)*sf_wj_4j

                                        wj_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_wj_4j
                                        wj_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_wj_4j

                                        wj_h1_mass[i] += h1_mass.GetBinContent(i)*sf_wj_4j
                                        wj_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_wj_4j
                                        wj_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_wj_4j
                                        wj_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_wj_4j

                                        wj_h2_mass[i] += h2_mass.GetBinContent(i)*sf_wj_4j
                                        wj_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_wj_4j
                                        wj_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_wj_4j
                                        wj_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_wj_4j

                                        wj_hh_mass[i] += hh_mass.GetBinContent(i)*sf_wj_4j
                                        wj_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_wj_4j
                                        wj_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_wj_4j
                                        wj_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_wj_4j


		# ST
                sf_st_schannel_incl = 122068.8/29561763.7269
                sf_st_tW_anti = 1300614/6933094.0
                sf_st_tW_top  = 1300614/6952830.0
		sf_st_tchannel_anti = 2942730/38811017.0
		sf_st_tchannel_top  = 4940880/58403420.0
                if (dr == 'st_schannel_incl' or dr == "st_tW_anti" or dr == "st_tW_top" or dr == "st_tW_ll" or dr == "st_tchannel_anti" or dr == "st_tchannel_top"):

                        if (dr == 'st_schannel_incl'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

                                        st_j1_mass[i] += j1_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_st_schannel_incl
                                        st_j1_btag[i] += j1_btag.GetBinContent(i)*sf_st_schannel_incl

                                        st_j2_mass[i] += j2_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_st_schannel_incl
                                        st_j2_btag[i] += j2_btag.GetBinContent(i)*sf_st_schannel_incl

                                        st_l1_mass[i] += l1_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_st_schannel_incl
                                        st_l1_charge[i] += l1_charge.GetBinContent(i)*sf_st_schannel_incl

                                        st_l2_mass[i] += l2_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_st_schannel_incl
                                        st_l2_charge[i] += l2_charge.GetBinContent(i)*sf_st_schannel_incl

                                        st_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_st_schannel_incl
                                        st_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_st_schannel_incl

                                        st_h1_mass[i] += h1_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_st_schannel_incl

					st_h2_mass[i] += h2_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_st_schannel_incl

					st_hh_mass[i] += hh_mass.GetBinContent(i)*sf_st_schannel_incl
                                        st_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_st_schannel_incl
                                        st_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_st_schannel_incl
                                        st_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_st_schannel_incl

			if (dr == 'st_tW_anti'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

                                        st_j1_mass[i] += j1_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_st_tW_anti
                                        st_j1_btag[i] += j1_btag.GetBinContent(i)*sf_st_tW_anti

                                        st_j2_mass[i] += j2_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_st_tW_anti
                                        st_j2_btag[i] += j2_btag.GetBinContent(i)*sf_st_tW_anti

                                        st_l1_mass[i] += l1_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_st_tW_anti
                                        st_l1_charge[i] += l1_charge.GetBinContent(i)*sf_st_tW_anti

                                        st_l2_mass[i] += l2_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_st_tW_anti
                                        st_l2_charge[i] += l2_charge.GetBinContent(i)*sf_st_tW_anti

                                        st_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_st_tW_anti
                                        st_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_st_tW_anti

                                        st_h1_mass[i] += h1_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_st_tW_anti

                                        st_h2_mass[i] += h2_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_st_tW_anti

                                        st_hh_mass[i] += hh_mass.GetBinContent(i)*sf_st_tW_anti
                                        st_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_st_tW_anti
                                        st_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_st_tW_anti
                                        st_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_st_tW_anti

			if (dr == 'st_tW_top'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

                                        st_j1_mass[i] += j1_mass.GetBinContent(i)*sf_st_tW_top
                                        st_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_st_tW_top
                                        st_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_st_tW_top
                                        st_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_st_tW_top
                                        st_j1_btag[i] += j1_btag.GetBinContent(i)*sf_st_tW_top

                                        st_j2_mass[i] += j2_mass.GetBinContent(i)*sf_st_tW_top
                                        st_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_st_tW_top
                                        st_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_st_tW_top
                                        st_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_st_tW_top
                                        st_j2_btag[i] += j2_btag.GetBinContent(i)*sf_st_tW_top

                                        st_l1_mass[i] += l1_mass.GetBinContent(i)*sf_st_tW_top
                                        st_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_st_tW_top
                                        st_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_st_tW_top
                                        st_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_st_tW_top
                                        st_l1_charge[i] += l1_charge.GetBinContent(i)*sf_st_tW_top

                                        st_l2_mass[i] += l2_mass.GetBinContent(i)*sf_st_tW_top
                                        st_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_st_tW_top
                                        st_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_st_tW_top
                                        st_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_st_tW_top
                                        st_l2_charge[i] += l2_charge.GetBinContent(i)*sf_st_tW_top

                                        st_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_st_tW_top
                                        st_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_st_tW_top

                                        st_h1_mass[i] += h1_mass.GetBinContent(i)*sf_st_tW_top
                                        st_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_st_tW_top
                                        st_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_st_tW_top
                                        st_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_st_tW_top

                                        st_h2_mass[i] += h2_mass.GetBinContent(i)*sf_st_tW_top
                                        st_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_st_tW_top
                                        st_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_st_tW_top
                                        st_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_st_tW_top

                                        st_hh_mass[i] += hh_mass.GetBinContent(i)*sf_st_tW_top
                                        st_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_st_tW_top
                                        st_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_st_tW_top
                                        st_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_st_tW_top

			if (dr == 'st_tchannel_anti'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

                                        st_j1_mass[i] += j1_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j1_btag[i] += j1_btag.GetBinContent(i)*sf_st_tchannel_anti

                                        st_j2_mass[i] += j2_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_st_tchannel_anti
                                        st_j2_btag[i] += j2_btag.GetBinContent(i)*sf_st_tchannel_anti

                                        st_l1_mass[i] += l1_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l1_charge[i] += l1_charge.GetBinContent(i)*sf_st_tchannel_anti

                                        st_l2_mass[i] += l2_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_st_tchannel_anti
                                        st_l2_charge[i] += l2_charge.GetBinContent(i)*sf_st_tchannel_anti

                                        st_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_st_tchannel_anti
                                        st_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_st_tchannel_anti

                                        st_h1_mass[i] += h1_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_st_tchannel_anti

                                        st_h2_mass[i] += h2_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_st_tchannel_anti

                                        st_hh_mass[i] += hh_mass.GetBinContent(i)*sf_st_tchannel_anti
                                        st_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_st_tchannel_anti
                                        st_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_st_tchannel_anti
                                        st_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_st_tchannel_anti

			if (dr == 'st_tchannel_top'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

                                        st_j1_mass[i] += j1_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_st_tchannel_top
                                        st_j1_btag[i] += j1_btag.GetBinContent(i)*sf_st_tchannel_top

                                        st_j2_mass[i] += j2_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_st_tchannel_top
                                        st_j2_btag[i] += j2_btag.GetBinContent(i)*sf_st_tchannel_top

                                        st_l1_mass[i] += l1_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_st_tchannel_top
                                        st_l1_charge[i] += l1_charge.GetBinContent(i)*sf_st_tchannel_top

                                        st_l2_mass[i] += l2_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_st_tchannel_top
                                        st_l2_charge[i] += l2_charge.GetBinContent(i)*sf_st_tchannel_top

                                        st_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_st_tchannel_top
                                        st_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_st_tchannel_top

                                        st_h1_mass[i] += h1_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_st_tchannel_top

                                        st_h2_mass[i] += h2_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_st_tchannel_top

                                        st_hh_mass[i] += hh_mass.GetBinContent(i)*sf_st_tchannel_top
                                        st_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_st_tchannel_top
                                        st_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_st_tchannel_top
                                        st_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_st_tchannel_top


		for i in range(1,h1_mass.GetNbinsX()+1):

                                hot_j1_mass.SetBinContent(i,wj_j1_mass[i]+st_j1_mass[i])
                                hot_j1_pt.SetBinContent(i,wj_j1_pt[i]+st_j1_pt[i])
                                hot_j1_eta.SetBinContent(i,wj_j1_eta[i]+st_j1_eta[i])
                                hot_j1_phi.SetBinContent(i,wj_j1_phi[i]+st_j1_phi[i])
                                hot_j1_btag.SetBinContent(i,wj_j1_btag[i]+st_j1_btag[i])

				hot_j2_mass.SetBinContent(i,wj_j2_mass[i]+st_j2_mass[i])
                                hot_j2_pt.SetBinContent(i,wj_j2_pt[i]+st_j2_pt[i])
                                hot_j2_eta.SetBinContent(i,wj_j2_eta[i]+st_j2_eta[i])
                                hot_j2_phi.SetBinContent(i,wj_j2_phi[i]+st_j2_phi[i])
                                hot_j2_btag.SetBinContent(i,wj_j2_btag[i]+st_j2_btag[i])

				hot_l1_mass.SetBinContent(i,wj_l1_mass[i]+st_l1_mass[i])
                                hot_l1_pt.SetBinContent(i,wj_l1_pt[i]+st_l1_pt[i])
                                hot_l1_eta.SetBinContent(i,wj_l1_eta[i]+st_l1_eta[i])
                                hot_l1_phi.SetBinContent(i,wj_l1_phi[i]+st_l1_phi[i])
                                hot_l1_charge.SetBinContent(i,wj_l1_charge[i]+st_l1_charge[i])

				hot_l2_mass.SetBinContent(i,wj_l2_mass[i]+st_l2_mass[i])
                                hot_l2_pt.SetBinContent(i,wj_l2_pt[i]+st_l2_pt[i])
                                hot_l2_eta.SetBinContent(i,wj_l2_eta[i]+st_l2_eta[i])
                                hot_l2_phi.SetBinContent(i,wj_l2_phi[i]+st_l2_phi[i])
                                hot_l2_charge.SetBinContent(i,wj_l2_charge[i]+st_l2_charge[i])

                                hot_dr_ll.SetBinContent(i,wj_dr_ll[i]+st_dr_ll[i])
                                hot_dr_jj.SetBinContent(i,wj_dr_jj[i]+st_dr_jj[i])

                                hot_h1_mass.SetBinContent(i,wj_h1_mass[i]+st_h1_mass[i])
                                hot_h1_pt.SetBinContent(i,wj_h1_pt[i]+st_h1_pt[i])
                                hot_h1_eta.SetBinContent(i,wj_h1_eta[i]+st_h1_eta[i])
                                hot_h1_phi.SetBinContent(i,wj_h1_phi[i]+st_h1_phi[i])

				hot_h2_mass.SetBinContent(i,wj_h2_mass[i]+st_h2_mass[i])
                                hot_h2_pt.SetBinContent(i,wj_h2_pt[i]+st_h2_pt[i])
                                hot_h2_eta.SetBinContent(i,wj_h2_eta[i]+st_h2_eta[i])
                                hot_h2_phi.SetBinContent(i,wj_h2_phi[i]+st_h2_phi[i])

				hot_hh_mass.SetBinContent(i,wj_hh_mass[i]+st_hh_mass[i])
                                hot_hh_pt.SetBinContent(i,wj_hh_pt[i]+st_hh_pt[i])
                                hot_hh_eta.SetBinContent(i,wj_hh_eta[i]+st_hh_eta[i])
                                hot_hh_phi.SetBinContent(i,wj_hh_phi[i]+st_hh_phi[i])

			#NOTE: Error should not be symmetric up and down
                                hot_h1_mass.SetBinError(i,hot_h1_mass[i]**0.5)
                                hot_h1_pt.SetBinError(i,  hot_h1_pt[i]**0.5)
                                hot_h1_eta.SetBinError(i, hot_h1_eta[i]**0.5)
                                hot_h1_phi.SetBinError(i, hot_h1_phi[i]**0.5)

                                hot_h2_mass.SetBinError(i,hot_h2_mass[i]**0.5)
                                hot_h2_pt.SetBinError(i,  hot_h2_pt[i]**0.5)
                                hot_h2_eta.SetBinError(i, hot_h2_eta[i]**0.5)
                                hot_h2_phi.SetBinError(i, hot_h2_phi[i]**0.5)

                                hot_hh_mass.SetBinError(i,hot_hh_mass[i]**0.5)
                                hot_hh_pt.SetBinError(i,  hot_hh_pt[i]**0.5)
                                hot_hh_eta.SetBinError(i, hot_hh_eta[i]**0.5)
                                hot_hh_phi.SetBinError(i, hot_hh_phi[i]**0.5)

		f.Close()


	hst_j1_mass = THStack("hst_j1_mass","Stacked Mass Distribution J1")
        hst_j1_pt   = THStack("hst_j1_pt","Stacked Pt Distribution J1")
        hst_j1_eta  = THStack("hst_j1_eta","Stacked Phi Distribution J1")
        hst_j1_phi  = THStack("hst_j1_phi","Stacked Eta Distribution J1")
        hst_j1_btag = THStack("hst_j1_btag","Stacked Btag Distribution J1")

        hst_j2_mass = THStack("hst_j2_mass","Stacked Mass Distribution J2")
        hst_j2_pt   = THStack("hst_j2_pt","Stacked Pt Distribution J2")
        hst_j2_eta  = THStack("hst_j2_eta","Stacked Phi Distribution J2")
        hst_j2_phi  = THStack("hst_j2_phi","Stacked Eta Distribution J2")
        hst_j2_btag = THStack("hst_j2_btag","Stacked Btag Distribution J2")

	hst_l1_mass = THStack("hst_l1_mass","Stacked Mass Distribution L1")
        hst_l1_pt   = THStack("hst_l1_pt","Stacked Pt Distribution L1")
        hst_l1_eta  = THStack("hst_l1_eta","Stacked Phi Distribution L1")
        hst_l1_phi  = THStack("hst_l1_phi","Stacked Eta Distribution L1")
        hst_l1_charge = THStack("hst_l1_charge","Stacked Charge Distribution L1")

        hst_l2_mass = THStack("hst_l2_mass","Stacked Mass Distribution L2")
        hst_l2_pt   = THStack("hst_l2_pt","Stacked Pt Distribution L2")
        hst_l2_eta  = THStack("hst_l2_eta","Stacked Phi Distribution L2")
        hst_l2_phi  = THStack("hst_l2_phi","Stacked Eta Distribution L2")
        hst_l2_charge = THStack("hst_l2_charge","Stacked Charge Distribution L2")

	hst_dr_ll   = THStack("hst_dr_ll","Stacked DeltaR Distribution LL")
	hst_dr_jj   = THStack("hst_dr_jj","Stacked DeltaR Distribution JJ")

	hst_h1_mass = THStack("hst_h1_mass","Stacked Mass Distribution H1")
        hst_h1_pt   = THStack("hst_h1_pt","Stacked Pt Distribution H1")
        hst_h1_eta  = THStack("hst_h1_eta","Stacked Phi Distribution H1")
        hst_h1_phi  = THStack("hst_h1_phi","Stacked Eta Distribution H1")
	
        hst_h2_mass = THStack("hst_h2_mass","Stacked Mass Distribution H2<-bb")
        hst_h2_pt   = THStack("hst_h2_pt","Stacked Pt Distribution H2<-bb")
        hst_h2_eta  = THStack("hst_h2_eta","Stacked Phi Distribution H2<-bb")
        hst_h2_phi  = THStack("hst_h2_phi","Stacked Eta Distribution H2<-bb")

	hst_hh_mass = THStack("hst_hh_mass","Stacked Mass Distribution HH")
	hst_hh_pt   = THStack("hst_hh_pt","Stacked Pt Distribution HH")
	hst_hh_eta  = THStack("hst_hh_eta","Stacked Eta Distribution HH")
	hst_hh_phi  = THStack("hst_hh_phi","Stacked Phi Distribution HH")


	htt_j1_mass.SetFillColor(31)
        hdy_j1_mass.SetFillColor(41)
        hot_j1_mass.SetFillColor(52)
        hst_j1_mass.Add(hot_j1_mass)
        hst_j1_mass.Add(hdy_j1_mass)
        hst_j1_mass.Add(htt_j1_mass)

        htt_j1_pt.SetFillColor(31)
        hdy_j1_pt.SetFillColor(41)
        hot_j1_pt.SetFillColor(52)
        hst_j1_pt.Add(hot_j1_pt)
        hst_j1_pt.Add(hdy_j1_pt)
        hst_j1_pt.Add(htt_j1_pt)

        htt_j1_eta.SetFillColor(31)
        hdy_j1_eta.SetFillColor(41)
        hot_j1_eta.SetFillColor(52)
        hst_j1_eta.Add(hot_j1_eta)
        hst_j1_eta.Add(hdy_j1_eta)
        hst_j1_eta.Add(htt_j1_eta)

        htt_j1_phi.SetFillColor(31)
        hdy_j1_phi.SetFillColor(41)
        hot_j1_phi.SetFillColor(52)
        hst_j1_phi.Add(hot_j1_phi)
        hst_j1_phi.Add(hdy_j1_phi)
        hst_j1_phi.Add(htt_j1_phi)

	htt_j1_btag.SetFillColor(31)
        hdy_j1_btag.SetFillColor(41)
        hot_j1_btag.SetFillColor(52)
        hst_j1_btag.Add(hot_j1_btag)
        hst_j1_btag.Add(hdy_j1_btag)
        hst_j1_btag.Add(htt_j1_btag)


	htt_j2_mass.SetFillColor(31)
        hdy_j2_mass.SetFillColor(41)
        hot_j2_mass.SetFillColor(52)
        hst_j2_mass.Add(hot_j2_mass)
        hst_j2_mass.Add(hdy_j2_mass)
        hst_j2_mass.Add(htt_j2_mass)

        htt_j2_pt.SetFillColor(31)
        hdy_j2_pt.SetFillColor(41)
        hot_j2_pt.SetFillColor(52)
        hst_j2_pt.Add(hot_j2_pt)
        hst_j2_pt.Add(hdy_j2_pt)
        hst_j2_pt.Add(htt_j2_pt)

        htt_j2_eta.SetFillColor(31)
        hdy_j2_eta.SetFillColor(41)
        hot_j2_eta.SetFillColor(52)
        hst_j2_eta.Add(hot_j2_eta)
        hst_j2_eta.Add(hdy_j2_eta)
        hst_j2_eta.Add(htt_j2_eta)

        htt_j2_phi.SetFillColor(31)
        hdy_j2_phi.SetFillColor(41)
        hot_j2_phi.SetFillColor(52)
        hst_j2_phi.Add(hot_j2_phi)
        hst_j2_phi.Add(hdy_j2_phi)
        hst_j2_phi.Add(htt_j2_phi)

	htt_j2_btag.SetFillColor(31)
        hdy_j2_btag.SetFillColor(41)
        hot_j2_btag.SetFillColor(52)
        hst_j2_btag.Add(hot_j2_btag)
        hst_j2_btag.Add(hdy_j2_btag)
        hst_j2_btag.Add(htt_j2_btag)


	htt_l1_mass.SetFillColor(31)
        hdy_l1_mass.SetFillColor(41)
        hot_l1_mass.SetFillColor(52)
        hst_l1_mass.Add(hot_l1_mass)
        hst_l1_mass.Add(hdy_l1_mass)
        hst_l1_mass.Add(htt_l1_mass)

        htt_l1_pt.SetFillColor(31)
        hdy_l1_pt.SetFillColor(41)
        hot_l1_pt.SetFillColor(52)
        hst_l1_pt.Add(hot_l1_pt)
        hst_l1_pt.Add(hdy_l1_pt)
        hst_l1_pt.Add(htt_l1_pt)

        htt_l1_eta.SetFillColor(31)
        hdy_l1_eta.SetFillColor(41)
        hot_l1_eta.SetFillColor(52)
        hst_l1_eta.Add(hot_l1_eta)
        hst_l1_eta.Add(hdy_l1_eta)
        hst_l1_eta.Add(htt_l1_eta)

        htt_l1_phi.SetFillColor(31)
        hdy_l1_phi.SetFillColor(41)
        hot_l1_phi.SetFillColor(52)
        hst_l1_phi.Add(hot_l1_phi)
        hst_l1_phi.Add(hdy_l1_phi)
        hst_l1_phi.Add(htt_l1_phi)

	htt_l1_charge.SetFillColor(31)
        hdy_l1_charge.SetFillColor(41)
        hot_l1_charge.SetFillColor(52)
        hst_l1_charge.Add(hot_l1_charge)
        hst_l1_charge.Add(hdy_l1_charge)
        hst_l1_charge.Add(htt_l1_charge)


	htt_l2_mass.SetFillColor(31)
        hdy_l2_mass.SetFillColor(41)
        hot_l2_mass.SetFillColor(52)
        hst_l2_mass.Add(hot_l2_mass)
        hst_l2_mass.Add(hdy_l2_mass)
        hst_l2_mass.Add(htt_l2_mass)

        htt_l2_pt.SetFillColor(31)
        hdy_l2_pt.SetFillColor(41)
        hot_l2_pt.SetFillColor(52)
        hst_l2_pt.Add(hot_l2_pt)
        hst_l2_pt.Add(hdy_l2_pt)
        hst_l2_pt.Add(htt_l2_pt)

        htt_l2_eta.SetFillColor(31)
        hdy_l2_eta.SetFillColor(41)
        hot_l2_eta.SetFillColor(52)
        hst_l2_eta.Add(hot_l2_eta)
        hst_l2_eta.Add(hdy_l2_eta)
        hst_l2_eta.Add(htt_l2_eta)

        htt_l2_phi.SetFillColor(31)
        hdy_l2_phi.SetFillColor(41)
        hot_l2_phi.SetFillColor(52)
        hst_l2_phi.Add(hot_l2_phi)
        hst_l2_phi.Add(hdy_l2_phi)
        hst_l2_phi.Add(htt_l2_phi)

	htt_l2_charge.SetFillColor(31)
        hdy_l2_charge.SetFillColor(41)
        hot_l2_charge.SetFillColor(52)
        hst_l2_charge.Add(hot_l2_charge)
        hst_l2_charge.Add(hdy_l2_charge)
        hst_l2_charge.Add(htt_l2_charge)


	htt_dr_ll.SetFillColor(31)
	hdy_dr_ll.SetFillColor(41)
	hot_dr_ll.SetFillColor(52)
	hst_dr_ll.Add(hot_dr_ll)
	hst_dr_ll.Add(hdy_dr_ll)
	hst_dr_ll.Add(htt_dr_ll)

	htt_dr_jj.SetFillColor(31)
        hdy_dr_jj.SetFillColor(41)
        hot_dr_jj.SetFillColor(52)
        hst_dr_jj.Add(hot_dr_jj)
        hst_dr_jj.Add(hdy_dr_jj)
        hst_dr_jj.Add(htt_dr_jj)


	htt_h1_mass.SetFillColor(31)
	hdy_h1_mass.SetFillColor(41)
	hot_h1_mass.SetFillColor(52)
	hst_h1_mass.Add(hot_h1_mass)
	hst_h1_mass.Add(hdy_h1_mass)
	hst_h1_mass.Add(htt_h1_mass)

	htt_h1_pt.SetFillColor(31)
        hdy_h1_pt.SetFillColor(41)
        hot_h1_pt.SetFillColor(52)
        hst_h1_pt.Add(hot_h1_pt)
        hst_h1_pt.Add(hdy_h1_pt)
        hst_h1_pt.Add(htt_h1_pt)

	htt_h1_eta.SetFillColor(31)
        hdy_h1_eta.SetFillColor(41)
        hot_h1_eta.SetFillColor(52)
        hst_h1_eta.Add(hot_h1_eta)
        hst_h1_eta.Add(hdy_h1_eta)
        hst_h1_eta.Add(htt_h1_eta)

	htt_h1_phi.SetFillColor(31)
        hdy_h1_phi.SetFillColor(41)
        hot_h1_phi.SetFillColor(52)
        hst_h1_phi.Add(hot_h1_phi)
        hst_h1_phi.Add(hdy_h1_phi)
        hst_h1_phi.Add(htt_h1_phi)


        htt_h2_mass.SetFillColor(31)
        hdy_h2_mass.SetFillColor(41)
	hot_h2_mass.SetFillColor(52)
	hst_h2_mass.Add(hot_h2_mass)
        hst_h2_mass.Add(hdy_h2_mass)
        hst_h2_mass.Add(htt_h2_mass)

        htt_h2_pt.SetFillColor(31)
        hdy_h2_pt.SetFillColor(41)
        hot_h2_pt.SetFillColor(52)
        hst_h2_pt.Add(hot_h2_pt)
        hst_h2_pt.Add(hdy_h2_pt)
        hst_h2_pt.Add(htt_h2_pt)

	htt_h2_eta.SetFillColor(31)
        hdy_h2_eta.SetFillColor(41)
        hot_h2_eta.SetFillColor(52)
        hst_h2_eta.Add(hot_h2_eta)
        hst_h2_eta.Add(hdy_h2_eta)
        hst_h2_eta.Add(htt_h2_eta)

	htt_h2_phi.SetFillColor(31)
        hdy_h2_phi.SetFillColor(41)
        hot_h2_phi.SetFillColor(52)
        hst_h2_phi.Add(hot_h2_phi)
        hst_h2_phi.Add(hdy_h2_phi)
        hst_h2_phi.Add(htt_h2_phi)


	htt_hh_mass.SetFillColor(31)
        hdy_hh_mass.SetFillColor(41)
        hot_hh_mass.SetFillColor(52)
        hst_hh_mass.Add(hot_hh_mass)
        hst_hh_mass.Add(hdy_hh_mass)
        hst_hh_mass.Add(htt_hh_mass)

	htt_hh_pt.SetFillColor(31)
        hdy_hh_pt.SetFillColor(41)
        hot_hh_pt.SetFillColor(52)
        hst_hh_pt.Add(hot_hh_pt)
        hst_hh_pt.Add(hdy_hh_pt)
        hst_hh_pt.Add(htt_hh_pt)

	htt_hh_eta.SetFillColor(31)
        hdy_hh_eta.SetFillColor(41)
        hot_hh_eta.SetFillColor(52)
        hst_hh_eta.Add(hot_hh_eta)
        hst_hh_eta.Add(hdy_hh_eta)
        hst_hh_eta.Add(htt_hh_eta)

	htt_hh_phi.SetFillColor(31)
        hdy_hh_phi.SetFillColor(41)
        hot_hh_phi.SetFillColor(52)
        hst_hh_phi.Add(hot_hh_phi)
        hst_hh_phi.Add(hdy_hh_phi)
        hst_hh_phi.Add(htt_hh_phi)

	for i in range(1,htot_j1_mass.GetNbinsX()+1):

		htot_j1_mass.SetBinContent(i,htt_j1_mass.GetBinContent(i)+hdy_j1_mass.GetBinContent(i)+hot_j1_mass.GetBinContent(i))
                htot_j1_pt.SetBinContent(i,htt_j1_pt.GetBinContent(i)+hdy_j1_pt.GetBinContent(i)+hot_j1_pt.GetBinContent(i))
                htot_j1_eta.SetBinContent(i,htt_j1_eta.GetBinContent(i)+hdy_j1_eta.GetBinContent(i)+hot_j1_eta.GetBinContent(i))
                htot_j1_phi.SetBinContent(i,htt_j1_phi.GetBinContent(i)+hdy_j1_phi.GetBinContent(i)+hot_j1_phi.GetBinContent(i))
                htot_j1_btag.SetBinContent(i,htt_j1_btag.GetBinContent(i)+hdy_j1_btag.GetBinContent(i)+hot_j1_btag.GetBinContent(i))

                htot_j2_mass.SetBinContent(i,htt_j2_mass.GetBinContent(i)+hdy_j2_mass.GetBinContent(i)+hot_j2_mass.GetBinContent(i))
                htot_j2_pt.SetBinContent(i,htt_j2_pt.GetBinContent(i)+hdy_j2_pt.GetBinContent(i)+hot_j2_pt.GetBinContent(i))
                htot_j2_eta.SetBinContent(i,htt_j2_eta.GetBinContent(i)+hdy_j2_eta.GetBinContent(i)+hot_j2_eta.GetBinContent(i))
                htot_j2_phi.SetBinContent(i,htt_j2_phi.GetBinContent(i)+hdy_j2_phi.GetBinContent(i)+hot_j2_phi.GetBinContent(i))
                htot_j2_btag.SetBinContent(i,htt_j2_btag.GetBinContent(i)+hdy_j2_btag.GetBinContent(i)+hot_j2_btag.GetBinContent(i))

                htot_l1_mass.SetBinContent(i,htt_l1_mass.GetBinContent(i)+hdy_l1_mass.GetBinContent(i)+hot_l1_mass.GetBinContent(i))
                htot_l1_pt.SetBinContent(i,htt_l1_pt.GetBinContent(i)+hdy_l1_pt.GetBinContent(i)+hot_l1_pt.GetBinContent(i))
                htot_l1_eta.SetBinContent(i,htt_l1_eta.GetBinContent(i)+hdy_l1_eta.GetBinContent(i)+hot_l1_eta.GetBinContent(i))
                htot_l1_phi.SetBinContent(i,htt_l1_phi.GetBinContent(i)+hdy_l1_phi.GetBinContent(i)+hot_l1_phi.GetBinContent(i))
                htot_l1_charge.SetBinContent(i,htt_l1_charge.GetBinContent(i)+hdy_l1_charge.GetBinContent(i)+hot_l1_charge.GetBinContent(i))

                htot_l2_mass.SetBinContent(i,htt_l2_mass.GetBinContent(i)+hdy_l2_mass.GetBinContent(i)+hot_l2_mass.GetBinContent(i))
                htot_l2_pt.SetBinContent(i,htt_l2_pt.GetBinContent(i)+hdy_l2_pt.GetBinContent(i)+hot_l2_pt.GetBinContent(i))
                htot_l2_eta.SetBinContent(i,htt_l2_eta.GetBinContent(i)+hdy_l2_eta.GetBinContent(i)+hot_l2_eta.GetBinContent(i))
                htot_l2_phi.SetBinContent(i,htt_l2_phi.GetBinContent(i)+hdy_l2_phi.GetBinContent(i)+hot_l2_phi.GetBinContent(i))
                htot_l2_charge.SetBinContent(i,htt_l2_charge.GetBinContent(i)+hdy_l2_charge.GetBinContent(i)+hot_l2_charge.GetBinContent(i))

                htot_dr_ll.SetBinContent(i,htt_dr_ll.GetBinContent(i)+hdy_dr_ll.GetBinContent(i)+hot_dr_ll.GetBinContent(i))
                htot_dr_jj.SetBinContent(i,htt_dr_jj.GetBinContent(i)+hdy_dr_jj.GetBinContent(i)+hot_dr_jj.GetBinContent(i))

                htot_h1_mass.SetBinContent(i,htt_h1_mass.GetBinContent(i)+hdy_h1_mass.GetBinContent(i)+hot_h1_mass.GetBinContent(i))
                htot_h1_pt.SetBinContent(i,htt_h1_pt.GetBinContent(i)+hdy_h1_pt.GetBinContent(i)+hot_h1_pt.GetBinContent(i))
                htot_h1_eta.SetBinContent(i,htt_h1_eta.GetBinContent(i)+hdy_h1_eta.GetBinContent(i)+hot_h1_eta.GetBinContent(i))
                htot_h1_phi.SetBinContent(i,htt_h1_phi.GetBinContent(i)+hdy_h1_phi.GetBinContent(i)+hot_h1_phi.GetBinContent(i))

                htot_h2_mass.SetBinContent(i,htt_h2_mass.GetBinContent(i)+hdy_h2_mass.GetBinContent(i)+hot_h2_mass.GetBinContent(i))
                htot_h2_pt.SetBinContent(i,htt_h2_pt.GetBinContent(i)+hdy_h2_pt.GetBinContent(i)+hot_h2_pt.GetBinContent(i))
                htot_h2_eta.SetBinContent(i,htt_h2_eta.GetBinContent(i)+hdy_h2_eta.GetBinContent(i)+hot_h2_eta.GetBinContent(i))
                htot_h2_phi.SetBinContent(i,htt_h2_phi.GetBinContent(i)+hdy_h2_phi.GetBinContent(i)+hot_h2_phi.GetBinContent(i))

                htot_hh_mass.SetBinContent(i,htt_hh_mass.GetBinContent(i)+hdy_hh_mass.GetBinContent(i)+hot_hh_mass.GetBinContent(i))
                htot_hh_pt.SetBinContent(i,htt_hh_pt.GetBinContent(i)+hdy_hh_pt.GetBinContent(i)+hot_hh_pt.GetBinContent(i))
                htot_hh_eta.SetBinContent(i,htt_hh_eta.GetBinContent(i)+hdy_hh_eta.GetBinContent(i)+hot_hh_eta.GetBinContent(i))
                htot_hh_phi.SetBinContent(i,htt_hh_phi.GetBinContent(i)+hdy_hh_phi.GetBinContent(i)+hot_hh_phi.GetBinContent(i))


	for i in range(1,htot_j1_mass.GetNbinsX()+1):

		htot_j1_mass.SetBinError(i,(htt_j1_mass.GetBinError(i)**2+hdy_j1_mass.GetBinError(i)**2+hot_j1_mass.GetBinError(i)**2)**.5)
		htot_j1_pt.SetBinError(i,(htt_j1_pt.GetBinError(i)**2+hdy_j1_pt.GetBinError(i)**2+hot_j1_pt.GetBinError(i)**2)**.5)
		htot_j1_eta.SetBinError(i,(htt_j1_eta.GetBinError(i)**2+hdy_j1_eta.GetBinError(i)**2+hot_j1_eta.GetBinError(i)**2)**.5)
		htot_j1_phi.SetBinError(i,(htt_j1_phi.GetBinError(i)**2+hdy_j1_phi.GetBinError(i)**2+hot_j1_phi.GetBinError(i)**2)**.5)
		htot_j1_btag.SetBinError(i,(htt_j1_btag.GetBinError(i)**2+hdy_j1_btag.GetBinError(i)**2+hot_j1_btag.GetBinError(i)**2)**.5)

		htot_j2_mass.SetBinError(i,(htt_j2_mass.GetBinError(i)**2+hdy_j2_mass.GetBinError(i)**2+hot_j2_mass.GetBinError(i)**2)**.5)
                htot_j2_pt.SetBinError(i,(htt_j2_pt.GetBinError(i)**2+hdy_j2_pt.GetBinError(i)**2+hot_j2_pt.GetBinError(i)**2)**.5)
                htot_j2_eta.SetBinError(i,(htt_j2_eta.GetBinError(i)**2+hdy_j2_eta.GetBinError(i)**2+hot_j2_eta.GetBinError(i)**2)**.5)
                htot_j2_phi.SetBinError(i,(htt_j2_phi.GetBinError(i)**2+hdy_j2_phi.GetBinError(i)**2+hot_j2_phi.GetBinError(i)**2)**.5)
                htot_j2_btag.SetBinError(i,(htt_j2_btag.GetBinError(i)**2+hdy_j2_btag.GetBinError(i)**2+hot_j2_btag.GetBinError(i)**2)**.5)

		htot_l1_mass.SetBinError(i,(htt_l1_mass.GetBinError(i)**2+hdy_l1_mass.GetBinError(i)**2+hot_l1_mass.GetBinError(i)**2)**.5)
                htot_l1_pt.SetBinError(i,(htt_l1_pt.GetBinError(i)**2+hdy_l1_pt.GetBinError(i)**2+hot_l1_pt.GetBinError(i)**2)**.5)
                htot_l1_eta.SetBinError(i,(htt_l1_eta.GetBinError(i)**2+hdy_l1_eta.GetBinError(i)**2+hot_l1_eta.GetBinError(i)**2)**.5)
                htot_l1_phi.SetBinError(i,(htt_l1_phi.GetBinError(i)**2+hdy_l1_phi.GetBinError(i)**2+hot_l1_phi.GetBinError(i)**2)**.5)
                htot_l1_charge.SetBinError(i,(htt_l1_charge.GetBinError(i)**2+hdy_l1_charge.GetBinError(i)**2+hot_l1_charge.GetBinError(i)**2)**.5)

		htot_l2_mass.SetBinError(i,(htt_l2_mass.GetBinError(i)**2+hdy_l2_mass.GetBinError(i)**2+hot_l2_mass.GetBinError(i)**2)**.5)
                htot_l2_pt.SetBinError(i,(htt_l2_pt.GetBinError(i)**2+hdy_l2_pt.GetBinError(i)**2+hot_l2_pt.GetBinError(i)**2)**.5)
                htot_l2_eta.SetBinError(i,(htt_l2_eta.GetBinError(i)**2+hdy_l2_eta.GetBinError(i)**2+hot_l2_eta.GetBinError(i)**2)**.5)
                htot_l2_phi.SetBinError(i,(htt_l2_phi.GetBinError(i)**2+hdy_l2_phi.GetBinError(i)**2+hot_l2_phi.GetBinError(i)**2)**.5)
                htot_l2_charge.SetBinError(i,(htt_l2_charge.GetBinError(i)**2+hdy_l2_charge.GetBinError(i)**2+hot_l2_charge.GetBinError(i)**2)**.5)

		htot_dr_ll.SetBinError(i,(htt_dr_ll.GetBinError(i)**2+hdy_dr_ll.GetBinError(i)**2+hot_dr_ll.GetBinError(i)**2)**.5)
		htot_dr_jj.SetBinError(i,(htt_dr_jj.GetBinError(i)**2+hdy_dr_jj.GetBinError(i)**2+hot_dr_jj.GetBinError(i)**2)**.5)

		htot_h1_mass.SetBinError(i,(htt_h1_mass.GetBinError(i)**2+hdy_h1_mass.GetBinError(i)**2+hot_h1_mass.GetBinError(i)**2)**.5)
                htot_h1_pt.SetBinError(i,(htt_h1_pt.GetBinError(i)**2+hdy_h1_pt.GetBinError(i)**2+hot_h1_pt.GetBinError(i)**2)**.5)
                htot_h1_eta.SetBinError(i,(htt_h1_eta.GetBinError(i)**2+hdy_h1_eta.GetBinError(i)**2+hot_h1_eta.GetBinError(i)**2)**.5)
                htot_h1_phi.SetBinError(i,(htt_h1_phi.GetBinError(i)**2+hdy_h1_phi.GetBinError(i)**2+hot_h1_phi.GetBinError(i)**2)**.5)

		htot_h2_mass.SetBinError(i,(htt_h2_mass.GetBinError(i)**2+hdy_h2_mass.GetBinError(i)**2+hot_h2_mass.GetBinError(i)**2)**.5)
                htot_h2_pt.SetBinError(i,(htt_h2_pt.GetBinError(i)**2+hdy_h2_pt.GetBinError(i)**2+hot_h2_pt.GetBinError(i)**2)**.5)
                htot_h2_eta.SetBinError(i,(htt_h2_eta.GetBinError(i)**2+hdy_h2_eta.GetBinError(i)**2+hot_h2_eta.GetBinError(i)**2)**.5)
                htot_h2_phi.SetBinError(i,(htt_h2_phi.GetBinError(i)**2+hdy_h2_phi.GetBinError(i)**2+hot_h2_phi.GetBinError(i)**2)**.5)

		htot_hh_mass.SetBinError(i,(htt_hh_mass.GetBinError(i)**2+hdy_hh_mass.GetBinError(i)**2+hot_hh_mass.GetBinError(i)**2)**.5)
                htot_hh_pt.SetBinError(i,(htt_hh_pt.GetBinError(i)**2+hdy_hh_pt.GetBinError(i)**2+hot_hh_pt.GetBinError(i)**2)**.5)
                htot_hh_eta.SetBinError(i,(htt_hh_eta.GetBinError(i)**2+hdy_hh_eta.GetBinError(i)**2+hot_hh_eta.GetBinError(i)**2)**.5)
                htot_hh_phi.SetBinError(i,(htt_hh_phi.GetBinError(i)**2+hdy_hh_phi.GetBinError(i)**2+hot_hh_phi.GetBinError(i)**2)**.5)


	lst_j1_mass = TLegend(0.75,0.75,0.95,0.90)
        lst_j1_mass.AddEntry(htt_j1_mass,"TT","f")
        lst_j1_mass.AddEntry(hdy_j1_mass,"DY","f")
        lst_j1_mass.AddEntry(hot_j1_mass,"Others","f")
        lst_j1_mass.AddEntry(hsig_j1_mass,"Signal","l")
        lst_j1_mass.SetTextAlign(13)

        lst_j1_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_j1_pt.AddEntry(htt_j1_pt,"TT","f")
        lst_j1_pt.AddEntry(hdy_j1_pt,"DY","f")
        lst_j1_pt.AddEntry(hot_j1_pt,"Others","f")
        lst_j1_pt.AddEntry(hsig_j1_pt,"Signal","l")
        lst_j1_pt.SetTextAlign(13)

        lst_j1_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_j1_eta.AddEntry(htt_j1_eta,"TT","f")
        lst_j1_eta.AddEntry(hdy_j1_eta,"DY","f")
        lst_j1_eta.AddEntry(hot_j1_eta,"Others","f")
        lst_j1_eta.AddEntry(hsig_j1_eta,"Signal","l")
        lst_j1_eta.SetTextAlign(13)

        lst_j1_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_j1_phi.AddEntry(htt_j1_phi,"TT","f")
        lst_j1_phi.AddEntry(hdy_j1_phi,"DY","f")
        lst_j1_phi.AddEntry(hot_j1_phi,"Others","f")
        lst_j1_phi.AddEntry(hsig_j1_phi,"Signal","l")
        lst_j1_phi.SetTextAlign(13)

	lst_j1_btag = TLegend(0.75,0.75,0.95,0.90)
        lst_j1_btag.AddEntry(htt_j1_btag,"TT","f")
        lst_j1_btag.AddEntry(hdy_j1_btag,"DY","f")
        lst_j1_btag.AddEntry(hot_j1_btag,"Others","f")
        lst_j1_btag.AddEntry(hsig_j1_btag,"Signal","l")
        lst_j1_btag.SetTextAlign(13)


	lst_j2_mass = TLegend(0.75,0.75,0.95,0.90)
        lst_j2_mass.AddEntry(htt_j2_mass,"TT","f")
        lst_j2_mass.AddEntry(hdy_j2_mass,"DY","f")
        lst_j2_mass.AddEntry(hot_j2_mass,"Others","f")
        lst_j2_mass.AddEntry(hsig_j2_mass,"Signal","l")
        lst_j2_mass.SetTextAlign(13)

        lst_j2_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_j2_pt.AddEntry(htt_j2_pt,"TT","f")
        lst_j2_pt.AddEntry(hdy_j2_pt,"DY","f")
        lst_j2_pt.AddEntry(hot_j2_pt,"Others","f")
        lst_j2_pt.AddEntry(hsig_j2_pt,"Signal","l")
        lst_j2_pt.SetTextAlign(13)

        lst_j2_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_j2_eta.AddEntry(htt_j2_eta,"TT","f")
        lst_j2_eta.AddEntry(hdy_j2_eta,"DY","f")
        lst_j2_eta.AddEntry(hot_j2_eta,"Others","f")
        lst_j2_eta.AddEntry(hsig_j2_eta,"Signal","l")
        lst_j2_eta.SetTextAlign(13)

        lst_j2_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_j2_phi.AddEntry(htt_j2_phi,"TT","f")
        lst_j2_phi.AddEntry(hdy_j2_phi,"DY","f")
        lst_j2_phi.AddEntry(hot_j2_phi,"Others","f")
        lst_j2_phi.AddEntry(hsig_j2_phi,"Signal","l")
        lst_j2_phi.SetTextAlign(13)

        lst_j2_btag = TLegend(0.75,0.75,0.95,0.90)
        lst_j2_btag.AddEntry(htt_j2_btag,"TT","f")
        lst_j2_btag.AddEntry(hdy_j2_btag,"DY","f")
        lst_j2_btag.AddEntry(hot_j2_btag,"Others","f")
        lst_j2_btag.AddEntry(hsig_j2_btag,"Signal","l")
        lst_j2_btag.SetTextAlign(13)


	lst_l1_mass = TLegend(0.75,0.75,0.95,0.90)
        lst_l1_mass.AddEntry(htt_l1_mass,"TT","f")
        lst_l1_mass.AddEntry(hdy_l1_mass,"DY","f")
        lst_l1_mass.AddEntry(hot_l1_mass,"Others","f")
        lst_l1_mass.AddEntry(hsig_l1_mass,"Signal","l")
        lst_l1_mass.SetTextAlign(13)

        lst_l1_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_l1_pt.AddEntry(htt_l1_pt,"TT","f")
        lst_l1_pt.AddEntry(hdy_l1_pt,"DY","f")
        lst_l1_pt.AddEntry(hot_l1_pt,"Others","f")
        lst_l1_pt.AddEntry(hsig_l1_pt,"Signal","l")
        lst_l1_pt.SetTextAlign(13)

        lst_l1_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_l1_eta.AddEntry(htt_l1_eta,"TT","f")
        lst_l1_eta.AddEntry(hdy_l1_eta,"DY","f")
        lst_l1_eta.AddEntry(hot_l1_eta,"Others","f")
        lst_l1_eta.AddEntry(hsig_l1_eta,"Signal","l")
        lst_l1_eta.SetTextAlign(13)

        lst_l1_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_l1_phi.AddEntry(htt_l1_phi,"TT","f")
        lst_l1_phi.AddEntry(hdy_l1_phi,"DY","f")
        lst_l1_phi.AddEntry(hot_l1_phi,"Others","f")
        lst_l1_phi.AddEntry(hsig_l1_phi,"Signal","l")
        lst_l1_phi.SetTextAlign(13)

        lst_l1_charge = TLegend(0.75,0.75,0.95,0.90)
        lst_l1_charge.AddEntry(htt_l1_charge,"TT","f")
        lst_l1_charge.AddEntry(hdy_l1_charge,"DY","f")
        lst_l1_charge.AddEntry(hot_l1_charge,"Others","f")
        lst_l1_charge.AddEntry(hsig_l1_charge,"Signal","l")
        lst_l1_charge.SetTextAlign(13)


	lst_l2_mass = TLegend(0.75,0.75,0.95,0.90)
        lst_l2_mass.AddEntry(htt_l2_mass,"TT","f")
        lst_l2_mass.AddEntry(hdy_l2_mass,"DY","f")
        lst_l2_mass.AddEntry(hot_l2_mass,"Others","f")
        lst_l2_mass.AddEntry(hsig_l2_mass,"Signal","l")
        lst_l2_mass.SetTextAlign(13)

        lst_l2_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_l2_pt.AddEntry(htt_l2_pt,"TT","f")
        lst_l2_pt.AddEntry(hdy_l2_pt,"DY","f")
        lst_l2_pt.AddEntry(hot_l2_pt,"Others","f")
        lst_l2_pt.AddEntry(hsig_l2_pt,"Signal","l")
        lst_l2_pt.SetTextAlign(13)

        lst_l2_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_l2_eta.AddEntry(htt_l2_eta,"TT","f")
        lst_l2_eta.AddEntry(hdy_l2_eta,"DY","f")
        lst_l2_eta.AddEntry(hot_l2_eta,"Others","f")
        lst_l2_eta.AddEntry(hsig_l2_eta,"Signal","l")
        lst_l2_eta.SetTextAlign(13)

        lst_l2_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_l2_phi.AddEntry(htt_l2_phi,"TT","f")
        lst_l2_phi.AddEntry(hdy_l2_phi,"DY","f")
        lst_l2_phi.AddEntry(hot_l2_phi,"Others","f")
        lst_l2_phi.AddEntry(hsig_l2_phi,"Signal","l")
        lst_l2_phi.SetTextAlign(13)

        lst_l2_charge = TLegend(0.75,0.75,0.95,0.90)
        lst_l2_charge.AddEntry(htt_l2_charge,"TT","f")
        lst_l2_charge.AddEntry(hdy_l2_charge,"DY","f")
        lst_l2_charge.AddEntry(hot_l2_charge,"Others","f")
        lst_l2_charge.AddEntry(hsig_l2_charge,"Signal","l")
        lst_l2_charge.SetTextAlign(13)


	lst_dr_ll = TLegend(0.75,0.75,0.95,0.90)
	lst_dr_ll.AddEntry(htt_dr_ll,"TT","f")
	lst_dr_ll.AddEntry(hdy_dr_ll,"DY","f")
	lst_dr_ll.AddEntry(hot_dr_ll,"Others","f")
	lst_dr_ll.AddEntry(hsig_dr_ll,"Signal","l")
	lst_dr_ll.SetTextAlign(13)

	lst_dr_jj = TLegend(0.75,0.75,0.95,0.90)
        lst_dr_jj.AddEntry(htt_dr_jj,"TT","f")
        lst_dr_jj.AddEntry(hdy_dr_jj,"DY","f")
        lst_dr_jj.AddEntry(hot_dr_jj,"Others","f")
        lst_dr_jj.AddEntry(hsig_dr_jj,"Signal","l")
        lst_dr_jj.SetTextAlign(13)


	lst_h1_mass = TLegend(0.75,0.75,0.95,0.90)
	lst_h1_mass.AddEntry(htt_h1_mass,"TT","f")
	lst_h1_mass.AddEntry(hdy_h1_mass,"DY","f")
	lst_h1_mass.AddEntry(hot_h1_mass,"Others","f")
	lst_h1_mass.AddEntry(hsig_h1_mass,"Signal","l")
	lst_h1_mass.SetTextAlign(13)

	lst_h1_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_h1_pt.AddEntry(htt_h1_pt,"TT","f")
        lst_h1_pt.AddEntry(hdy_h1_pt,"DY","f")
        lst_h1_pt.AddEntry(hot_h1_pt,"Others","f")
        lst_h1_pt.AddEntry(hsig_h1_pt,"Signal","l")
        lst_h1_pt.SetTextAlign(13)

	lst_h1_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_h1_eta.AddEntry(htt_h1_eta,"TT","f")
        lst_h1_eta.AddEntry(hdy_h1_eta,"DY","f")
        lst_h1_eta.AddEntry(hot_h1_eta,"Others","f")
        lst_h1_eta.AddEntry(hsig_h1_eta,"Signal","l")
        lst_h1_eta.SetTextAlign(13)

	lst_h1_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_h1_phi.AddEntry(htt_h1_phi,"TT","f")
        lst_h1_phi.AddEntry(hdy_h1_phi,"DY","f")
        lst_h1_phi.AddEntry(hot_h1_phi,"Others","f")
        lst_h1_phi.AddEntry(hsig_h1_phi,"Signal","l")
        lst_h1_phi.SetTextAlign(13)


	lst_h2_mass = TLegend(0.75,0.75,0.95,0.90)
        lst_h2_mass.AddEntry(htt_h2_mass,"TT","f")
        lst_h2_mass.AddEntry(hdy_h2_mass,"DY","f")
        lst_h2_mass.AddEntry(hot_h2_mass,"Others","f")
        lst_h2_mass.AddEntry(hsig_h2_mass,"Signal","l")
	lst_h2_mass.SetTextAlign(13)

	lst_h2_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_h2_pt.AddEntry(htt_h2_pt,"TT","f")
        lst_h2_pt.AddEntry(hdy_h2_pt,"DY","f")
        lst_h2_pt.AddEntry(hot_h2_pt,"Others","f")
        lst_h2_pt.AddEntry(hsig_h2_pt,"Signal","l")
        lst_h2_pt.SetTextAlign(13)

	lst_h2_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_h2_eta.AddEntry(htt_h2_eta,"TT","f")
        lst_h2_eta.AddEntry(hdy_h2_eta,"DY","f")
        lst_h2_eta.AddEntry(hot_h2_eta,"Others","f")
        lst_h2_eta.AddEntry(hsig_h2_eta,"Signal","l")
        lst_h2_eta.SetTextAlign(13)

	lst_h2_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_h2_phi.AddEntry(htt_h2_phi,"TT","f")
        lst_h2_phi.AddEntry(hdy_h2_phi,"DY","f")
        lst_h2_phi.AddEntry(hot_h2_phi,"Others","f")
        lst_h2_phi.AddEntry(hsig_h2_phi,"Signal","l")
        lst_h2_phi.SetTextAlign(13)


	lst_hh_mass = TLegend(0.75,0.75,0.95,0.90)
	lst_hh_mass.AddEntry(htt_hh_mass,"TT","f")
	lst_hh_mass.AddEntry(hdy_hh_mass,"DY","f")
	lst_hh_mass.AddEntry(hot_hh_mass,"Others","f")
	lst_hh_mass.AddEntry(hsig_hh_mass,"Signal","l")
	lst_hh_mass.SetTextAlign(13)

	lst_hh_pt = TLegend(0.75,0.75,0.95,0.90)
        lst_hh_pt.AddEntry(htt_hh_pt,"TT","f")
        lst_hh_pt.AddEntry(hdy_hh_pt,"DY","f")
        lst_hh_pt.AddEntry(hot_hh_pt,"Others","f")
        lst_hh_pt.AddEntry(hsig_hh_pt,"Signal","l")
        lst_hh_pt.SetTextAlign(13)

	lst_hh_eta = TLegend(0.75,0.75,0.95,0.90)
        lst_hh_eta.AddEntry(htt_hh_eta,"TT","f")
        lst_hh_eta.AddEntry(hdy_hh_eta,"DY","f")
        lst_hh_eta.AddEntry(hot_hh_eta,"Others","f")
        lst_hh_eta.AddEntry(hsig_hh_eta,"Signal","l")
        lst_hh_eta.SetTextAlign(13)

	lst_hh_phi = TLegend(0.75,0.75,0.95,0.90)
        lst_hh_phi.AddEntry(htt_hh_phi,"TT","f")
        lst_hh_phi.AddEntry(hdy_hh_phi,"DY","f")
        lst_hh_phi.AddEntry(hot_hh_phi,"Others","f")
        lst_hh_phi.AddEntry(hsig_hh_phi,"Signal","l")
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


	c1 = TCanvas("c1","Stacked Histograms of Basic Kinematics h1", 900,600)

	c1.Divide(2,2)
        c1.cd(1)
        hst_h1_mass.Draw("hist")
	hst_h1_mass.GetXaxis().SetTitle("Invariant Mass[GeV]")
	hst_h1_mass.GetYaxis().SetTitle("Events/5GeV")
	hsig_h1_mass.SetLineColor(kBlack)
	hsig_h1_mass.Draw("hist same")
	htot_h1_mass.SetFillColor(kBlue)
	htot_h1_mass.SetFillStyle(3004)
	htot_h1_mass.Draw("E2 same")
	lst_h1_mass.Draw()

	c1.cd(2)
	hst_h1_pt.Draw("hist")
        hst_h1_pt.GetXaxis().SetTitle("Transverse Momentum[GeV]")
        hst_h1_pt.GetYaxis().SetTitle("Events/10GeV")
	hsig_h1_pt.SetLineColor(kBlack)
        hsig_h1_pt.Draw("hist same")
	htot_h1_pt.SetFillColor(kBlue)
	htot_h1_pt.SetFillStyle(3004)
	htot_h1_pt.Draw("E2 same")
	lst_h1_pt.Draw()

	c1.cd(3)
	hst_h1_eta.Draw("hist")
        hst_h1_eta.GetXaxis().SetTitle("Eta")
        hst_h1_eta.GetYaxis().SetTitle("Events/.17")
	hsig_h1_eta.SetLineColor(kBlack)
        hsig_h1_eta.Draw("hist same")
	htot_h1_eta.SetFillColor(kBlue)
	htot_h1_eta.SetFillStyle(3004)
	htot_h1_eta.Draw("E2 same")
	lst_h1_eta.Draw()

	c1.cd(4)
	hst_h1_phi.Draw("hist")
        hst_h1_phi.GetXaxis().SetTitle("Phi")
        hst_h1_phi.GetYaxis().SetTitle("Events/.157deg")
	hsig_h1_phi.SetLineColor(kBlack)
        hsig_h1_phi.Draw("hist same")
	htot_h1_phi.SetFillColor(kBlue)
	htot_h1_phi.SetFillStyle(3004)
	htot_h1_phi.Draw("E2 same")
	lst_h1_phi.Draw()

        c1.SetGrid()
	c1.Draw()
	c1.SaveAs("h1.png")

	input("Press Enter to continue...")


if __name__ == '__main__':
        main()
