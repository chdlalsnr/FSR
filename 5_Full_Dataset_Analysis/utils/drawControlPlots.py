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

        hreal_h1_mass = TH1F("hreal_h1_mass","Real H1 Mass Distribution", 40,0,400)
        hreal_h1_pt   = TH1F("hreal_h1_pt","Real H1 Pt Distribution", 40,0,400)
        hreal_h1_eta  = TH1F("hreal_h1_eta","Real H1 Eta Distribution", 40,-3.4,3.4)
        hreal_h1_phi  = TH1F("hreal_h1_phi","Real H1 Phi Distribution", 40,-3.14,3.14)

        hreal_h2_mass = TH1F("hreal_h2_mass","Real H2 Mass Distribution", 40,0,400)
        hreal_h2_pt   = TH1F("hreal_h2_pt","Real H2 Pt Distribution", 40,0,400)
        hreal_h2_eta  = TH1F("hreal_h2_eta","Real H2 Eta Distribution", 40,-3.4,3.4)
        hreal_h2_phi  = TH1F("hreal_h2_phi","Real H2 Phi Distribution", 40,-3.14,3.14)

        hreal_hh_mass = TH1F("hreal_hh_mass","Real HH Mass Distribution", 40,0,800)
        hreal_hh_pt   = TH1F("hreal_hh_pt","Real HH Pt Distribution", 40,0,600)
        hreal_hh_eta  = TH1F("hreal_hh_eta","Real HH Eta Distribution", 40,-3.4,3.4)
        hreal_hh_phi  = TH1F("hreal_hh_phi","Real HH Phi Distribution", 40,-3.14,3.14)


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

	hsig_h1_mass = TH1F("hsig_h1_mass","Signal H1 Mass Distribution", 40,0,400)
	hsig_h1_pt   = TH1F("hsig_h1_pt","Signal H1 Pt Distribution", 40,0,400)
	hsig_h1_eta  = TH1F("hsig_h1_eta","Signal H1 Eta Distribution", 40,-3.4,3.4)
	hsig_h1_phi  = TH1F("hsig_h1_phi","Signal H1 Phi Distribution", 40,-3.14,3.14)

	hsig_h2_mass = TH1F("hsig_h2_mass","Signal H2 Mass Distribution", 40,0,400)
	hsig_h2_pt   = TH1F("hsig_h2_pt","Signal H2 Pt Distribution", 40,0,400)
	hsig_h2_eta  = TH1F("hsig_h2_eta","Signal H2 Eta Distribution", 40,-3.4,3.4)
	hsig_h2_phi  = TH1F("hsig_h2_phi","Signal H2 Phi Distribution", 40,-3.14,3.14)

	hsig_hh_mass = TH1F("hsig_hh_mass","Signal HH Mass Distribution", 40,0,800)
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

	htt_h1_mass = TH1F("htt_h1_mass","TTbar H1 Mass Distribution", 40,0,400)
	htt_h1_pt   = TH1F("htt_h1_pt","TTbar H1 Pt Distribution", 40,0,400)
	htt_h1_eta  = TH1F("htt_h1_eta","TTbar H1 Eta Distribution", 40,-3.4,3.4)
	htt_h1_phi  = TH1F("htt_h1_phi","TTbar H1 Phi Distribution", 40,-3.14,3.14)

	htt_h2_mass = TH1F("htt_h2_mass","TTbar H2 Mass Distribution", 40,0,400)
	htt_h2_pt   = TH1F("htt_h2_pt","TTbar H2 Pt Distribution", 40,0,400)
	htt_h2_eta  = TH1F("htt_h2_eta","TTbar H2 Eta Distribution", 40,-3.4,3.4)
	htt_h2_phi  = TH1F("htt_h2_phi","TTbar H2 Phi Distribution", 40,-3.14,3.14)

	htt_hh_mass  = TH1F("htt_hh_mass","TTbar HH Mass Distribution", 40,0,800)
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

	hdy_h1_mass = TH1F("hdy_h1_mass","DY H1 Mass Distribution", 40,0,400)
	hdy_h1_pt   = TH1F("hdy_h1_pt","DY H1 Pt Distribution", 40,0,400)
	hdy_h1_eta  = TH1F("hdy_h1_eta","DY H1 Eta Distribution", 40,-3.4,3.4)
	hdy_h1_phi  = TH1F("hdy_h1_phi","DY H1 Phi Distribution", 40,-3.14,3.14)

	hdy_h2_mass = TH1F("hdy_h2_mass","DY H2 Mass Distribution", 40,0,400)
	hdy_h2_pt   = TH1F("hdy_h2_pt","DY H2 Pt Distribution", 40,0,400)
	hdy_h2_eta  = TH1F("hdy_h2_eta","DY H2 Eta Distribution", 40,-3.4,3.4)
	hdy_h2_phi  = TH1F("hdy_h2_phi","DY H2 Phi Distribution", 40,-3.14,3.14)

	hdy_hh_mass  = TH1F("hdy_hh_mass","DY HH Mass Distribution", 40,0,800)
	hdy_hh_pt    = TH1F("hdy_hh_pt","DY HH Pt Distribution", 40,0,600)
	hdy_hh_eta   = TH1F("hdy_hh_eta","DY HH Eta Distribution", 40,-3.4,3.4)
	hdy_hh_phi   = TH1F("hdy_hh_phi","DY HH Phi Distribution", 40,-3.14,3.14)


	hwj_j1_mass = TH1F("hwj_j1_mass","WJ J1 Mass Distribution", 40,0,200)
        hwj_j1_pt   = TH1F("hwj_j1_pt","WJ J1 Pt Distribution", 40,0,400)
        hwj_j1_eta  = TH1F("hwj_j1_eta","WJ J1 Eta Distribution", 40,-3.4,3.4)
        hwj_j1_phi  = TH1F("hwj_j1_phi","WJ J1 Phi Distribution", 40,-3.14,3.14)
        hwj_j1_btag = TH1F("hwj_j1_btag","WJ J1 Btag Distribution", 40,0,1.01)

        hwj_j2_mass = TH1F("hwj_j2_mass","WJ J2 Mass Distribution", 40,0,200)
        hwj_j2_pt   = TH1F("hwj_j2_pt","WJ J2 Pt Distribution", 40,0,400)
        hwj_j2_eta  = TH1F("hwj_j2_eta","WJ J2 Eta Distribution", 40,-3.4,3.4)
        hwj_j2_phi  = TH1F("hwj_j2_phi","WJ J2 Phi Distribution", 40,-3.14,3.14)
        hwj_j2_btag = TH1F("hwj_j2_btag","WJ J2 Btag Distribution", 40,0,1.01)

        hwj_l1_mass = TH1F("hwj_l1_mass","WJ L1 Mass Distribution", 40,0,200)
        hwj_l1_pt   = TH1F("hwj_l1_pt","WJ L1 Pt Distribution", 40,0,400)
        hwj_l1_eta  = TH1F("hwj_l1_eta","WJ L1 Eta Distribution", 40,-3.4,3.4)
        hwj_l1_phi  = TH1F("hwj_l1_phi","WJ L1 Phi Distribution", 40,-3.14,3.14)
        hwj_l1_charge = TH1F("hwj_l1_charge","WJ L1 Charge Distribution", 40,-1,1.01)

        hwj_l2_mass = TH1F("hwj_l2_mass","WJ L2 Mass Distribution", 40,0,200)
        hwj_l2_pt   = TH1F("hwj_l2_pt","WJ L2 Pt Distribution", 40,0,400)
        hwj_l2_eta  = TH1F("hwj_l2_eta","WJ L2 Eta Distribution", 40,-3.4,3.4)
        hwj_l2_phi  = TH1F("hwj_l2_phi","WJ L2 Phi Distribution", 40,-3.14,3.14)
        hwj_l2_charge = TH1F("hwj_l2_charge","WJ L2 Charge Distribution", 40,-1,1.01)

        hwj_dr_ll   = TH1F("hwj_dr_ll","WJ DeltaR Between LL", 40,0,10)
        hwj_dr_jj   = TH1F("hwj_dr_jj","WJ DeltaR Between JJ", 40,0,10)

        hwj_h1_mass = TH1F("hwj_h1_mass","WJ H1 Mass Distribution", 40,0,400)
        hwj_h1_pt   = TH1F("hwj_h1_pt","WJ H1 Pt Distribution", 40,0,400)
        hwj_h1_eta  = TH1F("hwj_h1_eta","WJ H1 Eta Distribution", 40,-3.4,3.4)
        hwj_h1_phi  = TH1F("hwj_h1_phi","WJ H1 Phi Distribution", 40,-3.14,3.14)

        hwj_h2_mass = TH1F("hwj_h2_mass","WJ H2 Mass Distribution", 40,0,400)
        hwj_h2_pt   = TH1F("hwj_h2_pt","WJ H2 Pt Distribution", 40,0,400)
        hwj_h2_eta  = TH1F("hwj_h2_eta","WJ H2 Eta Distribution", 40,-3.4,3.4)
        hwj_h2_phi  = TH1F("hwj_h2_phi","WJ H2 Phi Distribution", 40,-3.14,3.14)

        hwj_hh_mass = TH1F("hwj_hh_mass","WJ HH Mass Distribution", 40,0,800)
        hwj_hh_pt   = TH1F("hwj_hh_pt","WJ HH Pt Distribution", 40,0,600)
        hwj_hh_eta  = TH1F("hwj_hh_eta","WJ HH Eta Distribution", 40,-3.4,3.4)
        hwj_hh_phi  = TH1F("hwj_hh_phi","WJ HH Phi Distribution", 40,-3.14,3.14)


	hst_j1_mass = TH1F("hst_j1_mass","ST J1 Mass Distribution", 40,0,200)
        hst_j1_pt   = TH1F("hst_j1_pt","ST J1 Pt Distribution", 40,0,400)
        hst_j1_eta  = TH1F("hst_j1_eta","ST J1 Eta Distribution", 40,-3.4,3.4)
        hst_j1_phi  = TH1F("hst_j1_phi","ST J1 Phi Distribution", 40,-3.14,3.14)
        hst_j1_btag = TH1F("hst_j1_btag","ST J1 Btag Distribution", 40,0,1.01)

        hst_j2_mass = TH1F("hst_j2_mass","ST J2 Mass Distribution", 40,0,200)
        hst_j2_pt   = TH1F("hst_j2_pt","ST J2 Pt Distribution", 40,0,400)
        hst_j2_eta  = TH1F("hst_j2_eta","ST J2 Eta Distribution", 40,-3.4,3.4)
        hst_j2_phi  = TH1F("hst_j2_phi","ST J2 Phi Distribution", 40,-3.14,3.14)
        hst_j2_btag = TH1F("hst_j2_btag","ST J2 Btag Distribution", 40,0,1.01)

        hst_l1_mass = TH1F("hst_l1_mass","ST L1 Mass Distribution", 40,0,200)
        hst_l1_pt   = TH1F("hst_l1_pt","ST L1 Pt Distribution", 40,0,400)
        hst_l1_eta  = TH1F("hst_l1_eta","ST L1 Eta Distribution", 40,-3.4,3.4)
        hst_l1_phi  = TH1F("hst_l1_phi","ST L1 Phi Distribution", 40,-3.14,3.14)
        hst_l1_charge = TH1F("hst_l1_charge","ST L1 Charge Distribution", 40,-1,1.01)

        hst_l2_mass = TH1F("hst_l2_mass","ST L2 Mass Distribution", 40,0,200)
        hst_l2_pt   = TH1F("hst_l2_pt","ST L2 Pt Distribution", 40,0,400)
        hst_l2_eta  = TH1F("hst_l2_eta","ST L2 Eta Distribution", 40,-3.4,3.4)
        hst_l2_phi  = TH1F("hst_l2_phi","ST L2 Phi Distribution", 40,-3.14,3.14)
        hst_l2_charge = TH1F("hst_l2_charge","ST L2 Charge Distribution", 40,-1,1.01)

        hst_dr_ll   = TH1F("hst_dr_ll","ST DeltaR Between LL", 40,0,10)
        hst_dr_jj   = TH1F("hst_dr_jj","ST DeltaR Between JJ", 40,0,10)

	hst_h1_mass = TH1F("hst_h1_mass","ST H1 Mass Distribution", 40,0,400)
	hst_h1_pt   = TH1F("hst_h1_pt","ST H1 Pt Distribution", 40,0,400)
	hst_h1_eta  = TH1F("hst_h1_eta","ST H1 Eta Distribution", 40,-3.4,3.4)
	hst_h1_phi  = TH1F("hst_h1_phi","ST H1 Phi Distribution", 40,-3.14,3.14)

	hst_h2_mass = TH1F("hst_h2_mass","ST H2 Mass Distribution", 40,0,400)
        hst_h2_pt   = TH1F("hst_h2_pt","ST H2 Pt Distribution", 40,0,400)
        hst_h2_eta  = TH1F("hst_h2_eta","ST H2 Eta Distribution", 40,-3.4,3.4)
        hst_h2_phi  = TH1F("hst_h2_phi","ST H2 Phi Distribution", 40,-3.14,3.14)

	hst_hh_mass = TH1F("hst_hh_mass","ST HH Mass Distribution", 40,0,800)
        hst_hh_pt   = TH1F("hst_hh_pt","ST HH Pt Distribution", 40,0,600)
        hst_hh_eta  = TH1F("hst_hh_eta","ST HH Eta Distribution", 40,-3.4,3.4)
        hst_hh_phi  = TH1F("hst_hh_phi","ST HH Phi Distribution", 40,-3.14,3.14)


	stack_j1_mass = THStack("stack_j1_mass","Stacked Mass Distribution J1")
        stack_j1_pt   = THStack("stack_j1_pt","Stacked Pt Distribution J1")
        stack_j1_eta  = THStack("stack_j1_eta","Stacked Eta Distribution J1")
        stack_j1_phi  = THStack("stack_j1_phi","Stacked Phi Distribution J1")
        stack_j1_btag = THStack("stack_j1_btag","Stacked Btag Distribution J1")

        stack_j2_mass = THStack("stack_j2_mass","Stacked Mass Distribution J2")
        stack_j2_pt   = THStack("stack_j2_pt","Stacked Pt Distribution J2")
        stack_j2_eta  = THStack("stack_j2_eta","Stacked Eta Distribution J2")
        stack_j2_phi  = THStack("stack_j2_phi","Stacked Phi Distribution J2")
        stack_j2_btag = THStack("stack_j2_btag","Stacked Btag Distribution J2")

        stack_l1_mass = THStack("stack_l1_mass","Stacked Mass Distribution L1")
        stack_l1_pt   = THStack("stack_l1_pt","Stacked Pt Distribution L1")
        stack_l1_eta  = THStack("stack_l1_eta","Stacked Eta Distribution L1")
        stack_l1_phi  = THStack("stack_l1_phi","Stacked Phi Distribution L1")
        stack_l1_charge = THStack("stack_l1_charge","Stacked Charge Distribution L1")

        stack_l2_mass = THStack("stack_l2_mass","Stacked Mass Distribution L2")
        stack_l2_pt   = THStack("stack_l2_pt","Stacked Pt Distribution L2")
        stack_l2_eta  = THStack("stack_l2_eta","Stacked Eta Distribution L2")
        stack_l2_phi  = THStack("stack_l2_phi","Stacked Phi Distribution L2")
        stack_l2_charge = THStack("stack_l2_charge","Stacked Charge Distribution L2")

        stack_dr_ll   = THStack("stack_dr_ll","Stacked DeltaR Distribution LL")
        stack_dr_jj   = THStack("stack_dr_jj","Stacked DeltaR Distribution JJ")

        stack_h1_mass = THStack("stack_h1_mass","Stacked Mass Distribution H1<-ll")
        stack_h1_pt   = THStack("stack_h1_pt","Stacked Pt Distribution H1<-ll")
        stack_h1_eta  = THStack("stack_h1_eta","Stacked Eta Distribution H1<-ll")
        stack_h1_phi  = THStack("stack_h1_phi","Stacked Phi Distribution H1<-ll")

        stack_h2_mass = THStack("stack_h2_mass","Stacked Mass Distribution H2<-bb")
        stack_h2_pt   = THStack("stack_h2_pt","Stacked Pt Distribution H2<-bb")
        stack_h2_eta  = THStack("stack_h2_eta","Stacked Eta Distribution H2<-bb")
        stack_h2_phi  = THStack("stack_h2_phi","Stacked Phi Distribution H2<-bb")

        stack_hh_mass = THStack("stack_hh_mass","Stacked Mass Distribution HH")
        stack_hh_pt   = THStack("stack_hh_pt","Stacked Pt Distribution HH")
        stack_hh_eta  = THStack("stack_hh_eta","Stacked Eta Distribution HH")
        stack_hh_phi  = THStack("stack_hh_phi","Stacked Phi Distribution HH")


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


	wj_j1_mass = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j1_pt   = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j1_eta  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j1_phi  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j1_btag = np.zeros(hwj_j1_btag.GetNbinsX()+1)

        wj_j2_mass = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j2_pt   = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j2_eta  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j2_phi  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_j2_btag = np.zeros(hwj_h1_mass.GetNbinsX()+1)

        wj_l1_mass = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l1_pt   = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l1_eta  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l1_phi  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l1_charge = np.zeros(hwj_l1_charge.GetNbinsX()+1)

        wj_l2_mass = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l2_pt   = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l2_eta  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l2_phi  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_l2_charge = np.zeros(hwj_l1_charge.GetNbinsX()+1)

        wj_dr_ll   = np.zeros(hwj_dr_ll.GetNbinsX()+1)
        wj_dr_jj   = np.zeros(hwj_dr_ll.GetNbinsX()+1)

        wj_h1_mass = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_h1_pt   = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_h1_eta  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_h1_phi  = np.zeros(hwj_h1_mass.GetNbinsX()+1)

        wj_h2_mass = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_h2_pt   = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_h2_eta  = np.zeros(hwj_h1_mass.GetNbinsX()+1)
        wj_h2_phi  = np.zeros(hwj_h1_mass.GetNbinsX()+1)

        wj_hh_mass = np.zeros(hwj_hh_mass.GetNbinsX()+1)
        wj_hh_pt   = np.zeros(hwj_hh_mass.GetNbinsX()+1)
        wj_hh_eta  = np.zeros(hwj_hh_mass.GetNbinsX()+1)
        wj_hh_phi  = np.zeros(hwj_hh_mass.GetNbinsX()+1)


	st_j1_mass = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j1_pt   = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j1_eta  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j1_phi  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j1_btag = np.zeros(hst_j1_btag.GetNbinsX()+1)

        st_j2_mass = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j2_pt   = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j2_eta  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j2_phi  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_j2_btag = np.zeros(hst_h1_mass.GetNbinsX()+1)

        st_l1_mass = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l1_pt   = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l1_eta  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l1_phi  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l1_charge = np.zeros(hst_l1_charge.GetNbinsX()+1)

        st_l2_mass = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l2_pt   = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l2_eta  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l2_phi  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_l2_charge = np.zeros(hst_l1_charge.GetNbinsX()+1)

        st_dr_ll   = np.zeros(hst_dr_ll.GetNbinsX()+1)
        st_dr_jj   = np.zeros(hst_dr_ll.GetNbinsX()+1)

        st_h1_mass = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_h1_pt   = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_h1_eta  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_h1_phi  = np.zeros(hst_h1_mass.GetNbinsX()+1)

        st_h2_mass = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_h2_pt   = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_h2_eta  = np.zeros(hst_h1_mass.GetNbinsX()+1)
        st_h2_phi  = np.zeros(hst_h1_mass.GetNbinsX()+1)

        st_hh_mass = np.zeros(hst_hh_mass.GetNbinsX()+1)
        st_hh_pt   = np.zeros(hst_hh_mass.GetNbinsX()+1)
        st_hh_eta  = np.zeros(hst_hh_mass.GetNbinsX()+1)
        st_hh_phi  = np.zeros(hst_hh_mass.GetNbinsX()+1)

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

		rh1_yield,rh2_yield,rhh_yield = 0,0,0
		
		# Real
                for i in range(1,h1_mass.GetNbinsX()+1):
			if (dr == 'run2016B' or dr == 'run2016C' or dr == 'run2016D' or dr == 'run2016E' or dr == 'run2016F' or dr == 'run2016G' or dr == 'run2016H'):
				real_j1_mass[i] += j1_mass.GetBinContent(i)
				real_j1_pt[i]   += j1_pt.GetBinContent(i)
				real_j1_eta[i]  += j1_eta.GetBinContent(i)
				real_j1_phi[i]  += j1_phi.GetBinContent(i)
				real_j1_btag[i] += j1_btag.GetBinContent(i)

				real_j2_mass[i] += j2_mass.GetBinContent(i)
				real_j2_pt[i]   += j2_pt.GetBinContent(i)
				real_j2_eta[i]  += j2_eta.GetBinContent(i)
				real_j2_phi[i]  += j2_phi.GetBinContent(i)
				real_j2_btag[i] += j2_btag.GetBinContent(i)

				real_l1_mass[i] += l1_mass.GetBinContent(i)
				real_l1_pt[i]   += l1_pt.GetBinContent(i)
				real_l1_eta[i]  += l1_eta.GetBinContent(i)
				real_l1_phi[i]  += l1_phi.GetBinContent(i)
				real_l1_charge[i] += l1_charge.GetBinContent(i)

				real_l2_mass[i] += l2_mass.GetBinContent(i)
				real_l2_pt[i]   += l2_pt.GetBinContent(i)
				real_l2_eta[i]  += l2_eta.GetBinContent(i)
				real_l2_phi[i]  += l2_phi.GetBinContent(i)
				real_l2_charge[i] += l2_charge.GetBinContent(i)

				real_dr_ll[i]   += deltaR_ll.GetBinContent(i)
				real_dr_jj[i]   += deltaR_jj.GetBinContent(i)

				real_h1_mass[i] += h1_mass.GetBinContent(i)
				real_h1_pt[i]   += h1_pt.GetBinContent(i)
				real_h1_eta[i]  += h1_eta.GetBinContent(i)
				real_h1_phi[i]  += h1_phi.GetBinContent(i)

				real_h2_mass[i] += h2_mass.GetBinContent(i)
				real_h2_pt[i]   += h2_pt.GetBinContent(i)
				real_h2_eta[i]  += h2_eta.GetBinContent(i)
				real_h2_phi[i]  += h2_phi.GetBinContent(i)

				real_hh_mass[i] += hh_mass.GetBinContent(i)
				real_hh_pt[i]   += hh_pt.GetBinContent(i)
				real_hh_eta[i]  += hh_eta.GetBinContent(i)
				real_hh_phi[i]  += hh_phi.GetBinContent(i)

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

			rh1_yield += real_h1_mass[i]
			rh2_yield += real_h2_mass[i]
			rhh_yield += real_hh_mass[i]

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


		# Signal
		sf_signal = 1
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
		sf_tt_2l2nu = 3211572/4891620536.23#4891620336.1
		sf_tt_semi  = 13296780/32366940544.5#32366940321.3
		sf_tt_had   = 13732740/21432311062.9#21432310862.9
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
		sf_dy_m50 = 0 #220886400/1.89714160794e+12
		sf_dy_0j  = 175837200/2.83853141463e+11#2.83853150541e+11
		sf_dy_1j  = 32624340/2.24631956079e+11#2.24631958677e+11
		sf_dy_2j  = 12206880/1.02265928152e+11#1.02265928169e+11
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
		sf_wj_jet = 2234295000/86916455.0
		sf_wj_1j  = 342955200/45283121.0#45283121.0
		sf_wj_2j  = 118072500/60438768.0#60438768.0
		sf_wj_3j  = 41779500/59228199.0#59228199.0
		sf_wj_4j  = 23033220/27868119.0#27868119.0
		if (dr == 'wj_jet' or dr == 'wj_1j' or dr == 'wj_2j' or dr == 'wj_3j' or dr == 'wj_4j'):
			if (dr == 'wj_jet'):
                                for i in range(1,h1_mass.GetNbinsX()+1):

                                        wj_j1_mass[i] += j1_mass.GetBinContent(i)*sf_wj_jet
                                        wj_j1_pt[i]   += j1_pt.GetBinContent(i)*sf_wj_jet
                                        wj_j1_eta[i]  += j1_eta.GetBinContent(i)*sf_wj_jet
                                        wj_j1_phi[i]  += j1_phi.GetBinContent(i)*sf_wj_jet
                                        wj_j1_btag[i] += j1_btag.GetBinContent(i)*sf_wj_jet

                                        wj_j2_mass[i] += j2_mass.GetBinContent(i)*sf_wj_jet
                                        wj_j2_pt[i]   += j2_pt.GetBinContent(i)*sf_wj_jet
                                        wj_j2_eta[i]  += j2_eta.GetBinContent(i)*sf_wj_jet
                                        wj_j2_phi[i]  += j2_phi.GetBinContent(i)*sf_wj_jet
                                        wj_j2_btag[i] += j2_btag.GetBinContent(i)*sf_wj_jet

                                        wj_l1_mass[i] += l1_mass.GetBinContent(i)*sf_wj_jet
                                        wj_l1_pt[i]   += l1_pt.GetBinContent(i)*sf_wj_jet
                                        wj_l1_eta[i]  += l1_eta.GetBinContent(i)*sf_wj_jet
                                        wj_l1_phi[i]  += l1_phi.GetBinContent(i)*sf_wj_jet
                                        wj_l1_charge[i] += l1_charge.GetBinContent(i)*sf_wj_jet

                                        wj_l2_mass[i] += l2_mass.GetBinContent(i)*sf_wj_jet
                                        wj_l2_pt[i]   += l2_pt.GetBinContent(i)*sf_wj_jet
                                        wj_l2_eta[i]  += l2_eta.GetBinContent(i)*sf_wj_jet
                                        wj_l2_phi[i]  += l2_phi.GetBinContent(i)*sf_wj_jet
                                        wj_l2_charge[i] += l2_charge.GetBinContent(i)*sf_wj_jet

                                        wj_dr_ll[i]   += deltaR_ll.GetBinContent(i)*sf_wj_jet
                                        wj_dr_jj[i]   += deltaR_jj.GetBinContent(i)*sf_wj_jet

                                        wj_h1_mass[i] += h1_mass.GetBinContent(i)*sf_wj_jet
                                        wj_h1_pt[i]   += h1_pt.GetBinContent(i)*sf_wj_jet
                                        wj_h1_eta[i]  += h1_eta.GetBinContent(i)*sf_wj_jet
                                        wj_h1_phi[i]  += h1_phi.GetBinContent(i)*sf_wj_jet

                                        wj_h2_mass[i] += h2_mass.GetBinContent(i)*sf_wj_jet
                                        wj_h2_pt[i]   += h2_pt.GetBinContent(i)*sf_wj_jet
                                        wj_h2_eta[i]  += h2_eta.GetBinContent(i)*sf_wj_jet
                                        wj_h2_phi[i]  += h2_phi.GetBinContent(i)*sf_wj_jet

                                        wj_hh_mass[i] += hh_mass.GetBinContent(i)*sf_wj_jet
                                        wj_hh_pt[i]   += hh_pt.GetBinContent(i)*sf_wj_jet
                                        wj_hh_eta[i]  += hh_eta.GetBinContent(i)*sf_wj_jet
                                        wj_hh_phi[i]  += hh_phi.GetBinContent(i)*sf_wj_jet

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

		for i in range(1,h1_mass.GetNbinsX()+1):

                        hwj_j1_mass.SetBinContent(i,wj_j1_mass[i])
                        hwj_j1_pt.SetBinContent(i,wj_j1_pt[i])
                        hwj_j1_eta.SetBinContent(i,wj_j1_eta[i])
                        hwj_j1_phi.SetBinContent(i,wj_j1_phi[i])
                        hwj_j1_btag.SetBinContent(i,wj_j1_btag[i])

                        hwj_j2_mass.SetBinContent(i,wj_j2_mass[i])
                        hwj_j2_pt.SetBinContent(i,wj_j2_pt[i])
                        hwj_j2_eta.SetBinContent(i,wj_j2_eta[i])
                        hwj_j2_phi.SetBinContent(i,wj_j2_phi[i])
                        hwj_j2_btag.SetBinContent(i,wj_j2_btag[i])

                        hwj_l1_mass.SetBinContent(i,wj_l1_mass[i])
                        hwj_l1_pt.SetBinContent(i,wj_l1_pt[i])
                        hwj_l1_eta.SetBinContent(i,wj_l1_eta[i])
                        hwj_l1_phi.SetBinContent(i,wj_l1_phi[i])
                        hwj_l1_charge.SetBinContent(i,wj_l1_charge[i])

                        hwj_l2_mass.SetBinContent(i,wj_l2_mass[i])
                        hwj_l2_pt.SetBinContent(i,wj_l2_pt[i])
                        hwj_l2_eta.SetBinContent(i,wj_l2_eta[i])
                        hwj_l2_phi.SetBinContent(i,wj_l2_phi[i])
                        hwj_l2_charge.SetBinContent(i,wj_l2_charge[i])

                        hwj_dr_ll.SetBinContent(i,wj_dr_ll[i])
                        hwj_dr_jj.SetBinContent(i,wj_dr_jj[i])

                        hwj_h1_mass.SetBinContent(i,wj_h1_mass[i])
                        hwj_h1_pt.SetBinContent(i,wj_h1_pt[i])
                        hwj_h1_eta.SetBinContent(i,wj_h1_eta[i])
                        hwj_h1_phi.SetBinContent(i,wj_h1_phi[i])

                        hwj_h2_mass.SetBinContent(i,wj_h2_mass[i])
                        hwj_h2_pt.SetBinContent(i,wj_h2_pt[i])
                        hwj_h2_eta.SetBinContent(i,wj_h2_eta[i])
                        hwj_h2_phi.SetBinContent(i,wj_h2_phi[i])

                        hwj_hh_mass.SetBinContent(i,wj_hh_mass[i])
                        hwj_hh_pt.SetBinContent(i,wj_hh_pt[i])
                        hwj_hh_eta.SetBinContent(i,wj_hh_eta[i])
                        hwj_hh_phi.SetBinContent(i,wj_hh_phi[i])

			#NOTE: Error should not be symmetric up and down
                        hwj_h1_mass.SetBinError(i,wj_h1_mass[i]**0.5)
                        hwj_h1_pt.SetBinError(i,  wj_h1_pt[i]**0.5)
                        hwj_h1_eta.SetBinError(i, wj_h1_eta[i]**0.5)
                        hwj_h1_phi.SetBinError(i, wj_h1_phi[i]**0.5)

                        hwj_h2_mass.SetBinError(i,wj_h2_mass[i]**0.5)
                        hwj_h2_pt.SetBinError(i,  wj_h2_pt[i]**0.5)
                        hwj_h2_eta.SetBinError(i, wj_h2_eta[i]**0.5)
                        hwj_h2_phi.SetBinError(i, wj_h2_phi[i]**0.5)

                        hwj_hh_mass.SetBinError(i,wj_hh_mass[i]**0.5)
                        hwj_hh_pt.SetBinError(i,  wj_hh_pt[i]**0.5)
                        hwj_hh_eta.SetBinError(i, wj_hh_eta[i]**0.5)
                        hwj_hh_phi.SetBinError(i, wj_hh_phi[i]**0.5)


		# ST
                sf_st_schannel_incl = 122069/3370668.41375#3370668.5184
                sf_st_tW_anti = 1300614/6933094.0#6933094.0
                sf_st_tW_top  = 1300614/6952830.0#6952830.0
		sf_st_tchannel_anti = 2942730/38811017.0#38811017.0
		sf_st_tchannel_top  = 4940880/58403420.0#58403420.0
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

                        hst_j1_mass.SetBinContent(i,st_j1_mass[i])
                        hst_j1_pt.SetBinContent(i,st_j1_pt[i])
                        hst_j1_eta.SetBinContent(i,st_j1_eta[i])
                        hst_j1_phi.SetBinContent(i,st_j1_phi[i])
                        hst_j1_btag.SetBinContent(i,st_j1_btag[i])

                        hst_j2_mass.SetBinContent(i,st_j2_mass[i])
                        hst_j2_pt.SetBinContent(i,st_j2_pt[i])
                        hst_j2_eta.SetBinContent(i,st_j2_eta[i])
                        hst_j2_phi.SetBinContent(i,st_j2_phi[i])
                        hst_j2_btag.SetBinContent(i,st_j2_btag[i])

                        hst_l1_mass.SetBinContent(i,st_l1_mass[i])
                        hst_l1_pt.SetBinContent(i,st_l1_pt[i])
                        hst_l1_eta.SetBinContent(i,st_l1_eta[i])
                        hst_l1_phi.SetBinContent(i,st_l1_phi[i])
                        hst_l1_charge.SetBinContent(i,st_l1_charge[i])

                        hst_l2_mass.SetBinContent(i,st_l2_mass[i])
                        hst_l2_pt.SetBinContent(i,st_l2_pt[i])
                        hst_l2_eta.SetBinContent(i,st_l2_eta[i])
                        hst_l2_phi.SetBinContent(i,st_l2_phi[i])
                        hst_l2_charge.SetBinContent(i,st_l2_charge[i])

                        hst_dr_ll.SetBinContent(i,st_dr_ll[i])
                        hst_dr_jj.SetBinContent(i,st_dr_jj[i])

                        hst_h1_mass.SetBinContent(i,st_h1_mass[i])
                        hst_h1_pt.SetBinContent(i,st_h1_pt[i])
                        hst_h1_eta.SetBinContent(i,st_h1_eta[i])
                        hst_h1_phi.SetBinContent(i,st_h1_phi[i])

                        hst_h2_mass.SetBinContent(i,st_h2_mass[i])
                        hst_h2_pt.SetBinContent(i,st_h2_pt[i])
                        hst_h2_eta.SetBinContent(i,st_h2_eta[i])
                        hst_h2_phi.SetBinContent(i,st_h2_phi[i])

                        hst_hh_mass.SetBinContent(i,st_hh_mass[i])
                        hst_hh_pt.SetBinContent(i,st_hh_pt[i])
                        hst_hh_eta.SetBinContent(i,st_hh_eta[i])
                        hst_hh_phi.SetBinContent(i,st_hh_phi[i])

			#NOTE: Error should not be symmetric up and down
                        hst_h1_mass.SetBinError(i,st_h1_mass[i]**0.5)
                        hst_h1_pt.SetBinError(i,  st_h1_pt[i]**0.5)
                        hst_h1_eta.SetBinError(i, st_h1_eta[i]**0.5)
                        hst_h1_phi.SetBinError(i, st_h1_phi[i]**0.5)

                        hst_h2_mass.SetBinError(i,st_h2_mass[i]**0.5)
                        hst_h2_pt.SetBinError(i,  st_h2_pt[i]**0.5)
                        hst_h2_eta.SetBinError(i, st_h2_eta[i]**0.5)
                        hst_h2_phi.SetBinError(i, st_h2_phi[i]**0.5)

                        hst_hh_mass.SetBinError(i,st_hh_mass[i]**0.5)
                        hst_hh_pt.SetBinError(i,  st_hh_pt[i]**0.5)
                        hst_hh_eta.SetBinError(i, st_hh_eta[i]**0.5)
                        hst_hh_phi.SetBinError(i, st_hh_phi[i]**0.5)


		f.Close()


	htt_j1_mass.SetFillColor(31)
        hdy_j1_mass.SetFillColor(41)
        hwj_j1_mass.SetFillColor(47)
        hst_j1_mass.SetFillColor(48)
        stack_j1_mass.Add(htt_j1_mass)
        stack_j1_mass.Add(hdy_j1_mass)
        stack_j1_mass.Add(hwj_j1_mass)
        stack_j1_mass.Add(hst_j1_mass)

        htt_j1_pt.SetFillColor(31)
        hdy_j1_pt.SetFillColor(41)
        hwj_j1_pt.SetFillColor(47)
        hst_j1_pt.SetFillColor(48)
        stack_j1_pt.Add(htt_j1_pt)
        stack_j1_pt.Add(hdy_j1_pt)
        stack_j1_pt.Add(hwj_j1_pt)
        stack_j1_pt.Add(hst_j1_pt)

        htt_j1_eta.SetFillColor(31)
        hdy_j1_eta.SetFillColor(41)
        hwj_j1_eta.SetFillColor(47)
        hst_j1_eta.SetFillColor(48)
        stack_j1_eta.Add(htt_j1_eta)
        stack_j1_eta.Add(hdy_j1_eta)
        stack_j1_eta.Add(hwj_j1_eta)
        stack_j1_eta.Add(hst_j1_eta)

        htt_j1_phi.SetFillColor(31)
        hdy_j1_phi.SetFillColor(41)
        hwj_j1_phi.SetFillColor(47)
        hst_j1_phi.SetFillColor(48)
        stack_j1_phi.Add(htt_j1_phi)
        stack_j1_phi.Add(hdy_j1_phi)
        stack_j1_phi.Add(hwj_j1_phi)
        stack_j1_phi.Add(hst_j1_phi)

	htt_j1_btag.SetFillColor(31)
        hdy_j1_btag.SetFillColor(41)
        hwj_j1_btag.SetFillColor(47)
        hst_j1_btag.SetFillColor(48)
        stack_j1_btag.Add(htt_j1_btag)
        stack_j1_btag.Add(hdy_j1_btag)
        stack_j1_btag.Add(hwj_j1_btag)
        stack_j1_btag.Add(hst_j1_btag)

	htt_j2_mass.SetFillColor(31)
        hdy_j2_mass.SetFillColor(41)
        hwj_j2_mass.SetFillColor(47)
        hst_j2_mass.SetFillColor(48)
        stack_j2_mass.Add(htt_j2_mass)
        stack_j2_mass.Add(hdy_j2_mass)
        stack_j2_mass.Add(hwj_j2_mass)
        stack_j2_mass.Add(hst_j2_mass)

        htt_j2_pt.SetFillColor(31)
        hdy_j2_pt.SetFillColor(41)
        hwj_j2_pt.SetFillColor(47)
        hst_j2_pt.SetFillColor(48)
        stack_j2_pt.Add(htt_j2_pt)
        stack_j2_pt.Add(hdy_j2_pt)
        stack_j2_pt.Add(hwj_j2_pt)
        stack_j2_pt.Add(hst_j2_pt)

        htt_j2_eta.SetFillColor(31)
        hdy_j2_eta.SetFillColor(41)
        hwj_j2_eta.SetFillColor(47)
        hst_j2_eta.SetFillColor(48)
        stack_j2_eta.Add(htt_j2_eta)
        stack_j2_eta.Add(hdy_j2_eta)
        stack_j2_eta.Add(hwj_j2_eta)
        stack_j2_eta.Add(hst_j2_eta)

        htt_j2_phi.SetFillColor(31)
        hdy_j2_phi.SetFillColor(41)
        hwj_j2_phi.SetFillColor(47)
        hst_j2_phi.SetFillColor(48)
        stack_j2_phi.Add(htt_j2_phi)
        stack_j2_phi.Add(hdy_j2_phi)
        stack_j2_phi.Add(hwj_j2_phi)
        stack_j2_phi.Add(hst_j2_phi)

	htt_j2_btag.SetFillColor(31)
        hdy_j2_btag.SetFillColor(41)
        hwj_j2_btag.SetFillColor(47)
        hst_j2_btag.SetFillColor(48)
        stack_j2_btag.Add(htt_j2_btag)
        stack_j2_btag.Add(hdy_j2_btag)
        stack_j2_btag.Add(hwj_j2_btag)
        stack_j2_btag.Add(hst_j2_btag)


	htt_l1_mass.SetFillColor(31)
        hdy_l1_mass.SetFillColor(41)
        hwj_l1_mass.SetFillColor(47)
        hst_l1_mass.SetFillColor(48)
        stack_l1_mass.Add(htt_l1_mass)
        stack_l1_mass.Add(hdy_l1_mass)
        stack_l1_mass.Add(hwj_l1_mass)
        stack_l1_mass.Add(hst_l1_mass)

        htt_l1_pt.SetFillColor(31)
        hdy_l1_pt.SetFillColor(41)
        hwj_l1_pt.SetFillColor(47)
        hst_l1_pt.SetFillColor(48)
        stack_l1_pt.Add(htt_l1_pt)
        stack_l1_pt.Add(hdy_l1_pt)
        stack_l1_pt.Add(hwj_l1_pt)
        stack_l1_pt.Add(hst_l1_pt)

        htt_l1_eta.SetFillColor(31)
        hdy_l1_eta.SetFillColor(41)
        hwj_l1_eta.SetFillColor(47)
        hst_l1_eta.SetFillColor(48)
        stack_l1_eta.Add(htt_l1_eta)
        stack_l1_eta.Add(hdy_l1_eta)
        stack_l1_eta.Add(hwj_l1_eta)
        stack_l1_eta.Add(hst_l1_eta)

        htt_l1_phi.SetFillColor(31)
        hdy_l1_phi.SetFillColor(41)
        hwj_l1_phi.SetFillColor(47)
        hst_l1_phi.SetFillColor(48)
        stack_l1_phi.Add(htt_l1_phi)
        stack_l1_phi.Add(hdy_l1_phi)
        stack_l1_phi.Add(hwj_l1_phi)
        stack_l1_phi.Add(hst_l1_phi)

	htt_l1_charge.SetFillColor(31)
        hdy_l1_charge.SetFillColor(41)
        hwj_l1_charge.SetFillColor(47)
        hst_l1_charge.SetFillColor(48)
        stack_l1_charge.Add(htt_l1_charge)
        stack_l1_charge.Add(hdy_l1_charge)
        stack_l1_charge.Add(hwj_l1_charge)
        stack_l1_charge.Add(hst_l1_charge)


	htt_l2_mass.SetFillColor(31)
        hdy_l2_mass.SetFillColor(41)
        hwj_l2_mass.SetFillColor(47)
        hst_l2_mass.SetFillColor(48)
        stack_l2_mass.Add(htt_l2_mass)
        stack_l2_mass.Add(hdy_l2_mass)
        stack_l2_mass.Add(hwj_l2_mass)
        stack_l2_mass.Add(hst_l2_mass)

        htt_l2_pt.SetFillColor(31)
        hdy_l2_pt.SetFillColor(41)
        hwj_l2_pt.SetFillColor(47)
        hst_l2_pt.SetFillColor(48)
        stack_l2_pt.Add(htt_l2_pt)
        stack_l2_pt.Add(hdy_l2_pt)
        stack_l2_pt.Add(hwj_l2_pt)
        stack_l2_pt.Add(hst_l2_pt)

        htt_l2_eta.SetFillColor(31)
        hdy_l2_eta.SetFillColor(41)
        hwj_l2_eta.SetFillColor(47)
        hst_l2_eta.SetFillColor(48)
        stack_l2_eta.Add(htt_l2_eta)
        stack_l2_eta.Add(hdy_l2_eta)
        stack_l2_eta.Add(hwj_l2_eta)
        stack_l2_eta.Add(hst_l2_eta)

        htt_l2_phi.SetFillColor(31)
        hdy_l2_phi.SetFillColor(41)
        hwj_l2_phi.SetFillColor(47)
        hst_l2_phi.SetFillColor(48)
        stack_l2_phi.Add(htt_l2_phi)
        stack_l2_phi.Add(hdy_l2_phi)
        stack_l2_phi.Add(hwj_l2_phi)
        stack_l2_phi.Add(hst_l2_phi)

	htt_l2_charge.SetFillColor(31)
        hdy_l2_charge.SetFillColor(41)
        hwj_l2_charge.SetFillColor(47)
        hst_l2_charge.SetFillColor(48)
        stack_l2_charge.Add(htt_l2_charge)
        stack_l2_charge.Add(hdy_l2_charge)
        stack_l2_charge.Add(hwj_l2_charge)
        stack_l2_charge.Add(hst_l2_charge)

	htt_dr_ll.SetFillColor(31)
	hdy_dr_ll.SetFillColor(41)
	hwj_dr_ll.SetFillColor(47)
	hst_dr_ll.SetFillColor(48)
	stack_dr_ll.Add(htt_dr_ll)
	stack_dr_ll.Add(hdy_dr_ll)
	stack_dr_ll.Add(hwj_dr_ll)
	stack_dr_ll.Add(hst_dr_ll)

	htt_dr_jj.SetFillColor(31)
        hdy_dr_jj.SetFillColor(41)
        hwj_dr_jj.SetFillColor(47)
        hst_dr_jj.SetFillColor(48)
        stack_dr_jj.Add(htt_dr_jj)
        stack_dr_jj.Add(hdy_dr_jj)
        stack_dr_jj.Add(hwj_dr_jj)
        stack_dr_jj.Add(hst_dr_jj)

	htt_h1_mass.SetFillColor(31)
	hdy_h1_mass.SetFillColor(41)
	hwj_h1_mass.SetFillColor(47)
	hst_h1_mass.SetFillColor(48)
	stack_h1_mass.Add(htt_h1_mass)
	stack_h1_mass.Add(hdy_h1_mass)
	stack_h1_mass.Add(hwj_h1_mass)
	stack_h1_mass.Add(hst_h1_mass)

	htt_h1_pt.SetFillColor(31)
        hdy_h1_pt.SetFillColor(41)
        hwj_h1_pt.SetFillColor(47)
        hst_h1_pt.SetFillColor(48)
        stack_h1_pt.Add(htt_h1_pt)
        stack_h1_pt.Add(hdy_h1_pt)
        stack_h1_pt.Add(hwj_h1_pt)
        stack_h1_pt.Add(hst_h1_pt)

	htt_h1_eta.SetFillColor(31)
        hdy_h1_eta.SetFillColor(41)
        hwj_h1_eta.SetFillColor(47)
        hst_h1_eta.SetFillColor(48)
        stack_h1_eta.Add(htt_h1_eta)
        stack_h1_eta.Add(hdy_h1_eta)
        stack_h1_eta.Add(hwj_h1_eta)
        stack_h1_eta.Add(hst_h1_eta)

	htt_h1_phi.SetFillColor(31)
        hdy_h1_phi.SetFillColor(41)
        hwj_h1_phi.SetFillColor(47)
        hst_h1_phi.SetFillColor(48)
        stack_h1_phi.Add(htt_h1_phi)
        stack_h1_phi.Add(hdy_h1_phi)
        stack_h1_phi.Add(hwj_h1_phi)
        stack_h1_phi.Add(hst_h1_phi)


        htt_h2_mass.SetFillColor(31)
        hdy_h2_mass.SetFillColor(41)
	hwj_h2_mass.SetFillColor(47)
	hst_h2_mass.SetFillColor(48)
        stack_h2_mass.Add(htt_h2_mass)
        stack_h2_mass.Add(hdy_h2_mass)
	stack_h2_mass.Add(hwj_h2_mass)
	stack_h2_mass.Add(hst_h2_mass)

        htt_h2_pt.SetFillColor(31)
        hdy_h2_pt.SetFillColor(41)
        hwj_h2_pt.SetFillColor(47)
        hst_h2_pt.SetFillColor(48)
        stack_h2_pt.Add(htt_h2_pt)
        stack_h2_pt.Add(hdy_h2_pt)
        stack_h2_pt.Add(hwj_h2_pt)
        stack_h2_pt.Add(hst_h2_pt)

	htt_h2_eta.SetFillColor(31)
        hdy_h2_eta.SetFillColor(41)
        hwj_h2_eta.SetFillColor(47)
        hst_h2_eta.SetFillColor(48)
        stack_h2_eta.Add(htt_h2_eta)
        stack_h2_eta.Add(hdy_h2_eta)
        stack_h2_eta.Add(hwj_h2_eta)
        stack_h2_eta.Add(hst_h2_eta)

	htt_h2_phi.SetFillColor(31)
        hdy_h2_phi.SetFillColor(41)
        hwj_h2_phi.SetFillColor(47)
        hst_h2_phi.SetFillColor(48)
        stack_h2_phi.Add(htt_h2_phi)
        stack_h2_phi.Add(hdy_h2_phi)
        stack_h2_phi.Add(hwj_h2_phi)
        stack_h2_phi.Add(hst_h2_phi)


	htt_hh_mass.SetFillColor(31)
        hdy_hh_mass.SetFillColor(41)
        hwj_hh_mass.SetFillColor(47)
        hst_hh_mass.SetFillColor(48)
        stack_hh_mass.Add(htt_hh_mass)
        stack_hh_mass.Add(hdy_hh_mass)
        stack_hh_mass.Add(hwj_hh_mass)
        stack_hh_mass.Add(hst_hh_mass)

	htt_hh_pt.SetFillColor(31)
        hdy_hh_pt.SetFillColor(41)
        hwj_hh_pt.SetFillColor(47)
        hst_hh_pt.SetFillColor(48)
        stack_hh_pt.Add(htt_hh_pt)
        stack_hh_pt.Add(hdy_hh_pt)
        stack_hh_pt.Add(hwj_hh_pt)
        stack_hh_pt.Add(hst_hh_pt)

	htt_hh_eta.SetFillColor(31)
        hdy_hh_eta.SetFillColor(41)
        hwj_hh_eta.SetFillColor(47)
        hst_hh_eta.SetFillColor(48)
        stack_hh_eta.Add(htt_hh_eta)
        stack_hh_eta.Add(hdy_hh_eta)
        stack_hh_eta.Add(hwj_hh_eta)
        stack_hh_eta.Add(hst_hh_eta)

	htt_hh_phi.SetFillColor(31)
        hdy_hh_phi.SetFillColor(41)
        hwj_hh_phi.SetFillColor(47)
        hst_hh_phi.SetFillColor(48)
        stack_hh_phi.Add(htt_hh_phi)
        stack_hh_phi.Add(hdy_hh_phi)
        stack_hh_phi.Add(hwj_hh_phi)
        stack_hh_phi.Add(hst_hh_phi)


	err_j1_mass = TH1F("err_j1_mass","Error J1 Mass Distribution", 40,0,200)
        err_j1_pt   = TH1F("err_j1_pt","Error J1 Pt Distribution", 40,0,400)
        err_j1_eta  = TH1F("err_j1_eta","Error J1 Eta Distribution", 40,-3.4,3.4)
        err_j1_phi  = TH1F("err_j1_phi","Error J1 Phi Distribution", 40,-3.14,3.14)
        err_j1_btag = TH1F("err_j1_btag","Error J1 Btag Distribution", 40,0,1.01)

        err_j2_mass = TH1F("err_j2_mass","Error J2 Mass Distribution", 40,0,200)
        err_j2_pt   = TH1F("err_j2_pt","Error J2 Pt Distribution", 40,0,400)
        err_j2_eta  = TH1F("err_j2_eta","Error J2 Eta Distribution", 40,-3.4,3.4)
        err_j2_phi  = TH1F("err_j2_phi","Error J2 Phi Distribution", 40,-3.14,3.14)
        err_j2_btag = TH1F("err_j2_btag","Error J2 Btag Distribution", 40,0,1.01)

        err_l1_mass = TH1F("err_l1_mass","Error L1 Mass Distribution", 40,0,200)
        err_l1_pt   = TH1F("err_l1_pt","Error L1 Pt Distribution", 40,0,400)
        err_l1_eta  = TH1F("err_l1_eta","Error L1 Eta Distribution", 40,-3.4,3.4)
        err_l1_phi  = TH1F("err_l1_phi","Error L1 Phi Distribution", 40,-3.14,3.14)
        err_l1_charge = TH1F("err_l1_charge","Error L1 Charge Distribution", 40,-1,1.01)

        err_l2_mass = TH1F("err_l2_mass","Error L2 Mass Distribution", 40,0,200)
        err_l2_pt   = TH1F("err_l2_pt","Error L2 Pt Distribution", 40,0,400)
        err_l2_eta  = TH1F("err_l2_eta","Error L2 Eta Distribution", 40,-3.4,3.4)
        err_l2_phi  = TH1F("err_l2_phi","Error L2 Phi Distribution", 40,-3.14,3.14)
        err_l2_charge = TH1F("err_l2_charge","Error L2 Charge Distribution", 40,-1,1.01)

        err_dr_ll   = TH1F("err_dr_ll","Error DeltaR Between LL", 40,0,10)
        err_dr_jj   = TH1F("err_dr_jj","Error DeltaR Between JJ", 40,0,10)

        err_h1_mass = TH1F("err_h1_mass","Error H1 Mass Distribution", 40,0,400)
        err_h1_pt   = TH1F("err_h1_pt","Error H1 Pt Distribution", 40,0,400)
        err_h1_eta  = TH1F("err_h1_eta","Error H1 Eta Distribution", 40,-3.4,3.4)
        err_h1_phi  = TH1F("err_h1_phi","Error H1 Phi Distribution", 40,-3.14,3.14)

        err_h2_mass = TH1F("err_h2_mass","Error H2 Mass Distribution", 40,0,400)
        err_h2_pt   = TH1F("err_h2_pt","Error H2 Pt Distribution", 40,0,400)
        err_h2_eta  = TH1F("err_h2_eta","Error H2 Eta Distribution", 40,-3.4,3.4)
        err_h2_phi  = TH1F("err_h2_phi","Error H2 Phi Distribution", 40,-3.14,3.14)

        err_hh_mass = TH1F("err_hh_mass","Error HH Mass Distribution", 40,0,800)
        err_hh_pt   = TH1F("err_hh_pt","Error HH Pt Distribution", 40,0,600)
        err_hh_eta  = TH1F("err_hh_eta","Error HH Eta Distribution", 40,-3.4,3.4)
        err_hh_phi  = TH1F("err_hh_phi","Error HH Phi Distribution", 40,-3.14,3.14)

	h1_yield,h2_yield,hh_yield = 0,0,0
	for i in range(1,err_j1_mass.GetNbinsX()+1):

		err_j1_mass.SetBinContent(i,htt_j1_mass.GetBinContent(i)+hdy_j1_mass.GetBinContent(i)+hwj_j1_mass.GetBinContent(i)+hst_j1_mass.GetBinContent(i))
                err_j1_pt.SetBinContent(i,htt_j1_pt.GetBinContent(i)+hdy_j1_pt.GetBinContent(i)+hwj_j1_pt.GetBinContent(i)+hst_j1_pt.GetBinContent(i))
                err_j1_eta.SetBinContent(i,htt_j1_eta.GetBinContent(i)+hdy_j1_eta.GetBinContent(i)+hwj_j1_eta.GetBinContent(i)+hst_j1_eta.GetBinContent(i))
                err_j1_phi.SetBinContent(i,htt_j1_phi.GetBinContent(i)+hdy_j1_phi.GetBinContent(i)+hwj_j1_phi.GetBinContent(i)+hst_j1_phi.GetBinContent(i))
                err_j1_btag.SetBinContent(i,htt_j1_btag.GetBinContent(i)+hdy_j1_btag.GetBinContent(i)+hwj_j1_btag.GetBinContent(i)+hst_j1_btag.GetBinContent(i))

                err_j2_mass.SetBinContent(i,htt_j2_mass.GetBinContent(i)+hdy_j2_mass.GetBinContent(i)+hwj_j2_mass.GetBinContent(i)+hst_j2_mass.GetBinContent(i))
                err_j2_pt.SetBinContent(i,htt_j2_pt.GetBinContent(i)+hdy_j2_pt.GetBinContent(i)+hwj_j2_pt.GetBinContent(i)+hst_j2_pt.GetBinContent(i))
                err_j2_eta.SetBinContent(i,htt_j2_eta.GetBinContent(i)+hdy_j2_eta.GetBinContent(i)+hwj_j2_eta.GetBinContent(i)+hst_j2_eta.GetBinContent(i))
                err_j2_phi.SetBinContent(i,htt_j2_phi.GetBinContent(i)+hdy_j2_phi.GetBinContent(i)+hwj_j2_phi.GetBinContent(i)+hst_j2_phi.GetBinContent(i))
                err_j2_btag.SetBinContent(i,htt_j2_btag.GetBinContent(i)+hdy_j2_btag.GetBinContent(i)+hwj_j2_btag.GetBinContent(i)+hst_j2_btag.GetBinContent(i))

                err_l1_mass.SetBinContent(i,htt_l1_mass.GetBinContent(i)+hdy_l1_mass.GetBinContent(i)+hwj_l1_mass.GetBinContent(i)+hst_l1_mass.GetBinContent(i))
                err_l1_pt.SetBinContent(i,htt_l1_pt.GetBinContent(i)+hdy_l1_pt.GetBinContent(i)+hwj_l1_pt.GetBinContent(i)+hst_l1_pt.GetBinContent(i))
                err_l1_eta.SetBinContent(i,htt_l1_eta.GetBinContent(i)+hdy_l1_eta.GetBinContent(i)+hwj_l1_eta.GetBinContent(i)+hst_l1_eta.GetBinContent(i))
                err_l1_phi.SetBinContent(i,htt_l1_phi.GetBinContent(i)+hdy_l1_phi.GetBinContent(i)+hwj_l1_phi.GetBinContent(i)+hst_l1_phi.GetBinContent(i))
                err_l1_charge.SetBinContent(i,htt_l1_charge.GetBinContent(i)+hdy_l1_charge.GetBinContent(i)+hwj_l1_charge.GetBinContent(i)+hst_l1_charge.GetBinContent(i))

                err_l2_mass.SetBinContent(i,htt_l2_mass.GetBinContent(i)+hdy_l2_mass.GetBinContent(i)+hwj_l2_mass.GetBinContent(i)+hst_l2_mass.GetBinContent(i))
                err_l2_pt.SetBinContent(i,htt_l2_pt.GetBinContent(i)+hdy_l2_pt.GetBinContent(i)+hwj_l2_pt.GetBinContent(i)+hst_l2_pt.GetBinContent(i))
                err_l2_eta.SetBinContent(i,htt_l2_eta.GetBinContent(i)+hdy_l2_eta.GetBinContent(i)+hwj_l2_eta.GetBinContent(i)+hst_l2_eta.GetBinContent(i))
                err_l2_phi.SetBinContent(i,htt_l2_phi.GetBinContent(i)+hdy_l2_phi.GetBinContent(i)+hwj_l2_phi.GetBinContent(i)+hst_l2_phi.GetBinContent(i))
                err_l2_charge.SetBinContent(i,htt_l2_charge.GetBinContent(i)+hdy_l2_charge.GetBinContent(i)+hwj_l2_charge.GetBinContent(i)+hst_l2_charge.GetBinContent(i))

                err_dr_ll.SetBinContent(i,htt_dr_ll.GetBinContent(i)+hdy_dr_ll.GetBinContent(i)+hwj_dr_ll.GetBinContent(i)+hst_dr_ll.GetBinContent(i))
                err_dr_jj.SetBinContent(i,htt_dr_jj.GetBinContent(i)+hdy_dr_jj.GetBinContent(i)+hwj_dr_jj.GetBinContent(i)+hst_dr_jj.GetBinContent(i))

                err_h1_mass.SetBinContent(i,htt_h1_mass.GetBinContent(i)+hdy_h1_mass.GetBinContent(i)+hwj_h1_mass.GetBinContent(i)+hst_h1_mass.GetBinContent(i))
                err_h1_pt.SetBinContent(i,htt_h1_pt.GetBinContent(i)+hdy_h1_pt.GetBinContent(i)+hwj_h1_pt.GetBinContent(i)+hst_h1_pt.GetBinContent(i))
                err_h1_eta.SetBinContent(i,htt_h1_eta.GetBinContent(i)+hdy_h1_eta.GetBinContent(i)+hwj_h1_eta.GetBinContent(i)+hst_h1_eta.GetBinContent(i))
                err_h1_phi.SetBinContent(i,htt_h1_phi.GetBinContent(i)+hdy_h1_phi.GetBinContent(i)+hwj_h1_phi.GetBinContent(i)+hst_h1_phi.GetBinContent(i))

                err_h2_mass.SetBinContent(i,htt_h2_mass.GetBinContent(i)+hdy_h2_mass.GetBinContent(i)+hwj_h2_mass.GetBinContent(i)+hst_h2_mass.GetBinContent(i))
                err_h2_pt.SetBinContent(i,htt_h2_pt.GetBinContent(i)+hdy_h2_pt.GetBinContent(i)+hwj_h2_pt.GetBinContent(i)+hst_h2_pt.GetBinContent(i))
                err_h2_eta.SetBinContent(i,htt_h2_eta.GetBinContent(i)+hdy_h2_eta.GetBinContent(i)+hwj_h2_eta.GetBinContent(i)+hst_h2_eta.GetBinContent(i))
                err_h2_phi.SetBinContent(i,htt_h2_phi.GetBinContent(i)+hdy_h2_phi.GetBinContent(i)+hwj_h2_phi.GetBinContent(i)+hst_h2_phi.GetBinContent(i))

                err_hh_mass.SetBinContent(i,htt_hh_mass.GetBinContent(i)+hdy_hh_mass.GetBinContent(i)+hwj_hh_mass.GetBinContent(i)+hst_hh_mass.GetBinContent(i))
                err_hh_pt.SetBinContent(i,htt_hh_pt.GetBinContent(i)+hdy_hh_pt.GetBinContent(i)+hwj_hh_pt.GetBinContent(i)+hst_hh_pt.GetBinContent(i))
                err_hh_eta.SetBinContent(i,htt_hh_eta.GetBinContent(i)+hdy_hh_eta.GetBinContent(i)+hwj_hh_eta.GetBinContent(i)+hst_hh_eta.GetBinContent(i))
                err_hh_phi.SetBinContent(i,htt_hh_phi.GetBinContent(i)+hdy_hh_phi.GetBinContent(i)+hwj_hh_phi.GetBinContent(i)+hst_hh_phi.GetBinContent(i))

		h1_yield += htt_h1_mass.GetBinContent(i)+hdy_h1_mass.GetBinContent(i)+hwj_h1_mass.GetBinContent(i)+hst_h1_mass.GetBinContent(i)
		h2_yield += htt_h2_mass.GetBinContent(i)+hdy_h2_mass.GetBinContent(i)+hwj_h2_mass.GetBinContent(i)+hst_h2_mass.GetBinContent(i)
		hh_yield += htt_hh_mass.GetBinContent(i)+hdy_hh_mass.GetBinContent(i)+hwj_hh_mass.GetBinContent(i)+hst_hh_mass.GetBinContent(i)

	for i in range(1,err_j1_mass.GetNbinsX()+1):

		err_j1_mass.SetBinError(i,(htt_j1_mass.GetBinError(i)**2+hdy_j1_mass.GetBinError(i)**2+hwj_j1_mass.GetBinError(i)**2+hst_j1_mass.GetBinError(i)**2)**.5)
		err_j1_pt.SetBinError(i,(htt_j1_pt.GetBinError(i)**2+hdy_j1_pt.GetBinError(i)**2+hwj_j1_pt.GetBinError(i)**2+hst_j1_pt.GetBinError(i)**2)**.5)
		err_j1_eta.SetBinError(i,(htt_j1_eta.GetBinError(i)**2+hdy_j1_eta.GetBinError(i)**2+hwj_j1_eta.GetBinError(i)**2+hst_j1_eta.GetBinError(i)**2)**.5)
		err_j1_phi.SetBinError(i,(htt_j1_phi.GetBinError(i)**2+hdy_j1_phi.GetBinError(i)**2+hwj_j1_phi.GetBinError(i)**2+hst_j1_phi.GetBinError(i)**2)**.5)
		err_j1_btag.SetBinError(i,(htt_j1_btag.GetBinError(i)**2+hdy_j1_btag.GetBinError(i)**2+hwj_j1_btag.GetBinError(i)**2+hst_j1_btag.GetBinError(i)**2)**.5)

		err_j2_mass.SetBinError(i,(htt_j2_mass.GetBinError(i)**2+hdy_j2_mass.GetBinError(i)**2+hwj_j2_mass.GetBinError(i)**2+hst_j2_mass.GetBinError(i)**2)**.5)
                err_j2_pt.SetBinError(i,(htt_j2_pt.GetBinError(i)**2+hdy_j2_pt.GetBinError(i)**2+hwj_j2_pt.GetBinError(i)**2+hst_j2_pt.GetBinError(i)**2)**.5)
                err_j2_eta.SetBinError(i,(htt_j2_eta.GetBinError(i)**2+hdy_j2_eta.GetBinError(i)**2+hwj_j2_eta.GetBinError(i)**2+hst_j2_eta.GetBinError(i)**2)**.5)
                err_j2_phi.SetBinError(i,(htt_j2_phi.GetBinError(i)**2+hdy_j2_phi.GetBinError(i)**2+hwj_j2_phi.GetBinError(i)**2+hst_j2_phi.GetBinError(i)**2)**.5)
                err_j2_btag.SetBinError(i,(htt_j2_btag.GetBinError(i)**2+hdy_j2_btag.GetBinError(i)**2+hwj_j2_btag.GetBinError(i)**2+hst_j2_btag.GetBinError(i)**2)**.5)

		err_l1_mass.SetBinError(i,(htt_l1_mass.GetBinError(i)**2+hdy_l1_mass.GetBinError(i)**2+hwj_l1_mass.GetBinError(i)**2+hst_l1_mass.GetBinError(i)**2)**.5)
                err_l1_pt.SetBinError(i,(htt_l1_pt.GetBinError(i)**2+hdy_l1_pt.GetBinError(i)**2+hwj_l1_pt.GetBinError(i)**2+hst_l1_pt.GetBinError(i)**2)**.5)
                err_l1_eta.SetBinError(i,(htt_l1_eta.GetBinError(i)**2+hdy_l1_eta.GetBinError(i)**2+hwj_l1_eta.GetBinError(i)**2+hst_l1_eta.GetBinError(i)**2)**.5)
                err_l1_phi.SetBinError(i,(htt_l1_phi.GetBinError(i)**2+hdy_l1_phi.GetBinError(i)**2+hwj_l1_phi.GetBinError(i)**2+hst_l1_phi.GetBinError(i)**2)**.5)
                err_l1_charge.SetBinError(i,(htt_l1_charge.GetBinError(i)**2+hdy_l1_charge.GetBinError(i)**2+hwj_l1_charge.GetBinError(i)**2+hst_l1_charge.GetBinError(i)**2)**.5)

		err_l2_mass.SetBinError(i,(htt_l2_mass.GetBinError(i)**2+hdy_l2_mass.GetBinError(i)**2+hwj_l2_mass.GetBinError(i)**2+hst_l2_mass.GetBinError(i)**2)**.5)
                err_l2_pt.SetBinError(i,(htt_l2_pt.GetBinError(i)**2+hdy_l2_pt.GetBinError(i)**2+hwj_l2_pt.GetBinError(i)**2+hst_l2_pt.GetBinError(i)**2)**.5)
                err_l2_eta.SetBinError(i,(htt_l2_eta.GetBinError(i)**2+hdy_l2_eta.GetBinError(i)**2+hwj_l2_eta.GetBinError(i)**2+hst_l2_eta.GetBinError(i)**2)**.5)
                err_l2_phi.SetBinError(i,(htt_l2_phi.GetBinError(i)**2+hdy_l2_phi.GetBinError(i)**2+hwj_l2_phi.GetBinError(i)**2+hst_l2_phi.GetBinError(i)**2)**.5)
                err_l2_charge.SetBinError(i,(htt_l2_charge.GetBinError(i)**2+hdy_l2_charge.GetBinError(i)**2+hwj_l2_charge.GetBinError(i)**2+hst_l2_charge.GetBinError(i)**2)**.5)

		err_dr_ll.SetBinError(i,(htt_dr_ll.GetBinError(i)**2+hdy_dr_ll.GetBinError(i)**2+hwj_dr_ll.GetBinError(i)**2+hst_dr_ll.GetBinError(i)**2)**.5)
		err_dr_jj.SetBinError(i,(htt_dr_jj.GetBinError(i)**2+hdy_dr_jj.GetBinError(i)**2+hwj_dr_jj.GetBinError(i)**2+hst_dr_jj.GetBinError(i)**2)**.5)

		err_h1_mass.SetBinError(i,(htt_h1_mass.GetBinError(i)**2+hdy_h1_mass.GetBinError(i)**2+hwj_h1_mass.GetBinError(i)**2+hst_h1_mass.GetBinError(i)**2)**.5)
                err_h1_pt.SetBinError(i,(htt_h1_pt.GetBinError(i)**2+hdy_h1_pt.GetBinError(i)**2+hwj_h1_pt.GetBinError(i)**2+hst_h1_pt.GetBinError(i)**2)**.5)
                err_h1_eta.SetBinError(i,(htt_h1_eta.GetBinError(i)**2+hdy_h1_eta.GetBinError(i)**2+hwj_h1_eta.GetBinError(i)**2+hst_h1_eta.GetBinError(i)**2)**.5)
                err_h1_phi.SetBinError(i,(htt_h1_phi.GetBinError(i)**2+hdy_h1_phi.GetBinError(i)**2+hwj_h1_phi.GetBinError(i)**2+hst_h1_phi.GetBinError(i)**2)**.5)

		err_h2_mass.SetBinError(i,(htt_h2_mass.GetBinError(i)**2+hdy_h2_mass.GetBinError(i)**2+hwj_h2_mass.GetBinError(i)**2+hst_h2_mass.GetBinError(i)**2)**.5)
                err_h2_pt.SetBinError(i,(htt_h2_pt.GetBinError(i)**2+hdy_h2_pt.GetBinError(i)**2+hwj_h2_pt.GetBinError(i)**2+hst_h2_pt.GetBinError(i)**2)**.5)
                err_h2_eta.SetBinError(i,(htt_h2_eta.GetBinError(i)**2+hdy_h2_eta.GetBinError(i)**2+hwj_h2_eta.GetBinError(i)**2+hst_h2_eta.GetBinError(i)**2)**.5)
                err_h2_phi.SetBinError(i,(htt_h2_phi.GetBinError(i)**2+hdy_h2_phi.GetBinError(i)**2+hwj_h2_phi.GetBinError(i)**2+hst_h2_phi.GetBinError(i)**2)**.5)

		err_hh_mass.SetBinError(i,(htt_hh_mass.GetBinError(i)**2+hdy_hh_mass.GetBinError(i)**2+hwj_hh_mass.GetBinError(i)**2+hst_hh_mass.GetBinError(i)**2)**.5)
                err_hh_pt.SetBinError(i,(htt_hh_pt.GetBinError(i)**2+hdy_hh_pt.GetBinError(i)**2+hwj_hh_pt.GetBinError(i)**2+hst_hh_pt.GetBinError(i)**2)**.5)
                err_hh_eta.SetBinError(i,(htt_hh_eta.GetBinError(i)**2+hdy_hh_eta.GetBinError(i)**2+hwj_hh_eta.GetBinError(i)**2+hst_hh_eta.GetBinError(i)**2)**.5)
                err_hh_phi.SetBinError(i,(htt_hh_phi.GetBinError(i)**2+hdy_hh_phi.GetBinError(i)**2+hwj_hh_phi.GetBinError(i)**2+hst_hh_phi.GetBinError(i)**2)**.5)


	lst_j1_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_j1_mass.AddEntry(htt_j1_mass,"TT","f")
        lst_j1_mass.AddEntry(hdy_j1_mass,"DY","f")
        lst_j1_mass.AddEntry(hwj_j1_mass,"WJ","f")
        lst_j1_mass.AddEntry(hst_j1_mass,"ST","f")
        lst_j1_mass.AddEntry(hsig_j1_mass,"Signal","l")
        lst_j1_mass.AddEntry(hreal_j1_mass,"Real","l")
        lst_j1_mass.SetTextAlign(13)

        lst_j1_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_j1_pt.AddEntry(htt_j1_pt,"TT","f")
        lst_j1_pt.AddEntry(hdy_j1_pt,"DY","f")
        lst_j1_pt.AddEntry(hwj_j1_pt,"WJ","f")
        lst_j1_pt.AddEntry(hst_j1_pt,"ST","f")
        lst_j1_pt.AddEntry(hsig_j1_pt,"Signal","l")
        lst_j1_pt.AddEntry(hreal_j1_pt,"Real","l")
        lst_j1_pt.SetTextAlign(13)

        lst_j1_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_j1_eta.AddEntry(htt_j1_eta,"TT","f")
        lst_j1_eta.AddEntry(hdy_j1_eta,"DY","f")
        lst_j1_eta.AddEntry(hwj_j1_eta,"WJ","f")
        lst_j1_eta.AddEntry(hst_j1_eta,"ST","f")
        lst_j1_eta.AddEntry(hsig_j1_eta,"Signal","l")
        lst_j1_eta.AddEntry(hreal_j1_eta,"Real","l")
        lst_j1_eta.SetTextAlign(13)

        lst_j1_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_j1_phi.AddEntry(htt_j1_phi,"TT","f")
        lst_j1_phi.AddEntry(hdy_j1_phi,"DY","f")
        lst_j1_phi.AddEntry(hwj_j1_phi,"WJ","f")
        lst_j1_phi.AddEntry(hst_j1_phi,"ST","f")
        lst_j1_phi.AddEntry(hsig_j1_phi,"Signal","l")
        lst_j1_phi.AddEntry(hreal_j1_phi,"Real","l")
        lst_j1_phi.SetTextAlign(13)

	lst_j1_btag = TLegend(0.80,0.70,0.95,0.90)
        lst_j1_btag.AddEntry(htt_j1_btag,"TT","f")
        lst_j1_btag.AddEntry(hdy_j1_btag,"DY","f")
        lst_j1_btag.AddEntry(hwj_j1_btag,"WJ","f")
        lst_j1_btag.AddEntry(hst_j1_btag,"ST","f")
        lst_j1_btag.AddEntry(hsig_j1_btag,"Signal","l")
        lst_j1_btag.AddEntry(hreal_j1_btag,"Real","l")
        lst_j1_btag.SetTextAlign(13)


	lst_j2_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_mass.AddEntry(htt_j2_mass,"TT","f")
        lst_j2_mass.AddEntry(hdy_j2_mass,"DY","f")
        lst_j2_mass.AddEntry(hwj_j2_mass,"WJ","f")
        lst_j2_mass.AddEntry(hst_j2_mass,"ST","f")
        lst_j2_mass.AddEntry(hsig_j2_mass,"Signal","l")
        lst_j2_mass.AddEntry(hreal_j2_mass,"Real","l")
        lst_j2_mass.SetTextAlign(13)

        lst_j2_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_pt.AddEntry(htt_j2_pt,"TT","f")
        lst_j2_pt.AddEntry(hdy_j2_pt,"DY","f")
        lst_j2_pt.AddEntry(hwj_j2_pt,"WJ","f")
        lst_j2_pt.AddEntry(hst_j2_pt,"ST","f")
        lst_j2_pt.AddEntry(hsig_j2_pt,"Signal","l")
        lst_j2_pt.AddEntry(hreal_j2_pt,"Real","l")
        lst_j2_pt.SetTextAlign(13)

        lst_j2_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_eta.AddEntry(htt_j2_eta,"TT","f")
        lst_j2_eta.AddEntry(hdy_j2_eta,"DY","f")
        lst_j2_eta.AddEntry(hwj_j2_eta,"WJ","f")
        lst_j2_eta.AddEntry(hst_j2_eta,"ST","f")
        lst_j2_eta.AddEntry(hsig_j2_eta,"Signal","l")
        lst_j2_eta.AddEntry(hreal_j2_eta,"Real","l")
        lst_j2_eta.SetTextAlign(13)

        lst_j2_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_phi.AddEntry(htt_j2_phi,"TT","f")
        lst_j2_phi.AddEntry(hdy_j2_phi,"DY","f")
        lst_j2_phi.AddEntry(hwj_j2_phi,"WJ","f")
        lst_j2_phi.AddEntry(hst_j2_phi,"ST","f")
        lst_j2_phi.AddEntry(hsig_j2_phi,"Signal","l")
        lst_j2_phi.AddEntry(hreal_j2_phi,"Real","l")
        lst_j2_phi.SetTextAlign(13)

        lst_j2_btag = TLegend(0.80,0.70,0.95,0.90)
        lst_j2_btag.AddEntry(htt_j2_btag,"TT","f")
        lst_j2_btag.AddEntry(hdy_j2_btag,"DY","f")
        lst_j2_btag.AddEntry(hwj_j2_btag,"WJ","f")
        lst_j2_btag.AddEntry(hst_j2_btag,"ST","f")
        lst_j2_btag.AddEntry(hsig_j2_btag,"Signal","l")
        lst_j2_btag.AddEntry(hreal_j2_btag,"Real","l")
        lst_j2_btag.SetTextAlign(13)


	lst_l1_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_mass.AddEntry(htt_l1_mass,"TT","f")
        lst_l1_mass.AddEntry(hdy_l1_mass,"DY","f")
        lst_l1_mass.AddEntry(hwj_l1_mass,"WJ","f")
        lst_l1_mass.AddEntry(hst_l1_mass,"ST","f")
        lst_l1_mass.AddEntry(hsig_l1_mass,"Signal","l")
        lst_l1_mass.AddEntry(hreal_l1_mass,"Real","l")
        lst_l1_mass.SetTextAlign(13)

        lst_l1_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_pt.AddEntry(htt_l1_pt,"TT","f")
        lst_l1_pt.AddEntry(hdy_l1_pt,"DY","f")
        lst_l1_pt.AddEntry(hwj_l1_pt,"WJ","f")
        lst_l1_pt.AddEntry(hst_l1_pt,"ST","f")
        lst_l1_pt.AddEntry(hsig_l1_pt,"Signal","l")
        lst_l1_pt.AddEntry(hreal_l1_pt,"Real","l")
        lst_l1_pt.SetTextAlign(13)

        lst_l1_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_eta.AddEntry(htt_l1_eta,"TT","f")
        lst_l1_eta.AddEntry(hdy_l1_eta,"DY","f")
        lst_l1_eta.AddEntry(hwj_l1_eta,"WJ","f")
        lst_l1_eta.AddEntry(hst_l1_eta,"ST","f")
        lst_l1_eta.AddEntry(hsig_l1_eta,"Signal","l")
        lst_l1_eta.AddEntry(hreal_l1_eta,"Real","l")
        lst_l1_eta.SetTextAlign(13)

        lst_l1_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_phi.AddEntry(htt_l1_phi,"TT","f")
        lst_l1_phi.AddEntry(hdy_l1_phi,"DY","f")
        lst_l1_phi.AddEntry(hwj_l1_phi,"WJ","f")
        lst_l1_phi.AddEntry(hst_l1_phi,"ST","f")
        lst_l1_phi.AddEntry(hsig_l1_phi,"Signal","l")
        lst_l1_phi.AddEntry(hreal_l1_phi,"Real","l")
        lst_l1_phi.SetTextAlign(13)

        lst_l1_charge = TLegend(0.80,0.70,0.95,0.90)
        lst_l1_charge.AddEntry(htt_l1_charge,"TT","f")
        lst_l1_charge.AddEntry(hdy_l1_charge,"DY","f")
        lst_l1_charge.AddEntry(hwj_l1_charge,"WJ","f")
        lst_l1_charge.AddEntry(hst_l1_charge,"ST","f")
        lst_l1_charge.AddEntry(hsig_l1_charge,"Signal","l")
        lst_l1_charge.AddEntry(hreal_l1_charge,"Real","l")
        lst_l1_charge.SetTextAlign(13)


	lst_l2_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_mass.AddEntry(htt_l2_mass,"TT","f")
        lst_l2_mass.AddEntry(hdy_l2_mass,"DY","f")
        lst_l2_mass.AddEntry(hwj_l2_mass,"WJ","f")
        lst_l2_mass.AddEntry(hst_l2_mass,"ST","f")
        lst_l2_mass.AddEntry(hsig_l2_mass,"Signal","l")
        lst_l2_mass.AddEntry(hreal_l2_mass,"Real","l")
        lst_l2_mass.SetTextAlign(13)

        lst_l2_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_pt.AddEntry(htt_l2_pt,"TT","f")
        lst_l2_pt.AddEntry(hdy_l2_pt,"DY","f")
        lst_l2_pt.AddEntry(hwj_l2_pt,"WJ","f")
        lst_l2_pt.AddEntry(hst_l2_pt,"ST","f")
        lst_l2_pt.AddEntry(hsig_l2_pt,"Signal","l")
        lst_l2_pt.AddEntry(hreal_l2_pt,"Real","l")
        lst_l2_pt.SetTextAlign(13)

        lst_l2_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_eta.AddEntry(htt_l2_eta,"TT","f")
        lst_l2_eta.AddEntry(hdy_l2_eta,"DY","f")
        lst_l2_eta.AddEntry(hwj_l2_eta,"WJ","f")
        lst_l2_eta.AddEntry(hst_l2_eta,"ST","f")
        lst_l2_eta.AddEntry(hsig_l2_eta,"Signal","l")
        lst_l2_eta.AddEntry(hreal_l2_eta,"Real","l")
        lst_l2_eta.SetTextAlign(13)

        lst_l2_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_phi.AddEntry(htt_l2_phi,"TT","f")
        lst_l2_phi.AddEntry(hdy_l2_phi,"DY","f")
        lst_l2_phi.AddEntry(hwj_l2_phi,"WJ","f")
        lst_l2_phi.AddEntry(hst_l2_phi,"ST","f")
        lst_l2_phi.AddEntry(hsig_l2_phi,"Signal","l")
        lst_l2_phi.AddEntry(hreal_l2_phi,"Real","l")
        lst_l2_phi.SetTextAlign(13)

        lst_l2_charge = TLegend(0.80,0.70,0.95,0.90)
        lst_l2_charge.AddEntry(htt_l2_charge,"TT","f")
        lst_l2_charge.AddEntry(hdy_l2_charge,"DY","f")
        lst_l2_charge.AddEntry(hwj_l2_charge,"WJ","f")
        lst_l2_charge.AddEntry(hst_l2_charge,"ST","f")
        lst_l2_charge.AddEntry(hsig_l2_charge,"Signal","l")
        lst_l2_charge.AddEntry(hreal_l2_charge,"Real","l")
        lst_l2_charge.SetTextAlign(13)


	lst_dr_ll = TLegend(0.80,0.70,0.95,0.90)
	lst_dr_ll.AddEntry(htt_dr_ll,"TT","f")
	lst_dr_ll.AddEntry(hdy_dr_ll,"DY","f")
	lst_dr_ll.AddEntry(hwj_dr_ll,"WJ","f")
	lst_dr_ll.AddEntry(hst_dr_ll,"ST","f")
	lst_dr_ll.AddEntry(hsig_dr_ll,"Signal","l")
	lst_dr_ll.AddEntry(hreal_dr_ll,"Real","l")
	lst_dr_ll.SetTextAlign(13)

	lst_dr_jj = TLegend(0.80,0.70,0.95,0.90)
        lst_dr_jj.AddEntry(htt_dr_jj,"TT","f")
        lst_dr_jj.AddEntry(hdy_dr_jj,"DY","f")
        lst_dr_jj.AddEntry(hwj_dr_jj,"WJ","f")
        lst_dr_jj.AddEntry(hst_dr_jj,"ST","f")
        lst_dr_jj.AddEntry(hsig_dr_jj,"Signal","l")
        lst_dr_jj.AddEntry(hreal_dr_jj,"Real","l")
        lst_dr_jj.SetTextAlign(13)


	lst_h1_mass = TLegend(0.80,0.70,0.95,0.90)
	lst_h1_mass.AddEntry(htt_h1_mass,"TT","f")
	lst_h1_mass.AddEntry(hdy_h1_mass,"DY","f")
	lst_h1_mass.AddEntry(hwj_h1_mass,"WJ","f")
	lst_h1_mass.AddEntry(hst_h1_mass,"ST","f")
	lst_h1_mass.AddEntry(hsig_h1_mass,"Signal","l")
	lst_h1_mass.AddEntry(hreal_h1_mass,"Real","l")
	lst_h1_mass.SetTextAlign(13)

	lst_h1_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_h1_pt.AddEntry(htt_h1_pt,"TT","f")
        lst_h1_pt.AddEntry(hdy_h1_pt,"DY","f")
        lst_h1_pt.AddEntry(hwj_h1_pt,"WJ","f")
        lst_h1_pt.AddEntry(hst_h1_pt,"ST","f")
        lst_h1_pt.AddEntry(hsig_h1_pt,"Signal","l")
        lst_h1_pt.AddEntry(hreal_h1_pt,"Real","l")
        lst_h1_pt.SetTextAlign(13)

	lst_h1_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_h1_eta.AddEntry(htt_h1_eta,"TT","f")
        lst_h1_eta.AddEntry(hdy_h1_eta,"DY","f")
        lst_h1_eta.AddEntry(hwj_h1_eta,"WJ","f")
        lst_h1_eta.AddEntry(hst_h1_eta,"ST","f")
        lst_h1_eta.AddEntry(hsig_h1_eta,"Signal","l")
        lst_h1_eta.AddEntry(hreal_h1_eta,"Real","l")
        lst_h1_eta.SetTextAlign(13)

	lst_h1_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_h1_phi.AddEntry(htt_h1_phi,"TT","f")
        lst_h1_phi.AddEntry(hdy_h1_phi,"DY","f")
        lst_h1_phi.AddEntry(hwj_h1_phi,"WJ","f")
        lst_h1_phi.AddEntry(hst_h1_phi,"ST","f")
        lst_h1_phi.AddEntry(hsig_h1_phi,"Signal","l")
        lst_h1_phi.AddEntry(hreal_h1_phi,"Real","l")
        lst_h1_phi.SetTextAlign(13)


	lst_h2_mass = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_mass.AddEntry(htt_h2_mass,"TT","f")
        lst_h2_mass.AddEntry(hdy_h2_mass,"DY","f")
        lst_h2_mass.AddEntry(hwj_h2_mass,"WJ","f")
        lst_h2_mass.AddEntry(hst_h2_mass,"ST","f")
        lst_h2_mass.AddEntry(hsig_h2_mass,"Signal","l")
        lst_h2_mass.AddEntry(hreal_h2_mass,"Real","l")
	lst_h2_mass.SetTextAlign(13)

	lst_h2_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_pt.AddEntry(htt_h2_pt,"TT","f")
        lst_h2_pt.AddEntry(hdy_h2_pt,"DY","f")
        lst_h2_pt.AddEntry(hwj_h2_pt,"WJ","f")
        lst_h2_pt.AddEntry(hst_h2_pt,"ST","f")
        lst_h2_pt.AddEntry(hsig_h2_pt,"Signal","l")
        lst_h2_pt.AddEntry(hreal_h2_pt,"Real","l")
        lst_h2_pt.SetTextAlign(13)

	lst_h2_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_eta.AddEntry(htt_h2_eta,"TT","f")
        lst_h2_eta.AddEntry(hdy_h2_eta,"DY","f")
        lst_h2_eta.AddEntry(hwj_h2_eta,"WJ","f")
        lst_h2_eta.AddEntry(hst_h2_eta,"ST","f")
        lst_h2_eta.AddEntry(hsig_h2_eta,"Signal","l")
        lst_h2_eta.AddEntry(hreal_h2_eta,"Real","l")
        lst_h2_eta.SetTextAlign(13)

	lst_h2_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_h2_phi.AddEntry(htt_h2_phi,"TT","f")
        lst_h2_phi.AddEntry(hdy_h2_phi,"DY","f")
        lst_h2_phi.AddEntry(hwj_h2_phi,"WJ","f")
        lst_h2_phi.AddEntry(hst_h2_phi,"ST","f")
        lst_h2_phi.AddEntry(hsig_h2_phi,"Signal","l")
        lst_h2_phi.AddEntry(hreal_h2_phi,"Real","l")
        lst_h2_phi.SetTextAlign(13)


	lst_hh_mass = TLegend(0.80,0.70,0.95,0.90)
	lst_hh_mass.AddEntry(htt_hh_mass,"TT","f")
	lst_hh_mass.AddEntry(hdy_hh_mass,"DY","f")
	lst_hh_mass.AddEntry(hwj_hh_mass,"WJ","f")
	lst_hh_mass.AddEntry(hst_hh_mass,"ST","f")
	lst_hh_mass.AddEntry(hsig_hh_mass,"Signal","l")
	lst_hh_mass.AddEntry(hreal_hh_mass,"Real","l")
	lst_hh_mass.SetTextAlign(13)

	lst_hh_pt = TLegend(0.80,0.70,0.95,0.90)
        lst_hh_pt.AddEntry(htt_hh_pt,"TT","f")
        lst_hh_pt.AddEntry(hdy_hh_pt,"DY","f")
        lst_hh_pt.AddEntry(hwj_hh_pt,"WJ","f")
        lst_hh_pt.AddEntry(hst_hh_pt,"ST","f")
        lst_hh_pt.AddEntry(hsig_hh_pt,"Signal","l")
        lst_hh_pt.AddEntry(hreal_hh_pt,"Real","l")
        lst_hh_pt.SetTextAlign(13)

	lst_hh_eta = TLegend(0.80,0.70,0.95,0.90)
        lst_hh_eta.AddEntry(htt_hh_eta,"TT","f")
        lst_hh_eta.AddEntry(hdy_hh_eta,"DY","f")
        lst_hh_eta.AddEntry(hwj_hh_eta,"WJ","f")
        lst_hh_eta.AddEntry(hst_hh_eta,"ST","f")
        lst_hh_eta.AddEntry(hsig_hh_eta,"Signal","l")
        lst_hh_eta.AddEntry(hreal_hh_eta,"Real","l")
        lst_hh_eta.SetTextAlign(13)

	lst_hh_phi = TLegend(0.80,0.70,0.95,0.90)
        lst_hh_phi.AddEntry(htt_hh_phi,"TT","f")
        lst_hh_phi.AddEntry(hdy_hh_phi,"DY","f")
        lst_hh_phi.AddEntry(hwj_hh_phi,"WJ","f")
        lst_hh_phi.AddEntry(hst_hh_phi,"ST","f")
        lst_hh_phi.AddEntry(hsig_hh_phi,"Signal","l")
        lst_hh_phi.AddEntry(hreal_hh_phi,"Real","l")
        lst_hh_phi.SetTextAlign(13)


#	dr = 'combine'
#	ff = TFile("%s.root"%(dr),"recreate")
#        ff.cd()
#        gDirectory.mkdir(dr)
#        ff.cd(dr)
#
#	stack_j1_mass.Write()
#        stack_j1_pt.Write()
#        stack_j1_eta.Write()
#        stack_j1_phi.Write()
#        stack_j1_btag.Write()
#
#	stack_j2_mass.Write()
#        stack_j2_pt.Write()
#        stack_j2_eta.Write()
#        stack_j2_phi.Write()
#        stack_j2_btag.Write()
#
#	stack_l1_mass.Write()
#        stack_l1_pt.Write()
#        stack_l1_eta.Write()
#        stack_l1_phi.Write()
#        stack_l1_charge.Write()
#
#        stack_l2_mass.Write()
#        stack_l2_pt.Write()
#        stack_l2_eta.Write()
#        stack_l2_phi.Write()
#        stack_l2_charge.Write()
#
#	stack_dr_ll.Write()
#	stack_dr_jj.Write()
#
#	stack_h1_mass.Write()
#        stack_h1_pt.Write()
#        stack_h1_eta.Write()
#        stack_h1_phi.Write()
#
#        stack_h2_mass.Write()
#        stack_h2_pt.Write()
#        stack_h2_eta.Write()
#        stack_h2_phi.Write()
#
#        stack_hh_mass.Write()
#        stack_hh_pt.Write()
#        stack_hh_eta.Write()
#        stack_hh_phi.Write()
#
#	ff.Close()

	print("Real data yield is as follows...:")
	print(rh1_yield)
	print(rh2_yield)
	print(rhh_yield)

	print("MC data yield is as follows...:")
	print(h1_yield)
        print(h2_yield)
        print(hh_yield)

	c1 = TCanvas("c1","Stacked Histograms of Basic Kinematics hh", 900,600)

	c1.Divide(2,2)
        c1.cd(1)
	stack_hh_mass.Draw("hist")
#	stack_hh_mass.SetMinimum(0)
#	stack_hh_mass.SetMaximum(1600)
	stack_hh_mass.GetXaxis().SetTitle("Invariant Mass[GeV]")
	stack_hh_mass.GetYaxis().SetTitle("Events/10GeV")
	hsig_hh_mass.SetLineColor(kRed)
	hsig_hh_mass.Draw("hist same")
	hreal_hh_mass.SetLineColor(kBlack)
	hreal_hh_mass.Draw("hist e1 same")
#	hreal_hh_mass.SetStats(0)
#	err_hh_mass.SetFillColor(kBlue)
#	err_hh_mass.SetFillStyle(3004)
#	err_hh_mass.Draw("E2 same")
	lst_hh_mass.Draw()

	c1.cd(2)
	stack_hh_pt.Draw("hist")
#	stack_hh_pt.SetMinimum(0)
#        stack_hh_pt.SetMaximum(2000)
        stack_hh_pt.GetXaxis().SetTitle("Transverse Momentum[GeV]")
        stack_hh_pt.GetYaxis().SetTitle("Events/10GeV")
        hsig_hh_pt.SetLineColor(kRed)
        hsig_hh_pt.Draw("hist same")
	hreal_hh_pt.SetLineColor(kBlack)
        hreal_hh_pt.Draw("hist e1 same")
#        err_hh_pt.SetFillColor(kBlue)
#        err_hh_pt.SetFillStyle(3004)
#        err_hh_pt.Draw("E2 same")
        lst_hh_pt.Draw()

	c1.cd(3)
	stack_hh_eta.Draw("hist")
#	stack_hh_eta.SetMinimum(0)
#        stack_hh_eta.SetMaximum(400)
        stack_hh_eta.GetXaxis().SetTitle("Eta")
        stack_hh_eta.GetYaxis().SetTitle("Events/.17")
        hsig_hh_eta.SetLineColor(kRed)
        hsig_hh_eta.Draw("hist same")
	hreal_hh_eta.SetLineColor(kBlack)
	hreal_hh_eta.Draw("hist e1 same")
#        err_hh_eta.SetFillColor(kBlue)
#        err_hh_eta.SetFillStyle(3004)
#        err_hh_eta.Draw("E2 same")
        lst_hh_eta.Draw()

	c1.cd(4)
	stack_hh_phi.Draw("hist")
	stack_hh_phi.SetMinimum(0)
        stack_hh_phi.SetMaximum(500)
        stack_hh_phi.GetXaxis().SetTitle("Phi")
        stack_hh_phi.GetYaxis().SetTitle("Events/.157rad")
        hsig_hh_phi.SetLineColor(kRed)
        hsig_hh_phi.Draw("hist same")
	hreal_hh_phi.SetLineColor(kBlack)
        hreal_hh_phi.Draw("hist e1 same")
#        err_hh_phi.SetFillColor(kBlue)
#        err_hh_phi.SetFillStyle(3004)
#        err_hh_phi.Draw("E2 same")
        lst_hh_phi.Draw()

        c1.SetGrid()
	c1.Draw()
	c1.SaveAs("hh.png")

	input("Press Enter to continue...")


if __name__ == '__main__':
        main()
