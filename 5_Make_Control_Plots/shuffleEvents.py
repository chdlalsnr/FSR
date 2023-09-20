# -*- coding: utf-8 -*-
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():

	columns_name = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet1_btag',
                        'jet2_mass','jet2_pt','jet2_eta','jet2_phi','jet2_btag',
                        'lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep1_charge',
                        'lep2_mass','lep2_pt','lep2_eta','lep2_phi','lep2_charge',
                        'dr_ll','dr_jj','met_pt','met_phi',
                        'H1_mass','H1_pt','H1_eta','H1_phi',
                        'H2_mass','H2_pt','H2_eta','H2_phi',
                        'HH_mass','HH_pt','HH_eta','HH_phi','genWeight']

        data = pd.read_csv("background.csv")
        data = data.sample(frac=1)
	data = data.values.tolist()

        df = pd.DataFrame(data, columns=columns_name)

        jet1_mass = df.jet1_mass.tolist()
        jet1_pt = df.jet1_pt.tolist()
        jet1_eta = df.jet1_eta.tolist()
        jet1_phi = df.jet1_phi.tolist()
	jet1_btag = df.jet1_btag.tolist()

        jet2_mass = df.jet2_mass.tolist()
        jet2_pt = df.jet2_pt.tolist()
        jet2_eta = df.jet2_eta.tolist()
        jet2_phi = df.jet2_phi.tolist()
	jet2_btag = df.jet2_btag.tolist()

        lep1_mass = df.lep1_mass.tolist()
        lep1_pt = df.lep1_pt.tolist()
        lep1_eta = df.lep1_eta.tolist()
        lep1_phi = df.lep1_phi.tolist()
	lep1_charge = df.lep1_charge.tolist()

        lep2_mass = df.lep2_mass.tolist()
        lep2_pt = df.lep2_pt.tolist()
        lep2_eta = df.lep2_eta.tolist()
        lep2_phi = df.lep2_phi.tolist()
	lep2_charge = df.lep2_charge.tolist()

        dr_ll = df.dr_ll.tolist()
        dr_jj = df.dr_jj.tolist()
        met_pt = df.met_pt.tolist()
        met_phi = df.met_phi.tolist()

        h1_mass = df.H1_mass.tolist()
        h1_pt = df.H1_pt.tolist()
        h1_eta = df.H1_eta.tolist()
        h1_phi = df.H1_phi.tolist()

        h2_mass = df.H2_mass.tolist()
        h2_pt = df.H2_pt.tolist()
        h2_eta = df.H2_eta.tolist()
        h2_phi = df.H2_phi.tolist()

	hh_mass = df.HH_mass.tolist()
        hh_pt = df.HH_pt.tolist()
        hh_eta = df.HH_eta.tolist()
        hh_phi = df.HH_phi.tolist()

	weight = df.genWeight.tolist()

	data = []
        for i in range(len(jet1_mass)):
                data.append([jet1_mass[i],jet1_pt[i],jet1_eta[i],jet1_phi[i],jet1_btag[i],
                             jet2_mass[i],jet2_pt[i],jet2_eta[i],jet2_phi[i],jet2_btag[i],
                             lep1_mass[i],lep1_pt[i],lep1_eta[i],lep1_phi[i],lep1_charge[i],
                             lep2_mass[i],lep2_pt[i],lep2_eta[i],lep2_phi[i],lep2_charge[i],
                             dr_ll[i],dr_jj[i],met_pt[i],met_phi[i],
                             h1_mass[i],h1_pt[i],h1_eta[i],h1_phi[i],
                             h2_mass[i],h2_pt[i],h2_eta[i],h2_phi[i],
                             hh_mass[i],hh_pt[i],hh_eta[i],hh_phi[i],weight[i]])

        df = pd.DataFrame(data, columns=columns_name)

        df.to_csv("shuffled.csv", header=False, index=False)

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()
