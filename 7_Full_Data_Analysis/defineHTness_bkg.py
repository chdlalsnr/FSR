from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
from scipy.optimize import minimize,NonlinearConstraint


def defineHTness():

	h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876

	bkg = pd.read_csv('csv/background.csv')
	bkg = bkg.values.tolist()

	columns_name = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet1_btag',
                        'jet2_mass','jet2_pt','jet2_eta','jet2_phi','jet2_btag',
                        'lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep1_charge',
                        'lep2_mass','lep2_pt','lep2_eta','lep2_phi','lep2_charge',
                        'dr_ll','dr_jj','met_pt','met_phi',
                        'H1_mass','H1_pt','H1_eta','H1_phi',
                        'H2_mass','H2_pt','H2_eta','H2_phi',
                        'HH_mass','HH_pt','HH_eta','HH_phi','genWeight']

	background = pd.DataFrame(bkg, columns=columns_name)

	jet1_mass = background.jet1_mass.tolist()
	jet1_pt = background.jet1_pt.tolist()
	jet1_eta = background.jet1_eta.tolist()
	jet1_phi = background.jet1_phi.tolist()

	jet2_mass = background.jet2_mass.tolist()
        jet2_pt = background.jet2_pt.tolist()
        jet2_eta = background.jet2_eta.tolist()
        jet2_phi = background.jet2_phi.tolist()

	lep1_mass = background.lep1_mass.tolist()
        lep1_pt = background.lep1_pt.tolist()
        lep1_eta = background.lep1_eta.tolist()
        lep1_phi = background.lep1_phi.tolist()

        lep2_mass = background.lep2_mass.tolist()
        lep2_pt = background.lep2_pt.tolist()
        lep2_eta = background.lep2_eta.tolist()
        lep2_phi = background.lep2_phi.tolist()

	dr_ll = background.dr_ll.tolist()
        dr_jj = background.dr_jj.tolist()

	met_pt = background.met_pt.tolist()
	met_phi = background.met_phi.tolist()

	h1_mass = background.H1_mass.tolist()
	h1_pt = background.H1_pt.tolist()
	h1_eta = background.H1_eta.tolist()
	h1_phi = background.H1_phi.tolist()

	h2_mass = background.H2_mass.tolist()
	h2_pt = background.H2_pt.tolist()
	h2_eta = background.H2_eta.tolist()
	h2_phi = background.H2_phi.tolist()

	data = []
	# =================================================================
	for i in range(len(jet1_mass)):

                lep1_p4 = TLorentzVector()
                lep2_p4 = TLorentzVector()

                lep1_p4.SetPtEtaPhiM(lep1_pt[i],lep1_eta[i],lep1_phi[i],lep1_mass[i])
                lep2_p4.SetPtEtaPhiM(lep2_pt[i],lep2_eta[i],lep2_phi[i],lep2_mass[i])

                jet1_p4 = TLorentzVector()
                jet2_p4 = TLorentzVector()

                jet1_p4.SetPtEtaPhiM(jet1_pt[i],jet1_eta[i],jet1_phi[i],jet1_mass[i])
                jet2_p4.SetPtEtaPhiM(jet2_pt[i],jet2_eta[i],jet2_phi[i],jet2_mass[i])

                nu1_p4 = TLorentzVector()
                nu2_p4 = TLorentzVector()

                def defineHiggsness(x):

                        nu1_p4.SetPtEtaPhiM(x[0], x[1], x[2], 0)
                        nu2_p4.SetPtEtaPhiM(x[3], x[4], x[5], 0)

                        w1_p4 = lep1_p4 + nu1_p4
                        w2_p4 = lep2_p4 + nu2_p4

                        h1_p4 = w1_p4 + w2_p4
                        dilep_p4 = lep1_p4 + lep2_p4

			sig_h = 34.9
                        sig_on,sig_off = 18.5,17.6
                        sig_lep = 16.0
                        peak_dilep = 30
                        peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )

                        inner = [((w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4), ((w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4)]

                        higgsness = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (dilep_p4.M()**2-peak_dilep**2)**2/sig_lep**4  + min(inner)

                        return higgsness

                ini = np.array([0,0,0, 0,0,0])
                bnd = [(0,100),(-2.4,2.4),(-3.14,3.14), (0,100),(-2.4,2.4),(-3.14,3.14)]

                con_eq = lambda x: x[0] + x[3]
                cons = NonlinearConstraint(con_eq, met_pt[i], np.inf)

                opt = minimize(defineHiggsness, ini, constraints=cons, bounds=bnd, method='SLSQP', options={'ftol':20}) #'disp':True})

		if not opt.success: continue

                nu1_p4.SetPtEtaPhiM(opt.x[0], opt.x[1], opt.x[2], 0)
                nu2_p4.SetPtEtaPhiM(opt.x[3], opt.x[4], opt.x[5], 0)

                w1_p4 = lep1_p4 + nu1_p4
                w2_p4 = lep2_p4 + nu2_p4

                h1_p4 = w1_p4 + w2_p4

                def defineTopness():

                        chi2 = []

                        # (t1, t2) and (t5, t6) have anti-particle relations
                        t1_p4 = w1_p4 + jet1_p4
                        t2_p4 = w1_p4 + jet2_p4
                        t3_p4 = w2_p4 + jet1_p4
                        t4_p4 = w2_p4 + jet2_p4

                        t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                        w_cand = [w1_p4.M(),w2_p4.M()]

			sig_t = 41.0
                        sig_w = 19.0

                        chi2.append( ((t_cand[0]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[3]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )
                        chi2.append( ((t_cand[1]**2-t_mass**2)**2/sig_t**4) + ((w_cand[0]**2-w_mass**2)**2/sig_w**4) + ((t_cand[2]**2-t_mass**2)**2/sig_t**4) + ((w_cand[1]**2-w_mass**2)**2/sig_w**4) )

                        topness = min(chi2)

                        return topness

		# =================================================================
		# Save CSV file
		data.append([jet1_mass[i],jet1_pt[i],jet1_eta[i],jet1_phi[i],
			     jet2_mass[i],jet2_pt[i],jet2_eta[i],jet2_phi[i],
			     lep1_mass[i],lep1_pt[i],lep1_eta[i],lep1_phi[i],
			     lep2_mass[i],lep2_pt[i],lep2_eta[i],lep2_phi[i],
			     dr_ll[i],dr_jj[i],met_pt[i],met_phi[i],
			     h1_mass[i],h1_pt[i],h1_eta[i],h1_phi[i],
			     h2_mass[i],h2_pt[i],h2_eta[i],h2_phi[i],
			     opt.fun,defineTopness(),0])

        df = pd.DataFrame(data, columns=
                       ['jet1_mass','jet1_pt','jet1_eta','jet1_phi',
                        'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
                        'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
                        'lep2_mass','lep2_pt','lep2_eta','lep2_phi',
                        'dr_ll','dr_jj','met_pt','met_phi',
			'h1_mass','h1_pt','h1_eta','h1_phi',
			'h2_mass','h2_pt','h2_eta','h2_phi',
                        'higgsness','topness','target'])

        df.to_csv("csv/htness_bkg.csv", header=True, index=False)

	input("Press Enter to continue...")


if __name__ == "__main__":
	defineHTness()
