from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
from scipy.optimize import minimize, NonlinearConstraint
import matplotlib.pyplot as plt


def defineHTness():


    h_mass,t_mass,w_mass,z_mass = 125.18,172.26,80.379,91.1876

    sig = pd.read_csv('htness_hh.csv')

    sig = sig.values.tolist()

    columns_name = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi',
                    'jet2_mass','jet2_pt','jet2_eta','jet2_phi',
                    'lep1_mass','lep1_pt','lep1_eta','lep1_phi',
                    'lep2_mass','lep2_pt','lep2_eta','lep2_phi',
                    'dr_jj','dr_ll','met_pt','met_phi',
                    'higgsness','topness','target']

    signal = pd.DataFrame(sig, columns=columns_name)

    jet1_mass = signal.jet1_mass.tolist()
    jet1_pt = signal.jet1_pt.tolist()
    jet1_eta = signal.jet1_eta.tolist()
    jet1_phi = signal.jet1_phi.tolist()

    jet2_mass = signal.jet2_mass.tolist()
    jet2_pt = signal.jet2_pt.tolist()
    jet2_eta = signal.jet2_eta.tolist()
    jet2_phi = signal.jet2_phi.tolist()

    lep1_mass = signal.lep1_mass.tolist()
    lep1_pt = signal.lep1_pt.tolist()
    lep1_eta = signal.lep1_eta.tolist()
    lep1_phi = signal.lep1_phi.tolist()

    lep2_mass = signal.lep2_mass.tolist()
    lep2_pt = signal.lep2_pt.tolist()
    lep2_eta = signal.lep2_eta.tolist()
    lep2_phi = signal.lep2_phi.tolist()


    # ========================================
    # Draw Objective Function Plane
    lep1_p4 = TLorentzVector()
    lep2_p4 = TLorentzVector()

    lep1_p4.SetPtEtaPhiM(lep1_pt[1],lep1_eta[1],lep1_phi[1],lep1_mass[1])
    lep2_p4.SetPtEtaPhiM(lep2_pt[1],lep2_eta[1],lep2_phi[1],lep2_mass[1])

    jet1_p4 = TLorentzVector()
    jet2_p4 = TLorentzVector()

    jet1_p4.SetPtEtaPhiM(jet1_pt[1],jet1_eta[1],jet1_phi[1],jet1_mass[1])
    jet2_p4.SetPtEtaPhiM(jet1_pt[1],jet2_eta[1],jet2_phi[1],jet2_mass[1])

    nu1_p4 = TLorentzVector()
    nu2_p4 = TLorentzVector()

    def defineHiggsness(x):

            nu1_p4.SetPtEtaPhiM(x[0], x[1], x[2], 0)
            nu2_p4.SetPtEtaPhiM(x[3], x[4], x[5], 0)

            w1_p4 = lep1_p4 + nu1_p4
            w2_p4 = lep2_p4 + nu2_p4

            dilep_p4 = lep1_p4 + lep2_p4
            h1_p4 = w1_p4 + w2_p4

            sig_h = 34.0
            sig_on,sig_off = 16.7,22.7
            sig_lep = 16.1
            peak_off = (1/np.sqrt(3)) * np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*(h_mass**2)*(w_mass**2) + w_mass**4) )
            peak_dilep = 30

            inner = [((w1_p4.M()**2-w_mass**2)**2/sig_on**4 + (w2_p4.M()**2-peak_off**2)**2/sig_off**4), ((w2_p4.M()**2-w_mass**2)**2/sig_on**4 + (w1_p4.M()**2-peak_off**2)**2/sig_off**4)]

            higgsness = (h1_p4.M()**2-h_mass**2)**2/sig_h**4 + (dilep_p4.M()**2-peak_dilep**2)**2/sig_lep**4  + min(inner)

            return higgsness

    ini = np.array([0,0,0, 0,0,0])
    bnd = [(0,100),(-2.4,2.4),(-3.14,3.14), (0,100),(-2.4,2.4),(-3.14,3.14)]

    con_eq = lambda x: x[0] + x[3]
    cons = NonlinearConstraint(con_eq, met_pt, np.inf)
    
    opt = minimize(defineHiggsness, ini, constraints=cons, bounds=bnd, method='SLSQP', options={'ftol':20}) #'disp':True})

    x0,x1 = np.linspace(0,100,100),np.linspace(0,100,100)
    z = np.array([])
    for i in range(len(x0)):
            for j in range(len(x1)):
                    z = np.append(z,defineHiggsness([x0[i],opt.x[1],opt.x[2], x1[j],opt.x[4],opt.x[5]]))
    z.resize(100,100)
    x,y = np.meshgrid(x0,x1)

    min_x0,min_x1 = opt.x[0],opt.x[3]
    min_z = defineHiggsness([opt.x[0],opt.x[1],opt.x[2], opt.x[3],opt.x[4],opt.x[5]])

    fig = plt.figure(figsize=(15,10))

    from mpl_toolkits.mplot3d import axes3d
    ax = fig.add_subplot(131, projection='3d')
    ax.contour3D(x, y, z, 70, cmap='plasma')
    ax.scatter(min_x0, min_x1, min_z, marker='o',color='red',linewidth=10)
    ax.set_xlabel('$x_{0}$')
    ax.set_ylabel('$x_{1}$')
    ax.set_zlabel('$f(x)$')

    ax = fig.add_subplot(132, projection='3d')
    ax.contour3D(x, y, z, 70, cmap='plasma')
    ax.scatter(min_x0, min_x1, min_z, marker='o',color='red',linewidth=10)
    ax.set_xlabel('$x_{0}$')
    ax.set_ylabel('$x_{1}$')
    ax.set_zlabel('$f(x)$')
    ax.axes.yaxis.set_ticklabels([])
    ax.view_init(0,80)

    ax = fig.add_subplot(133, projection='3d')
    ax.contour3D(x, y, z, 70, cmap='plasma')
    ax.scatter(min_x0, min_x1, min_z, marker='o',color='red',linewidth=10)
    ax.set_xlabel('$x_{0}$')
    ax.set_ylabel('$x_{1}$')
    ax.set_zlabel('$f(x)$')
    ax.axes.zaxis.set_ticklabels([])
    ax.view_init(89,-90)

    plt.show()


if __name__ == "__main__":
	defineHTness()
