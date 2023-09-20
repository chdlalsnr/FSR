from ROOT import *
import numpy as np
import json
import sys
import pickle
import glob
import string


def main(eta,pt,rho,area):

        fname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/jec/jetl1corr_2016.txt" 
        f = open(fname,"read")

	lines = f.read().split('\n')
	lines.pop(-1)

	jecl1 = []
	for line in lines:

		word = line.split(',')
		jecl1.append(word)

	f.close()

	if i in range(len(jecl1)):

		if (eta > jecl1[i][1]): continue

		corr = 1 - (( jecl1[i][2] + rho*jecl1[i][3]*(1 + jecl1[i][4]*np.log10(pt)) )*area)/pt

		pt_corr = pt*corr

	return pt_corr


        input("Press Enter to continue...")


if __name__ == "__main__":
        main()
