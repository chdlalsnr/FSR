from ROOT import *
import numpy as np
import json
import sys
import pickle
import glob
import string


def main():

        fname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/jec/Summer16_07Aug2017_V11_MC_L1FastJet_AK4PF.txt" 
        f = open(fname,"read")

	lines = f.read().split('\n')
	lines.pop(-1)

	cont = []
	for line in lines:

		if 'Jet' in line: continue
		word = line.split(' ')

		i = 0
		while True:
			try:
				if (word[i] == ''):
					word.pop(i)
					i -= 1
			except: break
			i += 1

		start = word[0]
		end   = word[1]

		alpha = word[-3]
		beta  = word[-2]
		gamma = word[-1]

		ss = start + "," + end + "," + alpha + "," + beta + "," + gamma
		ss = str(ss)
		cont.append(ss)

#	if i in range(len(cont)):
#
#		if (eta > cont[i][0]): continue
#
#		corr = 1 - (( cont[i][2] + rho*cont[i][3]*(1 + cont[i][4]*np.log10(pt)) )*area)/pt
#
#		pt_corr = pt*corr

	ff = open('/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/jec/jetl1corr_2016.txt', 'w')
        for i in range(len(cont)):
                ff.write(cont[i]+'\n')

	ff.close()
	f.close()

        input("Press Enter to continue...")


if __name__ == "__main__":
        main()
