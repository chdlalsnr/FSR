from ROOT import *
import numpy as np
import json
import sys
import pickle
import glob
import string


def main():

        fname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt" 
        f = open(fname,"read")

#	data = f.read()
#	print(data)

	cont = []
	lines = f.read().split('"')
	for line in lines:

		line = line.strip()
		if (line == ''): continue
		if (line[0] == '2'):
			ss = line[0]+line[1]+line[2]+line[3]+line[4]+line[5]

		start,end,sp = 0,0,[]
		if (line[0] == ':'):
			for s in range(len(line)):
				if (line[s] == '['):
					start = s+2
					break
			for s in range(len(line)):
				if (line[s] == ']'): end = s-1

			cc = ''
			for i in range(start,end):
				cc += line[i]
			cc = cc.strip()

			if ('\n    ' not in cc):
				sp = cc.split('], [')

			if ('\n    ' in cc):
				sp = cc.split('\n    ')
				for j in range(len(sp)): 
					sp[j] = sp[j].strip('[')
					sp[j] = sp[j].strip('],')

		if (sp == []): continue
		for i in range(len(sp)):
			res = ss + ', ' + sp[i]
			cont.append(res)

	f.close()

	ff = open('/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/json/runlumi_2016.txt', 'w')
	for i in range(len(cont)):
		ff.write(cont[i]+'\n')

	ff.close()

        input("Press Enter to continue...")

if __name__ == "__main__":
        main()
