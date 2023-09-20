from ROOT import *
import numpy as np
import json
import sys
import pickle
import glob

def main(index):

	with open('dataset_2016_data.json','r') as js:
		jj = json.load(js)

	keys = ['run2016B_mm','run2016C_mm','run2016D_mm','run2016E_mm','run2016F_mm','run2016G_mm','run2016H_mm',
		'run2016B_em','run2016C_em','run2016D_em','run2016E_em','run2016F_em','run2016G_em','run2016H_em',
		'run2016B_ee','run2016C_ee','run2016D_ee','run2016E_ee','run2016F_ee','run2016G_ee','run2016H_ee',
		'run2016B_mu','run2016C_mu','run2016D_mu','run2016E_mu','run2016F_mu','run2016G_mu','run2016H_mu',
		'run2016B_el','run2016C_el','run2016D_el','run2016E_el','run2016F_el','run2016G_el','run2016H_el']

	jsonArray = []
	for key in keys:
		jsonArray.append(jj[key]['dataset'])

	index = int(index)
	flist = glob.glob(jsonArray[index])
	key = keys[index]

	totalEventNum = 0
	for fname in flist:

		runs,start,end = [],[],[]
                readRunLumi(runs,start,end)

                f = TFile.Open(fname,"read")
                t = f.Get("Events")
                ll = f.Get("LuminosityBlocks")
                ll.GetEntry()

                run  = int(ll.run)
                lumi = int(ll.luminosityBlock)

                good_run,good_lumi = False,False
                if run in runs:
                        good_run = True
                        for jj in range(len(start)):
                                if (runs[jj] == run):
                                        if (start[jj] < lumi and end[jj] > lumi):
                                                good_lumi = True

                if (good_run != True or good_lumi != True): continue

		f = TFile.Open(fname,"read")
		t = f.Get("Events")
		t.GetEntry()

		totalEventNum += t.GetEntriesFast()

	xf = open("%s.txt"%(key),'w')
	xf.write("%s"%(totalEventNum))
	xf.close()

def readRunLumi(runs,start,end):

	fname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Make_Control_Plots/json/runlumi_2016.txt"
        f = open(fname,"read")

        cont = []
        lines = f.read().split('\n')
        for line in lines:

                line = line.strip()
                sp = line.split(',')

                cont.append(sp)
        cont.pop(-1)

        for i in range(len(cont)):

                runs.append(int(cont[i][0]))
                start.append(int(cont[i][1]))
                end.append(int(cont[i][2]))

        f.close()

        return runs, start, end


if __name__ == "__main__":
	main(sys.argv[1])

