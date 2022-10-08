from ROOT import *
import numpy as np
import glob
import sys


def main():

        flist = glob.glob("/xrootd/store/data/Run2016B_ver2/DoubleMuon/NANOAOD/Nano25Oct2019_ver2-v1/40000/*")

	ii = 0
	for fname in flist:

		runs,start,end = [],[],[]
                readRunLumi(runs,start,end)

                f = TFile(fname,"read")
                ll = f.Get("LuminosityBlocks")
		ll.GetEntry()

		run  = ll.run
		lumi = ll.luminosityBlock

		print(run, lumi)

		good_run, good_lumi = False,False 
                if run in runs:
			good_run = True
			for jj in range(len(start)):
				if (runs[jj] == run):
					if not (start[jj] > lumi or end[jj] < lumi):
						good_lumi = True

		if (good_run == True and good_lumi == True): print("good")
		else: print("bad")

		ii += 1

def readRunLumi(runs,start,end):

        fname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis/json/runlumi_2016.txt"
        f = open(fname,"read")

        lines = f.read().split('\n')
	lines.pop(-1)
        for line in lines:

                line = line.strip()
                sp = line.split(',')

                runs.append(int(sp[0]))
                start.append(int(sp[1]))
                end.append(int(sp[2]))

        return runs,start,end


if __name__ == "__main__":
        main()
