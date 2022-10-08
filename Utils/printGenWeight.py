from ROOT import *
import numpy as np
import glob
import sys

def main():

        flist = glob.glob("/xrootd/store/mc/RunIISummer16NanoAODv7/GluGluHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/*/*.root")

	for fname in flist:

                f = TFile(fname,"read")
                t = f.Get("Events")

                event_num = 50
		for i in range(event_num):

                        t.GetEntry(i)
			print(t.genWeight)


if __name__ == "__main__":
        main()
