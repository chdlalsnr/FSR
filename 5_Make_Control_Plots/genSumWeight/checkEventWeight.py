from ROOT import *
import numpy as np
import json
import sys
import pickle
import glob

def main(index):

	with open('dataset_2016.json','r') as js:
		jj = json.load(js)

	keys = ['signal','tt_2l2nu','tt_semi','tt_had','dy_m50','dy_0j','dy_1j','dy_2j','wj_jets','wj_1j','wj_2j','wj_3j','wj_4j','st_schannel_incl','st_tW_anti','st_tW_top','st_tchannel_anti','st_tchannel_top']

	jsonArray = []
	for key in keys:
		jsonArray.append(jj[key]['dataset'])

	index = int(index)
	flist = glob.glob(jsonArray[index])
	key = keys[index]

	genweightSum = 0
	for fname in flist:

		f = TFile.Open(fname,"read")
		t = f.Get("Events")
		r = f.Get("Runs")
		r.GetEntry()

		genweightSum += r.genEventSumw
#		genweightSum += t.GetEntriesFast()

	xf = open("%s.txt"%(key),'w')
	xf.write("%s"%(genweightSum))
	xf.close()


if __name__ == "__main__":
	main(sys.argv[1])

