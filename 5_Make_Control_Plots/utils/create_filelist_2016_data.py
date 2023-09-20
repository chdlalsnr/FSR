#!/usr/bin/env python
import json
import glob
import pickle

with open("dataset_2016_data.json","r") as js:

	jj= json.load(js)
	datasets= {}
        prefix = "root://cms-xrdr.private.lo:2094/"
	for key in jj.keys():
		glob_cmd = jj[key]['dataset']
		flist = glob.glob(glob_cmd)
		xrd_flist = map(lambda x: prefix+x.replace("/xrootd","/xrd"), flist)

                with open("dataset_2016_data/%s.dat"%(key),"wb") as output:
			pickle.dump(xrd_flist, output)

signal_filelist = pickle.load(open("dataset_2016_data/run2016B_mumu1.dat","rb"))
print(signal_filelist)
	
