#!/usr/bin/env python
import json
import glob
import pickle

with open("dataset_2016_el.json","r") as js:

	jj= json.load(js)
	datasets= {}
        prefix = "root://cms-xrdr.private.lo:2094/"
	subkey = { 'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z', 'aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al','am','an','ao','ap','aq','ar','as','at','au','av','aw','ax','ay','az', 'ba','bb','bc','bd','be','bf','bg','bh','bi','bj','bk','bl','bm','bn','bo','bp','bq','br','bs','bt','bu','bv','bw','bx','by','bz' }
	for key in jj.keys():
		glob_cmd = jj[key]['dataset']
                flist = glob.glob(glob_cmd)
		for skey in subkey:
			if len(flist) >= 2:
				fls = flist[:2]
				xrd_fls = map(lambda x: prefix+x.replace("/xrootd","/xrd"), fls)
				with open("dataset_2016_el/%s%s.dat"%(key,skey),"wb") as output:
					pickle.dump(xrd_fls, output)
				flist = [i for i in flist if i not in fls]
			if len(flist) < 2:
				fls = flist[:]
				xrd_fls = map(lambda x: prefix+x.replace("/xrootd","/xrd"), fls)
                                with open("dataset_2016_el/%s%s.dat"%(key,skey),"wb") as output:
					pickle.dump(xrd_fls, output)
				break

signal_filelist = pickle.load(open("dataset_2016_el/run2016B_el1aa.dat","rb"))
print(signal_filelist)
	
