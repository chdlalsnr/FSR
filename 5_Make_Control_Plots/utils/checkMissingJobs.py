from ROOT import *
import numpy as np
import glob
import re

def main():

	keys = 'signal,tt_2l2nu1,tt_2l2nu2,tt_2l2nu3,tt_2l2nu4,tt_2l2nu5,tt_2l2nu6,tt_2l2nu7,tt_2l2nu8,tt_semi1,tt_semi2,tt_semi3,tt_semi4,tt_semi5,tt_semi6,tt_semi7,tt_semi8,tt_had1,tt_had2,tt_had3,tt_had4,tt_had5,tt_had6,tt_had7,tt_had8,dy_0j1,dy_0j2,dy_0j3,dy_0j4,dy_0j5,dy_0j6,dy_1j1,dy_1j2,dy_1j3,dy_1j4,dy_1j5,dy_1j6,dy_2j1,dy_2j2,dy_2j3,dy_2j4,dy_2j5,dy_2j6,wj_1j1,wj_1j2,wj_2j1,wj_2j2,wj_2j3,wj_2j4,wj_3j1,wj_3j2,wj_3j3,wj_3j4,wj_4j1,wj_4j2,st_schannel_incl1,st_tW_anti1,st_tW_anti2,st_tW_anti3,st_tW_top1,st_tW_top2,st_tchannel_anti1,st_tchannel_anti2,st_tchannel_anti3,st_tchannel_anti4,st_tchannel_top1,st_tchannel_top2'

	keys = re.split(',',keys)

	flist = glob.glob('./*.out')

	for fname in flist:
		for key in keys:
			if key in fname: keys.remove(key)

	remained = ''
	for key in keys:
		remained += key+','

	print(remained[:-1])


if __name__ == '__main__':
        main()
