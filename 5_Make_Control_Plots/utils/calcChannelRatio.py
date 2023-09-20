import numpy as np
import sys


def main(key):


	keys = ['tt_2l2nu','tt_semi','tt_had','st_schannel_incl','st_tchannel_top','st_tchannel_anti','st_tW_top','st_tW_anti','dy_0j','dy_1j','dy_2j']
	xsec = [88.4      ,366      ,378     ,3.36              ,136              ,81                ,35.8       ,35.8        ,4840   ,898    ,336]

	summ = sum(xsec)

	ratio = 0
	for i in range(len(keys)):
		if keys[i] == key:
			ratio = xsec[i]/summ


	print(ratio)
			

if __name__ == "__main__":
        main(sys.argv[1])
