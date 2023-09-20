from ROOT import *
import numpy as np
import json
import sys
import pickle
import glob
import string


def readRunLumi():
#def readRunLumi(runs,start,end):

	runs,start,end = [],[],[]
        fname = "/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Make_Control_Plots/json/runlumi_2017.txt"
        f = open(fname,"read")

        cont = []
        lines = f.read().split('\n')
        for line in lines:

                line = line.strip()
                sp = line.split(',')

                cont.append(sp)
	cont.pop(-1)

	print(cont)

	for i in range(len(cont)):

		runs.append(int(cont[i][0]))
		start.append(int(cont[i][1]))
		end.append(int(cont[i][2]))

	for i in range(len(runs)):
		print(runs[i],start[i],end[i])

	f.close()
        input("Press Enter to continue...")


if __name__ == "__main__":
        readRunLumi()
