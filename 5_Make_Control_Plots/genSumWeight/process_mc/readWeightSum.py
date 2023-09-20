from ROOT import *
import sys
import glob
import string

def readWeightSum():

	flist = glob.glob("./*.txt")

	weightSum = 0
	for fname in flist:

		f = open(fname,"read")
		data = float(f.read())

		weightSum += data

	print(weightSum)

if __name__ == "__main__":
	readWeightSum()

		
