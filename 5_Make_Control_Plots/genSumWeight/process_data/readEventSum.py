from ROOT import *
import sys
import glob
import string

def readEventSum():

	flist = glob.glob("./*.txt")

	eventSum = 0
	for fname in flist:

		f = open(fname,"read")
		data = float(f.read())

		eventSum += data

	print(eventSum)

if __name__ == "__main__":
	readEventSum()

		
