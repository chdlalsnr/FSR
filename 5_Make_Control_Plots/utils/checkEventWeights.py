import numpy as np
import pandas as pd 
import glob
import sys


def main(key):


	df = pd.read_csv('%s.csv'%(key))

#	flist = glob.glob('%s.csv'%(key))
	genweightSum = 0

	genweight = sum(df.iloc[:,-1])
	genweightSum += genweight

	print(genweightSum)


if __name__ == "__main__":
        main(sys.argv[1])
