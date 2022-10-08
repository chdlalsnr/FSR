import pandas as pd
import glob
import os


def main():

	input_file  = "/cms/ldap_home/chdlalsnr/HH/CMSSW_12_0_0_pre2/src/HHAnalysis/HH/python/Analysis/4/csv" 

	flist = ['signal.csv','background.csv']

	dd = []
	for f in flist:
		df = pd.read_csv(f,index_col=0)
		dd.append(df)

	comb = pd.concat(dd, axis=0)
	comb.to_csv("total.csv")

#	tt_list = glob.glob(os.path.join(input_file,'tt*.csv'))
#	dy_list = glob.glob(os.path.join(input_file,'dy*.csv'))
#	wj_list  = glob.glob(os.path.join(input_file,'w*.csv'))
#
#	df_sig = pd.read_csv("signal.csv",index_col=0)
#
#	tt,dy,wj = [],[],[]
#	for f in tt_list:
#		df = pd.read_csv(f,index_col=0)
#		tt.append(df)
#
#	for f in dy_list:
#		df = pd.read_csv(f,index_col=0)
#                dy.append(df)
#
#	for f in wj_list:
#		df = pd.read_csv(f,index_col=0)
#                wj.append(df)
#
#	df_sig.loc[:,'flag'] = 0
#	df_sig.to_csv("signal+.csv")
#
#	comb_tt = pd.concat(tt, axis=0)
#	comb_tt.loc[:,'flag'] = 1
#        comb_tt.to_csv("tt.csv")
#
#	comb_dy = pd.concat(dy, axis=0)
#	comb_dy.loc[:,'flag'] = 2
#        comb_dy.to_csv("dy.csv")
#
#        comb_wj = pd.concat(wj, axis=0)
#	comb_wj.loc[:,'flag'] = 3
#        comb_wj.to_csv("wj.csv")

	input("Press Enter to continue...")


if __name__ == "__main__":
        main()

