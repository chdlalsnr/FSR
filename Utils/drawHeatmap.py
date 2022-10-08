from ROOT import *
import numpy as np
import math
import glob
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def main():

	df = pd.read_csv("signal.csv")
	sel = df.loc[:, ['jet1_eta','jet1_btag','jet2_eta','jet2_phi','jet2_btag','lep1_eta','lep2_eta','dr_jj','dr_ll','H1_mass','H1_eta','H2_eta','H2_phi','HH_mass','HH_pt','HH_eta','HH_phi']]

	sns.heatmap(sel.corr(), linewidth=0.05, vmax=0.5, cmap=plt.cm.hot, linecolor='white', annot=False)

	plt.show()

if __name__ == "__main__":
        main()

