# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns


def main():

        target = 'target'
        feature = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet2_mass','jet2_pt','jet2_eta','jet2_phi','lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep2_mass','lep2_pt','lep2_eta','lep2_phi','dr_jj','dr_ll','met_pt','met_phi','h1_mass','h1_pt','h1_eta','h1_phi','h2_mass','h2_pt','h2_eta','h2_phi','higgsness','topness']

        df_hh = pd.read_csv("htness_sig.csv")
        df_tt = pd.read_csv("htness_bkg.csv")

        data = pd.concat([df_hh, df_tt])

        sel = data.loc[:,feature]

        plt.figure(figsize=(12,8))

        sns.heatmap(sel.corr(),linewidth=0.05,vmax=0.5,cmap=plt.cm.hot,linecolor='white',annot=True,annot_kws={"size": 5})

        plt.title("Confusion matrix between contributing features")
        plt.xticks(range(30),feature)
        plt.xticks(fontsize=8)
        plt.xticks(rotation=45)
        plt.yticks(range(30),feature)
        plt.yticks(fontsize=8)

        plt.show()


if __name__ == "__main__":
        main()

