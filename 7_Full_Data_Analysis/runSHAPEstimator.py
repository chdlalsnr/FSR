# -*- coding: utf-8 -*-
#from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf

from xgboost import XGBRegressor, XGBClassifier

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression, SelectKBest

from sklearn.model_selection import GridSearchCV

import random

def main():

	target = 'target'
	feature = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet2_mass','jet2_pt','jet2_eta','jet2_phi','lep1_pt','lep1_eta','lep1_phi','lep2_pt','lep2_eta','lep2_phi','dr_jj','dr_ll','met_pt','met_phi','h1_mass','h1_pt','h1_eta','h1_phi','h2_mass','h2_pt','h2_eta','h2_phi','higgsness','topness']
#	feature = ['jet1_mass','jet1_pt','jet2_mass','jet2_pt','lep1_pt','lep2_pt','dr_ll','dr_jj','met_pt','h1_mass','h1_pt','h2_mass','h2_pt','h2_phi','higgsness','topness'] # SelectKBest

	df_hh = pd.read_csv("csv/htness_sig.csv")
	df_tt = pd.read_csv("csv/htness_bkg.csv")

	df_hh = df_hh.sample(frac=1)
	df_tt = df_tt.sample(frac=1)

	data = pd.concat([df_hh[:15000], df_tt[:15000]])

	train,valid = train_test_split(data, test_size=0.3, random_state=62)

	X_train = train[feature]
	y_train = train[target]
	X_valid = valid[feature]
	y_valid = valid[target]

	scaler = StandardScaler()

	X_train = scaler.fit_transform(X_train)
	X_valid = scaler.transform(X_valid)

	model = XGBClassifier(n_estimators=400,
                        learning_rate=0.1,
                        max_depth=6)

	model.fit(X_train, y_train)

	y_pred = model.predict(X_valid)
	y_real = y_valid.values.tolist()

	# SHAP
	import shap

	shap.initjs()
	explainer = shap.Explainer(model, X_valid, feature_names=feature, seed=62)
	shap_values = explainer(X_valid[:500])

#	shap.plots.bar(shap_values, max_display=15)
	shap.summary_plot(shap_values)

	input("Press Enter to continue...")


if __name__ == "__main__":
        main() 
