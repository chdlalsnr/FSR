# -*- coding: utf-8 -*-
from ROOT import *
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
#	feature = ['jet1_mass','jet1_pt','jet2_mass','jet2_pt','lep1_pt','lep2_pt','dr_ll','dr_jj','met_pt','h1_mass','h1_pt','h2_mass','h2_pt','h2_phi','higgsness','topness']

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

#	estimator = XGBRegressor()
#
#	param_grid = {  'n_estimators' : [100,200,300,400],
#			'learning_rate': [0.1,0.2,0.3,0.4],
#			'max_depth'    : [3,6,9,12] }
#
#	grid_search = GridSearchCV( estimator=estimator,
#					param_grid=param_grid,
#					n_jobs=16,
#					scoring="roc_auc" )
#
#	grid_search.fit(X_train,y_train)
#
#	print(grid_search.best_params_)
	# {'n_estimators': 400, 'learning_rate': 0.1, 'max_depth': 6}

	model = XGBRegressor(n_estimators=400,
			learning_rate=0.1,
			max_depth=6)

	model.fit(X_train, y_train)

	y_pred = model.predict(X_valid)
	y_real = y_valid.values.tolist()

	# ======================== Visualization Part ========================
	# Learning curve
#	model.fit(X_train, y_train, 
#		eval_set=[(X_train, y_train),(X_valid, y_valid)], 
#		early_stopping_rounds=20)
#
#	results = model.evals_result()
#
#	plt.figure(figsize=(9,6))
#	plt.plot(results["validation_0"]["rmse"], label="Training loss")
#	plt.plot(results["validation_1"]["rmse"], label="Validation loss")
#	plt.xlabel("Number of trees")
#	plt.ylabel("Loss")
#	plt.legend(loc='best')
#	plt.title("Learning Curve of XGB Classifier")
#
#	plt.grid(True, axis='y', alpha=0.5, linestyle='--')
#
#	plt.show()

	# XGBRegressor
#	c1 = TCanvas("c1","XGB Regression Model Prediction_DL", 900,600)
#
#	prev_hh = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)
#	prev_tt = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)
#
#	for i in range(len(y_real)):
#		if y_real[i] == 1:
#			prev_hh.Fill(y_pred.item(i))
#		if y_real[i] == 0:
#			prev_tt.Fill(y_pred.item(i))
#
#	prev_hh.Scale(1/prev_hh.GetEntries())
#	prev_tt.Scale(1/prev_tt.GetEntries())
#
#	hist_hh = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)
#	hist_tt = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)
#
#	for i in range(1,hist_hh.GetNbinsX()+1):
#		hist_hh.SetBinContent(i,prev_hh.GetBinContent(i))
#		hist_tt.SetBinContent(i,prev_tt.GetBinContent(i))
#
#	hist_hh.GetYaxis().SetRangeUser(0, 0.3)
#
#	hist_hh.SetLineColor(2)
#	hist_hh.SetFillColor(2)
#	hist_hh.SetFillStyle(3004)
#	hist_tt.SetFillColor(4)
#	hist_tt.SetFillStyle(3005)
#
#	hist_tt.GetXaxis().SetTitle("Prediction Score")
#	hist_tt.GetYaxis().SetTitle("Entry")
#
#	c1.cd()
#	hist_hh.Draw()
#	hist_tt.Draw("same")
#
#	pred_ld = y_pred.flatten()
#	y_pred = np.where(pred_ld > 0.5, 1, 0)
#
#	print(accuracy_score(y_real, y_pred))

	from sklearn.metrics import roc_curve
	y_pred = model.predict(X_valid).ravel()
	fpr,tpr,threshold = roc_curve(y_valid, y_pred)

	from sklearn.metrics import auc
	auc = auc(fpr, tpr)
	print(auc)

	f = open("xgb.csv",'a')
	f.write(str(auc)+",")

	f.close()

#	plt.figure(figsize=(9,6))
#	# Plot ROC curve
#	plt.plot([0, 1],[0, 1], 'k--',alpha=0.6)
#	plt.plot(fpr, tpr, c="red", label='AUC = {:.3f})'.format(auc))
#
#	plt.xlabel('False Positive Rate')
#	plt.ylabel('True Positive Rate')
#	plt.title('XGB Classifier ROC Curve')
#
#	plt.grid(True, axis='y', alpha=0.5, linestyle='--')
#
#	plt.legend(loc='best')
#        plt.show()

	# XGBClassifier
#	import eli5
#	from eli5.sklearn import PermutationImportance
#
#	perm = PermutationImportance(model, scoring="f1", random_state=62).fit(X_valid, y_valid)
#	ww = eli5.show_weights(perm, top=15, feature_names=feature)
#
#	with open('feature_importance.html','wb') as f:
#		f.write(ww.data.encode("UTF-8"))
#
#	input("Press Enter to continue...")


if __name__ == "__main__":
        main() 
