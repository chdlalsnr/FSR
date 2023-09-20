# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from ROOT import *

from xgboost import XGBRegressor, XGBClassifier

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression, SelectKBest

from sklearn.model_selection import GridSearchCV


def main():

	target = 'target'
#	feature = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet2_mass','jet2_pt','jet2_eta','jet2_phi','lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep2_mass','lep2_pt','lep2_eta','lep2_phi','dr_jj','dr_ll','met_pt','met_phi','h1_mass','h1_pt','h1_eta','h1_phi','h2_mass','h2_pt','h2_eta','h2_phi','higgsness','topness']
	feature = ['jet1_mass','jet1_pt','jet2_mass','jet2_pt','lep1_pt','lep2_pt','dr_jj','dr_ll','met_pt','h1_mass','h1_pt','h2_mass','h2_pt','higgsness','topness']

	df_hh = pd.read_csv("results/htness_sig.csv")
	df_tt = pd.read_csv("results/htness_bkg.csv")

	data = pd.concat([df_hh, df_tt])

	train,valid = train_test_split(data, test_size=0.3, random_state=42)

	X_train = train[feature]
	y_train = train[target]
	X_valid = valid[feature]
	y_valid = valid[target]

	scaler = StandardScaler()

	X_train = scaler.fit_transform(X_train)
	X_valid = scaler.fit_transform(X_valid)

#	estimator = XGBRegressor()
#
#	param_grid = {  'n_estimators' : [100,300,500],
#			'learning_rate': [0.1,0.2,0.3],
#			'max_depth'    : [3,9,15] }
#
#	grid_search = GridSearchCV( estimator=estimator,
#					param_grid=param_grid,
#					n_jobs=32,
#					verbose=2 )
#
#	grid_search.fit(X_train,y_train)
#
#	print(grid_search.best_params_)
	#{'learning_rate': 0.1, 'max_depth': 9, 'n_estimators': 500}

	model = XGBRegressor(n_estimators = 500,
			learning_rate = 0.1,
			max_depth = 9)

	model.fit(X_train, y_train)

	y_pred = model.predict(X_valid)
	y_real = y_valid.values.tolist()

	# ======================== Visualization Part ========================

	c1 = TCanvas("c1","XGB Regression Model Prediction_DL", 900,600)

	prev_hh = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)
	prev_tt = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)

	for i in range(len(y_real)):
		if y_real[i] == 1:
			prev_hh.Fill(y_pred.item(i))
		if y_real[i] == 0:
			prev_tt.Fill(y_pred.item(i))

	prev_hh.Scale(1/prev_hh.GetEntries())
	prev_tt.Scale(1/prev_tt.GetEntries())

	hist_hh = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)
	hist_tt = TH1D("h1","XGB Regression Prediction Score_DL", 40,0,1.01)

	for i in range(1,hist_hh.GetNbinsX()+1):
		hist_hh.SetBinContent(i,prev_hh.GetBinContent(i))
		hist_tt.SetBinContent(i,prev_tt.GetBinContent(i))

	hist_hh.SetLineColor(2)
	hist_hh.SetFillColor(2)
	hist_hh.SetFillStyle(3004)
	hist_tt.SetFillColor(4)
	hist_tt.SetFillStyle(3005)

	hist_tt.GetXaxis().SetTitle("Prediction Score")
	hist_tt.GetYaxis().SetTitle("Entry")

	c1.cd()
	hist_hh.Draw()
	hist_tt.Draw("same")

	pred_ld = y_pred.flatten()
	y_pred = np.where(pred_ld > 0.5, 1, 0)

	print(accuracy_score(y_real, y_pred))


#	from sklearn.metrics import roc_curve
#	y_pred = model.predict(X_valid).ravel()
#	fpr,tpr,threshold = roc_curve(y_valid, y_pred)
#
#	from sklearn.metrics import auc
#	auc1 = auc(fpr, tpr)
#
#	plt.figure(figsize=(9,6))
#
#	# Plot ROC curve
#	plt.plot([0, 1],[0, 1], 'k--',alpha=0.6)
#	plt.plot(fpr, tpr, c="red", label='AUC = {:.3f})'.format(auc1))
#
#	plt.xlabel('False Positive Rate')
#	plt.ylabel('True Positive Rate')
#	plt.title('XGB Classifier ROC Curve')
#	plt.grid(True, axis='y', alpha=0.5, linestyle='--')
#
#	plt.legend(loc='best')
#
#	plt.show()

#	# Draw accuracy and validation loss
#	y_vloss = history.history['val_loss']
#	y_acc  = history.history['acc']
#
#	x_len = np.arange(len(y_acc))
#	plt.plot(x_len, y_vloss, "o", c="red",  markersize=3, label="Validation Loss")
#	plt.plot(x_len, y_vloss, c="red", alpha=0.5, linewidth=1)
#	plt.plot(x_len, y_acc  , "o", c="blue", markersize=3, label="Accuracy")
#	plt.plot(x_len, y_acc, c="blue", alpha=0.5, linewidth=1)
#	plt.legend()
#
#	plt.title("Test Model Accuracy and Validation Loss")
#	plt.xlabel("Epoch")
#	plt.ylabel("Accuracy")
#	plt.grid(True, axis='y', alpha=0.5, linestyle='--')
#
#	plt.show()

#	import eli5
#	from eli5.sklearn import PermutationImportance
#
#	perm = PermutationImportance(model, scoring="f1", random_state=32).fit(X_valid, y_valid)
#	ww = eli5.show_weights(perm, top=15, feature_names=feature)
#
#	with open('feature_importance.html','wb') as f:
#		f.write(ww.data.encode("UTF-8"))


	input("Press Enter to continue...")


if __name__ == "__main__":
        main() 
