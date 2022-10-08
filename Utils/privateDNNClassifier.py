from ROOT import *

from keras.models import Sequential, load_model
from keras.layers import Dense
from keras.callbacks import ModelCheckpoint, EarlyStopping

import numpy as np
import pandas as pd
import seaborn as sns
import tensorflow as tf
import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split

import pickle
import math
import glob
import sys
import os

from keras import backend as k

def main():

	df = pd.read_csv("total.csv")
	sel = df.loc[:, ['jet1_eta','jet1_btag','jet2_eta','jet2_phi','jet2_btag','lep1_eta','lep2_eta','dr_jj','dr_ll','H1_mass','H1_eta','H2_eta','H2_phi','HH_mass','HH_pt','HH_eta','HH_phi']]
	dataset = df.values
	flag = dataset[:,-1]

	x_train,x_test,y_train,y_test = train_test_split(sel,flag, test_size=0.5)

#	# TMVA BDT Training
#	from xgboost import XGBClassifier
#	bdt = XGBClassifier(max_depth=3, n_estimators=500)
#	w = np.ones(2)
#	bdt.fit(x_train, y_train, w)
#
#	ROOT.TMVA.Experimental.SaveXGBoost(bdt,"myBDT","tmva.root")

	model = Sequential()
	model.add(Dense(36, input_dim=17, activation='relu'))
	model.add(Dense(12, activation='relu'))
	model.add(Dense(6, activation='relu'))
	model.add(Dense(1, activation='sigmoid'))

	model.compile(loss='binary_crossentropy',
		      optimizer='adam',
		      metrics=['accuracy'])
#		      metrics={'acc':'accuracy', 'rec':recall, 'prec':precision, 'f1':f1score})
#			       tf.keras.metrics.Precision(name='prec'),
#		               tf.keras.metrics.Recall(name='rec')])
#			       tf.keras.metrics.FalsePositives(name='fp'),
#			       tf.keras.metrics.FalseNegatives(name='fn')])

	MODEL_DIR = './output/model/'
	if not os.path.exists(MODEL_DIR):
		os.mkdir(MODEL_DIR)

	modelpath="./output/model/{epoch:02d}-{val_loss:.4f}.hdf5"

	checkpointer = ModelCheckpoint(filepath=modelpath, monitor='val_loss', verbose=1, save_best_only=True)
	early_stopping_callback = EarlyStopping(monitor='val_loss', patience=200)

	#history = model.fit(sel, flag, validation_split=0.5, epochs=1000, batch_size=100, verbose=0, callbacks=[early_stopping_callback,checkpointer])
	history = model.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=1000, batch_size=100, verbose=0, callbacks=[early_stopping_callback,checkpointer])

	# Draw ROC curve
#	from sklearn.metrics import roc_curve
#	y_pred = model.predict(x_test).ravel()
#	fpr,tpr,threshold = roc_curve(y_test, y_pred)
#	from sklearn.metrics import auc
#	auc = auc(fpr, tpr)
#
#	# Plot ROC curve
#	plt.plot([0, 1],[0, 1], 'k--',alpha=0.6)
#	plt.plot(fpr, tpr, c="red", label='AUC = {:.3f})'.format(auc))
#
#	plt.xlabel('False Positive Rate')
#	plt.ylabel('True Positive Rate')
#	plt.title('HH Classifier Model ROC Curve')
#	plt.legend(loc='best')
#	plt.grid(True, axis='y', alpha=0.5, linestyle='--')
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

	#model = load_model('final.hdf5')

	pred = model.predict(x_test)

	cvs = TCanvas("c","HH Classifier DNN Score", 900, 600)
	h_pred = TH1D("h_pred", "HH Classifier DNN Score", 40, 0, 1)
        h_pred.GetXaxis().SetTitle("DNN Score")
	h_pred.GetYaxis().SetTitle("Entry")
	h_pred.SetAxisRange(0,1000,"Y")

	for i in range(len(pred)):
		h_pred.Fill(pred[i])

	cvs.cd()
	cvs.SetGrid()
	h_pred.Draw()
	h_sig.Draw("same")
	cvs.SaveAs("dnn.png")
#	cvs.SetLogy()
#	cvs.SaveAs("dnn_log.png")

	print("End of DNN model training...")
	input("Press Enter to continue...")


# User-defined functions
def recall(y_target,y_pred):
	y_target_yn = k.round(k.clip(y_target,0,1))
	y_pred_yn   = k.round(k.clip(y_pred,  0,1))

	count_true_positive = k.sum(y_target_yn*y_pred_yn)
	count_true_positive_false_negative = k.sum(y_target_yn)

	recall = count_true_positive/(count_true_positive_false_negative + k.epsilon())
	return recall

def precision(y_target, y_pred):
	y_pred_yn   = k.round(k.clip(y_pred,  0,1))
	y_target_yn = k.round(k.clip(y_target,0,1))

	count_true_positive = k.sum(y_pred_yn)
	count_true_positive_false_positive = k.sum(y_pred_yn)

	precision = count_true_positive/(count_true_positive_false_positive + k.epsilon())
	return precision

def f1score(y_target, y_pred):
	_recall = recall(y_target, y_pred)
	_precision = precision(y_target, y_pred)
	_f1score = (2*_recall*_precision)/(_recall+_precision+k.epsilon())

	return _f1score


if __name__ == "__main__":
        main()

