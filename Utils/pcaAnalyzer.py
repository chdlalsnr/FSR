# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression, SelectKBest


def main():

	target = 'target'
	feature = ['jet1_mass','jet1_pt','jet1_eta','jet1_phi','jet2_mass','jet2_pt','jet2_eta','jet2_phi','lep1_mass','lep1_pt','lep1_eta','lep1_phi','lep2_mass','lep2_pt','lep2_eta','lep2_phi','dr_jj','dr_ll','met_pt','met_phi','h1_mass','h1_pt','h1_eta','h1_phi','h2_mass','h2_pt','h2_eta','h2_phi','higgsness','topness']

	df_hh = pd.read_csv("results/htness_sig.csv")
	df_tt = pd.read_csv("results/htness_bkg.csv")

	data = pd.concat([df_hh, df_tt])

	train,valid = train_test_split(data, test_size=0.3, random_state=42)

	X_train = train[feature]
	y_train = train[target]
	X_valid = valid[feature]
	y_valid = valid[target]

	selector = SelectKBest(score_func=f_regression, k=15)

	X_train_selected = selector.fit_transform(X_train, y_train)

	X_valid_selected = selector.transform(X_valid)
	X_train_selected.shape, X_valid_selected.shape

	all_names = X_train.columns
	selected_mask = selector.get_support()
	selected_names = all_names[selected_mask]

	print(selected_names)


#	df = data[feature]
#	df_scaled = scaler.fit_transform(df)
#
#	pca = PCA(n_components=30)
#
#	pca_arr = pca.fit_transform(df_scaled)
#	pca_df = pd.DataFrame(pca_arr, index=data.index,
#				columns=[f"PC{num+1}" for num in range(df.shape[1])])
#
#	pc_df = pd.DataFrame({'explained variances(eigenvalues)':pca.explained_variance_,
#				'explained_variance_ratio':pca.explained_variance_ratio_},
#				index=pca_df.columns.tolist())
#	pc_df['cumulative contribution'] = pc_df['explained_variance_ratio'].cumsum()
#
#	fig,ax1 = plt.subplots(figsize=(10,6))
#
#	ax1.set_title("Scree plot")
#	ax1.plot(pca_df.columns.tolist(), pc_df['explained variances(eigenvalues)'].tolist(), color='red',alpha=0.5)
#	ax1.plot(pca_df.columns.tolist(), pc_df['explained variances(eigenvalues)'].tolist(), color='red',marker='o',markersize=5)
#	ax1.set_xlabel("Principal component number")
#	ax1.set_ylabel("Explained variance number (eigenvalue)",color='red')
#	ax1.axvline(x='PC5', color='red', linestyle='dashed',alpha=0.7)
#	ax1.tick_params(axis='y', labelcolor='red')
#
#	ax1.set_xticklabels(pca_df.columns.tolist(),rotation=45)
#
#	ax2 = ax1.twinx()
#	ax2.plot(pca_df.columns.tolist(), pc_df['cumulative contribution'].tolist(), color='gray',alpha=0.5)
#	ax2.plot(pca_df.columns.tolist(), pc_df['cumulative contribution'].tolist(), color='gray',marker='^',markersize=5)
#	ax2.set_ylabel("Cumulative explained variance ratio", color='gray')
#	ax2.axvline(x='PC15', color='gray', linestyle='dashed',alpha=0.7)
#
#	plt.grid(True,axis='y')
#	plt.show()

	input("Press Enter to continue...")


if __name__ == "__main__":
        main() 
