# -*- coding: utf-8 -*-
from glob import glob
import numpy as np
import csv
import sys
from scipy.optimize import minimize


def main():


	t1_hh,t2_hh,t3_hh,t4_hh,h1_hh,h2_hh,h3_hh = [],[],[],[],[],[],[]
	t1_tt,t2_tt,t3_tt,t4_tt,h1_tt,h2_tt,h3_tt = [],[],[],[],[],[],[]

	hh = open('htness_hh.csv','r')
	hhline = csv.reader(hh)

	tt = open('htness_tt.csv','r')
	ttline = csv.reader(tt)

	j = 0
	for line in hhline:
		if 'jet1_mass' in line: continue
		if j > 5000: break
		h3_hh.append(float(line[-2]))
		h2_hh.append(float(line[-3]))
		h1_hh.append(float(line[-4]))
		t4_hh.append(float(line[-5]))
		t3_hh.append(float(line[-6]))
		t2_hh.append(float(line[-7]))
		t1_hh.append(float(line[-8]))
		j += 1

	j = 0
	for line in ttline:
		if 'jet1_mass' in line: continue
		if j > 5000: break
		h3_tt.append(float(line[-2]))
                h2_tt.append(float(line[-3]))
                h1_tt.append(float(line[-4]))
                t4_tt.append(float(line[-5]))
                t3_tt.append(float(line[-6]))
                t2_tt.append(float(line[-7]))
                t1_tt.append(float(line[-8]))
		j += 1

	ww = []
	for i in range(j):

		def optimizeHTTerm(x):

			w1,w2,w3,w4,w5,w6,w7 = x[0],x[1],x[2],x[3],x[4],x[5],x[6]

			hh_t1 = t1_hh[i]*w1
			hh_t2 = t2_hh[i]*w2
			hh_t3 = t3_hh[i]*w3
			hh_t4 = t4_hh[i]*w4
			hh_h1 = h1_hh[i]*w5
			hh_h2 = h2_hh[i]*w6
			hh_h3 = h3_hh[i]*w7

			tt_t1 = t1_tt[i]*w1
			tt_t2 = t2_tt[i]*w2
			tt_t3 = t3_tt[i]*w3
			tt_t4 = t4_tt[i]*w4
			tt_h1 = h1_tt[i]*w5
			tt_h2 = h2_tt[i]*w6
			tt_h3 = h3_tt[i]*w7

			res = -np.sqrt( ((hh_t1 + hh_t2 + hh_t3 + hh_t4) - (tt_t1 + tt_t2 + tt_t3 + tt_t4))**2 + ((hh_h1 + hh_h2 + hh_h3) - (tt_h1 + tt_h2 + tt_h3))**2 )
			# np.sqrt((hh_topness - tt_topness)**2 + (hh_higgsness - tt_higgsness)**2)

			return res

		ini = np.array([0,0,0,0, 0,0,0])
		bnd = [(0,1),(0,1),(0,1),(0,1), (0,1),(0,1),(0,1)]

		cons = ({'type':'eq',
			 'fun': lambda x: x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6] - 1})

		opt = minimize(optimizeHTTerm, ini, bounds=bnd, method='SLSQP', constraints=cons, options={'maxiter':1000})

		ww.append([opt.x[0],opt.x[1],opt.x[2],opt.x[3], opt.x[4],opt.x[5],opt.x[6]])


	w1,w2,w3,w4,w5,w6,w7 = [],[],[],[],[],[],[]
	for i in range(len(ww)):

		w1.append(ww[i][0])
		w2.append(ww[i][1])
		w3.append(ww[i][2])
		w4.append(ww[i][3])
		w5.append(ww[i][4])
		w6.append(ww[i][5])
		w7.append(ww[i][6])

	weight = [np.mean(w1),np.mean(w2),np.mean(w3),np.mean(w4),np.mean(w5),np.mean(w6),np.mean(w7)]

	print(weight)

	input("Press Enter to continue...")



if __name__ == "__main__":
        main()


