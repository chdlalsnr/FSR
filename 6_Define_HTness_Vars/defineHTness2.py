def defineHTness2():

	eta_cand = np.arange(-3.4,3.4,0.4)
        for i in range(len(jet1_mass)):
        	if (lep1_charge[i] == lep2_charge[i]): continue

                chi2_cand,H_cand = [],[]
                for j in range(len(eta_cand)):
                        for k in range(j,len(eta_cand)):

                                lep1_p4 = TLorentzVector()
                                lep2_p4 = TLorentzVector()
                                lep3_p4 = TLorentzVector()

                                jet1_p4 = TLorentzVector()
                                jet2_p4 = TLorentzVector()

                                nu1_p4 = TLorentzVector()
                                nu2_p4 = TLorentzVector()

                                if (lep1_charge < 0): # w1 is w_minus

                                        nu1_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[j],c_met[i][2],c_met[i][0])
                                        nu2_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[k],c_met[i][2],c_met[i][0])

                                        lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
                                        lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])
                                        lep3_p4.SetPtEtaPhiM(c_dilep[i][11],c_dilep[i][12],c_dilep[i][13],c_dilep[i][10])

                                        w1_p4 = lep1_p4 + nu1_p4
                                        w2_p4 = lep2_p4 + nu2_p4
                                        w3_p4 = lep3_p4 + nu1_p4

                                        jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][4])
                                        jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])

                                        t1_p4 = w1_p4 + jet1_p4 # because we don't know charges of jets, I tried to combine w1, w2 with both jet1 and jet2
                                        t2_p4 = w1_p4 + jet2_p4 # top quark
                                        t3_p4 = w2_p4 + jet1_p4 # t-bar
                                        t4_p4 = w2_p4 + jet2_p4 # t-bar

                                        t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                                        w_cand = [w1_p4.M(),w2_p4.M(),w3_p4.M()]
                                        sig_t = np.std(t_cand)
                                        sig_w = np.std(w_cand)

					# 4 tops exist, but 2 terms of 4 terms have same top charges
                                        # there are 2 of W-boson, so when the first top is assigned, its matched W-boson is also assigned
                                        # by the definition, only a t-bar can come to the first term and only the W+-boson can be assigned to the second term
                                        todular1 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                        todular2 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                        todular3 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                        todular4 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)

                                        chi2_cand.append(todular1)
                                        chi2_cand.append(todular2)
                                        chi2_cand.append(todular3)
                                        chi2_cand.append(todular4)

                                        h1_p4 = TLorentzVector()
                                        h2_p4 = TLorentzVector()
                                        h3_p4 = TLorentzVector()

                                        h1_p4 = w1_p4 + w2_p4
                                        h2_p4 = w1_p4 + w3_p4
                                        h3_p4 = w2_p4 + w3_p4

                                        h_cand = [h1_p4.M(),h2_p4.M(),h3_p4.M()]

                                        nunu_p4 = TLorentzVector()
                                        nunu_p4 = nu1_p4 + nu2_p4

                                        sig_h = np.std(h_cand)
                                        sig_nu = 1
                                        #sig_wask = 1
                                        nunu_peak = 30
                                        wask_peak = (1/np.sqrt(3))*np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*h_mass**2*w_mass**2 + w_mass**4) )
                                        # because the condition is lep1_charge is positive, the lep1 goes to only 4th & 5th term
                                        hodular_prev = ((lep2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep1_p4.M()**2-wask_peak**2)**2/sig_w**4) + ((lep1_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep2_p4.M()**2-wask_peak**2)**2/sig_w**4)

                                        hodular1 = ((h1_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev
                                        hodular2 = ((h2_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev
                                        hodular3 = ((h3_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev

                                        H_cand.append(hodular1)
                                        H_cand.append(hodular2)
                                        H_cand.append(hodular3)


				if (lep1_charge > 0): # w1 is w_minus

                                                nu1_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[j],c_met[i][2],c_met[i][0])
                                                nu2_p4.SetPtEtaPhiM(c_met[i][1],eta_cand[k],c_met[i][2],c_met[i][0])

                                                lep1_p4.SetPtEtaPhiM(c_dilep[i][1],c_dilep[i][2],c_dilep[i][3],c_dilep[i][0])
                                                lep2_p4.SetPtEtaPhiM(c_dilep[i][6],c_dilep[i][7],c_dilep[i][8],c_dilep[i][5])
                                                lep3_p4.SetPtEtaPhiM(c_dilep[i][11],c_dilep[i][12],c_dilep[i][13],c_dilep[i][10])

                                                w1_p4 = lep2_p4 + nu2_p4
                                                w2_p4 = lep1_p4 + nu1_p4
                                                w3_p4 = lep3_p4 + nu1_p4

                                                jet1_p4.SetPtEtaPhiM(c_dijet[i][1],c_dijet[i][2],c_dijet[i][3],c_dijet[i][4])
                                                jet2_p4.SetPtEtaPhiM(c_dijet[i][6],c_dijet[i][7],c_dijet[i][8],c_dijet[i][5])

                                                t1_p4 = w1_p4 + jet1_p4 # because we don't know charges of jets, I tried to combine w1, w2 with both jet1 and jet2
                                                t2_p4 = w1_p4 + jet2_p4 # top quark

                                                t3_p4 = w2_p4 + jet1_p4 # t-bar
                                                t4_p4 = w2_p4 + jet2_p4 # t-bar

                                                t_cand = [t1_p4.M(),t2_p4.M(),t3_p4.M(),t4_p4.M()]
                                                w_cand = [w1_p4.M(),w2_p4.M(),w3_p4.M()]
                                                sig_t = np.std(t_cand)
                                                sig_w = np.std(w_cand)

                                                # 4 tops exist, but 2 terms of 4 terms have same top charges
                                                # there are 2 of W-boson, so when the first top is assigned, its matched W-boson is also assigned
                                                # from the definition, only a t-bar can come to the first term and only the W+-boson can be assigned to the second term
                                                todular1 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                                todular2 = ((t3_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                                todular3 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t1_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)
                                                todular4 = ((t4_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((t2_p4.M()**2-t_mass**2)**2/sig_t**4) + ((w1_p4.M()**2-w_mass**2)**2/sig_w**4)

                                                chi2_cand.append(todular1)
                                                chi2_cand.append(todular2)
                                                chi2_cand.append(todular3)
                                                chi2_cand.append(todular4)

                                                h1_p4 = TLorentzVector()
                                                h2_p4 = TLorentzVector()
                                                h3_p4 = TLorentzVector()

						h_cand = [h1_p4.M(),h2_p4.M(),h3_p4.M()]

                                                nunu_p4 = TLorentzVector()
                                                nunu_p4 = nu1_p4 + nu2_p4

                                                sig_h = np.std(h_cand)
                                                sig_nu = 1
                                                #sig_wask = 1
                                                nunu_peak = 30
                                                wask_peak = (1/np.sqrt(3))*np.sqrt( 2*(h_mass**2 + w_mass**2) - np.sqrt(h_mass**4 + 14*h_mass**2*w_mass**2 + w_mass**4) )
                                                # because the condition is lep1_charge is positive, the lep1 goes to only 3rd & 6th term
                                                hodular_prev = ((lep1_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep2_p4.M()**2-wask_peak**2)**2/sig_w**4) + ((lep2_p4.M()**2-w_mass**2)**2/sig_w**4) + ((lep1_p4.M()**2-wask_peak**2)**2/sig_w**4)

                                                hodular1 = ((h1_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev
                                                hodular2 = ((h2_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev
                                                hodular3 = ((h3_p4.M()**2-h_mass**2)**2/sig_h**4) + ((nunu_p4.M()**2-nunu_peak**2)**2/sig_nu**4) + hodular_prev

                                                H_cand.append(hodular1)
                                                H_cand.append(hodular2)
                                                H_cand.append(hodular3)

		topness.append(np.log(min(chi2_cand)))
                higgsness.append(np.log(min(H_cand)))

	return topness,higgsness
