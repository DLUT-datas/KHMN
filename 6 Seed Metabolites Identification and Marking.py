# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 2021
Part 6: Seed Metabolites Identification and Marking
@author: Li Zimu
"""

import Cosine

# Peak filter
def window_filter_peaks(peaks, window_size, top_peaks):
    peak_masses_to_keep = []
    for i in range(len(peaks)):
        peak_window=[]
        mass=peaks[i][0]
        for j in range(len(peaks)):
            candidate_mass=peaks[j][0]
            if candidate_mass-mass<-window_size:
                continue
            if abs(candidate_mass-mass)<=window_size:
                peak_window.append(peaks[j])
            if candidate_mass-mass>50:
                break
        peak_window=sorted(peak_window, key=lambda peak: peak[1],reverse=True)
        for peak in peak_window[:top_peaks]:
            if mass==peak[0]:
                peak_masses_to_keep.append(peaks[i])
                break

    return peak_masses_to_keep

def filter_precursor_peaks(peaks, tolerance_to_precursor, mz):
    new_peaks = []
    for peak in peaks:
        if abs(peak[0] - mz) > tolerance_to_precursor:
            new_peaks.append(peak)
    return new_peaks

def MS2List2Spec(all_MS2List,min_peak_no):
    spec2=[]
    ID = []
    final_node = []
    for i in range(len(all_MS2List)):
        ID.append(i)
        tR=all_MS2List[i][3]
        pm=all_MS2List[i][4]
        spec=all_MS2List[i][5]
        if len(all_MS2List[i])>6:

            final_node.append([i,pm,tR,all_MS2List[i][0],all_MS2List[i][1],all_MS2List[i][2],all_MS2List[i][6],all_MS2List[i][7]])
        else:
            final_node.append([i,pm,tR,all_MS2List[i][0],all_MS2List[i][1],all_MS2List[i][2],'',''])

        spec=filter_precursor_peaks(spec, 17, pm)
        spec=window_filter_peaks(spec, 50, min_peak_no)

        spec2.append(spec)

    return spec2,final_node,ID

def identify_seed(MS2List_seed,MS2List,ms1_tolerance,ms2_tolerance):

    spec1,final_node,ID = MS2List2Spec(MS2List_seed,6)
    spec2,final_node,ID = MS2List2Spec(MS2List,6)

    for i in range(len(MS2List_seed)):
         pm1=MS2List_seed[i][4]
         rt1=MS2List_seed[i][3]

         for j in range(len(MS2List)):
            pm2=MS2List[j][4]
            rt2=MS2List[j][3]
            if abs(pm1 - pm2)/pm1*1000000<=15 and abs(rt1-rt2)/rt1*100<=20:
                MS2_similarity,reported_alignments=Cosine.score_alignment(spec1[i],spec2[j],pm1,pm2,ms1_tolerance,ms2_tolerance,max_charge_consideration=1)
                if len(reported_alignments)>=6 and MS2_similarity>=0.6:
                    MS2List[j][1] = 'seed_metabolite_from_experiment'
                    MS2List[j][2] = MS2List_seed[i][2]
                    print('seed_metabolite_from_experiment')
    return MS2List




# Mark seed netabolite
ms1_tolerance_threshold = 0.01
ms2_tolerance_threshold = 0.02
MS2List = identify_seed(MS2List_seed,MS2List,ms1_tolerance_threshold,ms2_tolerance_threshold)

print('Mark seed form experiment data Finished')
