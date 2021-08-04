# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 2021
Part 2: Feature Table Input and Alignment
The third column m/z,
The fourth column rt,
Feature align：m/z, 0.005Da; rt, 6s
@author: Li Zimu
"""
import os
import pandas as pd



# Read feature files
def read_xlsx_file(file_name):
    if '.csv' in file_name:
        df = pd.read_csv(file_name)
#        print('running--------')
    else:
        df = pd.read_excel(file_name)
    data=df.values.tolist()
    return data


def output_file(filename,output_data):
#    filename = 'output_list_network.xlsx'
    df = pd.DataFrame(output_data, columns=['Align_ID','ID', 'Index','Peak_Name','m/z','RT','Group','Use','Tissue'])
    df.to_excel(filename, index=False)


filePath = 'D:\\Python_code\\script_code\\feature_file\\'

peak_all_name = os.listdir(filePath)
peak_data = []
for peak_name in peak_all_name:
    peak_data += read_xlsx_file(filePath+peak_name)

peak_data=sorted(peak_data,key=lambda x:x[3])
mz_threshold = 0.005
rt_threshold = 6
align_id = [-1]*len(peak_data)
temp_id = 1

for i in range(len(peak_data)):
    if align_id[i] == -1:
        align_id[i] = temp_id
        temp_id+=1

    peak_1 = peak_data[i]
    rt1 = peak_1[4]
    mz1 = peak_1[3]
    for j in range(len(peak_data)):
        if i>=j:
            continue
        else:
            peak_2 = peak_data[j]
            rt2 = peak_2[4]
            mz2 = peak_2[3]
            if abs(mz1-mz2)<=mz_threshold and abs(rt1-rt2)*60<=rt_threshold:
                if align_id[j] == -1:
                    align_id[j] = align_id[i]
                    print(align_id[i])


output_data = []
for i in range(len(peak_data)):
#    print(['Align_'+str(align_id[i])]+peak_data[i])
    output_data.append(['Align_'+str(align_id[i])]+peak_data[i])
# output all information
output_filename = 'D:\\Python_code\\script_code\\POS_tissue_all_features_Align.xlsx'
output_file(output_filename,output_data)

output_data_silk = []
output_data_leave = []
output_data_seed = []
for i in range(len(peak_data)):
#    print(['Align_'+str(align_id[i])]+peak_data[i])
    if peak_data[i][7] == 'silk':
        output_data_silk.append(['Align_'+str(align_id[i])]+peak_data[i])
    if peak_data[i][7] == 'Leave':
        output_data_leave.append(['Align_'+str(align_id[i])]+peak_data[i])
    if peak_data[i][7] == 'Seed':
        output_data_seed.append(['Align_'+str(align_id[i])]+peak_data[i])

output_silk_filename = 'D:\\Python_code\\script_code\\POS_silk_all_features_Align.xlsx'
output_file(output_silk_filename,output_data_silk)

output_Leave_filename = 'D:\\Python_code\\script_code\\POS_Leave_all_features_Align.xlsx'
output_file(output_Leave_filename,output_data_leave)

output_Seed_filename = 'D:\\Python_code\\script_code\\POS_Seed_all_features_Align.xlsx'
output_file(output_Seed_filename,output_data_seed)

feature_file = "D:\\Python_code\\script_code\\POS_tissue_all_features_Align.xlsx"
feature_file_data = read_xlsx_file(feature_file)
feature_mz=[]
feature_rt=[]
feature_id= []
feature_align= []
for line in feature_file_data:
     feature_id.append(line[1])
     feature_align.append(line[0])
     feature_mz.append(float(line[4]))
     feature_rt.append(float(line[5])*60)
print('Read peak table Finished')
