# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:06:15 2021
Part 3: MS2 Input and Filtering
3.1: Simple Conversion
@author: Li Zimu
"""
import os
import re


def check(str):

  my_re = re.compile(r'[A-Za-z]')

  res = re.findall(my_re,str)
  if len(res):

      return False
  else:
      return True
# Read 'mgf' file, transform to list
def seed_MS2toList(file_path,mgf_name):
#    file_path = mgf_file_path_seed
#    mgf_name =  'CGFs_Database_Positive.mgf'
    MS2List=[]
    seed_name = ["HCAAs" ,"CGFs","BOAs"]
    if any(e in mgf_name for e in seed_name):
        mark = "seed_metabolite"
    else:
        mark = 0
    with open(file_path+mgf_name)as f:
        for line in f:
            line=line.rstrip()
            if len(line) < 4:
                continue
#            print(line)
            if line== "BEGIN IONS":
                oneMS2List=[0,0,0,0,0]
                spec=[]
                continue

            if line[:5]=='NAME=':
                oneMS2List[2]=line[5:]
                continue
#                oneMS2List[0]=int(line[1])
            # seed file RT name: Predicted_RT
            if line[:13] == "Predicted_RT=":
                oneMS2List[3]=float(line[13:])*60
                continue
            # experimental data RT name: RTINSECONDS
            if line[:12] == "RTINSECONDS=":
                oneMS2List[3]=float(line[12:])
                continue

            if line[:8]=='PEPMASS=':
                line=re.split('=| ',line)
                oneMS2List[4]=float(line[1])
                continue

            if line=="END IONS":
                oneMS2List[0]=mgf_name
                oneMS2List[1]=mark
                oneMS2List.append(spec)
                MS2List.append(oneMS2List)
                continue

#            if line.find("=") == -1:
            if check(line):
                mz_intensity=[]
                line=re.split('[ |\t]+',line)
                mz_intensity.append(float(line[0]))
                mz_intensity.append(float(line[1]))
                spec.append(mz_intensity)
                continue

    return MS2List

# Read '.mgf' file and transform to list
def seed_MS2toList_align(file_path,mgf_name):
#    mgf_name =  seed_file_name
    MS2List=[]
    seed_name = ["HCAAs" ,"CGFs","BOAs"]
    if any(e in mgf_name for e in seed_name):
        mark = "seed_metabolite"
    else:
        mark = 0
    with open(file_path+mgf_name)as f:
        for line in f:
            line=line.rstrip()
            if len(line) < 4:
                continue
#            print(line)
            if line== "BEGIN IONS":
                oneMS2List=[0,0,0,0,0]
                id_list_1 = ''
                id_list_2 = ''
                spec=[]
                continue
            if line[:3]=='ID=':
                id_list_1=line[3:]
                continue
            if line[:9]=='Align_ID=':
                id_list_2=line[9:]
                continue
            if line[:5]=='NAME=':
                oneMS2List[2]=line[5:]
                continue
#                oneMS2List[0]=int(line[1])
            if line[:12] == "RTINSECONDS=":
                oneMS2List[3]=float(line[12:])
                continue
            if line[:8]=='PEPMASS=':
                line=re.split('=| ',line)
                oneMS2List[4]=float(line[1])
                continue

            if line=="END IONS":
                oneMS2List[0]=mgf_name
                oneMS2List[1]=mark
                oneMS2List.append(spec)
                oneMS2List.append(id_list_1)
                oneMS2List.append(id_list_2)
                MS2List.append(oneMS2List)
                continue

#            if line.find("=") == -1:
            if check(line):
                mz_intensity=[]
                line=re.split('[ |\t]+',line)
                mz_intensity.append(float(line[0]))
                mz_intensity.append(float(line[1]))
                spec.append(mz_intensity)
                continue

    return MS2List

# Find the feature ID corresponding to the MS/MS
def seekMS2id_ppm_tissue(feature_id,weight,rt,MS2weight,MS2rt,MS2_mgf,delata_weight,rt_threshold):

    candidate_MS2id=[]
    MS2List=[]
    for i in range(len(weight)):
        candidate_MS2id.append([])
    
#   delata_weight=15                    #mass tolerance: ppm
    last_index=0
    
    for index_weight in range(len(weight)):
        knownmz=weight[index_weight]
        knownrt=rt[index_weight]
        knownid = feature_id[index_weight]
        Time=0
        index=last_index
        if (MS2weight[index]-knownmz)/knownmz*1000000>delata_weight:
            continue
        
        for index_MS2weight in range(index,len(MS2weight)+1):
            if knownid[0:4].lower() not in MS2_mgf[index].lower():
                break
            if index_MS2weight==len(MS2weight):
                break
            if (MS2weight[index_MS2weight]-knownmz)/knownmz*1000000>delata_weight: # MS interval
                break
            if abs(knownmz-MS2weight[index_MS2weight])/knownmz*1000000<=delata_weight:
                Time+=1
                if abs(knownrt-MS2rt[index_MS2weight])<=rt_threshold:             # RT interval
                    
                    candidate_MS2id[index_weight].append(index_MS2weight)
                if Time==1:
                    last_index=index_MS2weight
                    
        if last_index==len(MS2weight) and (knownmz-MS2weight[last_index])/knownmz*1000000>delata_weight:
            break
    
    for i in range(len(candidate_MS2id)):
        
        error=[]
        if len(candidate_MS2id[i]):
            for j in candidate_MS2id[i]:
                error.append(abs(weight[i]-MS2weight[j])+abs(rt[i]-MS2rt[j]))
            MS2List.append(candidate_MS2id[i][error.index(min(error))])
        else:
            MS2List.append(-1)
    
    return MS2List

# Output '.mgf' file
def get_mgf_string(spectrum):
    output_lines = []
    output_lines.append("BEGIN IONS")
#    output_lines.append("BEGIN IONS")

    output_lines.append("RTINSECONDS=" + str(spectrum[3]))
    output_lines.append("PEPMASS=" + str(spectrum[4]))
    output_lines.append(get_mgf_peak_string(spectrum[5]))
    output_lines.append("END IONS")

    return "\n".join(output_lines)

def get_mgf_peak_string(peaks):
    output_string = ""
    for peak in peaks:
        output_string += str(peak[0]) + "\t" + str(peak[1]) + "\n"

    return output_string
# Output to an new MGF
def save_to_mgf(MS2List, output_mgf):
    output_file = open(output_mgf, "w")
    for spectrum in MS2List:
        if spectrum != None:
            output_file.write(str(get_mgf_string(spectrum))+"\n")
    output_file.close()
def get_mgf_string_align(spectrum):
    output_lines = []
    output_lines.append("BEGIN IONS")
#    output_lines.append("BEGIN IONS")
    output_lines.append("ID=" + str(spectrum[6]))
    output_lines.append("Align_ID=" + str(spectrum[7]))
    output_lines.append("RTINSECONDS=" + str(spectrum[3]))
    output_lines.append("PEPMASS=" + str(spectrum[4]))
    output_lines.append(get_mgf_peak_string(spectrum[5]))
    output_lines.append("END IONS")

    return "\n".join(output_lines)
# output to an mgf
def save_to_mgf_align(MS2List, output_mgf):
    output_file = open(output_mgf, "w")
    for spectrum in MS2List:
        if spectrum != None:
            output_file.write(str(get_mgf_string_align(spectrum))+"\n")
    output_file.close()

def del_value(list_i,value):
   j = []
   j_loc = []
   for i in range(len(list_i)):
       if list_i[i] == value:
           continue
       else:
           j.append(list_i[i])
           j_loc.append(i)
   return j,j_loc



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Import all '.mgf' files, Feature matching and filter, output new '.mgf' file

filePath = 'D:\\Python_code\\script_code\\Condition-A_mgf_file\\'
output_mgf_file_path = 'D:\\Python_code\\script_code\\new_Condition-A_mgf_file\\' 
mfg_all_name = os.listdir(filePath)

mz_threshold = 15
rt_threshold = 18
MS2List = []
for name in mfg_all_name:
#    name = mfg_all_name[0]
    temp_MS2 = seed_MS2toList(filePath,name)
    for i in temp_MS2:
        MS2List.append(i)
MS2List=sorted(MS2List,key=lambda x:x[4])

MS2List_rt=[]
MS2List_mz=[]
MS2List_mgf=[]
for one_MS2List in MS2List:
    MS2List_rt.append(one_MS2List[3])
    MS2List_mz.append(one_MS2List[4])
    MS2List_mgf.append(one_MS2List[0])
unknown_MS2id = []
# MS2id

unknown_MS2id=seekMS2id_ppm_tissue(feature_id,feature_mz,feature_rt,MS2List_mz,MS2List_rt,MS2List_mgf,mz_threshold,rt_threshold)

del_unknown_MS2id, del_unknown_MS2id_loc= del_value(unknown_MS2id,-1)

output_mgf_file_name = output_mgf_file_path+'new_' + name
if len(del_unknown_MS2id) == 0:
    print('None')
else:
    MS2=[]
    for i in range(len(del_unknown_MS2id)):
        temp_MS2List = []
        id_1 = del_unknown_MS2id[i]
        id_2 = del_unknown_MS2id_loc[i]
        temp_MS2List = MS2List[id_1].copy()
        temp_MS2List.append(feature_id[id_2])
        temp_MS2List.append(feature_align[id_2])
        MS2.append(temp_MS2List)

    save_to_mgf_align(MS2, output_mgf_file_name)
print('Feature matching Finished')
print('Generate new .mgf files Finish')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
