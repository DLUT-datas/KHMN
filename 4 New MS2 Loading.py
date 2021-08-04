# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 2021
Part 4: New MS2 Loading
@author: Li Zimu
"""
import os
import pandas as pd
import re



# read feature files
def read_xlsx_file(file_name):
# file_name=Characteristic_ions_file_name
    if '.csv' in file_name:
        df = pd.read_csv(file_name)
# print('running--------')
    else:
        df = pd.read_excel(file_name)
    data=df.values.tolist()
    return data

def check(str):

  my_re = re.compile(r'[A-Za-z]')

  res = re.findall(my_re,str)
  if len(res):

      return False
  else:
      return True


# Read '.mgf' file, transform to list
def seed_MS2toList_align(file_path,mgf_name):
# mgf_name =  seed_file_name
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
                # oneMS2List[0]=int(line[1])
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

            # if line.find("=") == -1:
            if check(line):
                mz_intensity=[]
                line=re.split('[ |\t]+',line)
                mz_intensity.append(float(line[0]))
                mz_intensity.append(float(line[1]))
                spec.append(mz_intensity)
                continue

    return MS2List



# Import all the new '.mgf' files
output_mgf_file_path = 'D:\\Python_code\\script_code\\new_mgf_file\\'

mfg_all_name = os.listdir(output_mgf_file_path)
MS2List = []
for name in mfg_all_name:

    MS2List+=seed_MS2toList_align(output_mgf_file_path,name)
print('Read new generated mgf files Finished')
