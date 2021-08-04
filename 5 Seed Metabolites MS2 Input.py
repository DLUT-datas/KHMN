# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 2021
Part 5: Seed Metabolites MS2 Input
@author: Li Zimu
"""
import os
import pandas as pd
import re


# Read feature files
def read_xlsx_file(file_name):

    if '.csv' in file_name:
        df = pd.read_csv(file_name)

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
def seed_MS2toList(file_path,mgf_name):
# file_path = mgf_file_path_seed
# mgf_name =  'CGFs_Database_Positive.mgf'
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
                spec=[]
                continue

            if line[:5]=='NAME=':
                oneMS2List[2]=line[5:]
                continue
                # oneMS2List[0]=int(line[1])
                # Predicted_RT: RT in Seed file
            if line[:13] == "Predicted_RT=":
                oneMS2List[3]=float(line[13:])*60
                continue
            # RTINSECONDS: RT in experimental data
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

            # if line.find("=") == -1:
            if check(line):
                mz_intensity=[]
                line=re.split('[ |\t]+',line)
                mz_intensity.append(float(line[0]))
                mz_intensity.append(float(line[1]))
                spec.append(mz_intensity)
                continue

    return MS2List






# import seeds '.mgf' file
mgf_file_path_seed = 'D:\\Python_code\\script_code\\seed_mgf_file\\'
mfg_all_name_seed = os.listdir(mgf_file_path_seed)
MS2List_seed = []
for name in mfg_all_name_seed:
    MS2List_seed+=seed_MS2toList(mgf_file_path_seed,name)
print('Raed seeds mgf files Finished')
