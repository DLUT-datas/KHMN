# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 2021
Part 9: Nodes and Edges Annotation
9.3 Edges Annotation of BOAs
@author: Li Zimu
"""
import os
import pandas as pd
import re


def list_2_str(list_name):
    str = '//'
    if len(list_name):
        return str.join(list_name)
    else:
        return ''

def output_edge_annotation_result_align(filename,final_node):
    # filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_node, columns=['node1_ID', 'node2_ID','diff_ms','cosine_score','group','ID_1','ID_2','Align_1','Align_2', 'node2_precursor_mass','node2_RT','node2_mgf','node2_mark', 'node2_precursor_mass','node2_RT','node2_mgf','node2_mark','Enzyme_catalyzed_reactions'])
    df.to_excel(filename, index=False)

# Annotation of edges
def edge_annotation_align(output_filename,edge_data,Enzyme_data,mz_threshold=0.005):
    # output_filename = 'output_edge_annotation_result.xlsx'
    num_edge = len(edge_data)
    num_enzyme = len(Enzyme_data)

    for i in range(num_edge):
        mass_diff = edge_data[i][2]
        # print(mass_diff)
        temp_result=[]
        for j in range(num_enzyme):
            mass_diff_Enzyme = Enzyme_data[j][4]
            # print(mass_diff_Enzyme)
            if abs(abs(mass_diff)-mass_diff_Enzyme)<=mz_threshold:
                temp_result.append(Enzyme_data[j][2])

        edge_data[i].append(list_2_str(temp_result))
    output_edge_annotation_result_align(output_filename,edge_data)




# Read annotation files
def read_xlsx_file(file_name):
#    file_name=Characteristic_ions_file_name
    if '.csv' in file_name:
        df = pd.read_csv(file_name)
#        print('running--------')
    else:
        df = pd.read_excel(file_name)
    data=df.values.tolist()
    return data



# Read edges file
# Enzyme-catalyzed reactions list.xlsx
# Enzyme_catalyzed reactions list.csv
Enzyme_file_name = r'D:\Python_code\script_code\annotation_file\Enzyme_catalyzed_reactions_list.csv'
Enzyme_data = read_xlsx_file(Enzyme_file_name)
# Read edges of network
edge_file_name = 'output_seed_filter_network_info_align.xlsx'
edge_data = read_xlsx_file(edge_file_name)
# Output edges annotation result
# output_edge_annotation_result_filename = 'output_network_edge_annotation_result.xlsx'
# edge_data = edge_annotation(output_edge_annotation_result_filename,edge_data,Enzyme_data,mz_threshold=0.005)
output_edge_annotation_result_filename = 'output_network_edge_annotation_result_align.xlsx'
edge_data = edge_annotation_align(output_edge_annotation_result_filename,edge_data,Enzyme_data,mz_threshold=0.005)
print('Network edges annotation Finished')
