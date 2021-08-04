# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 16:09:37 2021
Part 10: Redundancy Removal
@author: Li Zimu
"""

import pandas as pd



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

def del_redundancy_info(annotation_data):
    output_result = []
    for i in range(len(annotation_data)):

        one_data = annotation_data[i].copy()
        if 'Database' in one_data[11]:
            one_data.append(one_data[7])
        elif '15V' in one_data[11]:
            one_data.append(one_data[7]+'_B')

        elif '30V' in one_data[11]:
            one_data.append(one_data[7]+'_C')
        elif '45V' in one_data[11]:
            one_data.append(one_data[7]+'_D')
        else:
            one_data.append(one_data[7]+'_A')

        if 'Database' in one_data[15]:
            one_data.append(one_data[8])
        elif '15V' in one_data[15]:
            one_data.append(one_data[8]+'_B')

        elif '30V' in one_data[15]:
            one_data.append(one_data[8]+'_C')
        elif '45V' in one_data[15]:
            one_data.append(one_data[8]+'_D')
        else:
            one_data.append(one_data[8]+'_A')

        output_result.append(one_data)

    return output_result

def output_HCCAs_annotation_del(final_edge_data,filename):

    df = pd.DataFrame(final_edge_data, columns=['node_ID1', 'node_ID2','diff_ms','cosine_score','group', 'Tissue_ID1','Tissue_ID2','Align_ID1','Align_ID2',\
                                                      'node1_precursor_mass','node1_RT','node1_mgf','node1_mark', \
                                                      'node2_precursor_mass','node2_RT','node2_mgf','node2_mark','Enzyme_catalyzed_reactions',\
                                                      'node1_Characteristic_ion','node1_Neutral_loss_of_Amine','node1_Neutral_loss_of_HCA',\
                                                      'node2_Characteristic_ion','node2_Neutral_loss_of_Amine','node2_Neutral_loss_of_HCA',\
                                                      'Align_ID_Platform1','Align_ID_Platform2'])
    df.to_excel(filename, index=False)

def output_CGFs_annotation_del(final_edge_data,filename):

    df = pd.DataFrame(final_edge_data, columns=['node_ID1', 'node_ID2','diff_ms','cosine_score','group',  'Tissue_ID1','Tissue_ID2','Align_ID1','Align_ID2',\
                                                      'node1_precursor_mass','node1_RT','node1_mgf','node1_mark', \
                                                      'node2_precursor_mass','node2_RT','node2_mgf','node2_mark','Enzyme_catalyzed_reactions',\
                                                      'node1_Sub_motif','node1_Glycosyl_neutral_loss',\
                                                      'node2_Sub_motif','node2_Glycosyl_neutral_loss',\
                                                      'Align_ID_Platform1','Align_ID_Platform2'])
    df.to_excel(filename, index=False)

def output_BOAs_annotation_del(final_edge_data,filename):

    df = pd.DataFrame(final_edge_data, columns=['node_ID1', 'node_ID2','diff_ms','cosine_score','group',  'Tissue_ID1','Tissue_ID2','Align_ID1','Align_ID2',\
                                                      'node1_precursor_mass','node1_RT','node1_mgf','node1_mark', \
                                                      'node2_precursor_mass','node2_RT','node2_mgf','node2_mark','Enzyme_catalyzed_reactions',\
                                                      'Align_ID_Platform1','Align_ID_Platform2'])
    df.to_excel(filename, index=False)


# BOAs only edges annotated
#BOAs_annotation_file = 'output_network_edge_annotation_result_align.xlsx'
#BOAs_annotation_data = read_xlsx_file(BOAs_annotation_file)
#BOAs_output_result = del_redundancy_info(BOAs_annotation_data)
#output_BOAs_annotation_del(BOAs_output_result,'final_output_BOAs_annotation_result_del.xlsx')


#HCCAs_annotation_file = 'final_output_HCCAs_annotation_result.xlsx'
#HCAAs_annotation_data = read_xlsx_file(HCCAs_annotation_file)
#HCAAs_output_result = del_redundancy_info(HCAAs_annotation_data)
#output_HCCAs_annotation_del(HCAAs_output_result,'final_output_HCCAs_annotation_result_del.xlsx')


CGF_annotation_file = 'final_output_CGFs_annotation_result.xlsx'
CGFs_annotation_data = read_xlsx_file(CGF_annotation_file)
CGFs_output_result = del_redundancy_info(CGFs_annotation_data)
output_CGFs_annotation_del(CGFs_output_result,'final_output_CGFs_annotation_result_del.xlsx')
