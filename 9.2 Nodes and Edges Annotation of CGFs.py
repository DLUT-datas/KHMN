# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:06:15 2021
Part 9: Nodes and Edges Annotation
9.2 Nodes and Edges Annotation of CGFs
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
#    output_filename = 'output_edge_annotation_result.xlsx'
    num_edge = len(edge_data)
    num_enzyme = len(Enzyme_data)

    for i in range(num_edge):
        mass_diff = edge_data[i][2]
#        print(mass_diff)
        temp_result=[]
        for j in range(num_enzyme):
            mass_diff_Enzyme = Enzyme_data[j][4]
#            print(mass_diff_Enzyme)
            if abs(abs(mass_diff)-mass_diff_Enzyme)<=mz_threshold:
                temp_result.append(Enzyme_data[j][2])

        edge_data[i].append(list_2_str(temp_result))
    output_edge_annotation_result_align(output_filename,edge_data)

# Normalize intensity for mgf files
def max_normalize_spectrum(oneMS2List):

    intermediate_output_spectrum = []
    intermediate_output_spectrum_mz = []
    one_spectrum = []
    sum_ms = 0.0
    for msms in oneMS2List[5]:
        sum_ms+= msms[1]
        intermediate_output_spectrum.append(msms[1])
        intermediate_output_spectrum_mz.append(msms[0])
    # for s in intermediate_output_spectrum:
        # print(s)
    max_ms = max(intermediate_output_spectrum)
    new_spectrum = [i/max_ms for i in intermediate_output_spectrum]
    one_spectrum.append(oneMS2List)
    one_spectrum.append(intermediate_output_spectrum_mz)
    one_spectrum.append(new_spectrum)

    return one_spectrum

def all_spectrum_normalize(MS2List):

    output_spectrum = []
    for oneMS2List in MS2List:
        # print(max_normalize_spectrum(oneMS2List))
        output_spectrum.append(max_normalize_spectrum(oneMS2List))
    return output_spectrum



def output_node_ions_annotation(filename,norm_MS2List,fianl_ion_result,fianl_Amine_result,fianl_HCA_result):
    # filename = 'output_list_network.xlsx'
    final_node = []
    for i in range(len(norm_MS2List)):
        one_MS2List = norm_MS2List[i][0]
        tR=one_MS2List[3]
        pm=one_MS2List[4]
        # print([i,pm,tR,one_MS2List[0],one_MS2List[1],one_MS2List[2],motifs_result_list_final[i],gly_result_list_final[i]])
        final_node.append([i,pm,tR,one_MS2List[0],one_MS2List[1],one_MS2List[2],fianl_ion_result[i],fianl_Amine_result[i],fianl_HCA_result[i]])

    df = pd.DataFrame(final_node, columns=['node_ID', 'node_precursor_mass','node_RT','node_mgf','mark','metabolite_name','Characteristic_ion','Neutral_loss_of_Amine','Neutral_loss_of_HCA'])
    df.to_excel(filename, index=False)

def match_Glycosyl(Glycosyl_data,diff_pepmass_M_H,mz_threshold):
    gly_result_list = []
    gly_result = ''
    for gly in Glycosyl_data:
        if abs(diff_pepmass_M_H-gly[3])<=mz_threshold:
            gly_result_list.append(gly[1])
            if gly_result == '':
                gly_result = gly[1]
            else:
                if gly[1] not in gly_result:
                    gly_result =  gly_result+'//'+ gly[1]
    return gly_result_list,gly_result
def match_one_motifs(Motifs_data,one_list,Glycosyl_data):

    motifs_name = []
    for motifs in Motifs_data:
        motifs_name.append(motifs[1])

    motifs_name = list(set(motifs_name))
    motifs_result=''
    motifs_result_list = []
    gly_result_list = []
    gly_result = ''
    for kk in range(len(motifs_name)):
        one_motifs_mz = []
        for one_motifs in Motifs_data:
            if one_motifs[1] ==motifs_name[kk]:
                one_motifs_mz.append(one_motifs[5])

                if one_motifs[3] =='[M+H]+':
                    M_H = one_motifs[5]
                    # print(M_H)
        # print(one_motifs_mz)
        # Determine whether it contains motifs
        if len(one_list[1])>=len(one_motifs_mz):

            motifs_loc = []
            for i in range(len(one_list[1])):
                for j in range(len(one_motifs_mz)):
                    if abs(one_list[1][i]-one_motifs_mz[j])<=mz_threshold:
                        motifs_loc.append(j)
            motifs_loc = list(set(motifs_loc))

            if len(motifs_loc) == len(one_motifs_mz):

                motifs_result_list.append(motifs_name[kk])
                if motifs_result == '':
                    motifs_result = motifs_name[kk]
                else:
                    if motifs_name[kk] not in motifs_result:
                        motifs_result =  motifs_result+'//'+ motifs_name[kk]

                PEPMASS = one_list[0][4]
                diff_pepmass_M_H = abs(PEPMASS-M_H)
                gly_list,gly = match_Glycosyl(Glycosyl_data,diff_pepmass_M_H,mz_threshold)
                if gly =='':
                    continue
                else:

                    gly_result_list.append(gly_list)

                    if gly_result == '':
                        gly_result = gly
                    else:
                        if gly not in gly_result:
                            gly_result =  gly_result+'//'+ gly

    return motifs_result_list, motifs_result,gly_result_list,gly_result
def output_node_motifs_annotation(filename,norm_MS2List,motifs_result_list_final,gly_result_list_final):
    # filename = 'output_list_network.xlsx'
    final_node = []
    for i in range(len(norm_MS2List)):
        one_MS2List = norm_MS2List[i][0]
        tR=one_MS2List[3]
        pm=one_MS2List[4]
        # print([i,pm,tR,one_MS2List[0],one_MS2List[1],one_MS2List[2],motifs_result_list_final[i],gly_result_list_final[i]])
        final_node.append([i,pm,tR,one_MS2List[0],one_MS2List[1],one_MS2List[2],motifs_result_list_final[i],gly_result_list_final[i]])

    df = pd.DataFrame(final_node, columns=['node_ID', 'node_precursor_mass','node_RT','node_mgf','mark','metabolite_name','Sub_motif','Glycosyl_neutral_loss'])
    df.to_excel(filename, index=False)




def output_CGFs_annotation(final_edge_data_HCCAs,filename):

    df = pd.DataFrame(final_edge_data_HCCAs, columns=['node1_ID', 'node2_ID','diff_ms','cosine_score','group', 'ID_1','ID_2','Align_1','Align_2',\
                                                      'node1_precursor_mass','node1_RT','node1_mgf','node1_mark', \
                                                      'node2_precursor_mass','node2_RT','node2_mgf','node2_mark','Enzyme_catalyzed_reactions',\
                                                      'node1_Sub_motif','node1_Glycosyl_neutral_loss',\
                                                      'node2_Sub_motif','node2_Glycosyl_neutral_loss'])
    df.to_excel(filename, index=False)

def merge_CGFs_edge(CGFs_data,edge_data):

    final_edge_data = []
    for i in range(len(edge_data)):
        node1_id = edge_data[i][0]
        node2_id = edge_data[i][1]
        one_edge_data = []
        # one_edge_data.append( edge_data[i])
        one_edge_data.append(CGFs_data[node1_id][6])
        one_edge_data.append(CGFs_data[node1_id][7])

        one_edge_data.append(CGFs_data[node2_id][6])
        one_edge_data.append(CGFs_data[node2_id][7])


        final_edge_data.append(edge_data[i]+one_edge_data)
    return final_edge_data



# Read annotation files
def read_xlsx_file(file_name):
    # file_name=Characteristic_ions_file_name
    if '.csv' in file_name:
        df = pd.read_csv(file_name)
        # print('running--------')
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# CGFs nodes annotation
Motifs_file_name= r'D:\Python_code\script_code\annotation_file\Motifs_list.csv'
Motifs_data = read_xlsx_file(Motifs_file_name)
Glycosyl_file_name= r'D:\Python_code\script_code\annotation_file\Glycosyl_neutral_losses_list.csv'
Glycosyl_data = read_xlsx_file(Glycosyl_file_name)

# Normalize intensity
norm_MS2List = all_spectrum_normalize(all_MS2List)
motifs_result_list_final = []
gly_result_list_final = []
mz_threshold = 0.005
int_threshold = 0.2
for i in range(len(norm_MS2List)) :
    one_list = norm_MS2List[i]
    motifs_result_list, motifs_result,gly_result_list,gly_result = match_one_motifs(Motifs_data,one_list,Glycosyl_data)
    # if len(motifs_result)>0:
        # print(motifs_result)
    # if len(gly_result)>0:
        # print(gly_result)
    motifs_result_list_final.append(motifs_result)
    gly_result_list_final.append(gly_result)

node_motifs_annotation_file = 'output_node_motifs_annotation_info.xlsx'
output_node_motifs_annotation(node_motifs_annotation_file,norm_MS2List,motifs_result_list_final,gly_result_list_final)

# Read annotated edges result
edge_file_name = 'output_network_edge_annotation_result_align.xlsx'
edge_data = read_xlsx_file(edge_file_name)

# Read annotated nodes result
CGFs_file_name = 'output_node_motifs_annotation_info.xlsx'
CGFs_data = read_xlsx_file(CGFs_file_name)

# Merge annotated edges and nodes files
final_edge_data_CGFs = merge_CGFs_edge(CGFs_data,edge_data)

output_CGFs_annotation_filename = 'final_output_CGFs_annotation_result.xlsx'
output_CGFs_annotation(final_edge_data_CGFs,output_CGFs_annotation_filename)


print('CGFs network annotation Finished')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
