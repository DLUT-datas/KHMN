# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 2021
Part 8: Candidates Features Screening
@author: Li Zimu
"""
import pandas as pd
# Filter network by seed nodes
def filter_network_by_seed(graph_info_final):
    # Get the unique group
    group= []
    for info in graph_info_final:
        group.append(info[4])
	unque_group = list(set(group))

    filter_group = []
    for un_info in unque_group:

        loc = [i for i in range(len(group)) if group[i] == un_info]

        for j in loc:
            if 'seed_metabolite' in graph_info_final[j] or 'seed_metabolite_from_experiment' in graph_info_final[j]:
                filter_group.append(un_info)
    unque_filter_group = list(set(filter_group))

    return unque_filter_group
def output_node_edge_info_align(filename,final_node):
    # filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_node, columns=['node1_ID', 'node2_ID','diff_ms','cosine_score','group','ID_1','ID_2','Align_1','Align_2', 'node1_precursor_mass','node1_RT','node1_mgf','node1_mark', 'node2_precursor_mass','node2_RT','node2_mgf','node2_mark'])
    df.to_excel(filename, index=False)
def match_node_info(MS2List,loc):

    node1_mgf = MS2List[int(loc)][0]
    mark = MS2List[int(loc)][1]
    tR=MS2List[int(loc)][3]
    pm=MS2List[int(loc)][4]

    node_info = [pm,tR,node1_mgf,mark]
    return node_info
def output_filter_network_node_edge_info_align(graph_info, unque_filter_group,MS2List,filename):

    num_graph = len(graph_info)
    graph_info_final = []
    for i in range(num_graph):
        if graph_info[i][4] in unque_filter_group:
            graph_info_final.append(graph_info[i]+match_node_info(MS2List,int(graph_info[i][0]))+match_node_info(MS2List,int(graph_info[i][1])))

    output_node_edge_info_align(filename,graph_info_final)
    return graph_info_final






# Extract group information
unque_filter_group=filter_network_by_seed(graph_info_final)

# Output network with align_ID
output_filter_graph_info_final_name = 'output_seed_filter_network_info_align.xlsx'
filter_graph_info_final=output_filter_network_node_edge_info_align(graph_info,unque_filter_group,all_MS2List, output_filter_graph_info_final_name)
len(filter_graph_info_final)
print('Network filter by seeds Finished')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
