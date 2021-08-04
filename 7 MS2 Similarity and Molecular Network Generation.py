# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 2021
Part 7: MS2 Similarity and Molecular Network Generation
@author: Li Zimu
"""

import pandas as pd

import Cosine

import networkx as nx

# Calculate MS/MS similarity
def MS2_similarity_score_align(spec2,MS2List,similarity_threshold,ms1_tolerance,ms2_tolerance,ID,num_match_ion):
    final_network = []
    num_spec = len(MS2List)

    for i in range(0,num_spec):
        for j in range(0,num_spec):
            if i<=j:
                continue
            else:
                MS2_ID_1 = i
                MS2_ID_2 = j
                pm1=MS2List[MS2_ID_1][4]
                pm2=MS2List[MS2_ID_2][4]

                MS2_similarity,reported_alignments=Cosine.score_alignment(spec2[MS2_ID_1],spec2[MS2_ID_2],pm1,pm2,ms1_tolerance,ms2_tolerance,max_charge_consideration=1) #二级质谱对齐
#
#                print(MS2_similarity)
#                print(len(reported_alignments))
#                reported_alignments
                if len(reported_alignments)>=num_match_ion:
                    if MS2_similarity>=similarity_threshold:
                        output = [str(ID[MS2_ID_1]),str(ID[MS2_ID_2]),str(pm1-pm2),MS2_similarity]

                        final_network.append(output)
    return final_network

# Filter features
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


def output_node_info(filename,final_node):
#    filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_node, columns=['node_ID', 'node_precursor_mass','node_RT','node_mgf','mark','metabolite_name'])
    df.to_excel(filename, index=False)
def output_node_info_Align(filename,final_node):
#    filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_node, columns=['node_ID', 'node_precursor_mass','node_RT','node_mgf','mark','metabolite_name','ID','Align_ID'])
    df.to_excel(filename, index=False)
def output_int_network(filename,final_network):
#    filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_network, columns=['node1_ID', 'node2_ID','diff_ms','cosine_score'])
    df.to_excel(filename, index=False)

def output_int_network_align(filename,final_network,all_MS2List):

    for i in range(len(final_network)):
        network_edge = final_network[i]
        if len(all_MS2List[int(network_edge[0])])>6:
            ID_1 = all_MS2List[int(network_edge[0])][6]
            Align_1 = all_MS2List[int(network_edge[0])][7]
        else:
            ID_1 = network_edge[0]
            Align_1 = network_edge[0]
        if len(all_MS2List[int(network_edge[1])])>6:
            ID_2 = all_MS2List[int(network_edge[1])][6]
            Align_2 = all_MS2List[int(network_edge[1])][7]
        else:
             ID_2 = network_edge[1]
             Align_2 = network_edge[1]
        final_network[i].append(ID_1)
        final_network[i].append(ID_2)
        final_network[i].append(Align_1)
        final_network[i].append(Align_2)

#    filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_network, columns=['node1_ID', 'node2_ID','diff_ms','cosine_score','ID_1','ID_2','Align_1','Align_2'])
    df.to_excel(filename, index=False)

    return final_network



def filter_top_k(G, top_k):
    print("Filter Top_K", top_k)
    #Keeping only the top K scoring edges per node
    print("Starting Numer of Edges", len(G.edges()))

    node_cutoff_score = {}
    for node in G.nodes():
        node_edges = G.edges((node), data=True)
        node_edges = sorted(node_edges, key=lambda edge: edge[2]["cosine_score"], reverse=True)

        edges_to_delete = node_edges[top_k:]
        edges_to_keep = node_edges[:top_k]

        if len(edges_to_keep) == 0:
            continue

        node_cutoff_score[node] = edges_to_keep[-1][2]["cosine_score"]
    edge_to_remove = []
    for edge in G.edges(data=True):
        cosine_score = edge[2]["cosine_score"]
        threshold1 = node_cutoff_score[edge[0]]
        threshold2 = node_cutoff_score[edge[1]]

        if cosine_score < threshold1 or cosine_score < threshold2:
            edge_to_remove.append(edge)

    for edge in edge_to_remove:
        G.remove_edge(edge[0], edge[1])

    print("After Top K Mutual", len(G.edges()))


def filter_component_additive(G, max_component_size):
    if max_component_size == 0:
        return

    all_edges = list(G.edges(data=True))
    G.remove_edges_from(list(G.edges))

    all_edges = sorted(all_edges, key=lambda x: x[2]["cosine_score"], reverse=True)

    for edge in all_edges:
        G.add_edges_from([edge])
        largest_cc = max(nx.connected_components(G), key=len)

        if len(largest_cc) > max_component_size:
            G.remove_edge(edge[0], edge[1])
    return G
def get_edges_of_component(G, component):
    component_edges = {}
    for node in component:
        node_edges = G.edges((node), data=True)
        for edge in node_edges:
            if edge[0] < edge[1]:
                key = edge[0] + "-" + edge[1]
            else:
                key = edge[1] + "-" + edge[0]
#            key = edge[0] + "-" + edge[1]
            component_edges[key] = edge

    component_edges = component_edges.values()
    return component_edges



def output_graph_align(G, filename):
    output_file = open(filename, "w")
    #Outputting the graph
    component_index = 0
    graph_info = []
    for component in nx.connected_components(G):
        component_index += 1
        for edge in get_edges_of_component(G, component):
            output_list = []
#            
#            if int(edge[0]) < int(edge[1]):
#                output_list.append(edge[0])
#                output_list.append(edge[1])
#            else:
#                output_list.append(edge[1])
#                output_list.append(edge[0])
            output_list.append(edge[2]["node1"])
            output_list.append(edge[2]["node2"])
            output_list.append(edge[2]["mass_difference"])
#            output_list.append(edge[2]["property1"])
            output_list.append(str(edge[2]["cosine_score"]))
#            output_list.append(str(edge[2]["explained_intensity"]))
            output_list.append(str(component_index))
            output_list.append(edge[2]["ID1"])
            output_list.append(edge[2]["ID2"])
            output_list.append(edge[2]["Align1"])
            output_list.append(edge[2]["Align2"])
            output_file.write("\t".join(output_list) + "\n")
            graph_info.append(output_list)
    output_file.close()
    return graph_info



# Output nodes information
def match_node_info(MS2List,loc):

    node1_mgf = MS2List[int(loc)][0]
    mark = MS2List[int(loc)][1]
    tR=MS2List[int(loc)][3]
    pm=MS2List[int(loc)][4]

    node_info = [pm,tR,node1_mgf,mark]
    return node_info
def output_node_edge_info_align(filename,final_node):
#    filename = 'output_list_network.xlsx'
    df = pd.DataFrame(final_node, columns=['node1_ID', 'node2_ID','diff_ms','cosine_score','group','ID_1','ID_2','Align_1','Align_2', 'node1_precursor_mass','node1_RT','node1_mgf','node1_mark', 'node2_precursor_mass','node2_RT','node2_mgf','node2_mark'])
    df.to_excel(filename, index=False)

def output_network_node_edge_info_align(graph_info, MS2List,filename):

    num_graph = len(graph_info)
    graph_info_final = []
    for i in range(num_graph):
        graph_info_final.append(graph_info[i]+match_node_info(MS2List,int(graph_info[i][0]))+match_node_info(MS2List,int(graph_info[i][1])))

    output_node_edge_info_align(filename,graph_info_final)
    return graph_info_final

def loading_network_align_new(final_network):
    row_count = len(final_network)
    edge_property_map = {}
    edge_object_list = []
    intermediate_graph_nodes = set()
    intermediate_edges_to_add = []
    for i in range(row_count):
        if final_network[i][3]<0.6:
            continue
        else:

            edge_object = {}
            edge_object["node1"] = final_network[i][0]
            edge_object["node2"] = final_network[i][1]
            edge_object["mass_difference"] = final_network[i][2]
            edge_object["property1"] = []
            edge_object["cosine_score"] = float(final_network[i][3])
            edge_object["explained_intensity"] = []
            edge_object["component"] = -1
            edge_object["EdgeType"] = "Cosine"
            edge_object["EdgeAnnotation"] = []
            edge_object["EdgeScore"] = float(final_network[i][3])
            edge_object["ID1"] = final_network[i][4]
            edge_object["ID2"] = final_network[i][5]
            edge_object["Align1"] = final_network[i][6]
            edge_object["Align2"] = final_network[i][7]
            edge_key = final_network[i][0] + "-" + final_network[i][1]

            edge_property_map[edge_key] = edge_object

            intermediate_graph_nodes.add(edge_object["node1"])
            intermediate_graph_nodes.add(edge_object["node2"])

            intermediate_edges_to_add.append((edge_object["node1"], edge_object["node2"], edge_object))
    G=nx.MultiGraph()
    G.add_nodes_from(intermediate_graph_nodes)
    G.add_edges_from(intermediate_edges_to_add)

    return G




# Caculate MS/MS similarity
# mgf file processing
all_MS2List = MS2List_seed+MS2List
spec2,final_node,ID = MS2List2Spec(all_MS2List,6)

# Output all nodes information
output_node_info_Align('output_mgf_node_info_Silk_CGFs_POS.xlsx',final_node)
print('Merge the database file and the experiment file, Generate and output the nodes information of the network, Finished')


# Caculate MS/MS similarity
similarity_threshold = 0.85
ms1_tolerance_threshold = 0.01
ms2_tolerance_threshold = 0.02
num_match_ion = 7

import time
start =time.clock()
print('Start calculate the similarity of each pair MS/MS')

final_network = MS2_similarity_score_align(spec2,all_MS2List,similarity_threshold,ms1_tolerance_threshold,ms2_tolerance_threshold,ID,num_match_ion)
print('Calculate the similarity of each pair MS/MS Finished')
end = time.clock()
print('Running time: %s min'%((end-start)/60))
# Output initial network information
output_int_network('output_initiall_network_info.xlsx',final_network)
print('Output network information Finished')

# Output network information, add Align_ID
final_network = output_int_network_align('output_initiall_network_info_align.xlsx',final_network,all_MS2List)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Generat network and filter
G =  loading_network_align_new(final_network)
top_k = 50


filter_top_k(G, top_k)
print('filter_top_k Finished')
max_component_size = 500
G = filter_component_additive(G, max_component_size)
print('filter_component_additive Finished')
# Output filtered network information
graph_info = output_graph_align(G, 'output_filter_network_info_align.xlsx')

graph_info_final = output_network_node_edge_info_align(graph_info,all_MS2List, 'output_filter_network_node_edge_info_align.xlsx')

print('Network nodes filter Finished')
len(graph_info_final)
