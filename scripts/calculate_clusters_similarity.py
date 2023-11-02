#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
from IO_lib import write_csv, read_csv_to_list
import itertools

def transform_clust_tab(clust_tab):

    new_tab=list()
    new_tab.append(['Cluster', 'Intensity', 'Condition', 'Organ', 'Method', 'Group'])
    clust_res_csv=read_csv_to_list(clust_tab, headless=True)
    clust_dict=defaultdict(list)
 
    for row in clust_res_csv:
        key=':'.join([row[1]] + row[6:]).replace(' ', '|')
        clust_dict[key].append([float(x) for x in row[2:5]])

    for key in  clust_dict:
        control_intensities=[]
        anoxia_intensities=[]
        reaeration_intensities=[]
        print(key)
        for key_list in clust_dict[key]: #3_Shoot_h_clust_All
            control_intensities.append(key_list[0])
            anoxia_intensities.append(key_list[1])
            reaeration_intensities.append(key_list[2])

            control_mean=sum(control_intensities)/(len(control_intensities))
            anoxia_mean=sum(anoxia_intensities)/(len(anoxia_intensities))
            reaeration_mean=sum(reaeration_intensities)/(len(reaeration_intensities))

        new_tab.append([key.split(':')[0]] + [control_mean]+['Control']+key.split(':')[1:3]+[key.split(':')[3].replace('|',' ')])
        new_tab.append([key.split(':')[0]] + [anoxia_mean]+['Anoxia']+key.split(':')[1:3]+[key.split(':')[3].replace('|',' ')])
        new_tab.append([key.split(':')[0]] + [reaeration_mean]+['Reaeration']+key.split(':')[1:3]+[key.split(':')[3].replace('|',' ')])
    write_csv(new_tab, 'Clusters_intensity.tsv')
    
def create_clust_dict( organ, cluster, clust_tab):
    spots_dict=defaultdict(dict)
    for row in clust_tab:
        if row[6]==organ and row[7]==cluster:
            if row[8] not in spots_dict:
                spots_dict[row[8]]=dict()
            if int(row[1]) not in spots_dict[row[8]]:
                spots_dict[row[8]][int(row[1])]=list()
             
            spots_dict[row[8]][int(row[1])].append(row[0])
    return(spots_dict)
    
def clust_comparision(list1, list2, comp_type):
     if comp_type=='Min':
         return(len(set(list1).intersection(set(list2)))/min(len(list1), len(list2)))
     else:
         return(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2))))

def make_clust_comparision(clust_pair, spots_clust, comp_fun):
    max_count=max([len(spots_clust[clust_pair[0]]), len(spots_clust[clust_pair[1]])])
    indexes=[p for p in itertools.product(list(range(1, max_count+1)), repeat=2)]
    comp_list=list()
    for ind_tuple in indexes:
        try:
            clust_inter=clust_comparision(spots_clust[clust_pair[0]][ind_tuple[0]],spots_clust[clust_pair[1]][ind_tuple[1]], comp_type=comp_fun)
            comp_list.append(clust_inter)
        except:
            pass      
    return(sum([x for x in comp_list if x!=0])/max_count)

def make_similarity_heatmap(spot_dict):
    clust_pairs=list(itertools.combinations(spot_dict.keys(), 2))
    clust_heatmap=[['', 'All', 'Annotated', 'Significant', 'Annotated significant'],['All', 1, '','','' ],['Annotated',0,1,'',''],[ 'Significant',0,0,1,''],['Annotated significant',0,0,0,1]]
    ind_dict={'All|Annotated':(2,1), 'All|Significant':(3,1),'All|Annotated significant':(4,1),'Annotated|Significant':(3,2),'Annotated|Annotated significant':(4,2),'Significant|Annotated significant':(4,3)}
    for clust_pair in clust_pairs:
        similarity=make_clust_comparision(clust_pair, spot_dict, 'Min')
        clust_key=clust_pair[0]+'|'+clust_pair[1]
        clust_heatmap[ind_dict[clust_key][0]][ind_dict[clust_key][1]]=similarity
    return(clust_heatmap)
    
def make_inter_clust_comparision(group, spots_clust1, spots_clust2, comp_fun):
    max_count=max([len(spots_clust1[group]), len(spots_clust2[group])])
    indexes=[p for p in itertools.product(list(range(1, max_count+1)), repeat=2)]
    comp_list=list()
    for ind_tuple in indexes:
        try:
            clust_inter=clust_comparision(spots_clust1[group][ind_tuple[0]],spots_clust1[group][ind_tuple[1]], comp_type=comp_fun)
            comp_list.append(clust_inter)
        except:
            pass      
    return(sum([x for x in comp_list if x!=0])/max_count)

def compare_all_clusters(clust_tab):
    clust_res_csv=read_csv_to_list(clust_tab, headless=True)
    spots_dict_shoot_kmeans=create_clust_dict('Shoot', 'K-means', clust_res_csv)
    spots_dict_shoot_hclust=create_clust_dict('Shoot', 'h_clust', clust_res_csv)
    spots_dict_root_kmeans=create_clust_dict('Root', 'K-means', clust_res_csv)
    spots_dict_root_hclust=create_clust_dict('Root', 'h_clust', clust_res_csv)

    #write_csv(make_similarity_heatmap(spots_dict_shoot_kmeans), 'Shoot_k-means_clusters_similarity.tsv')
    #write_csv(make_similarity_heatmap(spots_dict_shoot_hclust), 'Shoot_hclust_clusters_similarity.tsv')
    #write_csv(make_similarity_heatmap(spots_dict_root_kmeans), 'Root_k-means_clusters_similarity.tsv')
    #write_csv(make_similarity_heatmap(spots_dict_root_hclust), 'Root_hclust_clusters_similarity.tsv')
    inter_clust_list=[['Shoot','All',make_inter_clust_comparision('All', spots_dict_shoot_kmeans, spots_dict_shoot_hclust, comp_fun='Jaccard')],
                      ['Shoot','Annotated',make_inter_clust_comparision('Annotated', spots_dict_shoot_kmeans, spots_dict_shoot_hclust, comp_fun='Jaccard')],
                      ['Shoot','Significant',make_inter_clust_comparision('Significant', spots_dict_shoot_kmeans, spots_dict_shoot_hclust, comp_fun='Jaccard')],
                     ['Shoot','Annotated significant',make_inter_clust_comparision('Annotated significant', spots_dict_shoot_kmeans, spots_dict_shoot_hclust, comp_fun='Jaccard')],
                      ['Root','All',make_inter_clust_comparision('All', spots_dict_root_kmeans, spots_dict_root_hclust, comp_fun='Jaccard')],
                      ['Root','Annotated',make_inter_clust_comparision('Annotated', spots_dict_root_kmeans, spots_dict_root_hclust, comp_fun='Jaccard')],
                      ['Root','Significant',make_inter_clust_comparision('Significant', spots_dict_root_kmeans, spots_dict_root_hclust, comp_fun='Jaccard')],
                      ['Root','Annotated significant',make_inter_clust_comparision('Annotated significant', spots_dict_root_kmeans, spots_dict_root_hclust, comp_fun='Jaccard')]]
    write_csv(inter_clust_list, 'Inter_clusters_similarity.tsv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluates clusters similarity for different clusterization')
    parser.add_argument('-c', '--clust', dest='clust_tab', help='the file with clusterization results',
                        type=str)
    args = parser.parse_args()

    #transform_clust_tab(args.clust_tab)
    compare_all_clusters(args.clust_tab)
