import sys
import csv
from collections import defaultdict
import os
import argparse
from IO_lib import write_csv, read_csv_to_list
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pygraphviz
import statistics as st
import math

# root - "#CD853F" - brown
# shoot - '#EEE8AA' - light yellow

def mean_in_list(calc_list):
    return(sum(calc_list)/len(calc_list))


def make_spots_dict(spot_tab, organ):

    ret_spots_dict=defaultdict(list)
    spot_tab_csv=read_csv_to_list(spot_tab, headless=True)
    pre_spot_dict=defaultdict(dict)

    for row in spot_tab_csv:
        if organ=='root':
            ind=13
        else:
            ind=11
        if row[ind]!='-':
            row[ind]=row[ind].replace(' ','')

            if organ=='root':
                if row[ind] not in pre_spot_dict:
                    pre_spot_dict[row[ind]]['anoxia']=[]
                    pre_spot_dict[row[ind]]['control']=[]
                    pre_spot_dict[row[ind]]['reaeration']=[]
                pre_spot_dict[row[ind]]['anoxia'].extend([float(x) for x in row[2:5]])
                pre_spot_dict[row[ind]]['control'].extend([float(x) for x in row[5:8]])
                pre_spot_dict[row[ind]]['reaeration'].extend([float(x) for x in row[8:11]])
                pre_spot_dict[row[ind]]['annotation']=row[12]

            elif organ=='shoot':
                if row[ind] not in pre_spot_dict:
                    pre_spot_dict[row[ind]]['anoxia']=[]
                    pre_spot_dict[row[ind]]['control']=[]
                    pre_spot_dict[row[ind]]['reaeration']=[]
                pre_spot_dict[row[ind]]['anoxia'].extend([float(x) for x in row[2:5]])
                pre_spot_dict[row[ind]]['control'].extend([float(x) for x in row[5:7]])
                pre_spot_dict[row[ind]]['reaeration'].extend([float(x) for x in row[7:9]])
                pre_spot_dict[row[ind]]['annotation']=row[10]
     
    for accession in pre_spot_dict:
        accession=accession.replace(' ', '') 

        ret_spots_dict[accession]=[mean_in_list(pre_spot_dict[accession]['anoxia']), 
                                   mean_in_list(pre_spot_dict[accession]['control']),
                                   mean_in_list(pre_spot_dict[accession]['reaeration']),
                                   pre_spot_dict[accession]['annotation'],
                                   st.stdev(pre_spot_dict[accession]['anoxia'])/math.sqrt(len(pre_spot_dict[accession]['anoxia'])),
                                   st.stdev(pre_spot_dict[accession]['control'])/math.sqrt(len(pre_spot_dict[accession]['anoxia'])),
                                   st.stdev(pre_spot_dict[accession]['reaeration'])/math.sqrt(len(pre_spot_dict[accession]['anoxia']))]
    return(ret_spots_dict)

def process_TFs_names(TF):
    if TF=='Myb/SANT; MYB':
        TF='Myb/SANT'
    elif TF=='Dof;GATA':
        TF='Dof'
    elif TF=='Homeodomain; HD-ZIP':
        TF='HD-ZIP'
    elif TF=='Homeodomain; WOX':
        TF='WOX'
    elif TF=='Homeodomain; ZF-HD':
        TF='ZF-HD'
    elif TF=='NAC; NAM':
        TF='NAM'
    elif TF=='AP2; B3; RAV':
        TF='AP2/B3'
    elif TF=='Homeodomain; TALE':
        TF='TALE'
    elif TF=='B3; bZIP;bZIP':
        TF='bZIP'
    elif TF=='MADS box; MIKC; M-type':
        TF='MIKC'
    elif TF=='MADS box; MIKC':
        TF='MIKC'
    elif TF=='Myb/SANT; G2-like':
        TF='Myb/SANT'
    elif TF=='AP2; ERF;ERF':
        TF='ERF'
    elif TF=='AP2; ERF':
        TF='ERF'
    elif TF=='B3; ARF;ARF':
        TF='ARF'
    elif TF=='B3; ARF':
        TF='ARF'
    elif TF=='bZIP; Homeodomain; HD-ZIP':
        TF='bZIP'

    return(TF)

def signif_dict_to_annotations(signif_dict):
    resulting_csv=[['Accession', 'Anoxia', 'Control', 'Re-aeration', 'Annotation', 'err_Anoxia', 'err_Control', 'err_Re-aeration','Type']]
    res_dict=dict()
    for key in signif_dict:
        if signif_dict[key][0]> signif_dict[key][1] and signif_dict[key][0]> signif_dict[key][2]:
            resulting_csv.append([key]+ signif_dict[key]+['Anoxia'])
            res_dict[key]='Anoxia'
        elif signif_dict[key][1]> signif_dict[key][0] and signif_dict[key][1]> signif_dict[key][2]:
            resulting_csv.append([key]+ signif_dict[key]+['Control'])
            res_dict[key]='Control'
        elif signif_dict[key][2]> signif_dict[key][0] and signif_dict[key][2]> signif_dict[key][1]:
            resulting_csv.append([key]+ signif_dict[key]+['Re-aeration'])
            res_dict[key]='Re-aeration'
        else:
            print('unknown')

    #Change accessions for alternative equal names
    if 'XP_015617521.1' in res_dict:
        res_dict['XP_015617521.1']='Re-aeration'
    if 'BAC79528.1' in res_dict:
        res_dict['BAC79528.1']='Re-aeration'

    write_csv(resulting_csv, 'mean_err_signif_points_root.csv')
    return(res_dict)
    
def check_accession(accession, spots_dict, condition):
    if accession in spots_dict:
        if spots_dict[accession]==condition:
            return(True)
        else:
            return(False)
    elif accession+'_dub' in spots_dict:
        if spots_dict[accession+'_dub']==condition:
            return(True)
        else:
            return(False)
     
def summarize_sites_signals(root_dir, shoot_dir, spots_dict_root, spots_dict_shoot, annot_flag=False):

    TF_drawing_dict=defaultdict(dict)

    for prom_tab in os.listdir(root_dir):
        if '~' not in prom_tab:
            if annot_flag:
                if prom_tab.replace('.csv','') not in spots_dict_root.keys():
                    continue
            tab_rel_name = os.path.join(os.path.realpath(root_dir),prom_tab)
            TF_tab_csv=read_csv_to_list(tab_rel_name, headless=True)
            prom_tab=prom_tab.replace('.csv','')
            for row in TF_tab_csv:
                if row[1]!='(Motif sequence only)' and row[1]!='(Others)' and row[1]!='Alpha-amylase':
                    TF_name = process_TFs_names(row[1])
                    if TF_name not in TF_drawing_dict:
                        TF_drawing_dict[TF_name]=dict()
                    if row[2] not in TF_drawing_dict[TF_name]:
                        TF_drawing_dict[TF_name][row[2]]=defaultdict(list)
                        TF_drawing_dict[TF_name][row[2]]['root_genes']=[]
                        TF_drawing_dict[TF_name][row[2]]['root_calls']=[]
                        TF_drawing_dict[TF_name][row[2]]['shoot_genes']=[]
                        TF_drawing_dict[TF_name][row[2]]['shoot_calls']=[]
                    TF_drawing_dict[TF_name][row[2]]['root_genes']=list(set(list(TF_drawing_dict[TF_name][row[2]]['root_genes']) +[prom_tab]))
                    TF_drawing_dict[TF_name][row[2]]['root_calls'].append(prom_tab)


    for prom_tab in os.listdir(shoot_dir):
        if '~' not in prom_tab:
            if annot_flag:
                if prom_tab.replace('.csv','') not in spots_dict_shoot.keys():
                    continue



            tab_rel_name = os.path.join(os.path.realpath(shoot_dir),prom_tab)
            TF_tab_csv=read_csv_to_list(tab_rel_name, headless=True)
            prom_tab=prom_tab.replace('.csv','')

            for row in TF_tab_csv:
                if row[1]!='(Motif sequence only)' and row[1]!='(Others)' and row[1]!='Alpha-amylase':
                    TF_name = process_TFs_names(row[1])
                    if TF_name not in TF_drawing_dict:
                        TF_drawing_dict[TF_name]=dict()
                    if row[2] not in TF_drawing_dict[TF_name]:
                        TF_drawing_dict[TF_name][row[2]]=defaultdict(list)
                        TF_drawing_dict[TF_name][row[2]]['root_genes']=[]
                        TF_drawing_dict[TF_name][row[2]]['root_calls']=[]
                        TF_drawing_dict[TF_name][row[2]]['shoot_genes']=[]
                        TF_drawing_dict[TF_name][row[2]]['shoot_calls']=[]
                    TF_drawing_dict[TF_name][row[2]]['shoot_genes']=list(set(list(TF_drawing_dict[TF_name][row[2]]['shoot_genes']) +[prom_tab]))
                    TF_drawing_dict[TF_name][row[2]]['shoot_calls'].append(prom_tab)

    TF_coords_tab=[['TF','Site', 'Root_genes', 'Root_calls', 'Shoot_genes', 'Shoot_calls']]

    for TF in TF_drawing_dict:
        for site in TF_drawing_dict[TF]:
            TF_coords_tab.append([TF, site, len(TF_drawing_dict[TF][site]['root_genes']),
                  len(TF_drawing_dict[TF][site]['root_calls']),len(TF_drawing_dict[TF][site]['shoot_genes']), len(TF_drawing_dict[TF][site]['shoot_calls'])])


    write_csv(TF_coords_tab, 'rice_promotor_sites_coords_signif.tsv')

    TF_sites_genes_num=[['TF', 'Root_genes', 'Root_sites', 'Shoot_genes','Shoot_sites', 'All_genes', 'All_sites']]
    TF_sites_signif_num=[['TF','Root_genes', 'Shoot_genes','All_genes', 'Type', 'Root_list', 'Shoot_list', 'Root_sites', 'Shoot_sites', 'All_sites']]

    for TF in TF_drawing_dict:
        TF_sites_roots=0
        TF_genes_roots=[]

        TF_sites_shoots=0
        TF_genes_shoots=[]

        TF_sites_roots_anoxia=0
        TF_sites_roots_control=0
        TF_sites_roots_reaeration=0

        TF_sites_shoots_anoxia=0
        TF_sites_shoots_control=0
        TF_sites_shoots_reaeration=0
        
        for site in TF_drawing_dict[TF]:        
            for gene in TF_drawing_dict[TF][site]['shoot_calls']:
                TF_sites_shoots_anoxia+=int(gene in [key for key in spots_dict_shoot if spots_dict_shoot[key]=='Anoxia'])
                TF_sites_shoots_control+=int(gene in [key for key in spots_dict_shoot if spots_dict_shoot[key]=='Control'])
                TF_sites_shoots_reaeration+=int(gene in [key for key in spots_dict_shoot if spots_dict_shoot[key]=='Re-aeration'])

            for gene in TF_drawing_dict[TF][site]['root_calls']:
                TF_sites_roots_anoxia+=int(gene in [key for key in spots_dict_root if spots_dict_root[key]=='Anoxia'])
                TF_sites_roots_control+=int(gene in [key for key in spots_dict_root if spots_dict_root[key]=='Control'])
                TF_sites_roots_reaeration+=int(gene in [key for key in spots_dict_root if spots_dict_root[key]=='Re-aeration'])

            
            TF_genes_roots.extend(TF_drawing_dict[TF][site]['root_genes']) 
            TF_genes_shoots.extend(TF_drawing_dict[TF][site]['shoot_genes']) 

            TF_sites_roots+=len(TF_drawing_dict[TF][site]['root_calls'])
            TF_sites_shoots+=len(TF_drawing_dict[TF][site]['shoot_calls'])
        TF_genes_roots=list(set(TF_genes_roots))
        TF_genes_shoots=list(set(TF_genes_shoots))
        TF_signif_roots=[el.replace('.csv', '') for el in TF_genes_roots if el.replace('.csv', '') in list(spots_dict_root.keys())]
        TF_signif_shoots=[el.replace('.csv', '') for el in TF_genes_shoots if el.replace('.csv', '') in list(spots_dict_shoot.keys())]

        if len(TF_signif_shoots) + len(TF_signif_roots)>0:
            anoxia_root=[el for el in TF_signif_roots if check_accession(el,spots_dict_root,'Anoxia')]
            control_root=[el for el in TF_signif_roots if check_accession(el,spots_dict_root,'Control')]
            reaeration_root=[el for el in TF_signif_roots if check_accession(el,spots_dict_root,'Re-aeration')]

            anoxia_shoot=[el for el in TF_signif_shoots if check_accession(el,spots_dict_shoot,'Anoxia')]
            control_shoot=[el for el in TF_signif_shoots if check_accession(el,spots_dict_shoot,'Control')]
            reaeration_shoot=[el for el in TF_signif_shoots if check_accession(el,spots_dict_shoot,'Re-aeration')] 

            if len(anoxia_root) + len(anoxia_shoot)>0:
                TF_sites_signif_num.append([TF,len(anoxia_root), len(anoxia_shoot), len(anoxia_root)+len(anoxia_shoot), 'Anoxia', '; '.join(anoxia_root), '; '.join(anoxia_shoot),
                                           TF_sites_roots_anoxia, TF_sites_shoots_anoxia, TF_sites_shoots_anoxia+TF_sites_roots_anoxia])
            if len(control_root) + len(control_shoot)>0:
                TF_sites_signif_num.append([TF,len(control_root), len(control_shoot), len(control_root)+len(control_shoot), 'Control', '; '.join(control_root), '; '.join(control_shoot),
                                           TF_sites_roots_control, TF_sites_shoots_control, TF_sites_shoots_control+TF_sites_roots_control])
            if len(reaeration_root) + len(reaeration_shoot)>0:
                TF_sites_signif_num.append([TF,len(reaeration_root), len(reaeration_shoot), len(reaeration_root)+len(reaeration_shoot),'Re-aeration', '; '.join(reaeration_root), '; '.join(reaeration_shoot),TF_sites_roots_reaeration, TF_sites_shoots_reaeration, TF_sites_shoots_reaeration+TF_sites_roots_reaeration])


        TF_sites_genes_num.append([TF, len(TF_genes_roots), TF_sites_roots, len(TF_genes_shoots), TF_sites_shoots,len(TF_genes_roots)+ len(TF_genes_shoots), TF_sites_roots+TF_sites_shoots ])


    write_csv(TF_sites_signif_num, 'rice_promotor_sites_change_signif.tsv')
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A script for summarizing promotors results')
    parser.add_argument('-r', '--root', dest='root_tab', help='the path to the root spots table',
                        type=str)
    parser.add_argument('-s', '--shoot', dest='shoot_tab', help='the path to the shoot spots table',
                        type=str)
    parser.add_argument('-sr', '--sigroot', dest='sigroot_tab', help='the path to the root signifficant spots table',
                        type=str)
    parser.add_argument('-ss', '--sigshoot', dest='sigshoot_tab', help='the path to the shoot signifficant spots table',
                        type=str)
    parser.add_argument('-rd', '--root_d', dest='root_dir', help='the path to the root TF sites',
                        type=str)
    parser.add_argument('-sd', '--shoot_d', dest='shoot_dir', help='the path to the shoot TF sites',
                        type=str)
    args = parser.parse_args()

    spots_dict_root = make_spots_dict(args.sigroot_tab, 'root')
    spots_dict_shoot = make_spots_dict(args.sigshoot_tab, 'shoot')

    signif_dict_annot_root=signif_dict_to_annotations(spots_dict_root)
    signif_dict_annot_shoot=signif_dict_to_annotations(spots_dict_shoot)

    spots_dict_root_all = make_spots_dict(args.root_tab, 'root')
    spots_dict_shoot_all = make_spots_dict(args.shoot_tab, 'shoot')


    summarize_sites_signals(args.root_dir, args.shoot_dir, signif_dict_annot_root, signif_dict_annot_shoot, annot_flag=True)



