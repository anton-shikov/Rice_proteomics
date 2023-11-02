#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, Counter
import csv
from Bio.SeqRecord import SeqRecord
from IO_lib import write_csv, read_csv_to_list

def AA_dist(in_file):
    mearged_str=''
    AA_perc_dict=defaultdict(list)
    for record in SeqIO.parse(in_file,"fasta"):
        seq_str=str(record.seq).upper()
        AA_sub_dict = Counter(seq_str)
        for AA in AA_sub_dict:
            AA_perc_dict[AA].append(AA_sub_dict[AA]/sum([AA_sub_dict[AA] for AA in AA_sub_dict]))
            #if AA=='K':
                #print(AA, AA_sub_dict[AA]/sum([AA_sub_dict[AA] for AA in AA_sub_dict]), record.description, sep='\t')
                #print(AA_sub_dict)
        for index, rec_element in enumerate(seq_str):
            if rec_element=='K':
                rec_name = ' '.join(record.description.split()[1:]).split(', chloroplastic')[0].split('[')[0].split(', mitochondrial ')[0].split(', cytosolic')[0].split(', chloroplast')[0].split('(chloroplast)')[0]
                #print(rec_name)
                print(index, rec_element, rec_name, len(seq_str), sep='\t')
    #print(*AA_perc_dict['E'], sep='\n')
    #print(AA_perc_dict)
    #write_csv([[key, sum(AA_perc_dict[key])/len(AA_perc_dict[key])] for key in AA_perc_dict],'lysin_percent_root_dist.tsv')
        
              

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculates the distribution of amino_acids for the fasta alignment')
    parser.add_argument('-f', '--fas', dest='fas_file', help='the file with the aligned sequences',
                        type=str)
    args = parser.parse_args()

    AA_dist(args.fas_file)
