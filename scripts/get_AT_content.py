#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
from Bio.SeqRecord import SeqRecord
from IO_lib import write_csv, read_csv_to_list

def calculate_AT(in_file):
    for record in SeqIO.parse(in_file,"fasta"):
        seq_str=str(record.seq).upper()
        AT_content=(seq_str.count('A')+seq_str.count('T'))/len(seq_str)
        print(record.id.split(':')[0], AT_content, sep='\t')

def make_coverage_hist(in_file):
    out_csv=[['Site','AT_content']]
    AT_dict=defaultdict(dict)
    for i in range(500):
        AT_dict[i]=''

    for record in SeqIO.parse(in_file,"fasta"):
        seq_str=str(record.seq).upper()

        for i in range(500):
            AT_dict[i]+=(seq_str[i])
    for site in AT_dict:
        AT_content=(AT_dict[site].count('A')+AT_dict[site].count('T'))/len(AT_dict[site])
        out_csv.append([site, AT_content])
    write_csv(out_csv,'root_AT_content_per_site')
        
              

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculates AT-content for the fasta file')
    parser.add_argument('-f', '--fas', dest='fas_file', help='the file with the aligned sequences',
                        type=str)
    args = parser.parse_args()

    #calculate_AT(args.fas_file)
    make_coverage_hist(args.fas_file)
