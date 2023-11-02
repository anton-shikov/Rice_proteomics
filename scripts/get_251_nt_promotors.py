#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
from Bio.SeqRecord import SeqRecord

def fix_fasta(in_file):
    new_recs=[]
    for record in SeqIO.parse(in_file,"fasta"):
        rec1=(SeqRecord(record.seq[0:251],id=record.id +'.1', description=record.description))
        rec2=(SeqRecord(record.seq[249:500],id=record.id +'.2', description=record.description))
        new_recs.append(rec1)
        new_recs.append(rec2)
    SeqIO.write(new_recs, in_file+'.251nt.fasta', 'fasta')
        
            
              

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='deletes gaps for the fasta file')
    parser.add_argument('-f', '--fas', dest='fas_file', help='the file with the aligned sequences',
                        type=str)
    args = parser.parse_args()

    fix_fasta(args.fas_file)
