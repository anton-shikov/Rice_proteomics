#/usr/bin/python3.7

import argparse
import csv
from collections import defaultdict


def create_dict_from_list(parse_list, key_ind, *val_inds):
    """
        Creates list dict from parsing list
        Key refers to the asserted index (key_ind), value - by default takes all the indexes (including the key index)
    """
    parse_dict=defaultdict(list)
    for string in parse_list:
        if not val_inds:
            parse_dict[string[key_ind]]=string
        else:
            parse_dict[string[key_ind]]=[string[i] for i in range(len(string)) if i in val_inds]
    return(parse_dict)


def read_csv_to_list(in_file, headless=True, delim='\t'):
    """
        Reads csv file and returns list without header by default 
        If headless argument is false, parses the whole file
    """
    ret_list=list()
    with open(in_file,'r') as csv_file:
        my_reader = csv.reader(csv_file, delimiter=delim) 
        if headless:
            next(my_reader)
        for row in my_reader:
            ret_list.append(list(row))
    return(ret_list)


def write_csv(row_list,out_name,*header_strings : str, sep='\t'):
    """
       A universal function for writing lists to csv-files
       If input is list of lists uses writerows function else iteratively writes 
       If strings for header are specified, writes header to the output, saves headless table otherwise
    """
    with open(out_name,'w',newline='') as result_file:
        wr = csv.writer(result_file, delimiter=sep)
        if header_strings:
            wr.writerow([name for name in header_strings])
        if type(row_list[0]) is list:
            wr.writerows(row_list)
        else:
            for row in row_list:
                wr.writerow([row])

