import sys
import csv

trypsinfile = sys.argv[1]
trypsinlist=list()
parsefile = sys.argv[2]

with open (trypsinfile, 'r',newline='') as csvfile2:
    my_reader3 = csv.reader(csvfile2, delimiter='\t')
    for row in my_reader3:
        trypsinlist.append(row[0])


ret_rows=[]
with open (parsefile, 'r',newline='') as csvfile2:
    my_reader3 = csv.reader(csvfile2, delimiter='\t')
    for row in my_reader3:
        if row[0].split('.')[0]+'.' not in trypsinlist: 
            ret_rows.append(row)

with open (parsefile.split('.')[0]+'_filtered.mgf', 'w',newline='') as csvfile2:
    my_writer3 = csv.writer(csvfile2, delimiter='\t')
    for row in ret_rows:
        my_writer3.writerow(row)
