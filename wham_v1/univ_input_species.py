#!/usr/bin/python
import sys
import csv
import re
import glob

acc_check = ''
count_species=0
taxa_final = ''
old_tail = 0
header = ''
data = ''
acc=''
taxa=''
arg=''
output=''
gene=''

f_in= open(sys.argv[1],'r')
#f_out= open('gfam_collapsed.tsv', 'w')

for line in f_in.readlines():
    spline = line.split('\t')
    heading = spline[0:2]
    data = spline[3:]
    data1 = '\t'.join(data)
    data2 = data1.replace('\n', '')
    if 'Gene_Family' in line:
	print('Acc Number' + '\t' + 'Gene Family' + '\t' + 'Species' + '\t' + "Filler" + '\t' + data2)
    else:	
	acc = spline[0]
	gene = spline[1]
	taxa = spline[2]
	print(acc + '\t' + gene + '\t' + taxa + '\t' + acc + '\t' + data2)
#f_out.close()
f_in.close()

