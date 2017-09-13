#!/usr/bin/python
import sys
import csv
import re
import glob
import pandas as pd
import numpy as np

gene_check = ''
acc_final = ''
taxa_final = ''
count_species = 0

gfam_list = list()

# index = [0]*82
# numDF = pd.DataFrame(index)
# numDF = numDF.fillna(0)
# numPrint = pd.DataFrame(index)
# numPrint = numPrint.fillna(0)
# sample_num = 0


f_in= open(sys.argv[1],'r')
#f_out= open('gfam_collapsed2.tsv', 'w')

def f5(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

fileDF = list()
for line in f_in.readlines():
    spline = line.split('\t')
    header = spline[0:2]
    expr = spline[3:]
    expr1 = '\t'.join(expr)
    expr2 = expr1.replace('\n', '')
    if "Gene_Family" in header:
        print("Acc" + '\t' "Gene_Family" + '\t' + "Species_num" + '\t' + "Species" + '\t' + expr2)
        #samp_list = header.split('\t')
        sample_num = len(expr)
    #full = spline[1:]
    else:
    	things = line.split('\t')
        gfam = things[1]
        gfam_list.append(gfam)
        fileDF.append(line)

new_gfam = f5(gfam_list)


for i in range(0, len(new_gfam)):
    search_term = str(new_gfam[i]+"\t")
    matching = [s for s in fileDF if search_term in s]
    #matching = [s for s in fileDF if search_term in s]
    acc_list = list()
    taxa_list = list()
    orig_vec = np.array([0]*sample_num)
    if len(matching) > 1:
        for item in matching:
            items = item.strip('\n')
            stuff = items.split('\t')
            acc = stuff[0]
            gfam = stuff[1]
            taxa = stuff[2]
            nums = stuff[3:]
            if acc not in acc_list:
                acc_list.append(acc)
            if taxa not in taxa_list:
                taxa_list.append(taxa)
            data_num = [float(i) for i in nums]
            new_vec = np.array(data_num)
            sum_vec = orig_vec + new_vec
            orig_vec = sum_vec
        sum_list = list(sum_vec)
        sum_str = str(sum_list)
        sum_str = sum_str.replace(", ", "\t")
        sum_str = sum_str.replace("[", "")
        sum_str = sum_str.replace("]", "")
        print(','.join(acc_list) + '\t' + gfam + '\t' + str(len(taxa_list)) + '\t' + ','.join(taxa_list) + '\t' + sum_str)

    else:
        sing_item = '\t'.join(matching)
        sing_items = sing_item.strip('\n')
        sing_stuff = sing_items.split('\t')
        sing_acc = sing_stuff[0]
        sing_gfam = sing_stuff[1]
        sing_taxa = sing_stuff[2]
        sing_nums = sing_stuff[3:]
        print(sing_acc + '\t' + sing_gfam + '\t' + "1" + '\t' + sing_taxa + '\t' + '\t'.join(sing_nums))

f_in.close()
#f_out.close()

#### duplication issue in the matching step!!!!!!!!!!!!