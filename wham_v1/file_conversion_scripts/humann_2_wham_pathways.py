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

f_in= open(sys.argv[1],'r')

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
    spline = line.split('\r')
    header = spline[0]
    if "#" in header:
        print('Acc' + '\t' + "Gene_Family" + '\t' + "Species" + header[9:len(header)-1])
        samp_list = header.split('\t')
        sample_num = len(samp_list)-1
    if "|" in line:
    	  things = line.split('\t')
        heading = things[0]
        expr = things[1:]
        expr_list = ','.join(expr)
        expr_vals1 = expr_list.replace(',', '\t')
        expr_vals = expr_vals1.replace('\n', '')
        expr_vals_print = expr_vals[0:len(expr_vals)]
	ex=str.split(expr_vals, '\t')
	ex_float=[]
	for item in ex:
		aa=float(item)
		ex_float.append(aa)
	if sum(item==0 for item in ex_float) > 1*len(ex_float):
		next
	else:
        	if "NO_NAME" not in line:
          		if "unclassified" not in line:
              			if "unknown" not in line:
                  			R=re.search('(.*)\:(.*)\|(.*)', heading)
                  			acc=R.group(1)
                  			gfam=R.group(2)
                  			taxa=R.group(3)
                  			print(acc + '\t' + gfam + '\t' + taxa + '\t' + expr_vals_print)
          		else:
              			if "unknown" not in line:
                  			things = line.split('\t')
                  			heading = things[0]
                  			R=re.search('(.*)\:(.*)\|(.*)', heading)
                  			acc=R.group(1)
                  			gfam=R.group(2)
                  			taxa=R.group(3)
                  			print(acc + '\t' + gfam + '\t' + taxa + '\t' + expr_vals_print)
            				#gfam_list.append(gfam)

f_in.close()
