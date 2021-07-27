# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 18:41:11 2020

@author: pjotr
"""
import numpy as np
from kpal.klib import Profile
from Bio import SeqIO

fname = 'C:/Users/pjotr/Documents/Python Scripts/Seq/64s/64s.fasta' # 25s/assembly.fa'

# labels from scaffolds
sc = [record.id for record in SeqIO.parse(fname, 'fasta')]

#%% kmer profiles
profs = Profile.from_fasta_by_record(open(fname), 4)
p_list = [ p for p in profs]
n = len(p_list)
p_mat = np.zeros((n,256))

for i in range(n):     
    p = p_list[i]
    p.balance()
    p_mat[i] = p.counts
    
i_vec = -np.ones(256)
for i in range(256): 
    j = p.reverse_complement(i)
    if i_vec[j]==-1: 
        i_vec[i] = i

inds = i_vec[i_vec > -1]       
x_raw_mat = p_mat[:, i_vec > -1]   

#%% labels from fasta headers (provisional ...) 

#import csv
#
#fname ='C:/Users/pjotr/Documents/Python Scripts/Seq/25s/member.txt'
# 
#with open(fname) as f:
#    reader = csv.reader(f, delimiter="\t")
#    d = list(reader)
##    
#print(d[0][1]) 
#n = len(d)
#
#dt = [([int(s) for s in d[i][0].split('-') if s.isdigit()][0],d[i][1]) for i in range(n)]
#dt.sort(key=lambda tup: tup[0])
#
#lab_list = [dt[i][1] for i in range(n)]
#
#lab_set = set(lab_list)
#  
#dic_lab = dict(zip(lab_set, range(0,len(lab_set))))
#
#lab_vec = np.zeros(n, dtype=np.int)
#for i in range(n):
#    lab_vec[i] = dic_lab[lab_list[i]]
#
#fname ='C:/Users/pjotr/Documents/Python Scripts/Seq/25s/My.depth'
#
#labs_ref = lab_vec

