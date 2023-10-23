# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 2023

@author: Maria
"""

from Bio import Entrez
import pandas as pd

#import re
#import subprocess
###########################################################################
    
Entrez.email = "m.kogadeeva@gmail.com"

###########################################################################
# read the strain table
geneTable = pd.read_table('table_selected_enzymes.csv', sep=',')

###########################################################################
# load B.theta genome and feature table from NCBI
protein_names = geneTable["proteinID"]
gene_locustag = geneTable["locus_tag"]
gene_ec = geneTable["EC.1"]


protein_info = []

#############################################
# get protein sequences from NCBI by protein ID
# concatenate names with gene locustag and EC
# write all proteins to one sequence file
proteinfilename = 'diet_selected_enzymes.fasta'
for i in range(len(protein_names)):
    record = Entrez.efetch(db="protein", id=protein_names[i], rettype="fasta", retmode="text")
    curseq = record.read()
    # add locus tag and EC to the name
    # from curseq reove first symbol (>) and last EOL
    curec = gene_ec[i]
    if type(curec)==float: # if curec is not string set it to empty string
        curec=''
    curseq = '>'+gene_locustag[i] + ' ' + curec + ' ' + curseq[1:-1]
    if i==0:
        with open(proteinfilename, 'w') as f: #it is first file, rewrite previous
            f.write(curseq) 
    else: #it is not first file, append
        with open(proteinfilename, 'a') as f:
            f.write(curseq) # write sequence without last EOL symbol
        
#############################################
