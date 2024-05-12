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
#geneTable = pd.read_table('.\ProcessedData\\example_output\\selected_genes_ann.csv', sep=',')
geneTable = pd.read_table('.\ProcessedData\\example_output\\selected_genes_from_ec_ann.csv', sep=',')
###########################################################################
# load B.theta genome and feature table from NCBI
protein_names = geneTable["proteinID"]
gene_locustag = geneTable["locus_tag"]
gene_ec = geneTable["EC"]
# remove empty protein names
protein_empty = [type(i)==float for i in protein_names]
protein_names = [protein_names for (protein_names, remove) 
                 in zip(protein_names, protein_empty) if not remove]
gene_locustag = [gene_locustag for (gene_locustag, remove) 
                 in zip(gene_locustag, protein_empty) if not remove]
gene_ec = [gene_ec for (gene_ec, remove) 
                 in zip(gene_ec, protein_empty) if not remove]

# replace brackets in protein names
protein_names = [x.replace("['","") for x in protein_names]
protein_names = [x.replace("']","") for x in protein_names]


#############################################
# get protein sequences from NCBI by protein ID
# concatenate names with gene locustag and EC
# write all proteins to one sequence file
#proteinfilename = '.\ProcessedData\\example_output\\selected_genes_proteins.fasta'
proteinfilename = '.\ProcessedData\\example_output\\selected_genes_from_ec_proteins.fasta'
for i in range(len(protein_names)):
    record = Entrez.efetch(db="protein", id=protein_names[i], rettype="fasta", retmode="text")
    curseq = record.read()
    # add locus tag and EC to the name
    # from curseq reove first symbol (>) and last EOL
    curec = gene_ec[i]
    if type(curec)==float: # if curec is not string set it to empty string
        curec=''
    curseq = '>'+gene_locustag[i] + ' ' + curec + ' ' + curseq[1:-1]
    curseq_lines = curseq.splitlines()
    curseq_id = curseq_lines[0]
    curseq_seq = ''.join(curseq_lines[1:])
    if i==0:
        with open(proteinfilename, 'w') as f: #it is first file, rewrite previous
            f.write(curseq_id)
            f.write(curseq_seq)
            f.write('\n')
    else: #it is not first file, append
        with open(proteinfilename, 'a') as f:
            f.write(curseq_id) # write sequence without last EOL symbol
            f.write(curseq_seq)
            f.write('\n')
#############################################
