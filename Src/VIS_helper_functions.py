#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:10:28 2023

@author: weichan
"""

import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

   
def chunks(lst, n):
    """
    Returns list of n-sized lists in the nested format: [[],[]]
    """
    collection = []
    for i in range(0, len(lst), n):
        collection.append(lst[i:i + n])
    return collection

def fragmentation_fasta(fasta, equal_fragments, outfilename):
    """
    Splits FASTA file into equal_sized fragments and stores it as a single FASTA file with numbered fragments as ids
    """
    filename= outfilename #str(fasta) + "_fragmentation_" + str(equal_fragments) #the name of the output file
    with open(filename, 'w') as output_file:
        # opening given fasta file using the file path
        with open(fasta, 'r') as fasta_file:
            # extracting multiple data in single fasta file using biopython
            for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
                #split FASTA sequence into equally sized sub-lists with different id
                record_id = record.id
                allchunks = chunks(record.seq, equal_fragments)
                for i,chunk in enumerate(allchunks):
                    #print(">"+str(record_id)+"\n"+str(chunk) + "\n")
                    record_id_chunk = str(record_id) + '_%i' % (i)
                    output_file.write(">"+str(record_id_chunk)+"\n"+str(chunk) + "\n")            


#Hardcoded visualisations
                    
def read_length_distribution(fasta):
    """
    Plots read lengths of a FASTA file
    """
    lengths = map(len, SeqIO.parse(fasta, 'fasta'))
    plt.hist(np.log10(pd.Series(lengths)), color='black', bins=250)
    plt.ylabel('log10(frequency)')
    plt.xlabel("log10(bp)")
    plt.title('Read length distribution')
    plt.savefig('read_length_distribution.pdf')

#read_length_distribution(sys.argv[1])
#fragmentation_fasta(sys.argv[1], int(sys.argv[2]))
#print(chunks("EASTBEJDNEDENDNLKEDLKEDLEDLKMDLKMDKMD", 10))
def fragmentation_match_distribution(data):
    """
    Takes the vector fragments and plots histogram of their frequency in the alignment
    """
    blasted = pd.read_csv(data, sep="\t")
    blasted[['Vector', 'Fragment']] = blasted['SubjectID'].str.split('_', n=1, expand=True)
    print(blasted.head())
    blasted["Fragment"] = pd.to_numeric(blasted["Fragment"])
    plt.hist(blasted["Fragment"].sort_values(), color='black')
    #plt.xticks(range(0,106,10))
    plt.ylabel('Alignment Frequency')
    plt.xlabel("Vector Fragment")
    plt.title('Dummy 500 bp fragment distribution')
    plt.savefig('Dummy_500_fragmentation_distribution.pdf')

#fragmentation_match_distribution(sys.argv[1])
#plt.close()
def fragmentation_read_match_distribution(data):
    """
    Takes the read ids with matches and plots histogram of their frequency
    """
    blasted = pd.read_csv(data, sep="\t")
    plt.hist(blasted["QueryID"], color='black')
    plt.xticks(rotation=90)
    plt.ylabel('Read match Frequency')
    plt.xlabel("Read")
    plt.title('Dummy 500 bp read match fragment distribution')
    plt.savefig('Dummy_500_read_match_fragmentation_distribution.pdf', bbox_inches='tight')

#fragmentation_read_match_distribution(sys.argv[1])