#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:10:28 2023

@author: weichan
"""

import sys
import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import collections
import seaborn as sns
#from pybedtools import BedTool
   
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

def insertion_normalisation(insertions,bases, n, outpath):
    """
    Uses the number of aligned bases, number of insertions, and scaling factor n for normalisation.
    Results in insertions per n aligned bases. 
    E.g. A Normalised count of 3 with the scaling factor 10000 means that there were on average 3 insertions detected for each 1kb aligned bases. 
    """
    n_insertions = sum(1 for _ in open(insertions)) #bed entries are by definition 1 per line
    print(n_insertions)
    with open(bases, 'r') as file:
        # Read the second line of the file
        next(file)
        second_line = file.readline()
        # Convert the content to an integer
        n_bases = int(second_line.strip())
        
    normalised = (n/n_bases) * n_insertions
    
    with open(outpath, 'w') as output_file:
        output_file.write("Number of aligned bases: " + str(n_bases) + "\n") 
        output_file.write("Number of Insertions: " + str(n_insertions) + "\n")
        output_file.write("Normalisation factor: " + str(n) + "\n")
        output_file.write("Insertions per " + str(n) + " bases: " + str(normalised))


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
def fragmentation_match_distribution(data, fragment_specifier, outpath):
    """
    Takes the vector fragments and plots histogram of their frequency in the alignment
    """
    blasted = pd.read_csv(data, sep="\t")
    blasted[['Vector', 'Fragment']] = blasted['SubjectID'].str.split('_', n=1, expand=True)
    blasted["Fragment"] = pd.to_numeric(blasted["Fragment"])
    freq = collections.Counter(blasted["Fragment"].sort_values())
    plt.bar(freq.keys(), freq.values(), color='black')
    #print(blasted["Fragment"].sort_values())
    if (10000/fragment_specifier) < 50:    
        plt.xticks(np.arange(0, (10000/fragment_specifier)+1))
    else:
        plt.xticks(np.arange(0, (10000/fragment_specifier)+1, step=(10000/fragment_specifier)/10))
    plt.ylabel('Alignment Frequency')
    plt.xlabel("Vector Fragment")
    plt.title(f'{fragment_specifier} bp fragment distribution')
    outfile = outpath + str("/") + f'{fragment_specifier}_fragmentation_distribution.png'
    plt.savefig(outfile)
    plt.close()

#fragmentation_match_distribution(sys.argv[1])
#plt.close()
def fragmentation_read_match_distribution(data, fragment_specifier, outpath):
    """
    Takes the read ids with matches and plots histogram of their frequency
    """
    blasted = pd.read_csv(data, sep="\t")
    freq = collections.Counter(blasted["QueryID"])
    plt.bar(freq.keys(), freq.values(), color='black')
    plt.xticks(rotation=90)
    plt.ylabel('Read match Frequency')
    plt.xlabel("Read")
    plt.title(f'{fragment_specifier} bp read match fragment distribution')
    outfile = outpath + str("/") + f'{fragment_specifier}_read_match_fragmentation_distribution.png'
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()

#fragmentation_read_match_distribution(sys.argv[1])
#under construction: Not callable from snakemake so far 
def group_read_venn(data1, data2, data3):
    """
    Manual use (not pipeline optimized yet: Creates Venn diagram showing the chromosome-specific density across the samples
    """
    set1 = set(pd.read_table(data1, header =None).iloc[:,0].tolist()) #for just the reads, header=None needs to be removed
    set2 = set(pd.read_table(data2, header=None).iloc[:,0].tolist())
    set3 = set(pd.read_table(data3, header=None).iloc[:,0].tolist())
    venn3([set1, set2, set3], ('ADS', 'CD123+', 'UTD'))
    plt.title('Chromosomal overlap in matches between samples')
    plt.savefig("Venn_100_reads_chromosome.png")

#group_read_venn(sys.argv[1], sys.argv[2], sys.argv[3])

def plot_bed_files_as_heatmap(bed_files, outfile):
    """
    Creates heatmap the chromosome-specific density of matches in BED across the samples
    """
    data = {}
    # Process each BED file
    for file_id, bed_file in enumerate(bed_files):
        # Read the BED file into a DataFrame
        df = pd.read_csv(bed_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Read_ID', 'X','XX'])

        # Count chromosome occurrences
        chromosome_counts = df['Chromosome'].value_counts()

        # Store the counts in the data dictionary
        head, tail = os.path.split(bed_file)
        data[tail.split(".")[0]] = chromosome_counts

    # Create a DataFrame from the data
    counts_df = pd.DataFrame(data).fillna(0)

    # Create a Seaborn heatmap
    plt.figure(figsize=(16, 9))
    sns.clustermap(counts_df, cmap="YlGnBu", annot=True, cbar_pos=(0, .2, .03, .4))
    plt.xlabel('File ID')
    plt.ylabel('Chromosome')
    plt.title('Chromosome Occurrences \n  in BED Files')
    plt.savefig(outfile, bbox_inches="tight")

def insertion_proximity(bedfile, binsize, outfile):
    """
    Takes a list of insertion sites (BED file) and adds a user-defined number of bp to each site (start and end). Returns new BED file.
    """
    df = pd.read_csv(bedfile, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "read"])
    df_modified = df.copy()
    df_modified['start'] = df['start'] - binsize
    df_modified['end'] = df['end'] + binsize
    df_modified.to_csv(outfile, sep = '\t', header = False, index = False)

def add_sequence_column(bed_file_path, fasta_file_path, output_bed_path):
    """
    Adds the fasta sequence to the respective bed gaps. This makes the following methylation calling step easier.
    """
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, usecols=[0,1,2], names=['Chromosome', 'Start', 'End'])
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Add a new column to the DataFrame with the corresponding sequence from the FASTA file
    bed_df['Sequence'] = [fasta_sequences[i].seq for i,n in enumerate(fasta_sequences)]
    
    bed_df.to_csv(output_bed_path, sep='\t', header=False, index=False)

'''
def methylation_in_insertion_proximity(meth_bed, insertion_bed, window_size, max_distance, outfile):
    # Read the genomic coordinates BED file
    insertion_bed = pd.read_csv(insertion_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "seq"])

    # Read the target BED file
    meth_bed = pd.read_csv(meth_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "mod"])
    #max_distance = int(max_distance)
    #window_size = int(window_size)
    # Iterate through each genomic coordinate
    for index, row in insertion_bed.iterrows():
        start = int(row[1])
        end = int(row[2])
        print(start)
        i=0
        modifications={}
        while i < max_distance +1:
            start_interval = start + i
            count=0
            for index,entry in meth_bed.iterrows(): #if entry on chr Z is within interval, modification counter +1 #pack into function later that returns number of mods
                if row['Chr'] == entry['Chr'] and \
                   start_interval <= entry['start'] and \
                   start + max_distance >= entry['end']:
                    count += 1
            modifications[start_interval] = count
            i = i + window_size
            #interval_entries = meth_bed[start_interval : start_interval + window_size] #all entries in certain interval
            #
            #end_interval = end - i 
            # Calculate the mean of the entries' values (assuming a single-column BED file)
            #mean_value = sum(int(entry.fields[3]) for entry in interval_entries) / len(interval_entries)
            print(modifications)
'''