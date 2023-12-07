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
    Adds the fasta sequence to the respective bed gaps. This makes the following methylation counting step easier.
    """
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, usecols=[0,1,2], names=['Chromosome', 'Start', 'End'])
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Add a new column to the DataFrame with the corresponding sequence from the FASTA file
    bed_df['Sequence'] = [fasta_sequences[i].seq for i,n in enumerate(fasta_sequences)]
    
    bed_df.to_csv(output_bed_path, sep='\t', header=False, index=False)

def bed_in_bed(chromosome, start, end, bed_data):
    """
    Uses start and stop coordinates and checks how many entries of bed_data are contained within start and stop. 
    """
    count=0
    for index,entry in bed_data.iterrows(): #if entry on chr Z is within interval, modification counter +1 #pack into function later that returns number of mods
        if chromosome == entry['Chr'] and start <= entry['start'] and end >= entry['end']:
                    count += 1
    return count
 
def C_in_range(fasta, start, stop):
    """
    ´Uses start and stop coordinates and checks how many Cs are contained within start and stop. 10k from each end in 500 steps. Only read have more than 500 Cs
    """
    count=0
    subfasta = fasta[start:stop]
    count = subfasta.lower().count('c')
    if count != 0:
        return count
    else:
        return 1
def methylation_in_insertion_proximity(meth_bed, insertion_bed, window_size, max_distance, outfile):
    # Read the genomic coordinates BED file
    insertion_bed = pd.read_csv(insertion_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "seq"])
    print("This will take a while...")

    # Read the target BED file
    meth_bed = pd.read_csv(meth_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "mod"])

    #output BED
    final_bed = insertion_bed.copy()
    # Iterate through each genomic coordinate
    for index, row in insertion_bed.iterrows():
        chromosome = row.iloc[0]
        start = int(row.iloc[1]) + max_distance #to end up with the original coordinates again
        end = int(row.iloc[2]) - max_distance
        intervalsize = end - start
        orig_insertion= str(start)+"_"+str(end)
        final_bed.loc[index, "Insertion"] = str(orig_insertion)
        #methylation at insertion itself: Good overview of basic principle
        modifications = bed_in_bed(chromosome, start, end, meth_bed) #counts number of entries in modification bed in range start-end
        bases = C_in_range(row.iloc[3], max_distance, len(row.iloc[3]) - max_distance) #the proximity bed file with fasta has the structure: max_dist-start-stop-max_dist
        final_bed.loc[index, str("Insertion_point")] = modifications/(bases) *100 #puts key (0)-value(methylation-ratio) pair into dataframe in the respective line (index)
        #methylation in 3' direction in window-size steps up to max_distance
        i=0
        while i < max_distance: #3' direction
            start_interval = end + i
            stop_interval = start_interval + window_size
            modifications = bed_in_bed(chromosome, start_interval, stop_interval, meth_bed) #still bad quality code, but better
            column_name="+" + str(i+window_size)
            bases = C_in_range(row.iloc[3], max_distance+intervalsize+i, max_distance+intervalsize+i+window_size)
            print(i, chromosome, orig_insertion, modifications, bases)
            final_bed.loc[index, column_name] = (modifications/bases) *100
            i = i + window_size
        #methylation in 5' direction in window-size steps up to max_distance
        i= 0
        while abs(i) < max_distance: #5' direction
            start_interval = start - i
            stop_interval = start_interval - window_size
            modifications = bed_in_bed(chromosome, stop_interval, start_interval, meth_bed) #still bad quality code, but better; start and stop exchanged!
            column_name=str(i-window_size)
            bases = C_in_range(row.iloc[3], max_distance-i-window_size, max_distance-i)
            print(i, chromosome, orig_insertion, modifications, bases)
            final_bed.loc[index, column_name] = (modifications/bases) *100
            i = i - window_size 
            
    print(list(final_bed.columns))
    print(np.array(final_bed.iloc[:,5:47]).sum(axis=0))
    final_bed.to_csv(outfile, sep='\t', header=True, index=False)

def plot_modification_proximity(bedfile, outfile):
    """
    Creates heatmap
    """
    # Process each BED file
    df = pd.read_csv(bedfile, sep='\t')
    df["Concatenated"] = df["Chr"] + "_" + df["Insertion"]
    df = df.set_index('Concatenated')
    df = df.drop(columns=["Chr","start","end","seq","Insertion"])
    print(df.head())

    # Create a Seaborn heatmap
    plt.figure(figsize=(16, 9))
    sns.clustermap(df, cmap="YlGnBu", annot=True, cbar_pos=(0, .2, .03, .4))
    plt.xlabel('Distance from Insertion site')
    plt.ylabel('Read with insertion')
    plt.title('')
    plt.savefig(outfile, bbox_inches="tight")