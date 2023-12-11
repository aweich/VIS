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
from pybedtools import BedTool
   
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

def fragmentation_match_distribution(data, fragment_specifier, outpath):
    """
    Takes the vector fragments and plots histogram of their frequency in the alignment
    """
    blasted = pd.read_csv(data, sep="\t")
    if blasted['QueryID'][0].startswith("V"):
        blasted[['Vector', 'Fragment']] = blasted['QueryID'].str.split('_', n=1, expand=True)
        blasted["Fragment"] = pd.to_numeric(blasted["Fragment"])
        freq = collections.Counter(blasted["Fragment"].sort_values())
    else:
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
    if blasted["QueryID"][0].startswith("V"):
        freq = collections.Counter(blasted["SubjectID"])
    else:
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
    Adds the fasta sequence to the respective bed gaps. This makes the counting of Cytosines easier for a later step.
    """
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, usecols=[0,1,2], names=['Chromosome', 'Start', 'End'])
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Add a new column to the DataFrame with the corresponding sequence from the FASTA file
    bed_df['Sequence'] = [fasta_sequences[i].seq for i,n in enumerate(fasta_sequences)]
    
    bed_df.to_csv(output_bed_path, sep='\t', header=False, index=False)

####all below are just for the mean methylation in the equal sized bin in proximity of the read with the insertion
def bed_intersect_count(chromosome, start, stop, meth_bed):
    """
    Python wrapper of bed intersect function that returns the number of intersections. 
    Takes in a single chr start stop combination and a whole file of possible intersections. 
    """
    # Create a temporary BED file from the list of coordinates
    temp_bed = BedTool(f"{chromosome}\t{start}\t{stop}", from_string=True)
    # Perform the intersect operation and count the number of overlaps (-c option)
    intersected = temp_bed.intersect(meth_bed, c=True)
    # Calculate and return the number of overlaps
    num_overlaps = intersected[0][3]
    return int(num_overlaps)
'''
def count_equal_entries(bed_file_path):
    """
    Counts how often a chr start stop entry occurs in a bed file. 
    """
    # Read the BED file using pybedtools
    bed = BedTool(bed_file_path)

    # Initialize a dictionary to store counts
    entry_counts = {}

    # Count occurrences of each entry using only the first three columns
    for entry in bed:
        entry_str = '\t'.join(map(str, entry.fields[:3]))
        entry_counts[entry_str] = entry_counts.get(entry_str, 0) + 1

    return entry_counts
'''
def collapse_equal_entries(bed_file):
    """
    Collapses duplicates chr start stop entries into a single representative.
    """
    bed = BedTool(bed_file)
    # Initialize a set to keep track of unique entries
    unique_entries = set()
    # Create a list to store representative entries
    representative_entries = []
    # Iterate through each entry
    for entry in bed:
        entry_str = '\t'.join(map(str, entry.fields[:3]))

        # Check if the entry is unique
        if entry_str not in unique_entries:
            unique_entries.add(entry_str)
            representative_entries.append(entry)

    # Create a new BedTool with representative entries
    collapsed_entries = BedTool(representative_entries)

    return collapsed_entries


def C_in_range(fasta, start, stop):
    """
    Uses start and stop coordinates and checks how many Cs are contained within start and stop. 
    """
    count=0
    subfasta = fasta[start:stop]
    count = subfasta.lower().count('c')
    if count != 0:
        return count
    else:
        return 1

def methylation_in_insertion_proximity(meth_bed, insertion_bed, window_size, max_distance, outfile):
    """
    For each entry in insertion bed, a interval of window_size will be scanned for occurring Cytosines (=N total Cs) up to max_distance in positive and negative direction.
    For each of these intervals, the meth_bed file is scanned for entries contained within the genomic ranges (=N modified Cs). 
    The modification ratio in % for each interval is calculated and a dataframe is returned.
    """
    # Read the genomic coordinates BED file
    insertion_bed = pd.read_csv(insertion_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "seq"])
    print("This will take a while...")

    # Read the target BED file
    #meth_bed = pd.read_csv(meth_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "mod"])
    with_duplicates = BedTool(meth_bed) #if I at some point figure out why one entry can be in this stupid file up to twelve times, I might change this
    meth_bed = collapse_equal_entries(with_duplicates)
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
        modifications = bed_intersect_count(chromosome, start, end, meth_bed)
        bases = C_in_range(row.iloc[3], max_distance, len(row.iloc[3]) - max_distance) #the proximity bed file with fasta has the structure: max_dist-start-stop-max_dist
        final_bed.loc[index, str("Insertion_point")] = (modifications/bases) *100 #puts key (0)-value(methylation-ratio) pair into dataframe in the respective line (index)
        #methylation in 3' direction in window-size steps up to max_distance
        i=0
        while i < max_distance+1: #3' direction
            start_interval = end + i
            stop_interval = start_interval + window_size
            modifications = bed_intersect_count(chromosome, start_interval, stop_interval, meth_bed)
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
            modifications = bed_intersect_count(chromosome, stop_interval, start_interval, meth_bed)
            column_name=str(i-window_size)
            bases = C_in_range(row.iloc[3], max_distance-i-window_size, max_distance-i)
            print(i, chromosome, orig_insertion, modifications, bases)
            final_bed.loc[index, column_name] = (modifications/bases) *100
            i = i - window_size 
            
    print(list(final_bed.columns))
    #print(np.array(final_bed.iloc[:,5:47]).sum(axis=0))
    #result = count_equal_entries(meth_bed)
    #for entry, count in result.items():
    #    print(f"Entry: {entry}, Count: {count}")

    final_bed.to_csv(outfile, sep='\t', header=True, index=False)

def flatten_comprehension(matrix):
    """
    Flattens a nested list and makes each element a str.
    """
    return [str(item) for row in matrix for item in row]

def plot_modification_proximity(bedfile, outfile):
    """
    Creates heatmap of interval and methylation mean based on the data created in 'methylation_in_insertion_proximity()'
    """
    # Process each BED file
    df = pd.read_csv(bedfile, sep='\t')
    df["Read with insertion"] = df["Chr"] + "_" + df["Insertion"]
    df = df.set_index('Read with insertion')
    df = df.drop(columns=["Chr","start","end","seq","Insertion"])
    
    #for x axis
    pos=list(range(0, 10001, 500))
    pos = [f'+{num}' if num > 0 else num for num in pos]
    column_order = flatten_comprehension([list(range(-10000, 0, 500)), ["Insertion_point"], pos])
    column_order.remove('0')
    # Reorder the DataFrame columns
    df_reordered = df[column_order]

    # Create a Seaborn heatmap
    plt.figure(figsize=(16, 9))
    sns.clustermap(df_reordered, cmap="YlGnBu", annot=False, row_cluster=True, linewidths=.5, col_cluster=False,\
                   cbar_pos=(0.25, .9, .5, 0.01), cbar_kws={'orientation': 'horizontal'},\
                   dendrogram_ratio=(.15))
    
    #plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45) # For x axis #add g =sns...
    plt.xlabel("%C with modification")
    plt.ylabel('')
    plt.title('')
    plt.savefig(outfile, bbox_inches="tight")

#sniffles
def preprocess_sniffles(vcf_path, outfile, outfile2, outfile3):
    """
    Reads sniffles output and re-formats it into something more useful. 
    Deletion gets also a start coordinate although this does not make any biological sense! This is just for the plotting later!
    """
    # Read VCF file into a DataFrame
    vcf = pd.read_csv(vcf_path, comment='#', header=None, sep='\t', dtype=str)
    # Extract relevant columns
    relevant_columns = pd.DataFrame(vcf.iloc[:, [0, 1, 2]])

    # Extract SVType from V3
    relevant_columns['SVType'] = vcf.iloc[:, 2].str.split('.').str[1]
    vcf["Start"] = vcf.iloc[:, 7].str.extract(r'SVLEN=(\d+)')
    vcf['Start'] = vcf['Start'].fillna(0)
    # Calculate start coordinates (subtract SVLEN from V2)
    relevant_columns['Start'] = vcf.iloc[:, 1].astype(int) - vcf["Start"].astype(int)
    # Rename and reshape columns
    relevant_columns = relevant_columns.drop(2, axis=1)
    relevant_columns = relevant_columns.rename(columns={0: 'Chromosome', 1: 'Stop'})
    #relevant_columns.columns = ['Chromosome', 'Stop', 'SVType', 'Start']
    relevant_columns = relevant_columns[['Chromosome', 'Start', 'Stop', 'SVType']]
    # Filter out chromosomes
    relevant_columns = relevant_columns[relevant_columns['Chromosome'].str.match(r'^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$')]
    print(relevant_columns["Chromosome"].value_counts())
    relevant_columns.to_csv(outfile, sep='\t', header=True, index=False)
    #split into INS and DEL
    mask = relevant_columns.SVType.str.contains("INS")
    INS = relevant_columns[mask]
    OTHER = relevant_columns[~mask]
    INS.to_csv(outfile2, sep='\t', header=False, index=False)
    OTHER.to_csv(outfile3, sep='\t', header=False, index=False)

def blastn_bed_merger(blast, bed, outfile):
    """
    Merge script for blastn output and bamtobed output to further narrow insertion from read level to coordinate level.
    """
    blast = pd.read_csv(blast, sep='\t')
    bed = pd.read_csv(bed, sep='\t', header=None, usecols=[0,1,2,3])
    bed["ReadSize"] = bed[2] - bed[1]
    bed = bed.drop(bed[bed.ReadSize < 80].index) #gets rid of matches that can't be the insertion site (since they are shorted than the smallest fragment (80bp/VO99))
    # Merge based on read names
    #merged_bed = pd.merge(blast, bed, left_on=["QueryID", "QueryLength"], right_on=[3, "ReadSize"])
    merged_bed = pd.merge(blast, bed, left_on=["QueryID"], right_on=[3])
    #reshape
    #start = merged_bed[1] +
    #stop = merged_bed[2] - 
    #insertions = pd.DataFrame({"Chromosome": merged_bed[0], "Start": merged_bed[1] , "Stop": merged_bed[2], "Fragment": merged_bed["SubjectID"], "Read": merged_bed["QueryID"]})
    #merged_bed = merged_bed[['Chromosome', 'Start', 'Stop', 'SVType']]
    print(merged_bed.head())
    # Save the merged BED file
    merged_bed.to_csv(outfile, sep='\t', index=False)
