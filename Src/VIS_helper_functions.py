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
import collections
import seaborn as sns
from pybedtools import BedTool
import json
import subprocess
from Bio.Align.Applications import ClustalOmegaCommandline
   
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

def insertion_normalisation(insertions, n, fasta, outpath):
    """
    Uses the number of aligned bases, number of insertions, and scaling factor n for normalization.
    Results in insertions per n aligned bases. 
    E.g. A Normalised count of 3 with the scaling factor 10000 means that there were on average 3 insertions detected for each 1kb aligned bases. 
    """
    N50 = get_N50(fasta)
    
    n_insertions = sum(1 for _ in open(insertions))  # bed entries are by definition 1 per line
    print(n_insertions)
    
    # Calculate the total number of bases in the FASTA file
    n_bases = 0
    read_lengths = []
    with open(fasta, 'r') as fasta_file:
        for line in fasta_file:
            if not line.startswith('>'):
                n_bases += len(line.strip())
                read_lengths.append(len(line.strip()))

    mean_read_length = sum(read_lengths) / len(read_lengths) if read_lengths else 0

    
    if n_bases != 0:
        normalised = (n / n_bases) * n_insertions
        with open(outpath, 'w') as output_file:
            output_file.write("N50: " + str(N50) + "\n")
            output_file.write("Mean Read Length: " + str(mean_read_length) + "\n")
            output_file.write("Number of Bases: " + str(n_bases) + "\n") 
            output_file.write("Number of Insertions: " + str(n_insertions) + "\n")
            output_file.write("Normalisation factor: " + str(n) + "\n")
            output_file.write("Insertions per " + str(n) + " Bases: " + str(normalised))
    else:
        with open(outpath, 'w') as output_file:
            output_file.write("N50: " + str(N50) + "\n")
            output_file.write("Mean Read Length: 0\n")
            output_file.write("Number of Bases: " + str(n_bases) + "\n") 
            output_file.write("Number of Insertions: " + str(n_insertions) + "\n")
            output_file.write("Normalisation factor: " + str(n) + "\n")
            output_file.write("No reads found. Check the input!")
        

'''                    
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
'''
def get_read_lengths_from_fasta(fasta_file):
    """
    returns fasta read lengths
    """
    read_lengths = []
    with open(fasta_file, 'r') as file:
        sequence_length = 0
        for line in file:
            if line.startswith('>'):
                if sequence_length > 0:
                    read_lengths.append(sequence_length)
                sequence_length = 0
            else:
                sequence_length += len(line.strip())

        if sequence_length > 0:
            read_lengths.append(sequence_length)
    return read_lengths

def get_N50(fasta):
    """Calculate read length N50.
    Based on https://github.com/PapenfussLab/Mungo/blob/master/bin/fasta_stats.py
    """
    lengths = get_read_lengths_from_fasta(fasta)
    return lengths[np.where(np.cumsum(lengths) >= 0.5 * np.sum(lengths))[0][0]]


def fragmentation_match_distribution(data, fragment_specifier, outpath):
    """
    Takes the vector fragments and plots histogram of their frequency in the alignment
    """
    blasted = pd.read_csv(data, sep="\t")

    # Check if the data is empty
    if blasted.empty:
        plt.figure()
        plt.title('Empty Data')
        plt.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=12)
        plt.xticks([])
        plt.yticks([])
        outfile = outpath + str("/") + f'{fragment_specifier}_fragmentation_distribution.png'
        plt.savefig(outfile)
        plt.close()
        return
    
    if any(x.isupper() for x in blasted['QueryID'][0]): #to make sure the right column is used for plotting. Reads do not have any uppercase letters
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
    
    #check if empty
    if blasted.empty:
        plt.figure()
        plt.title('Empty Data')
        plt.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=12)
        plt.xticks([])
        plt.yticks([])
        outfile = outpath + str("/") + f'{fragment_specifier}_fragmentation_distribution.png'
        plt.savefig(outfile)
        plt.close()
        return

    if any(x.isupper() for x in blasted['QueryID'][0]):
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
    if len(counts_df.columns) > 1:
        sns.clustermap(counts_df, cmap="YlGnBu", annot=True, cbar_pos=(0, .2, .03, .4))
        plt.xlabel('File ID')
        plt.ylabel('Chromosome')
        plt.title('Chromosome Occurrences \n  in BED Files')
        plt.savefig(outfile, bbox_inches="tight")
    else:    
        sns.heatmap(counts_df, annot=True,cmap="YlGnBu")
        plt.xlabel('File ID')
        plt.ylabel('Chromosome')
        plt.title('Chromosome Occurrences \n  in BED Files')
        plt.savefig(outfile, bbox_inches="tight")

def insertion_proximity(bedfile, binsize, outfile):
    """
    Takes a list of insertion sites (BED file) and adds a user-defined number of bp to each site (start and end). Returns new BED file.
    """
    df = pd.read_csv(bedfile, sep='\t', lineterminator='\n', usecols=[0,1,2,3,4,5], names = ["Chr","start","end", "read","CleavageSites","strand"])
    df_modified = df.copy()
    df_modified['start'] = df['start'] - binsize
    df_modified['end'] = df['end'] + binsize
    df_modified.to_csv(outfile, sep = '\t', header = False, index = False)

def add_sequence_column(bed_file_path, fasta_file_path, output_bed_path):
    """
    Adds the fasta sequence to the respective bed gaps. This makes the counting of Cytosines easier for a later step.
    """
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, usecols=[0,1,2,3,4,5], names=['Chromosome', 'Start', 'End','Read',"CleavageSites",'strand'])
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))
    # Add a new column to the DataFrame with the corresponding sequence from the FASTA file
    bed_df['Sequence'] = [fasta_sequences[i].seq for i,n in enumerate(fasta_sequences)]
    
    bed_df.to_csv(output_bed_path, sep='\t', header=False, index=False)

def add_insertion_sequence(bed_file_path, fasta_file_path, output_bed_path):
    """
    Adds the fasta sequence to the insertion itself. This makes the counting of Cytosines easier for a later step.
    """
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, usecols=[0,1,2,3,4,5], names=['Chromosome', 'Start', 'End', 'Read','CleavageSites','strand'])
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))
    # Add a new column to the DataFrame with the corresponding sequence from the FASTA file
    seqs=[]

    for entry in bed_df["Read"]: #part before _Insertion
        entry = entry.split("_")[0]
        for i,n in enumerate(fasta_sequences):  
            if fasta_sequences[i].id.split("_")[0] == entry: #part before _0/_1 etc.     
                seqs.append(fasta_sequences[i].seq)
    
    bed_df["Sequence"] = seqs
    bed_df = bed_df[['Chromosome', 'Start', 'End', 'Sequence', 'Read','CleavageSites','strand']]
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
        return 1 #otherwise division 0/0 -> breaks

def methylation_in_insertion_proximity(meth_bed, insertion_bed, window_size, max_distance, outfile):
    """
    For each entry in insertion bed, a interval of window_size will be scanned for occurring Cytosines (=N total Cs) up to max_distance in positive and negative direction.
    For each of these intervals, the meth_bed file is scanned for entries contained within the genomic ranges (=N modified Cs). 
    The modification ratio in % for each interval is calculated and a dataframe is returned.
    window_size !=0  and max_distance =0: Modification ratio over FASTA and BED cooridnates
    wind_size !=0 and max_distance != 0: Modification ratio over FASTA and proximity/max_distance +/- BED coordinates.

    """
    # Read the genomic coordinates BED file
    insertion_bed = pd.read_csv(insertion_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3,4,5], names=['Chromosome', 'Start', 'End', 'Read','CleavageSites','strand'])
    print("This will take a while...")

    # Read the target BED file
    #meth_bed = pd.read_csv(meth_bed, sep='\t', lineterminator='\n', usecols=[0,1,2,3], names = ["Chr","start","end", "mod"])
    with_duplicates = BedTool(meth_bed) #if I at some point figure out why one entry can be in this stupid file up to twelve times, I might change this
    meth_bed = collapse_equal_entries(with_duplicates)
    #output BED
    final_bed = insertion_bed.copy()
    print(insertion_bed.head())
    #percent_dict={} #collection for the insertions
    # Iterate through each genomic coordinate
    for index, row in insertion_bed.iterrows():
        chromosome = row.iloc[0]
        start = int(row.iloc[1]) + max_distance #to end up with the original coordinates again; only relevant in case of buffer
        end = int(row.iloc[2]) - max_distance
        intervalsize = end - start #FASTA length
        orig_insertion= str(start)+"_"+str(end)
        final_bed.loc[index, "Insertion"] = str(orig_insertion)
        
        #methylation on insertion read itself: Good overview of basic principle
        #modifications = bed_intersect_count(chromosome, start, end, meth_bed)
        #bases = C_in_range(row.iloc[3], max_distance, len(row.iloc[3])+1 - max_distance) #the proximity bed file with fasta has the structure: max_dist-start-stop-max_dist/ for max_distance =0, the whole FASTA is used
        #final_bed.loc[index, str("Mean_Mod_Percentage")] = (modifications/bases) *100 #puts key (0)-value(methylation-ratio) pair into dataframe in the respective line (index)
        #methylation in 3' direction in window-size steps up to max_distance
    
        #settings for the proximity option (= not rad-level only)
        if max_distance != 0:
            i=0
            while i < max_distance+1: #3' direction
                start_interval = end + i
                stop_interval = start_interval + window_size
                modifications = bed_intersect_count(chromosome, start_interval, stop_interval, meth_bed)
                column_name="+" + str(i+window_size)
                bases = C_in_range(row.iloc[3], max_distance+intervalsize+i, max_distance+intervalsize+i+window_size)
                #print(i, chromosome, orig_insertion, modifications, bases)
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
                #print(i, chromosome, orig_insertion, modifications, bases)
                final_bed.loc[index, column_name] = (modifications/bases) *100
                i = i - window_size 
        
            final_bed.to_csv(outfile, sep='\t', header=True, index=False)    
        
    #settings for the insertion itself:
    #this is - and + dependent! if the FASTA mapped on the - strand, then the C numbers must be caluculated from FASTA back to front 
        for n,i in enumerate(range(start, end, window_size)):
            
            if row[4] == '-':
                reversed_seq = row.iloc[3][::-1]    
                bases = C_in_range(reversed_seq, n*window_size, (n+1)*window_size)
                modifications = bed_intersect_count(chromosome, i, i+window_size, meth_bed)
                final_bed.loc[index, (n+1)*window_size] = (modifications/bases) *100
            
            bases = C_in_range(row.iloc[3], n*window_size, (n+1)*window_size)
            modifications = bed_intersect_count(chromosome, i, i+window_size, meth_bed)
            final_bed.loc[index, (n+1)*window_size] = (modifications/bases) *100
    
    #final_bed['Modification'] = final_bed['seq'].map(percent_dict)
    #print(percent_dict)
    final_bed.to_csv(outfile, sep='\t', header=True, index=False)

def flatten_comprehension(matrix):
    """
    Flattens a nested list and makes each element a str.
    """
    return [str(item) for row in matrix for item in row]

def plot_modification_proximity(bedfile,window_size, max_distance, outfile): #needs arguments for the range and step-size
    """
    Creates heatmap of interval and methylation mean based on the data created in 'methylation_in_insertion_proximity()'
    """ 
    if "ID" in pd.read_csv(bedfile, sep='\t'):
        print("Can use existing IDs...")
        df = pd.read_csv(bedfile, sep='\t')
    elif type(bedfile) is not str: #checks if we have snakemake list input or regular str input (i.e. multiple files (snakemake list) or just one)
        bedlist=[]
        print("Create IDs...")
        for filename in bedfile:
            df = pd.read_csv(filename, sep='\t')
            df["ID"] = str(filename).split("/")[-1]
            bedlist.append(df)
        df = pd.concat(bedlist, axis=0, ignore_index=True)
    else:
        print("IDs from filename...")
        df = pd.read_csv(bedfile, sep='\t')
        df["ID"] = str(bedfile).split("/")[-1]
    
    df["Read with insertion"] = df["Chr"] + "_" + df["Insertion"]
    df = df.set_index('Read with insertion')
    df = df.drop(columns=["Chr","start","end","seq","Insertion"])
    df['Mean_Mod_Percentage'] = 100 #replace insertion site values with max to create a border
    
    #ID column
    ID=df["ID"]
    palette = sns.color_palette("hls", len(ID.unique()))
    lut = dict(zip(ID.unique(), palette)) #needs to be adjusted for number of samples of course #'rgb'
    row_color = ID.map(lut)
    #for x axis
    pos=list(range(0, max_distance+1, window_size))
    pos = [f'+{num}' if num > 0 else num for num in pos]
    column_order = flatten_comprehension([list(range(-max_distance, 0, window_size)),["Mean_Mod_Percentage"], pos]) #["Mean_Mod_Percentage"], before pos
    column_order.remove('0')
    # Reorder the DataFrame columns
    df_reordered = df[column_order]
    # Create a Seaborn heatmap
    plt.figure(figsize=(16, 9))
    sns.clustermap(df_reordered, row_colors=row_color,yticklabels=True, cmap="YlGnBu", annot=False, row_cluster=True, linewidths=.75, col_cluster=False,\
                   cbar_pos=(0.25, .9, .5, 0.01), cbar_kws={'orientation': 'horizontal'},\
                   dendrogram_ratio=(.15))
    
    #legend for ID column
    for label in df["ID"].unique():
        plt.bar(0, 0, color=lut[label], label=label, linewidth=0)
    plt.legend(loc='lower center', bbox_to_anchor=[0.5, 0.75])

    #plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45) # For x axis #add g =sns...
    plt.xlabel("%C with modification")
    plt.ylabel('')
    plt.title('')
    plt.savefig(outfile, bbox_inches="tight")



####this part here is dedicated to the splitting of blast-match including fasta reads
def merge_intervals(intervals, overlap, filtering, filtervalue):
    # Sort intervals by start coordinates
    #sorted_intervals = sorted(intervals, key=lambda x: x[0])
    intervals.sort() #sorts all intervals in ascending order: Overlaps are possible! #sorts inplace
    print(intervals)

    merged_intervals = []
    #the following part does not work yet: It needs to take our different list properties into account!
    
    #first coordinates
    merged_intervals.append(intervals[0]) #start that needs to be there always
    
    #intermediate coordinates if needed
    for i in range(0, len(intervals)-1): #iterating over start and stops with 1 offset! this means that we take a look at stop_interval1 and start_interval2
        #print(str(intervals[i]),str(intervals[i + 1]))
        if abs(intervals[i] - intervals[i + 1]) >= overlap: #strictness filter
            merged_intervals.append(intervals[i]) #adds intermediate coordinate as end of interval
            merged_intervals.append(intervals[i+1]) # = new start
    
    #last element
    merged_intervals.append(intervals[-1]) #last element    

    print(merged_intervals)
    
    #iterate throgh list and substract current from previous for the filter
    if filtering:
        for i in [y - x for x,y in zip(merged_intervals,merged_intervals[1:])]:
            if abs(i) > filtervalue:
                return merged_intervals
            else:
                return None
    #if no filter, just return all of them
    print("no filter...")
    
    return merged_intervals

   
def splitting_borders(blast_file, filteroption, filtervalue,  overlap, outfile1, outfile2):
    """
    This function takes in a BLASTn output file and returns a file with the reads and the respective borders for fasta splitting.
    Output format: Read_ID cut1, cut2, cutN 
    """
    # Read the BLAST output into a DataFrame
    blast_data = pd.read_csv(blast_file, sep='\t')

    # Group the data by QueryID
    grouped_data = blast_data.groupby('QueryID')

    # Dictionary to store intervals for each QueryID
    intervals_dict = {}

    # Iterate over groups and extract intervals
    for query_id, group in grouped_data:
        intervals = []
        for _, row in group.iterrows():
            start, end = row['QueryStart'], row['QueryEnd']
            intervals.append(start)
            intervals.append(end)
        # Store the intervals for each QueryID
        intervals_dict[query_id] = intervals

    #output if something is found
    result_dict = {key: merge_intervals(intervals, overlap, filteroption,filtervalue) for key, intervals in intervals_dict.items()}
    result_dict = {k: v for k, v in result_dict.items() if v is not None}
    json.dump(result_dict, open(outfile1,'w'))
    with open(outfile2, 'w') as output_file:
        for key in result_dict.keys():
            output_file.write(f"{key}\n")

def split_fasta_with_breakpoints(fasta_string, breakpoints):
    """
    Takes fasta string and corresponding list of breakpoints as arguments and returns list of sequences :breakpoint1,breakpoint2:breakpoint3, ... breakpointn:
    """
    sequences = []
    
    # Add the sequence before the first breakpoint
    sequences.append(fasta_string[:breakpoints[0]])
    #print(fasta_string)
    print("Breakpoints: " + str(breakpoints))

    # Iterate over pairs of breakpoints
    for i in range(len(breakpoints) - 1):
        start = breakpoints[i]
        end = breakpoints[i + 1]
        sequences.append(fasta_string[start:end])

    # Add the sequence after the last breakpoint
    sequences.append(fasta_string[breakpoints[-1]:])
    #drop every second interval if there are more than two snippets of the fasta: -> To save only the intervals without vector
    if len(sequences) > 2:
        del sequences[1::2]  
    
    #extraction of the vector sequence
    vector=[]
    
    for i in range(len(breakpoints)):
        if i % 2 == 0 or i == 0:
            vector.append(fasta_string[breakpoints[i]:breakpoints[i+1]])
    
    return sequences, vector


def split_fasta_by_borders(border_dict, fasta, mode, outfasta, outvector):
    """
    Uses previously created breakpoints in border dict to cut out insertions from fasta file
    """
    border_dict = json.load(open(border_dict)) 
    with open(outfasta, 'w') as output_file:
        with open(outvector, 'w') as output_vector_file:
            # opening given fasta file using the file path
            with open(fasta, 'r') as fasta_file:
                # extracting multiple data in single fasta file using biopython
                for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
                    #split FASTA sequence into equally sized sub-lists with different id
                    record_id = record.id
                    record_seq = record.seq
                    if record_id in border_dict:
                        record_list, vector_list = split_fasta_with_breakpoints(record_seq, border_dict[record_id])
                        # write out inserted vector sequence
                        for i,entry in enumerate(vector_list):
                            vector_seq_new = vector_list[i]
                            vector_id_new = str(record_id) + '_%i' % (i)
                            print("Split " + str(record_id) + "into " + str(vector_id_new))
                            output_vector_file.write(">"+str(vector_id_new)+"\n"+str(vector_seq_new) + "\n")
                        # write out fasta without insertion        
                        if mode == "Separated":
                            #this part only if the non-insertion fragments should NOT be combined
                            for i,entry in enumerate(record_list):
                                record_seq_new = record_list[i]
                                record_id_new = str(record_id) + '_%i' % (i)
                                print("Split " + str(record_id) + "into " + str(record_id_new))
                                output_file.write(">"+str(record_id_new)+"\n"+str(record_seq_new) + "\n")
                        elif mode == "Buffer":
                            n_buffer = 100
                            buffer = n_buffer*'N'
                            record_seq_new = buffer.join([str(n) for n in record_list])
                            print(record_seq_new)
                            record_id_new = str(record_id) + '_Buffer'+str(n_buffer)+'Insertion'
                            output_file.write(">"+str(record_id_new)+"\n"+str(record_seq_new) + "\n")
                        else: 
                            record_seq_new = ''.join([str(n) for n in record_list])
                            record_id_new = str(record_id) + '_Insertion'
                            output_file.write(">"+str(record_id_new)+"\n"+str(record_seq_new) + "\n")
                    else:
                        output_file.write(">"+str(record_id)+"\n"+str(record_seq) + "\n")     


def adjust_coordinates_for_strand(start, stop, strand, insertion_values):
    """
    Adjusts the insertion coordinates based on strand directionality.
    
    Parameters:
        start (int): Start coordinate of the read from the BED file.
        stop (int): Stop coordinate of the read from the BED file.
        strand (str): Strand information ('+' or '-').
        insertion_values (list): List of insertion coordinates from the border_dict.
    
    Returns:
        list: Adjusted coordinates as a list of (new_start, new_stop) tuples.
    """
    read_length = stop - start
    adjusted_coords = []

    if strand == "+":
        # Use insertion values as-is for the + strand
        for i in range(0, len(insertion_values), 2):
            new_start = start + insertion_values[i]
            new_stop = start + insertion_values[i + 1]
            adjusted_coords.append((new_start, new_stop))
    elif strand == "-":
        # Reverse and flip insertion values for the - strand
        for i in range(0, len(insertion_values), 2):
            new_start = start + (read_length - insertion_values[i + 1])
            new_stop = start + (read_length - insertion_values[i])
            adjusted_coords.append((new_start, new_stop))

    return adjusted_coords

def exact_insertion_coordinates(border_dict, bed, outfile):
    """
    Updates BED file coordinates based on border_dict values.
    Adjusts for strand directionality and creates new BED entries for each insertion.
    """
    bed = pd.read_csv(bed, sep='\t', header=None, usecols=[0, 1, 2, 3, 5])
    bed['BaseRead'] = bed[3].str.split("_").str[0]  # Get the original read name

    # Dict with nsertion intervals
    border_dict = json.load(open(border_dict))

    # Filter BED entries to only include those present in the border_dict
    matching_entries = bed[bed["BaseRead"].isin(border_dict.keys())].copy()

    # Check if there are matching entries
    if matching_entries.empty:
        print("No matching reads found in the BED file.")
        return

    # Process matching entries by updating cooridnates
    results = []
    for _, row in matching_entries.iterrows():
        read_name = row["BaseRead"]
        old_start = row[1]
        old_stop = row[2]
        strand = row[5]
        chromosome = row[0]

        # Get insertion coordinates from the dictionary
        insertion_values = border_dict[read_name]

        # Ensure insertion_values are in pairs (start, stop)
        if len(insertion_values) % 2 != 0:
            print(f"Error: Uneven number of coordinates for read {read_name}")
            continue

        # Adjust coordinates based on strand
        adjusted_coords = adjust_coordinates_for_strand(old_start, old_stop, strand, insertion_values)
        
        old_coords = [old_start, old_stop]
        for new_start, new_stop in adjusted_coords:
            results.append([chromosome, new_start, new_stop, read_name, old_coords, strand])

    updated_bed = pd.DataFrame(results, columns=[0, 1, 2, 3, 5, 6])
    updated_bed.to_csv(outfile, sep='\t', index=False, header=False)
    print(f"Insertion BED file saved to {outfile}")


def combine_beds_add_ID(beds, outfile): #maybe at some point add to the heatmap plotting function
    """
    Combines BEDs and gives them an ID column
    """
    bedlist=[]
    for filename in beds:
        df = pd.read_csv(filename, sep='\t', header=None)
        df[len(df.columns)] = str(filename).split("/")[-1] #add column to dataframe without header
        bedlist.append(df)
    df = pd.concat(bedlist, ignore_index=True)
    df.to_csv(outfile, sep='\t', index=False, header=False)

def add_annotation_column_bed(bed1, bed2, outfile):
    """
    Adds column of one bed to another bed if coordinates match
    """
    bed1 = pd.read_csv(bed1, sep='\t')
    bed2 = pd.read_csv(bed2, sep='\t', usecols=[0,1,2,3,4], names = ["Chr","start","end","seq","ID"])
    bed1 = bed1.merge(bed2, on=["Chr","start","end","seq"])
    bed1.to_csv(outfile, sep='\t', index=False, header=True)

def plot_modification_per_vectorlength(meanmodbed, window_size, outpath):
    """
    PLots how many bases are modified for the length of the sequence
    """
    print(outpath)
    mod = pd.read_csv(meanmodbed, sep='\t')
    mod["ID"] = mod["Chr"] + "_" + mod["Insertion"]
    mod = mod.drop(columns=['Chr', 'start','end','seq', 'Insertion','strand'])
    mod = pd.melt(mod, id_vars=['ID'])
    for i in mod["ID"].unique():
        print(i)
        print(mod["ID"])
        modx = mod[mod["ID"] == i]
        plt.figure(figsize=(16, 9))
        sns.lineplot(data=modx, x="variable", y="value", hue="ID", alpha=0.7, linewidth=4)
        plt.xticks(rotation=90)
        plt.xlabel("Read")
        plt.ylabel('%C modified')
        plt.title('')
        name = f'ReadInsertion_{i}_Modifications.png'
        plt.savefig(outpath +"/"+name, bbox_inches="tight")

def variant_bed_to_fasta(bed, outfasta):
    """
    Uses the variant callers sequence column output to generate a fasta from it.
    """
    variants = pd.read_csv(bed, sep='\t')
    variants = variants.iloc[:, [3,6]]
    with open(outfasta, 'w') as fasta:
         for i in range(len(variants)):
             if len(variants.iloc[i, 1]) >= 100:
                 fasta.write(">"+str(variants.iloc[i, 0])+"\n"+str(variants.iloc[i, 1]) + "\n")            

def blast2gff(blast,outfile):
    """
    Transforms blast output into gff format.
    """
    with open(blast, 'r') as blast_file, open(outfile, 'w') as gff_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            sequence_id, subject_id, start, end = fields[0], fields[1], fields[8], fields[9]
            gff_file.write(f"{subject_id}\t{start}\t{end}\tBLAST\tfeature\t.\t.\t{sequence_id}\n")

def reversevector(fastain, fastaout):
    with open(fastaout, 'w') as output_file:
        # opening given fasta file using the file path
        with open(fastain, 'r') as fasta_file:
            # extracting multiple data in single fasta file using biopython
            for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
                #split FASTA sequence into equally sized sub-lists with different id
                record_id = record.id
                record_seq = record.seq
                reversed_record_seq = record.seq[::-1]
                record_id_reversed = str(record_id) + '_reversed'
                output_file.write(">"+str(record_id)+"\n"+str(record_seq) + "\n")
                output_file.write(">"+str(record_id_reversed)+"\n"+str(reversed_record_seq) + "\n") 

def plot_insertion_length(bed, outfile):
    """
    Takes the full coordinates bed files and plots them by their length.
    """
    plt.figure(figsize=(16, 9))
    
    dfs = list()
    for f in bed:
        data = pd.read_csv(f, sep='\t', names=["Chr", "Start", "End", "Read"], usecols=[0,1,2,3])
        #add number to reads with multiple insertions so they don't overlap in the plot
        mask = data['Read'].duplicated(keep=False)
        data.loc[mask, 'Read'] += data.groupby('Read').cumcount().add(1).astype(str) #adds 1/2/3 respectively
        data['Read'] = data['Read'] + "_" + data['Chr']
        data["Length"] = data["End"] - data["Start"]
        head, tail = os.path.split(f)
        data['ID'] = tail.split(".")[0]
        dfs.append(data)
    
    df = pd.concat(dfs, ignore_index=True)
    df2 = df.copy()
    df2["Length"] = 0
    df_all = pd.concat([df,df2], ignore_index=True)
    sns.pointplot(data=df_all, x="Length", y="Read", hue='ID', linestyle="None", marker="_", legend=None, linewidth=2)
    sns.scatterplot(data=df, x="Length", y="Read", hue='ID', linestyle="None", marker="o")
    plt.legend(loc='lower center', bbox_to_anchor=[0.5, 1])
    plt.savefig(outfile, bbox_inches="tight")
            
def proximity_generator_for_bed_file(input_bed, output_bed, offsets):
    # Read the BED file into a DataFrame
    bed_df = pd.read_csv(input_bed, sep='\t', header=None, usecols=[0,1,2,3], names=['chrom', 'start', 'end', 'read'])
    print(bed_df)
    bed_df = bed_df[bed_df['chrom'].str.contains('chr')]
    bed_df = bed_df[~bed_df['chrom'].str.contains('_')]
    bed_df2 = bed_df.copy()

    # Iterate through the list of offsets and modify the genomic locations
    for offset in offsets:
        #downstream
        bed_df['start'] = bed_df['start'] - offset
        bed_df['end'] = bed_df['end']
        bed_df["offset"] = "-" + str(offset)
        # Write the modified DataFrame back to the original BED file
        bed_df.to_csv(output_bed, sep='\t', header=False, index=False, mode='a')  # 'a' for append
    for offset in offsets:
        #upstream
        bed_df2['start'] = bed_df2['start']
        bed_df2['end'] = bed_df2['end'] + offset
        bed_df2["offset"] = "+" + str(offset)
        # Write the modified DataFrame back to the original BED file
        bed_df2.to_csv(output_bed, sep='\t', header=False, index=False, mode='a')  # 'a' for append

def reshape_functional_tables(input_bed, output_bed):
    """
    Reformat bed-like output from TF and GENE overlays to get better table
    """
    bed_df = pd.read_csv(input_bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'read','distance','symbol'])
    bed_df=bed_df.drop_duplicates()
    bed_df['abs_distance'] = abs(bed_df["distance"])
    bed_df = bed_df.loc[bed_df.groupby(['read','symbol']).abs_distance.idxmin()]
    print(bed_df.head())
    bed_df['entry'] = range(len(bed_df))
    pivoted = bed_df.pivot(index=["entry","chrom","start","end","read"], columns="distance", values="symbol")
    pivoted.to_csv(output_bed, sep='\t')

def reads_with_insertions(infasta,border_dict,outfasta):
    '''
    Opens fasta and returns a fasta with only those reads that have an id that matches the border dict keys
    '''
    border_dict = json.load(open(border_dict))
    with open(outfasta, 'w') as output_file: #append mode, if it does not work, just store the orfs in list an then paste into file
        # opening given fasta file using the file path
        with open(infasta, 'r') as fasta_file:
            # extracting multiple data in single fasta file using biopython
            for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
                if record.id.split("_")[0] in border_dict:
                    output_file.write(">"+str(record.id)+"\n"+str(record.seq) + "\n")

def summary_reads_with_insertions(fullvectorsequence, vector, precutfasta, postcutfasta, border_dict, outpath):
    '''
    Opens vector fasta, precut fasta, and postcut fasta and combines read-sepcific the contents into single files, that can be used as an input for clustal omega 
    '''
    with open(fullvectorsequence, 'r') as vectorseq:
        for record in SeqIO.parse(vectorseq, 'fasta'):
            border_dict = json.load(open(border_dict))
            inputfiles = [vector, precutfasta, postcutfasta]
            processed = []
            for i in inputfiles:
                i = SeqIO.to_dict(SeqIO.parse(i, "fasta"))
                i = {k.split("_")[0]: v for k, v in i.items()} #dict comprehension
                keyfiltered= {key: i[key] for key in border_dict.keys()}
                processed.append(keyfiltered)
            
            for k in border_dict.keys():
                outfasta = outpath + "/" + str(k) +".fasta"  
                with open(outfasta, 'w') as output_file: #append mode, if it does not work, just store the orfs in list an then paste into file
                    output_file.write(">"+str(record.id)+"\n"+str(record.seq) + "\n")
                    output_file.write(">"+str(processed[0][k].id)+"\n"+str(processed[0][k].seq) + "\n")
                    output_file.write(">"+str(processed[1][k].id)+"\n"+str(processed[1][k].seq) + "\n")
                    output_file.write(">"+str(processed[2][k].id)+"\n"+str(processed[2][k].seq) + "\n")


def orf_reshape(orf_fasta, filename):
    """
    Uses FASTA-like output from orf finder and reshapes it into bed-like file for plotting. - Strandedness is added if stop > start in the ORF. 
    """
    with open(filename, 'w') as output_file:
        with open(orf_fasta, 'r') as fasta_file:
            # extracting multiple data in single fasta file using biopython
            for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
                #split FASTA sequence into equally sized sub-lists with different id
                orf_id = record.description.split(" ")[1]
                read_id = orf_id.split(":")[0].split("_")[1]
                orf_number = orf_id.split("_")[0]
                orf_start=orf_id.split(":")[1]
                orf_stop=orf_id.split(":")[2]
                orf_seq=record.seq
                if int(orf_start) > int(orf_stop):
                    print(orf_start, orf_stop)
                    output_file.write(str(read_id)+"\t"+str(orf_stop)+"\t"+str(orf_start)+"\t"+"-"+"\t"+str(orf_seq)+"\t"+str(orf_number)+"\n")
                else:
                    output_file.write(str(read_id)+"\t"+str(orf_start)+"\t"+str(orf_stop)+"\t"+"+"+"\t"+str(orf_seq)+"\t"+str(orf_number)+"\n") 

def full_coordinates_bed(border_dict, global_bed, outbed):
    """
    Uses the global bed file and the coordinates from the insertion dictionary and combined its content
    """
    bed = pd.read_csv(global_bed, sep='\t', header=None, usecols=[0,1,2,3])
    border_dict = json.load(open(border_dict))
    keys_to_filter = list(border_dict.keys())
    filtered_bed = bed[bed[3].str.split("_", expand=True)[0].isin(keys_to_filter)]
    filtered_bed['CleavageSites'] = filtered_bed[3].str.split("_", expand=True)[0].map(border_dict)
    filtered_bed.to_csv(outbed, sep='\t', index=False, header=False)


def find_longest_interval(matches, buffer):
    """
    Uses the list-transformed blast entries for each query id and returns the longest consecutive interval. Buffer defines a number of bases than can be jumped. without the break of the interval!
    """
    # Sort matches based on query_start
    matches.sort(key=lambda x: x['query_start'])
    
    longest_start, longest_end = None, None
    current_start, current_end = matches[0]['query_start'], matches[0]['query_end']
    longest_length = current_end - current_start
    longest_subject_ids = [matches[0]['subject_id'].split("_")[-1]]
    current_subject_ids = [matches[0]['subject_id'].split("_")[-1]]
    
    for match in matches[1:]:
        if match['query_start'] <= current_end + buffer:
            current_end = max(current_end, match['query_end'])
            current_subject_ids.append(match['subject_id'].split("_")[-1])
        else:
            current_length = current_end - current_start
            if current_length > longest_length:
                longest_start, longest_end = current_start, current_end
                longest_length = current_length
                longest_subject_ids = current_subject_ids[:]
            current_start, current_end = match['query_start'], match['query_end']
            current_subject_ids = [match['subject_id'].split("_")[-1]]
    
    current_length = current_end - current_start
    if current_length > longest_length:
        longest_start, longest_end = current_start, current_end
        longest_length = current_length
        longest_subject_ids = current_subject_ids[:]

    return longest_start, longest_end, longest_subject_ids

def plot_longest_interval(matches, longest_start, longest_end, longest_subject_ids, outfile):
    """
    PLotting of the longest interval.
    """
    plt.figure(figsize=(16, 2))
    # Plot the underlying sequence as a line
    interval=longest_end - longest_start
    plt.plot([0, interval], [0, 0], color='#E9967A', linewidth=2, label=matches[0]['query_id'])
    #plt.scatter(interval + interval*0.005 ,0,marker=">", s=100, color='red')
    # Plot the subject IDs as markers along the line
    for match in matches:
        if longest_start <= match['query_start'] <= longest_end:
            x_position = match['query_start'] - longest_start
            plt.scatter(x_position, 0, color='black',s=100, marker="|")
            plt.text(
                x_position, 
                0.01, 
                match['subject_id'].split("_")[-1], 
                ha='center', 
                va='bottom', 
                fontsize=5, 
                rotation=90
            )
    
    
    #invisible y axis
    plt.box(False)
    plt1 = plt.gca()
    plt1.axes.get_yaxis().set_visible(False)
    # Set plot limits and labels
    plt.xlim(-interval*0.1, interval + interval*0.1)
    plt.ylim(-0.05, 0.1)
    plt.ylabel('Matches')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(outfile, dpi=600)
    plt.show()

def find_and_plot_longest_blast_interval(blastn, buffer, threshold, output_dir):
    """
    Uses the blast result table and creates an illustration of the longest blast interval, including the description of the matching part of the vector. 
    """
    df = pd.read_csv(blastn, sep='\t')
    # Extract matches and group by QueryID
    query_groups = df.groupby('QueryID')

    # Iterate over each QueryID group to find the longest interval and plot
    for query_id, group in query_groups:
        matches = []
        for idx, row in group.iterrows():
            matches.append({
                'query_id': row['QueryID'],
                'subject_id': row['SubjectID'],
                'query_start': row['QueryStart'],
                'query_end': row['QueryEnd'],
                'subject_start': row['SubjectStart'],
                'subject_end': row['SubjectEnd']
            })
        
        # Find the longest interval
        longest_start, longest_end, longest_subject_ids = find_longest_interval(matches, buffer)
        print(f"QueryID: {query_id} - Longest interval: {longest_start} - {longest_end}")
        print(f"Ordered list of matching parts: {longest_subject_ids}")
        

        # Plot the longest interval
        if (longest_start is not None) and (longest_end is not None):
            if longest_end - longest_start > threshold:
                print(f"Longest interval for {query_id} is longer than {threshold}. Proceeding...")
                out_file = f'{output_dir}/Longest_interval_{query_id}.png'
                plot_longest_interval(matches, longest_start, longest_end, longest_subject_ids, out_file)
        else:
            print(f"Not enough matches for {query_id}.")

def plot_element_distance(bed, distances, distance_threshold, output_path):
    """
    Uses the bed file from the distance calculations and returns a plot to visualize the respective elements with their distance.
    Entries farther than the defined threshold are excluded.
    
    Parameters:
        bed (str): Path to the input BED file.
        distances (list): List of distances to define x-axis ticks.
        output_path (str): Path to save the plot.
        distance_threshold (int, optional): Maximum distance to include in the plot.
    """
    # Read the table
    df = pd.read_csv(
        bed,
        sep='\t',
        header=None,
        names=["chr", "start_insertion", "stop_insertion", "read", "source", "element_name", "distance"],
    )
    
    # Apply threshold filtering if provided
    print(distance_threshold)
    print(type(distance_threshold))
    if distance_threshold is not None:
        df = df[df['distance'].abs() <= int(distance_threshold)]

    # Ensure absolute distance and sort by absolute distance within groups
    df['abs_distance'] = df['distance'].abs()
    df = df.sort_values(by=['read', 'abs_distance']).drop_duplicates(subset=['read', 'element_name'], keep='first').reset_index(drop=True)
    
    # Prepare data for plotting
    plt.figure(figsize=(12, 6))
    
    # Create scatter plot
    sns.scatterplot(
        data=df,
        x='distance',
        y='element_name',
        hue='read',
        palette='tab10',
        s=100,
        style='source'
    )
    
    # Highlight genes with zero distance
    zero_distance_genes = df[df['distance'] == 0]['element_name']
    for gene in zero_distance_genes:
        plt.axhline(y=df[df['element_name'] == gene].index[0], color='gray', linestyle='--', linewidth=0.8)
    
    # Plot binned rugplot for distances
    bin_size = 100  # Bin size for grouping distances
    df['distance_bin'] = (df['distance'] // bin_size) * bin_size
    sns.rugplot(x=df['distance_bin'], color='black', height=0.05, linewidth=1)
    
    # Configure plot aesthetics
    plt.xticks(sorted({x for n in distances for x in (n, -n)}), rotation=45)
    plt.xlabel("Distance (bp)")
    plt.ylabel("Element Name")
    plt.title("Distance Distribution to Elements")
    sns.despine()
    plt.legend(title="",  bbox_to_anchor=(0.5, -0.2),  loc='upper center', fontsize=8)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_path}")

def join_read_mapq(file_list, prefixes, output_file):
    """
    Combine multiple text files with specified column prefixes, handling different row counts.

    Args:
        file_list (list): List of file paths to combine.
        prefixes (list): List of prefixes for the column titles.
        output_file (str): Path to save the combined output.

    Returns:
        pd.DataFrame: Combined DataFrame.
    """

    if len(file_list) != len(prefixes):
        raise ValueError("The number of files and prefixes must be the same.")
    
    combined_df = None
    
    for file, prefix in zip(file_list, prefixes):
        print(file, prefix)
        # Load the file and set custom column names with the prefix
        df = pd.read_csv(file, sep=" ", header=None, names=["Read", f"{prefix}_InsertedChr", f"{prefix}_MAPQ"])
        
        # Truncate the Read column to the part before '_'
        df["Read"] = df["Read"].str.split('_').str[0]
        
        # Merge into the combined dataframe
        if combined_df is None:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on="Read", how="outer")
    
    # Save the combined DataFrame to a file
    combined_df.to_csv(output_file, sep="\t", index=False)

def plot_mapq_changes(input_file, output_file):
    """
    Generates a line plot showing the changes in MAPQ values across the three stages 
    (Precut, Postcut, and Postcut_filtered) for each read individually.
    """
    # Load the data into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t")

    # Prepare the data for plotting
    plot_data = df[['Read', 'Precut_MAPQ', 'Postcut_MAPQ', 'Postcut_filtered_MAPQ']]
    
    # Reshape the data to have one row per Read and one column per MAPQ stage
    plot_data_melted = plot_data.melt(id_vars=['Read'], var_name='Stage', value_name='MAPQ')

    # Set the order for the stages on the x-axis
    stage_order = ['Precut_MAPQ', 'Postcut_MAPQ', 'Postcut_filtered_MAPQ']
    plot_data_melted['Stage'] = pd.Categorical(plot_data_melted['Stage'], categories=stage_order, ordered=True)

    # Plotting the lineplot
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=plot_data_melted, x='Stage', y='MAPQ', hue='Read', marker='o', dashes=False, markersize=4, alpha=0.5, legend=None)

    # Add jittered scatter points on top of the lineplot
    sns.scatterplot(data=plot_data_melted, x='Stage', y='MAPQ', hue='Read', alpha=0.5, marker='o', s=40)

    plt.title('MAPQ Changes Across Stages for Each Read')
    plt.xlabel('Stage')
    plt.ylabel('Mapping Quality (MAPQ)')
    plt.xticks(rotation=45)
    plt.legend(title="Read", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save the plot to the output file
    plt.savefig(output_file)
    plt.close()