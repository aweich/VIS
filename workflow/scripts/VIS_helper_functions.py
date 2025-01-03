#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:10:28 2023

@author: weichan
"""

import sys
import os
import re
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import collections
import seaborn as sns
import pybedtools
import json
import subprocess
#from Bio.Align.Applications import ClustalOmegaCommandline
from matplotlib.patches import Patch

#wrapper
from functools import wraps
import inspect

#wrapper for logging
def redirect_logging(logfile_param="logfile"):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Inspect the arguments of the function
            func_signature = inspect.signature(func)
            bound_args = func_signature.bind(*args, **kwargs)
            bound_args.apply_defaults()
            
            # Get the logfile path from the function arguments
            logfile = bound_args.arguments.get(logfile_param)
            if not logfile:
                raise ValueError(f"The parameter '{logfile_param}' must be provided with a valid file path.")
            
            # Ensure the directory for the logfile exists
            os.makedirs(os.path.dirname(logfile), exist_ok=True)
            
            # Open the logfile for writing
            with open(logfile, 'w') as log:
                # Redirect stdout and stderr
                original_stdout = sys.stdout
                original_stderr = sys.stderr
                sys.stdout = log
                sys.stderr = log
                
                try:
                    # Execute the wrapped function
                    result = func(*args, **kwargs)
                finally:
                    # Restore stdout and stderr
                    sys.stdout = original_stdout
                    sys.stderr = original_stderr
                
                return result
        return wrapper
    return decorator


#main
def chunks(lst, n):
    """
    Returns list of n-sized lists in the nested format: [[],[]]
    """
    collection = []
    for i in range(0, len(lst), n):
        collection.append(lst[i:i + n])
    return collection

@redirect_logging(logfile_param="logfile")
def fragmentation_fasta(fasta, equal_fragments, outfilename, logfile):
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

@redirect_logging(logfile_param="logfile")
def plot_bed_files_as_heatmap(bed_files, outfile, logfile):
    """
    Creates heatmap for chromosome-specific density of matches in BED files across samples.
    """
    data = {}
    # Process each BED file
    for file_id, bed_file in enumerate(bed_files):
        # Read the BED file into a DataFrame
        df = pd.read_csv(bed_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Read_ID', 'X', 'XX'])

        # Count chromosome occurrences
        chromosome_counts = df['Chromosome'].value_counts()

        # Store the counts in the data dictionary
        head, tail = os.path.split(bed_file)
        data[tail.split(".")[0]] = chromosome_counts

    # Create a DataFrame from the data
    counts_df = pd.DataFrame(data).fillna(0)

    # Debugging: Print the DataFrame to ensure correctness
    print("Counts DataFrame:")
    print(counts_df)

    if counts_df.empty:
        print("No data to plot. The counts DataFrame is empty.")
        return

    # Create the heatmap
    try:
        if len(counts_df.columns) > 1:
            # Attempt clustering heatmap
            plt.figure(figsize=(16, 9))
            sns.clustermap(counts_df, cmap="YlGnBu", annot=True, cbar_pos=(0, .2, .03, .4))
            plt.xlabel('File ID')
            plt.ylabel('Chromosome')
            plt.title('Chromosome Occurrences \n in BED Files')
        else:
            raise ValueError("Single column or insufficient data for clustering.")
    except ValueError as e:
        print(f"Clustering failed: {e}. Falling back to a simple heatmap.")
        plt.figure(figsize=(16, 9))
        sns.heatmap(counts_df, annot=True, cmap="YlGnBu")
        plt.xlabel('File ID')
        plt.ylabel('Chromosome')
        plt.title('Chromosome Occurrences \n in BED Files')

    # Save the heatmap
    plt.savefig(outfile, bbox_inches="tight")

####this part here is dedicated to the splitting of blast-match including fasta reads
def merge_intervals(intervals, overlap, filtering, filtervalue):
    # Sort intervals by start coordinates
    #sorted_intervals = sorted(intervals, key=lambda x: x[0])
    intervals.sort() #sorts all intervals in ascending order: Overlaps are possible! #sorts inplace
    print(intervals)

    merged_intervals = []
    
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

@redirect_logging(logfile_param="logfile")   
def splitting_borders(blast_file, filteroption, filtervalue,  overlap, outfile1, outfile2, logfile):
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
    #drop every second interval if there are more than two snippets of the fasta: -> To save only the intervals without insertion
    if len(sequences) > 2:
        del sequences[1::2]  
    
    #extraction of the insertion sequence
    insertion=[]
    
    for i in range(len(breakpoints)):
        if i % 2 == 0 or i == 0:
            insertion.append(fasta_string[breakpoints[i]:breakpoints[i+1]])
    
    return sequences, insertion

@redirect_logging(logfile_param="logfile")
def split_fasta_by_borders(border_dict, fasta, mode, outfasta, outinsertion, logfile):
    """
    Uses previously created breakpoints in border dict to cut out insertions from fasta file
    """
    border_dict = json.load(open(border_dict)) 
    with open(outfasta, 'w') as output_file:
        with open(outinsertion, 'w') as output_insertion_file:
            # opening given fasta file using the file path
            with open(fasta, 'r') as fasta_file:
                # extracting multiple data in single fasta file using biopython
                for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
                    #split FASTA sequence into equally sized sub-lists with different id
                    record_id = record.id
                    record_seq = record.seq
                    if record_id in border_dict:
                        record_list, insertion_list = split_fasta_with_breakpoints(record_seq, border_dict[record_id])
                        # write out insertion sequence
                        for i,entry in enumerate(insertion_list):
                            insertion_seq_new = insertion_list[i]
                            insertion_id_new = str(record_id) + '_%i' % (i)
                            print("Split " + str(record_id) + " into " + str(insertion_id_new))
                            output_insertion_file.write(">"+str(insertion_id_new)+"\n"+str(insertion_seq_new) + "\n")
                        # write out fasta without insertion        
                        if mode == "Buffer":
                            print("Buffer mode selected: Insertion sequences are replaced by N.")
                            print("Insertion coordinates will be reported as global coordinates relative to the reference genome.")
                            bufferlength =len(record) - len(''.join([str(n) for n in record_list]))
                            buffer = bufferlength*'N'
                            record_seq_new = buffer.join([str(n) for n in record_list])
                            print(record_seq_new)
                            record_id_new = str(record_id)+ '_Buffer'+str(len(buffer))+'Insertion'
                            output_file.write(">"+str(record_id_new)+"\n"+str(record_seq_new) + "\n")
                        elif mode == "Split":
                            print("Split mode selected: Insertion sequences are cut from the read and the read is split at the borders.")
                            print("Insertion coordinates will be reported as surrounding read coordinates.")
                            #this part only if the non-insertion fragments should NOT be combined
                            for i,entry in enumerate(record_list):
                                record_seq_new = record_list[i]
                                record_id_new = str(record_id) + '_%i' % (i)
                                print("Split " + str(record_id) + "into " + str(record_id_new))
                                output_file.write(">"+str(record_id_new)+"\n"+str(record_seq_new) + "\n")
                        else: 
                            print("No mode or unknown mode selected.")
                            print("Default: Insertions are cut from respective reads and the loose ends are joined together.")
                            print("Insertion coordinates will be reported as read coordinates")
                            record_seq_new = ''.join([str(n) for n in record_list])
                            record_id_new = str(record_id) + '_Insertion'
                            output_file.write(">"+str(record_id_new)+"\n"+str(record_seq_new) + "\n")
                    else:
                        output_file.write(">"+str(record_id)+"\n"+str(record_seq) + "\n")     

#cigar for exact cooridinates
def parse_cigar(cigar, start, strand):
    """
    Parses a CIGAR string to reconstruct the original FASTA from the alignment. 
    Soft-clipping is used to "correct" the start coordinates based on the real read length and the optimal alignment.
    This new adjusted_start can then be used to really localize the insertion in the reference genome on the base level. 
    
    Some addiitonal reasoning for the CIGAR caluclation: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
    Parameters:
        cigar (str): The CIGAR string (e.g., "50M5I30M10S").
        start (int): The genomic start position of the read (1-based).
        strand (str): The strand orientation ('+' or '-').

    Returns:
        list: A list of (genomic_position, operation) tuples.
        int: Adjusted genomic start position after considering soft clipping.
    """
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)  # Splits numbers and letters in CIGAR string
    genomic_pos = start
    positions = []

    # Adjust for soft-clipping to compute true start position: Import for insertions a the beginning or end of reads!
    left_soft_clip = 0
    right_soft_clip = 0

    if operations[0][1] == "S":  # Check for left soft-clipping
        left_soft_clip = int(operations[0][0])
    if operations[-1][1] == "S":  # Check for right soft-clipping
        right_soft_clip = int(operations[-1][0])

    adjusted_start = start - left_soft_clip

    print(operations)

    # Convert counts to integers and tally totals by operation type
    operation_totals = {}
    for length, operation in operations:
        length = int(length)
        if operation not in operation_totals:
            operation_totals[operation] = 0
        operation_totals[operation] += length

    # Display the results
    for op, total in operation_totals.items():
        print(f"Total {op}: {total}")



    for length, op in operations:
        length = int(length)
        if op in "M=X":  # Matches or alignment regions
            for _ in range(length):
                positions.append((genomic_pos, "M"))
                genomic_pos += 1
        elif op == "I":  # Insertions relative to the reference
            for _ in range(length):
                positions.append((None, "I"))  # No reference position
        elif op == "N":  # Skipped regions (splice junctions)
            for _ in range(length):
                #skipped regions do not "consume" query
                continue
        elif op == "D":  # Deletions in the reference
            for _ in range(length):
                #deletions do not "consume" query
                continue
        elif op == "S":  # Soft clipping
            for _ in range(length):
                positions.append((None, "S"))  # Clipped bases
        elif op == "H":  # Hard clipping
            for _ in range(length):
                #hard clipped should not exist in the FASTA 
                print("The CIGAR contains hard-clipped information although this should not exist at this point.")
                print("Make sure your FASTAs with the insertions can be reconstructed from their CIGARs after the alignment.")
                continue

    return positions, adjusted_start, start

@redirect_logging(logfile_param="logfile")
def reconstruct_coordinates(bed_with_cigar, fasta_coordinates, splitmode, output, logfile):
    """
    Reconstructs global coordinates of reads by translating FASTA-based positions using CIGAR strings.
    If there was not "Buffer" mode chosen for the split fasta function, the coordinates can not be re-constructed on base level accuracy.
    In case of "Separated", the coordinates of the split reads are reported.
    In case of "Join", the coordinates of the combined read are reported. This is the least accurate representation and is only recommended for developmental puproses so far.   

    Parameters:
        bed_with_cigar (str): Path to the BED file with CIGAR strings (from `bamtobed -cigar`).
        fasta_coordinates (dict): Dictionary of FASTA-based coordinates for reads.
        splitmode (str): Either Buffer or default to reporting the read coordinates. 
        output (str): Path to save the reconstructed BED file.

    """
    if logfile is None:
        raise ValueError("Logfile must be provided to log output.")

    bed = pd.read_csv(bed_with_cigar, sep='\t', header=None)
    bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'cigar']

    # Load FASTA coordinates
    fasta_coordinates = json.load(open(fasta_coordinates))
    insertionreads = set(fasta_coordinates.keys())
    
    # Normalize 'name' in BED to match the normalized FASTA keys
    bed['normalized_name'] = bed['name'].str.split("_").str[0] 

    # Subset BED to include only rows with matching normalized names
    bed_subset = bed[bed['normalized_name'].isin(insertionreads)]

    results=[]

    if splitmode == "Buffer":
        print("Reporting the coordinates of the insertions. Tracing the FASTA length from the CIGAR string...")
        # add condition for the Separated and options
        for _, row in bed_subset.iterrows():
            read_name = row['normalized_name']
            start = row['start']
            stop = row['end']
            strand = row['strand']
            cigar = row['cigar']
            chrom = row['chrom']
            origcoord = [start, stop]


            # Parse CIGAR string and adjust start position
            positions, adjusted_start, orig_start = parse_cigar(cigar, start, strand)
            # Map FASTA coordinates to BED coordinates
            fasta_coords = fasta_coordinates[read_name]
            fastalength = len(positions)
            
            for i in range(0, len(fasta_coords), 2):
                fasta_start = fasta_coords[i]
                fasta_end = fasta_coords[i + 1]

                # Translate FASTA to BED using CIGAR-based mapping
                if positions[0][1] == 'S':  # Handle soft-clipped regions
                    if strand == "+":
                        genomic_start = adjusted_start + fasta_start
                    else:  # Reverse strand handling
                        genomic_start = adjusted_start + fastalength - fasta_end
                else:  # Aligned regions
                    print("No soft clipping dected; No adjustment of start necessary")
                    genomic_start = orig_start + fasta_start

                if positions[0][1] == 'S':  # Handle soft-clipped regions
                    if strand == "+":
                        genomic_end = adjusted_start + fasta_end
                    else:  # Reverse strand handling
                        genomic_end = adjusted_start + fastalength - fasta_start
                else:  # Aligned regions
                    print("No soft clipping dected; No adjustment of end necessary")
                    genomic_end = orig_start + fasta_end
                
                # Debugging information
                print(f"Read: {read_name}, Strand: {strand}, Length: {fastalength}")
                print(f"FASTA Start: {fasta_start}, FASTA End: {fasta_end}")
                print(f"Genomic Start: {genomic_start}, Genomic End: {genomic_end}")
                print(f"Original coordinates: {origcoord}")
                results.append([chrom, genomic_start, genomic_end, read_name, origcoord, strand])

    elif splitmode == "Separated":
        print("Reporting surrounding-read coordinates.")
        for _, row in bed_subset.iterrows():
            read_name = row['normalized_name']
            start = row['start']
            stop = row['end']
            strand = row['strand']
            chrom = row['chrom']
            results.append([chrom, start, stop, read_name, "coordinates_of_surrounding_read", strand])
    else:
        print("Reporting read coordinates...")
        for _, row in bed_subset.iterrows():
            read_name = row['normalized_name']
            start = row['start']
            stop = row['end']
            strand = row['strand']
            chrom = row['chrom']
            results.append([chrom, start, stop, read_name, "coordinates_of_read_with_insertion", strand])
    

    # Save results to a new BED file
    reconstructed_bed = pd.DataFrame(results)
    reconstructed_bed.to_csv(output, sep='\t', index=False, header=False)


def blast2gff(blast,outfile):
    """
    Transforms blast output into gff format.
    """
    with open(blast, 'r') as blast_file, open(outfile, 'w') as gff_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            sequence_id, subject_id, start, end = fields[0], fields[1], fields[8], fields[9]
            gff_file.write(f"{subject_id}\t{start}\t{end}\tBLAST\tfeature\t.\t.\t{sequence_id}\n")

@redirect_logging(logfile_param="logfile")
def reverseinsertion(fastain, fastaout, logfile):
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

@redirect_logging(logfile_param="logfile")
def plot_insertion_length(bed, outfile, logfile):
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
    print("dataframe overview...")
    print(df.head())
    df2 = df.copy()
    df2["Length"] = 0
    df_all = pd.concat([df,df2], ignore_index=True)
    sns.pointplot(data=df_all, x="Length", y="Read", hue='ID', linestyle="None", marker="_", legend=None, linewidth=2)
    sns.scatterplot(data=df, x="Length", y="Read", hue='ID', linestyle="None", marker="o")
    plt.legend(loc='lower center', bbox_to_anchor=[0.5, 1])
    plt.savefig(outfile, bbox_inches="tight", dpi=600)

### this is for plotting each single read with a match and the matching parts            
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
    plt.figure(figsize=(8, 2))
    # Plot the underlying sequence as a line
    interval=longest_end - longest_start
    plt.plot([0, interval], [0, 0], color='red', linewidth=2, label=matches[0]['query_id'])
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
                fontsize=6, 
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

@redirect_logging(logfile_param="logfile")
def find_and_plot_longest_blast_interval(blastn, buffer, threshold, output_dir, logfile):
    """
    Uses the blast result table and creates an illustration of the longest blast interval, including the description of the matching part of the insertion. 
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

#functional
@redirect_logging(logfile_param="logfile")
def calculate_element_distance(insertions_bed, output_bed, logfile, annotation_files):
    """
    Calculates distances between insertion sites and genomic annotations using bedtools closest.
    
    Parameters:
    - insertions_bed (str): Path to the BED file containing the insertion sites.
    - output_bed (str): Path to save the output BED file.
    - annotation_files (str): List of paths to annotation BED files.
    """

    # At least one annotation input necessary
    if not annotation_files:
        raise ValueError("At least one annotation file must be provided.")
    # Create DataFrame for combined annotations
    combined_df = pd.DataFrame()

    for file in annotation_files:
        try:
            tag = os.path.basename(file).split(".")[0]  # Use file name without extension as tag
            df = pd.read_csv(file, sep="\t", header=None)  # Load as DataFrame
            df["source"] = tag  # Add source column
            combined_df = pd.concat([combined_df, df], ignore_index=True)
        except:
            continue
    
    # Convert combined DataFrame bed object
    combined_bed = pybedtools.BedTool.from_dataframe(combined_df)
    sorted_annotations = combined_bed.sort()

    insertions = pybedtools.BedTool(insertions_bed)

    #bedtools closest operation
    closest = insertions.closest(sorted_annotations, D="a", filenames=True)

    # Convert BedTool output to DataFrame
    closest_df = closest.to_dataframe(
        names=[
            "InsertionChromosome",
            "InsertionStart",
            "InsertionEnd",
            "InsertionRead",
            "InsertionOrig",
            "InsertionStrand",
            "AnnotationChromosome",
            "AnnotationStart",
            "AnnotationEnd",
            "AnnotationID",
            "AnnotationScore",
            "AnnotationStrand",
            "AnnotationSource",
            "Distance",
        ]
    )

    # Save DataFrame to a file with headers
    closest_df.to_csv(output_bed, sep="\t", index=False)
    print(f"Distances calculated and saved to {output_bed}")

@redirect_logging(logfile_param="logfile")
def plot_element_distance(bed, distances, distance_threshold, output_path, logfile):
    """
    Uses the bed file from the distance calculations and provides a plot to visualize the respective elements with their distance. Entries that are further away than the defined threshold value are excluded.
    """
    # Read the table
    df = pd.read_csv(
        bed,
        sep='\t',
    )
    
    print(df.head())
    # Apply threshold filtering if provided
    if distance_threshold is not None:
        df = df[df['Distance'].abs() <= int(distance_threshold)]
  
    # Ensure absolute distance and sort by absolute distance within groups
    df['abs_distance'] = df['Distance'].abs()
    df = df.sort_values(by=['InsertionRead', 'abs_distance']).drop_duplicates(subset=['InsertionRead', 'AnnotationID', "AnnotationSource"], keep='first').reset_index()
    
    # Create scatter plot
    sns.scatterplot(
        data=df,
        x='Distance',
        y='AnnotationID',
        hue='InsertionRead',
        palette='tab10',
        s=100,
        style='AnnotationSource'
    )
    
    # Binned rugplot for distances
    bin_size = 100  # Bin size grouping distances
    df['distance_bin'] = (df['Distance'] // bin_size) * bin_size
    sns.rugplot(x=df['distance_bin'], color='black', height=0.05, linewidth=1)
    
    # Configure plot aesthetics
    plt.xticks(sorted({x for n in distances for x in (n, -n)}), rotation=45)
    plt.xlabel("Distance (bp)")
    plt.ylabel("Element Name")
    plt.title("Distance Distribution to Elements")
    sns.despine()
    plt.legend(title="", fontsize=8)
    
    # Save plot
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_path}")

@redirect_logging(logfile_param="logfile")
def plot_element_distance_violin(bed, distances, distance_threshold, output_path, logfile):
    # Read the table
    df = pd.read_csv(
        bed,
        sep='\t',
    )

    # Apply threshold filtering if provided
    if distance_threshold is not None:
        df = df[df['Distance'].abs() <= int(distance_threshold)]

    # Create a count of how many times each source appears at each distance
    distance_counts = df.groupby(['Distance', 'AnnotationSource', 'InsertionRead']).size().reset_index(name='count')

    # Create the bar plot
    print(distance_counts.head())
    plt.figure(figsize=(10, 6))
    sns.displot(
        data=distance_counts,
        x='Distance', y='count', hue='InsertionRead', col='AnnotationSource',
        palette='Set2'
    )

    # Customize the plot
    plt.xticks(rotation=45)
    plt.xlabel("Distance (bp)")
    plt.ylabel("Count of Sources")
    plt.title("Distribution of Sources at Different Distances")
    sns.despine()

    # Save the plot
    plt.tight_layout()  # To ensure everything fits without overlap
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Plot saved to {output_path}")


@redirect_logging(logfile_param="logfile")
def scoring_insertions(data, output_file, logfile):
    """
    Uses custom conditions to visualize the entries of the annotated insertion summary table.
    """
    colnames = ["chr", "start_insertion", "stop_insertion", "read", "origcoord", "strand","chr_element","start_element", "stop_element","element_name", "element_score", "element_strand", "source", "distance"]
    df = pd.read_csv(
        data,
        sep='\t',
    )
    
    # Drop duplicate entries
    df = df.drop_duplicates().reset_index(drop=True)

    # Define conditions
    conditions = [
        ("COSMIC", 0), ("TF", 0), ("GENCODEV44", 0), ("HiC_Tcells", 0), ("Exons", 0),
        ("COSMIC", 10_000), ("TF", 10_000), ("GENCODEV44", 10_000), ("HiC_Tcells", 10_000),
        ("COSMIC", 50_000), ("TF", 50_000), ("GENCODEV44", 50_000), ("HiC_Tcells", 50_000),
        ("COSMIC", None), ("TF", None), ("GENCODEV44", None), ("HiC_Tcells", None),
    ]

    for source, distance in conditions:
        if distance == 0:
            df[f"{source}_0"] = df["source"].str.contains(source) & (df["distance"] == 0)
        elif distance == 10_000:
            df[f"{source}_10kb"] = df["source"].str.contains(source) & ((10_000 > abs(df["distance"])) & (abs(df["distance"]) > 0))
        elif distance == 50_000:
            df[f"{source}_50kb"] = df["source"].str.contains(source) & ((50_000 > abs(df["distance"])) & (abs(df["distance"]) > 10_000))
        else:
            df[f"{source}_Safe"] = df["source"].str.contains(source) & (abs(df["distance"]) > 50_000)

    # Aggregate data
    heatmap_data = df.groupby(["read", "chr", "start_insertion", "stop_insertion"]).agg("sum").reset_index()

    # Select numeric columns
    heatmap_matrix = heatmap_data.drop(columns=colnames)
    heatmap_matrix.index = heatmap_data["chr"] + "_" + \
                           heatmap_data["start_insertion"].astype(str) + "_" + \
                           heatmap_data["stop_insertion"].astype(str)

    # Calculate Final Score
    def calculate_score(row):
        if row["COSMIC_0"] > 0 or row["Exons_0"] > 0 or (row["TF_0"] + row["GENCODEV44_0"] > 1):
            return "Dangerous"
        elif (row["TF_0"] <= 1 or row["GENCODEV44_0"] <= 1 or row["COSMIC_10kb"] > 0 or row["HiC_Tcells_10kb"] > 0):
            return "Likely Dangerous"
        elif all(value > 0 for col, value in row.items() if "_Safe" in col) and row.sum() == row[[col for col in row.index if "_Safe" in col]].sum():
            return "Safe"
        else:
            return "Intermediate"

    heatmap_matrix["Risk"] = heatmap_matrix.apply(calculate_score, axis=1)

    # Map scores to colors
    score_colors = {
        "Dangerous": "red",
        "Likely Dangerous": "orange",
        "Intermediate": "yellow",
        "Safe": "green"
    }
    row_colors = heatmap_matrix["Risk"].map(score_colors)

    # Prepare row color map
    row_color_cmap = pd.DataFrame({
        "Risk": row_colors
    })

    # Plot clustermap with annotated row colors and custom colormap
    cluster_grid = sns.clustermap(
        heatmap_matrix.drop(columns=["Risk"]),
        cmap="Greys",
        row_colors=row_colors,
        figsize=(10, 5),
        dendrogram_ratio=(0.1, 0.1),
        linewidths=0.5,
        linecolor="grey",
        annot=True,
        col_cluster=False,
        row_cluster=False,
        clip_on=False,
        cbar_kws={
            "label": "Counts",
            "ticks": [0, 1, 2, 3, 4, 5],
            "shrink": 0.3,
            "orientation": "horizontal", 
        },
        vmin=0, vmax=5
    )

    # Add a custom legend
    legend_handles = [Patch(color=color, label=label) for label, color in score_colors.items()]
    cluster_grid.ax_heatmap.legend(
        handles=legend_handles,
        title="Risk Assessment",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),  # Position the legend below the plot
        ncol=4,  # Number of columns in the legend
        frameon=True
    )


    # Save the plot
    plt.savefig(output_file, bbox_inches="tight")

    # Optionally, print the heatmap matrix with Final Scores
    print(heatmap_matrix)

# qc
@redirect_logging(logfile_param="logfile")
def join_read_mapq(file_list, prefixes, output_file, logfile):
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
    print(combined_df.head())
    combined_df.columns = ["Read", "PrecutChr", "PrecutMAPQ", "PostcutChr", "PostcutMAPQ", "FilteredChr", "FilteredMAPQ"]
    combined_df.to_csv(output_file, sep="\t", index=False)

@redirect_logging(logfile_param="logfile")
def plot_mapq_changes(input_file, output_file, logfile):
    """
    Generates a line plot showing the changes in MAPQ values across the three stages 
    (Precut, Postcut, and Postcut_filtered) for each read individually.
    """
    # Load the data into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t")

    print(df.head())
    # Prepare the data for plotting
    plot_data = df[['Read', 'PrecutMAPQ', 'PostcutMAPQ', 'FilteredMAPQ']]
    
    # Reshape the data to have one row per Read and one column per MAPQ stage
    plot_data_melted = plot_data.melt(id_vars=['Read'], var_name='Stage', value_name='MAPQ')

    # Set the order for the stages on the x-axis
    stage_order = ['PrecutMAPQ', 'PostcutMAPQ', 'FilteredMAPQ']
    plot_data_melted['Stage'] = pd.Categorical(plot_data_melted['Stage'], categories=stage_order, ordered=True)

    # Plotting the lineplot
    plt.figure(figsize=(6, 6))
    sns.lineplot(data=plot_data_melted, x='Stage', y='MAPQ', hue='Read', marker='o', dashes=False, markersize=6, alpha=0.5)

    #plt.title('MAPQ Changes Across Stages for Each Read')
    plt.xlabel('Stage')
    plt.ylabel('Mapping Quality (MAPQ)')
    plt.xticks(rotation=45)
    plt.legend(title="Read", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save the plot to the output file
    plt.savefig(output_file, dpi=300)
    plt.close()

@redirect_logging(logfile_param="logfile")
def fragmentation_match_distribution(data, fragment_specifier, outpath, logfile):
    """
    Takes the insertion fragments and plots histogram of their frequency in the alignment
    """
    blasted = pd.read_csv(data, sep="\t")

    try: 
        if any(x.isupper() for x in blasted['QueryID'][0]) and "Read" not in blasted['QueryID'][0]: #to make sure the right column is used for plotting. Reads do not have any uppercase letters
            blasted[['Insertion', 'Fragment']] = blasted['QueryID'].str.split('_', n=1, expand=True)
            blasted["Fragment"] = pd.to_numeric(blasted["Fragment"])
            freq = collections.Counter(blasted["Fragment"].sort_values())
        else:
            blasted[['Insertion', 'Fragment']] = blasted['SubjectID'].str.split('_', n=1, expand=True)
            blasted["Fragment"] = pd.to_numeric(blasted["Fragment"])
            freq = collections.Counter(blasted["Fragment"].sort_values())
        
        plt.figure(figsize=(5, 4))
        plt.bar(freq.keys(), freq.values(), color='black')
        upperlimit = max(blasted["Fragment"])
        plt.xticks(np.arange(0, upperlimit+1, step=round(upperlimit/10)), fontsize=10)
        
        plt.ylabel('Count')
        plt.xlabel("Insertion Fragment")
        plt.title(f'Combined distribution of all {fragment_specifier} bp fragments')
        outfile = outpath + str("/") + f'{fragment_specifier}_fragmentation_distribution.png'
        plt.savefig(outfile, bbox_inches='tight', dpi=600)
        plt.close()
    except:
        print("The provided input could not be processed.")
        plt.figure()
        plt.title('Empty Data.')
        plt.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=12)
        plt.xlabel("No BLAST matches or read names not all lowercase")
        plt.xticks([])
        plt.yticks([])
        outfile = outpath + str("/") + f'{fragment_specifier}_fragmentation_distribution.png'
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
        return

@redirect_logging(logfile_param="logfile")
def fragmentation_read_match_distribution(data, fragment_specifier, outpath, logfile):
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

    plt.figure(figsize=(5, 4))
    if any(x.isupper() for x in blasted['QueryID'][0]) and "Read" not in blasted['QueryID'][0]:
        freq = collections.Counter(blasted["SubjectID"])
    else:
        freq = collections.Counter(blasted["QueryID"])
    plt.bar(freq.keys(), freq.values(), color='black')
    plt.xticks(rotation=90, fontsize=10)
    plt.ylabel('Fragment Count')
    plt.xlabel("Read")
    plt.title(f'Contribution of reads to the total count of {fragment_specifier} bp fragments')
    outfile = outpath + str("/") + f'{fragment_specifier}_read_match_fragmentation_distribution.png'
    plt.savefig(outfile, bbox_inches='tight', dpi=600)
    plt.close()
