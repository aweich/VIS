#!/usr/bin/env python3

import pandas as pd
import pybedtools

#wrapper
from VIS_helper_functions import redirect_logging

#functional
@redirect_logging(logfile_param="logfile")
def calculate_element_distance(insertions_bed, output_bed, logfile, annotation_files):
    """
    Calculates distances between insertion sites and genomic annotations using bedtools closest.
    
    Parameters:
    - insertions_bed (str): Path to the BED file containing the insertion sites.
    - output_bed (str): Path to save the output BED file.
    - annotation_files (dict): Dict with config_entry as keys and pathways as values.
    """
    
    # At least one annotation input necessary
    if not annotation_files:
        raise ValueError("At least one annotation file must be provided.")
    
    # Create DataFrame for combined annotations
    combined_df = pd.DataFrame()

    for tag, file in annotation_files.items():
        try:
            df = pd.read_csv(file, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5])
            if df.iloc[0, 0].startswith("chr"): 
                df["source"] = tag  # Add source column
                print(f"Loaded {tag}: {df.head()}")
                combined_df = pd.concat([combined_df, df], ignore_index=True)
        except:
            print(f"Error reading {file})")
            print("BED files are expected to follow BED6 format convention. If more than 6 columns are provided, the first 6 will be used.")
            continue
    
    # Convert combined DataFrame bed object
    combined_bed = pybedtools.BedTool.from_dataframe(combined_df)
    sorted_annotations = combined_bed.sort()

    insertions = pybedtools.BedTool(insertions_bed)
    
    #bedtools closest operation
    closest = insertions.closest(sorted_annotations, D="a", filenames=True)

    print(type(closest))
    print(closest)
    closest.saveas(output_bed)
    print(f"Distances calculated and saved to {output_bed}")

@redirect_logging(logfile_param="logfile")
def run_bedtools_intersect(insertions_bed, output_files, logfile, annotations):
    """
    Runs bedtools intersect for each annotation using pybedtools.
    """
    insertions = pybedtools.BedTool(insertions_bed)
    
    for key, annotation_file in annotations.items():
        output_file = output_files[key]
        annotation = pybedtools.BedTool(annotation_file)
        
        intersected = insertions.intersect(annotation, wb=True)
        
        intersected.saveas(output_file)
        print(f"Completed: {key}")
