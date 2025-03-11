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
    - annotation_files (dict): Dict with config_entry as keys and file paths as values.
    """

    # Ensure at least one annotation file is provided
    if not annotation_files:
        raise ValueError("At least one annotation file must be provided.")

    # Load insertions as a BedTool object
    insertions = pybedtools.BedTool(insertions_bed)

    # Store processed closest distances for all annotations
    closest_results = []

    for tag, file in annotation_files.items():
        try:
            df = pd.read_csv(file, sep="\t", header=None)

            # Check if the file has at least 4 columns
            if df.shape[1] < 4:
                raise ValueError(f"{file} has fewer than 4 columns. This is unexpected.")

            # If the file has only 5 columns, add a dummy strand column (+)
            while df.shape[1] < 6:
                df[df.shape[[1]]] = "."
            
            #makes sure no other irrelevant columns are introduced since the downstream scripts depend on the specific size
            df = df.iloc[:,:6] 
            
            # Ensure first column has "chr" to confirm it's a BED file
            if df.iloc[0, 0].startswith("chr"): 
                df["source"] = tag  # Add annotation source
                print(f"Loaded {tag}: {df.head()}")

                # df to bed and sort
                bed = pybedtools.BedTool.from_dataframe(df).sort()
            
                # bedtools closest 
                closest = insertions.closest(bed, D="a", t='first')
            
                # Convert result to df and add to results
                closest_results.append(closest.to_dataframe(header=None))

        except Exception as e:
            print(f"Error reading {file}: {e}")
            with open(logfile, "a") as log:
                log.write(f"Error reading {file}: {e}\n")
            continue
    
    print(closest_results)
    # Combine all results and save to output
    if closest_results:
        final_df = pd.concat(closest_results, ignore_index=True)
        final_df.to_csv(output_bed, sep="\t", index=False, header=None)
        print(f"Distances calculated and saved to {output_bed}")
    else:
        print("No valid annotations processed. Output file not created.")


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
