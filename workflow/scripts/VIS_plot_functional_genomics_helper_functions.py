#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

#wrapper
from VIS_helper_functions import redirect_logging

@redirect_logging(logfile_param="logfile")
def plot_element_distance(bed, distances, distance_threshold, output_path, logfile):
    """
    Uses the bed file from the distance calculations and provides a plot to visualize the respective elements with their distance. Entries that are further away than the defined threshold value are excluded.
    """
    # Read the table
    df = pd.read_csv(
        bed,
        sep='\t',
        header=None
    )
    
    print(df.head())
    df.columns = [
            "InsertionChromosome",
            "InsertionStart",
            "InsertionEnd",
            "InsertionRead",
            "InsertionOrigStart",
            "InsertionOrigEnd",
            "InsertionStrand",
            "AnnotationChromosome",
            "AnnotationStart",
            "AnnotationEnd",
            "AnnotationID",
            "AnnotationScore",
            "AnnotationStrand",
            "AnnotationSource",
            "Distance"
        ]	
   
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
	
    print(df.head())
    df.columns = [
            "InsertionChromosome",
            "InsertionStart",
            "InsertionEnd",
            "InsertionRead",
            "InsertionOrigStart",
            "InsertionOrigEnd",
            "InsertionStrand",
            "AnnotationChromosome",
            "AnnotationStart",
            "AnnotationEnd",
            "AnnotationID",
            "AnnotationScore",
            "AnnotationStrand",
            "AnnotationSource",
            "Distance"
        ]	
	
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
def scoring_insertions(data, output_plot, output_file, logfile):
    """
    Uses custom conditions to visualize the entries of the annotated insertion summary table.
    """
    colnames =[
            "InsertionChromosome",
            "InsertionStart",
            "InsertionEnd",
            "InsertionRead",
            "InsertionOrigStart",
            "InsertionOrigEnd",
            "InsertionStrand",
            "AnnotationChromosome",
            "AnnotationStart",
            "AnnotationEnd",
            "AnnotationID",
            "AnnotationScore",
            "AnnotationStrand",
            "AnnotationSource",
            "Distance"
        ]

    df = pd.read_csv(
        data,
        sep='\t',
        names=colnames
    )
    # Drop duplicate entries
    df = df.drop_duplicates().reset_index(drop=True)

    # Only distance = 0 matters here
    #df = df[df["Distance"] == 0]

    # Define conditions
    conditions = [
        ("cosmic", 0), ("tf", 0), ("intron", 0), ("hic", 0), ("exon", 0), ("promoter", 0)
    ]

    for source, distance in conditions:
        df[f"{source}_0"] = df["AnnotationSource"].str.contains(source) & (df["Distance"] == 0)

    # Aggregate data
    heatmap_data = df.groupby(["InsertionRead", "InsertionChromosome", "InsertionStart", "InsertionEnd"]).agg("sum").reset_index()

    print(heatmap_data)
    # Select numeric columns
    heatmap_matrix = heatmap_data.drop(columns=colnames)
    heatmap_matrix.index = heatmap_data["InsertionChromosome"] + "_" + \
                           heatmap_data["InsertionStart"].astype(str) + "_" + \
                           heatmap_data["InsertionEnd"].astype(str)


    # Calculate Final Score
    def calculate_score(row):
        if row.loc["intron_0"] < 1 and row.loc["exon_0"] < 1 and row.loc["promoter_0"] < 1:
            if row.loc["tf_0"] < 1 and row.loc["hic_0"] < 1:
                return "0"
            else:
                return "1"    
        elif row.loc["intron_0"] >= 1 and row.loc["exon_0"] < 1:
            if row.loc["tf_0"] >= 1 or row.loc["hic_0"] >= 1:
                if row.loc["cosmic_0"] < 1:
                    return "2"
                else:
                    return "4"
            else:
                if row.loc["cosmic_0"] < 1:
                    return "1"
                else:
                    return "4"
        elif row.loc["exon_0"] >= 1:
            if row.loc["tf_0"] > 1 or row.loc["hic_0"] > 1:
                if row.loc["cosmic_0"] < 1:
                    return "3"
                else:
                    return "4"
            else:
                if row.loc["cosmic_0"] < 1:
                    return "2"
                else:
                    return "4"
        elif row.loc["promoter_0"] >= 1:
            if row.loc["tf_0"] > 1 or row.loc["hic_0"] > 1:
                if row.loc["cosmic_0"] < 1:
                    return "3"
                else:
                    return "4"
            else:
                if row.loc["cosmic_0"] < 1:
                    return "3"
                else:
                    return "4"
        else:
            return "undefined"

    heatmap_matrix["Risk"] = heatmap_matrix.apply(calculate_score, axis=1)

    # Map scores to colors
    score_colors = {
        "4": "red",
        "3": "orange",
        "2": "yellow",
        "1": "lightgrey",
        "0": "green"
    }
    row_colors = heatmap_matrix["Risk"].map(score_colors)

    # Prepare row color map
    row_color_cmap = pd.DataFrame({
        "Risk": row_colors
    })

    #pretty names
    heatmap_matrix.columns = ["Cancer gene", "TF", "Intron", "HiC", "Exon", "Promoter", "Risk"]
    
    # Plot clustermap with annotated row colors and custom colormap
    cluster_grid = sns.clustermap(
        heatmap_matrix.drop(columns=["Risk"]),
        cmap="Greys",
        row_colors=row_colors,
        figsize=(10,10),
        dendrogram_ratio=(0.1, 0.1),
        linewidths=0.5,
        linecolor="grey",
        annot=True,
        col_cluster=False,
        row_cluster=False,
        clip_on=False,
        label=False,
        cbar_kws={
            "label": "Counts",
            "ticks": [0, 1, 2, 3, 4, 5],
            "shrink": 0.5,
            "orientation": "horizontal", 
        },
        vmin=0, vmax=5
    )

    # Add a "5+" label to the last tick
    cbar = cluster_grid.ax_cbar  # Access the color bar
    cbar.set_xticks([0, 1, 2, 3, 4, 5])  # Set tick positions
    cbar.set_xticklabels(["0", "1", "2", "3", "4", "5+"]) 


    # Add a custom legend
    legend_handles = [Patch(color=color, label=label) for label, color in score_colors.items()]
    cluster_grid.ax_heatmap.legend(
        handles=legend_handles,
        title="Risk Assessment",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),  # Position the legend below the plot
        ncol=5,  # Number of columns in the legend
        frameon=True
    )

    # Save the plot
    plt.savefig(output_plot, format="svg", bbox_inches="tight")

    # Heatmap matrix with Final Scores
    print(heatmap_matrix)
    heatmap_matrix.to_csv(output_file, sep="\t")
