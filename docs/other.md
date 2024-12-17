# Output directory structure
```bash
├── config_settings.yml
├── final
│   ├── functional_genomics
│   │   ├── BarPlot_Distance_to_Genes_100_Sim1.png
│   │   ├── Functional_distances_to_Insertions_Sim1.bed
│   │   ├── Insertion_Scoring_Sim1.png
│   │   └── Plot_Distance_to_Genes_100_Sim1.png
│   ├── localization
│   │   ├── ExactInsertions_Sim1.bed
│   │   ├── Heatmap_Insertion_Chr.png
│   │   └── Insertion_length.png
│   └── qc
│       ├── Fragmentation
│       │   ├── Insertions
│       │   │   └── insertions_100_Sim1
│       │   │       ├── 100_fragmentation_distribution.png
│       │   │       └── 100_read_match_fragmentation_distribution.png
│       │   ├── Longest_Interval
│       │   │   └── Sim1
│       │   │       ├── Longest_interval_Read-1008.png
│       │   │       ├── Longest_interval_Read-1009.png
│       │   │       ├── Longest_interval_Read-1059.png
│       │   │       ├── Longest_interval_Read-1081.png
│       │   │       ├── Longest_interval_Read-1387.png
│       │   │       ├── Longest_interval_Read-1568.png
│       │   │       ├── Longest_interval_Read-1571.png
│       │   │       ├── Longest_interval_Read-1589.png
│       │   │       ├── Longest_interval_Read-1717.png
│       │   │       ├── Longest_interval_Read-1719.png
│       │   │       ├── Longest_interval_Read-1915.png
│       │   │       ├── Longest_interval_Read-365.png
│       │   │       ├── Longest_interval_Read-402.png
│       │   │       ├── Longest_interval_Read-41.png
│       │   │       ├── Longest_interval_Read-423.png
│       │   │       ├── Longest_interval_Read-424.png
│       │   │       ├── Longest_interval_Read-463.png
│       │   │       ├── Longest_interval_Read-641.png
│       │   │       ├── Longest_interval_Read-650.png
│       │   │       └── Longest_interval_Read-822.png
│       │   └── Reference
│       │       └── reference_100_Sim1
│       │           └── 100_fragmentation_distribution.png
│       └── multiqc_report.html
└── intermediate
    ├── blastn
    │   ├── 100_VectorMatches_Sim1.blastn
    │   ├── 100_VectorMatches_Sim1.gff
    │   ├── Annotated_100_VectorMatches_Sim1.blastn
    │   ├── CleavageSites_100_VectorMatches_Sim1.blastn
    │   ├── Filtered_Annotated_100_VectorMatches_Sim1.blastn
    │   ├── humanref
    │   │   ├── 100_VectorMatches_Sim1.gff
    │   │   ├── Annotated_100_VectorMatches_Sim1.blastn
    │   │   └── Filtered_Annotated_100_VectorMatches_Sim1.blastn
    │   └── Readnames_100_VectorMatches_Sim1.txt
    ├── fasta
    │   ├── Cleaved_Sim1_noVector.fa
    │   ├── fragments
    │   │   ├── 100_Vector_fragments.fa
    │   │   ├── 100_Vector_fragments.fa.ndb
    │   │   ├── 100_Vector_fragments.fa.nhr
    │   │   ├── 100_Vector_fragments.fa.nin
    │   │   ├── 100_Vector_fragments.fa.njs
    │   │   ├── 100_Vector_fragments.fa.not
    │   │   ├── 100_Vector_fragments.fa.nsq
    │   │   ├── 100_Vector_fragments.fa.ntf
    │   │   ├── 100_Vector_fragments.fa.nto
    │   │   └── Forward_Backward_Vector.fa
    │   ├── Full_Sim1.fa
    │   └── Insertion_Sim1_Vector.fa
    ├── localization
    │   ├── ExactInsertions_Sim1.bed
    │   └── Sorted_ExactInsertions_Sim1.bed
    ├── log
    │   ├── detection
    │   │   ├── BAM_to_BED
    │   │   │   ├── Postcut_Sim1.log
    │   │   │   └── Precut_Sim1.log
    │   │   ├── basic_insertion_plots
    │   │   │   ├── heat.log
    │   │   │   └── length.log
    │   │   ├── build_insertion_reference
    │   │   │   └── out.log
    │   │   ├── calculate_exact_insertion_coordinates
    │   │   │   └── Sim1.log
    │   │   ├── clean_postcut_by_maping_quality
    │   │   │   └── Sim1.log
    │   │   ├── collect_outputs
    │   │   │   └── Sim1.log
    │   │   ├── copy_config_version
    │   │   │   └── out.log
    │   │   ├── extract_by_length
    │   │   │   └── Sim1.log
    │   │   ├── find_vector_BLASTn
    │   │   │   └── Sim1.log
    │   │   ├── find_vector_BLASTn_in_humanRef
    │   │   │   └── Sim1.log
    │   │   ├── get_cleavage_sites_for_fasta
    │   │   │   └── Sim1.log
    │   │   ├── hardcode_blast_header
    │   │   │   └── Sim1.log
    │   │   ├── insertion_mapping
    │   │   │   └── Sim1.log
    │   │   ├── make_blastn_DB
    │   │   │   └── out.log
    │   │   ├── make_fasta_without_tags
    │   │   │   └── Sim1.log
    │   │   ├── minimap_index
    │   │   │   └── out.log
    │   │   ├── Non_insertion_mapping
    │   │   │   └── Sim1.log
    │   │   ├── prepare_vector
    │   │   │   └── out.log
    │   │   ├── split_fasta_by_borders
    │   │   │   └── Sim1.log
    │   │   └── vector_fragmentation
    │   │       └── out.log
    │   ├── functional_genomics
    │   │   ├── calc_distance_to_elements
    │   │   │   └── Sim1.log
    │   │   ├── plot_distance_to_elements
    │   │   │   ├── scatter_Sim1.log
    │   │   │   └── violin_Sim1.log
    │   │   ├── plot_scoring
    │   │   │   └── Sim1.log
    │   │   └── sort_insertion_file
    │   │       └── Sim1.log
    │   └── qc
    │       ├── detailed_fragmentation_length_plot
    │       │   └── Sim1.log
    │       ├── extract_fastq_insertions
    │       │   └── Sim1.log
    │       ├── extract_mapping_quality
    │       │   └── Sim1.log
    │       ├── finalize_mapping_quality
    │       │   └── Sim1.log
    │       ├── fragmentation_distribution_plots
    │       │   ├── fragmentation_match_distribution_Sim1.log
    │       │   └── fragmentation_read_match_distribution_Sim1.log
    │       ├── generate_mapq_heatmap
    │       │   └── Sim1.log
    │       ├── multiqc
    │       │   └── out.log
    │       ├── nanoplot
    │       │   └── Sim1.log
    │       └── read_level_fastqc
    │           └── Sim1.log
    ├── mapping
    │   ├── Postcut_Sim1.bed
    │   ├── Postcut_Sim1_sorted.bam
    │   ├── Postcut_Sim1_sorted.bam.bai
    │   ├── Postcut_Sim1_unfiltered_sorted.bam
    │   ├── Postcut_Sim1_unfiltered_sorted.bam.bai
    │   ├── Precut_Sim1.bed
    │   ├── Precut_Sim1_sorted.bam
    │   ├── Precut_Sim1_sorted.bam.bai
    │   └── vector_ref_genome.fa
    └── qc
        ├── fastqc
        │   ├── readlevel_Sim1
        │   │   ├── Sim1_read_Read-1008.fastq
        │   │   ├── Sim1_read_Read-1008_fastqc.html
        │   │   ├── Sim1_read_Read-1008_fastqc.zip
        │   │   ├── Sim1_read_Read-1009.fastq
        │   │   ├── Sim1_read_Read-1009_fastqc.html
        │   │   ├── Sim1_read_Read-1009_fastqc.zip
        │   │   ├── Sim1_read_Read-1059.fastq
        │   │   ├── Sim1_read_Read-1059_fastqc.html
        │   │   ├── Sim1_read_Read-1059_fastqc.zip
        │   │   ├── Sim1_read_Read-1081.fastq
        │   │   ├── Sim1_read_Read-1081_fastqc.html
        │   │   ├── Sim1_read_Read-1081_fastqc.zip
        │   │   ├── Sim1_read_Read-1387.fastq
        │   │   ├── Sim1_read_Read-1387_fastqc.html
        │   │   ├── Sim1_read_Read-1387_fastqc.zip
        │   │   ├── Sim1_read_Read-1568.fastq
        │   │   ├── Sim1_read_Read-1568_fastqc.html
        │   │   ├── Sim1_read_Read-1568_fastqc.zip
        │   │   ├── Sim1_read_Read-1571.fastq
        │   │   ├── Sim1_read_Read-1571_fastqc.html
        │   │   ├── Sim1_read_Read-1571_fastqc.zip
        │   │   ├── Sim1_read_Read-1589.fastq
        │   │   ├── Sim1_read_Read-1589_fastqc.html
        │   │   ├── Sim1_read_Read-1589_fastqc.zip
        │   │   ├── Sim1_read_Read-1717.fastq
        │   │   ├── Sim1_read_Read-1717_fastqc.html
        │   │   ├── Sim1_read_Read-1717_fastqc.zip
        │   │   ├── Sim1_read_Read-1719.fastq
        │   │   ├── Sim1_read_Read-1719_fastqc.html
        │   │   ├── Sim1_read_Read-1719_fastqc.zip
        │   │   ├── Sim1_read_Read-1915.fastq
        │   │   ├── Sim1_read_Read-1915_fastqc.html
        │   │   ├── Sim1_read_Read-1915_fastqc.zip
        │   │   ├── Sim1_read_Read-365.fastq
        │   │   ├── Sim1_read_Read-365_fastqc.html
        │   │   ├── Sim1_read_Read-365_fastqc.zip
        │   │   ├── Sim1_read_Read-402.fastq
        │   │   ├── Sim1_read_Read-402_fastqc.html
        │   │   ├── Sim1_read_Read-402_fastqc.zip
        │   │   ├── Sim1_read_Read-41.fastq
        │   │   ├── Sim1_read_Read-41_fastqc.html
        │   │   ├── Sim1_read_Read-41_fastqc.zip
        │   │   ├── Sim1_read_Read-423.fastq
        │   │   ├── Sim1_read_Read-423_fastqc.html
        │   │   ├── Sim1_read_Read-423_fastqc.zip
        │   │   ├── Sim1_read_Read-424.fastq
        │   │   ├── Sim1_read_Read-424_fastqc.html
        │   │   ├── Sim1_read_Read-424_fastqc.zip
        │   │   ├── Sim1_read_Read-463.fastq
        │   │   ├── Sim1_read_Read-463_fastqc.html
        │   │   ├── Sim1_read_Read-463_fastqc.zip
        │   │   ├── Sim1_read_Read-641.fastq
        │   │   ├── Sim1_read_Read-641_fastqc.html
        │   │   ├── Sim1_read_Read-641_fastqc.zip
        │   │   ├── Sim1_read_Read-650.fastq
        │   │   ├── Sim1_read_Read-650_fastqc.html
        │   │   ├── Sim1_read_Read-650_fastqc.zip
        │   │   ├── Sim1_read_Read-822.fastq
        │   │   ├── Sim1_read_Read-822_fastqc.html
        │   │   └── Sim1_read_Read-822_fastqc.zip
        │   └── Sim1_filtered.fastq
        ├── mapq
        │   ├── Insertions_Sim1_mapq.txt
        │   └── Sim1_mapq_heatmap_image.png
        ├── multiqc_data
        │   ├── fastqc_adapter_content_plot.txt
        │   ├── fastqc_overrepresented_sequences_plot.txt
        │   ├── fastqc_per_base_n_content_plot.txt
        │   ├── fastqc_per_base_sequence_quality_plot.txt
        │   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
        │   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
        │   ├── fastqc_per_sequence_quality_scores_plot.txt
        │   ├── fastqc_sequence_counts_plot.txt
        │   ├── fastqc_sequence_duplication_levels_plot.txt
        │   ├── fastqc-status-check-heatmap.txt
        │   ├── fastqc_top_overrepresented_sequences_table.txt
        │   ├── multiqc_citations.txt
        │   ├── multiqc_data.json
        │   ├── multiqc_fastqc.txt
        │   ├── multiqc_general_stats.txt
        │   ├── multiqc.log
        │   ├── multiqc_nanostat.txt
        │   ├── multiqc_software_versions.txt
        │   ├── multiqc_sources.txt
        │   ├── nanostat_aligned_stats_table.txt
        │   └── nanostat_quality_dist.txt
        ├── multiqc_report.html
        └── nanoplot
            └── Sim1
                ├── AlignedReadlengthvsSequencedReadLength_dot.html
                ├── AlignedReadlengthvsSequencedReadLength_dot.png
                ├── AlignedReadlengthvsSequencedReadLength_kde.html
                ├── AlignedReadlengthvsSequencedReadLength_kde.png
                ├── MappingQualityvsReadLength_dot.html
                ├── MappingQualityvsReadLength_dot.png
                ├── MappingQualityvsReadLength_kde.html
                ├── MappingQualityvsReadLength_kde.png
                ├── NanoPlot_20241217_1605.log
                ├── NanoPlot-report.html
                ├── NanoStats.txt
                ├── Non_weightedHistogramReadlength.html
                ├── Non_weightedHistogramReadlength.png
                ├── Non_weightedLogTransformed_HistogramReadlength.html
                ├── Non_weightedLogTransformed_HistogramReadlength.png
                ├── PercentIdentityHistogramDynamic_Histogram_percent_identity.html
                ├── PercentIdentityHistogramDynamic_Histogram_percent_identity.png
                ├── PercentIdentityvsAlignedReadLength_dot.html
                ├── PercentIdentityvsAlignedReadLength_dot.png
                ├── PercentIdentityvsAlignedReadLength_kde.html
                ├── PercentIdentityvsAlignedReadLength_kde.png
                ├── WeightedHistogramReadlength.html
                ├── WeightedHistogramReadlength.png
                ├── WeightedLogTransformed_HistogramReadlength.html
                ├── WeightedLogTransformed_HistogramReadlength.png
                ├── Yield_By_Length.html
                └── Yield_By_Length.png

63 directories, 214 files
```