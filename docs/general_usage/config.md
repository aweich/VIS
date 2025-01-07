
# Configuration file

The configuration file is necessary to specify where the pipeline can find the input files that it needs for a proper execution. Below is a table with all currently implemented options, exemplary parameters, and a corresponding description for each field.  

An example for a ready-to-use `config.yml` is presented during the [tutorial](../tutorial/tutorial.md/#before-running-the-pipeline) of the pipeline. 

| **Parameter**            | **Exemplary value**                                                                                      | **Description**                                                                                              |
|---------------------------|------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| `experiment`             | `tutorial`                                                                                  | Name of the experiment.                                                                                 |
| `samples`                | `S1: /path/to/S1.bam` <br> `S2: /path/to/S2.bam` | Paths to BAM files for the samples. Each sample should have its file path defined.                      |
| `processing_dir`         | `/path/to/outdir`                                                                           | Directory where output files will be saved.                                                             |
| `insertion_fasta`           | `/path/to/insertionseq.fa`                                                      | Path to the insertion sequence file in FASTA format.                                                       |
| `blastn_db`              | `/path/to/blastNdb`                                                             | BLASTN database for reference nucleotides to check matches with the insertion sequence. If the insertion sequence contains reference genome parts, the matches here are needed to evaluate the true positive detections.                          |
| `splitmode`              | `Buffer`, `Split`, or `Join`                                                                                      | Mode used for processing reads based on the insertion: <ul><li><strong>Buffer:</strong> Replaces the insertion with "N" and retains the full length of each read. Recommended for the most accurate insertion location (CIGAR-based).</li><li><strong>Split:</strong> Cuts the insertion from the read and creates individual reads from the remaining sequence. The locations of the insertions are reported as the locations of each of the individual reads. Recommended for exploratory search for insertions that might fuse otherwise non-neighboring parts of the genome together.</li><li><strong>Join:</strong> Cuts the insertion from the read and joins the remaining sequence together. The locations of the insertions are reported as the location of the joined read. Recommended for debugging of otherwise unmappable reads.</li></ul> |
| `fragment_size`          | `100`                                                                                         | Size of fragments (in bp) for splitting sequences. Fragments of this size will be used to construct the BLASTN database of the insertion sequence.                                                              |
| `bridging_size`          | `300`                                                                                         | Acceptable gap size (in bp) before splitting the longest consecutive interval. It is recommended to customize this parameter according to the underlying read quality and insertion length.                                                               |
| `MinLength`              | `1`                                                                                           | Minimum read length (in bp) for BLASTN matches processing.                                                      |
| `MAPQ`                   | `10`                                                                                          | Minimum mapping quality score for reads. Filtering is applied after the modification of the reads with insertions.                                               |
| `MinInsertionLength`     | `500`                                                                                         | Minimum length (in bp) of insertions to be detected. This is dependent on the respective insertion and potentially its matches with the reference genome.                                                            |
| `ref_genome_ctrl`        | `/path/to/ref.fa`                                                     | Reference genome file in FASTA format.                                                      |
| `annotation_1`             | `/path/to/anno1.bed`                 | Sorted annotation file in [BED6 format](https://samtools.github.io/hts-specs/BEDv1.pdf).
| `annotation_2`             | `/path/to/anno2.bed`                 | Sorted annotation file in [BED6 format](https://samtools.github.io/hts-specs/BEDv1.pdf).
| `annotation_3`             | `/path/to/anno3.bed`                 | Sorted annotation file in [BED6 format](https://samtools.github.io/hts-specs/BEDv1.pdf).
| `annotation_4`             | `/path/to/anno4.bed`                 | Sorted annotation file in [BED6 format](https://samtools.github.io/hts-specs/BEDv1.pdf).                                                                                                                                                                                               |
| `detection`              | `rules/detection.smk`                                                                         | Snakemake rule file for the detection and localization of insertions.                                                                      |
| `quality_control`        | `rules/qc.smk`                                                                                | Snakemake rule file for the quality control of reads and insertions.                                                                |
| `functional_genomics`    | `rules/functional_genomics.smk`                                                               | Snakemake rule file for the functional annotation of insertions on the genome level.                                                            |


!!! info

    The parameters `blastn_db` and `annotation_2`/`annotation_3`/`annotation_4` are optional. 


