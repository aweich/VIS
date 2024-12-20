## Simulation data for tutorial
#### Reference Genome

The first 500kb of chromosome 1 from the `hg38` assembly was used as the reference genome. The extracted sequence was saved as `chr1` in the FASTA format, which was then used for subsequent simulations and as reference genome for the [tutorial](tutorial.md/#prepare-configyml).

#### Insertion Sequence

The insertion sequence for the simulation was based on the [FTCAR2:pFlagCMV-mCAR-TVV vector construct](http://n2t.net/addgene:87239), obtained from Addgene (Plasmid # 87239). This construct comprises a vector system for the transduction of CXADR (Coxsackievirus-Adenovirus Receptor) into mammalian cells.

#### Simulation Process Overview

1. **Reference Genome Sampling:** 1000 reads were sampled from the reference genome with a mean read length of 10,000 bp, generated from a log-normal distribution.
2. **Insertion Introduction:** 5 of the 1000 reads were randomly chosen to receive an insertion, either the full-length or part of the construct, with random insertion directionality (`+` or `-`).
3. **Generation of Samples:** The process was repeated twice to generate two samples (`S1` and `S2`), each with a summary of the expected insertion locations.
4. **BAM generation:** Both samples were mapped against the reference genome to generate BAM files using:
 
```bash
minimap2 -ax map-ont chr1_1_50000_ref.fa S1.fa | samtools sort | samtools view -F 2304 -o S1.bam
```

This simulation provides a set of reads with varying insertion locations for evaluating the pipeline's detection and analysis capabilities.


### Annotation data

<details>
  <summary>Full Simulation Code:</summary>

```python
#!/usr/bin/env python3

import random
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq  # For reverse complement

# Paths to input files
reference_genome_path = "chr1_1_50000_ref.fa" # Reference genome
vector_sequence_path = "vectorseq.fa"   # Vector sequence

random.seed(2079)  # Seed for reproducibility

def collapse_fasta(fastapath):
    """Collapses multi-line FASTA sequences into a single string."""
    with open(fastapath, 'r') as fasta_file:
        seq_list = []
        for record in SeqIO.parse(fasta_file, 'fasta'):
            seq_list.append(str(record.seq))
    return ''.join(seq_list)

def generate_reads(reference_genome, mean_read_length=10000, num_reads=1000):
    """Generates random reads from the reference genome."""
    reads = []
    # Create read lengths following a log-normal distribution
    read_length_distribution = np.random.lognormal(mean=np.log(mean_read_length) - 0.5, sigma=1.0, size=num_reads)
    for read_length in read_length_distribution:
        read_length = int(read_length)  # Ensure lengths are integers
        start_position = random.randint(0, len(reference_genome) - read_length)
        read = reference_genome[start_position:start_position + read_length]
        reads.append(read)
    return reads

def add_insertions_to_reads(reads, insertion_sequence, num_insertions):
    """Adds insertions with strandedness to randomly selected reads."""
    insertion_summary = []  # To store insertion details for the summary
    for _ in range(num_insertions):
        # Randomly select a read to modify
        read_index = random.randint(0, len(reads) - 1)
        read = reads[read_index]

        # Decide whether to use the full insertion or a partial sequence
        if random.choice([True, False]):  # 50% chance
            insertion = insertion_sequence  # Full insertion
        else:
            start = random.randint(0, len(insertion_sequence) - 2000)
            end = start + random.randint(2000, min(5000, len(insertion_sequence) - start))
            insertion = insertion_sequence[start:end]  # Partial insertion

        # Decide the strandedness of the insertion
        if random.choice([True, False]):  # 50% chance for negative strand
            insertion = str(Seq(insertion).reverse_complement())  # Reverse complement for negative strand

        # Insert the sequence at a random position in the read
        insert_position = random.randint(0, len(read))
        modified_read = read[:insert_position] + insertion + read[insert_position:]

        # Replace the original read with the modified read
        reads[read_index] = modified_read

        # Store the details of the insertion
        strand = "-" if insertion != insertion_sequence else "+"
        insertion_summary.append((f'Read-{read_index + 1}', len(insertion), strand))

    return reads, insertion_summary

# Collapse the FASTA files into single sequences
reference_genome = collapse_fasta(reference_genome_path)
vector_sequence = collapse_fasta(vector_sequence_path)

# Generate reads and add insertions
non_insertion_reads = generate_reads(reference_genome)  # Reads without insertions
insertion_reads, insertion_summary = add_insertions_to_reads(non_insertion_reads.copy(), vector_sequence, 5)

# Combine reads and save them to a file
total_reads = insertion_reads + non_insertion_reads

output_file_path = "S2.fa"
with open(output_file_path, 'w') as output_file:
    for i, read in enumerate(total_reads):
        output_file.write(f'>Read-{i+1}\n{read}\n')

# Print the summary of insertions
print("Summary of Inserted Reads:")
for read_name, insertion_length, strand in insertion_summary:
    print(f"{read_name}: Insertion length = {insertion_length}, Strand = {strand}")

# Save summary to a file
summary_file_path = "S2_InsertionSummary.txt"
with open(summary_file_path, 'w') as summary_file:
    summary_file.write("Summary of Inserted Reads:\n")
    for read_name, insertion_length, strand in insertion_summary:
        summary_file.write(f"{read_name}: Insertion length = {insertion_length}, Strand = {strand}\n")

# Debugging outputs
print(f"Total reads: {len(total_reads)}")
```
</details>

<details>
  <summary>Summary of Inserted Reads for S1:</summary>
```
Read-221: Insertion length = 5711, Strand = +
Read-628: Insertion length = 2139, Strand = -
Read-536: Insertion length = 2018, Strand = -
Read-399: Insertion length = 5711, Strand = +
Read-46: Insertion length = 2137, Strand = -
```
</details>
<details>
  <summary>Summary of Inserted Reads for S2:</summary>

```
Read-389: Insertion length = 3328, Strand = -
Read-920: Insertion length = 2039, Strand = -
Read-532: Insertion length = 5711, Strand = +
Read-328: Insertion length = 5711, Strand = +
Read-347: Insertion length = 4109, Strand = -
```
</details>
<br> 

# show the vector sequence and how it is fragmented into 100bp sized-fragments