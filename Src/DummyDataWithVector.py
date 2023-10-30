#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:57:35 2023

@author: weichan
"""
import random
import sys

# Define the parameters
read_length = 5000  # Adjust the read length as needed
coverage = 40  # Adjust the coverage as needed

dummy=''.join(random.choice('ACGT') for _ in range(10000))
with open("/home/weichan/permanent/Projects/VIS/dev/DummyVector.fa", 'w') as output_file:
         output_file.write("> Dummy\n" + dummy)

# Generate random long-read sequences
reads = []
for _ in range(coverage):
    # Create a random long read without the insert sequence
    read = ''.join(random.choice('ACGT') for _ in range(read_length))

    # Randomly choose the insertion point
    #insertion_point = random.randint(0, read_length - 1)

    # Insert the custom sequence at the chosen point
    #read = read[:insertion_point] + dummy + read[insertion_point:]

    reads.append(read)
# 20 20 20; 40 with full insertion, 20 with half insertion, 40 without insertion
for read in reads[:20]: #full vector insertion
     # Randomly choose the insertion point
    insertion_point = random.randint(0, read_length - 1)

    # Insert the custom sequence at the chosen point
    read = read[:insertion_point] + dummy + read[insertion_point:]

    reads.append(read)

for read in reads[20:40]: #first half vector insertion
     # Randomly choose the insertion point
    insertion_point = random.randint(0, read_length - 1)

    # Insert the custom sequence at the chosen point
    read = read[:insertion_point] + dummy[:len(dummy)//2] + read[insertion_point:]

    reads.append(read)
    
# Save the generated reads
with open("/home/weichan/permanent/Projects/VIS/dev/DummyData.fa", 'w') as output_file:
    for i, read in enumerate(reads):
         output_file.write(f'>Read_{i+1}\n{read}\n')