#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:10:28 2023

@author: weichan
"""

import sys
from Bio import SeqIO


def parse_FASTQ(fastq_read):
    id = fastq_read.id
    seq = fastq_read.seq
    #scan for matches
    #cut out match
    #create new ids based on the match _1, _2, ... _n
    #create new read files
    print(id, seq)
    return
