#!/usr/bin/env python3

# This script takes all multi-record GenBank files in a directory, and outputs
# concatenated single-entry files of the same format. By default, the entries are combined into
# a single sequence with 50 "N" characters separating each of the original entries, but this
# spacer can be adjusted or omitted.

# The script was adapted from a ResearchGate post reply by Sanjay Kumar:
# https://www.researchgate.net/post/How-can-I-load-a-multi-record-genbank-file-into-Artemis

import sys
import os
import glob
from os.path import join, basename, splitext
from Bio import SeqIO

path = sys.argv[1]
files = join(path, '*.gb*')
out_dir = join(path, "concatenated_seqs")
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for filename in glob.glob(files):
    file_base = splitext(basename(filename))[0]
    already_merged = file_base+".gb"

    if not os.path.exists(already_merged):
        merged_rec = ''

        for rec in SeqIO.parse(open(filename,"r"), "genbank") :
            # Insert a spacer to preserve contig boundaries:
            merged_rec += rec + ("N" * 50)
            #Or to merge contigs with no spacers:
            #merged_rec += rec

        merged_rec.id = file_base[:15] #limit id to first 15 characters to prevent error
        merged_rec.description = "merged seq"
        out_file = join(out_dir, file_base+'_merge.gb')
        SeqIO.write(merged_rec, out_file, "genbank")
