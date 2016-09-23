#!/usr/bin/python

from Bio import SeqIO
import sys

#Usage: python filter_contigs_by_size.py [infile] [size (bp)] [outfile]

fasta = SeqIO.parse(sys.argv[1], 'fasta')

min_size = int(sys.argv[2])

to_keep = []
for f in fasta:
	if len(f.seq) > min_size:
		to_keep.append(f)


SeqIO.write(to_keep, sys.argv[3], 'fasta')
