#!/usr/local/bin/python

# Calculate GC content of each read in a fastq file, and output results to a tab-separated text file

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import argparse
import sys

def GC_content(read):
    seq = str(read.seq)

    GC = len([base for base in seq if base == 'G' or base == 'C'])

    GC_percent = float(GC)/len(seq) * 100

    return GC_percent


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help="input fastq file")
    args = parser.parse_args()

    if args.fastq.endswith('.gz'):
        fastq = gzip.open(args.fastq)
    else:
        fastq = args.fastq

    fq_reads = SeqIO.parse(fastq, 'fastq')

    for read in fq_reads:
        GC = GC_content(read)
        GC_line = '%s\t%.2f\n' % (read.id, GC)
        sys.stdout.write(GC_line)