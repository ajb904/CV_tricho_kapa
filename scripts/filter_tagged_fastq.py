#!/usr/local/bin/python

# Filter fastq files that have been tagged by the Fastq_screen program to:
# a) keep reads that don't align to any contaminant genomes, while
# b) preserving read pairing

import argparse