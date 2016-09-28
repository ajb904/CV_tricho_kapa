#!/usr/local/bin/python

import subprocess
import os
from Bio import SeqIO
import argparse

# Cluster contigs from an assembly using CD-HIT-EST at a variety of different identity levels. Shorter contigs
# should match longer contigs over 100% of their length

CD_HIT_EST = '~/Downloads/cd-hit-v4.6.6-2016-0711/cd-hit-est'

def run_CDHIT(fasta, identity, word_size, outdir, cov=1.0):
    # identity should be given as a decimal between 0 and 1, up to 2 dp
    percent_id = float(identity)*100

    outfile = get_outfile(fasta, percent_id, outdir)

    cdhit_cline = [CD_HIT_EST,
                   '-i', fasta,
                   '-o', outfile,
                   '-c', '%.2f' % float(identity),
                   '-n', str(word_size),
                   '-aS', '1.0',
                   '-d', '0',
                   '-T', '0',
                   '-M', '32000',
                   '>', '%s.log' % outfile]

    cdhit_cline = ' '.join(cdhit_cline)

    print "Running: %s" % cdhit_cline

    #Only run clustering if result file doesn't exist
    if not os.path.exists(outfile):
        subprocess.call(cdhit_cline, shell=True)

    return outfile

def count_clusters(clustered_file, size_threshold):
    #Return number of clusters in a file that are greater than a specified size)

    fasta = SeqIO.parse(clustered_file, 'fasta')

    seq_count = 0
    for seq in fasta:
        if len(seq.seq) >= int(size_threshold):
            seq_count += 1

    return seq_count

def get_outfile(fasta, percent_id, outdir):
    infilename = os.path.basename(fasta)
    outfile = infilename.replace('.fasta', '%d.fasta' % percent_id)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, outfile)

    return outfile

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fasta", help="Fasta file containing sequences to be clustered", required=True)
    parser.add_argument("-t", "--tmpdir", help="Directory for clustering output", required=True)
    parser.add_argument("-o", "--outfile", help="Result file", required=True)

    args = parser.parse_args()

    percents = range(80, 101, 5)
    word_sizes = range(5, 11)

    input_fasta = args.fasta
    outdir = args.tmpdir

    result_file = open(args.outfile, 'w')

    thresholds_to_check = [0, 200, 500, 1000, 2000, 5000, 10000]

    result_file.write('Identity' + ('\tContigs.gt.%dbp' * len(thresholds_to_check)) % tuple(thresholds_to_check) + '\n')

    for p, ws in zip(percents, word_sizes):
        identity = float(p)/100

        clustered_file = run_CDHIT(input_fasta, identity, ws, outdir)

        cluster_counts = []
        for threshold in thresholds_to_check:
            cluster_count = count_clusters(clustered_file, threshold)
            cluster_counts.append(cluster_count)

        outline = str(p) + '\t%d' * len(cluster_counts) % tuple(cluster_counts) + '\n' 

        result_file.write(outline)

    result_file.close()
            
