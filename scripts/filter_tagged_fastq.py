#!/usr/local/bin/python

# Filter fastq files that have been tagged by the Fastq_screen program to:
# a) keep reads that don't align to any contaminant genomes, while
# b) preserving read pairing

import argparse
from Bio import SeqIO
import gzip
import os
import subprocess

def read_fastq_gz(gzfile):
    handle = gzip.open(gzfile, 'rb')

    fastq = SeqIO.parse(handle, 'fastq')

    return fastq

def get_tag(read):
    # Fastq_screen tag is appended to end of read description. Should be of the form '#FQST:<tag>', except for
    # the first in the file, where the genomes used for filtering are included (e.g. #FQST:PhiX:Pfluor:00)

    try:
        fq_screen_tag = read.description.split('#FQST')[1]
        filter_tag = fq_screen_tag.split(':')[-1]

        return filter_tag

    except IndexError:
        print 'Read not correctly tagged'
        raise


def is_not_contaminant(read):
    # Fastq_screen tags use '0' to denote that a read did not align to a particular contaminant genome
    # Therefore, a read that doesn't align to any contaminants will have a tag that evaluates to 0 when
    # converted to an integer, regardless of how many genomes are used for filtering

    tag = get_tag(read)

    if int(tag) == 0:
        return True

    else:
        return False


def tag_reads(fastq, out_dir):
    # Run fastq_screen --tag on each input file and return the filename of the tagged file


    fastq_screen_cline = 'fastq_screen --tag --outdir %s %s' % (out_dir, fastq)
    print 'Running: %s' % fastq_screen_cline
    subprocess.call(fastq_screen_cline, shell=True)

    tagged_fastq_file = os.path.basename(fastq).replace('.fastq.gz', '.tagged.fastq.gz')

    tagged_fastq_filepath = os.path.join(out_dir, tagged_fastq_file)

    return tagged_fastq_filepath



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--fastq1', help="forward reads", required=True)
    parser.add_argument('-2', '--fastq2', help="reverse reads", required=True)
    parser.add_argument('-t', '--tagged', help="are reads already tagged?", action='store_true')
    parser.add_argument('-o', '--out_dir', help="output directory")
    args = parser.parse_args()

    fq_file1 = args.fastq1
    fq_file2 = args.fastq2

    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_dir = os.path.dirname(fq_file1)

    # If reads are not already tagged, run fastq_screen to tag them. Intermediate 'tagged' files will have
    # '.tagged' inserted into the file name before '.fastq.gz'

    if not args.tagged:
        fq_file1 = tag_reads(fq_file1, out_dir)
        fq_file2 = tag_reads(fq_file2, out_dir)

    fq1 = read_fastq_gz(fq_file1)
    fq2 = read_fastq_gz(fq_file2)

    keep_r1 = []
    keep_r2 = []

    read_counter = 0

    while 1:
        try:
            r1 = fq1.next()
            r2 = fq2.next()

            assert r1.id == r2.id, 'Fastq files must be correctly paired'

            if is_not_contaminant(r1) and is_not_contaminant(r2):
                keep_r1.append(r1)
                keep_r2.append(r2)

            read_counter += 1

        except StopIteration:
            break


    keep = len(keep_r1)
    total = read_counter
    keep_percent = float(keep)/total * 100

    print 'Done filtering. Keeping %d out of %d read pairs (%.2f %%)' % (len(keep_r1),
                                                                         read_counter,
                                                                         keep_percent)

    outfile1 = fq_file1.replace('.tagged', '.filtered')
    outfile2 = fq_file2.replace('.tagged', '.filtered')

    out1 = gzip.open(outfile1, mode='wb')
    out2 = gzip.open(outfile2, mode='wb')

    SeqIO.write(keep_r1, out1, 'fastq')
    SeqIO.write(keep_r2, out2, 'fastq')

    out1.close()
    out2.close()