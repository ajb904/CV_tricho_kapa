#!/usr/local/bin/python

import argparse
import subprocess
import os
import sys



# Assembly wrapper script. Usage: assemble.py <method> <R1.fastq.gz> <R2.fastq.gz>

def zcat_ps(filename):
    # bash command to pipe a gzipped file in place of a regular fastq (using process substitution)
    if sys.platform == 'darwin':
        zcat = 'gzcat'
    else:
        zcat = 'zcat'

    return '<(%s %s)' % (zcat, filename)

def paste_fq_read(filename):
    # bash paste command for joining together four lines of fastq
    return 'paste - - - - < %s' % filename

def interleave_reads_cline(r1, r2, outfile):
    # Make a bash command line for interleaving paired end fastq files - it's much faster than any
    # python solution
    interleave_cline = "paste <(%s) <(%s) | tr '\\t' '\\n' | gzip > %s" % (paste_fq_read( zcat_ps( r1 ) ),
                                                                           paste_fq_read( zcat_ps( r2 ) ),
                                                                           outfile)
    return interleave_cline

def fastq2fasta_cline(fastq, fasta):
    return 'fq2fa %s %s' % (fastq, fasta)

def Velvet_cline(r1, r2, mink = 31, maxk = 151, step = 20):
    pass

def Megahit_cline():
    pass

def IDBA_UD_cline(fq1, fq2, out_dir):
    # IDBA requires paired end reads in a single, interleaved fasta file. Check if this exists and create if not.
    # By default, fasta file will have same stem as individual read files, and will be stored in main assembly
    # directory, one level above the run directory.

    prefix = os.path.basename(fq1).split('_R1')[0]
    interleaved_fastq = prefix + '_interleaved.fastq.gz'
    interleaved_fasta = prefix + '_interleaved.fasta'

    #Generate interleaved fastq file
    interleave_cline = interleave_reads_cline(fq1, fq2, interleaved_fastq)
    print 'Running "%s"' % interleave_cline
    subprocess.call( interleave_cline, shell=True, executable='/bin/bash' )

    #Convert fastq to fasta
    fq2fa_cline = fastq2fasta_cline(interleaved_fastq, interleaved_fasta)
    print fq2fa_cline

    #Make IDBA_UD command line
    idba_ud_cline = 'idba_ud -r %s -o %s --pre_correction' % (interleaved_fasta, out_dir)

    print idba_ud_cline


def SPAdes_cline():
        pass

def run_assembly(cline):
        pass


if __name__=='__main__':

    # Get assembly method and read files from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--method", help="Assembly method to use. Must be one of 'spades', 'megahit', 'velvet', 'idba_ud'", choices=['spades', 'megahit', 'velvet', 'idba_ud'], required=True)
        
    # For the moment, only deal with paired end reads. Require both read 1 and read 2
    # Need a better check here (are fastq files properly formatted?)
    parser.add_argument("-f", "--forward", help="Fastq file containing first reads of pair", required=True)
    parser.add_argument("-r", "--reverse", help="Fastq file containing second reads of pair", required=True)
        
    parser.add_argument("-o", "--out_prefix", help="Prefix to add to run output folder (after method name)", default="")
    args = parser.parse_args()
        
        
    print args.method
    print args.forward, args.reverse
    run_dir = "%s_%s" % (args.method, args.out_prefix)
    print run_dir

        
    if args.method == 'idba_ud':
        print IDBA_UD_cline(args.forward, args.reverse, run_dir)