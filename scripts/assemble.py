#!/usr/local/bin/python

import argparse
import subprocess
import os
import sys

# Assembly wrapper script. Usage: assemble.py <method> <R1.fastq.gz> <R2.fastq.gz>

class AssemblerWrapper(object):

    def __init__(self, r1, r2, out_dir):
        self.fastq1 = r1
        self.fastq2 = r2
        self.out_dir = out_dir

    def prep(self):
        pass

    def run(self):
        # Take the assembly cline made by prep method, and run it
        print 'Assembling'
        print self.assembly_cline
        subprocess.call(self.assembly_cline, shell=True)
        print 'Done!'


class IDBAUDWrapper(AssemblerWrapper):

    def __init__(self, r1, r2, out_dir):
        super(IDBAUDWrapper, self).__init__(r1, r2, out_dir)

        self.basename = os.path.basename(self.fastq1).split('_R1')[0]
        self.interleaved_fastq = self.basename + "_interleaved.fastq.gz"
        self.fasta = self.basename + "_interleaved.fasta"

    def prep(self, update_fasta = False):
        # Check whether the interleaved fasta file exists. If not, create it.

        if not os.path.exists(self.fasta):
            fastq2fasta(self.fastq1, self.fastq2, self.fasta)

        self.assembly_cline = 'idba_ud -r %s -o %s --pre_correction' % (self.fasta, self.out_dir)


class SPAdesWrapper(AssemblerWrapper):

    def prep(self):
        self.assembly_cline = "spades.py -1 %s -2 %s -o %s --careful -t 8" % (self.fastq1,
                                                                              self.fastq2,
                                                                              self.out_dir)


class VelvetWrapper(AssemblerWrapper):

    def __init__(self, r1, r2, out_dir, mink = 31, maxk = 151, step = 10):
        super(VelvetWrapper, self).__init__(r1, r2, out_dir)

        self.mink = mink
        self.maxk = maxk
        self.step = step

        self.optimiser_exe = "~/Downloads/VelvetOptimiser-2.2.5/VelvetOptimiser.pl"

    def prep(self):

        velveth_string = "-fastq.gz -shortPaired -separate %s %s" % (self.fastq1, self.fastq2)

        self.assembly_cline = "%s -t 8 -s %d -e %d -x %d -f '%s' -p %s" % (self.optimiser_exe,
                                                                           self.mink,
                                                                           self.maxk,
                                                                           self.step,
                                                                           velveth_string,
                                                                           self.out_dir)


class MegahitWrapper(AssemblerWrapper):

    def prep(self):
        self.assembly_cline = 'megahit -1 %s -2 %s -o %s --presets meta' % (self.fastq1,
                                                                            self.fastq2,
                                                                            self.out_dir)



def zcat_ps(filename):
    # bash command to pipe a gzipped file in place of a regular fastq (using process substitution)
    if sys.platform == 'darwin':
        zcat = 'gzcat'
    else:
        zcat = 'zcat'

    return '<(%s %s)' % (zcat, filename)


def fastq2fasta(r1, r2, fasta):

    fq2fa_cline = 'fq2fa --merge %s %s %s' % ( zcat_ps( r1 ),
                                               zcat_ps( r2 ),
                                               fasta)
    print "merging reads with '%s'" % fq2fa_cline

    subprocess.call(fq2fa_cline, shell=True, executable='/bin/bash')


if __name__=='__main__':

    # Get assembly method and read files from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--method", help="Assembly method to use. Must be one of 'spades', 'megahit', 'velvet', 'idba_ud'", choices=['spades', 'megahit', 'velvet', 'idba_ud'], required=True)
        
    # For the moment, only deal with paired end reads. Require both read 1 and read 2
    # Need a better check here (are fastq files properly formatted?)
    parser.add_argument("-f", "--forward", help="Fastq file containing first reads of pair", required=True)
    parser.add_argument("-r", "--reverse", help="Fastq file containing second reads of pair", required=True)
        
    parser.add_argument("-o", "--out_dir", help="Assembly directory", default="assembly")
    parser.add_argument("-p", "--prefix", help="prefix to add to output directory (after method name)", default = "test")
    args = parser.parse_args()
        
        
    print args.method
    print args.forward, args.reverse
    run_dir = os.path.join(args.out_dir, "%s_%s" % (args.method, args.prefix))
    print run_dir

        
    if args.method == 'idba_ud':
        assembler = IDBAUDWrapper(args.forward, args.reverse, run_dir)

    elif args.method == 'spades':
        assembler = SPAdesWrapper(args.forward, args.reverse, run_dir)

    elif args.method == 'megahit':
        assembler = MegahitWrapper(args.forward, args.reverse, run_dir)

    elif args.method == 'velvet':
        assembler = VelvetWrapper(args.forward, args.reverse, run_dir)

    assembler.prep()
    assembler.run()
