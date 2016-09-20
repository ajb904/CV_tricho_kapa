#!/usr/local/bin/python

import argparse
import os
import subprocess


class Bowtie2Wrapper(object):

    def __init__(self, fastq1, fastq2, index, aligned_dir, unaligned_dir=False):

        self.fastq1 = fastq1
        self.fastq2 = fastq2

        self.sample = os.path.basename(self.fastq1).split('_R1')[0]

        self.index = index
        self.reference = os.path.basename(index)

        self.out_prefix = '%s_v_%s' % (self.sample, self.reference)

        self.aligned_dir = aligned_dir

        if unaligned_dir:
            self.unaligned_dir = unaligned_dir


    def prep(self, aligned_only=True, keep_unaligned=True, bam_convert=True, sort=True):

        self.bam_convert = bam_convert
        self.sort = sort

        #Check if index exists. Create if not

        if not os.path.exists(self.index + '.1.bt2'):
            bt2_build_cline = 'bowtie2-build %s %s' % (self.index + '.fasta', self.index)
            print bt2_build_cline

            subprocess.call(bt2_build_cline, shell=True, executable='/bin/bash')

        self.bt2_cline = 'bowtie2 -tp 8 -X 1000 -x %s -1 %s -2 %s' % (self.index,
                                                                  self.fastq1,
                                                                  self.fastq2)

        if aligned_only:
            self.bt2_cline += ' --no-unal'

        if keep_unaligned:
            unaligned_file = os.path.join(self.unaligned_dir,
                                          '%s_unaligned_R%%.fastq.gz' % self.out_prefix)
            self.bt2_cline += ' --un-conc-gz %s' % unaligned_file

        if not self.sort and not self.bam_convert:
            samfile = os.path.join(self.aligned_dir, '%s.sam' % self.out_prefix)
            self.bt2_cline += ' -S %s' % samfile

        elif self.bam_convert:
            bamfile = os.path.join(self.aligned_dir, self.out_prefix + '.bam')
            self.bt2_cline += ' | samtools view -uS'

            if not self.sort:
                self.bt2_cline += ' -o %s' % bamfile

            else:
                bamfile = os.path.join(self.aligned_dir, self.out_prefix)
                self.bt2_cline += ' - | samtools sort - %s' % bamfile


    def run(self):

        print self.bt2_cline
        subprocess.call(self.bt2_cline, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--forward", help="Fastq file containing first reads of pair", required=True)
    parser.add_argument("-r", "--reverse", help="Fastq file containing second reads of pair", required=True)
    parser.add_argument("-i", "--index", help="Bowtie2 index for reference sequence", required=True)
    parser.add_argument('-a', "--align_dir", help="directory for aligned reads", required=True)
    parser.add_argument('-u', "--unalign_dir", help="directory for unaligned reads")

    args = parser.parse_args()

    if not os.path.exists(args.align_dir):
        os.makedirs(args.align_dir)

    if args.unalign_dir:
        if not os.path.exists(args.unalign_dir):
            os.makedirs(args.unalign_dir)

    bt2 = Bowtie2Wrapper(args.forward, args.reverse, args.index, args.align_dir, args.unalign_dir)
    bt2.prep()
    bt2.run()