#!/usr/local/bin/python

import subprocess
import os

# Cluster contigs from an assembly using CD-HIT-EST at a variety of different identity levels. Shorter contigs
# should match longer contigs over 100% of their length

CD_HIT_EST = '~/Downloads/cd-hit-v4.6.6-2016-0711/cd-hit-est'

def run_CDHIT(fasta, identity, outdir, cov=1.0):
    # identity should be given as a decimal between 0 and 1, up to 2 dp
    percent_id = float(identity)*100
    outfile_suffix = ('%.0f' % percent_id).replace('.','_') + '.fasta'
    outfile = fasta.replace('.fasta', outfile_suffix)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, outfile)

    cdhit_cline = [CD_HIT_EST,
                   '-i', fasta,
                   '-o', outfile,
                   '-c', '%.1f' % float(identity),
                   '-aS', '1.0',
                   '-d', '0',
                   '-T', '0',
                   '-M', '32000']

    cdhit_cline = ' '.join(cdhit_cline)
    print cdhit_cline

    proc = subprocess.Popen(cdhit_cline, shell=True, stdout=subprocess.PIPE)

    stdout = proc.communicate()[0]

    num_clusters=0
    for line in stdout:
        print line
        if line.endswith('clusters'):
            num_clusters = line.split[3]

    return outfile, num_clusters

if __name__=='__main__':
    percents = range(80, 101, 5)

    for p in percents:
        identity = float(p)/100
        s = run_CDHIT('temp.fasta', identity, 'test')
        print s

