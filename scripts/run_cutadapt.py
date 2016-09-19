#!/usr/local/bin/python

#Usage: python run_cutadapt.py <in_R1> <in_R2> <out_R1> <out_R2> <log>

import sys
import subprocess

QUALITY_THRESHOLD = 15
MAX_TRIMMED_SIZE = 250
FORWARD_ADAPTOR = 'AGATCGGAAGAGC'
REVERSE_ADAPTOR = 'AGATCGGAAGAGC'

rawR1, rawR2, trimmedR1, trimmedR2, logfile = tuple(sys.argv[1:])

cline = "cutadapt -a %s -A %s -q %d -m %d -o %s -p %s %s %s" % (FORWARD_ADAPTOR,
                                                                REVERSE_ADAPTOR,
                                                                QUALITY_THRESHOLD,
                                                                MAX_TRIMMED_SIZE,
                                                                trimmedR1,
                                                                trimmedR2,
                                                                rawR1,
                                                                rawR2)

print '\tRunning cutadapt with command line: \n\t%s' % cline

proc = subprocess.Popen(cline, shell=True, stdout=subprocess.PIPE)

stdout = proc.communicate()[0]

log = open(logfile, 'w')
log.write(stdout)
log.close()

