import sys
import os
import tempfile
import argparse
import time
from subprocess import call

from Bio import SeqIO
from mpi4py import MPI

def countSeqs(fastafile):
    nbseq = 0
    fasta = SeqIO.parse(fastafile, 'fasta')
    for seq in fasta:
        nbseq+=1
    
    return nbseq

def splitfasta(fastafile, nbpart, outdir):
    '''split fasta file into nbpart chunks, and save split files in temporary directory'''

    total_seqs = countSeqs(fastafile)

    modulus = total_seqs % nbpart

    if modulus == 0:
        seqs_per_file = total_seqs//nbpart
    else:
        seqs_per_file = (total_seqs//nbpart) + 1

    fasta = SeqIO.parse(fastafile, 'fasta')
    fastaname = os.path.basename(fastafile)

    for part in range(nbpart):
        seqnum = 0
        seqlist = []
        outfile = '%s/%s_%d.fasta' % (outdir, fastaname, part+1)

        while seqnum < seqs_per_file:
            try:
                seqlist.append(fasta.next())
                seqnum += 1
            except StopIteration:
                break

        SeqIO.write(seqlist, outfile, 'fasta')

def replace_extension(fasta, out_ext, db):
    # Replace '.fa' or '.fasta' with required extension, and also include the database that was used.

    fasta_base = os.path.splitext(os.path.basename(fasta))[0]
    
    new_ext = '_v_%s.%s' % (db, out_ext)
    
    return fasta_base + new_ext


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nbpart = comm.Get_size()

if comm.rank == 0:
    print '#########################'
    print 'MPI wrapper for blast+'
    print 'This job will use %d cores' % nbpart
    print 'Blast run started at ' + time.strftime('%I:%M:%S %p', time.localtime())
    print '#########################\n'

    # User arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('-q', '--query', help='query sequences in fasta format', required=True)
    parser.add_argument('-o', '--out', help='result file', required=True)
    parser.add_argument('-x', '--blast_exe', help='blast program to run', required=True)
    parser.add_argument('-d', '--db', help='blast database to search', required=True)
    parser.add_argument('-b', '--blast_opt', help='string specifying extra blast options')
    parser.add_argument('-f', '--out_format', help='number representing blast output format (5 for xml, 6 for tab)', required=True)

    args=parser.parse_args()

    fastafile = args.query
    outfile = args.out
    outdir = os.path.dirname(outfile)
    blast_exe = args.blast_exe
    blastdb = args.db
    if args.blast_opt:
        blastOpt = args.blast_opt
    else:
        blastOpt = ''

    outfmt = int(args.out_format)
    if outfmt == 5:
        out_ext = 'xml'
    elif outfmt == 6:
        out_ext = 'tab'
    else:
        out_ext = 'out'

    print 'outfmt = ' + str(outfmt)

    #create temp directory
    tmpdir = tempfile.mkdtemp(prefix=os.path.basename(blast_exe), dir=outdir)
    print '#########################'
    print 'writing temporary files into ' + tmpdir
    print '#########################\n'

    splitfasta(fastafile, nbpart, tmpdir)

    print 'Temporary files written (' + time.strftime('%I:%M:%S %p', time.localtime()) + ')\n'

    blast_files = os.listdir(tmpdir)
    out_files = [replace_extension(f, out_ext, os.path.basename(blastdb)) for f in blast_files]

else:
    fastafile = None
    blastdb = None
    blast_exe = None
    blastOpt = None
    tmpdir = None
    blast_files = None
    out_files = None
    outfmt = None


c_blast_exe = comm.bcast(blast_exe, root=0)
c_blastOpt = comm.bcast(blastOpt, root=0)
c_tmpdir = comm.bcast(tmpdir, root=0)
c_blast_files = comm.bcast(blast_files, root=0)
c_out_files = comm.bcast(out_files, root=0)
c_blastdb = comm.bcast(blastdb, root=0)
c_outfmt = comm.bcast(outfmt, root=0)

##### for tests

print'\n'
print 'prog ',c_blast_exe , ' and my rank is ', rank
print 'opt ',c_blastOpt, ' and my rank is ', rank
print 'tmp ',c_tmpdir, ' and my rank is ', rank
print 'query ',c_blast_files[rank], ' and my rank is ', rank
print 'out ', c_out_files[rank], ' and my rank is ', rank

#Run blast

cmd = [c_blast_exe,
       '-db', c_blastdb,
       '-query', os.path.join(c_tmpdir, c_blast_files[rank]),
       '-out', os.path.join(c_tmpdir, c_out_files[rank]),
       '-outfmt', str(c_outfmt),
       c_blastOpt
       ]

cline = ' '.join(cmd)

print '\n'
print 'Command is:\n', cline, '\nRank: ', rank, '\n'

blast_run = call( cline, shell=True)

#Wait for all processes
comm.Barrier()

#Finish up after all processes return
if comm.rank == 0:
    concat = 'cat %s/*.%s > %s' % (tmpdir, out_ext, outfile) 
    call(concat, shell=True)

    print 'Done! Blast results in %s' % outfile
    print 'End time: ' + time.strftime('%I:%M:%S %p', time.localtime())