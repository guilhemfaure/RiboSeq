__author__ = 'Guilhem Faure'
__description__ = '''Ribosome Occupancy reads bam index and extract the 3 prime position
of each well aligned reads'''

import os
import sys
import configparser
import subprocess
import optparse
import numpy as np


from Bio import SeqIO
import pysam


if __name__ == '__main__':

    usage = "usage: %prog -b file.bam [-g <Genbank file> or -f <Genome file>] -o root.wig"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-b', '--bam', dest='p_bam', help='Bam file', default = None)
    parser.add_option('-f', '--fasta', dest='p_fasta', help='Fasta file', default = None)
    parser.add_option('-o', '--output', dest='p_output', help='Root of wig files', default = None)
    (options, args) = parser.parse_args()


    f_bam = options.p_bam


    if options.p_fasta:
        genome_fasta = SeqIO.read(open(options.p_fasta), "fasta")
        gi = genome_fasta.id
        size_genome = len(genome_fasta.seq)


    p_output = options.p_output

    sam_file = pysam.AlignmentFile(f_bam, "rb")



    with open(p_output, 'w') as fout:
        fout.write('#command {0}\n # read.id start stop strand\n'.format(' '.join(sys.argv)))
        # REF start is 0 in BAM
        # gi|16127994|ref|NC_000913.1|
        # gi|49175990|ref|NC_000913.2|
        # gi|556503834|ref|NC_000913.3|
        for read in sam_file.fetch('gi|556503834|ref|NC_000913.3|', 0, size_genome):
            if read.get_tag('NM') != 0:
                continue

            if read.query_alignment_length < 25:
                continue

            l_pos = read.get_reference_positions()

            strand = '-' if read.is_reverse else '+'

            fout.write('{0}\t{1}\t{2}\t{3}\n'.format(read.query_name, l_pos[0], l_pos[-1], strand))

