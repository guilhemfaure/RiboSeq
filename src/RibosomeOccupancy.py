__author__ = 'Guilhem Faure'
__description__ = '''Ribosome Occupancy reads bam index and extract the 3 prime position
of each well aligned reads'''

import os
import sys
import configparser
import subprocess
import optparse
import numpy as np

# try:
#     import pysam
# except:
#     from VirtualEnvOnDemand import enableOnDemandImporter, ensureImportGlobal
#     enableOnDemandImporter()   # Activate the hook on the "import" keyword.
#     pysam = ensureImportGlobal('pysam', 'pysam')
#     Bio = ensureImportGlobal('Bio', 'biopython')
#     import Bio
#     import pysam

from Bio import SeqIO
import pysam


if __name__ == '__main__':

    usage = "usage: %prog -b file.bam [-g <Genbank file> or -f <Genome file>] -o root.wig"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-b', '--bam', dest='p_bam', help='Bam file', default = None)
    parser.add_option('-g', '--genbank', dest='p_genbank', help='Genbank file', default = None)
    parser.add_option('-f', '--fasta', dest='p_fasta', help='Fasta file', default = None)
    parser.add_option('-o', '--output', dest='p_output', help='Root of wig files', default = None)
    (options, args) = parser.parse_args()


    f_bam = options.p_bam

    if options.p_genbank:
        genome_gbk = SeqIO.read(options.p_genbank,'genbank')
        size_genome = len(genome_gbk.seq)
        # Find how to get id from genbank with Biopython
    if options.p_fasta:
        genome_fasta = SeqIO.read(open(options.p_fasta), "fasta")
        gi = genome_fasta.id
        size_genome = len(genome_fasta.seq)


    p_output = options.p_output

    sam_file = pysam.AlignmentFile(f_bam, "rb")
    sam_file = pysam.AlignmentFile(f_bam, "rb")

    print (size_genome)

    d_rf = {'minus': [0]*size_genome,
            'plus': [0]*size_genome }



    # REF start is 0 in BAM
    # gi|16127994|ref|NC_000913.1|
    # gi|49175990|ref|NC_000913.2|
    # gi|556503834|ref|NC_000913.3|
    for read in sam_file.fetch('gi|255767013|ref|NC_000964.3|', 0, size_genome):
        if read.get_tag('NM') != 0:
            continue

        if read.query_alignment_length < 25:
            continue

        l_pos = read.get_reference_positions()

        # Minus strand
        if read.is_reverse:
            position_3p = l_pos[0]
            d_rf['minus'][position_3p] += 1

        # Plus strand
        else:
            position_3p = l_pos[-1]
            d_rf['plus'][position_3p] += 1


    nb_read = sum(d_rf['minus']) + sum(d_rf['plus'])
    m = 1000000 # to get RPM
    print ('Nb Read:', nb_read)
    with open(p_output+'_minus.wig', 'w') as fout:
        fout.write('Nb read {0}\n3 prime position \n'.format(nb_read))
        fout.write('\n'.join( map(str, list( np.array(d_rf['minus']) * m / float(nb_read) ) )))

    with open(p_output+'_plus.wig', 'w') as fout:
        fout.write('Nb read {0}\n3 prime position \n'.format(nb_read))
        fout.write('\n'.join( map(str, list( np.array(d_rf['plus']) * m / float(nb_read) ) )))


    print (size_genome)


