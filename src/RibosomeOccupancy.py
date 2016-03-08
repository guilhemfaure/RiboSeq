__author__ = 'Guilhem Faure'
__description__ = '''ExtractSeq.py extract sequence from different format file
  Already implemented Genbank and extraction of tRNA and rRNA sequence'''

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

    f_bam = sys.argv[1]
    p_genbank = sys.argv[2]
    p_output = sys.argv[3]

    genome_gbk = SeqIO.read(p_genbank,'genbank')
    size_genome = len(genome_gbk.seq)

    sam_file = pysam.AlignmentFile(f_bam, "rb")

    print (size_genome)

    d_rf = {'minus': [0]*size_genome,
            'plus': [0]*size_genome }



    # REF start is 0 in BAM
    for read in sam_file.fetch('gi|49175990|ref|NC_000913.2|', 0, size_genome):
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
    print ('Nb Read:', nb_read)
    with open(p_output+'_minus.wig', 'w') as fout:
        fout.write('Nb read {0}\n3 prime position \n'.format(nb_read))
        fout.write('\n'.join( map(str, list( np.array(d_rf['minus']) / float(nb_read) ) )))

    with open(p_output+'_plus.wig', 'w') as fout:
        fout.write('Nb read {0}\n3 prime position \n'.format(nb_read))
        fout.write('\n'.join( map(str, list( np.array(d_rf['plus']) / float(nb_read) ) )))


    print (size_genome)


