__author__ = 'Guilhem Faure'
__description__ = '''ExtractSeq.py extract sequence from different format file
  Already implemented Genbank and extraction of tRNA and rRNA sequence'''

import os
import sys
import configparser
import subprocess
import optparse

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

    genome_gbk = SeqIO.read(p_genbank,'genbank')

    read_bam = pysam.AlignmentFile(f_bam, "rb")


    for ite_gene, gene in enumerate(genome_gbk.features):

        if gene.type in ['CDS']:

            d_info = {'type' : gene.type, 'genome':genome_gbk.id}
            d_info['start']     = gene.location.start.position
            d_info['stop']      = gene.location.end.position
            d_info['strand']    = gene.location.strand
            d_info['seq']       = gene.extract(genome_gbk.seq)

            print ('{seq} {strand}'.format(**d_info))

            print ( d_info['start'],  d_info['stop'])
            for ite, pilup in enumerate(read_bam.pileup('gi|49175990|ref|NC_000913.2|',
                                                          d_info['start'],
                                                          d_info['stop'])):
                print ('#####################')
                print ( pilup.pos, pilup.n)

                if ite > 50 :
                    break

                continue

            break
                # for i in pilup.pileups:
                #     print (i.alignment.query_sequence)
                #     print (i.alignment.query_alignment_start)
                #     print (i.alignment.query_alignment_end)
                #     print (i.alignment.query_length)
                #
                #     print (i.alignment.reference_start)
                #     print (i.alignment.reference_end)






    # plus = 0
    # minus = 0
    # for ite, i in enumerate(file.fetch()):
    #     seq_length = len(i.seq)
    #     seq_3p_position = i.positions[-1]
    #
    #     strand = '-' if i.is_reverse else '+'
    #
    #     print(seq_length, seq_3p_position, strand)
    #     break



