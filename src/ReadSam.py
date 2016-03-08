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

    genome_gbk = SeqIO.read(p_genbank,'genbank')

    read_bam = pysam.AlignmentFile(f_bam, "rb")



    # CDS coords
    for ite_gene, gene in enumerate(genome_gbk.features):

        if gene.type in ['CDS']:

            d_info = {'type' : gene.type, 'genome':genome_gbk.id}
            d_info['start']     = gene.location.start.position
            d_info['stop']      = gene.location.end.position
            d_info['strand']    = gene.location.strand
            d_info['seq']       = gene.extract(genome_gbk.seq)
            d_info['gene']      = gene.qualifiers['gene']
            d_info['locus']     = ','.join(gene.qualifiers['locus_tag'])

            print ('{seq} {strand}'.format(**d_info))


            # Coord CDS
            print ( d_info['start'],  d_info['stop'])

            # Plus strand
            if d_info['strand'] != 1:
                print (d_info['strand'])
                continue


            for ite, pileupcolumn in enumerate(read_bam.pileup('gi|49175990|ref|NC_000913.2|',
                                                          d_info['start'],
                                                          d_info['stop'])):

                for pileupread in pileupcolumn.pileups:
                    if pileupread.alignment.get_tag('NM') != 0:
                        continue



                    start = int(pileupread.alignment.reference_start)
                    stop = int(pileupread.alignment.reference_end)

                    position_3p_on_gene = stop-d_info['start']-1



            #print (l_read / sum(l_read))
            for position, n_read in enumerate(l_read / sum(l_read)):
                f_out.write ('{gene} {position} {read}\n'.
                       format(gene = d_info['locus'],
                              position = position,
                              read = n_read))
    f_out.close()