__author__ = 'Guilhem Faure'
__description__ = '''ExtractSeq.py extract sequence from different format file
  Already implemented Genbank and extraction of tRNA and rRNA sequence'''

import os
import configparser
import subprocess
import optparse

try:
    from Bio import SeqIO
except:
    from VirtualEnvOnDemand import enableOnDemandImporter, ensureImportGlobal
    enableOnDemandImporter()   # Activate the hook on the "import" keyword.
    Bio = ensureImportGlobal('Bio', 'biopython')
    from Bio import SeqIO


def genebank_extract_tRNA_rRNA(p_genbank, p_output = None):
    '''
    Extract tRNA and rRNA from genbank file
    :param: p_genbank path to the genbank file
    :param: p_output path to the fasta output file containing tRNA and rRNA

    :return:
    '''

    if p_output == None:
        p_output = os.path.basename(p_genbank)+'.tRNA_rRNA'


    genome=SeqIO.read(p_genbank,'genbank')

    fasta_format = '>{type}|{genome}|position={start}-{stop}:{strand}|locus={locus}|gene={gene}|product={product}\n{seq}\n'

    fout = open(p_output, 'w')
    nb_sequence = 0
    for gene in genome.features:
        if gene.type in ['tRNA', 'rRNA']:

            d_info = {'type' : gene.type, 'genome':genome.id}

            d_info['seq']       = gene.extract(genome.seq)
            d_info['start']     = gene.location.start.position
            d_info['stop']      = gene.location.end.position
            d_info['strand']    = gene.location.strand

            d_info['product']   = ','.join(gene.qualifiers['product'])
            d_info['gene']      = ','.join(gene.qualifiers['gene'])
            d_info['locus']     = ','.join(gene.qualifiers['locus_tag'])

            fout.write(fasta_format.format(**d_info))

            nb_sequence += 1

    fout.close()
    print('%d tRNA and rRNA have been extracted in %s'%(nb_sequence ,p_output))

    return None

def genebank_extract_exon(p_genbank, p_output = None):
    '''
    Extract exons from genbank file
    :param: p_genbank path to the genbank file
    :param: p_output path to the fasta output file containing tRNA and rRNA

    :return:
    '''

    if p_output == None:
        p_output = os.path.basename(p_genbank)+'.exon'


    genome=SeqIO.read(p_genbank,'genbank')

    fasta_format = '>{type}|{genome}|position={start}-{stop}:{strand}|locus={locus}|gene={gene}|product={product}\n{seq}\n'

    fout = open(p_output, 'w')
    nb_sequence = 0
    for gene in genome.features:
        if gene.type in ['CDS']:
            d_info = {'type' : gene.type, 'genome':genome.id}

            d_info['seq']       = gene.extract(genome.seq)
            d_info['start']     = gene.location.start.position
            d_info['stop']      = gene.location.end.position
            d_info['strand']    = gene.location.strand

            # some gene, like pseudo gene, transposon do not have product
            # We do not take pseudogene
            if 'product' not in gene.qualifiers:
                if 'pseudogene' not in ''.join(gene.qualifiers['note']):
                    d_info['product'] = ''
                else:
                    continue # we do not take pseudogene

            else:
                d_info['product']  = ','.join(gene.qualifiers['product'])



            d_info['gene']      = ','.join(gene.qualifiers['gene'])
            d_info['locus']     = ','.join(gene.qualifiers['locus_tag'])

            fout.write(fasta_format.format(**d_info))

            nb_sequence += 1

    fout.close()
    print('%d exons have been extracted in %s'%(nb_sequence ,p_output))

    return None




if __name__ == '__main__':

    usage = "usage: %prog [-h -w <workdir>] -g <Genbank file> [-r tRNArRNA] [-e exon]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--genbank', dest='p_genbank', help='Genbank file full size', default = None)
    parser.add_option('-r', '--rna', dest='p_rna', help='Extract tRNA and rRNA into a fasta file', default = None)
    parser.add_option('-e', '--exon', dest='p_exon', help='Extract exons into a fasta file', default = None)
    (options, args) = parser.parse_args()

    if options.p_genbank:

        if options.p_rna:
            genebank_extract_tRNA_rRNA(options.p_genbank, options.p_rna)
        if options.p_exon:
            genebank_extract_exon(options.p_genbank, options.p_exon)
    else:
        parser.error("You should provide a genbank file with -g or --genbank, -h for help")





