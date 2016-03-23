__author__ = 'Guilhem Faure'
__description__ = '''Gather Genbank from the genome assembly id'''

import optparse
import sys
from Bio import Entrez
from Bio import SeqIO


class Genbank:
    def __init__(self, **kwargs):
        '''

        :param id:
        :return:
        '''

        # From an id we load from internet
        if kwargs.get('id'):
            self.id = kwargs.get('id')
            self._handle = self._download_gb()

        # If dump, we dump in the output
        if kwargs.get('dump'):
            self.dump_gb(kwargs.get('dump'))
            return None

        # Create a _gb variable
        self._load_gb()

        # Extract annotated 16S
        self.extract_16S()

    def _download_gb(self):
        '''
        Download genbank from the web
        :return:
        '''

        Entrez.email = 'xxx@gmail.com'
        handle = Entrez.efetch(db='nucleotide', id=options.genbank, rettype='gbwithparts', retmode="text")

        return handle


    def dump_gb(self, p_output):
        '''
        Write the Genbank file
        :param p_output: output path for genbank
        :return:
        '''

        with open(p_output, 'w') as fout:
            fout.write(self._handle.read())

        #self._handle.close()

        return

    def _load_gb(self):
        '''
        Load GB into memory
        :return:
        '''

        self._gb = SeqIO.read(self._handle, 'genbank')

    def extract_16S(self):
        '''

        :return:
        '''

        fasta_format = '>{type}|{genome}|' \
                       'position={start}-{stop}:{strand}' \
                       '|locus={locus}' \
                       '|gene={gene}' \
                       '|product={product}\n' \
                       '{seq}\n'
        is_16S = lambda x:True if sum([True for i in x if '16S ribosomal RNA' in i]) else False

        for gene in self._gb.features:

            if 'rRNA' in gene.type:
                if 'product' in gene.qualifiers:

                    # Look for 16S annotation
                    if is_16S(gene.qualifiers['product']):
                        d_info = {'genome':self.id+' '+self._gb.id, 'type':gene.type}

                        d_info['seq']  = gene.extract(self._gb.seq)
                        d_info['name'] = ' '.join(gene.qualifiers['product'])
                        d_info['gene'] = ','.join(gene.qualifiers['gene'])
                        d_info['locus']= ','.join(gene.qualifiers['locus_tag'])
                        d_info['start']     = gene.location.start.position
                        d_info['stop']      = gene.location.end.position
                        d_info['strand']    = gene.location.strand
                        d_info['product']   = ','.join(gene.qualifiers['product'])

                        print (fasta_format.format(**d_info))

        print (self._gb.annotations)
        print (self._gb.dbxrefs)
        print (self._gb.description)
        print (self._gb.features)
        print (self._gb.id)
        print (self._gb.name)




if __name__ == '__main__':

    usage = "usage: %prog [-h ] -g <Genome id assembly NC_xxxx> [-o output.gb]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--genbank', dest='genbank', help='Genbank file full size', default = None)
    parser.add_option('-o', '--output', dest='output', help='output file', default = None)
    (options, args) = parser.parse_args()


    if options.genbank is None:
        parser.error('You should provide a name assembly i.e. NC_000913 see -h')

    if options.genbank and options.output:
        MyGB = Genbank(id = options.genbank, dump = options.output)
        sys.exit('Genbank:'+options.genbank+' available at:'+options.output)

    if options.genbank:
        MyGB = Genbank(id = options.genbank)

