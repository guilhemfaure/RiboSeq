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

        self._gb = SeqIO.parse(self._handle, 'genbank')

    def extract_16S(self):
        '''

        :return:
        '''

        print (self._gb)
        for gene in self._gb:
            print (dir(gene))
            break






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

