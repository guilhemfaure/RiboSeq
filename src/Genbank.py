__author__ = 'Guilhem Faure'
__description__ = '''Managing Genbank from the genome assembly id
Biopython and python3 are required

<Script mode>
Get genbank file
python3 Genbank.py -g NC_000913.3 [-o NC_000913.3.gb]

Get 16s
python3 Genbank.py -g NC_000913.3 -m 16s.fasta

<Class mode>
myGB = Genbank(id = 'NC_000913.3')
myGB.   extract_16s('output.file')
        save_genbank('output.gb')

'''

import optparse
import sys
from Bio import Entrez
from Bio import SeqIO


class Genbank:
    def __init__(self, **kwargs):
        '''

        :param kwargs: id -> chromosome version (required)

                        <Script Mode>
                       dump -> dump genbank file in script mode
                       rRNAS -> dump 16s fasta sequence in script mode
        :return:
        '''

        # From an id we load from internet
        if kwargs.get('id'):
            self.id = kwargs.get('id')
            self._handle = self._download_gb()


        ## Script MODE
        # If dump, we dump in the output
        if kwargs.get('dump'):
            self.dump_gb(kwargs.get('dump'))
            return None

        # Extract 16S into a file
        if kwargs.get('rRNA16s'):
            self._load_gb()
            self.extract_16S(kwargs.get('rRNA16s'))
            return None


        ## Class mode
        # Create a _gb variable
        self._load_gb()

        # Extract annotated 16S
        #self.extract_16S(output = None)

        # Save genbank
        #self.save_genbank(output = 'output.genbank')

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
        <Script Mode>
        Write the Genbank file from the script mode
        :param p_output: output path for genbank
        :return:
        '''

        with open(p_output, 'w') as fout:
            fout.write(self._handle.read())


        print ('Genbank available at: ', p_output)
        #self._handle.close()

        return

    def _load_gb(self):
        '''
        Load GB into memory
        :return:
        '''

        self._gb = SeqIO.read(self._handle, 'genbank')

        return


    def extract_16S(self, output = None):
        '''
        <Class mode or Script mode>
        Extract 16S fasta from the genbank annotation
        if output is specified it dumps the sequences into a file
        :return: d_seq [header] -> sequence
        '''

        if output:
            fout = open(output, 'w')
        d_seq = {}
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

                        str_fasta = fasta_format.format(**d_info)

                        header, seq, _ = str_fasta.split('\n')
                        d_seq[header] = seq

                        if output:
                            fout.write(str_fasta)

        if output:
            fout.close()

        return d_seq
        # print (self._gb.annotations)
        # print (self._gb.dbxrefs)
        # print (self._gb.description)
        # print (self._gb.features)
        # print (self._gb.id)
        # print (self._gb.name)

    def save_genbank(self, output):
        '''
        <Class mode>
        Save genbank into a file
        :return:
        '''


        with open("short_seqs.fasta", "w") as output:
            SeqIO.write(self._gb, output, "genbank")

        return


if __name__ == '__main__':

    usage = "usage: %prog [-h ] -g <Genome id assembly NC_xxxx> [-o output.gb]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--genbank', dest='genbank', help='Genbank file full size', default = None)
    parser.add_option('-o', '--output', dest='output', help='output file', default = None)
    parser.add_option('-m', '--r16s', dest='rRNA16s', help='dump 16S fasta', default = None)
    (options, args) = parser.parse_args()


    if options.genbank is None:
        parser.error('You should provide a name assembly i.e. NC_000913.3 see -h')

    if options.genbank and options.output:
        MyGB = Genbank(id = options.genbank, dump = options.output)
        sys.exit('Genbank:'+options.genbank+' available at:'+options.output)

    if options.genbank and options.rRNA16s:
        MyGB = Genbank(id = options.genbank, rRNA16s = options.rRNA16s)
        sys.exit('Genbank:'+options.genbank+' available at:'+options.rRNA16s)

    if options.genbank:
        MyGB = Genbank(id = options.genbank, dump = options.genbank+'.gb')

