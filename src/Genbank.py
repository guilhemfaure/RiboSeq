__author__ = 'Guilhem Faure'
__description__ = '''Managing Genbank from the genome assembly id
Biopython and python3 are required

<Script mode>
Get genbank file
python3 Genbank.py -g NC_000913.3 [-o NC_000913.3.gb]

Get 16s
python3 Genbank.py -g NC_000913.3 -f 16s.fasta -m

Get Exon
python3 Genbank.py -g NC_000913.3 -f exon.fasta -e

Get tRNA and rRNA
python3 Genbank.py -g NC_000913.3 -f exon.fasta -r

from a genbank file:
python3 Genbank.py -i NC_000913.3.gb -f exon.fasta -r

<Class mode>
myGB = Genbank(id = 'NC_000913.3')
or
myGB = Genbank(input = 'NC_000913.3.gb')
myGB.   extract(sequence = ['CDS', 'rRNA', 'tRNA', 'r16s'])
        save_genbank('output.gb')
        _gb object

'''

import optparse
import sys
from Bio import Entrez
from Bio import SeqIO


class Genbank:
    def __init__(self, **kwargs):
        '''

        :param kwargs: id -> chromosome version (required)
                       input -> genbank file

                        <Script Mode>
                       extract -> list of sequence type to extract ['CDS', 'rRNA', 'tRNA', 'r16s']
                       fasta -> output file to record the sequence
                       output -> output file to record the genbank file
        :return:
        '''

        # From an id we load from internet
        if kwargs.get('id'):
            self.id = kwargs.get('id')
            self._handle = self._download_gb()
            self._load_gb()

        # From a local file
        if kwargs.get('input'):
            self._handle = kwargs.get('input')
            self._load_gb()
            self.id = self._gb.id


        ## Script MODE
        if kwargs.get('fasta'):
            self.extract(sequence = kwargs.get('extract'), output = kwargs.get('fasta'))

        if kwargs.get('output'):
            self.save_genbank(kwargs.get('output'))


        ## Class mode
        # Genbank(id = 'NC_000913.3') or
        # Genbank(input = 'NC_000913.3.gbfile')
        # Create a _gb variable
        #self._load_gb()

        # Extract annotated sequence
        #self.extract(sequence = ['CDS', 'tRNA', 'rRNA', 'r16s'], output = 'file')

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
        <Script Mode> Obsolete ==> But normally faster than the current save_genbank,
        I switched to save_genbank since it is more conveninent and I case used loading data
        Once this function is launched, we lost data from memory
        Write the Genbank file from the script mode
        :param p_output: output path for genbank
        :return:
        '''

        with open(p_output, 'w') as fout:
            fout.write(self._handle.read())

        print ('Genbank available at: ', p_output)


        return

    def _load_gb(self):
        '''
        <Class mode or Script mode>
        Load GB into memory
        :return:
        '''

        self._gb = SeqIO.read(self._handle, 'genbank')

        return


    def extract(self, sequence = [], output = None):
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

        print (sequence, output)

        for gene in self._gb.features:

            if gene.type in sequence:
                if 'product' in gene.qualifiers:


                    d_info = {'genome':self._gb.id, 'type':gene.type}

                    d_info['seq']  = gene.extract(self._gb.seq)
                    d_info['name'] = ' '.join(gene.qualifiers['product'])
                    d_info['gene'] = ','.join(gene.qualifiers['gene'])
                    d_info['locus']= ','.join(gene.qualifiers['locus_tag'])
                    d_info['start']     = gene.location.start.position
                    d_info['stop']      = gene.location.end.position
                    d_info['strand']    = gene.location.strand
                    d_info['product']   = ','.join(gene.qualifiers['product'])

                    str_fasta = fasta_format.format(**d_info)

                    if 'r16s' in sequence and is_16S(gene.qualifiers['product']) is False:
                        continue

                    header, seq, _ = str_fasta.split('\n')
                    d_seq[header] = seq

                    if output:
                        fout.write(str_fasta)

        if output:
            fout.close()

        return d_seq



    def save_genbank(self, output):
        '''
        <Class mode or Script mode>
        Save genbank into a file
        :return:
        '''


        with open(output, "w") as output:
            SeqIO.write(self._gb, output, "genbank")

        return


if __name__ == '__main__':

    usage = "usage: %prog [-h ] -g <Genome id assembly NC_xxxx> or -i gb.file\n" \
            "[-o output.gb] download GB file \n" \
            "[-f] extract sequence use -e (for exon) -m (for 16S) -r (for tRNA and rRNA)\n"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--genbank', dest='genbank', help='Input genome assembly name', default = None)
    parser.add_option('-i', '--input', dest='input', help='Input Genbank file', default = None)
    parser.add_option('-o', '--output', dest='output', help='dump Genbankfile', default = None)
    parser.add_option('-f', '--fasta', dest='fasta', help='dump fasta sequence', default = None)

    parser.add_option('-e', '--exon', dest='exon', action='store_true', help='extract exon sequence')
    parser.add_option('-m', '--r16s', dest='r16s', action='store_true', help='extract 16S')
    parser.add_option('-r', '--rRNA', dest='rRNA', action='store_true', help='extract rRNA and tRNA')

    (options, args) = parser.parse_args()

    if options.genbank is None and options.input is None:
        parser.error('You should provide an assembly with -g NC_000913.3 or a genbank file with -i gb.file see -h')

    l_extract = []
    if options.exon:
        l_extract.append('CDS')
    if options.r16s:
        l_extract.append('r16s')
        l_extract.append('rRNA')
    if options.rRNA:
        l_extract.append('rRNA')
        l_extract.append('tRNA')

    MyGB = Genbank(id = options.genbank,
                   extract = l_extract,
                   fasta = options.fasta,
                   output = options.output,
                   input = options.input)
    #sys.exit('Genbank:'+options.genbank+' sequences dumped:'+' '.join(l_extract)+' available at:'+options.fasta)
