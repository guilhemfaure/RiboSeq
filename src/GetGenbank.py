__author__ = 'Guilhem Faure'
__description__ = '''Gather Genbank from the genome assembly id'''

import optparse
from Bio import Entrez


if __name__ == '__main__':

    usage = "usage: %prog [-h ] -g <Genome id assembly NC_xxxx> [-o output.gb]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--genbank', dest='genbank', help='Genbank file full size', default = None)
    parser.add_option('-o', '--output', dest='p_output', help='output file', default = None)
    (options, args) = parser.parse_args()


    if options.genbank is None:
        parser.error('You should provide a name assembly i.e. NC_000913 see -h')


    Entrez.email = 'xxx@gmail.com'
    handle = Entrez.efetch(db='nucleotide', id=options.genbank, rettype='gbfull', retmod='text')


    if options.p_output is None:
        print (handle.read())

    else:
        with open(options.p_output, 'w') as fout:
            fout.write(handle.read())

