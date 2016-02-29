__author__ = 'Guilhem Faure'
__description__ = ''''''

import os
import sys



if __name__ == '__main__':

    f_fasta = sys.argv[1]
    f_blast = sys.argv[2]


    seq = ''
    with open(f_fasta) as f:
        for line in f:
            if line.startswith('>'):
                header = line
                continue
            seq += line.strip()

    print (seq)


