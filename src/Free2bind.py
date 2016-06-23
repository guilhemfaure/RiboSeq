__author__ = 'Guilhem Faure'
__description__ = '''ExtractSeq.py extract sequence from different format file
  Already implemented Genbank and extraction of tRNA and rRNA sequence'''

import os
import configparser
import subprocess
import optparse
import sys
from Bio import SeqIO
import itertools



def grep_dg(res):

    for i in res.split('\n'):
        if i.startswith('Delta'):
            sp = i.split()
            dg1 = sp[-5]
            dg2 = sp[-3]
            dg3 = sp[-1]
            break
    return dg1, dg2, dg3


if __name__ == '__main__':

    nb_nt = 8
    l_nt = [''.join(i) for i in itertools.product('AUCG', repeat=nb_nt)]

    sd = 'ACCTCC'[::-1].replace('T', 'U')

    print (sd)

    command_free2bind = 'perl free_align.pl -f {seq} {sd}'

    print (os.getcwd())

    for ite, seq in enumerate(l_nt):
        res  = subprocess.getoutput(command_free2bind.format(seq=seq, sd=sd))

        print (sd, seq, ' '.join(grep_dg(res)))
