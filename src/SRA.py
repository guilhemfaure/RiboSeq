__author__ = 'Guilhem Faure'
__description__ = '''SRA.py download SRR files from GSM or SRX list  '''

import os
import configparser
import subprocess
import optparse


def download_sra(command, workdir):
    '''
    Download SRR file in workdir
    :param command: command to download all SRR from GSM number
    :param workdir: work directory
    :return: list of SRR files downloaded
    '''


    p_ori = os.getcwd()
    if not os.path.exists(sample_gsm):
        os.mkdir(sample_gsm)
    os.chdir(sample_gsm)

    print ('Downloading', sample_gsm)

    subprocess.run(command, check = True, shell = True)

    l_srr = os.listdir('.')
    os.chdir(p_ori)

    return l_srr

if __name__ == '__main__':

    usage = "usage: %prog [-h -w <workdir>] -g <GSEsample>"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--gse', dest='gse', help='GSE sample generate by GSE.py')
    parser.add_option('-w', '--workdir', dest='workdir', help='Workdir output', default = None)
    (options, args) = parser.parse_args()

    if options.gse == None:
        parser.error("You should provide a GSEsample file with -g or --gse, -h for help")



    ## Loading dependency and command lines
    p_d_script = os.path.dirname(os.path.realpath(__file__))
    Config_dependancy = configparser.ConfigParser()
    Config_dependancy.read(os.path.join( p_d_script, "dependency.config"))

    Config_command = configparser.ConfigParser()
    Config_command.read(os.path.join( p_d_script, "command.config"))

    d_program = dict(dict(Config_dependancy['EUTILS']), **dict(Config_dependancy['SRATOOLKIT']))

    ## Output path
    p_out = os.path.basename(options.gse).split('.')[0]+'.sra'

    # Read GSE.sample and download selected samples
    with open(options.gse) as f:

        # if workdir is enter by the user
        if options.workdir and not os.path.exists(options.workdir):
            os.mkdir(options.workdir)
            os.chdir(options.workdir)

        with open(p_out, 'w') as fout:
            fout.write('#GSM and list of SRA files associated\n')

            for line in f:
                if line.startswith('#'):
                    continue
                sp = line.strip().split('\t')
                sample_gsm = sp[0]

                command_download_sra = Config_command['SRA']['download_sra'].format(**dict(d_program, **{'gsm':sample_gsm}))

                l_srr = download_sra(command_download_sra, sample_gsm)

                fout.write('{gsm}\t{nbsrr}\t{srr}\n'.format(gsm=sample_gsm,\
                                                          nbsrr = len(l_srr),
                                                          srr = ' '.join(l_srr)))

    print ( 'Downloaded SRA files', p_out   )


