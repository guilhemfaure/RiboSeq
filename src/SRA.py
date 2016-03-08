__author__ = 'Guilhem Faure'
__description__ = '''SRA.py download SRR files from GSM or SRX list  '''

import os
import configparser
import subprocess
import optparse
import copy
import sys


def download_sra(command):
    '''
    Download SRR file in workdir
    :param command: command to download all SRR from GSM number
    :return: stdout of download_sra
    '''


    print ('Downloading', sample_gsm)
    log_download_sra = subprocess.getoutput(command)

    return log_download_sra

def gsm2srr(command):
    '''

    :param command_gsm2srr: command to grab the SRR from the GSM
    :return: srr_id
    '''

    print ('Getting the SRR identification')
    srr_id = subprocess.getoutput(command)

    return srr_id

def clip_adapter(command):
    '''
    Clip adapter
    :param command: command to clip adapter from gzip fastq
    :return: log clip adapter
    '''


    print ('Clipping adapter', sample_gsm)

    log_clip_adapter = subprocess.getoutput(command)

    return log_clip_adapter


if __name__ == '__main__':

    usage = "usage: %prog [-h -w <workdir>] -g <GSEsample>"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--gse', dest='gse', help='GSE sample generate by GSE.py')
    parser.add_option('-w', '--workdir', dest='workdir', help='Workdir output', default = None)
    parser.add_option('-a', '--adapter', dest='adapter', help='adapter', default = 'CTGTAGGCACCATCAAT')
    parser.add_option('-c', '--clip', dest='clip', action='store_true',
                      help='clip adapter; use -a to specify adapter (default CTGTAGGCACCATCAAT)')
    (options, args) = parser.parse_args()

    if options.gse == None:
        parser.error("You should provide a GSEsample file with -g or --gse, -h for help")


    ## Loading dependency and command lines
    p_d_script = os.path.dirname(os.path.realpath(__file__))
    Config_dependancy = configparser.ConfigParser()
    Config_dependancy.read(os.path.join( p_d_script, "dependency.config"))

    Config_command = configparser.ConfigParser()
    Config_command.read(os.path.join( p_d_script, "command.config"))


    d_options = copy.copy(Config_dependancy['EUTILS'])
    d_options.update(Config_dependancy['SRATOOLKIT'])



    ## Output path
    p_out = os.path.basename(options.gse).split('.')[0]+'.log'


    p_root = os.getcwd()

    # Read GSE.sample and download selected samples
    with open(options.gse) as f:

        # if workdir is enter by the user
        if options.workdir and not os.path.exists(options.workdir):
            os.mkdir(options.workdir)
            os.chdir(options.workdir)



        # Read sample to download
        for line in f:
            if line.startswith('#'):
                continue
            sp = line.strip().split('\t')
            sample_gsm = sp[0]
            d_options['gsm'] = sample_gsm


            # Create GSM directory
            if not os.path.exists(sample_gsm):
                os.mkdir(sample_gsm)
            os.chdir(sample_gsm)

            # Log file in each GSM directory
            with open(sample_gsm, 'w') as fout:
                fout.write('#Log '+ ' '.join(sys.argv)+'\n')

                # Get SRR number
                command_gsm2srr = Config_command['SRA']['gsm2srr'].format(**d_options)
                fout.write(command_gsm2srr+'\n')
                srr_id = gsm2srr(command_gsm2srr)
                d_options['srr'] = srr_id
                fout.write('GSM2SRR {0} {1}\n'.format(sample_gsm, srr_id))


                # Download fastq.gz
                command_download_sra = Config_command['SRA']['download_sra'].format(**d_options)
                fout.write(command_download_sra+'\n')
                log_download_sra = download_sra(command_download_sra)
                fout.write(log_download_sra+'\n')

                if options.clip is True:
                    d_options.update(Config_dependancy['FASTX'])
                    d_options['min_read_size'] = '25'
                    d_options['adapter'] = options.adapter
                    d_options['input_fastqgz'] = srr_id+'.fastq.gz'
                    d_options['output_clip'] = srr_id+'_clip.fastq.gz'
                    command_clip_adapter = Config_command['FASTX']['clip_adapter'].format(**d_options)
                    fout.write(command_clip_adapter+'\n')
                    log_clip_adapter = clip_adapter(command_clip_adapter)
                    fout.write(log_clip_adapter+'\n')

            # Coming back to root directory
            os.chdir(p_root)


    print ( 'Log file', p_root )




