'''
Dependencies
http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software


Execute this line to install eutils unix command line
http://www.ncbi.nlm.nih.gov/books/NBK179288/

  cd ~
  perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
     $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");'
  unzip -u -q edirect.zip
  rm edirect.zip
  export PATH=$PATH:$HOME/edirect
  ./edirect/setup.sh


Aspera
http://downloads.asperasoft.com/connect2///


esearch -db sra -query SRX403935| efetch --format runinfo | cut -d ',' -f 1 | grep SRR |  xargs fastq-dump -X 10 --split-files

Remove cache for SRA toolkit
./vdb-config -i
uncross Cache

'''
import os
import configparser
import subprocess
import optparse



class SRA:
    def __init__(self):
        '''

        :return:
        '''

        return None

def download_sra(command, workdir):
    '''
    Download SRR file in workdir
    :param command: command to download all SRR from GSM number
    :param workdir: work directory
    :return:
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

    usage = "usage: %prog [-h] -g <GSEsample>"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--gse', dest='gse', help='GSE sample generate by GSE.py')
    (options, args) = parser.parse_args()

    if options.gse == None:
        parser.error("You should provide a GSEsample file with -g or --gse, -h for help")

    p_d_script = os.path.dirname(os.path.realpath(__file__))
    Config_dependancy = configparser.ConfigParser()
    Config_dependancy.read(os.path.join( p_d_script, "dependency.config"))

    Config_command = configparser.ConfigParser()
    Config_command.read(os.path.join( p_d_script, "command.config"))

    d_program = dict(dict(Config_dependancy['EUTILS']), **dict(Config_dependancy['SRATOOLKIT']))


    with open(options.gse) as f:
        for line in f:
            if line.startswith('#'):
                continue
            sp = line.strip().split('\t')
            sample_gsm = sp[0]

            command_download_sra = Config_command['SRA']['download_sra'].format(**dict(d_program, **{'gsm':sample_gsm}))

            l_srr = download_sra(command_download_sra, sample_gsm)
            print (l_srr)
            break

