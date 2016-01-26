__author__ = 'Guilhem Faure'
__description__ = '''GSE.py seeks for RNA-seq and Ribo-seq experiment form a GSE number '''


import os
import xml.etree.ElementTree as etree
from urllib.request import urlopen
import tarfile
import optparse
import sys


class GSE:

    def __init__(self, gse, workdir = os.getcwd()):
        '''

        :param gse: GSE identification
        :param workdir: path to work directory
        :return:
        '''

        # Get variables
        self.gse = gse
        self.workdir = os.getcwd() if workdir == None else workdir
        self.path_xml = self.gse+'.xml'

        # Create Workdir
        if not os.path.exists(self.workdir):
            self.setup_dir()

        # Download xml
        if not os.path.exists(self.path_xml):
            self.get_xml()

        # Get sample name
        self.l_sample = self.get_sample_from_xml()

        # Stats and sample dumping
        self.dump_display()

        return None

    def dump_display(self):
        '''
        Display number of sample per species
        Dump the data into a .sample file to edit
        :return:
        '''

        d_species = {}
        with open(self.gse+'.sample', 'w') as fout:
            fout.write('# Comment the sample you do not want to use\n')

            for sample_name,\
                sample_id, \
                sample_species, \
                sample_genome_assembly, \
                sample_URL\
                    in self.l_sample:
                if sample_species not in d_species:
                    d_species[sample_species] = [0,sample_genome_assembly]
                d_species[sample_species][0] += 1
                fout.write('\t'.join((sample_species, sample_name, sample_id, sample_genome_assembly, ' '.join(sample_URL)))+'\n')

        for species in d_species:
            print(species, d_species[species][0], d_species[species][1])

        print('Total number of sample', len(self.l_sample))
        print('Edit your sample file, comment the samples', self.gse+'.sample')


    def setup_dir(self):
        '''
        Create workdir if non existing
        Workdir name is the GSE number
        :return:
        '''

        path_root = os.path.join(self.workdir, self.gse)
        if not os.path.isdir(path_root):
            os.makedirs(path_root)
            print('Create workdir', path_root)
        else:
            print('Workdir existing', path_root)

        os.chdir(path_root)
        print('Moving to', os.getcwd())

        return None

    def get_xml(self):
        '''
        Download xml file of the experiment
        :return: path of the xml file in gz
        '''

        xml_url_template = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/miniml/{accession}_family.xml.tgz'
        xml_url = xml_url_template.format(prefix = self.gse[:5],
                                          accession = self.gse)

        f_xml = urlopen(xml_url, os.path.basename(xml_url))

        # Dump the file
        with open(self.gse+'.xml.gz', 'wb') as f: f.write(f_xml.read())
        print('Download', self.gse+'.xml')

        # Unzip and untar the file
        with tarfile.open(self.gse+'.xml.gz') as f: f.extractall()
        print('Untar', self.gse+'.xml')

        # Rename the xml file
        os.rename(self.gse+'_family.xml', self.gse+'.xml')

        # Remove the tar file
        os.remove(self.gse+'.xml.gz')

        return self.gse+'.xml'

    def get_sample_from_xml(self):
        '''Adapted from JHussmann
        Get data from th xml
        Returns a list of (sample_name, sample_id, sample_species, sample_genome_assembly, sample_URL) tuples.
        '''


        l_sample = []

        tree = etree.parse(self.gse+'.xml')
        root = tree.getroot()
        for child in root:

            if child.tag.endswith('Sample'):
                sample_id = child.attrib['iid']
                sample_URL = []
                sample_name = None
                sample_species = None
                sample_genome_assembly = None
                for grand in child:

                    if grand.tag.endswith('Title'):
                        sample_name = grand.text

                        if not type(sample_name) == str:
                            # Some sample names (knockout strains) have a unicode
                            # delta in them. Replace this with a spelled out
                            # 'delta' and coerce to a non-unicode string.
                            deltas = [u'\u0394', u'\u2206']
                            for delta in deltas:
                                sample_name = sample_name.replace(delta, 'delta_')
                            sample_name = str(sample_name)


                    elif grand.tag.endswith('Supplementary-Data') and grand.attrib['type'] == 'SRA Experiment':
                        sample_URL.append(grand.text.strip())

                    # Grab the species
                    elif grand.tag.endswith('Channel'):
                        for great in grand:
                            if great.tag.endswith('Organism'):
                                sample_species = great.text.strip()

                    elif grand.tag.endswith('Data-Processing'):
                        for line in grand.text.split('\n'):
                            if line.startswith('Genome_build'):
                                sample_genome_assembly =  line.split('Genome_build:')[1].strip()

                l_sample.append((sample_name, sample_id, sample_species, sample_genome_assembly, sample_URL))


        return l_sample



if __name__ == '__main__':

    usage = "usage: %prog [-h -w <workdir>] -g <GSEnumber>"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-g', '--gse', dest='gse', help='GSE identification')
    parser.add_option('-w', '--workdir', dest='workdir', help='Workdir, where .xml and .sample will be generated', default=None)
    (options, args) = parser.parse_args()



    if options.gse == None:
        parser.error("You should provide a GSE number with -g or --gse, -h for help")

    GSE(options.gse, options.workdir)