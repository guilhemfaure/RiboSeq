import os
import subprocess
import urllib2

class GSE:

    def __init__(self, gse, d_workdir):
        '''

        :param gse:
        :return:
        '''

        # Get variables
        self.gse = gse
        self.d_workdir = d_workdir

        # Create Workdir
        self.setup_dir()

        self.path_xml = self.get_xml()

    def setup_dir(self):
        '''
        Create workdir if non existing
        Workdir name is the GSE number
        :return:
        '''

        path_root = os.path.join(self.d_workdir, self.gse)
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

        f_xml = urllib2.urlopen(xml_url, os.path.basename(xml_url))
        with open(self.gse+'.xml.gz', 'wb') as f: f.write(f_xml.read())
        print('Download', self.gse+'.xml.gz')

        return self.gse+'.xml.gz'

if __name__ == '__main__':

    GSE('GSE67387', '../../../RiboSeq_dev')