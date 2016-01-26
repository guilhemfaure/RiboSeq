import os
import subprocess

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
        self.setup_dir

        self.get_xml()

    def setup_dir(self):
        '''
        Create workdir if non existing
        Workdir name is the GSE number
        :return:
        '''

        p_root = os.path.join(self.d_workdir, self.gse)
        if not os.path.isdir(p_root):
            os.makedirs(p_root)
            print('Create workdir', p_root)
        else:
            print('Workdir existing', p_root)

        os.chdir(p_root)

        return None

    def get_xml(self):
        '''

        :return:
        '''

        xml_url_template = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/miniml/{accession}_family.xml.tgz'
        xml_url = xml_url_template.format(prefix = self.gse[:5],
                                          accession = self.gse)

        wget_command = ['wget', '--quiet', '-P', paper_dir, xml_url]
        subprocess.check_call(wget_command)

if __name__ == '__main__':

    GSE('GSE67387', '../../RiboSeq_dev')