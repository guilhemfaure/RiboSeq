import os
import subprocess

import xml.etree.ElementTree as etree
import string

import gzip
from gzip import GzipFile
from urllib.request import urlopen
import sys
import base64
import tarfile
import string


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

        # Download xml
        self.path_xml = self.get_xml()

        # Get sample name
        self.sample = self.get_sample_from_xml()



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

    def get_sample_from_xml(self, condition=lambda x: True):
        '''Adapted from JHussmann
        Parse an xml file describing a GEO accession number to extract samples.

        Sanitizes the names of samples by replacing spaces, slashes, and brackets
        with underscores, replacing unicode deltas with the string 'delta_', and
        removing parentheses.

        condition -- a filtering function to be applied to sample names

        Returns a list of (sample name, list of URLs) tuples.
        '''


        samples = []

        tree = etree.parse(self.gse+'.xml')
        root = tree.getroot()
        for child in root:

            if child.tag.endswith('Sample'):

                for grand in child:
                    if grand.tag.endswith('Title'):
                        sample_name = grand.text
                        sample_URLs = []
                        if not type(sample_name) == str:
                            # Some sample names (knockout strains) have a unicode
                            # delta in them. Replace this with a spelled out
                            # 'delta' and coerce to a non-unicode string.
                            deltas = [u'\u0394', u'\u2206']
                            for delta in deltas:
                                sample_name = sample_name.replace(delta, 'delta_')
                            sample_name = str(sample_name)

                    elif grand.tag.endswith('Supplementary-Data') and grand.attrib['type'] == 'SRA Experiment':
                        sample_URL = grand.text.strip()
                        sample_URLs.append(sample_URL)

                samples.append((sample_name, sample_URLs))

        print('Number of sample', len(samples))
        return samples



if __name__ == '__main__':

    GSE('GSE67387', '../../../RiboSeq_dev')