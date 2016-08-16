"""
Assemble genomic data from various sources (e.g. UCSC, Ensembl) and organize them in a standard way.
"""

import os
import sys
import urllib.parse
import urllib.request
import urllib.error
import gzip
import configparser
import atg.config
import atg.data

DATA_SOURCE_PATH = os.path.join(os.path.dirname(__file__), 'data_sources.ini')
GENOME_FILES = [
    'chrom.sizes',                     # UCSC chromosome sizes
    'genome.fa',                       # Ensembl genome sequence
    'genes.gtf',                       # Ensembl gene annotation
    'gene_go.csv',                     # Ensembl gene <-> GO accession
    'go_definition.csv',               # information for all GO terms
    'ensembl_gene.csv',                # Ensembl gene information
    'ensembl_gene_transcript.csv',     # Ensembl gene <-> transcript ID
    'ensembl_transcript.csv'           # Ensembl transcript information
]


def fetch_ensembl(xml_string, output_filename):
    """
    retrieve data from Ensembl using XML specifications (e.g. for gene annotation or Gene Ontology sets)

    :param xml_string:
    :param output_filename:
    :return:
    """
    query = urllib.parse.quote(xml_string)
    url = 'http://www.ensembl.org/biomart/martservice?query=' + query
    return fetch_url(url, output_filename)


def fetch_url(url, path, overwrite=False):
    """
    Download specified url to destination path, returning True if successful. Will not overwrite existing files by
    default.
    :param url:
    :param path:
    :param overwrite:
    :return:
    """
    if not overwrite and os.path.exists(path):
        return False

    # for compressed files, make sure to decompress before writing to disk
    if url.endswith('.gz'):
        response = urllib.request.urlopen(url)
        with open(path, 'wb') as output_file:
            output_file.write(gzip.decompress(response.read()))

    # download uncompressed files directly
    else:
        try:
            urllib.request.urlretrieve(url, filename=path)

        except urllib.error.HTTPError as error:
            print("Could not retrieve %s: %s [HTTP Error %s]" % (url, error.reason, error.code), file=sys.stderr)
            return False

    return True


class ATGDataTracker:
    """

    """

    def __init__(self):
        self.data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
        self.config = configparser.ConfigParser(default_section="Common",
                                                interpolation=configparser.ExtendedInterpolation())
        self.config.read(DATA_SOURCE_PATH)
        self.genome_path = {}

        for genome_id in self.config.sections():
            organism = self.config.get(genome_id, 'organism')
            current_genome_path = os.path.join(self.data_root, organism, 'Current', genome_id)
            self.genome_path[genome_id] = current_genome_path
            os.makedirs(current_genome_path, exist_ok=True)

    def retrieve_data(self, genome, overwrite=False):
        """
        Fetch all data files for a specified genome

        :param genome:
        :param overwrite:
        :return:
        """

        if genome not in self.genome_path:
            print('%s is not an available genome.' % genome, file=sys.stderr)
            return False

        for filename, current_url in self.config.items(genome):
            if filename not in GENOME_FILES:
                continue

            current_path = os.path.join(self.genome_path[genome], filename)
            if current_url.startswith('<?xml'):
                fetch_ensembl(current_url, current_path)
            else:
                fetch_url(current_url, current_path, overwrite=overwrite)


if __name__ == '__main__':
    tracker = ATGDataTracker()
    tracker.retrieve_data('hg38', overwrite=False)
    tracker.retrieve_data('mm10', overwrite=False)
