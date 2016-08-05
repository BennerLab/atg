"""
Assemble genomic data from various sources (e.g. UCSC, ENSEMBL) and organize them in a standard way.
"""

import os
import sys
import urllib.parse
import urllib.request
import urllib.error
import atg.config
import atg.data


def fetch_ensembl(xml_filename, output_filename):
    """
    retrieve data from Ensembl using XML specifications (e.g. for gene annotation or Gene Ontology sets)

    :param xml_filename:
    :param output_filename:
    :return:
    """
    input_file = open(xml_filename)
    definition = ''.join(input_file.readlines())
    query = urllib.parse.quote(definition)
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
        self.data_root = atg.config.settings['Data']['Root']
        self.genome_path = {}
        for genome, organism in atg.data.genome_hierarchy.items():
            current_genome_path = os.path.join(self.data_root, organism, 'Current', genome)
            self.genome_path[genome] = current_genome_path
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

        for filename, base_url in atg.data.genome_data.items():
            current_url = base_url.format(genome=genome)
            current_path = os.path.join(self.genome_path[genome], filename)
            fetch_status = fetch_url(current_url, current_path, overwrite=overwrite)


if __name__ == '__main__':
    tracker = ATGDataTracker()
    tracker.retrieve_data('hg38')
