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
import atg.data.ontology

DATA_SOURCE_PATH = os.path.join(os.path.dirname(__file__), 'data_sources.ini')
GENOME_FILES = [
    'chrom.sizes',                      # UCSC chromosome sizes
    'genome.fa',                        # Ensembl genome sequence
    'genes.gtf',                        # Ensembl gene annotation
    'gene_go.csv',                      # Ensembl gene <-> GO accession
    'go_definition.csv',                # information for all GO terms
    'ensembl_gene.csv',                 # Ensembl gene information
    'ensembl_gene_homology_human.csv',  # Homology information with Human
    'ensembl_gene_transcript.csv',      # Ensembl gene <-> transcript ID
    'ensembl_transcript.csv'            # Ensembl transcript information
    # TODO: get APPRIS transcript annotation from BioMart (canonical transcript info)
    # TODO: get NCBI ID (Entrez Gene ID) for translation
]
BINARY_CHUNKSIZE = 2**30


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


def fetch_url(url, path, overwrite=False, verbose=False):
    """
    Download specified url to destination path, returning True if successful. Will not overwrite existing files by
    default.
    :param url:
    :param path:
    :param overwrite:
    :param verbose:
    :return:
    """
    if not overwrite and os.path.exists(path):
        if verbose:
            print('The file %s already exists, so %s was not downloaded.' % (path, url))
        return False

    # skip empty URLs, e.g. human homologs for human genome data
    if url == '':
        pass

    # for compressed files, make sure to decompress before writing to disk
    elif url.endswith('.gz'):
        response = urllib.request.urlopen(url)
        with open(path, 'wb') as output_file:
            if sys.platform == 'darwin':
                # write output in chunks to avoid bug in MacOS [https://bugs.python.org/issue24658]
                binary_result = gzip.decompress(response.read())
                result_length = len(binary_result)
                chunk_start = 0

                while chunk_start < result_length:
                    chunk_end = min(result_length, chunk_start+BINARY_CHUNKSIZE)
                    output_file.write(binary_result[chunk_start:chunk_end])
                    chunk_start = chunk_end
            else:
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

    def __init__(self, data_root=None):
        if data_root:
            self.data_root = data_root
        else:
            self.data_root = os.path.expanduser(atg.config.settings['Data']['Root'])

        self.config = configparser.ConfigParser(default_section="Common",
                                                interpolation=configparser.ExtendedInterpolation())
        self.config.read(DATA_SOURCE_PATH)
        self.genome_path = {}

        for organism in self.config.sections():
            ensembl_genome = self.config.get(organism, 'ensembl_genome')
            current_genome_path = os.path.join(self.data_root, organism, ensembl_genome)
            self.genome_path[organism] = current_genome_path
            os.makedirs(current_genome_path, exist_ok=True)

    def check_annotation(self, organism):
        return organism in self.genome_path

    def retrieve_data(self, genome, overwrite=False, selected_files=GENOME_FILES):
        """
        Fetch all or selected data files for a specified genome

        :param genome:
        :param overwrite:
        :param selected_files: a restricted list of files to retrieve
        :return:
        """

        if not self.check_annotation(genome):
            print('%s is not an available genome.' % genome, file=sys.stderr)
            return False

        for filename, current_url in self.config.items(genome):
            if filename not in selected_files:
                continue

            current_path = os.path.join(self.genome_path[genome], filename)
            if current_url.startswith('<?xml'):
                fetch_ensembl(current_url, current_path)
            else:
                fetch_url(current_url, current_path, overwrite=overwrite)

    def derive_data(self, genome, overwrite=False):
        """
        Generate additional data files from raw retrieved data, e.g. GO biological process info.

        :param genome:
        :param overwrite:
        :return:
        """

        if not self.check_annotation(genome):
            print('%s is not an available genome.' % genome, file=sys.stderr)
            return False

        go_term_filename = os.path.join(self.genome_path[genome], "gene_go.csv")
        go_biological_process_filename = os.path.join(self.genome_path[genome], "go_biological_process.csv")

        if not overwrite and os.path.exists(go_biological_process_filename):
            return True

        go_biological_process_df = atg.data.ontology.process_ontology(go_term_filename)
        go_biological_process_df.to_csv(go_biological_process_filename, index=False)


def retrieve_data(namespace):
    tracker = ATGDataTracker()

    if not tracker.check_annotation(namespace.organism):
        print('%s genome information is not available.\n' % namespace.organism, file=sys.stderr)
        print("The following organisms are available:", file=sys.stderr)
        for organism in sorted(tracker.genome_path.keys()):
            print("\t%s" % organism, file=sys.stderr)
        print('', file=sys.stderr)
    else:
        tracker.retrieve_data(namespace.organism, namespace.overwrite)
        tracker.derive_data(namespace.organism, namespace.overwrite)


def setup_subparsers(subparsers):
    retrieval_parser = subparsers.add_parser('ensembl', help="Use Ensembl's main site")
    retrieval_parser.add_argument('organism', help="a species common name, e.g. human or mouse")
    retrieval_parser.add_argument('-o', '--overwrite', action="store_true", help="Overwrite existing files")
    retrieval_parser.set_defaults(func=retrieve_data)
