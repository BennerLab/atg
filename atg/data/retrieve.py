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
import pandas
import progress.bar
import atg.config
import atg.data
import atg.data.ontology

ENSEMBL_PERMANENT_SITE = 'Apr2019.archive.ensembl.org'
DATA_SOURCE_PATH = os.path.join(os.path.dirname(__file__), 'data_sources.ini')
GENOME_FILES = [
    'chrom.sizes',                      # UCSC chromosome sizes
    'genome.fa',                        # Ensembl genome sequence
    'genes.gtf',                        # Ensembl gene annotation
    'gene_go.csv',                      # Ensembl gene <-> GO accession
    'go_definition.csv',                # information for all GO terms
    'ensembl_gene.csv',                 # Ensembl gene information
    'ensembl_gene_homology_human.csv',  # Homology information with Human
    'ensembl_transcript.csv'            # Ensembl transcript information
    # TODO: get APPRIS transcript annotation from BioMart (canonical transcript info)
    # TODO: get NCBI ID (Entrez Gene ID) for translation
]
BINARY_CHUNKSIZE = 2**30


def fetch_ensembl(xml_string, output_filename, overwrite=False):
    """
    retrieve data from Ensembl using XML specifications (e.g. for gene annotation or Gene Ontology sets)

    :param xml_string:
    :param output_filename:
    :param overwrite:
    :return:
    """
    query = urllib.parse.quote(xml_string)
    url = f'http://{ENSEMBL_PERMANENT_SITE}/biomart/martservice?query={query}'
    return fetch_url(url, output_filename, overwrite)


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
        return True

    # create empty file for blank URLs, e.g. human homologs for human genome data
    if url == '':
        open(path, 'a').close()

    # for compressed files, make sure to decompress before writing to disk
    elif url.endswith('.gz'):
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.HTTPError as error:
            print("Could not retrieve %s: %s [HTTP Error %s]" % (url, error.reason, error.code), file=sys.stderr)
            return False
        except urllib.error.URLError as error:
            print("Could not retrieve %s: %s" % (url, error.reason), file=sys.stderr)
            return False

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
        except urllib.error.URLError as error:
            print("Could not retrieve %s: %s" % (url, error.reason), file=sys.stderr)
            return False

    return True


class ATGDataTracker:
    """

    """

    def __init__(self, data_root=None):
        """

        :param data_root:
        """
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

            if os.path.exists(current_genome_path):
                self.genome_path[organism] = current_genome_path

    def check_available(self, organism):
        return self.config.has_section(organism)

    def report_available(self, *args):
        """
        Parse the config file for available organisms and print a table to stdout.
        :return: None
        """
        records = []
        section_list = self.config.sections()
        for organism in section_list:
            section = self.config[organism]
            records.append(list(map(section.get, ('genus_species', 'ensembl_genome', 'release'))))

        available_df = pandas.DataFrame(records, index=section_list,
                                        columns=['Scientific name', 'Assembly', 'Ensembl release'])
        print(available_df.rename_axis('Organism').reset_index().to_string(index=False))

    def list_installed(self, *args):
        """
        Look through ATG data directory for installed organisms. Print a table to stdout.

        :return:
        """
        records = []
        section_list = self.config.sections()
        # check if each organism has all downloaded files and get Ensembl release number
        for organism in section_list:
            installation_status = 'OK'
            # check for directory existing
            if organism in self.genome_path:
                current_genome_path = self.genome_path[organism]
                # check presence of individual files
                for filename in GENOME_FILES:
                    if not os.path.exists(os.path.join(current_genome_path, filename)):
                        installation_status = 'broken'
                # read Ensembl release if present
                genome_config_path = os.path.join(current_genome_path, 'info.ini')
                try:
                    genome_config = configparser.ConfigParser()
                    genome_config.read(genome_config_path)
                    release = genome_config.get('Ensembl', 'release')
                except:
                    release = 'unknown'

                records.append([organism, installation_status, release])

        installed_df = pandas.DataFrame(records, columns=['Organism', 'Installation status', 'Ensembl release'])

        print(installed_df.to_string(index=False))

    def retrieve_data(self, organism, overwrite=False, selected_files=GENOME_FILES):
        """
        Fetch all or selected data files for a specified genome. Record Ensembl version of downloaded files.

        :param organism:
        :param overwrite:
        :param selected_files: a restricted list of files to retrieve
        :return: True if files were successfully retrieved.
        """

        if not self.check_available(organism):
            print('%s is not an available genome.' % organism, file=sys.stderr)
            return False

        # create directory if necessary
        ensembl_genome = self.config.get(organism, 'ensembl_genome')
        current_genome_path = os.path.join(self.data_root, organism, ensembl_genome)
        os.makedirs(current_genome_path, exist_ok=True)

        # set up a progress bar
        progress_bar = progress.bar.Bar('Downloading', max=len(selected_files))
        progress_bar.suffix = '%(index)d/%(max)d %(filename)s'

        # fetch each file
        for filename, current_url in self.config.items(organism):
            if filename not in selected_files:
                continue

            progress_bar.filename = filename
            progress_bar.next()

            current_path = os.path.join(current_genome_path, filename)
            if current_url.startswith('<?xml'):
                successful_retrieval = fetch_ensembl(current_url, current_path, overwrite=overwrite)
            else:
                successful_retrieval = fetch_url(current_url, current_path, overwrite=overwrite)

            if not successful_retrieval:
                print('\nDid not successfully retrieve data for %s.\nFailed to get %s.' %
                      (organism, filename), file=sys.stderr)
                return False

        progress_bar.filename = ''
        progress_bar.update()
        progress_bar.finish()

        self.genome_path[organism] = current_genome_path

        # write a config file with Ensembl release number
        genome_config_path = os.path.join(current_genome_path, 'info.ini')
        genome_config = configparser.ConfigParser()
        genome_config.read(genome_config_path)
        if not genome_config.has_section('Ensembl'):
            genome_config.add_section('Ensembl')
        genome_config['Ensembl']['Release'] = self.config.get(organism, 'release')
        with open(genome_config_path, 'w') as genome_config_output:
            genome_config.write(genome_config_output)

        return True

    def derive_data(self, organism, overwrite=False):
        """
        Generate additional data files from raw retrieved data, e.g. GO biological process info.

        :param organism:
        :param overwrite:
        :return:
        """

        if not self.check_available(organism):
            print('%s is not an available genome.' % organism, file=sys.stderr)
            return False

        go_term_filename = os.path.join(self.genome_path[organism], "gene_go.csv")
        go_biological_process_filename = os.path.join(self.genome_path[organism], "go_biological_process.csv")

        if not overwrite and os.path.exists(go_biological_process_filename):
            return True

        go_biological_process_df = atg.data.ontology.process_ontology(go_term_filename)
        go_biological_process_df.to_csv(go_biological_process_filename, index=False)


def retrieve_data(namespace):
    tracker = ATGDataTracker()

    if not tracker.check_available(namespace.organism):
        print('%s genome information is not available.\n' % namespace.organism, file=sys.stderr)
        print("The following organisms are available:", file=sys.stderr)
        for organism in sorted(tracker.genome_path.keys()):
            print("\t%s" % organism, file=sys.stderr)
        print('', file=sys.stderr)
    else:
        retrieval_successful = tracker.retrieve_data(namespace.organism, namespace.overwrite)
        if retrieval_successful:
            tracker.derive_data(namespace.organism, namespace.overwrite)


def setup_subparsers(subparsers):
    retrieval_parser = subparsers.add_parser('install', help="Download data from Ensembl's main site")
    retrieval_parser.add_argument('organism', help="a species common name, e.g. human or mouse")
    retrieval_parser.add_argument('-o', '--overwrite', action="store_true", help="Overwrite existing files")
    retrieval_parser.set_defaults(func=retrieve_data)

    available_parser = subparsers.add_parser('available', help='Show organisms available for download')
    available_parser.set_defaults(func=ATGDataTracker().report_available)

    installed_parser = subparsers.add_parser('list', help='List installed organisms')
    installed_parser.set_defaults(func=ATGDataTracker().list_installed)
