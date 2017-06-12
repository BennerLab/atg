"""
Find species data in Ensembl, recording genome and annotation URLs.
"""

import os
import sys
import pandas
import ftplib
import string
import atg.config
import atg.data.retrieve

ENSEMBL_SPECIES_INFORMATION = 'ftp://ftp.ensemblgenomes.org/pub/current/species.txt'
ENSEMBL_DNA_BASE_LOCATION = string.Template('pub/current/$division/fasta$collection/$species/dna/')
ENSEMBL_GTF_BASE_LOCATION = string.Template('pub/current/$division/gtf$collection/$species/$assembly.'
                                            '$version.gtf.gz')


class EnsemblSpecies:
    """
    A class for fetching and managing species data from Ensembl Genomes, which include many organisms not found on
    the main Ensembl site. Files for these organisms are stored in individual subfolders in e.g.
    ~/ATGData/ensemblgenomes/.

    """

    def __init__(self):
        self.data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
        ensembl_genome_file = os.path.join(self.data_root, 'ensembl_species.txt')
        if not os.path.exists(ensembl_genome_file):
            atg.data.retrieve.fetch_url(ENSEMBL_SPECIES_INFORMATION, ensembl_genome_file)
        self.ensembl_species_df = pandas.read_table(ensembl_genome_file, index_col=False)

    def get_species_information(self, species):
        """

        :param species: genus and species (as named by Ensembl), e.g. zea_mays
        :return: dictionary containing URLs to genome fasta and gene annotation (GTF), if found
        """

        if sum(self.ensembl_species_df.species.isin([species])) == 0:
            return {'species': species}

        # pull out first matching record
        ensembl_record = self.ensembl_species_df.ix[self.ensembl_species_df['species'] == species].iloc[0]
        ensembl_division = ensembl_record.ix['division'].lstrip('Ensembl').lower()
        # could access assembly ID or accession from record, but the Ensembl files don't use one consistently

        ensembl_core_db = ensembl_record.ix['core_db']
        if "collection" in ensembl_core_db:
            collection_path = '/' + ensembl_core_db.split('_core_')[0]
        else:
            collection_path = ''

        with ftplib.FTP('ftp.ensemblgenomes.org') as ftp:
            ftp.login()
            genome_listing = ftp.nlst(ENSEMBL_DNA_BASE_LOCATION.safe_substitute(division=ensembl_division,
                                                                                species=species,
                                                                                collection=collection_path))
            genome_location = ''
            annotation_location = ''
            genome_assembly_version = ''

            # find toplevel unmasked genome
            for filename in genome_listing:
                if 'dna.toplevel' in filename:
                    genome_location = filename
                    break

            if genome_location != '':
                genome_filename = genome_location.split('/')[-1]
                genome_assembly = genome_filename.rstrip('.dna.toplevel.fa.gz')
                genome_assembly_version = genome_assembly.split('.', maxsplit=1)[1]

                annotation_listing = ftp.nlst(ENSEMBL_GTF_BASE_LOCATION.safe_substitute(division=ensembl_division,
                                                                                        species=species,
                                                                                        assembly=genome_assembly,
                                                                                        collection=collection_path,
                                                                                        version=35))

                if len(annotation_listing) == 0:
                    annotation_location = ''
                elif len(annotation_listing) == 1:
                    annotation_location = annotation_listing[0]
                else:
                    annotation_location = 'multiple'

            ftp.close()

        return {'species': species, 'genome': genome_location, 'annotation': annotation_location,
                'version': genome_assembly_version}

    def collect_species_information(self, species_list):
        """
        Given a list of species names, create a dataframe containing all information
        :param species_list:
        :return: dataframe
        """

        record_list = []

        for species in species_list:
            record_list.append(self.get_species_information(species))

        return pandas.DataFrame.from_records(record_list)

    def retrieve_species_data(self, species):
        """
        Download data from Ensembl.
        :param species:
        :return: True if successful
        """

        species_information = self.get_species_information(species)
        if len(species_information) == 1:
            return False

        ensembl_species_path = os.path.join(self.data_root, 'ensemblgenomes', species)
        os.makedirs(ensembl_species_path, exist_ok=True)
        for filetype in ('genome', 'annotation'):
            filename = os.path.split(species_information[filetype])[-1]
            ensembl_url = 'ftp://ftp.ensemblgenomes.org/' + species_information[filetype]
            output_filename = os.path.join(ensembl_species_path, filename)
            atg.data.retrieve.fetch_url(ensembl_url, output_filename)

        return True


def retrieve_ensembl_species(namespace):
    # get list of species from file or namespace
    if namespace.list:
        species_list = pandas.read_csv(namespace.species_name[0], index_col=False, header=None).iloc[:, 0].tolist()
    else:
        species_list = namespace.species_name

    tracker = EnsemblSpecies()
    # output species information as table, or download
    if namespace.table:
        species_df = tracker.collect_species_information(species_list)
        species_df.to_csv(sys.stdout, sep="\t", index=False, columns=['species', 'genome', 'annotation', 'version'])
    else:
        for species in species_list:
            retrieval_success = tracker.retrieve_species_data(species)
            if retrieval_success:
                print('%s retrieved successfully.' % species)
            else:
                print('%s information not retrieved.' % species)


def setup_subparsers(subparsers):
    retrieval_parser = subparsers.add_parser('species', help="Retrieve genome sequences and related information")
    retrieval_parser.add_argument('species_name', nargs="+",
                                  help="one or more genus/species for an organism in Ensembl, e.g. zea_mays")
    retrieval_parser.add_argument('-l', '--list', action="store_true", help="species are provided in a text file given"
                                                                            "as the only argument")
    retrieval_parser.add_argument('-t', '--table', action="store_true",
                                  help="instead of downloading data, write the species information to stdout")
    # retrieval_parser.add_argument('-o', '--overwrite', action="store_true", help="Overwrite existing files")
    retrieval_parser.set_defaults(func=retrieve_ensembl_species)
