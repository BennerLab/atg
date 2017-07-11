import pandas
import os
import atg.data

ID_COLUMN_NAME = ['ensembl', 'chromosome', 'symbol']
ID_COLUMN_INDEX = [0, 2, 6]
VALID_ID_TYPES = set(ID_COLUMN_NAME)  # entrez not supported yet


def guess_identifier_type(input_series):
    # ENSEMBL identifiers all begin with ENS
    if all([str(x).startswith('ENS') for x in input_series]):
        return 'ensembl'
    # assume integer values are Entrez IDs
    elif all([str(x).isdigit() for x in input_series]):
        return 'entrez'
    # otherwise, gene symbols
    else:
        return 'symbol'


class GeneIDTranslator:
    def __init__(self, species):
        self.species = species

        if species not in atg.data.genome_path.keys():
            supported_species = '\n'.join(sorted(atg.data.genome_path.keys()))
            raise ValueError('\n\n%s is not a supported species.\n\nTry one of the following:\n%s'
                             % (species, supported_species))

        gene_annotation_location = os.path.join(atg.data.genome_path[species], 'ensembl_gene.csv')
        gene_annotation = pandas.read_csv(gene_annotation_location, names=ID_COLUMN_NAME, usecols=ID_COLUMN_INDEX,
                                          skiprows=1)

        # sort by chromosome string length and chromosome to prioritize main chromosomes
        self.gene_annotation = (gene_annotation.assign(chr_strlen = lambda x: x.chromosome.str.len())
                                               .sort_values(['chr_strlen', 'chromosome'])
                                               .drop(['chr_strlen', 'chromosome'], axis=1))

    def translate_identifiers(self, input_series, input_type, output_type):
        """
        get a Series of translated IDs
        - input IDs that don't map to anything are ignored
        - input IDs are limited to at most one output ID
        :param input_series:
        :param input_type:
        :param output_type:
        :return:
        """

        mapping_df = self.map_identifiers(input_series, input_type, output_type)
        deduplicated_series = mapping_df.drop_duplicates(mapping_df.columns[0])[output_type].drop_duplicates()
        return deduplicated_series

    def map_identifiers(self, input_series, input_type, output_type):
        """
        make a DataFrame containing all mappings between input and translated IDs
        :param input_series:
        :param input_type:
        :param output_type:
        :return:
        """
        if not input_type:
            input_type = guess_identifier_type(input_series)

        if input_type == output_type:
            return input_series

        if input_type not in VALID_ID_TYPES:
            raise ValueError('Input type %s is not supported.' % input_type)
        if output_type not in VALID_ID_TYPES:
            raise ValueError('Output type %s is not supported.' % output_type)

        input_dataframe = input_series.to_frame(input_type)
        merged_dataframe = input_dataframe.merge(self.gene_annotation, how="left", on=input_type)

        return merged_dataframe
