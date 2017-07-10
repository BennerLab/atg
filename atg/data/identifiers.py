import pandas
import os
import atg.data

ID_COLUMN_NAME = ['ensembl', 'symbol']
ID_COLUMN_INDEX = [0, 6]
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
        self.gene_annotation = pandas.read_csv(gene_annotation_location, names=ID_COLUMN_NAME, usecols=ID_COLUMN_INDEX,
                                               skiprows=1)

    def translate_identifiers(self, input_series, input_type, output_type):
        if not input_type:
            input_type = guess_identifier_type(input_series)

        if input_type not in VALID_ID_TYPES:
            raise ValueError('Input type %s is not supported.' % input_type)
        if output_type not in VALID_ID_TYPES:
            raise ValueError('Output type %s is not supported.' % output_type)

        input_dataframe = input_series.to_frame(input_type)
        merged_dataframe = input_dataframe.merge(self.gene_annotation, on=input_type)
        return merged_dataframe[output_type].drop_duplicates()
