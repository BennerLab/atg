import os
import configparser
import numpy
import pandas
import atg.config
import atg.data.retrieve
import atg.data.ensembl

DATA_SOURCE_PATH = os.path.join(os.path.dirname(__file__), 'data_sources.ini')

data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
config = configparser.ConfigParser(default_section="Common",
                                   interpolation=configparser.ExtendedInterpolation())
config.read(DATA_SOURCE_PATH)
genome_path = {}

for organism in config.sections():
    ensembl_genome = config.get(organism, 'ensembl_genome')
    current_genome_path = os.path.join(data_root, organism, ensembl_genome)
    genome_path[organism] = current_genome_path


def get_gene_bed6(organism_common_name):
    """
    read Ensembl gene annotation file, reformatting strand column and returning only BED6-like columns.

    :param organism_common_name:
    :return: a BED6-like DataFrame with gene intervals, including a seventh column for gene type
    """
    gene_annotation_path = os.path.join(genome_path[organism_common_name], 'ensembl_gene.csv')
    gene_df = pandas.read_csv(gene_annotation_path)
    bed_df = gene_df.loc[:, ['Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)', 'Transcript count',
                             'Gene name', 'Strand', 'Gene type', 'Gene stable ID']]
    bed_df.replace({'Strand': {-1: '-', 1: '+'}}, inplace=True)
    bed_df.rename(index=str, columns={'Chromosome/scaffold name': 'chrom',
                                      'Gene start (bp)': 'start',
                                      'Gene end (bp)': 'end'}, inplace=True)
    return bed_df


def get_tss_bed6(organism_common_name, level="transcript"):
    """
    read Ensembl transcript annotation file, reformatting strand column and returning only BED6-like columns. Adjusts
    given start positions based on strand. Includes gene-level information when specified.

    :param organism_common_name:
    :param level: for BED field, use either ensembl gene ID ('gene'), transcript id ('transcript'), or
    gene symbol ('symbol'). Duplicate removal will take this value into account.
    :return: a BED6-like DataFrame with TSS intervals, including a seventh column for transcript type
    """

    transcript_annotation_path = os.path.join(genome_path[organism_common_name], 'ensembl_transcript.csv')
    transcript_df = pandas.read_csv(transcript_annotation_path)
    transcript_df.replace({'Strand': {-1: '-', 1: '+'}}, inplace=True)
    transcript_df.rename(index=str, columns={'Chromosome/scaffold name': 'chrom',
                                             'Strand': 'strand'}, inplace=True)

    # assign TSS start/end based on strand
    transcript_bed = transcript_df.assign(start=numpy.where(transcript_df['strand'] == '+',
                                                            transcript_df['Transcription start site (TSS)'],
                                                            transcript_df['Transcription start site (TSS)'] - 1),
                                          end=numpy.where(transcript_df['strand'] == '+',
                                                          transcript_df['Transcription start site (TSS)'] + 1,
                                                          transcript_df['Transcription start site (TSS)']),
                                          score=0,
                                          name=transcript_df['Transcript stable ID'])

    # assign gene/symbol if specified
    if level == 'symbol':
        gene_df = get_gene_bed6(organism_common_name)
        transcript_bed = transcript_bed.merge(gene_df.loc[:, ['Gene stable ID', 'Gene name']])
        transcript_bed = transcript_bed.assign(name=transcript_bed['Gene name'])
    elif level == 'gene':
        transcript_bed = transcript_bed.assign(name=transcript_bed['Gene stable ID'])

    # remove redundancy
    transcript_bed = (transcript_bed.loc[:, ['chrom', 'start', 'end', 'score', 'name', 'strand', 'Transcript type']]
                                    .drop_duplicates())

    return transcript_bed


def setup_subparsers(subparsers):
    data_parser = subparsers.add_parser('data', help='Retrieve genomes and annotation')
    data_subparsers = data_parser.add_subparsers(title='Retrieve genomes and annotation')
    data_parser.set_defaults(func=lambda x: data_parser.print_help())

    atg.data.retrieve.setup_subparsers(data_subparsers)
    atg.data.ensembl.setup_subparsers(data_subparsers)
