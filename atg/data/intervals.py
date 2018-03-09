"""

"""

import os
import sys
import pandas
import atg.config
import atg.data


def generate_gene_bed(organism, gene_type=None, ucsc_names=False):
    gene_annotation_path = os.path.join(atg.data.genome_path[organism], 'ensembl_gene.csv')
    gene_df = pandas.read_csv(gene_annotation_path)
    bed_df = gene_df.loc[:, ['Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)', 'Transcript count',
                             'Gene name', 'Strand', 'Gene type']]
    bed_df.replace({'Strand': {-1: '-', 1: '+'}}, inplace=True)
    bed_df.rename(index=str, columns={'Chromosome/scaffold name': 'chrom', 'Gene start (bp)': 'start',
                                      'Gene end (bp)': 'end'}, inplace=True)

    if gene_type:
        bed_df = bed_df.loc[bed_df['Gene type'] == gene_type]

    if ucsc_names:
        bed_df = bed_df.loc[bed_df['chrom'].str.len() <= 2, :]
        bed_df.replace({'chrom': {'MT': 'M'}}, inplace=True)
        bed_df['chrom'] = 'chr' + bed_df['chrom']

    bed_df.iloc[:, 0:6].sort_values(['chrom', 'start', 'end']).to_csv(sys.stdout, header=None, sep='\t', index=False)


def generate_tss_bed(organism, gene_type=None, ucsc_names=False):
    pass


def generate_intervals(namespace):
    if namespace.organism not in atg.data.genome_path:
        print('%s is not a valid organism choice.\nTry one of the following: %s' %
              (namespace.organism, ' '.join(sorted(atg.data.genome_path.keys()))))

    if namespace.type == 'gene':
        generate_gene_bed(namespace.organism, 'protein_coding', True)
        # generate_gene_bed(namespace.organism, None, True)
        # generate_gene_bed(namespace.organism, None, False)

    elif namespace.type == 'tss':
        generate_tss_bed(namespace.organism, None, namespace.ucsc)
    else:
        raise ValueError('Invalid interval type selected.')


def setup_subparsers(subparsers):
    interval_parser = subparsers.add_parser('interval', help="Generate BED files from genome annotation")
    interval_parser.add_argument('organism', help="a UCSC genome code, e.g. hg38\n"
                                                  "in the future, this will likely change to a common name")
    interval_parser.add_argument('type', choices=['gene', 'tss'], help='BED file for gene body or TSS')
    interval_parser.add_argument('--ucsc', help="Use UCSC chromosome names", action="store_true")
    interval_parser.add_argument('-p', '--protein', help="Generate intervals only for protein-coding genes")
    interval_parser.set_defaults(func=generate_intervals)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)
