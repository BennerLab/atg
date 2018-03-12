"""

"""

import sys
import atg.config
import atg.data


def ucsc_chromosomes(interval_df):
    return (interval_df.loc[interval_df['chrom'].str.len() <= 2, :]
                       .replace({'chrom': {'MT': 'M'}})
                       .assign(chrom='chr' + interval_df['chrom']))


def generate_gene_bed(organism, gene_type=None, ucsc_names=False):
    bed_df = atg.data.get_gene_bed6(organism)

    if gene_type:
        bed_df = bed_df.loc[bed_df['Gene type'] == gene_type]

    if ucsc_names:
        bed_df = ucsc_chromosomes(bed_df)

    bed_df.iloc[:, 0:6].sort_values(['chrom', 'start', 'end']).to_csv(sys.stdout, header=None, sep='\t', index=False)


def generate_tss_bed(organism, transcript_type=None, ucsc_names=False):
    transcript_df = atg.data.get_tss_bed6(organism, level='symbol')

    if transcript_type:
        transcript_df = transcript_df.loc[transcript_df['Transcript type'] == transcript_type]

    if ucsc_names:
        transcript_df = ucsc_chromosomes(transcript_df)

    transcript_df.iloc[:, 0:6].sort_values(['chrom', 'start', 'end']).drop_duplicates().to_csv(sys.stdout, header=None,
                                                                                               sep='\t', index=False)


def generate_intervals(namespace):
    if namespace.organism not in atg.data.genome_path:
        print('%s is not a valid organism choice.\nTry one of the following: %s' %
              (namespace.organism, ' '.join(sorted(atg.data.genome_path.keys()))))

    if namespace.type == 'gene':
        generate_gene_bed(namespace.organism, namespace.protein, namespace.ucsc)
    elif namespace.type == 'tss':
        generate_tss_bed(namespace.organism, namespace.protein, namespace.ucsc)
    else:
        raise ValueError('Invalid interval type selected.')


def setup_subparsers(subparsers):
    interval_parser = subparsers.add_parser('interval', help="Generate BED files from genome annotation")
    interval_parser.add_argument('organism', help="a UCSC genome code, e.g. hg38\n"
                                                  "in the future, this will likely change to a common name")
    interval_parser.add_argument('type', choices=['gene', 'tss'], help='BED file for gene body or TSS')
    interval_parser.add_argument('--ucsc', help="Use UCSC chromosome names", action="store_true")
    interval_parser.add_argument('-p', '--protein', help="Generate intervals only for protein-coding genes",
                                 action='store_const', const='protein_coding')
    interval_parser.set_defaults(func=generate_intervals)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)
