#!/usr/bin/env python

import pandas

BEDGRAPH_COLUMNS = ['chr', 'start', 'end', 'score']
BED6_COLUMNS = ['chr', 'start', 'end', 'name', 'score', 'strand']
ENSEMBL_CHR = [str(x) for x in range(1, 24)] + ['X', 'Y']
UCSC_CHR = [('chr' + x) for x in ENSEMBL_CHR]
CHR_DICT = dict(zip(ENSEMBL_CHR, UCSC_CHR))
CHR_DICT['MT'] = 'chrM'


def merge_bedgraph(positive_bedgraph, negative_bedgraph, convert_ucsc):
    pos_df = (pandas.read_csv(positive_bedgraph, sep='\t', names=BEDGRAPH_COLUMNS, dtype={'chr': 'str'})
                    .assign(strand='+', name='.'))
    neg_df = (pandas.read_csv(negative_bedgraph, sep='\t', names=BEDGRAPH_COLUMNS, dtype={'chr': 'str'})
                    .assign(strand='-', name='.'))

    merged_df: pandas.DataFrame = pandas.concat([pos_df, neg_df])
    if convert_ucsc:
        basic_chromosomes = merged_df.chr.isin(CHR_DICT.keys())
        merged_df = merged_df.loc[basic_chromosomes].replace(CHR_DICT)

    merged_df = merged_df.sort_values(['chr', 'start', 'end']).loc[:, BED6_COLUMNS]
    return merged_df


if __name__ == '__main__':
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('bedgraph', nargs='+', 
                        help='')
    parser.add_argument('-u', '--ucsc', action='store_true',
                        help='Naively change chromosome names to UCSC-style names')

    args = parser.parse_args()
    merge_bedgraph(args.bedgraph[0], args.bedgraph[1], args.ucsc).to_csv(sys.stdout, sep='\t', index=False,
                                                                         header=False)
