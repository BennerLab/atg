#!/usr/bin/env python3

import os
import subprocess
import pandas


def merge_fastq(filename_list, execute=False):
    """
    concatenate fastq.gz files that have been split across multiple sequencing lanes.
    """
    fastq_df = pandas.DataFrame(filename_list, columns=['filename'])
    fastq_merge_df = fastq_df.filename.apply(os.path.basename).str.extract(r'(.*)_S\d+_L\d\d\d_(R\d)_\d\d\d', expand=True).join(fastq_df)
    fastq_merge_df['output'] = fastq_merge_df[0] + '_' + fastq_merge_df[1] + '.fastq.gz'

    for output_name, group in fastq_merge_df.groupby('output'):
        if execute:
            command_args = ['cat'] + group['filename'].sort_values().tolist()
            with open(output_name, 'w') as output:
                subprocess.run(command_args, stdout=output, check=True)
        else:
            pandas.set_option("display.max_colwidth", 120)
            print(output_name)
            print("\t" + group.to_string(columns=['filename'], header=False, index=False, formatters={'filename': lambda x: "\t" + x}) + "\n")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename_list', nargs='+')
    parser.add_argument('-x', '--execute', action='store_true', help="execute merge (just prints commands by default)")

    args = parser.parse_args()
    merge_fastq(args.filename_list, args.execute)

