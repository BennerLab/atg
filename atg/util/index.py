import os
import pysam
import string
import shlex
import subprocess
import argparse
import math


def get_genome_stats(genome_fasta):
    reference_fasta_index = genome_fasta + '.fai'
    if not os.path.exists(reference_fasta_index):
        print("\nIndexing %s\n" % os.path.abspath(genome_fasta))
        pysam.faidx(genome_fasta)

    reference_genome = pysam.FastaFile(genome_fasta)
    total_bases = sum(reference_genome.lengths)

    return reference_genome.nreferences, total_bases


class SequenceIndexer:
    """
    A base class for NGS read aligner indexing.
    """
    name = 'test'

    def __init__(self):
        self.command_template = string.Template('echo [executable] --input $input --output $output --threads $threads')
        self.index_dir = 'ExampleIndex'

    @classmethod
    def get_argument_parser(cls, indexer_argument_parser=None):
        """

        :param indexer_argument_parser:
        :return:
        """
        if not indexer_argument_parser:
            indexer_argument_parser = argparse.ArgumentParser()

        default_threads = max(1, subprocess.os.cpu_count()//2)

        indexer_argument_parser.set_defaults(Indexer=cls)
        indexer_argument_parser.add_argument('input', help="location of the sequence file")
        indexer_argument_parser.add_argument('-o', '--output', help="directory to store indexes")
        indexer_argument_parser.add_argument('-c', '--check', help="print commands to stdout instead of running them",
                                             action="store_true")
        indexer_argument_parser.add_argument('-t', '--threads', type=int,
                                             help="number of threads (default: %d)" % default_threads,
                                             default=default_threads)
        return indexer_argument_parser

    def index_sequences(self, **kwargs):
        # write indexes to a a default directory if nothing is specified
        if kwargs['output']:
            index_dir = kwargs['output']
        else:
            genome_path = os.path.abspath(kwargs['input'])
            genome_dir, genome_basename = os.path.split(genome_path)
            index_dir = os.path.join(genome_dir, self.index_dir)

        # replace specified input filename with absolute path
        command = self.command_template.safe_substitute(kwargs, input=os.path.abspath(kwargs['input']),
                                                        output=index_dir)

        if kwargs['check']:
            print(command)
            return True

        if not os.path.exists(index_dir):
            os.makedirs(index_dir)

        subprocess.run(args=shlex.split(command), env=os.environ.copy(), cwd=index_dir)
        return True


class STARIndexer(SequenceIndexer):
    # TODO: add support for gene annotation
    name = 'star'
    DEFAULT_SAINDEXNBASES = 14
    DEFAULT_CHRBINNBITS = 18

    def __init__(self):
        self.command_template = string.Template('STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $output '
                                                '--genomeFastaFiles $input --genomeSAindexNbases $nbases '
                                                '--genomeChrBinNbits $nbits')
        self.index_dir = 'STARIndex'

    def index_sequences(self, **kwargs):
        # assess sequence files to determine indexing parameters for STAR
        num_seq, total_bases = get_genome_stats(kwargs['input'])

        nbases = min(self.DEFAULT_SAINDEXNBASES, int(round(math.log(total_bases, 2))) - 1)
        nbits = min(self.DEFAULT_CHRBINNBITS, int(round(math.log(total_bases / num_seq, 2))))

        super().index_sequences(**kwargs, nbases=nbases, nbits=nbits)


class HISAT2Indexer(SequenceIndexer):
    # TODO: add support for sequence variants
    # TODO: add support for gene annotation
    name = 'hisat2'

    def __init__(self):
        self.command_template = string.Template('hisat2-build -p $threads $input $output/$index_name')
        self.index_dir = 'HISAT2Index'

    @classmethod
    def get_argument_parser(cls, indexer_argument_parser=None):
        indexer_argument_parser = super(HISAT2Indexer, cls).get_argument_parser(indexer_argument_parser)
        indexer_argument_parser.add_argument('-n', '--index_name', default='genome',
                                             help='base name for index files (default: genome)')
        return indexer_argument_parser


def run_index(namespace):
    indexer = namespace.Indexer()
    indexer.index_sequences(**vars(namespace))


def setup_subparsers(subparsers):
    indexing_parser = subparsers.add_parser('index', help='Index sequences for read alignment')
    indexing_subparser = indexing_parser.add_subparsers(title="indexer", dest="indexer",
                                                        description="available alignment programs")

    for indexer in [SequenceIndexer, STARIndexer, HISAT2Indexer]:
        current_subparser = indexing_subparser.add_parser(indexer.name)
        indexer.get_argument_parser(current_subparser)

    indexing_parser.set_defaults(func=run_index)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)
