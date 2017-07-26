import os
import argparse
import subprocess
import shlex
import pandas
import string

DEFAULT_THREADS = max(1, subprocess.os.cpu_count() // 2)


class ReadFileMismatch(Exception):
    pass


def group_files_from_sample_sheet(sample_sheet, force_pair=False):
    """

    :param sample_sheet: a delimited file containing filenames and sample IDs
    :param force_pair: interpret _1 _2 at end of file basename as R1 and R2
    :return: a dictionary containing sampleID -> string containing formatted file names
    """

    # read in sample sheet, splitting on whitespace or comma
    read_table = pandas.read_table(sample_sheet, header=None, sep="[ \t,]", names=['file', 'sampleID'], engine="python")
    read_files_strings = {}
    for sampleID, file_list in read_table.groupby("sampleID"):
        r1_samples = []
        r2_samples = []

        for filename in sorted(list(file_list.file)):
            if "R2" in filename or "mate2" in filename:
                r2_samples.append(filename)
            elif force_pair and os.path.basename(filename).split(os.path.extsep)[0].endswith('_2'):
                r2_samples.append(filename)
            else:
                r1_samples.append(filename)

        read_files_strings[sampleID] = (','.join(r1_samples) + ' ' + ','.join(r2_samples)).strip()

    return read_files_strings


def group_files_from_list(read_list, delimiter, force_pair=False):
    """

    :param read_list: list of fastq files
    :param delimiter:
    :param force_pair: interpret _1 _2 at end of file basename as R1 and R2
    :return: a dictionary containing sampleID -> string containing formatted file names
    """
    sample_group = {}
    for read_file in read_list:
        sample_name = os.path.basename(read_file).split(delimiter)[0]
        sample_group.setdefault(sample_name, []).append(read_file)

    read_files_strings = {}

    for sample_name, grouped_samples in sample_group.items():
        r1_samples = []
        r2_samples = []
        for sample in sorted(grouped_samples):
            if 'R2' in sample or 'mate2' in sample:
                r2_samples.append(sample)
            elif force_pair and os.path.basename(sample).split(os.path.extsep)[0].endswith('_2'):
                r2_samples.append(sample)
            else:
                r1_samples.append(sample)

        # check pairing of reads
        if len(r2_samples) > 0 and len(r2_samples) != len(r1_samples):
            raise ReadFileMismatch('Different number of paired-end read files found.')

        read_files_strings[sample_name] = (' '.join([','.join(r1_samples), ','.join(r2_samples)])).strip()

    return read_files_strings


class ReadAlignmentHelper:
    """
    A base class for NGS read aligners. Common functions are:
        - aligning reads from a list of file names
        - aligning reads from a given file (sample sheet)
        - summarizing results
    """

    name = 'test'

    def __init__(self):
        """

        """
        self.command_template = string.Template('echo [executable] --index $index --output $output --threads $threads '
                                                '--input $input')

    @classmethod
    def get_argument_parser(cls, aligner_argument_parser=None):
        """
        Common parameters:
            - an index (usually genome or transcriptome)
            - output directory location
            - threads/processes to run in parallel
            -
        """
        if not aligner_argument_parser:
            aligner_argument_parser = argparse.ArgumentParser()

        aligner_argument_parser.set_defaults(Aligner=cls)
        aligner_argument_parser.add_argument('index', help="location of the reference index")
        aligner_argument_parser.add_argument('output', help="directory to store all output")
        aligner_argument_parser.add_argument('fastq', nargs='+', help="input files or sample sheet")
        aligner_argument_parser.add_argument('-s', '--sample_sheet', action="store_true",
                                             help="use a sample sheet for read files and sample IDs",)
        aligner_argument_parser.add_argument('-c', '--check', help="output commands to stdout instead of running them",
                                             action="store_true")
        aligner_argument_parser.add_argument('-t', '--threads', help="number of threads", type=int,
                                             default=DEFAULT_THREADS)
        aligner_argument_parser.add_argument("-d", "--delimiter", help="delimiter for sample ID extraction (default _)",
                                             default="_")
        aligner_argument_parser.add_argument('-p', '--force_pair', action="store_true",
                                             help='interpret "_2".fastq.gz as a paired-end read')

        return aligner_argument_parser

    def format_output_name(self, output_path, sample_name):
        return os.path.join(output_path, sample_name)

    def set_params(self):
        pass

    @staticmethod
    def group_reads(**kwargs):
        # group reads based on input as sample sheet or filename list
        if kwargs['sample_sheet']:
            read_files_dict = group_files_from_sample_sheet(kwargs['fastq'][0], kwargs['force_pair'])
        else:
            read_files_dict = group_files_from_list(kwargs['fastq'], kwargs['delimiter'], kwargs['force_pair'])

        return read_files_dict

    def get_command_list(self, read_files_dict, **kwargs):
        command_list = []
        for sample_name in sorted(read_files_dict.keys()):
            output_file = os.path.join(kwargs['output'], sample_name)
            command_list.append(self.command_template.safe_substitute(kwargs, read_files=read_files_dict[sample_name],
                                                                      basename=output_file))
        return command_list

    def align_reads(self, **kwargs):
        """
        primary method, which takes care of:
        1. grouping read files into a dictionary of samples
        2. generating a list of commands to be run
        3. executing the commands via subprocess
        
        :return:
        """

        read_files_dict = self.group_reads(**kwargs)
        command_list = self.get_command_list(read_files_dict, **kwargs)

        # check commands or run them
        if kwargs['check']:
            for command in command_list:
                print(command)

        else:
            for command in command_list:
                subprocess.run(args=shlex.split(command), env=os.environ.copy())

        return True

    def summarize_results(self):
        pass


class STARAligner(ReadAlignmentHelper):
    name = 'star'

    def __init__(self):
        """
        run spliced alignment with STAR
            - assume BAM output
        """

        self.command_template = string.Template(
            "STAR --genomeDir $index --genomeLoad NoSharedMemory --runThreadN $threads "
            "--readFilesIn $read_files --outFileNamePrefix $basename. "
            "--outSAMtype BAM SortedByCoordinate")

        # NoSharedMemory
        # --genomeLoad LoadAndKeep
        # --limitBAMsortRAM 50000000000

    @classmethod
    def get_argument_parser(cls, aligner_argument_parser=None):
        aligner_argument_parser = super(STARAligner, cls).get_argument_parser(aligner_argument_parser)
        aligner_argument_parser.add_argument('-k', '--keep_unmapped', action='store_true')
        return aligner_argument_parser

    def get_command_list(self, read_files_dict, **kwargs):
        command_list = []
        for sample_name in sorted(read_files_dict.keys()):
            output_file = os.path.join(kwargs['output'], sample_name)
            # compressed files need some additional handling
            # assume that all grouped files are compressed with the same format
            input_string = read_files_dict[sample_name]
            file_extension = os.path.splitext(input_string)[1]
            if file_extension == '.gz':
                additional_options = ' --readFilesCommand gunzip -c'
            elif file_extension == '.bz2':
                additional_options = ' --readFilesCommand bunzip2 -c'
            else:
                additional_options = ''

            command_list.append(self.command_template.safe_substitute(kwargs, read_files=read_files_dict[sample_name],
                                                                      basename=output_file) + additional_options)
        return command_list


class HISAT2Aligner(ReadAlignmentHelper):
    name = 'hisat2'

    def __init__(self):
        super().__init__()
        self.command_template = string.Template("hisat2 --new-summary -x $index -p $threads $processed_input")

    @classmethod
    def get_argument_parser(cls, aligner_argument_parser=None):
        aligner_argument_parser = super(HISAT2Aligner, cls).get_argument_parser(aligner_argument_parser)
        aligner_argument_parser.add_argument('--bam', action='store_true', help="write sorted BAM rather than "
                                                                                "unsorted SAM output")
        aligner_argument_parser.add_argument('--dta', action='store_true', help="add information for downstream"
                                                                                "transcriptome assembly")
        return aligner_argument_parser

    def get_command_list(self, read_files_dict, **kwargs):
        command_list = []
        for sample_name in sorted(read_files_dict.keys()):
            output_file = os.path.join(kwargs['output'], sample_name)

            # manipulate read file input for HISAT2 command
            read_list = read_files_dict[sample_name].split(' ')
            if len(read_list) == 1:
                processed_input = '-U %s' % read_list[0]
            elif len(read_list) == 2:
                processed_input = '-1 %s -2 %s' % (read_list[0], read_list[1])
            else:
                raise ReadFileMismatch('Improperly formatted read file input:\n%s\n' % read_files_dict[sample_name])

            additional_options = ''
            if kwargs['dta']:
                additional_options += ' --dta'

            additional_options += ' -S %s.sam' % output_file

            command_list.append(self.command_template.safe_substitute(kwargs, processed_input=processed_input,
                                                                      output_basename=output_file) +
                                additional_options)
        return command_list

    def align_reads(self, **kwargs):
        """
        primary method, which takes care of:
        1. grouping read files into a dictionary of samples
        2. generating a list of commands to be run
        3. executing the commands via subprocess

        :return:
        """

        read_files_dict = self.group_reads(**kwargs)
        command_list = self.get_command_list(read_files_dict, **kwargs)

        if kwargs['bam']:
            samtools_bam_command = shlex.split('samtools view -u')

            for alignment_command in command_list:
                alignment_command_list = shlex.split(alignment_command)
                # remove the last two elements from the list (e.g. '-S' and 'output.sam')
                # the first element removed should be the SAM output filename
                output_basename = alignment_command_list.pop()[:-4]
                output_filename = output_basename + '.bam'
                log_filename = output_basename + '.log'
                alignment_command_list.pop()

                samtools_sort_command = shlex.split('samtools sort -@ %d -o %s' % (DEFAULT_THREADS, output_filename))

                if kwargs['check']:
                    command_pipe = ' | '.join([
                        ' '.join(alignment_command_list),
                        ' '.join(samtools_bam_command),
                        ' '.join(samtools_sort_command)
                    ])
                    print(command_pipe)
                else:
                    with open(log_filename, 'w') as log:
                        aligner = subprocess.Popen(alignment_command_list, stdout=subprocess.PIPE, stderr=log)
                        samtools_bam = subprocess.Popen(samtools_bam_command, stdin=aligner.stdout,
                                                        stdout=subprocess.PIPE)
                        samtools_sort = subprocess.Popen(samtools_sort_command, stdin=samtools_bam.stdout)
                        samtools_sort.communicate()

        else:
            for alignment_command in command_list:
                if kwargs['check']:
                    print(alignment_command)
                else:
                    output_basename = alignment_command.split()[-1][:-4]
                    log_filename = output_basename + '.log'
                    with open(log_filename, 'w') as log:
                        subprocess.run(args=shlex.split(alignment_command), env=os.environ.copy(), stderr=log)

        return True


def run_aligner(namespace):
    aligner = namespace.Aligner()
    aligner.align_reads(**vars(namespace))


def setup_subparsers(subparsers):
    aligner_parser = subparsers.add_parser('align', help='Align reads to reference genome/transcriptome')
    aligner_subparser = aligner_parser.add_subparsers(title="aligner", dest="aligner",
                                                        description="available alignment programs")

    for aligner in [ReadAlignmentHelper, STARAligner, HISAT2Aligner]:
        current_subparser = aligner_subparser.add_parser(aligner.name)
        aligner.get_argument_parser(current_subparser)

    aligner_parser.set_defaults(func=run_aligner)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)
