import os
import argparse
import subprocess
import shlex
import pandas
import string
import progress.bar

# ENHANCEMENT: use functionality from MultiQC for parsing. MultiQC is not currently written to allow convenient access.

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


class ReadAlignmentBase:
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
                                                '--input $read_files $extra')
        self.alignment_list = []
        self.log_list = []

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
        aligner_argument_parser.add_argument('extra', nargs=argparse.REMAINDER,
                                             help='additional options to be passed through')

        return aligner_argument_parser

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
            extra_options = ' '.join(kwargs['extra'])
            command_list.append(self.command_template.safe_substitute(kwargs, read_files=read_files_dict[sample_name],
                                                                      basename=output_file, extra=extra_options))
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
            # create output directory if it doesn't already exist; ignore for this base class
            if self.name != 'test':
                os.makedirs(kwargs['output'], exist_ok=TRUE)

            progress_bar = progress.bar.Bar('Processed %(index)d/%(max)d',
                                            suffix='Remaining: %(eta_td)s', max=len(command_list))
            progress_bar.start()
            for command in command_list:
                subprocess.run(args=shlex.split(command), env=os.environ.copy(), stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
                progress_bar.next()
            progress_bar.finish()
            print('\nCompleted in %s.\n' % progress_bar.elapsed_td)

        return True

    def process_alignments(self):
        pass

    def summarize_results(self):
        pass

    @classmethod
    def parse_log(cls, filename):
        pass


class STARAligner(ReadAlignmentBase):
    name = 'star'
    MAX_SORTING_BINS = 50
    MAX_OPEN_FILES = 500

    def __init__(self):
        """
        run spliced alignment with STAR
            - assume BAM output
        """
        super().__init__()
        self.command_template = string.Template(
            "STAR --genomeDir $index --runThreadN $threads "
            "--readFilesIn $read_files --outFileNamePrefix $basename. "
            "--outSAMtype BAM SortedByCoordinate $extra")

    @classmethod
    def get_argument_parser(cls, aligner_argument_parser=None):
        aligner_argument_parser = super(STARAligner, cls).get_argument_parser(aligner_argument_parser)
        aligner_argument_parser.add_argument('-q', '--quantify', action='store_true', help="count reads per gene")
        aligner_argument_parser.add_argument('-k', '--keep_unmapped', action='store_true',
                                             help="write unmapped reads to file")
        aligner_argument_parser.add_argument('--low_resource', action='store_true',
                                             help="avoid using shared memory")
        return aligner_argument_parser

    def get_command_list(self, read_files_dict, **kwargs):
        extra_options = ' '.join(kwargs['extra'])

        command_list = []
        for sample_name in sorted(read_files_dict.keys()):
            output_file = os.path.join(kwargs['output'], sample_name)
            # compressed files need some additional handling
            # assume that all grouped files are compressed with the same format
            input_string = read_files_dict[sample_name]
            file_extension = os.path.splitext(input_string)[1]

            additional_options = ''
            if file_extension == '.gz':
                additional_options += ' --readFilesCommand gunzip -c'
            elif file_extension == '.bz2':
                additional_options += ' --readFilesCommand bunzip2 -c'

            if not kwargs['low_resource']:
                sorting_bins = min(STARAligner.MAX_SORTING_BINS, STARAligner.MAX_OPEN_FILES//kwargs['threads'])
                # STAR will create temp files ~ threads*bins, and the system will limit this to approximately 1,000.
                # default is 50, which creates a problem when using 48 threads

                additional_options += ' --genomeLoad LoadAndKeep'
                additional_options += ' --limitBAMsortRAM 50000000000'
                additional_options += f' --outBAMsortingBinsN {sorting_bins}'
            if kwargs['keep_unmapped']:
                additional_options += ' --outReadsUnmapped Fastx'
            if kwargs['quantify']:
                additional_options += ' --quantMode GeneCounts'

            command_list.append(self.command_template.safe_substitute(kwargs, read_files=read_files_dict[sample_name],
                                                                      basename=output_file, extra=extra_options) +
                                additional_options)

            if not kwargs['check']:
                self.alignment_list.append(output_file + '.Aligned.sortedByCoord.out.bam')
                self.log_list.append(output_file + '.Log.final.out')
        return command_list

    def process_alignments(self):
        for alignment_file in self.alignment_list:
            new_alignment_filename = alignment_file.replace('.Aligned.sortedByCoord.out.bam', '.bam')
            os.rename(alignment_file, new_alignment_filename)

    def summarize_results(self):
        if len(self.log_list) == 0:
            return

        summary_list = []
        for log_file in self.log_list:
            summary_list.append(self.parse_log(log_file))

        # output a full result file containing all values from original Log.final.out
        result_df = pandas.DataFrame(summary_list)
        result_df.to_csv('star_alignment_complete.txt', sep='\t', index=False)

        # output simplified results too
        # total reads, uniquely mapped, uniquely mapped %, multi-mapped %, unmapped %
        simplified_df = result_df.loc[:, ['Sample', 'Number of input reads', 'Uniquely mapped reads number',
                                          'Uniquely mapped reads %', '% of reads mapped to multiple loci']]
        simplified_df.rename(index=str, columns={'Number of input reads': 'Total reads',
                                                 'Uniquely mapped reads number': 'Uniquely mapped',
                                                 'Uniquely mapped reads %': 'Unique %',
                                                 '% of reads mapped to multiple loci': 'Multimapped %'},
                             inplace=True)
        simplified_df['Unmapped %'] = 1 - simplified_df['Unique %'] - simplified_df['Multimapped %']
        simplified_df.to_csv('star_alignment_summary.txt', sep='\t', index=False)

        print(simplified_df.to_string(index=False, formatters={'Total reads': '{:,}'.format,
                                                               'Uniquely mapped': '{:,}'.format,
                                                               'Unique %': '{:.1%}'.format,
                                                               'Multimapped %': '{:.1%}'.format,
                                                               'Unmapped %': '{:.1%}'.format}))

    @classmethod
    def parse_log(cls, filename):
        """
        parse a STAR log, returning a pandas.Series
        :return: a Series
        """

        log_entries = ['Sample']
        log_values = [os.path.split(filename)[-1].replace('.Log.final.out', '')]

        with open(filename) as log_input:
            for line in log_input:
                if 'Started' in line or 'Finished' in line:
                    continue
                elif '|' in line:
                    name = line.split('|')[0].strip()
                    value = line.split('|')[1].strip()
                    if '%' in value or '.' in value:
                        value = float(value[:-1]) / 100.0  # change % to float
                    else:
                        value = int(value)
                    log_entries.append(name)
                    log_values.append(value)
        return pandas.Series(log_values, log_entries)


class SalzbergAligner(ReadAlignmentBase):
    """
    Base class for aligners from the Salzberg lab, which have similar options and do not natively output sorted BAM
    files. (Bowtie2 and HISAT2)
    """

    @classmethod
    def get_argument_parser(cls, aligner_argument_parser=None):
        aligner_argument_parser = super(SalzbergAligner, cls).get_argument_parser(aligner_argument_parser)
        aligner_argument_parser.add_argument('--bam', action='store_true', help="write sorted BAM rather than "
                                                                                "unsorted SAM output")

        return aligner_argument_parser

    def get_command_list(self, read_files_dict, **kwargs):
        extra_options = ' '.join(kwargs['extra'])

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

            additional_options = ' -S %s.sam' % output_file

            command_list.append(self.command_template.safe_substitute(kwargs, processed_input=processed_input,
                                                                      output_basename=output_file, extra=extra_options)
                                + additional_options)
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
        progress_bar = progress.bar.Bar('Processed %(index)d/%(max)d', suffix='Remaining: %(eta_td)s',
                                        max=len(command_list))

        # create output directory if it doesn't already exist; ignore for this base class
        if not kwargs['check'] and self.name != 'test':
            os.makedirs(kwargs['output'], exist_ok=True)

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
                    progress_bar.update()
                    with open(log_filename, 'w') as log:
                        aligner = subprocess.Popen(alignment_command_list, stdout=subprocess.PIPE, stderr=log)
                        samtools_bam = subprocess.Popen(samtools_bam_command, stdin=aligner.stdout,
                                                        stdout=subprocess.PIPE)
                        samtools_sort = subprocess.Popen(samtools_sort_command, stdin=samtools_bam.stdout)
                        samtools_sort.communicate()
                    self.log_list.append(log_filename)
                    progress_bar.next()

        else:
            for alignment_command in command_list:
                if kwargs['check']:
                    print(alignment_command)
                else:
                    progress_bar.update()
                    output_basename = alignment_command.split()[-1][:-4]
                    log_filename = output_basename + '.log'
                    with open(log_filename, 'w') as log:
                        subprocess.run(args=shlex.split(alignment_command), env=os.environ.copy(), stderr=log)
                    self.log_list.append(log_filename)
                    progress_bar.next()

        progress_bar.finish()
        print('\nCompleted in %s.\n' % progress_bar.elapsed_td)

        return True


class Bowtie2Aligner(SalzbergAligner):
    name = 'bowtie2'

    def __init__(self):
        super().__init__()
        self.command_template = string.Template("bowtie2 $extra -x $index -p $threads $processed_input")

    def summarize_results(self):
        if len(self.log_list) == 0:
            return

        summary_list = []
        for log_file in self.log_list:
            summary_list.append(self.parse_log(log_file))

        # output a full result file containing all values from original Log.final.out
        result_df = pandas.DataFrame(summary_list)
        result_df.to_csv('bowtie2_alignment_complete.txt', sep='\t', index=False)

        # output simplified results too
        # total reads, uniquely mapped, uniquely mapped %, multi-mapped %, unmapped %
        result_df['Unique %'] = result_df['Uniquely mapped']/result_df['Total reads']
        result_df['Multimapped %'] = result_df['Multimapped'] / result_df['Total reads']
        result_df['Unmapped %'] = result_df['Unmapped'] / result_df['Total reads']
        result_df.drop(columns=['Multimapped', 'Unmapped'], inplace=True)

        result_df.to_csv('bowtie2_alignment_summary.txt', sep='\t', index=False)

        print(result_df.to_string(index=False, formatters={'Total reads': '{:,}'.format,
                                                           'Uniquely mapped': '{:,}'.format,
                                                           'Unique %': '{:.1%}'.format,
                                                           'Multimapped %': '{:.1%}'.format,
                                                           'Unmapped %': '{:.1%}'.format}))

    @classmethod
    def parse_log(cls, filename):
        """
         parse a Bowtie2 log, returning a pandas.Series
         :return: a Series
         """
        log_entries = dict()
        log_entries['Sample'] = os.path.split(filename)[-1].replace('.log', '')

        with open(filename) as bowtie2_log:
            log_entries['Total reads'] = int(bowtie2_log.readline().strip().split()[0])
            bowtie2_log.readline()
            log_entries['Unmapped'] = int(bowtie2_log.readline().strip().split()[0])
            log_entries['Uniquely mapped'] = int(bowtie2_log.readline().strip().split()[0])
            log_entries['Multimapped'] = int(bowtie2_log.readline().strip().split()[0])

        return pandas.Series(log_entries)


class HISAT2Aligner(SalzbergAligner):
    name = 'hisat2'
    SINGLE_END_RESULT_COLUMNS = 5

    def __init__(self):
        super().__init__()
        self.command_template = string.Template("hisat2 $extra --new-summary -x $index -p $threads $processed_input")

    @classmethod
    def get_argument_parser(cls, aligner_argument_parser=None):
        aligner_argument_parser = super(HISAT2Aligner, cls).get_argument_parser(aligner_argument_parser)
        aligner_argument_parser.add_argument('--dta', action='store_true', help="add information for downstream"
                                                                                "transcriptome assembly")
        return aligner_argument_parser

    def summarize_results(self):
        if len(self.log_list) == 0:
            return

        summary_list = []
        for log_file in self.log_list:
            summary_list.append(self.parse_log(log_file))

        # output a full result file containing all values from original Log.final.out
        result_df = pandas.DataFrame(summary_list)
        result_df.to_csv('hisat2_alignment_complete.txt', sep='\t', index=False)

        # output simplified results too
        # total reads, uniquely mapped, uniquely mapped %, multi-mapped %, unmapped %

        # paired/single end results have different column names
        if result_df.shape[1] > HISAT2Aligner.SINGLE_END_RESULT_COLUMNS:
            simplified_df = (result_df.loc[:, ['Sample', 'Total pairs', 'Aligned concordantly 1 time',
                                               'Aligned concordantly >1 times',
                                               'Aligned concordantly or discordantly 0 time']]
                                      .rename(index=str,
                                              columns={'Total pairs': 'Total reads',
                                                       'Aligned concordantly 1 time': 'Uniquely mapped',
                                                       'Aligned concordantly >1 times': 'Multimapped',
                                                       'Aligned concordantly or discordantly 0 time': 'Unmapped'}))
        else:
            simplified_df = (result_df.loc[:, ['Sample', 'Total reads', 'Aligned 1 time', 'Aligned >1 times',
                                               'Aligned 0 time']]
                                      .rename(index=str,
                                              columns={'Aligned 1 time': 'Uniquely mapped',
                                                       'Aligned >1 times': 'Multimapped',
                                                       'Aligned 0 time': 'Unmapped'}))

        simplified_df['Unique %'] = simplified_df['Uniquely mapped']/simplified_df['Total reads']
        simplified_df['Multimapped %'] = simplified_df['Multimapped'] / simplified_df['Total reads']
        simplified_df['Unmapped %'] = simplified_df['Unmapped'] / simplified_df['Total reads']
        simplified_df.drop(columns=['Multimapped', 'Unmapped'], inplace=True)

        simplified_df.to_csv('hisat2_alignment_summary.txt', sep='\t', index=False)

        print(simplified_df.to_string(index=False, formatters={'Total reads': '{:,}'.format,
                                                               'Uniquely mapped': '{:,}'.format,
                                                               'Unique %': '{:.1%}'.format,
                                                               'Multimapped %': '{:.1%}'.format,
                                                               'Unmapped %': '{:.1%}'.format}))

    @classmethod
    def parse_log(cls, filename):
        """
         parse a HISAT2 log, returning a pandas.Series
         :return: a Series
         """

        log_entries = ['Sample']
        log_values = [os.path.split(filename)[-1].replace('.log', '')]

        log_df = pandas.read_csv(filename, sep=":", header=None, names=['entry', 'value'], skipinitialspace=True,
                                 skiprows=1, skipfooter=1, engine='python',
                                 converters={'entry': str.strip, 'value': lambda x: int(x.split(' ')[0])})

        log_entries += log_df.entry.tolist()
        log_values += log_df.value.tolist()

        return pandas.Series(log_values, log_entries)


class KallistoAligner(ReadAlignmentBase):
    name = 'kallisto'

    def __init__(self):
        super().__init__()
        self.command_template = string.Template("kallisto quant -i $index --bias --rf-stranded -t $threads "
                                                "-o $basename $read_files $extra")
        self.output_list = []

    def get_command_list(self, read_files_dict, **kwargs):
        """
        kallisto doesn't use comma-separated list of files
        kallisto outputs to directories,
        1. change sample output to
        2.

        :param read_files_dict:
        :param kwargs:
        :return:
        """
        command_list = []

        for sample_name in sorted(read_files_dict.keys()):
            additional_options = ''

            input_string = read_files_dict[sample_name]
            # paired-end if a space is present
            if ' ' in input_string:
                input_string = KallistoAligner.arrange_read_filenames(input_string)

            # single-end analysis requires fragment length and std. dev. to be specified
            else:
                input_string = input_string.replace(',', ' ')
                additional_options += ' --single -l 200 -s 30'

            output_dir = os.path.join(kwargs['output'], sample_name)
            self.output_list.append(output_dir)
            extra_options = ' '.join(kwargs['extra'])

            command_list.append(self.command_template.safe_substitute(kwargs, read_files=input_string,
                                                                      basename=output_dir, extra=extra_options) +
                                additional_options)
        return command_list

    def align_reads(self, **kwargs):
        """
        after running inherited method, build up a list of log files
        :param kwargs:
        :return:
        """
        super().align_reads(**kwargs)
        if not kwargs['check']:
            self.log_list = [os.path.join(output_dir, 'run_info.json') for output_dir in self.output_list]

    def summarize_results(self):
        if len(self.log_list) == 0:
            return

        summary_list = []
        for log_file in self.log_list:
            sample_name = os.path.split(os.path.dirname(log_file))[1]
            summary_series = self.parse_log(log_file)
            summary_series['Sample'] = sample_name
            summary_list.append(summary_series)

        # output a full result file containing all values from original log
        result_df = pandas.DataFrame(summary_list)
        result_df.to_csv('kallisto_quant_complete.txt', sep='\t', index=False)

        # output simplified results too
        # total reads, uniquely mapped, uniquely mapped %, multi-mapped %, unmapped %
        simplified_df = result_df.loc[:, ['Sample', 'n_processed', 'n_pseudoaligned',
                                          'n_unique', 'p_pseudoaligned', 'p_unique']]
        simplified_df.rename(index=str, columns={'n_processed': 'Total reads',
                                                 'n_unique': 'Uniquely aligned',
                                                 'p_unique': 'Unique %',
                                                 'n_pseudoaligned': 'Pseudoaligned',
                                                 'p_pseudoaligned': 'Pseudoaligned %'},
                             inplace=True)
        simplified_df['Unaligned %'] = 100 - simplified_df['Pseudoaligned %']
        simplified_df.to_csv('star_alignment_summary.txt', sep='\t', index=False)

        print(simplified_df.to_string(index=False, formatters={'Total reads': '{:,}'.format,
                                                               'Pseudoaligned': '{:,}'.format,
                                                               'Uniquely aligned': '{:,}'.format,
                                                               'Unique %': '{:.1f}'.format,
                                                               'Pseudoaligned %': '{:.1f}'.format,
                                                               'Unaligned %': '{:.1f}'.format}))

    @classmethod
    def parse_log(cls, filename):
        """
        extract useful values from the run_info.json file created by all kallisto quant runs.
        :param filename:
        :return:
        """
        return pandas.read_json(filename, typ='Series')

    @classmethod
    def arrange_read_filenames(cls, filename_string):
        """
        unlike other aligners, kallisto's input is an interleaved list of filenames
        :param filename_string: e.g., '1_1_R1.fastq,1_R1.fastq 1_1_R2.fastq,1_R2.fastq'
        :return: e.g. '1_1_R1.fastq  1_1_R2.fastq 1_R1.fastq 1_R2.fastq'
        """
        r1, r2 = filename_string.split(' ')
        interleaved_string = ' '.join([' '.join(x) for x in zip(r1.split(','), r2.split(','))])
        return interleaved_string


def run_aligner(namespace):
    aligner = namespace.Aligner()
    aligner.align_reads(**vars(namespace))
    aligner.process_alignments()
    aligner.summarize_results()


def setup_subparsers(subparsers):
    aligner_parser = subparsers.add_parser('align', help='Align reads to reference genome/transcriptome')
    aligner_subparser = aligner_parser.add_subparsers(title="aligner", dest="aligner",
                                                      description="available alignment programs")

    for aligner in [ReadAlignmentBase, STARAligner, HISAT2Aligner, Bowtie2Aligner, KallistoAligner]:
        current_subparser = aligner_subparser.add_parser(aligner.name)
        aligner.get_argument_parser(current_subparser)

    aligner_parser.set_defaults(func=run_aligner)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)
