import pytest
import os
import sys
import shlex
from atg.util import align

BOWTIE2_LOG = os.path.join(os.path.dirname(__file__), 'data', 'bowtie2_SE.log')
KALLISTO_LOG = os.path.join(os.path.dirname(__file__), 'data', 'kallisto_run.json')


class TestSTARCommand:
    def test_star_arg_parsing(self):
        arguments = shlex.split(' -k -t 16 /directory/genomeIndex output_directory a.fastq b.fastq')
        parser = align.STARAligner.get_argument_parser()
        star_argument_namespace = parser.parse_args(arguments)

        assert star_argument_namespace.index == '/directory/genomeIndex'
        assert star_argument_namespace.fastq == ['a.fastq', 'b.fastq']
        assert star_argument_namespace.keep_unmapped
        assert star_argument_namespace.threads == 16

    def test_star_command_construction(self, capsys):
        correct_output = 'STAR --genomeDir /directory/genomeIndex --runThreadN 16 ' \
                         '--readFilesIn a_1.fastq --outFileNamePrefix output_directory/a. --outSAMtype BAM ' \
                         'SortedByCoordinate  --genomeLoad LoadAndKeep --limitBAMsortRAM 50000000000 ' \
                         '--outBAMsortingBinsN 31\n' \
                         'STAR --genomeDir /directory/genomeIndex --runThreadN 16 ' \
                         '--readFilesIn b_1.fastq --outFileNamePrefix output_directory/b. --outSAMtype BAM ' \
                         'SortedByCoordinate  --genomeLoad LoadAndKeep --limitBAMsortRAM 50000000000 ' \
                         '--outBAMsortingBinsN 31'

        arguments = shlex.split('-t 16 -c /directory/genomeIndex output_directory a_1.fastq b_1.fastq')

        parser = align.STARAligner.get_argument_parser()
        args = parser.parse_args(arguments)
        aligner = align.STARAligner()
        aligner.align_reads(**vars(args))

        captured = capsys.readouterr()
        aligner_output = captured.out.strip()
        assert aligner_output == correct_output

    def test_star_command_construction_gzip(self, capsys):
        # self.maxDiff = None
        correct_output = 'STAR --genomeDir /directory/genomeIndex --runThreadN 16 ' \
                         '--readFilesIn a_1.fastq.gz --outFileNamePrefix output_directory/a. --outSAMtype BAM ' \
                         'SortedByCoordinate  --readFilesCommand gunzip -c --genomeLoad LoadAndKeep ' \
                         '--limitBAMsortRAM 50000000000 --outBAMsortingBinsN 31\n' \
                         'STAR --genomeDir /directory/genomeIndex --runThreadN 16 ' \
                         '--readFilesIn b_1.fastq.gz --outFileNamePrefix output_directory/b. --outSAMtype BAM ' \
                         'SortedByCoordinate  --readFilesCommand gunzip -c --genomeLoad LoadAndKeep ' \
                         '--limitBAMsortRAM 50000000000 --outBAMsortingBinsN 31'

        arguments = shlex.split('-c -t 16 /directory/genomeIndex output_directory a_1.fastq.gz b_1.fastq.gz')

        parser = align.STARAligner.get_argument_parser()
        args = parser.parse_args(arguments)
        aligner = align.STARAligner()
        aligner.align_reads(**vars(args))

        captured = capsys.readouterr()
        aligner_output = captured.out.strip()
        assert aligner_output == correct_output


class TestBowtie2Command:
    def test_log_parsing(self):
        log_series = align.Bowtie2Aligner.parse_log(BOWTIE2_LOG)
        assert log_series['Total reads'] == 12494316
        assert log_series['Unmapped'] == 72860
        assert log_series['Uniquely mapped'] == 7693710
        assert log_series['Multimapped'] == 4727746
        assert log_series[2:5].sum() == log_series['Total reads']


class TestKallistoCommand:
    def test_paired_end(self, capsys):
        correct_output = 'kallisto quant -i /directory/transcriptome.tdx --bias --rf-stranded -t 8 ' \
                         '-o output_directory/a a_1_R1.fastq.gz a_1_R2.fastq.gz a_2_R1.fastq.gz a_2_R2.fastq.gz'

        arguments = shlex.split('-c -t 8 /directory/transcriptome.tdx output_directory a_1_R1.fastq.gz a_1_R2.fastq.gz '
                                'a_2_R1.fastq.gz a_2_R2.fastq.gz')

        parser = align.KallistoAligner.get_argument_parser()
        args = parser.parse_args(arguments)
        aligner = align.KallistoAligner()
        aligner.align_reads(**vars(args))

        captured = capsys.readouterr()
        aligner_output = captured.out.strip()
        assert aligner_output == correct_output

    def test_single_end(self, capsys):
        correct_output = 'kallisto quant -i /directory/transcriptome.tdx --bias --rf-stranded -t 8 ' \
                         '-o output_directory/a a_1_R1.fastq.gz a_2_R1.fastq.gz  --single -l 200 -s 30'

        arguments = shlex.split('-c -t 8 /directory/transcriptome.tdx output_directory a_1_R1.fastq.gz '
                                'a_2_R1.fastq.gz')

        parser = align.KallistoAligner.get_argument_parser()
        args = parser.parse_args(arguments)
        aligner = align.KallistoAligner()
        aligner.align_reads(**vars(args))

        captured = capsys.readouterr()
        aligner_output = captured.out.strip()
        assert aligner_output == correct_output

    def test_log_parsing(self):
        log_series = align.KallistoAligner.parse_log(KALLISTO_LOG)
        assert log_series['n_unique'] == 16090035
        assert log_series['n_processed'] == 44891704
        assert log_series['n_pseudoaligned'] == 36521130
