import unittest
import sys
import shlex
from atg.util import align
from atg.util.align import Bowtie2Aligner


class STARCommandTestCase(unittest.TestCase):
    def test_star_arg_parsing(self):
        arguments = shlex.split(' -k -t 16 /directory/genomeIndex output_directory a.fastq b.fastq')
        parser = align.STARAligner.get_argument_parser()
        star_argument_namespace = parser.parse_args(arguments)

        self.assertEqual(star_argument_namespace.index, '/directory/genomeIndex')
        self.assertListEqual(star_argument_namespace.fastq, ['a.fastq', 'b.fastq'])
        self.assertEqual(star_argument_namespace.keep_unmapped, True)
        self.assertEqual(star_argument_namespace.threads, 16)

    def test_star_command_construction(self):
        CORRECT_OUTPUT = 'STAR --genomeDir /directory/genomeIndex --runThreadN 16 ' \
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

        aligner_output = sys.stdout.getvalue().strip()
        self.assertEqual(aligner_output, CORRECT_OUTPUT)

    def test_star_command_construction_gzip(self):
        # self.maxDiff = None
        CORRECT_OUTPUT = 'STAR --genomeDir /directory/genomeIndex --runThreadN 16 ' \
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

        aligner_output = sys.stdout.getvalue().strip()
        self.assertEqual(aligner_output, CORRECT_OUTPUT)


class Bowtie2CommandTestCase(unittest.TestCase):
    def test_log_parsing(self):
        log_series = Bowtie2Aligner.parse_log('data/bowtie2_log_SE.txt')
        self.assertEqual(log_series['Total reads'], 12494316)
        self.assertEqual(log_series['Unmapped'], 72860)
        self.assertEqual(log_series['Uniquely mapped'], 7693710)
        self.assertEqual(log_series['Multimapped'], 4727746)
        self.assertEqual(log_series[2:5].sum(), log_series['Total reads'])

    def test_output(self):
        bowtie2_aligner = Bowtie2Aligner()
        bowtie2_aligner.log_list = ['data/bowtie2_log_SE.txt']
        bowtie2_aligner.summarize_results()


if __name__ == '__main__':
    unittest.main()
