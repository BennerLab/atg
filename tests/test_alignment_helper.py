import unittest
import sys
import shlex
from atg.util import align


class MyTestCase(unittest.TestCase):
    def test_star_arg_parsing(self):
        arguments = shlex.split('/directory/genomeIndex output_directory a.fastq b.fastq -k -t 16')
        parser = align.STARAligner.get_argument_parser()
        star_argument_namespace = parser.parse_args(arguments)

        self.assertEqual(star_argument_namespace.index, '/directory/genomeIndex')
        self.assertListEqual(star_argument_namespace.fastq, ['a.fastq', 'b.fastq'])
        self.assertEqual(star_argument_namespace.keep_unmapped, True)
        self.assertEqual(star_argument_namespace.threads, 16)

    def test_star_command_construction(self):
        CORRECT_OUTPUT = 'STAR --genomeDir /directory/genomeIndex --genomeLoad LoadAndKeep --runThreadN 16 ' \
                         '--readFilesIn a_1.fastq --outFileNamePrefix output_directory/a. --outSAMtype BAM ' \
                         'SortedByCoordinate --limitBAMsortRAM 50000000000\n' \
                         'STAR --genomeDir /directory/genomeIndex --genomeLoad LoadAndKeep --runThreadN 16 ' \
                         '--readFilesIn b_1.fastq --outFileNamePrefix output_directory/b. --outSAMtype BAM ' \
                         'SortedByCoordinate --limitBAMsortRAM 50000000000'

        arguments = shlex.split('/directory/genomeIndex output_directory a_1.fastq b_1.fastq -k -t 16 -c')

        parser = align.STARAligner.get_argument_parser()
        args = parser.parse_args(arguments)
        aligner = align.STARAligner()
        aligner.align_reads(**vars(args))

        self.assertEqual(sys.stdout.getvalue().strip(), CORRECT_OUTPUT)

    def test_star_command_construction_gzip(self):
        # self.maxDiff = None
        CORRECT_OUTPUT = 'STAR --genomeDir /directory/genomeIndex --genomeLoad LoadAndKeep --runThreadN 16 ' \
                         '--readFilesIn a_1.fastq.gz --outFileNamePrefix output_directory/a. --outSAMtype BAM ' \
                         'SortedByCoordinate --limitBAMsortRAM 50000000000 --readFilesCommand gunzip -c\n' \
                         'STAR --genomeDir /directory/genomeIndex --genomeLoad LoadAndKeep --runThreadN 16 ' \
                         '--readFilesIn b_1.fastq.gz --outFileNamePrefix output_directory/b. --outSAMtype BAM ' \
                         'SortedByCoordinate --limitBAMsortRAM 50000000000 --readFilesCommand gunzip -c'

        arguments = shlex.split('/directory/genomeIndex output_directory a_1.fastq.gz b_1.fastq.gz -k -t 16 -c')

        parser = align.STARAligner.get_argument_parser()
        args = parser.parse_args(arguments)
        aligner = align.STARAligner()
        aligner.align_reads(**vars(args))

        self.assertEqual(sys.stdout.getvalue().strip(), CORRECT_OUTPUT)


if __name__ == '__main__':
    unittest.main()
