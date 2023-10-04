import os
import tempfile
import csv

from pybedtools.helpers import chromsizes_to_file

from atg.quantification import coverage

UNSTRANDED_SE_RNASEQ_FILE = os.path.join(os.path.dirname(__file__), 'data', 'rnaseq_unstranded_SE.bam')
STRANDED_PE_RNASEQ_FILE = os.path.join(os.path.dirname(__file__), 'data', 'rnaseq.bam')
CHROM_SIZES_HG38 = os.path.join(os.path.dirname(__file__), 'data', 'hg38_chrom.sizes')


class TestUnstrandedCoverageCalculator:
    def setup_method(self):
        self.coverage_calculator = coverage.UnstrandedCoverageCalculator(UNSTRANDED_SE_RNASEQ_FILE)

    def test_read_count(self):
        assert self.coverage_calculator.mapped_read_count == 1224

    def test_bedgraph(self):
        """
        all entries should be on chromosome 11
        forward file should have only positive entries in column 4, reverse file should have only negative
        """
        with tempfile.TemporaryDirectory() as working_directory:
            forward_filename = os.path.join(working_directory, "forward.bedgraph")
            self.coverage_calculator.write_bedgraph(forward_filename)

            for entry in csv.reader(open(forward_filename), delimiter='\t'):
                assert entry[0] == 'chr11'
                assert float(entry[3]) > 0

    def test_bigwig_with_genome(self):
        """
        just check that files were created
        TODO: use checksum?
        """
        with tempfile.TemporaryDirectory() as working_directory:
            forward_filename = os.path.join(working_directory, "forward.bigwig")
            self.coverage_calculator.write_bigwig(forward_filename, 'hg38')

            assert os.path.exists(forward_filename)

    def test_bigwig_with_chrom_sizes(self):
        """
        just check that files were created
        TODO: use checksum?
        """
        with tempfile.TemporaryDirectory() as working_directory:
            forward_filename = os.path.join(working_directory, "forward.bigwig")
            self.coverage_calculator.write_bigwig(forward_filename, chrom_sizes=CHROM_SIZES_HG38)

            assert os.path.exists(forward_filename)


class TestStrandedCoverageCalculator:
    def setup_method(self):
        self.coverage_calculator = coverage.StrandedCoverageCalculator(STRANDED_PE_RNASEQ_FILE)

    def test_read_count(self):
        """
        test file contains 198 uniquely mapped R1 reads
        """
        assert self.coverage_calculator.mapped_read_count == 196

    def test_bedgraph(self):
        """
        all entries should be on chromosome 11
        forward file should have only positive entries in column 4, reverse file should have only negative
        """
        with tempfile.TemporaryDirectory() as working_directory:
            forward_filename = os.path.join(working_directory, "forward.bedgraph")
            reverse_filename = os.path.join(working_directory, "reverse.bedgraph")
            self.coverage_calculator.write_bedgraph(forward_filename, reverse_filename)

            for entry in csv.reader(open(forward_filename), delimiter='\t'):
                assert entry[0] == 'chr11'
                assert float(entry[3]) > 0

            for entry in csv.reader(open(reverse_filename), delimiter='\t'):
                assert entry[0] == 'chr11'
                assert float(entry[3]) < 0

    def test_bigwig_with_genome(self):
        """
        just check that files were created
        TODO: use checksum?
        """
        with tempfile.TemporaryDirectory() as working_directory:
            forward_filename = os.path.join(working_directory, "forward.bigwig")
            reverse_filename = os.path.join(working_directory, "reverse.bigwig")
            self.coverage_calculator.write_bigwig(forward_filename, reverse_filename, 'hg38')

            assert os.path.exists(forward_filename)
            assert os.path.exists(reverse_filename)

    def test_bigwig_with_chrom_sizes(self):
        """
        just check that files were created
        TODO: use checksum?
        """
        with tempfile.TemporaryDirectory() as working_directory:
            forward_filename = os.path.join(working_directory, "forward.bigwig")
            reverse_filename = os.path.join(working_directory, "reverse.bigwig")
            self.coverage_calculator.write_bigwig(forward_filename, reverse_filename, chrom_sizes=CHROM_SIZES_HG38)

            assert os.path.exists(forward_filename)
            assert os.path.exists(reverse_filename)
