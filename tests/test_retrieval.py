import os
import unittest
import tempfile
import pandas
import atg.data.retrieve

UCSC_GZIPPED_FILE = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.agp.gz'
UCSC_CHROMOSOME_SIZE_FILE = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
ASSEMBLY_FILE_ENTRIES_PER_LINE = 9
UCSC_CHROMOSOME_ENTRIES = 455
BIOMART_XML_FILE = os.path.join(os.path.dirname(__file__), 'data', 'GRCh38_chrY.xml')


class RetrievalTest(unittest.TestCase):
    def test_xml_fetch(self):
        with open(BIOMART_XML_FILE, 'r') as biomart_input:
            with tempfile.TemporaryDirectory() as working_directory:
                biomart_output = os.path.join(working_directory, "gene.csv")
                xml_string = ''.join(biomart_input.readlines())
                atg.data.retrieve.fetch_ensembl(xml_string, biomart_output)

                gene_df = pandas.read_csv(biomart_output)
                # at least 10 genes should be present
                self.assertGreater(gene_df.shape[0], 10)
                # the XML file restricts retrievals to chromosome Y
                self.assertEqual(gene_df['Chromosome/scaffold name'][0], 'Y')

    def test_url_fetch(self):
        with tempfile.TemporaryDirectory() as working_directory:
            assembly_filename = os.path.join(working_directory, "test.agp")
            atg.data.retrieve.fetch_url(UCSC_GZIPPED_FILE, assembly_filename)

            with open(assembly_filename, 'r') as assembly:
                firstline = assembly.readline().split()
                self.assertEqual(firstline[0], 'chr1')
                self.assertNotEqual(firstline[0], 'chr2')
                self.assertEqual(len(firstline), ASSEMBLY_FILE_ENTRIES_PER_LINE)

            chromosome_filename = os.path.join(working_directory, "chrom.sizes")
            atg.data.retrieve.fetch_url(UCSC_CHROMOSOME_SIZE_FILE, chromosome_filename)
            chromosome_df = pandas.read_table(chromosome_filename, names=['chrom', 'size'])
            self.assertEqual(chromosome_df['chrom'][0], 'chr1')
            self.assertEqual(chromosome_df.shape[0], UCSC_CHROMOSOME_ENTRIES)
