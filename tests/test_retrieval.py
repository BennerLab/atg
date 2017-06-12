import os
import unittest
import tempfile
import pandas
import atg.data.retrieve
import atg.data.ensembl

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


class EnsemblGenomesTest(unittest.TestCase):
    def setUp(self):
        self.ensembl_genomes = atg.data.ensembl.EnsemblSpecies()

    def test_ensembl_table_read(self):
        corn_record = self.ensembl_genomes.ensembl_species_df.ix[self.ensembl_genomes.ensembl_species_df['species'] ==
                                                                 'zea_mays']
        self.assertEqual(corn_record.shape[0], 1)
        self.assertEqual(corn_record.iloc[0]['division'], 'EnsemblPlants')

    def test_species_information(self):
        bee_information = {'annotation': 'pub/current/metazoa/gtf/apis_mellifera/Apis_mellifera.'
                                         'GCA_000002195.1.35.gtf.gz',
                           'version': 'GCA_000002195.1',
                           'genome': 'pub/current/metazoa/fasta/apis_mellifera/dna/Apis_mellifera.'
                                     'GCA_000002195.1.dna.toplevel.fa.gz'}

        corn_information = {'annotation': 'pub/current/plants/gtf/zea_mays/Zea_mays.AGPv4.35.gtf.gz',
                            'version': 'AGPv4',
                            'genome': 'pub/current/plants/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz'}

        mushroom_information = {'genome': 'pub/current/fungi/fasta/fungi_basidiomycota1_collection/'
                                          'agaricus_bisporus_var_bisporus_h97/dna/Agaricus_bisporus_var_bisporus_h97.'
                                          'Agabi_varbisH97_2.dna.toplevel.fa.gz',
                                'version': 'Agabi_varbisH97_2',
                                'annotation': 'pub/current/fungi/gtf/fungi_basidiomycota1_collection/'
                                              'agaricus_bisporus_var_bisporus_h97/Agaricus_bisporus_var_bisporus_h97.'
                                              'Agabi_varbisH97_2.35.gtf.gz'}

        self.assertDictEqual(self.ensembl_genomes.get_species_information('zea_mays'), corn_information)
        self.assertDictEqual(self.ensembl_genomes.get_species_information('apis_mellifera'), bee_information)
        self.assertDictEqual(self.ensembl_genomes.get_species_information('agaricus_bisporus_var_bisporus_h97'),
                             mushroom_information)