import unittest
import os
import atg.stats.enrich
import numpy

# Genes present in GO:0007259
SAMPLE_ENSEMBL_GENE_LIST = ['ENSG00000159110', 'ENSG00000185338', 'ENSG00000183709', 'ENSG00000182393',
                            'ENSG00000126561', 'ENSG00000197110', 'ENSG00000170677', 'ENSG00000134470',
                            'ENSG00000008710', 'ENSG00000124356', 'ENSG00000113494', 'ENSG00000123609',
                            'ENSG00000120833', 'ENSG00000096968', 'ENSG00000168610', 'ENSG00000198730',
                            'ENSG00000115415', 'ENSG00000112964', 'ENSG00000121807', 'ENSG00000170581',
                            'ENSG00000173757', 'ENSG00000068078', 'ENSG00000175505', 'ENSG00000171150',
                            'ENSG00000142166', 'ENSG00000108691', 'ENSG00000259384', 'ENSG00000164509',
                            'ENSG00000138378', 'ENSG00000033800', 'ENSG00000118762']

SAMPLE_ENSEMBL_GENE_LIST2 = ['ENSG00000196365', 'ENSG00000140451', 'ENSG00000198836', 'ENSG00000117020',
                             'ENSG00000068305', 'ENSG00000068305', 'ENSG00000275199', 'ENSG00000151729',
                             'ENSG00000114120', 'ENSG00000154719']


def test_adjust_pvalue():
    sample_p_values = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216, 0.222, 0.251, 0.269,
                       0.275, 0.34, 0.341, 0.384, 0.569, 0.594, 0.696, 0.762, 0.94, 0.942, 0.975, 0.986]

    bh_corrected_p_values = [0.025, 0.1, 0.21, 0.21, 0.21, 0.25, 0.2643, 0.4911, 0.4911, 0.4911, 0.4911, 0.4911,
                             0.4911, 0.4911, 0.5328, 0.5328, 0.5647, 0.7816, 0.7816, 0.87, 0.9071, 0.986, 0.986,
                             0.986, 0.986]
    numpy.testing.assert_almost_equal(atg.stats.enrich.p_adjust_bh(sample_p_values), bh_corrected_p_values,
                                      decimal=4)

    bh_corrected_p_values_n100 = atg.stats.enrich.p_adjust_bh(sample_p_values, n=100)
    numpy.testing.assert_array_less(bh_corrected_p_values, bh_corrected_p_values_n100)

    bh_corrected_p_values_n25 = atg.stats.enrich.p_adjust_bh(sample_p_values, n=25)
    numpy.testing.assert_almost_equal(bh_corrected_p_values, bh_corrected_p_values_n25, decimal=4)


class EnrichmentTest(unittest.TestCase):
    def setUp(self):
        data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
        go_term_path = os.path.join(data_root, 'human', 'GRCh38', 'gene_go.csv')
        go_definition_path = os.path.join(data_root, 'human', 'GRCh38', 'go_definition.csv')
        self.calculator = atg.stats.enrich.EnrichmentCalculator(go_term_path, go_definition_path)
        self.all_genes = self.calculator.gene_term_df.ix[:, atg.stats.enrich.GENE_COLUMN_LABEL_INDEX].drop_duplicates()

    def test_single_term(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:5] + self.all_genes.sample(100).tolist()
        good_pvalue = self.calculator.get_single_enrichment(good_gene_list, 'GO:0007259')

        better_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:10] + self.all_genes.sample(100).tolist()
        better_pvalue = self.calculator.get_single_enrichment(better_gene_list, 'GO:0007259')
        self.assertLess(better_pvalue, good_pvalue)
        self.assertLess(good_pvalue, -5.0)

        bad_pvalue = self.calculator.get_single_enrichment(good_gene_list, 'GO:0000122')
        self.assertGreater(bad_pvalue, -5.0)

    def test_full_enrichment(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:5] + self.all_genes.sample(100).tolist()

        enrichment_df = self.calculator.get_all_enrichment(good_gene_list)
        self.assertIn('GO:0007259', enrichment_df.index)
        self.assertLess(sum(enrichment_df['log_pvalue'] < -5), 100)

    def test_plot(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:15] + SAMPLE_ENSEMBL_GENE_LIST2
        self.calculator.plot_enrichment_single(good_gene_list)

    def test_iterative_enrichment(self):
        self.calculator.iterative_enrichment(SAMPLE_ENSEMBL_GENE_LIST + SAMPLE_ENSEMBL_GENE_LIST2 +
                                             self.all_genes.sample(100).tolist())
