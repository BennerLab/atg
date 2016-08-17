import unittest
import os
import atg.stats.enrich

# Genes present in GO:0007259
SAMPLE_ENSEMBL_GENE_LIST = ['ENSG00000159110', 'ENSG00000185338', 'ENSG00000183709', 'ENSG00000182393',
                            'ENSG00000126561', 'ENSG00000197110', 'ENSG00000170677', 'ENSG00000134470',
                            'ENSG00000008710', 'ENSG00000124356', 'ENSG00000113494', 'ENSG00000123609',
                            'ENSG00000120833', 'ENSG00000096968', 'ENSG00000168610', 'ENSG00000198730',
                            'ENSG00000115415', 'ENSG00000112964', 'ENSG00000121807', 'ENSG00000170581',
                            'ENSG00000173757', 'ENSG00000068078', 'ENSG00000175505', 'ENSG00000171150',
                            'ENSG00000142166', 'ENSG00000108691', 'ENSG00000259384', 'ENSG00000164509',
                            'ENSG00000138378', 'ENSG00000033800', 'ENSG00000118762']


class EnrichmentTest(unittest.TestCase):
    def setUp(self):
        data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
        go_term_path = os.path.join(data_root, 'human', 'Current', 'hg38', 'gene_go.csv')
        self.calculator = atg.stats.enrich.EnrichmentCalculator(go_term_path)

    def test_single_term(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:5] + list(range(100))
        good_pvalue = self.calculator.get_single_enrichment(good_gene_list, 'GO:0007259')

        better_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:10] + list(range(100))
        better_pvalue = self.calculator.get_single_enrichment(better_gene_list, 'GO:0007259')
        self.assertLess(better_pvalue, good_pvalue)
        self.assertLess(good_pvalue, -5.0)

        bad_pvalue = self.calculator.get_single_enrichment(good_gene_list, 'GO:0000122')
        self.assertGreater(bad_pvalue, -5.0)

    def test_full_enrichment(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:5] + list(range(100))

        enrichment_df = self.calculator.get_all_enrichment(good_gene_list)
        self.assertIn('GO:0007259', enrichment_df.index)
        self.assertLess(enrichment_df.shape[0], 100)
