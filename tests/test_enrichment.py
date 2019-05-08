import pytest
import os
import atg.stats.enrich
import atg.data.retrieve
import numpy
import pandas
import argparse
import shlex

# Genes present in GO:0007259 (JAK-STAT cascade)
SAMPLE_ENSEMBL_GENE_LIST = ['ENSG00000159110', 'ENSG00000185338', 'ENSG00000183709', 'ENSG00000182393',
                            'ENSG00000126561', 'ENSG00000197110', 'ENSG00000170677', 'ENSG00000134470',
                            'ENSG00000008710', 'ENSG00000124356', 'ENSG00000113494', 'ENSG00000123609',
                            'ENSG00000120833', 'ENSG00000096968', 'ENSG00000168610', 'ENSG00000198730',
                            'ENSG00000115415', 'ENSG00000112964', 'ENSG00000121807', 'ENSG00000170581',
                            'ENSG00000173757', 'ENSG00000068078', 'ENSG00000175505', 'ENSG00000171150',
                            'ENSG00000142166', 'ENSG00000108691', 'ENSG00000259384', 'ENSG00000164509',
                            'ENSG00000138378', 'ENSG00000033800', 'ENSG00000118762']

# Genes present in GO:0000002 (mitochondrial genome maintenance)
SAMPLE_ENSEMBL_GENE_LIST2 = ['ENSG00000196365', 'ENSG00000140451', 'ENSG00000198836', 'ENSG00000117020',
                             'ENSG00000068305', 'ENSG00000068305', 'ENSG00000275199', 'ENSG00000151729',
                             'ENSG00000114120', 'ENSG00000154719']

# Genes present in GO:0050821 (protein stabilization)
SAMPLE_ENSEMBL_GENE_LIST3 = ['ENSG00000116288', 'ENSG00000124766', 'ENSG00000071537', 'ENSG00000107937',
                             'ENSG00000049246', 'ENSG00000134480', 'ENSG00000164134', 'ENSG00000224200',
                             'ENSG00000232804', 'ENSG00000124587', 'ENSG00000124762', 'ENSG00000115977',
                             'ENSG00000115484', 'ENSG00000132589', 'ENSG00000237724', 'ENSG00000162407',
                             'ENSG00000137312', 'ENSG00000080824', 'ENSG00000231555', 'ENSG00000108518',
                             'ENSG00000173928', 'ENSG00000005893', 'ENSG00000111640', 'ENSG00000148584',
                             'ENSG00000235941', 'ENSG00000180228', 'ENSG00000224501', 'ENSG00000169877',
                             'ENSG00000151929', 'ENSG00000171209']

GENE_LIST_LOCATION = os.path.join(os.path.dirname(__file__), 'data', 'genelists')
HALLMARK_APOPTOSIS = os.path.join(GENE_LIST_LOCATION, 'HALLMARK_APOPTOSIS.txt')
HALLMARK_P53_PATHWAY = os.path.join(GENE_LIST_LOCATION, 'HALLMARK_P53_PATHWAY.txt')


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


def test_enrichment_calculation():
    term_row = {"hit_count": 4,
                "universe": 18466,
                "term_count": 44,
                "list_size": 135}

    numpy.testing.assert_almost_equal(atg.stats.enrich.enrichment_significance(term_row), -3.52925402)


class TestEnrichment:
    def setup(self):
        data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
        go_term_path = os.path.join(data_root, 'human', 'GRCh38', 'go_biological_process.csv')
        go_definition_path = os.path.join(data_root, 'human', 'GRCh38', 'go_definition.csv')

        data_tracker = atg.data.retrieve.ATGDataTracker()
        data_tracker.retrieve_data("human", selected_files=['ensembl_gene.csv'])

        self.calculator = atg.stats.enrich.EnrichmentCalculator(go_term_path, go_definition_path)
        self.all_genes = self.calculator.gene_term_df.iloc[:, atg.stats.enrich.GENE_COLUMN_LABEL_INDEX].\
            drop_duplicates()

    def test_single_term(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:5] + self.all_genes.sample(100).tolist()
        good_pvalue = self.calculator.get_single_enrichment(good_gene_list, 'GO:0007259')

        better_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:10] + self.all_genes.iloc[0:100].tolist()
        better_pvalue = self.calculator.get_single_enrichment(better_gene_list, 'GO:0007259')
        assert better_pvalue < good_pvalue
        assert good_pvalue < -5.0

        bad_pvalue = self.calculator.get_single_enrichment(SAMPLE_ENSEMBL_GENE_LIST3, 'GO:0006367')
        assert bad_pvalue > -5.0

    def test_full_enrichment(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:5] + self.all_genes.sample(100).tolist()

        enrichment_df = self.calculator.get_all_enrichment(good_gene_list)
        assert 'GO:0007259' in enrichment_df.index
        assert sum(enrichment_df['log_pvalue'] < -5) < 100

    def test_plot(self):
        good_gene_list = SAMPLE_ENSEMBL_GENE_LIST[0:15] + SAMPLE_ENSEMBL_GENE_LIST2
        self.calculator.plot_enrichment_single(good_gene_list, output_filename='/tmp/enrich.pdf')
        self.calculator.plot_enrichment_single(good_gene_list, iterative=True,
                                               output_filename='/tmp/iterative_enrichment.pdf')

    def test_iterative_enrichment(self):
        iterative_result = self.calculator.iterative_enrichment(SAMPLE_ENSEMBL_GENE_LIST + SAMPLE_ENSEMBL_GENE_LIST2 +
                                                                self.all_genes.loc[0:500].tolist())
        assert 'GO:0007259' in iterative_result.index
        assert 'GO:0000002' in iterative_result.index

    def test_iterative_enrichment_multilist(self):
        multi_gene_list = {'A': SAMPLE_ENSEMBL_GENE_LIST, 'B': SAMPLE_ENSEMBL_GENE_LIST2,
                           'C': SAMPLE_ENSEMBL_GENE_LIST3}
        multi_enrichment_result = self.calculator.iterative_enrichment_multilist(multi_gene_list)
        assert multi_enrichment_result.index.equals(pandas.Index(['GO:0007259', 'GO:0050821', 'GO:0000002']))

    def test_plot_multiple_enrichment(self):
        multi_gene_list = {'A': SAMPLE_ENSEMBL_GENE_LIST, 'B': SAMPLE_ENSEMBL_GENE_LIST2,
                           'C': SAMPLE_ENSEMBL_GENE_LIST3}
        self.calculator.plot_enrichment_multiple(multi_gene_list, output_filename="/tmp/multi.pdf")
        assert True

    def test_namespace_run(self):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(dest="command", help='commands')
        atg.stats.enrich.setup_subparsers(subparsers)

        bad_species_command = shlex.split('enrich test.txt -s robot')

        with pytest.raises(ValueError):
            atg.stats.enrich.run_enrichment(parser.parse_args(bad_species_command))

        atg.stats.enrich.run_enrichment(parser.parse_args(['enrich', HALLMARK_APOPTOSIS]))
        atg.stats.enrich.run_enrichment(parser.parse_args(['enrich', HALLMARK_APOPTOSIS, HALLMARK_P53_PATHWAY]))
