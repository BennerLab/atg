import unittest
import pandas
import atg.data.retrieve
import atg.data.identifiers

human_symbol_series = pandas.Series(['CXCL9', 'ACTB'])
human_ensembl_series = pandas.Series(['ENSG00000138755', 'ENSG00000075624'])

mouse_symbol_series = pandas.Series(['Cxcl9', 'Actb'])
mouse_ensembl_series = pandas.Series(['ENSMUSG00000029417', 'ENSMUSG00000029580'])


class IDConversionTest(unittest.TestCase):
    def setUp(self):
        # make sure that relevant files are available
        data_tracker = atg.data.retrieve.ATGDataTracker()
        data_tracker.retrieve_data("human", selected_files=['ensembl_gene.csv'])
        data_tracker.retrieve_data("mouse", selected_files=['ensembl_gene.csv'])

        self.human_translator = atg.data.identifiers.GeneIDTranslator('human')
        self.mouse_translator = atg.data.identifiers.GeneIDTranslator('mouse')

        with self.assertRaises(ValueError):
            atg.data.identifiers.GeneIDTranslator('robot')

    def test_translate_symbol_ensembl(self):
        """
        test a couple of one-to-one translations both ways
        :return:
        """
        for symbol_series, ensembl_series, translator in [(human_symbol_series, human_ensembl_series, self.human_translator),
                                              (mouse_symbol_series, mouse_ensembl_series, self.mouse_translator)]:

            ensembl_result = translator.translate_identifiers(symbol_series, input_type=None, output_type='ensembl')
            pandas.util.testing.assert_series_equal(ensembl_series, ensembl_result, check_names=False)

            ensembl_result = translator.translate_identifiers(ensembl_series, input_type=None, output_type='symbol')
            pandas.util.testing.assert_series_equal(symbol_series, ensembl_result, check_names=False)

    def test_translation_fails(self):
        with self.assertRaises(ValueError):
            self.human_translator.translate_identifiers([], input_type='nothing', output_type="something")
        with self.assertRaises(ValueError):
            self.human_translator.translate_identifiers([], input_type=None, output_type="something")

    def test_guessing(self):
        self.assertEqual('ensembl', atg.data.identifiers.guess_identifier_type(human_ensembl_series))
        self.assertEqual('symbol', atg.data.identifiers.guess_identifier_type(human_symbol_series))
        self.assertEqual('entrez', atg.data.identifiers.guess_identifier_type([1, 2, 3]))
        self.assertEqual('entrez', atg.data.identifiers.guess_identifier_type(['1', '2', '3']))

        self.assertEqual('symbol', atg.data.identifiers.guess_identifier_type(['1', 1, 'a']))
