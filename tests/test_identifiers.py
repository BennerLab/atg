import pytest
import pandas
import atg.data.retrieve
import atg.data.identifiers
from pandas.testing import assert_frame_equal, assert_series_equal

human_symbol = pandas.Series(['CXCL9', 'ACTB'])
human_ensembl = pandas.Series(['ENSG00000138755', 'ENSG00000075624'])

mouse_symbol = pandas.Series(['Cxcl9', 'Actb'])
mouse_ensembl = pandas.Series(['ENSMUSG00000029417', 'ENSMUSG00000029580'])

# genes that appear multiple times in the genome (due to presence on haplotype blocks?)
human_symbol_multiple = pandas.Series(['CCL16', 'RAD17'])
mouse_symbol_multiple = pandas.Series(['Ccl19', 'Gatm'])
human_ensembl_first = pandas.Series(['ENSG00000275152', 'ENSG00000152942'])
mouse_ensembl_first = pandas.Series(['ENSMUSG00000071005', 'ENSMUSG00000027199'])


class TestIDConversion:
    def setup(self):
        # make sure that relevant files are available
        data_tracker = atg.data.retrieve.ATGDataTracker()
        data_tracker.retrieve_data("human", selected_files=['ensembl_gene.csv'])
        data_tracker.retrieve_data("mouse", selected_files=['ensembl_gene.csv'])

        self.human_translator = atg.data.identifiers.GeneIDTranslator('human')
        self.mouse_translator = atg.data.identifiers.GeneIDTranslator('mouse')

        with pytest.raises(ValueError):
            atg.data.identifiers.GeneIDTranslator('robot')

    def test_translate_symbol_ensembl(self):
        """
        test a couple of one-to-one translations both ways
        :return:
        """
        for symbol_series, ensembl_series, translator in [(human_symbol, human_ensembl, self.human_translator),
                                                          (mouse_symbol, mouse_ensembl, self.mouse_translator)]:

            ensembl_result = translator.translate_identifiers(symbol_series, input_type=None, output_type='ensembl')
            assert_series_equal(ensembl_series, ensembl_result, check_names=False)

            ensembl_result = translator.translate_identifiers(ensembl_series, input_type=None, output_type='symbol')
            assert_series_equal(symbol_series, ensembl_result, check_names=False)

    def test_translate_to_same(self):
        assert_series_equal(self.human_translator.translate_identifiers(human_ensembl, input_type=None,
                                                                        output_type="ensembl"),
                            human_ensembl, check_names=False)

    def test_translation_fails(self):
        with pytest.raises(ValueError):
            self.human_translator.translate_identifiers([], input_type='nothing', output_type="something")
        with pytest.raises(ValueError):
            self.human_translator.translate_identifiers([], input_type=None, output_type="something")

    def test_guessing(self):
        assert 'ensembl' == atg.data.identifiers.guess_identifier_type(human_ensembl)
        assert 'symbol' == atg.data.identifiers.guess_identifier_type(human_symbol)
        assert 'entrez' == atg.data.identifiers.guess_identifier_type([1, 2, 3])
        assert 'entrez' == atg.data.identifiers.guess_identifier_type(['1', '2', '3'])

        assert 'symbol' == atg.data.identifiers.guess_identifier_type(['1', 1, 'a'])

    def test_mapping(self):
        human_dataframe = pandas.DataFrame({'ensembl': human_ensembl, 'symbol': human_symbol})
        map_result = self.human_translator.map_identifiers(human_ensembl, input_type=None, output_type='symbol')
        assert_frame_equal(human_dataframe, map_result)

        mouse_dataframe = pandas.DataFrame({'ensembl': mouse_ensembl, 'symbol': mouse_symbol})
        map_result = self.mouse_translator.map_identifiers(mouse_ensembl, input_type=None, output_type='symbol')
        assert_frame_equal(mouse_dataframe, map_result)

    def test_multiple_mapping(self):
        for symbol_series, ensembl_series, translator in ([human_symbol_multiple, human_ensembl_first,
                                                           self.human_translator],
                                                          [mouse_symbol_multiple, mouse_ensembl_first,
                                                           self.mouse_translator]):
            map_result = translator.map_identifiers(symbol_series, input_type=None, output_type='ensembl')
            assert len(map_result) > len(symbol_series)
            simple_translation_result = translator.translate_identifiers(symbol_series, None, 'ensembl')
            assert len(simple_translation_result) == len(symbol_series)
            assert simple_translation_result.tolist() == ensembl_series.tolist()
