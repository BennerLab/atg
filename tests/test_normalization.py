import os
import numpy
import pandas
import unittest
import atg.quantification.normalization

RNASEQ_COUNT_FILE = os.path.join(os.path.dirname(__file__), 'data', 'rnaseq.csv')
# normalization factors calculated using edgeR's calcNormFactors
EDGER_NORM_FACTOR = [0.9632356, 0.9739227, 0.9621926, 0.9607130, 1.1352795, 1.0157446]
EDGER_TMM_VCP = [501.1333611, 540.4635384, 567.8398203, 345.3321018, 316.5981401, 377.8790519]
EDGER_CPM_VCP = [482.7128, 526.3641, 546.3750, 331.7619, 359.4104, 383.8492]
LIMMA_VOOM_VCP = [8.969148, 9.078111, 9.149407, 8.431966, 8.306551, 8.561956]

class CountNormalizationTest(unittest.TestCase):
    def setUp(self):
        self.read_count_df = pandas.read_csv(RNASEQ_COUNT_FILE, index_col=0)

    def test_tmm_normalization_factors(self):
        norm_factor = atg.quantification.normalization.tmm_norm_factor(self.read_count_df)
        numpy.testing.assert_almost_equal(EDGER_NORM_FACTOR, norm_factor, decimal=3)

    def test_tmm_normalization(self):
        normalized_count = atg.quantification.normalization.tmm_normalization(self.read_count_df)
        numpy.testing.assert_almost_equal(EDGER_TMM_VCP, normalized_count.ix['VCP', :], decimal=4)

    def test_cpm_normalization(self):
        normalized_count = atg.quantification.normalization.cpm_normalization(self.read_count_df)
        numpy.testing.assert_almost_equal(EDGER_CPM_VCP, normalized_count.ix['VCP', :], decimal=4)

    def test_voom_transformation(self):
        transformed_count = atg.quantification.normalization.voom(self.read_count_df)
        numpy.testing.assert_almost_equal(LIMMA_VOOM_VCP, transformed_count.ix['VCP', :], decimal=3)

if __name__ == '__main__':
    unittest.main()
