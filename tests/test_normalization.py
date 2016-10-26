import os
import numpy
import pandas
import unittest
from atg.quantification.normalization import tmm_norm_factor

RNASEQ_COUNT_FILE = os.path.join(os.path.dirname(__file__), 'data', 'rnaseq.csv')
# normalization factors calculated using edgeR's calcNormFactors
EDGER_NORM_FACTOR = [0.9632356, 0.9739227, 0.9621926, 0.9607130, 1.1352795, 1.0157446]


class CountNormalizationTest(unittest.TestCase):
    def setUp(self):
        self.read_count_df = pandas.read_csv(RNASEQ_COUNT_FILE, index_col=0)

    def test_tmm_normalization_factors(self):
        norm_factor = tmm_norm_factor(self.read_count_df)
        numpy.testing.assert_almost_equal(EDGER_NORM_FACTOR, norm_factor, decimal=3)


if __name__ == '__main__':
    unittest.main()
