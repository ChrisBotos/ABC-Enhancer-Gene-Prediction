import os
import sys
import unittest

import numpy as np
import pandas as pd

SCRIPTS_DIR = os.path.abspath("workflow/scripts")
sys.path.insert(0, SCRIPTS_DIR)
from predict import format_full_precision, get_predict_argument_parser


class TestFullPrecision(unittest.TestCase):
    def test_flag_is_off_by_default(self) -> None:
        required = ["--enhancers", "e", "--outdir", "o", "--chrom_sizes", "c"]
        required += ["--genes", "g", "--hic_pseudocount_distance", "5000"]
        args, _ = get_predict_argument_parser().parse_known_args(required)
        self.assertFalse(args.full_precision)

    def test_numerator_reads_back_exactly(self) -> None:
        values = pd.Series([1.234567891234e-07, 0.1 + 0.2, 12345.678901234, np.nan])
        written = format_full_precision(values)
        self.assertEqual(written.iloc[-1], "NaN")
        read_back = pd.to_numeric(written.replace("NaN", np.nan))
        pd.testing.assert_series_equal(read_back, values)
        self.assertNotEqual(f"{values.iloc[0]:.6f}", written.iloc[0])


if __name__ == "__main__":
    unittest.main()
