from pathlib import Path
import tempfile
import unittest

from opendss_tools import convert_reduced_dataset


CONVERTER=Path(__file__).resolve().parents[1]
REPO=CONVERTER.parent


class ReducedDatasetConverterIntegrationTest(unittest.TestCase):
    """The full round trip, against the datasets the package actually ships."""

    def test_ieee123_57pct_aggregates_and_solves(self):
        with tempfile.TemporaryDirectory() as directory:
            report=convert_reduced_dataset(
                CONVERTER/'opendss_cases/123Bus/IEEE123Master.dss',
                REPO/'data/ieee123',
                REPO/'reduced_networks/ieee123/57pct',
                Path(directory),
                5.0,
            )
        self.assertEqual(report['reduced_bus_count'],56)
        self.assertEqual(report['radial_edge_count'],55)
        self.assertEqual(report['aggregation_max_abs_power_error_pu'],0.0)
        self.assertEqual(report['target_vs_full_max_complex_voltage_error_pu'],0.0)
        self.assertLess(report['relative_synthesized_ybus_error'],1e-9)
        self.assertLess(report['max_voltage_magnitude_error_pu'],0.002)


if __name__=='__main__':
    unittest.main()

