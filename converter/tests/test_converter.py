from pathlib import Path
import tempfile
import unittest

from opendss_tools import convert


ROOT = Path(__file__).resolve().parent


class ConverterIntegrationTest(unittest.TestCase):
    def test_nontrivial_radial_reduction_matches_surviving_voltages(self):
        with tempfile.TemporaryDirectory() as directory:
            report = convert(
                ROOT / "fixture/full/Master.dss",
                ROOT / "fixture/assignment.csv",
                ROOT / "fixture/bus_order.csv",
                Path(directory),
            )
        self.assertEqual(report["full_bus_count"], 4)
        self.assertEqual(report["reduced_bus_count"], 3)
        self.assertLess(report["max_magnitude_error_pu"], 1e-10)


if __name__ == "__main__":
    unittest.main()
