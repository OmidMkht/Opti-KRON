from pathlib import Path
import csv
import tempfile
import unittest

from opendss_tools import export_dataset


ROOT = Path(__file__).resolve().parent


class DatasetExporterIntegrationTest(unittest.TestCase):
    def test_ordered_per_unit_export_passes_serialized_kcl(self):
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            report=export_dataset(ROOT/'fixture/full/Master.dss',out,5.0)
            with (out/'bus.csv').open() as f: buses=list(csv.DictReader(f))
            with (out/'voltage.csv').open() as f: voltages=list(csv.DictReader(f))
            with (out/'load.csv').open() as f: injections=list(csv.DictReader(f))
            with (out/'bus.csv').open() as f: bus_header=next(csv.reader(f))
            with (out/'ybus.csv').open() as f: ybus_header=next(csv.reader(f))
            with (out/'load.csv').open() as f: load_header=next(csv.reader(f))
            equipment_files=['bus_coordinates.csv','regulator.csv','capacitor_bank.csv','phase_shift_equipment.csv']
            equipment_files_exist=all((out/name).exists() for name in equipment_files)
        self.assertEqual(report['bus_count'],4)
        self.assertEqual(report['phase_node_count'],12)
        self.assertEqual(len(voltages),sum(len(x['phases']) for x in buses))
        self.assertEqual(len(injections),len(voltages))
        self.assertEqual(bus_header,['bus_id','phases','base_kv','type'])
        self.assertEqual(ybus_header,['row','col','g','b'])
        self.assertEqual(load_header,['bus_id','phase','scenario','p_pu','q_pu'])
        self.assertLessEqual(report['serialized_max_abs_yv_minus_i_pu'],1e-12)
        self.assertLess(report['independent_max_abs_kcl_pu'],1e-6)
        self.assertTrue(equipment_files_exist)

    def test_parked_open_switch_is_not_stamped_as_an_edge(self):
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            export_dataset(ROOT/'fixture/open_switch/Master.dss',out,5.0)
            with (out/'bus.csv').open() as f: buses=list(csv.DictReader(f))
            with (out/'switch.csv').open() as f: switches=list(csv.DictReader(f))
        self.assertNotIn('b3_open',{row['bus_id'] for row in buses})
        self.assertEqual(
            [(row['switch_id'],row['to_bus'],row['status']) for row in switches],
            [('sw_open','b3','open')],
        )


if __name__=='__main__':
    unittest.main()
