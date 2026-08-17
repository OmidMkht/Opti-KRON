import csv
from pathlib import Path
import tempfile
import unittest

import opendssdirect as dss

from opendss_tools import ConversionError, convert_reduced_dataset, export_dataset
from opendss_tools.reduced_dataset_converter import (
    _original_equipment_commands,
    _transformed_equipment_object,
)


CONVERTER=Path(__file__).resolve().parents[1]
REPO=CONVERTER.parent


class ReducedDatasetConverterIntegrationTest(unittest.TestCase):
    def test_location_only_serialization_covers_der_and_shunts(self):
        dss.Basic.ClearAll()
        for command in (
            'new Circuit.devices bus1=sourcebus basekv=12.47 phases=3',
            'new Line.feeder bus1=sourcebus bus2=b1 phases=3 r1=0.1 x1=0.2 length=1 units=km',
            'new Generator.g1 bus1=b1 phases=3 kv=7.2 kw=100 kvar=20 model=3',
            'new PVSystem.pv1 bus1=b1 phases=3 kv=12.47 kva=80 pmpp=70 model=1',
            'new Storage.bat1 bus1=b1 phases=3 kv=12.47 kwrated=50 kwhrated=200 %stored=60 state=idling',
            'new Capacitor.cap1 bus1=b1 phases=3 kv=12.47 kvar=300',
            'solve',
        ):
            dss.Text.Command(command)
        _,_,injections,shunts,counts,_,audit=_original_equipment_commands(
            {'sourcebus':'sourcebus','b1':'super'},
            {'sourcebus':7200.0,'b1':7200.0,'super':7200.0},
        )
        text=''.join(injections+shunts)
        self.assertEqual(counts,{'Generator':1,'PVSystem':1,'Storage':1,'Capacitor':1})
        self.assertIn('Generator.g1 Bus1="super"',text)
        self.assertIn('Model=3',text)
        self.assertIn('PVSystem.pv1',text)
        self.assertIn('Storage.bat1',text)
        self.assertIn('Capacitor.cap1 Bus1="super"',text)
        dss.Basic.ClearAll()
        dss.Text.Command('new Circuit.rebuilt bus1=sourcebus basekv=12.47 phases=3')
        dss.Text.Command(
            'new Line.feeder bus1=sourcebus bus2=super phases=3 r1=0.1 x1=0.2 length=1 units=km'
        )
        for command in injections+shunts:
            dss.Text.Command(command.strip())
        dss.Text.Command('solve')
        self.assertTrue(dss.Solution.Converged())
        self.assertEqual(set(dss.Generators.AllNames()),{'g1'})
        self.assertEqual(set(dss.PVsystems.AllNames()),{'pv1'})
        self.assertEqual(set(dss.Storages.AllNames()),{'bat1'})
        self.assertEqual(set(dss.Capacitors.AllNames()),{'cap1'})
        self.assertEqual(audit['transformed_equipment_count'],0)

    def test_regular_device_properties_are_rebased_in_per_unit(self):
        mapping={'low':'high','high':'high'}
        bases={'low':2400.0,'high':7200.0}
        load,ratio=_transformed_equipment_object(
            'Load',{'Name':'z','Bus1':'low.1','kV':2.4,'kW':10.0,'Model':2},
            mapping,bases,
        )
        source,_=_transformed_equipment_object(
            'Isource',{'Name':'i','Bus1':'low.1','Amps':30.0},mapping,bases
        )
        reactor,_=_transformed_equipment_object(
            'Reactor',{'Name':'r','Bus1':'low.1','Z':[3.0,4.0],'NormAmps':90.0},
            mapping,bases,
        )
        capacitor,_=_transformed_equipment_object(
            'Capacitor',{'Name':'c','Bus1':'low.1','kV':2.4,'Cuf':[9.0]},
            mapping,bases,
        )
        self.assertEqual(ratio,3.0)
        self.assertAlmostEqual(load['kV'],7.2)
        self.assertEqual(load['kW'],10.0)
        self.assertEqual(load['Model'],2)
        self.assertEqual(source['Amps'],10.0)
        self.assertEqual(reactor['Z'],[27.0,36.0])
        self.assertEqual(reactor['NormAmps'],30.0)
        self.assertAlmostEqual(capacitor['kV'],7.2)
        self.assertEqual(capacitor['Cuf'],[1.0])

    def test_default_voltage_rating_omitted_by_json_is_still_rebased(self):
        dss.Basic.ClearAll()
        dss.Text.Command('new Circuit.defaults bus1=source basekv=12.47 phases=3')
        dss.Text.Command('new Load.z bus1=low.1 phases=1 kw=10 model=2')
        _,loads,_,_,_,_,audit=_original_equipment_commands(
            {'source':'source','low':'high'},
            {'source':7200.0,'low':2400.0,'high':4800.0},
        )
        self.assertIn('kV=24.94',loads[0])
        self.assertEqual(audit['transformed_equipment_counts'],{'Load':1})

    def test_two_winding_wye_cross_voltage_conversion_solves(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);master=root/'Master.dss'
            master.write_text(
                'clear\n'
                'new Circuit.x bus1=source phases=3 basekv=12.47 pu=1\n'
                'new Transformer.t phases=3 windings=2 buses=[source low] '
                'conns=[wye wye] kvs=[12.47 4.16] kvas=[1000 1000] '
                '%rs=[0.5 0.5] xhl=2\n'
                'new Load.z bus1=low phases=3 conn=wye kv=4.16 kw=300 kvar=120 '
                'model=2 vminpu=0 vmaxpu=2\n'
                'set voltagebases=[12.47 4.16]\ncalcvoltagebases\nsolve\n',
                encoding='ascii',
            )
            full=root/'full';export_dataset(master,full,1.0)
            reduced=root/'reduced';reduced.mkdir()
            with (full/'bus.csv').open(newline='',encoding='utf-8') as stream:
                buses=list(csv.DictReader(stream))
            with (reduced/'assignment.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.writer(stream);writer.writerow(['bus_id','super_node','kept'])
                for row in buses:
                    writer.writerow([row['bus_id'],'source',int(row['bus_id'].lower()=='source')])
            source_bus=next(row for row in buses if row['bus_id'].lower()=='source')
            with (reduced/'bus.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=source_bus.keys())
                writer.writeheader();writer.writerow(source_bus)
            with (full/'load.csv').open(newline='',encoding='utf-8') as stream:
                loads=list(csv.DictReader(stream))
            with (reduced/'load.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=loads[0].keys());writer.writeheader()
                for phase in 'abc':
                    phase_rows=[row for row in loads if row['phase']==phase]
                    row=dict(phase_rows[0]);row['bus_id']='source'
                    row['p_pu']=sum(float(item['p_pu']) for item in phase_rows)
                    row['q_pu']=sum(float(item['q_pu']) for item in phase_rows)
                    writer.writerow(row)
            with (full/'voltage.csv').open(newline='',encoding='utf-8') as stream:
                voltages=[
                    row for row in csv.DictReader(stream) if row['bus_id'].lower()=='source'
                ]
            with (reduced/'voltage.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=voltages[0].keys())
                writer.writeheader();writer.writerows(voltages)
            output=root/'out'
            report=convert_reduced_dataset(master,full,reduced,output,1.0)
            load_text=(output/'Loads.dss').read_text(encoding='ascii')
            transformer_text=(output/'Transformers.dss').read_text(encoding='ascii')
        self.assertEqual(
            report['equipment_per_unit_transformation']['transformed_equipment_counts'],
            {'Load':1},
        )
        self.assertTrue(report['voltage_level_audit']['per_unit_transformed_equivalent_valid'])
        self.assertEqual(report['voltage_level_audit']['transformed_transformer_winding_count'],1)
        self.assertIn('kV=12.47',load_text)
        self.assertIn('Model=2',load_text)
        self.assertIn('kvs=[12.47 12.47]',transformer_text)
        self.assertLess(report['max_voltage_magnitude_error_pu'],1e-5)
        self.assertLess(report['relative_synthesized_ybus_error'],2e-8)

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

    def test_ieee34_preserves_mixed_wye_and_delta_loads(self):
        with tempfile.TemporaryDirectory() as directory:
            report=convert_reduced_dataset(
                CONVERTER/'opendss_cases/34Bus/ieee34Mod1.dss',
                REPO/'data/ieee34',
                REPO/'reduced_networks/ieee34/51pct',
                Path(directory),
                2.5,
            )
        self.assertEqual(report['load_representation'],'original-models')
        self.assertEqual(report['mapped_load_count'],68)
        self.assertEqual(report['relocated_equipment_counts'],{'Load':68,'Capacitor':2})
        self.assertEqual(report['original_load_model_counts'],{'1':38,'2':18,'4':2,'5':10})
        self.assertLess(report['relative_synthesized_ybus_error'],1e-12)
        self.assertLess(report['max_voltage_magnitude_error_pu'],0.0002)

    def test_ieee37_preserves_three_wire_delta_loads(self):
        with tempfile.TemporaryDirectory() as directory:
            report=convert_reduced_dataset(
                CONVERTER/'opendss_cases/37Bus/ieee37.dss',
                REPO/'data/ieee37',
                REPO/'reduced_networks/ieee37/67pct',
                Path(directory),
                2.5,
            )
        self.assertEqual(report['load_representation'],'original-models')
        self.assertEqual(report['mapped_load_count'],30)
        self.assertEqual(report['relocated_equipment_counts'],{'Load':30})
        self.assertEqual(report['original_load_model_counts'],{'1':15,'2':7,'4':8})
        self.assertLess(report['relative_synthesized_ybus_error'],1e-12)
        # This level reduces 39 -> 13 by crossing xfm1's 0.277/2.771 kV winding,
        # which the converter rebases. Budget is 0.003; the solved error is a
        # fifth of it.
        self.assertLess(report['max_voltage_magnitude_error_pu'],0.001)
        audit=report['voltage_level_audit']
        self.assertFalse(audit['location_only_voltage_level_valid'])
        self.assertTrue(audit['per_unit_transformed_equivalent_valid'])
        self.assertEqual(audit['incompatible_phase_domain_assignment_count'],0)
        self.assertEqual(audit['transformed_transformer_winding_count'],1)

    def test_ieee8500_pins_center_taps_so_no_domain_is_crossed(self):
        """The shipped ieee8500 levels reduce with `preserve = :required`.

        Center-tapped transformers are pinned, so the assignment never crosses a
        phase reference domain -- the one boundary no change of per-unit base can
        repair. Upstream's equivalent test uses a reduction that leaves them free
        and asserts the converter *rejects* it; here the property worth pinning is
        that there is nothing to reject.
        """
        with tempfile.TemporaryDirectory() as directory:
            report=convert_reduced_dataset(
                CONVERTER/'opendss_cases/8500-Node/Master.dss',
                REPO/'data/ieee8500',
                REPO/'reduced_networks/ieee8500/31pct',
                Path(directory),
                27.5,
            )
        self.assertEqual(report['reduced_bus_count'],3346)
        audit=report['voltage_level_audit']
        self.assertEqual(audit['incompatible_phase_domain_assignment_count'],0)
        self.assertTrue(audit['location_only_voltage_level_valid'])
        self.assertEqual(report['preserved_transformer_count'],1190)
        self.assertLess(report['max_energized_voltage_magnitude_error_pu'],0.005)


if __name__=='__main__':
    unittest.main()
