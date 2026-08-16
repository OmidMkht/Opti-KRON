from __future__ import annotations

import csv
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import opendssdirect as dss
from scipy import sparse

from .converter import (
    ConversionError,
    NODE_TO_PHASE,
    SKIP_CLASSES,
    Node,
    _complex_array,
    _element_nodes,
    _graph_edges,
    _source_command,
    _source_properties,
    _transformer_command,
    _validate_tree,
    decompose_radial_y,
    normalize_bus,
)


def _rows(path: Path) -> list[dict[str,str]]:
    with path.open(newline='',encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def _phase_nodes(bus_rows: list[dict[str,str]]) -> list[Node]:
    number={'a':1,'b':2,'c':3}
    return [Node(row['bus_id'].lower(),number[phase]) for row in bus_rows for phase in row['phases']]


def _read_complex(rows: list[dict[str,str]], real: str, imag: str) -> np.ndarray:
    return np.asarray([complex(float(row[real]),float(row[imag])) for row in rows])


def _read_ybus(path: Path, size: int) -> sparse.csr_matrix:
    rr=[];cc=[];vv=[]
    for row in _rows(path):
        rr.append(int(row['row'])-1);cc.append(int(row['col'])-1)
        vv.append(complex(float(row['g']),float(row['b'])))
    y=sparse.coo_matrix((vv,(rr,cc)),shape=(size,size),dtype=complex).tocsr()
    y.sum_duplicates();y.eliminate_zeros();y.sort_indices()
    return y


def _read_assignment(path: Path) -> tuple[dict[str,str],set[str]]:
    mapping: dict[str,str]={};kept:set[str]=set()
    for row in _rows(path):
        bus=normalize_bus(row['bus_id']);super_bus=normalize_bus(row['super_node'])
        if bus in mapping:
            raise ConversionError(f'Duplicate assignment for bus {bus!r}')
        mapping[bus]=super_bus
        if int(row['kept']):kept.add(bus)
    return mapping,kept


# The dataset exporter raises both limits to 100 before solving, and a feeder
# large enough to need that there needs it here too -- ieee8500 does not converge
# inside the OpenDSS defaults of 15 power-flow and 10 control iterations.
SOLVE_MAX_ITERATIONS = 100
SOLVE_MAX_CONTROL_ITERATIONS = 100


def _compile_and_solve(master: Path) -> None:
    dss.Basic.ClearAll()
    dss.Text.Command(f'compile "{master.resolve()}"')
    dss.Solution.MaxIterations(SOLVE_MAX_ITERATIONS)
    dss.Solution.MaxControlIterations(SOLVE_MAX_CONTROL_ITERATIONS)
    dss.Solution.Solve()


def _compile_full(master: Path, required_buses: set[str], reduced_nodes: list[Node]):
    _compile_and_solve(master)
    if not dss.Solution.Converged():
        raise ConversionError('The full OpenDSS circuit did not converge')
    source_name,source_bus,source_properties=_source_properties()
    bases={}
    for bus in sorted(required_buses):
        dss.Circuit.SetActiveBus(bus)
        if normalize_bus(dss.Bus.Name())!=bus:
            raise ConversionError(f'Assigned bus {bus!r} is absent from the compiled full circuit')
        value=1000.0*float(dss.Bus.kVBase())
        if not np.isfinite(value) or value<=0:
            raise ConversionError(f'Bus {bus!r} has invalid voltage base')
        bases[bus]=value
    node_index={node:i for i,node in enumerate(reduced_nodes)}
    preserved=np.zeros((len(reduced_nodes),len(reduced_nodes)),dtype=complex)
    transformer_commands=[]
    for element_name in dss.Circuit.AllElementNames():
        cls=element_name.split('.',1)[0].lower()
        if cls not in {'transformer','autotrans'}:continue
        dss.Circuit.SetActiveElement(element_name)
        if not dss.CktElement.Enabled():continue
        nterms=dss.CktElement.NumTerminals();nconds=dss.CktElement.NumConductors()
        local_nodes=_element_nodes(dss.CktElement.BusNames(),dss.CktElement.NodeOrder(),nconds)
        terminal_buses={node.bus for node in local_nodes if node is not None}
        missing=terminal_buses-set(required_buses)
        if missing:
            raise ConversionError(
                f'Transformer {element_name} has non-surviving terminal buses {sorted(missing)}'
            )
        yprim=_complex_array(dss.CktElement.YPrim());order=nterms*nconds
        if yprim.size!=order*order:
            raise ConversionError(f'Unexpected YPrim size for {element_name}')
        local_y=yprim.reshape((order,order),order='F')
        for a,node_a in enumerate(local_nodes):
            if node_a not in node_index:continue
            for b,node_b in enumerate(local_nodes):
                if node_b in node_index:preserved[node_index[node_a],node_index[node_b]]+=local_y[a,b]
        transformer_commands.append(_transformer_command(element_name))
    return (
        source_name,source_bus,source_properties,float(dss.Solution.Frequency()),bases,
        preserved,transformer_commands,
    )


def _validate_and_load(full_dataset: Path, reduced_dataset: Path, aggregation_tolerance: float):
    full_buses=_rows(full_dataset/'bus.csv');reduced_buses=_rows(reduced_dataset/'bus.csv')
    full_names={normalize_bus(row['bus_id']) for row in full_buses}
    mapping,kept=_read_assignment(reduced_dataset/'assignment.csv')
    if set(mapping)!=full_names:
        raise ConversionError(
            f'Assignment/full bus mismatch: missing={sorted(full_names-set(mapping))[:10]}, '
            f'extra={sorted(set(mapping)-full_names)[:10]}'
        )
    reduced_names=[normalize_bus(row['bus_id']) for row in reduced_buses]
    if set(reduced_names)!=kept or len(reduced_names)!=len(kept):
        raise ConversionError('Reduced bus list must equal the assignment rows marked kept')
    for bus,super_bus in mapping.items():
        if super_bus not in kept:
            raise ConversionError(f'Bus {bus!r} maps to non-surviving bus {super_bus!r}')
    phases={normalize_bus(row['bus_id']):set(row['phases']) for row in full_buses}
    for bus,super_bus in mapping.items():
        if not phases[bus].issubset(phases[super_bus]):
            raise ConversionError(f'Phase incompatibility: {bus} -> {super_bus}')

    nodes=_phase_nodes(reduced_buses)
    expected=[(node.bus,NODE_TO_PHASE[node.number]) for node in nodes]
    voltage_rows=_rows(reduced_dataset/'voltage.csv');load_rows=_rows(reduced_dataset/'load.csv')
    for name,rows in [('voltage.csv',voltage_rows),('load.csv',load_rows)]:
        actual=[(normalize_bus(row['bus_id']),row['phase'].lower()) for row in rows]
        if actual!=expected:
            raise ConversionError(f'{name} does not follow reduced bus/phase order')
    voltage=_read_complex(voltage_rows,'v_re_pu','v_im_pu')
    power=_read_complex(load_rows,'p_pu','q_pu')
    ypu=_read_ybus(reduced_dataset/'ybus.csv',len(nodes))

    aggregate=defaultdict(complex)
    for row in _rows(full_dataset/'load.csv'):
        aggregate[(mapping[normalize_bus(row['bus_id'])],row['phase'].lower())]+=complex(float(row['p_pu']),float(row['q_pu']))
    expected_power=np.asarray([aggregate[key] for key in expected])
    aggregation_error=float(np.max(np.abs(expected_power-power))) if power.size else 0.0
    if aggregation_error>aggregation_tolerance:
        raise ConversionError(
            f'Reduced load is not the phase-wise assignment aggregation; '
            f'maximum mismatch={aggregation_error:.6g} pu'
        )

    full_voltage={
        (normalize_bus(row['bus_id']),row['phase'].lower()):complex(float(row['v_re_pu']),float(row['v_im_pu']))
        for row in _rows(full_dataset/'voltage.csv')
    }
    target_error=float(np.max(np.abs(voltage-np.asarray([full_voltage[key] for key in expected]))))
    return mapping,kept,reduced_buses,nodes,ypu,voltage,power,aggregation_error,target_error


def _write_mapping(path: Path, mapping: dict[str,str], kept: set[str]) -> None:
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['original_bus','super_bus','survives'])
        for bus,super_bus in mapping.items():writer.writerow([bus,super_bus,int(bus in kept)])


def _write_equipment_mapping(path: Path, full_dataset: Path, mapping: dict[str,str]) -> int:
    specifications=[
        ('capacitor','capacitor_bank.csv','capacitor_id',[('bus','bus_id')]),
        ('transformer','transformer.csv','transformer_id',[('from','from_bus'),('to','to_bus')]),
        ('switch','switch.csv','switch_id',[('from','from_bus'),('to','to_bus')]),
        ('regulator','regulator.csv','regulator_id',[('from','from_bus'),('regulated','regulated_bus')]),
        ('phase_shift','phase_shift_equipment.csv','equipment_id',[('from','from_bus'),('to','to_bus')]),
    ]
    output=[]
    for equipment_class,filename,id_field,terminals in specifications:
        source=full_dataset/filename
        if not source.exists():continue
        for row in _rows(source):
            for terminal,bus_field in terminals:
                bus=normalize_bus(row[bus_field])
                if bus not in mapping:
                    raise ConversionError(f'{filename} references bus {bus!r} absent from assignment')
                output.append([equipment_class,row[id_field],terminal,bus,mapping[bus]])
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['equipment_class','equipment_id','terminal','original_bus','reduced_bus'])
        writer.writerows(output)
    return len(output)


def _write_loads(path: Path, nodes: list[Node], power_pu: np.ndarray, vbase: np.ndarray,
                 source_bus: str, sbase_va: float) -> None:
    lines=['! Phase-wise constant-PQ injections aggregated by assignment; positive DSS power is consumption.\n']
    for index,(node,s,base) in enumerate(zip(nodes,power_pu,vbase),1):
        if node.bus==source_bus or abs(s)<1e-15:continue
        lines.append(
            f'new Load.agg_{index} phases=1 bus1={node.bus}.{node.number} conn=wye '
            f'kv={base/1000.0:.16g} kw={-s.real*sbase_va/1000.0:.16g} '
            f'kvar={-s.imag*sbase_va/1000.0:.16g} model=1\n'
        )
    path.write_text(''.join(lines),encoding='ascii')


def _validate_output(master: Path, nodes: list[Node], target_pu: np.ndarray, target_ypu: sparse.csr_matrix,
                     vbase: np.ndarray, sbase_va: float, output_csv: Path) -> dict[str,float]:
    _compile_and_solve(master)
    if not dss.Solution.Converged():raise ConversionError('Generated reduced circuit did not converge')
    order=[Node(normalize_bus(x),int(x.rsplit('.',1)[1])) for x in dss.Circuit.YNodeOrder()]
    values=_complex_array(dss.Circuit.YNodeVArray());lookup=dict(zip(order,values))
    missing=[node for node in nodes if node not in lookup]
    if missing:raise ConversionError(f'Reduced circuit is missing nodes: {missing[:10]}')
    solved=np.asarray([lookup[node] for node in nodes])/vbase
    mag_error=np.abs(np.abs(solved)-np.abs(target_pu));complex_error=np.abs(solved-target_pu)
    node_index={node:i for i,node in enumerate(nodes)}
    passive=np.zeros((len(nodes),len(nodes)),dtype=complex)
    for element_name in dss.Circuit.AllElementNames():
        dss.Circuit.SetActiveElement(element_name);cls=element_name.split('.',1)[0].lower()
        if cls in SKIP_CLASSES or not dss.CktElement.Enabled():continue
        nterms=dss.CktElement.NumTerminals();nconds=dss.CktElement.NumConductors()
        local_nodes=_element_nodes(dss.CktElement.BusNames(),dss.CktElement.NodeOrder(),nconds)
        yprim=_complex_array(dss.CktElement.YPrim());order_local=nterms*nconds
        if yprim.size==0:continue
        if yprim.size!=order_local*order_local:
            raise ConversionError(f'Unexpected generated YPrim size for {element_name}')
        local_y=yprim.reshape((order_local,order_local),order='F')
        for a,node_a in enumerate(local_nodes):
            if node_a not in node_index:continue
            for b,node_b in enumerate(local_nodes):
                if node_b in node_index:passive[node_index[node_a],node_index[node_b]]+=local_y[a,b]
    actual_ypu=(vbase[:,None]/sbase_va)*passive*vbase[None,:]
    y_error=actual_ypu-target_ypu.toarray()
    with output_csv.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['bus_id','phase','target_v_mag_pu','solved_v_mag_pu','magnitude_error_pu','complex_error_pu'])
        for node,target,actual,me,ce in zip(nodes,target_pu,solved,mag_error,complex_error):
            writer.writerow([node.bus,NODE_TO_PHASE[node.number],f'{abs(target):.15f}',f'{abs(actual):.15f}',f'{me:.15f}',f'{ce:.15f}'])
    worst=int(np.argmax(mag_error))
    return {
        'max_voltage_magnitude_error_pu':float(mag_error[worst]),
        'max_voltage_magnitude_error_bus':nodes[worst].bus,
        'max_voltage_magnitude_error_phase':NODE_TO_PHASE[nodes[worst].number],
        'max_complex_voltage_error_pu':float(np.max(complex_error)),
        'max_abs_synthesized_ybus_error_pu':float(np.max(np.abs(y_error))),
        'relative_synthesized_ybus_error':float(
            np.max(np.abs(y_error))/max(np.max(np.abs(target_ypu.data)),np.finfo(float).tiny)
        ),
    }


def _convert_impl(full_master: Path, full_dataset: Path, reduced_dataset: Path, output: Path,
                  sbase_mva: float, matrix_tolerance: float, aggregation_tolerance: float):
    mapping,kept,bus_rows,nodes,ypu,target_pu,power_pu,aggregation_error,target_error=(
        _validate_and_load(full_dataset,reduced_dataset,aggregation_tolerance)
    )
    (
        source_name,source_bus,source_properties,frequency,bases,
        preserved_y_si,transformer_commands,
    )=_compile_full(full_master,set(kept),nodes)
    if source_bus not in kept:raise ConversionError(f'Source bus {source_bus!r} must survive')
    vbase=np.asarray([bases[node.bus] for node in nodes])
    sbase_va=sbase_mva*1e6
    ysi=sparse.diags(sbase_va/vbase)@ypu@sparse.diags(1.0/vbase)
    ysi=ysi.tocsr();ysi.eliminate_zeros();ysi.sort_indices()
    _validate_tree(ysi,nodes,matrix_tolerance)
    edges=_graph_edges(ysi.toarray(),nodes,matrix_tolerance)
    synthetic_y=ysi.toarray()-preserved_y_si
    branches,shunts,_=decompose_radial_y(synthetic_y,nodes,matrix_tolerance)

    output.mkdir(parents=True,exist_ok=True)
    source_data=type('SourceData',(),{'source_properties':source_properties})()
    (output/'Source.dss').write_text(_source_command(source_data),encoding='ascii')
    (output/'Transformers.dss').write_text(
        '! Preserved physical transformers and solved regulator taps\n'+''.join(transformer_commands),
        encoding='ascii',
    )
    (output/'Branches.dss').write_text('! Radial Kron-equivalent series elements\n'+''.join(branches),encoding='ascii')
    (output/'Shunts.dss').write_text('! Kron-equivalent grounded shunts\n'+''.join(shunts),encoding='ascii')
    _write_loads(output/'Loads.dss',nodes,power_pu,vbase,source_bus,sbase_va)
    source_kv=source_properties.get('basekv','1')
    master=(
        'clear\n'
        f'new Circuit.Reduced basekv={source_kv} phases=3 bus1={source_bus} frequency={frequency:g}\n'
        'redirect Source.dss\nredirect Transformers.dss\nredirect Branches.dss\nredirect Shunts.dss\nredirect Loads.dss\n'
        'set controlmode=off\nsolve mode=snapshot\n'
    )
    (output/'Master.dss').write_text(master,encoding='ascii')
    _write_mapping(output/'bus_mapping.csv',mapping,kept)
    equipment_mapping_rows=_write_equipment_mapping(output/'equipment_mapping.csv',full_dataset,mapping)
    metrics=_validate_output(
        output/'Master.dss',nodes,target_pu,ypu,vbase,sbase_va,output/'validation.csv'
    )
    report={
        'full_bus_count':len(mapping),'reduced_bus_count':len(kept),
        'reduced_phase_node_count':len(nodes),'radial_edge_count':len(edges),
        'preserved_transformer_count':len(transformer_commands),
        'equipment_mapping_row_count':equipment_mapping_rows,
        'sbase_mva':sbase_mva,'aggregation_max_abs_power_error_pu':aggregation_error,
        'target_vs_full_max_complex_voltage_error_pu':target_error,**metrics,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


def convert_reduced_dataset(full_master: Path, full_dataset: Path, reduced_dataset: Path,
                            output: Path, sbase_mva: float, matrix_tolerance: float=1e-9,
                            aggregation_tolerance: float=1e-12):
    original=Path.cwd()
    try:
        return _convert_impl(
            full_master.resolve(),full_dataset.resolve(),reduced_dataset.resolve(),output.resolve(),
            sbase_mva,matrix_tolerance,aggregation_tolerance,
        )
    finally:
        os.chdir(original)
