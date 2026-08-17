from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import opendssdirect as dss
from scipy import sparse


class DatasetError(RuntimeError):
    pass


INJECTION_CLASSES = {
    "load", "generator", "pvsystem", "storage", "isource", "vccs",
    "indmach012", "gicsource",
}
SOURCE_CLASSES = {"vsource"}
NON_ELECTRICAL_CLASSES = {
    "energymeter", "monitor", "sensor", "regcontrol", "capcontrol",
    "relay", "recloser", "fuse", "swtcontrol", "invcontrol", "expcontrol",
}
PHASE_NAMES = {1: "a", 2: "b", 3: "c"}
CSV_DECIMALS = 8
YBUS_DECIMALS = 15
STATE_DECIMALS = 15


@dataclass(frozen=True)
class PhaseNode:
    bus_id: str
    node_number: int


def _normalize_bus(spec: str) -> str:
    return spec.strip().strip('"').split(".")[0].lower()


def _logical_bus(spec: str) -> str:
    bus=_normalize_bus(spec)
    return bus[:-5] if bus.endswith('_open') else bus


def _complex(values) -> np.ndarray:
    raw = np.asarray(values, dtype=float)
    if raw.size == 1 and raw[0] == 0:
        return np.empty(0, dtype=complex)
    if raw.size % 2:
        raise DatasetError("OpenDSS returned an odd-length complex array")
    return raw[0::2] + 1j * raw[1::2]


def _format_float(value: float, decimals: int=CSV_DECIMALS) -> str:
    """Plain decimal output with bounded precision and no negative zero."""
    rounded=round(float(value),decimals)
    if rounded == 0:
        return '0'
    return f'{rounded:.{decimals}f}'.rstrip('0').rstrip('.')


def _formatted_rows(rows):
    for row in rows:
        yield [_format_float(value) if isinstance(value,(float,np.floating)) else value for value in row]


def _element_nodes(bus_specs: list[str], node_order: list[int], nconds: int) -> list[PhaseNode | None]:
    result: list[PhaseNode | None] = []
    for terminal, spec in enumerate(bus_specs):
        bus_id = _normalize_bus(spec)
        for conductor in range(nconds):
            local = terminal*nconds + conductor
            number = int(node_order[local]) if local < len(node_order) else conductor+1
            result.append(PhaseNode(bus_id, number) if number in (1,2,3) else None)
    return result


def _canonical_order(y_node_order: list[str]) -> tuple[list[str], list[PhaseNode]]:
    buses: list[str] = []
    phases_by_bus: dict[str,set[int]] = {}
    for entry in y_node_order:
        bus, raw_node = entry.rsplit(".",1)
        bus_id = _normalize_bus(bus)
        if bus_id.endswith('_open'):
            continue
        node = int(raw_node)
        if node not in (1,2,3):
            raise DatasetError(f"Unsupported non-phase node in Y order: {entry}")
        if bus_id not in phases_by_bus:
            buses.append(bus_id)
            phases_by_bus[bus_id] = set()
        phases_by_bus[bus_id].add(node)
    nodes = [PhaseNode(bus,node) for bus in buses for node in sorted(phases_by_bus[bus])]
    return buses,nodes


def _active_line_switch_state(name: str) -> tuple[bool,bool]:
    """Return (is_switch, is_open) and leave the corresponding line active."""
    dss.Lines.Name(name)
    dss.Circuit.SetActiveElement(f'Line.{name}')
    lname=name.lower()
    specs=dss.CktElement.BusNames()
    is_switch=bool(dss.Lines.IsSwitch() or lname.startswith('sw') or '_sw' in lname)
    conductor_open=any(
        dss.CktElement.IsOpen(term,cond)
        for term in range(1,dss.CktElement.NumTerminals()+1)
        for cond in range(1,dss.CktElement.NumConductors()+1)
    )
    parked_open=any(_normalize_bus(spec).endswith('_open') for spec in specs)
    return is_switch,bool(is_switch and (conductor_open or parked_open))


def _compile_and_solve(master: Path, max_iterations: int, max_control_iterations: int) -> None:
    dss.Basic.ClearAll()
    dss.Text.Command(f'compile "{master.resolve()}"')
    dss.Solution.MaxIterations(max_iterations)
    dss.Solution.MaxControlIterations(max_control_iterations)
    dss.Solution.Solve()
    if not dss.Solution.Converged():
        raise DatasetError(
            f"OpenDSS did not converge: iterations={dss.Solution.Iterations()}, "
            f"control_iterations={dss.Solution.ControlIterations()}"
        )


def _load_available_bus_coordinates(master: Path) -> None:
    """Load a colocated OpenDSS bus-coordinate file when the master omitted it."""
    names=list(dss.Circuit.AllBusNames())
    if any(dss.Circuit.SetActiveBus(name) and dss.Bus.Coorddefined() for name in names):
        return
    candidates=[
        path for path in master.parent.iterdir()
        if path.is_file()
        and any(token in path.stem.lower().replace("_","").replace("-","") for token in ("buscoord","busxy"))
        and path.suffix.lower() in {".csv",".dat",".txt",".dss"}
    ]
    if candidates:
        dss.Text.Command(f'BusCoords "{sorted(candidates,key=lambda p:p.name.lower())[0].resolve()}"')


def _assemble(master: Path, sbase_va: float, max_iterations: int, max_control_iterations: int):
    _compile_and_solve(master,max_iterations,max_control_iterations)
    _load_available_bus_coordinates(master)
    raw_order = list(dss.Circuit.YNodeOrder())
    buses,nodes = _canonical_order(raw_order)
    raw_nodes = [PhaseNode(_normalize_bus(x.rsplit('.',1)[0]),int(x.rsplit('.',1)[1])) for x in raw_order]
    raw_index = {node:i for i,node in enumerate(raw_nodes)}
    canonical_index = {node:i for i,node in enumerate(nodes)}
    permutation = np.asarray([raw_index[node] for node in nodes],dtype=int)
    voltage_raw = _complex(dss.Circuit.YNodeVArray())
    voltage = voltage_raw[permutation]

    n=len(nodes)
    rows: list[int]=[]; cols: list[int]=[]; values: list[complex]=[]
    element_injection = np.zeros(n,dtype=complex)
    ignored_with_yprim: list[str]=[]

    for element_name in dss.Circuit.AllElementNames():
        dss.Circuit.SetActiveElement(element_name)
        cls=element_name.split('.',1)[0].lower()
        if not dss.CktElement.Enabled():
            continue
        if cls == 'line':
            _,open_switch=_active_line_switch_state(element_name.split('.',1)[1])
            if open_switch:
                continue
        nterms=dss.CktElement.NumTerminals(); nconds=dss.CktElement.NumConductors()
        local_nodes=_element_nodes(dss.CktElement.BusNames(),dss.CktElement.NodeOrder(),nconds)
        if cls in INJECTION_CLASSES:
            currents=_complex(dss.CktElement.Currents())
            for k,node in enumerate(local_nodes):
                if node in canonical_index and k<currents.size:
                    element_injection[canonical_index[node]]-=currents[k]
            continue
        if cls in SOURCE_CLASSES or cls in NON_ELECTRICAL_CLASSES:
            continue
        yprim=_complex(dss.CktElement.YPrim())
        order=nterms*nconds
        if yprim.size==0:
            continue
        if yprim.size!=order*order:
            ignored_with_yprim.append(element_name)
            continue
        local_y=yprim.reshape((order,order),order='F')
        for a,node_a in enumerate(local_nodes):
            ia=canonical_index.get(node_a)
            if ia is None: continue
            for b,node_b in enumerate(local_nodes):
                ib=canonical_index.get(node_b)
                if ib is not None and local_y[a,b]!=0:
                    rows.append(ia); cols.append(ib); values.append(local_y[a,b])
    if ignored_with_yprim:
        raise DatasetError(f"Electrical primitives with unexpected dimensions: {ignored_with_yprim[:10]}")
    y_si=sparse.coo_matrix((values,(rows,cols)),shape=(n,n),dtype=complex).tocsr()
    y_si.sum_duplicates(); y_si.sort_indices()

    source_buses=set()
    for name in dss.Vsources.AllNames():
        dss.Circuit.SetActiveElement(f"Vsource.{name}")
        source_buses.add(_normalize_bus(dss.CktElement.BusNames()[0]))
    if not source_buses:
        raise DatasetError("Circuit has no voltage source")

    vbase=np.zeros(n)
    bus_metadata={}
    for bus_idx,bus in enumerate(buses,1):
        dss.Circuit.SetActiveBus(bus)
        kvbase=float(dss.Bus.kVBase())
        if not np.isfinite(kvbase) or kvbase<=0:
            raise DatasetError(f"Bus {bus!r} has invalid kVBase={kvbase}; establish voltage bases in DSS")
        coords=(float(dss.Bus.X()),float(dss.Bus.Y()))
        coord_defined=bool(dss.Bus.Coorddefined())
        bus_nodes=[node.node_number for node in nodes if node.bus_id==bus]
        bus_metadata[bus]={"bus_index":bus_idx,"phases":"".join(PHASE_NAMES[x] for x in bus_nodes),"base_kv_ln":kvbase,"x":coords[0],"y":coords[1],"coordinate_defined":coord_defined}
        for node in bus_nodes:
            vbase[canonical_index[PhaseNode(bus,node)]]=kvbase*1000.0

    # Phase-domain convention: common Sbase for every phase entry. This is
    # compatible with S_pu = V_pu*conj(I_pu) and sum(S_phase_pu)=S_total/Sbase.
    ibase=sbase_va/vbase
    y_pu=sparse.diags(1.0/ibase)@y_si@sparse.diags(vbase)
    y_pu=y_pu.tocsr(); y_pu.eliminate_zeros(); y_pu.sort_indices()
    v_pu=voltage/vbase
    i_pu=np.asarray(y_pu@v_pu).ravel()
    s_pu=v_pu*np.conj(i_pu)

    non_source=np.asarray([node.bus_id not in source_buses for node in nodes])
    independent_residual_si=y_si@voltage-element_injection
    independent_residual_pu=independent_residual_si/ibase
    return {
        "buses":buses,"nodes":nodes,"bus_metadata":bus_metadata,"source_buses":source_buses,
        "y_pu":y_pu,"v_pu":v_pu,"i_pu":i_pu,"s_pu":s_pu,
        "element_injection_pu":element_injection/ibase,"independent_residual_pu":independent_residual_pu,
        "non_source":non_source,"sbase_va":sbase_va,
    }


def _bus_phase_string(spec: str) -> str:
    parts=spec.strip().split('.')[1:]
    nums=[int(x) for x in parts if x.lstrip('-').isdigit() and int(x) in (1,2,3)]
    return ''.join(PHASE_NAMES[x] for x in nums)


def _bus_spec_nodes(spec: str) -> tuple[int,...]:
    return tuple(
        int(value) for value in spec.strip().split('.')[1:]
        if value.lstrip('-').isdigit()
    )


def _center_tap_pair(buses: list[str], winding_kvs: list[float], phase_count: int):
    """Return the two zero-based secondary winding indices of a center tap."""
    if phase_count!=1 or len(buses)<3:return None
    for left in range(1,len(buses)):
        for right in range(left+1,len(buses)):
            if _normalize_bus(buses[left])!=_normalize_bus(buses[right]):continue
            if abs(winding_kvs[left]-winding_kvs[right])>1e-8*max(
                1.0,abs(winding_kvs[left]),abs(winding_kvs[right])
            ):continue
            left_nodes=_bus_spec_nodes(buses[left]);right_nodes=_bus_spec_nodes(buses[right])
            if len(left_nodes)<2 or len(right_nodes)<2:continue
            forward=(
                left_nodes[-1]==right_nodes[0]
                and left_nodes[0]!=right_nodes[-1]
            )
            reverse=(
                left_nodes[0]==right_nodes[-1]
                and left_nodes[-1]!=right_nodes[0]
            )
            if forward or reverse:return left,right
    return None


def _write_equipment(output: Path, bus_metadata: dict, sbase_va: float):
    transformers=[]
    phase_shift_equipment=[]
    center_tapped_transformers=[]
    for name in dss.Transformers.AllNames():
        dss.Transformers.Name(name)
        dss.Circuit.SetActiveElement(f"Transformer.{name}")
        buses=dss.CktElement.BusNames()
        enabled=int(dss.CktElement.Enabled())
        phase_count=int(dss.Properties.Value('Phases'))
        lead_lag=dss.Properties.Value('LeadLag').strip().lower()
        winding_data=[]
        for winding,spec in enumerate(buses,1):
            dss.Transformers.Wdg(winding)
            bus=_normalize_bus(spec)
            winding_data.append({
                'winding':winding,'bus':bus,'bus_spec':spec.strip().lower(),
                'phases':_bus_phase_string(spec) or bus_metadata[bus]['phases'],
                'connection':'delta' if dss.Transformers.IsDelta() else 'wye',
                'kv':float(dss.Transformers.kV()),'kva':float(dss.Transformers.kVA()),
                'tap':float(dss.Transformers.Tap()),
            })
        primary=winding_data[0]
        from_bus=primary['bus'];from_phases=primary['phases']
        from_connection=primary['connection'];from_kv=primary['kv']
        for winding in range(2,len(buses)+1):
            detail=winding_data[winding-1]
            to_bus=detail['bus'];to_phases=detail['phases']
            to_connection=detail['connection'];to_kv=detail['kv']
            voltage_ratio=to_kv/from_kv if from_kv else ''
            transformers.append([name,winding-1,from_bus,from_phases,to_bus,to_phases,from_connection,to_connection,voltage_ratio,detail['tap'],enabled])
            if phase_count == 3 and from_connection != to_connection:
                # OpenDSS LeadLag specifies the LV displacement relative to HV:
                # Lag/ANSI=-30 degrees; Lead/Euro=+30 degrees. Report the
                # signed shifted-winding angle relative to winding 1.
                lv_relative_hv=30.0 if lead_lag in {'lead','euro'} else -30.0
                shift=lv_relative_hv if from_kv >= to_kv else -lv_relative_hv
                phase_shift_equipment.append([
                    name,from_bus,from_phases,to_bus,to_phases,
                    from_connection,to_connection,shift,enabled,
                ])
        center_pair=_center_tap_pair(buses,[item['kv'] for item in winding_data],phase_count)
        if center_pair is not None:
            for index,detail in enumerate(winding_data):
                if index==0:role='primary'
                elif index==center_pair[0]:role='secondary_1'
                elif index==center_pair[1]:role='secondary_2'
                else:role='other'
                center_tapped_transformers.append([
                    name,detail['winding'],role,detail['bus'],detail['bus_spec'],
                    detail['connection'],detail['kv']/from_kv if from_kv else '',
                    1000.0*detail['kva']/sbase_va,detail['tap'],enabled,
                ])
    with (output/'transformer.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['transformer_id','secondary_winding','from_bus','from_phases','to_bus','to_phases','from_connection','to_connection','voltage_ratio','tap_pu','enabled']); w.writerows(_formatted_rows(transformers))
    with (output/'phase_shift_equipment.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['equipment_id','from_bus','from_phases','to_bus','to_phases','from_connection','to_connection','phase_shift_deg','enabled'])
        w.writerows(_formatted_rows(phase_shift_equipment))
    with (output/'center_tapped_transformer.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow([
            'transformer_id','winding','role','bus_id','bus_spec','connection',
            'voltage_ratio','rated_s_pu','tap_pu','enabled',
        ])
        w.writerows(_formatted_rows(center_tapped_transformers))

    switches=[]
    for name in dss.Lines.AllNames():
        is_switch,opened=_active_line_switch_state(name)
        if not is_switch: continue
        specs=dss.CktElement.BusNames(); b1=_logical_bus(specs[0]); b2=_logical_bus(specs[1])
        switches.append([name,b1,_bus_phase_string(specs[0]) or bus_metadata[b1]['phases'],b2,_bus_phase_string(specs[1]) or bus_metadata[b2]['phases'],'open' if opened else 'closed',int(dss.CktElement.Enabled())])
    with (output/'switch.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['switch_id','from_bus','from_phases','to_bus','to_phases','status','enabled']); w.writerows(_formatted_rows(switches))

    regulators=[]
    for name in dss.RegControls.AllNames():
        dss.RegControls.Name(name)
        transformer=dss.RegControls.Transformer()
        winding=int(dss.RegControls.Winding())
        tap_winding=int(dss.RegControls.TapWinding())
        dss.Circuit.SetActiveElement(f'Transformer.{transformer}')
        specs=dss.CktElement.BusNames()
        from_bus=_normalize_bus(specs[0]) if specs else ''
        regulated_bus=_normalize_bus(specs[min(max(winding-1,0),len(specs)-1)]) if specs else ''
        if transformer in dss.Transformers.AllNames():
            dss.Transformers.Name(transformer); dss.Transformers.Wdg(tap_winding); tap=dss.Transformers.Tap()
        else: tap=''
        dss.Circuit.SetActiveElement(f'RegControl.{name}')
        vreg=dss.RegControls.ForwardVreg(); band=dss.RegControls.ForwardBand(); pt=dss.RegControls.PTRatio()
        base_v=1000.0*bus_metadata[regulated_bus]['base_kv_ln'] if regulated_bus in bus_metadata else 0.0
        vreg_pu=vreg*pt/base_v if base_v else ''
        band_pu=band*pt/base_v if base_v else ''
        regulators.append([name,transformer,from_bus,regulated_bus,winding,tap,vreg_pu,band_pu,int(dss.CktElement.Enabled())])
    with (output/'regulator.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['regulator_id','transformer_id','from_bus','regulated_bus','controlled_winding','tap_pu','vreg_pu','band_pu','enabled']); w.writerows(_formatted_rows(regulators))

    controls={}
    for name in dss.CapControls.AllNames():
        dss.CapControls.Name(name)
        capacitor=dss.CapControls.Capacitor().lower()
        controls.setdefault(capacitor,[]).append({
            'id':name,'mode':dss.CapControls.Mode(),'monitored_object':dss.CapControls.MonitoredObj(),
            'monitored_terminal':dss.CapControls.MonitoredTerm(),'on':dss.CapControls.ONSetting(),
            'off':dss.CapControls.OFFSetting(),'pt':dss.CapControls.PTRatio(),'ct':dss.CapControls.CTRatio(),
        })
    capacitors=[]
    for name in dss.Capacitors.AllNames():
        dss.Capacitors.Name(name); dss.Circuit.SetActiveElement(f'Capacitor.{name}')
        specs=dss.CktElement.BusNames(); bus=_normalize_bus(specs[0]); phases=_bus_phase_string(specs[0]) or bus_metadata[bus]['phases']
        associated=controls.get(name.lower(),[{}])
        for control in associated:
            kvar=dss.Capacitors.kvar()
            capacitors.append([name,bus,phases,'delta' if dss.Capacitors.IsDelta() else 'wye',1000.0*kvar/sbase_va,dss.Capacitors.NumSteps(),' '.join(str(x) for x in dss.Capacitors.States()),control.get('mode',''),int(dss.CktElement.Enabled())])
    with (output/'capacitor_bank.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['capacitor_id','bus_id','phases','connection','rated_q_pu','num_steps','states','control_mode','enabled']); w.writerows(_formatted_rows(capacitors))


def _write_dataset(output: Path, data: dict, scenario: str):
    output.mkdir(parents=True,exist_ok=True)
    buses=data['buses']; nodes=data['nodes']; meta=data['bus_metadata']; sources=data['source_buses']
    with (output/'bus.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['bus_id','phases','base_kv','type'])
        for bus in buses:
            m=meta[bus]; w.writerow([bus,m['phases'],'1.0','slack' if bus in sources else 'pq'])
    with (output/'bus_coordinates.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['bus_id','x','y'])
        for bus in buses:
            m=meta[bus]; w.writerow([bus,_format_float(m['x']) if m['coordinate_defined'] else '',_format_float(m['y']) if m['coordinate_defined'] else ''])
    original_y=data['y_pu'].tocoo()
    rounded_y_data=np.round(original_y.data.real,YBUS_DECIMALS)+1j*np.round(original_y.data.imag,YBUS_DECIMALS)
    nonzero=rounded_y_data != 0
    y=sparse.coo_matrix((rounded_y_data[nonzero],(original_y.row[nonzero],original_y.col[nonzero])),shape=original_y.shape)
    order=np.lexsort((y.col,y.row))
    with (output/'ybus.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['row','col','g','b'])
        for k in order: w.writerow([int(y.row[k])+1,int(y.col[k])+1,_format_float(y.data[k].real,YBUS_DECIMALS),_format_float(y.data[k].imag,YBUS_DECIMALS)])
    voltage=np.round(data['v_pu'].real,STATE_DECIMALS)+1j*np.round(data['v_pu'].imag,STATE_DECIMALS)
    current=np.asarray(y.tocsr()@voltage).ravel()
    power=voltage*np.conj(current)
    with (output/'voltage.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['bus_id','phase','scenario','v_re_pu','v_im_pu'])
        for node,v in zip(nodes,voltage):
            w.writerow([node.bus_id,PHASE_NAMES[node.node_number],scenario,_format_float(v.real,STATE_DECIMALS),_format_float(v.imag,STATE_DECIMALS)])
    with (output/'load.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['bus_id','phase','scenario','p_pu','q_pu'])
        for node,s in zip(nodes,power):
            w.writerow([node.bus_id,PHASE_NAMES[node.node_number],scenario,_format_float(s.real,STATE_DECIMALS),_format_float(s.imag,STATE_DECIMALS)])
    _write_equipment(output,meta,data['sbase_va'])


def _roundtrip_validate(output: Path, tolerance: float) -> dict:
    with (output/'bus.csv').open(encoding='utf-8') as f: buses=list(csv.DictReader(f))
    with (output/'voltage.csv').open(encoding='utf-8') as f: volts=list(csv.DictReader(f))
    with (output/'load.csv').open(encoding='utf-8') as f: injections=list(csv.DictReader(f))
    n=sum(len(x['phases']) for x in buses)
    if len(volts)!=n or len(injections)!=n:
        raise DatasetError('Serialized node-table dimensions do not match bus phase expansion')
    v=np.asarray([complex(float(x['v_re_pu']),float(x['v_im_pu'])) for x in volts])
    s=np.asarray([complex(float(x['p_pu']),float(x['q_pu'])) for x in injections])
    expected=np.conj(s/v)
    rr=[];cc=[];vv=[]
    with (output/'ybus.csv').open(encoding='utf-8') as f:
        for x in csv.DictReader(f): rr.append(int(x['row'])-1);cc.append(int(x['col'])-1);vv.append(complex(float(x['g']),float(x['b'])))
    y=sparse.coo_matrix((vv,(rr,cc)),shape=(n,n)).tocsr()
    residual=np.asarray(y@v).ravel()-expected
    maximum=float(np.max(np.abs(residual)))
    if maximum>tolerance:
        raise DatasetError(f'Serialized YV=I residual {maximum:.6g} pu exceeds {tolerance:.6g} pu')
    return {'serialized_max_abs_yv_minus_i_pu':maximum,'phase_node_count':n,'ybus_nnz':int(y.nnz)}


def export_dataset(master: Path, output: Path, sbase_mva: float, scenario: str='h001',
                   max_iterations: int=100, max_control_iterations: int=100,
                   serialized_tolerance: float=1e-12, independent_tolerance: float=1e-5) -> dict:
    original=Path.cwd(); master=master.resolve(); output=output.resolve()
    try:
        data=_assemble(master,sbase_mva*1e6,max_iterations,max_control_iterations)
        independent=float(np.max(np.abs(data['independent_residual_pu'][data['non_source']])))
        if independent>independent_tolerance:
            raise DatasetError(f'Independent element-current KCL residual {independent:.6g} pu exceeds {independent_tolerance:.6g} pu')
        _write_dataset(output,data,scenario)
        report=_roundtrip_validate(output,serialized_tolerance)
        magnitudes=np.abs(data['v_pu'])
        report.update({'bus_count':len(data['buses']),'sbase_mva':sbase_mva,'scenario':scenario,'independent_max_abs_kcl_pu':independent,'min_voltage_pu':float(np.min(magnitudes)),'max_voltage_pu':float(np.max(magnitudes)),'phase_nodes_below_0_5_pu':int(np.sum(magnitudes<0.5)),'converged':True})
        (output/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        return report
    finally:
        os.chdir(original)
