from __future__ import annotations

import csv
import copy
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import opendssdirect as dss
from scipy import sparse
from scipy.sparse.linalg import splu

from .converter import (
    ConversionError,
    NODE_TO_PHASE,
    SKIP_CLASSES,
    Node,
    _complex_array,
    _element_nodes,
    _matrix_text,
    _source_command,
    _source_properties,
    _transformer_command,
    _write_reactor,
    normalize_bus,
)


MODEL_CLASSES=('LoadShape','GrowthShape','TShape','XYcurve','Spectrum','XfmrCode')
INJECTION_DEVICE_CLASSES=('Load','Generator','PVSystem','Storage','Isource','VCCS')
PHYSICAL_SHUNT_CLASSES=('Capacitor','Reactor')
TRANSFORM_REFERENCE_PROPERTIES={
    'Load':('kV',),'Generator':('kV',),'PVSystem':('kV',),'Storage':('kV',),
    'Isource':('Amps',),'VCCS':('VRated',),
    'Capacitor':('kV','NormAmps','EmergAmps'),
    'Reactor':('kV','NormAmps','EmergAmps'),
}


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


def _write_ybus(path: Path, y: sparse.spmatrix) -> None:
    matrix=y.tocoo();order=np.lexsort((matrix.col,matrix.row))
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['row','col','g','b'])
        for k in order:
            value=matrix.data[k]
            writer.writerow([
                int(matrix.row[k])+1,int(matrix.col[k])+1,
                f'{value.real:.15f}'.rstrip('0').rstrip('.') or '0',
                f'{value.imag:.15f}'.rstrip('0').rstrip('.') or '0',
            ])


def _rebuild_reduced_ybus(full_dataset: Path, reduced_dataset: Path, reduced_nodes: list[Node],
                          relative_tolerance: float, block_size: int=64):
    full_buses=_rows(full_dataset/'bus.csv');full_nodes=_phase_nodes(full_buses)
    full_y=_read_ybus(full_dataset/'ybus.csv',len(full_nodes)).tocsc()
    full_index={node:index for index,node in enumerate(full_nodes)}
    try:
        kept=np.asarray([full_index[node] for node in reduced_nodes],dtype=int)
    except KeyError as exc:
        raise ConversionError(f'Reduced phase node {exc.args[0]} is absent from the full Y-bus') from exc
    retained_mask=np.zeros(len(full_nodes),dtype=bool);retained_mask[kept]=True
    eliminated=np.flatnonzero(~retained_mask)
    ykk=full_y[kept,:][:,kept].tocsc()
    if eliminated.size==0:
        reduced=ykk.tocsr()
    else:
        ykr=full_y[kept,:][:,eliminated].tocsr()
        yrk=full_y[eliminated,:][:,kept].tocsc()
        yrr=full_y[eliminated,:][:,eliminated].tocsc()
        try:
            factor=splu(yrr)
        except RuntimeError as exc:
            raise ConversionError('Eliminated Y-bus block is singular; cannot rebuild Kron Y') from exc

        maximum=0.0
        for start in range(0,len(kept),block_size):
            stop=min(start+block_size,len(kept))
            block=ykk[:,start:stop].toarray()-ykr@factor.solve(yrk[:,start:stop].toarray())
            maximum=max(maximum,float(np.max(np.abs(block))))
        threshold=relative_tolerance*maximum
        rr=[];cc=[];vv=[]
        for start in range(0,len(kept),block_size):
            stop=min(start+block_size,len(kept))
            block=ykk[:,start:stop].toarray()-ykr@factor.solve(yrk[:,start:stop].toarray())
            retain=np.abs(block)>threshold
            diagonal=np.arange(start,stop)
            local_diagonal=diagonal-start
            retain[diagonal,local_diagonal]=True
            local_rows,local_cols=np.nonzero(retain)
            values=block[local_rows,local_cols]
            rr.extend(local_rows.tolist());cc.extend((local_cols+start).tolist());vv.extend(values.tolist())
        reduced=sparse.coo_matrix((vv,(rr,cc)),shape=(len(kept),len(kept)),dtype=complex).tocsr()
        reduced.sum_duplicates()
        # The phase-domain network is complex symmetric. Averaging reciprocal
        # entries removes asymmetric solver round-off introduced by block solves.
        reduced=((reduced+reduced.T)*0.5).tocsr()
        reduced.eliminate_zeros();reduced.sort_indices()
    try:
        splu(reduced.tocsc())
    except RuntimeError as exc:
        raise ConversionError(
            f'Sparsified Kron Y-bus lost rank at relative tolerance {relative_tolerance:g}'
        ) from exc
    output=reduced_dataset/'ybus.csv';_write_ybus(output,reduced)
    serialized=_read_ybus(output,len(reduced_nodes))
    try:
        splu(serialized.tocsc())
    except RuntimeError as exc:
        raise ConversionError('Serialized rebuilt Kron Y-bus is rank deficient') from exc
    metadata={
        'source':'rebuilt-kron','relative_sparsification_tolerance':relative_tolerance,
        'phase_node_count':len(reduced_nodes),'nnz':int(serialized.nnz),'full_rank':True,
    }
    (reduced_dataset/'ybus_rebuild.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    return serialized,metadata


def _sparse_graph_edges(y: sparse.spmatrix, nodes: list[Node], tolerance: float):
    bus_names=list(dict.fromkeys(node.bus for node in nodes))
    bus_position={bus:index for index,bus in enumerate(bus_names)}
    node_bus=np.asarray([bus_position[node.bus] for node in nodes],dtype=int)
    matrix=y.tocoo();edges=set()
    for row,col,value in zip(matrix.row,matrix.col,matrix.data):
        left=int(node_bus[row]);right=int(node_bus[col])
        if left!=right and abs(value)>tolerance:
            edges.add((min(left,right),max(left,right)))
    return sorted(edges)


def _validate_sparse_tree(y: sparse.spmatrix, nodes: list[Node], tolerance: float):
    bus_names=list(dict.fromkeys(node.bus for node in nodes));edges=_sparse_graph_edges(y,nodes,tolerance)
    if len(edges)!=max(0,len(bus_names)-1):
        raise ConversionError(
            f'Reduced admittance graph is not a tree: {len(bus_names)} buses, {len(edges)} edges'
        )
    adjacency=[[] for _ in bus_names]
    for left,right in edges:
        adjacency[left].append(right);adjacency[right].append(left)
    seen={0} if bus_names else set();stack=list(seen)
    while stack:
        current=stack.pop()
        for neighbor in adjacency[current]:
            if neighbor not in seen:seen.add(neighbor);stack.append(neighbor)
    if len(seen)!=len(bus_names):
        raise ConversionError('Reduced admittance graph is disconnected')
    return edges


def _write_ground_padded_reactor(name: str, left: list[Node], right: list[Node],
                                 y_matrix: np.ndarray, left_positions: list[int],
                                 right_positions: list[int], tolerance: float) -> str:
    size=y_matrix.shape[0];z=np.linalg.inv(y_matrix)
    reconstruction=np.linalg.norm(y_matrix@z@y_matrix-y_matrix,ord=np.inf)
    if reconstruction>max(tolerance,1e-7*np.linalg.norm(y_matrix,ord=np.inf)):
        raise ConversionError(f'Equivalent {name} cannot be represented by a padded Reactor Z matrix')
    def terminal(active: list[Node],positions: list[int]):
        numbers=['0']*size
        for node,position in zip(active,positions):numbers[position]=str(node.number)
        return active[0].bus+'.'+'.'.join(numbers)
    return (
        f'new Reactor.{name} phases={size} bus1={terminal(left,left_positions)} '
        f'bus2={terminal(right,right_positions)} '
        f'rmatrix=[{_matrix_text(z.real)}] xmatrix=[{_matrix_text(z.imag)}]\n'
    )


def _rectangular_series_admittance(sub: np.ndarray):
    rows,cols=sub.shape;size=rows+cols;result=np.zeros((size,size),complex)
    left_positions=list(range(rows));right_positions=list(range(rows,size))
    for row in range(rows):
        for col in range(cols):
            value=-sub[row,col]
            left=left_positions[row];right=right_positions[col]
            result[left,right]=value;result[right,left]=value
    scale=max(float(np.max(np.abs(result))),np.finfo(float).tiny)
    regularization=scale*1e-6*(1+0.01j)
    for _ in range(6):
        candidate=result+regularization*np.eye(size)
        if np.linalg.cond(candidate)<1e10:
            result=candidate;break
        regularization*=10
    else:
        raise ConversionError('Could not construct a stable grounded-conductor branch completion')
    return result,left_positions,right_positions


def _decompose_sparse_radial_y(y: sparse.spmatrix, nodes: list[Node], tolerance: float):
    bus_names=list(dict.fromkeys(node.bus for node in nodes))
    by_bus={bus:[i for i,node in enumerate(nodes) if node.bus==bus] for bus in bus_names}
    edges=_sparse_graph_edges(y,nodes,tolerance);branch_lines=[]
    incident={bus:np.zeros((len(by_bus[bus]),len(by_bus[bus])),complex) for bus in bus_names}
    for i,j in edges:
        bus_i=bus_names[i];bus_j=bus_names[j];rows=by_bus[bus_i];cols=by_bus[bus_j]
        block=y[rows,:][:,cols].toarray();reverse=y[cols,:][:,rows].toarray()
        block_scale=float(np.max(np.abs(block))) if block.size else 0.0
        coupling_tolerance=max(tolerance,1e-7*block_scale)
        active_rows=np.flatnonzero(np.max(np.abs(block),axis=1)>coupling_tolerance)
        active_cols=np.flatnonzero(np.max(np.abs(block),axis=0)>coupling_tolerance)
        sub=block[np.ix_(active_rows,active_cols)]
        reverse_sub=reverse[np.ix_(active_cols,active_rows)]
        if np.linalg.norm(sub-reverse_sub.T,ord=np.inf)>max(
            tolerance,1e-7*np.linalg.norm(sub,ord=np.inf)
        ):
            raise ConversionError(f'Edge {bus_i}-{bus_j} is non-reciprocal and cannot be written as a Reactor')
        left=[nodes[rows[k]] for k in active_rows];right=[nodes[cols[k]] for k in active_cols]
        if active_rows.size==active_cols.size:
            ys=-sub
            command=_write_reactor(f'kr_{i+1}_{j+1}',left,right,ys,tolerance)
            left_incident=ys;right_incident=ys.T
        else:
            ys,left_positions,right_positions=_rectangular_series_admittance(sub)
            command=_write_ground_padded_reactor(
                f'kr_{i+1}_{j+1}',left,right,ys,left_positions,right_positions,tolerance
            )
            left_incident=ys[np.ix_(left_positions,left_positions)]
            right_incident=ys[np.ix_(right_positions,right_positions)]
        if command:branch_lines.append(command)
        incident[bus_i][np.ix_(active_rows,active_rows)]+=left_incident
        incident[bus_j][np.ix_(active_cols,active_cols)]+=right_incident
    shunt_lines=[]
    for i,bus in enumerate(bus_names):
        indices=by_bus[bus]
        shunt=y[indices,:][:,indices].toarray()-incident[bus]
        command=_write_reactor(f'sh_{i+1}',[nodes[k] for k in indices],None,shunt,tolerance)
        if command:shunt_lines.append(command)
    return branch_lines,shunt_lines,edges


def _read_assignment(path: Path) -> tuple[dict[str,str],set[str]]:
    mapping: dict[str,str]={};kept:set[str]=set()
    for row in _rows(path):
        bus=normalize_bus(row['bus_id']);super_bus=normalize_bus(row['super_node'])
        if bus in mapping:
            raise ConversionError(f'Duplicate assignment for bus {bus!r}')
        mapping[bus]=super_bus
        if int(row['kept']):kept.add(bus)
    return mapping,kept


def _dss_literal(value) -> str:
    if isinstance(value,bool):
        return 'yes' if value else 'no'
    if isinstance(value,(int,float)):
        return f'{value:.16g}'
    if isinstance(value,list):
        items=[]
        for item in value:
            if isinstance(item,str) and not any(char.isspace() for char in item):
                items.append(item)
            else:
                items.append(_dss_literal(item))
        return '['+' '.join(items)+']'
    escaped=str(value).replace('"','""')
    return f'"{escaped}"'


def _mapped_bus_spec(spec: str, mapping: dict[str,str]) -> str:
    raw=str(spec).strip().strip('"')
    bus=normalize_bus(raw)
    if bus not in mapping:
        raise ConversionError(f'Equipment references bus {bus!r} absent from assignment')
    prefix=raw.split('.',1)[0]
    return mapping[bus]+raw[len(prefix):]


def _voltage_ratio(bus: str, mapping: dict[str,str], bases: dict[str,float]) -> float:
    """Return Vbase(new)/Vbase(old) for a mapped physical bus."""
    return bases[mapping[bus]]/bases[bus]


def _mapped_transformer_command(element_name: str, mapping: dict[str,str],
                                bases: dict[str,float]) -> str:
    command=_transformer_command(element_name)
    replacements=[]
    original_buses=dss.CktElement.BusNames()
    for original in original_buses:
        replacements.append((original,_mapped_bus_spec(original,mapping)))
    for original,mapped in sorted(replacements,key=lambda pair:len(pair[0]),reverse=True):
        command=command.replace(original,mapped)
    raw_kvs=str(dss.Properties.Value('kvs')).strip().strip('[]()')
    kvs=[float(value) for value in raw_kvs.replace(',',' ').split()]
    if len(kvs)!=len(original_buses):
        raise ConversionError(
            f'{element_name} exposes {len(original_buses)} terminals but {len(kvs)} winding kVs'
        )
    transformed=[]
    for kv,spec in zip(kvs,original_buses):
        bus=normalize_bus(spec)
        transformed.append(kv*_voltage_ratio(bus,mapping,bases))
    command=command.rstrip()+f' kvs=[{" ".join(f"{value:.16g}" for value in transformed)}]\n'
    return command


def _bus_spec_has_phase_conductor(spec: str) -> bool:
    parts=str(spec).strip().strip('"').split('.')[1:]
    return any(part.lstrip('-').isdigit() and int(part) in (1,2,3) for part in parts)


def _json_command(cls: str, obj: dict, mapping: dict[str,str] | None=None) -> str:
    name=obj.get('Name')
    if not name:
        raise ConversionError(f'{cls} JSON object has no name')
    properties=[]
    for key,value in obj.items():
        if key=='Name':
            continue
        if mapping is not None and key.lower() in {'bus1','bus2'}:
            value=_mapped_bus_spec(value,mapping)
        elif mapping is not None and key.lower() in {'bus','buses'}:
            value=[_mapped_bus_spec(spec,mapping) for spec in value]
        dss_key='%'+key[3:] if key.startswith('pct') and len(key)>3 and key[3].isupper() else key
        properties.append(f'{dss_key}={_dss_literal(value)}')
    return f'new {cls}.{name} '+' '.join(properties)+'\n'


def _scaled_value(value, scale: float):
    if isinstance(value,bool):
        return value
    if isinstance(value,(int,float)):
        return value*scale
    if isinstance(value,list):
        return [_scaled_value(item,scale) for item in value]
    raise ConversionError(f'Cannot scale nonnumeric OpenDSS property value {value!r}')


def _object_bus_ratios(obj: dict, mapping: dict[str,str], bases: dict[str,float]) -> list[float]:
    ratios=[]
    for key,value in obj.items():
        if key.lower() in {'bus1','bus2'}:
            specifications=[value]
        elif key.lower() in {'bus','buses'}:
            specifications=value if isinstance(value,list) else [value]
        else:
            continue
        for spec in specifications:
            if not str(spec).strip():continue
            bus=normalize_bus(str(spec))
            if bus not in mapping:
                raise ConversionError(f'Equipment references bus {bus!r} absent from assignment')
            ratios.append(_voltage_ratio(bus,mapping,bases))
    return ratios


def _transformed_equipment_object(cls: str, obj: dict, mapping: dict[str,str],
                                  bases: dict[str,float]) -> tuple[dict,float]:
    """Rebase a regular device while preserving its voltage-per-unit behavior."""
    result=copy.deepcopy(obj)
    ratios=_object_bus_ratios(obj,mapping,bases)
    ratio=ratios[0] if ratios else 1.0
    if any(abs(value-ratio)>1e-9*max(1.0,abs(value),abs(ratio)) for value in ratios[1:]):
        raise ConversionError(
            f'{cls}.{obj.get("Name", "?")} has terminals with different voltage-base mappings; '
            'a regular scalar transformed equivalent is not possible'
        )
    voltage_keys={
        'Load':{'kv'},'Generator':{'kv'},'PVSystem':{'kv'},'Storage':{'kv'},
        'VCCS':{'vrated'},'Capacitor':{'kv'},'Reactor':{'kv'},
    }.get(cls,set())
    current_keys={
        'Isource':{'amps'},'Capacitor':{'normamps','emergamps'},
        'Reactor':{'normamps','emergamps'},
    }.get(cls,set())
    impedance_keys={
        'Capacitor':{'r','xl'},
        'Reactor':{'rmatrix','xmatrix','r','x','rp','z1','z2','z0','z','lmh'},
    }.get(cls,set())
    admittance_keys={'Capacitor':{'cmatrix','cuf'}}.get(cls,set())
    for key,value in list(result.items()):
        normalized=key.lower()
        if normalized in voltage_keys:
            result[key]=_scaled_value(value,ratio)
        elif normalized in current_keys:
            result[key]=_scaled_value(value,1.0/ratio)
        elif normalized in impedance_keys:
            result[key]=_scaled_value(value,ratio*ratio)
        elif normalized in admittance_keys:
            result[key]=_scaled_value(value,1.0/(ratio*ratio))
    return result,ratio


def _complete_transform_reference_properties(cls: str, obj: dict) -> dict:
    """Add transform-critical values that compact OpenDSS JSON may omit."""
    required=TRANSFORM_REFERENCE_PROPERTIES.get(cls,())
    missing=[key for key in required if key.lower() not in {item.lower() for item in obj}]
    if not missing:return obj
    name=obj.get('Name')
    if not name or not dss.Circuit.SetActiveElement(f'{cls}.{name}'):
        raise ConversionError(f'Cannot reactivate {cls}.{name} to read transformation properties')
    result=copy.deepcopy(obj)
    for key in missing:
        value=dss.Properties.Value(key)
        if value not in (None,''):
            try:result[key]=float(value)
            except ValueError as exc:
                raise ConversionError(
                    f'{cls}.{name} property {key}={value!r} is not a numeric scalar'
                ) from exc
    return result


def _original_equipment_commands(mapping: dict[str,str], bases: dict[str,float]):
    """Serialize original models at mapped buses, rebased to retain per-unit behavior."""
    circuit=json.loads(dss.Circuit.ToJSON())
    models=[];loads=[];injections=[];shunts=[];counts={};load_models=defaultdict(int)
    transformed_counts=defaultdict(int);transformed_examples=[]
    for cls in MODEL_CLASSES:
        for obj in circuit.get(cls,[]):
            models.append(_json_command(cls,obj))
    for cls in INJECTION_DEVICE_CLASSES:
        objects=circuit.get(cls,[])
        if objects:counts[cls]=len(objects)
        for obj in objects:
            obj=_complete_transform_reference_properties(cls,obj)
            transformed,ratio=_transformed_equipment_object(cls,obj,mapping,bases)
            command=_json_command(cls,transformed,mapping)
            if abs(ratio-1.0)>1e-9:
                transformed_counts[cls]+=1
                if len(transformed_examples)<12:
                    transformed_examples.append({
                        'class':cls,'element':str(obj.get('Name','')),'voltage_base_ratio':ratio,
                    })
            if cls=='Load':
                loads.append(command)
                load_models[str(obj.get('Model',1))]+=1
            else:
                injections.append(command)
    for cls in PHYSICAL_SHUNT_CLASSES:
        objects=circuit.get(cls,[])
        copied=0
        for obj in objects:
            # A two-terminal reactor is a series network element, not a shunt.
            if cls=='Reactor' and obj.get('Bus2'):
                bus1=normalize_bus(str(obj.get('Bus1','')))
                bus2=normalize_bus(str(obj['Bus2']))
                if bus1!=bus2 or _bus_spec_has_phase_conductor(obj['Bus2']):continue
            obj=_complete_transform_reference_properties(cls,obj)
            transformed,ratio=_transformed_equipment_object(cls,obj,mapping,bases)
            shunts.append(_json_command(cls,transformed,mapping));copied+=1
            if abs(ratio-1.0)>1e-9:
                transformed_counts[cls]+=1
                if len(transformed_examples)<12:
                    transformed_examples.append({
                        'class':cls,'element':str(obj.get('Name','')),'voltage_base_ratio':ratio,
                    })
        if copied:counts[cls]=copied
    transformation_audit={
        'transformed_equipment_counts':dict(sorted(transformed_counts.items())),
        'transformed_equipment_count':sum(transformed_counts.values()),
        'examples':transformed_examples,
    }
    return models,loads,injections,shunts,counts,dict(load_models),transformation_audit


def _mapped_constant_pq_loads(mapping: dict[str,str]) -> tuple[list[str],complex]:
    """Freeze solved OpenDSS loads as constant-PQ devices at their super-nodes.

    Keeping the original terminal specification is important for three-wire
    feeders: replacing a phase-to-phase delta load with two phase-to-ground
    nodal loads changes the common-mode voltage even at identical total power.
    """
    commands=[]
    total_kva=0j
    for index,name in enumerate(dss.Loads.AllNames(),1):
        dss.Loads.Name(name)
        dss.Circuit.SetActiveElement(f'Load.{name}')
        if not dss.CktElement.Enabled():
            continue
        specs=dss.CktElement.BusNames()
        if len(specs)!=1:
            raise ConversionError(f'Load {name!r} has {len(specs)} terminals; expected one')
        original_bus=normalize_bus(specs[0])
        if original_bus not in mapping:
            raise ConversionError(f'Load {name!r} references bus {original_bus!r} absent from assignment')
        raw_spec=specs[0].strip().strip('"')
        suffix=raw_spec[len(raw_spec.split('.',1)[0]):]
        mapped_spec=f'{mapping[original_bus]}{suffix}'
        terminal_power=_complex_array(dss.CktElement.TotalPowers())
        if terminal_power.size!=1:
            raise ConversionError(f'Load {name!r} returned unexpected terminal powers')
        power_kva=terminal_power[0]
        total_kva+=power_kva
        commands.append(
            f'new Load.mapped_{index} phases={dss.CktElement.NumPhases()} '
            f'bus1={mapped_spec} conn={"delta" if dss.Loads.IsDelta() else "wye"} '
            f'kv={dss.Loads.kV():.16g} kw={power_kva.real:.16g} '
            f'kvar={power_kva.imag:.16g} model=1 vminpu=0 vmaxpu=2\n'
        )
    return commands,total_kva


def _array_tokens(value: str) -> list[str]:
    return str(value).strip().strip('[]()').replace(',',' ').lower().split()


def _phase_domain_violations(mapping: dict[str,str], phase_shift_names: set[str]):
    """Find assignments crossing transformer connections that require a dense T.

    Regular two-winding transformers with identical winding connections join a
    phase-reference domain, even when their voltage bases differ. Multi-winding
    and connection-changing transformers form domain boundaries.
    """
    parent={bus:bus for bus in mapping}

    def find(bus: str) -> str:
        while parent[bus]!=bus:
            parent[bus]=parent[parent[bus]];bus=parent[bus]
        return bus

    def union(left: str, right: str) -> None:
        left=find(left);right=find(right)
        if left!=right:parent[right]=left

    incompatible={}
    for element_name in dss.Circuit.AllElementNames():
        cls,name=element_name.split('.',1);cls=cls.lower();name=name.lower()
        if cls not in {'transformer','autotrans'}:continue
        dss.Circuit.SetActiveElement(element_name)
        if not dss.CktElement.Enabled():continue
        windings=int(float(dss.Properties.Value('windings') or 0))
        conns=_array_tokens(dss.Properties.Value('conns'))
        reasons=[]
        if windings!=2:reasons.append(f'{windings}-winding')
        if len(conns)>=2 and len(set(conns))>1:reasons.append('connection-changing')
        if name in phase_shift_names:reasons.append('phase-shifting')
        if reasons:incompatible[element_name.lower()]=','.join(reasons)

    branch_classes={'line','reactor','transformer','autotrans'}
    for element_name in dss.Circuit.AllElementNames():
        cls=element_name.split('.',1)[0].lower()
        if cls not in branch_classes or element_name.lower() in incompatible:continue
        dss.Circuit.SetActiveElement(element_name)
        if not dss.CktElement.Enabled():continue
        buses=list(dict.fromkeys(
            normalize_bus(spec) for spec in dss.CktElement.BusNames()
            if normalize_bus(spec) in mapping
        ))
        for bus in buses[1:]:union(buses[0],bus)

    violations=[]
    for bus,super_bus in mapping.items():
        if find(bus)!=find(super_bus):violations.append((bus,super_bus))
    return violations,incompatible


def _compile_full(master: Path, required_buses: set[str], reduced_nodes: list[Node],
                  mapping: dict[str,str], load_representation: str,
                  phase_shift_names: set[str]):
    dss.Basic.ClearAll();dss.Text.Command(f'compile "{master.resolve()}"')
    dss.Solution.MaxIterations(100);dss.Solution.MaxControlIterations(100);dss.Solution.Solve()
    if not dss.Solution.Converged():
        raise ConversionError('The full OpenDSS circuit did not converge')
    source_name,source_bus,source_properties=_source_properties()
    all_bases={}
    for bus in sorted(mapping):
        dss.Circuit.SetActiveBus(bus)
        if normalize_bus(dss.Bus.Name())!=bus:
            raise ConversionError(f'Assigned bus {bus!r} is absent from the compiled full circuit')
        value=1000.0*float(dss.Bus.kVBase())
        if not np.isfinite(value) or value<=0:
            raise ConversionError(f'Bus {bus!r} has invalid voltage base')
        all_bases[bus]=value
    bases={bus:all_bases[bus] for bus in required_buses}
    domain_violations,incompatible_transformers=_phase_domain_violations(
        mapping,phase_shift_names
    )
    if domain_violations:
        examples=', '.join(f'{bus}->{super_bus}' for bus,super_bus in domain_violations[:8])
        barriers=', '.join(
            f'{name} ({reason})' for name,reason in list(incompatible_transformers.items())[:8]
        )
        raise ConversionError(
            f'Assignment crosses incompatible transformer phase domains at '
            f'{len(domain_violations)} buses: {examples}. Domain boundaries: {barriers}'
        )
    assignment_base_mismatches=[
        (bus,super_bus,all_bases[bus],all_bases[super_bus])
        for bus,super_bus in mapping.items()
        if abs(all_bases[bus]-all_bases[super_bus])>1e-6*max(all_bases[bus],all_bases[super_bus])
    ]
    equipment_base_mismatches=defaultdict(list)
    audited_classes={
        'load','generator','pvsystem','storage','isource','vccs','capacitor',
        'transformer','autotrans','reactor',
    }
    for element_name in dss.Circuit.AllElementNames():
        cls=element_name.split('.',1)[0].lower()
        if cls not in audited_classes:continue
        dss.Circuit.SetActiveElement(element_name)
        specs=dss.CktElement.BusNames();nterms=dss.CktElement.NumTerminals()
        if cls=='reactor' and nterms>1 and len({normalize_bus(spec) for spec in specs})>1:
            continue
        for spec in specs:
            bus=normalize_bus(spec);super_bus=mapping[bus]
            if abs(all_bases[bus]-all_bases[super_bus])>1e-6*max(
                all_bases[bus],all_bases[super_bus]
            ):
                record=(element_name,bus,super_bus,all_bases[bus],all_bases[super_bus])
                if record not in equipment_base_mismatches[cls]:
                    equipment_base_mismatches[cls].append(record)
    voltage_level_audit={
        'assignment_bus_mismatch_count':len(assignment_base_mismatches),
        'equipment_terminal_mismatch_counts':{
            cls:len(records) for cls,records in sorted(equipment_base_mismatches.items())
        },
        'location_only_voltage_level_valid':not any(equipment_base_mismatches.values()),
        'per_unit_transformed_equivalent_valid':True,
        'incompatible_phase_domain_assignment_count':0,
        'examples':[
            {
                'class':cls,'element':record[0],'original_bus':record[1],
                'super_bus':record[2],'original_kv_ln':record[3]/1000.0,
                'super_kv_ln':record[4]/1000.0,
            }
            for cls,records in sorted(equipment_base_mismatches.items()) for record in records[:3]
        ],
    }
    node_index={node:i for i,node in enumerate(reduced_nodes)}
    preserved_rows=[];preserved_cols=[];preserved_values=[];transformer_commands=[]
    transformed_transformer_windings=0
    physical_shunt_count=0
    for element_name in dss.Circuit.AllElementNames():
        cls=element_name.split('.',1)[0].lower()
        dss.Circuit.SetActiveElement(element_name)
        if not dss.CktElement.Enabled():continue
        nterms=dss.CktElement.NumTerminals();nconds=dss.CktElement.NumConductors()
        local_nodes=_element_nodes(dss.CktElement.BusNames(),dss.CktElement.NodeOrder(),nconds)
        is_transformer=cls in {'transformer','autotrans'}
        is_physical_shunt=(
            load_representation=='original-models'
            and (
                cls=='capacitor'
                or (
                    cls=='reactor'
                    and (nterms==1 or all(node is None for node in local_nodes[nconds:]))
                )
            )
        )
        if not is_transformer and not is_physical_shunt:continue
        terminal_buses={node.bus for node in local_nodes if node is not None}
        missing=terminal_buses-set(mapping)
        if missing:
            raise ConversionError(f'{element_name} references buses absent from assignment: {sorted(missing)}')
        if is_transformer or is_physical_shunt:
            local_nodes=[Node(mapping[node.bus],node.number) if node is not None else None for node in local_nodes]
        yprim=_complex_array(dss.CktElement.YPrim());order=nterms*nconds
        if yprim.size!=order*order:
            raise ConversionError(f'Unexpected YPrim size for {element_name}')
        local_y=yprim.reshape((order,order),order='F')
        # V_old = D V_new and I_new = D^H I_old. Voltage-base ratios are real,
        # hence Y_new = D Y_old D for a scalar/diagonal per-unit rebasing.
        local_scale=np.ones(order)
        for index,node in enumerate(_element_nodes(
            dss.CktElement.BusNames(),dss.CktElement.NodeOrder(),nconds
        )):
            if node is not None:
                local_scale[index]=all_bases[node.bus]/all_bases[mapping[node.bus]]
        local_y=local_scale[:,None]*local_y*local_scale[None,:]
        for a,node_a in enumerate(local_nodes):
            if node_a not in node_index:continue
            for b,node_b in enumerate(local_nodes):
                if node_b in node_index and local_y[a,b]!=0:
                    preserved_rows.append(node_index[node_a]);preserved_cols.append(node_index[node_b])
                    preserved_values.append(local_y[a,b])
        if is_transformer:
            transformed_transformer_windings+=sum(
                abs(_voltage_ratio(normalize_bus(spec),mapping,all_bases)-1.0)>1e-9
                for spec in dss.CktElement.BusNames()
            )
            transformer_commands.append(_mapped_transformer_command(element_name,mapping,all_bases))
        if is_physical_shunt:
            physical_shunt_count+=1
    original_equipment=_original_equipment_commands(mapping,all_bases)
    voltage_level_audit['transformed_transformer_winding_count']=(
        transformed_transformer_windings
    )
    preserved=sparse.coo_matrix(
        (preserved_values,(preserved_rows,preserved_cols)),
        shape=(len(reduced_nodes),len(reduced_nodes)),dtype=complex,
    ).tocsr()
    preserved.sum_duplicates();preserved.eliminate_zeros();preserved.sort_indices()
    load_commands,total_load_kva=_mapped_constant_pq_loads(mapping)
    unsupported_injections=[]
    for element_name in dss.Circuit.AllElementNames():
        cls=element_name.split('.',1)[0].lower()
        if cls in {'generator','pvsystem','storage','isource','vccs','indmach012','gicsource'}:
            dss.Circuit.SetActiveElement(element_name)
            if dss.CktElement.Enabled():unsupported_injections.append(element_name)
    if load_representation=='mapped-elements' and unsupported_injections:
        raise ConversionError(
            'Mapped-equipment load realization does not yet support non-load injection elements: '
            f'{unsupported_injections[:10]}'
        )
    return (
        source_name,source_bus,source_properties,float(dss.Solution.Frequency()),bases,
        preserved,transformer_commands,load_commands,total_load_kva,original_equipment,
        physical_shunt_count,float(dss.Solution.LoadMult()),voltage_level_audit,
    )


def _validate_and_load(full_dataset: Path, reduced_dataset: Path, aggregation_tolerance: float,
                       kron_relative_tolerance: float):
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
    for bus in kept:
        if mapping[bus]!=bus:
            raise ConversionError(f'Surviving bus {bus!r} must represent itself')
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
    ybus_path=reduced_dataset/'ybus.csv'
    if ybus_path.exists():
        ypu=_read_ybus(ybus_path,len(nodes))
        metadata_path=reduced_dataset/'ybus_rebuild.json'
        ybus_info=(json.loads(metadata_path.read_text(encoding='utf-8')) if metadata_path.exists()
                   else {'source':'csv','relative_sparsification_tolerance':None,
                         'phase_node_count':len(nodes),'nnz':int(ypu.nnz)})
    else:
        ypu,ybus_info=_rebuild_reduced_ybus(
            full_dataset,reduced_dataset,nodes,kron_relative_tolerance
        )

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
    return mapping,kept,reduced_buses,nodes,ypu,voltage,power,aggregation_error,target_error,ybus_info


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
    center_tap_source=full_dataset/'center_tapped_transformer.csv'
    if center_tap_source.exists():
        for row in _rows(center_tap_source):
            bus=normalize_bus(row['bus_id'])
            if bus not in mapping:
                raise ConversionError(
                    f'center_tapped_transformer.csv references bus {bus!r} absent from assignment'
                )
            output.append([
                'center_tapped_transformer',row['transformer_id'],
                f'winding_{row["winding"]}',bus,mapping[bus],
            ])
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
            f'kvar={-s.imag*sbase_va/1000.0:.16g} model=1 vminpu=0 vmaxpu=2\n'
        )
    path.write_text(''.join(lines),encoding='ascii')


def _write_mapped_loads(path: Path, commands: list[str]) -> None:
    path.write_text(
        '! Solved physical loads moved to super-nodes; connection and phase pairs preserved.\n'
        +''.join(commands),
        encoding='ascii',
    )


def _write_commands(path: Path, heading: str, commands: list[str]) -> None:
    path.write_text(f'! {heading}\n'+''.join(commands),encoding='ascii')


def _validate_output(master: Path, nodes: list[Node], target_pu: np.ndarray, target_ypu: sparse.csr_matrix,
                     vbase: np.ndarray, sbase_va: float, output_csv: Path) -> dict[str,float]:
    dss.Basic.ClearAll();dss.Text.Command(f'compile "{master.resolve()}"');dss.Solution.Solve()
    if not dss.Solution.Converged():raise ConversionError('Generated reduced circuit did not converge')
    order=[Node(normalize_bus(x),int(x.rsplit('.',1)[1])) for x in dss.Circuit.YNodeOrder()]
    values=_complex_array(dss.Circuit.YNodeVArray());lookup=dict(zip(order,values))
    missing=[node for node in nodes if node not in lookup]
    if missing:raise ConversionError(f'Reduced circuit is missing nodes: {missing[:10]}')
    solved=np.asarray([lookup[node] for node in nodes])/vbase
    mag_error=np.abs(np.abs(solved)-np.abs(target_pu));complex_error=np.abs(solved-target_pu)
    energized=np.abs(target_pu)>=0.5
    node_index={node:i for i,node in enumerate(nodes)}
    passive_rows=[];passive_cols=[];passive_values=[]
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
                if node_b in node_index and local_y[a,b]!=0:
                    passive_rows.append(node_index[node_a]);passive_cols.append(node_index[node_b])
                    passive_values.append(local_y[a,b])
    passive=sparse.coo_matrix(
        (passive_values,(passive_rows,passive_cols)),shape=(len(nodes),len(nodes)),dtype=complex
    ).tocsr()
    passive.sum_duplicates();passive.eliminate_zeros();passive.sort_indices()
    actual_ypu=sparse.diags(vbase/sbase_va)@passive@sparse.diags(vbase)
    y_error=(actual_ypu-target_ypu).tocsr();y_error.eliminate_zeros()
    max_y_error=float(np.max(np.abs(y_error.data))) if y_error.nnz else 0.0
    with output_csv.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['bus_id','phase','target_v_mag_pu','solved_v_mag_pu','magnitude_error_pu','complex_error_pu'])
        for node,target,actual,me,ce in zip(nodes,target_pu,solved,mag_error,complex_error):
            writer.writerow([node.bus,NODE_TO_PHASE[node.number],f'{abs(target):.15f}',f'{abs(actual):.15f}',f'{me:.15f}',f'{ce:.15f}'])
    worst=int(np.argmax(mag_error))
    metrics={
        'max_voltage_magnitude_error_pu':float(mag_error[worst]),
        'max_voltage_magnitude_error_bus':nodes[worst].bus,
        'max_voltage_magnitude_error_phase':NODE_TO_PHASE[nodes[worst].number],
        'max_complex_voltage_error_pu':float(np.max(complex_error)),
        'max_abs_synthesized_ybus_error_pu':max_y_error,
        'relative_synthesized_ybus_error':float(
            max_y_error/max(np.max(np.abs(target_ypu.data)),np.finfo(float).tiny)
        ),
        'target_phase_nodes_below_0_5_pu':int(np.sum(~energized)),
    }
    if np.any(energized):
        energized_indices=np.flatnonzero(energized)
        worst_energized=int(energized_indices[np.argmax(mag_error[energized])])
        metrics.update({
            'max_energized_voltage_magnitude_error_pu':float(mag_error[worst_energized]),
            'max_energized_voltage_magnitude_error_bus':nodes[worst_energized].bus,
            'max_energized_voltage_magnitude_error_phase':NODE_TO_PHASE[nodes[worst_energized].number],
        })
    return metrics


def _convert_impl(full_master: Path, full_dataset: Path, reduced_dataset: Path, output: Path,
                  sbase_mva: float, matrix_tolerance: float, aggregation_tolerance: float,
                  load_representation: str, kron_relative_tolerance: float):
    mapping,kept,bus_rows,nodes,ypu,target_pu,power_pu,aggregation_error,target_error,ybus_info=(
        _validate_and_load(
            full_dataset,reduced_dataset,aggregation_tolerance,kron_relative_tolerance
        )
    )
    phase_shift_path=full_dataset/'phase_shift_equipment.csv'
    phase_shift_names={
        str(row['equipment_id']).strip().lower() for row in _rows(phase_shift_path)
        if int(row.get('enabled','1'))
    } if phase_shift_path.exists() else set()
    (
        source_name,source_bus,source_properties,frequency,bases,
        preserved_y_si,transformer_commands,load_commands,total_load_kva,original_equipment,
        physical_shunt_count,load_multiplier,voltage_level_audit,
    )=_compile_full(
        full_master,set(kept),nodes,mapping,load_representation,phase_shift_names
    )
    if source_bus not in kept:raise ConversionError(f'Source bus {source_bus!r} must survive')
    vbase=np.asarray([bases[node.bus] for node in nodes])
    sbase_va=sbase_mva*1e6
    ysi=sparse.diags(sbase_va/vbase)@ypu@sparse.diags(1.0/vbase)
    ysi=ysi.tocsr();ysi.eliminate_zeros();ysi.sort_indices()
    edges=_validate_sparse_tree(ysi,nodes,matrix_tolerance)
    synthetic_y=(ysi-preserved_y_si).tocsr()
    synthetic_y.eliminate_zeros();synthetic_y.sort_indices()
    branches,shunts,_=_decompose_sparse_radial_y(synthetic_y,nodes,matrix_tolerance)

    output.mkdir(parents=True,exist_ok=True)
    source_data=type('SourceData',(),{'source_properties':source_properties})()
    (output/'Source.dss').write_text(_source_command(source_data),encoding='ascii')
    (output/'Transformers.dss').write_text(
        '! Preserved physical transformers and solved regulator taps\n'+''.join(transformer_commands),
        encoding='ascii',
    )
    (output/'Branches.dss').write_text('! Radial Kron-equivalent series elements\n'+''.join(branches),encoding='ascii')
    (output/'Shunts.dss').write_text('! Kron-equivalent grounded shunts\n'+''.join(shunts),encoding='ascii')
    (
        model_commands,original_loads,injection_commands,physical_shunts,equipment_counts,
        load_models,equipment_transformation_audit,
    )=(
        original_equipment
    )
    if load_representation=='original-models':
        _write_commands(output/'Models.dss','Original supporting model objects',model_commands)
        _write_commands(output/'Loads.dss','Original load models relocated and rebased in per unit',original_loads)
        _write_commands(output/'Injections.dss','Original injection models relocated and rebased in per unit',injection_commands)
        _write_commands(output/'PhysicalShunts.dss','Original shunt models relocated and rebased in per unit',physical_shunts)
    elif load_representation=='mapped-elements':
        _write_commands(output/'Models.dss','No supporting model objects in constant-PQ mode',[])
        _write_mapped_loads(output/'Loads.dss',load_commands)
        _write_commands(output/'Injections.dss','No additional injection elements in constant-PQ mode',[])
        _write_commands(output/'PhysicalShunts.dss','Physical shunts remain embedded in the Kron Y-bus',[])
    elif load_representation=='nodal-wye':
        _write_commands(output/'Models.dss','No supporting model objects in nodal-wye mode',[])
        _write_loads(output/'Loads.dss',nodes,power_pu,vbase,source_bus,sbase_va)
        _write_commands(output/'Injections.dss','No additional injection elements in nodal-wye mode',[])
        _write_commands(output/'PhysicalShunts.dss','Physical shunts remain embedded in the Kron Y-bus',[])
    else:
        raise ConversionError(f'Unknown load representation {load_representation!r}')
    source_kv=source_properties.get('basekv','1')
    master=(
        'clear\n'
        f'new Circuit.Reduced basekv={source_kv} phases=3 bus1={source_bus} frequency={frequency:g}\n'
        'redirect Models.dss\nredirect Source.dss\nredirect Transformers.dss\nredirect Branches.dss\n'
        'redirect Shunts.dss\nredirect PhysicalShunts.dss\nredirect Loads.dss\nredirect Injections.dss\n'
        f'set loadmult={load_multiplier:.16g}\nset controlmode=off\nsolve mode=snapshot\n'
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
        'load_representation':load_representation,
        'voltage_level_audit':voltage_level_audit,
        'equipment_per_unit_transformation':equipment_transformation_audit,
        'reduced_ybus':ybus_info,
        'relocated_equipment_counts':equipment_counts if load_representation=='original-models' else {},
        'original_load_model_counts':load_models if load_representation=='original-models' else {},
        'relocated_physical_shunt_count':physical_shunt_count,
        'mapped_load_count':len(load_commands),
        'mapped_solved_load_p_pu':float(total_load_kva.real*1000.0/sbase_va),
        'mapped_solved_load_q_pu':float(total_load_kva.imag*1000.0/sbase_va),
        'sbase_mva':sbase_mva,'aggregation_max_abs_power_error_pu':aggregation_error,
        'target_vs_full_max_complex_voltage_error_pu':target_error,**metrics,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


def convert_reduced_dataset(full_master: Path, full_dataset: Path, reduced_dataset: Path,
                            output: Path, sbase_mva: float, matrix_tolerance: float=1e-9,
                            aggregation_tolerance: float=1e-12,
                            load_representation: str='original-models',
                            kron_relative_tolerance: float=1e-12):
    original=Path.cwd()
    try:
        return _convert_impl(
            full_master.resolve(),full_dataset.resolve(),reduced_dataset.resolve(),output.resolve(),
            sbase_mva,matrix_tolerance,aggregation_tolerance,load_representation,
            kron_relative_tolerance,
        )
    finally:
        os.chdir(original)
