from __future__ import annotations

import csv
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import opendssdirect as dss
from scipy.linalg import pinv


class ConversionError(RuntimeError):
    pass


PC_CLASSES = {
    "load", "generator", "pvsystem", "storage", "isource", "vccs",
    "indmach012", "gicsource",
}
SOURCE_CLASSES = {"vsource"}
SKIP_CLASSES = PC_CLASSES | SOURCE_CLASSES | {"energymeter", "monitor", "sensor"}
PHASE_TO_NODE = {"a": 1, "b": 2, "c": 3}
NODE_TO_PHASE = {1: "a", 2: "b", 3: "c"}


@dataclass(frozen=True)
class Node:
    bus: str
    number: int


@dataclass
class CircuitData:
    buses: list[str]
    phases: dict[str, tuple[int, ...]]
    nodes: list[Node]
    y: np.ndarray
    voltage: np.ndarray
    injection: np.ndarray
    preserved_y: np.ndarray
    transformer_commands: list[str]
    transformer_buses: set[str]
    source_name: str
    source_bus: str
    source_properties: dict[str, str]
    frequency: float


def normalize_bus(value: str) -> str:
    return value.strip().strip('"').split(".")[0].lower()


def read_bus_order(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.reader(stream))
    if not rows:
        raise ConversionError(f"Empty bus-order file: {path}")
    if rows[0] and rows[0][0].strip().lower() in {"bus", "bus_id", "name"}:
        rows = rows[1:]
    buses = [normalize_bus(row[0]) for row in rows if row and row[0].strip()]
    if len(buses) != len(set(buses)):
        raise ConversionError("bus_order.csv contains duplicate bus names")
    return buses


def read_assignment(path: Path) -> np.ndarray:
    if path.suffix.lower() == ".npy":
        return np.asarray(np.load(path), dtype=float)
    rows: list[list[float]] = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for raw in csv.reader(stream):
            if not raw or all(not cell.strip() for cell in raw):
                continue
            try:
                rows.append([float(cell) for cell in raw])
            except ValueError:
                if not rows:
                    continue
                # Permit a row-label column.
                rows.append([float(cell) for cell in raw[1:]])
    matrix = np.asarray(rows, dtype=float)
    if matrix.ndim != 2:
        raise ConversionError("Assignment matrix must be two-dimensional")
    if matrix.shape[1] == matrix.shape[0] + 1:
        matrix = matrix[:, 1:]
    return matrix


def _complex_array(values: list[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.size == 1 and array[0] == 0:
        return np.asarray([], dtype=complex)
    if array.size % 2:
        raise ConversionError("OpenDSS returned an odd-length complex array")
    return array[0::2] + 1j * array[1::2]


def _element_nodes(bus_specs: list[str], node_order: list[int], nconds: int) -> list[Node | None]:
    result: list[Node | None] = []
    for terminal, spec in enumerate(bus_specs):
        parts = spec.strip().split(".")
        bus = normalize_bus(parts[0])
        explicit = [int(x) for x in parts[1:] if x.lstrip("-").isdigit()]
        for conductor in range(nconds):
            idx = terminal * nconds + conductor
            number = node_order[idx] if idx < len(node_order) else (explicit[conductor] if conductor < len(explicit) else conductor + 1)
            result.append(None if number <= 0 else Node(bus, number))
    return result


def _active_element_class() -> str:
    return dss.CktElement.Name().split(".", 1)[0].lower()


def _source_properties() -> tuple[str, str, dict[str, str]]:
    names = dss.Vsources.AllNames()
    if len(names) != 1:
        raise ConversionError(f"Exactly one Vsource is currently supported; found {len(names)}")
    dss.Circuit.SetActiveElement(f"Vsource.{names[0]}")
    props: dict[str, str] = {}
    for key in ("bus1", "bus2", "phases", "basekv", "pu", "angle", "frequency", "r1", "x1", "r0", "x0"):
        value = dss.Properties.Value(key)
        if value not in (None, ""):
            props[key] = str(value)
    return names[0], normalize_bus(props["bus1"]), props


def _transformer_command(element_name: str) -> str:
    cls, name = element_name.split(".", 1)
    keys = ("phases", "windings", "buses", "conns", "kvs", "kvas", "taps", "%rs", "xhl", "xht", "xlt", "%noloadloss", "%imag", "ppm_antifloat")
    values = []
    for key in keys:
        value = dss.Properties.Value(key)
        if value not in (None, ""):
            values.append(f"{key}={value}")
    return f"new {cls}.{_safe_name(name)} " + " ".join(values) + "\n"


def import_circuit(master: Path, ordered_buses: list[str]) -> CircuitData:
    dss.Basic.ClearAll()
    dss.Text.Command(f'compile "{master.resolve()}"')
    dss.Solution.Solve()
    if not dss.Solution.Converged():
        raise ConversionError("The full OpenDSS circuit did not converge")

    y_order = [entry.lower() for entry in dss.Circuit.YNodeOrder()]
    nodes = [Node(normalize_bus(entry), int(entry.rsplit(".", 1)[1])) for entry in y_order]
    node_index = {node: i for i, node in enumerate(nodes)}
    if len(node_index) != len(nodes):
        raise ConversionError("OpenDSS Y-node order contains duplicate node names")
    voltage = _complex_array(dss.Circuit.YNodeVArray())
    if voltage.size != len(nodes):
        raise ConversionError("Voltage and Y-node dimensions differ")

    phases: dict[str, tuple[int, ...]] = {}
    available = {node.bus for node in nodes}
    if set(ordered_buses) != available:
        missing = sorted(available - set(ordered_buses))
        extra = sorted(set(ordered_buses) - available)
        raise ConversionError(f"bus_order.csv mismatch; missing={missing[:10]}, extra={extra[:10]}")
    for bus in ordered_buses:
        phases[bus] = tuple(sorted(node.number for node in nodes if node.bus == bus and node.number in (1, 2, 3)))
        if not phases[bus]:
            raise ConversionError(f"Bus {bus!r} has no active phase nodes 1/2/3")

    n = len(nodes)
    y = np.zeros((n, n), dtype=complex)
    preserved_y = np.zeros((n, n), dtype=complex)
    injection = np.zeros(n, dtype=complex)
    transformer_commands: list[str] = []
    transformer_buses: set[str] = set()

    element_names = list(dss.Circuit.AllElementNames())
    for element_name in element_names:
        dss.Circuit.SetActiveElement(element_name)
        cls = _active_element_class()
        nterms = dss.CktElement.NumTerminals()
        nconds = dss.CktElement.NumConductors()
        local_nodes = _element_nodes(dss.CktElement.BusNames(), dss.CktElement.NodeOrder(), nconds)

        if cls in PC_CLASSES:
            currents = _complex_array(dss.CktElement.Currents())
            for terminal in range(nterms):
                for conductor in range(nconds):
                    local = terminal * nconds + conductor
                    node = local_nodes[local]
                    if node in node_index and local < currents.size:
                        injection[node_index[node]] -= currents[local]
            continue
        if cls in SKIP_CLASSES:
            continue

        yprim = _complex_array(dss.CktElement.YPrim())
        order = nterms * nconds
        if yprim.size != order * order:
            # Controls and other non-electrical objects may report no primitive.
            if yprim.size == 0:
                continue
            raise ConversionError(f"Unexpected YPrim size for {element_name}")
        local_y = yprim.reshape((order, order), order="F")
        element_y = np.zeros((n, n), dtype=complex)
        for a, node_a in enumerate(local_nodes):
            if node_a not in node_index:
                continue
            ia = node_index[node_a]
            for b, node_b in enumerate(local_nodes):
                if node_b in node_index:
                    element_y[ia, node_index[node_b]] += local_y[a, b]
        y += element_y
        if cls in {"transformer", "autotrans"}:
            preserved_y += element_y
            transformer_commands.append(_transformer_command(element_name))
            transformer_buses.update(normalize_bus(spec) for spec in dss.CktElement.BusNames())

    source_name, source_bus, source_properties = _source_properties()
    return CircuitData(
        buses=ordered_buses, phases=phases, nodes=nodes, y=y, voltage=voltage,
        injection=injection, preserved_y=preserved_y,
        transformer_commands=transformer_commands, transformer_buses=transformer_buses,
        source_name=source_name, source_bus=source_bus,
        source_properties=source_properties, frequency=float(dss.Solution.Frequency()),
    )


def validate_assignment(a: np.ndarray, data: CircuitData, tolerance: float = 1e-9) -> list[int]:
    n = len(data.buses)
    if a.shape != (n, n):
        raise ConversionError(f"A has shape {a.shape}; expected {(n, n)}")
    if np.max(np.abs(a - np.round(a))) > tolerance or np.any((a < -tolerance) | (a > 1+tolerance)):
        raise ConversionError("A must be binary")
    a[:] = np.round(a)
    if np.max(np.abs(a.sum(axis=0) - 1)) > tolerance:
        raise ConversionError("Every column of A must contain exactly one assignment")
    survivors = [i for i in range(n) if a[i, i] == 1]
    if not survivors:
        raise ConversionError("A has no surviving buses")
    for row, col in zip(*np.nonzero(a)):
        if row not in survivors:
            raise ConversionError(f"Bus {data.buses[col]} is assigned to non-survivor {data.buses[row]}")
        if not set(data.phases[data.buses[col]]).issubset(data.phases[data.buses[row]]):
            raise ConversionError(f"Phase incompatibility: {data.buses[col]} -> {data.buses[row]}")
    source_idx = data.buses.index(data.source_bus)
    if source_idx not in survivors:
        raise ConversionError(f"Source bus {data.source_bus!r} must survive")
    for bus in sorted(data.transformer_buses):
        idx = data.buses.index(bus)
        if idx not in survivors or a[idx, idx] != 1:
            raise ConversionError(f"Transformer/regulator terminal bus {bus!r} must survive and represent itself")
    return survivors


def lift_assignment(a: np.ndarray, data: CircuitData, survivors: list[int]) -> tuple[np.ndarray, list[Node], list[int]]:
    node_index = {node: i for i, node in enumerate(data.nodes)}
    kept_nodes = [Node(data.buses[i], phase) for i in survivors for phase in data.phases[data.buses[i]]]
    lift = np.zeros((len(kept_nodes), len(data.nodes)), dtype=float)
    kept_pos = {node: i for i, node in enumerate(kept_nodes)}
    for child_idx, child_bus in enumerate(data.buses):
        parent_idx = int(np.flatnonzero(a[:, child_idx])[0])
        parent_bus = data.buses[parent_idx]
        for phase in data.phases[child_bus]:
            lift[kept_pos[Node(parent_bus, phase)], node_index[Node(child_bus, phase)]] = 1.0
    kept_original = [node_index[node] for node in kept_nodes]
    return lift, kept_nodes, kept_original


def kron_reduce(y: np.ndarray, kept: list[int]) -> np.ndarray:
    reduced = [i for i in range(y.shape[0]) if i not in set(kept)]
    ykk = y[np.ix_(kept, kept)]
    if not reduced:
        return ykk.copy()
    ykr = y[np.ix_(kept, reduced)]
    yrk = y[np.ix_(reduced, kept)]
    yrr = y[np.ix_(reduced, reduced)]
    return ykk - ykr @ pinv(yrr) @ yrk


def _safe_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_]", "_", value)
    return clean or "bus"


def _matrix_text(matrix: np.ndarray) -> str:
    rows = []
    for i in range(matrix.shape[0]):
        rows.append(" ".join(f"{matrix[i,j]:.16g}" for j in range(i + 1)))
    return " | ".join(rows)


def _bus_spec(nodes: list[Node]) -> str:
    if not nodes:
        raise ConversionError("Cannot write an empty terminal")
    return nodes[0].bus + "".join(f".{node.number}" for node in nodes)


def _write_reactor(name: str, left: list[Node], right: list[Node] | None, y_matrix: np.ndarray, tolerance: float) -> str | None:
    if np.linalg.norm(y_matrix, ord=np.inf) <= tolerance:
        return None
    if y_matrix.shape[0] != y_matrix.shape[1] or y_matrix.shape[0] != len(left):
        raise ConversionError(f"Equivalent {name} has a non-square phase coupling")
    z = pinv(y_matrix)
    reconstruction = np.linalg.norm(y_matrix @ z @ y_matrix - y_matrix, ord=np.inf)
    if reconstruction > max(tolerance, 1e-7 * np.linalg.norm(y_matrix, ord=np.inf)):
        raise ConversionError(f"Equivalent {name} admittance cannot be represented by a Reactor Z matrix")
    bus1 = _bus_spec(left)
    bus2 = _bus_spec(right) if right else left[0].bus + "." + ".".join("0" for _ in left)
    return (
        f"new Reactor.{name} phases={len(left)} bus1={bus1} bus2={bus2} "
        f"rmatrix=[{_matrix_text(z.real)}] xmatrix=[{_matrix_text(z.imag)}]\n"
    )


def _graph_edges(y: np.ndarray, nodes: list[Node], tolerance: float) -> list[tuple[int,int]]:
    bus_names = list(dict.fromkeys(node.bus for node in nodes))
    by_bus = {bus: [i for i,node in enumerate(nodes) if node.bus == bus] for bus in bus_names}
    return [(i,j) for i in range(len(bus_names)) for j in range(i+1,len(bus_names))
            if np.max(np.abs(y[np.ix_(by_bus[bus_names[i]],by_bus[bus_names[j]])])) > tolerance]


def _validate_tree(y: np.ndarray, nodes: list[Node], tolerance: float) -> None:
    bus_names = list(dict.fromkeys(node.bus for node in nodes))
    edges = _graph_edges(y,nodes,tolerance)
    if len(edges) != max(0,len(bus_names)-1):
        raise ConversionError(f"Reduced admittance graph is not a tree: {len(bus_names)} buses, {len(edges)} edges")
    seen={0} if bus_names else set()
    changed=True
    while changed:
        changed=False
        for i,j in edges:
            if i in seen and j not in seen: seen.add(j); changed=True
            if j in seen and i not in seen: seen.add(i); changed=True
    if len(seen)!=len(bus_names):
        raise ConversionError("Reduced admittance graph is disconnected")


def decompose_radial_y(y: np.ndarray, nodes: list[Node], tolerance: float) -> tuple[list[str], list[str], list[tuple[int,int]]]:
    bus_names = list(dict.fromkeys(node.bus for node in nodes))
    by_bus = {bus: [i for i, node in enumerate(nodes) if node.bus == bus] for bus in bus_names}
    branch_lines: list[str] = []
    edges: list[tuple[int,int]] = []
    incident = {bus: np.zeros((len(by_bus[bus]), len(by_bus[bus])), complex) for bus in bus_names}

    for i, bus_i in enumerate(bus_names):
        rows = by_bus[bus_i]
        for j in range(i + 1, len(bus_names)):
            bus_j = bus_names[j]
            cols = by_bus[bus_j]
            block = y[np.ix_(rows, cols)]
            reverse = y[np.ix_(cols, rows)]
            # Kron arithmetic can leave tiny cross-phase terms beside a
            # dominant one-phase edge. OpenDSS Reactor terminals require the
            # same active conductor count at both ends, so classify active
            # conductors relative to the edge's dominant coupling.
            block_scale=float(np.max(np.abs(block))) if block.size else 0.0
            coupling_tolerance=max(tolerance,1e-7*block_scale)
            active_rows = np.flatnonzero(np.max(np.abs(block), axis=1) > coupling_tolerance)
            active_cols = np.flatnonzero(np.max(np.abs(block), axis=0) > coupling_tolerance)
            if not active_rows.size and not active_cols.size:
                continue
            if active_rows.size != active_cols.size:
                raise ConversionError(f"Edge {bus_i}-{bus_j} has incompatible phase dimensions")
            sub = block[np.ix_(active_rows, active_cols)]
            reverse_sub = reverse[np.ix_(active_cols, active_rows)]
            if np.linalg.norm(sub - reverse_sub.T, ord=np.inf) > max(tolerance, 1e-7*np.linalg.norm(sub,ord=np.inf)):
                raise ConversionError(f"Edge {bus_i}-{bus_j} is non-reciprocal and cannot be written as a Reactor")
            ys = -sub
            left = [nodes[rows[k]] for k in active_rows]
            right = [nodes[cols[k]] for k in active_cols]
            text = _write_reactor(f"kr_{i+1}_{j+1}", left, right, ys, tolerance)
            if text:
                branch_lines.append(text)
                edges.append((i, j))
                incident[bus_i][np.ix_(active_rows, active_rows)] += ys
                incident[bus_j][np.ix_(active_cols, active_cols)] += ys.T

    shunt_lines: list[str] = []
    for i, bus in enumerate(bus_names):
        idx = by_bus[bus]
        shunt = y[np.ix_(idx, idx)] - incident[bus]
        text = _write_reactor(f"sh_{i+1}", [nodes[k] for k in idx], None, shunt, tolerance)
        if text:
            shunt_lines.append(text)
    return branch_lines, shunt_lines, edges


def _source_command(data: CircuitData) -> str:
    p = data.source_properties
    allowed = ("bus1", "bus2", "phases", "basekv", "pu", "angle", "frequency", "r1", "x1", "r0", "x0")
    values = " ".join(f"{key}={p[key]}" for key in allowed if key in p)
    return f"edit Vsource.source {values}\n"


def _write_equivalent_injections(path: Path, nodes: list[Node], voltage: np.ndarray, injection: np.ndarray) -> None:
    lines = ["! Operating-point constant-PQ equivalents; positive kW/kvar means consumption.\n"]
    for i, (node, v, current) in enumerate(zip(nodes, voltage, injection), start=1):
        s_inj = v * np.conj(current)
        if abs(s_inj) < 1e-6:
            continue
        kv_ln = abs(v) / 1000.0
        lines.append(
            f"new Load.eq_{i} phases=1 bus1={node.bus}.{node.number} conn=wye "
            f"kv={kv_ln:.16g} kw={-s_inj.real/1000:.16g} kvar={-s_inj.imag/1000:.16g} model=1\n"
        )
    path.write_text("".join(lines), encoding="ascii")


def write_reduced_circuit(output: Path, data: CircuitData, kept_nodes: list[Node], yred: np.ndarray, ysynth: np.ndarray,
                          vred: np.ndarray, ired: np.ndarray, tolerance: float) -> None:
    output.mkdir(parents=True, exist_ok=True)
    _validate_tree(yred, kept_nodes, tolerance)
    branches, shunts, _ = decompose_radial_y(ysynth, kept_nodes, tolerance)
    (output / "Source.dss").write_text(_source_command(data), encoding="ascii")
    (output / "Transformers.dss").write_text("! Preserved transformers and fixed regulator taps\n" + "".join(data.transformer_commands), encoding="ascii")
    (output / "Branches.dss").write_text("! Kron-equivalent series elements\n" + "".join(branches), encoding="ascii")
    (output / "Shunts.dss").write_text("! Kron-equivalent grounded shunts\n" + "".join(shunts), encoding="ascii")
    _write_equivalent_injections(output / "Loads.dss", kept_nodes, vred, ired)
    source_kv = data.source_properties.get("basekv", "1")
    master = (
        "clear\n"
        f"new Circuit.Reduced basekv={source_kv} phases=3 bus1={data.source_bus} frequency={data.frequency:g}\n"
        "redirect Source.dss\nredirect Transformers.dss\nredirect Branches.dss\nredirect Shunts.dss\nredirect Loads.dss\n"
        "set controlmode=off\nsolve mode=snapshot\n"
    )
    (output / "Master.dss").write_text(master, encoding="ascii")


def _write_mapping(path: Path, a: np.ndarray, data: CircuitData) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["original_bus", "super_bus", "survives"])
        for j, bus in enumerate(data.buses):
            i = int(np.flatnonzero(a[:, j])[0])
            writer.writerow([bus, data.buses[i], int(a[j,j] == 1)])


def validate_output(master: Path, kept_nodes: list[Node], target: np.ndarray, output_csv: Path) -> dict[str,float]:
    dss.Basic.ClearAll()
    dss.Text.Command(f'compile "{master.resolve()}"')
    dss.Solution.Solve()
    if not dss.Solution.Converged():
        raise ConversionError("Generated reduced circuit did not converge")
    order = [Node(normalize_bus(x), int(x.rsplit(".",1)[1])) for x in dss.Circuit.YNodeOrder()]
    values = _complex_array(dss.Circuit.YNodeVArray())
    lookup = dict(zip(order, values))
    missing = [node for node in kept_nodes if node not in lookup]
    if missing:
        raise ConversionError(f"Reduced circuit is missing nodes: {missing[:10]}")
    actual = np.asarray([lookup[node] for node in kept_nodes])
    error = actual - target
    with output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["bus_id","phase","target_re_v","target_im_v","reduced_re_v","reduced_im_v","complex_error_v","magnitude_error_pu"])
        for node, expected, obtained, delta in zip(kept_nodes,target,actual,error):
            base = max(abs(expected), 1.0)
            writer.writerow([node.bus,NODE_TO_PHASE.get(node.number,str(node.number)),expected.real,expected.imag,obtained.real,obtained.imag,abs(delta),abs(abs(obtained)-abs(expected))/base])
    return {
        "max_complex_error_v": float(np.max(np.abs(error))),
        "max_magnitude_error_pu": float(np.max(np.abs(np.abs(actual)-np.abs(target))/np.maximum(np.abs(target),1.0))),
    }


def _convert_impl(master: Path, assignment_path: Path, bus_order_path: Path, output: Path,
            tolerance: float = 1e-9, max_voltage_error_pu: float = 1e-6) -> dict[str, object]:
    master = master.resolve()
    assignment_path = assignment_path.resolve()
    bus_order_path = bus_order_path.resolve()
    output = output.resolve()
    buses = read_bus_order(bus_order_path)
    assignment = read_assignment(assignment_path)
    data = import_circuit(master, buses)
    survivors = validate_assignment(assignment, data)
    lift, kept_nodes, kept_original = lift_assignment(assignment, data, survivors)
    ysynth = kron_reduce(data.y - data.preserved_y, kept_original)
    ypreserved = data.preserved_y[np.ix_(kept_original, kept_original)]
    yred = ysynth + ypreserved
    vred = data.voltage[kept_original]
    # Fitted operating-point injections make the original surviving voltage
    # vector an exact solution of the Kron network. Direct A*I aggregation is
    # the Opti-KRON approximation and generally does not preserve V_K exactly.
    ired = yred @ vred
    source_positions = [i for i,node in enumerate(kept_nodes) if node.bus == data.source_bus]
    for i in source_positions:
        ired[i] = 0.0  # the retained Vsource supplies the slack injection
    write_reduced_circuit(output, data, kept_nodes, yred, ysynth, vred, ired, tolerance)
    _write_mapping(output / "bus_mapping.csv", assignment, data)
    metrics = validate_output(output / "Master.dss", kept_nodes, vred, output / "validation.csv")
    report: dict[str, object] = {
        "full_bus_count": len(data.buses), "reduced_bus_count": len(survivors),
        "full_phase_node_count": len(data.nodes), "reduced_phase_node_count": len(kept_nodes),
        **metrics,
    }
    (output / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if metrics["max_magnitude_error_pu"] > max_voltage_error_pu:
        raise ConversionError(
            f"Reduced circuit voltage error {metrics['max_magnitude_error_pu']:.6g} pu exceeds "
            f"limit {max_voltage_error_pu:.6g} pu; output retained for diagnosis"
        )
    return report


def convert(master: Path, assignment_path: Path, bus_order_path: Path, output: Path,
            tolerance: float = 1e-9, max_voltage_error_pu: float = 1e-6) -> dict[str, object]:
    """Convert and validate a circuit while insulating the caller from DSS cwd changes."""
    original_directory = Path.cwd()
    try:
        return _convert_impl(master, assignment_path, bus_order_path, output, tolerance, max_voltage_error_pu)
    finally:
        os.chdir(original_directory)
