# --------------------------------------------------------------------------- #
# Opti-KRON -- optimal Kron-based network reduction for distribution feeders.
#
# Two solver families sit on one shared core:
#
#   src/optimization/  MILP formulation (JuMP). Exact, and the method of the
#                      journal paper. One phase-aware model covers both phase
#                      counts -- single-phase is its degenerate case.
#   src/search/        exhaustive search at q=1, CPU and GPU. The method of the
#                      PSCC three-phase paper. Handles both phase counts.
#
# Both consume a `Network` from src/io and nothing else, so a new input format
# is one constructor rather than an edit to either solver.
# --------------------------------------------------------------------------- #
module OptiKRON

using LinearAlgebra
using Printf
using SparseArrays
using CSV
using DataFrames
using Graphs
using JuMP
const MOI = JuMP.MOI
using HiGHS

"""
    ReductionSolution

What a reduction backend returned. `MilpSolution` and `SearchSolution` both
subtype this, so [`Reduction`](@ref) can hold either. Declared here rather than
beside one of them because neither owns the concept.
"""
abstract type ReductionSolution end

# ---- Input layer ---------------------------------------------------------- #
include("io/network.jl")
include("io/from_csv.jl")
include("io/from_ybus.jl")
include("io/devices.jl")
include("io/to_matpower.jl")
include("io/to_csv.jl")

export Network, Branch, Device
export nnodes, nphase_rows, nscenarios, is_three_phase, has_branch_data, node_rows
export read_network_csv, read_network_ybus, network_from_matrices
export write_matpower, write_reduced_csv, read_voltage, read_devices, preserved_buses

# ---- Core: topology, Kron reduction, radialization ------------------------ #
include("core/topology.jl")
include("core/kron.jl")
include("core/radialization.jl")
include("core/powerflow.jl")

export RadialTree, orient_radial, interior_path, hop_distance, admissible_pairs
export nedges, path_to_root, subtree_degrees, phases_of
export expand_assignment, aggregate_injections, lift_voltages
export super_nodes, reduction_ratio, kron_reduce, reduced_adjacency
export identity_assignment, assign!
export critical_nodes, radialize, is_radial, spanning_subtree
export powerflow, powerflow_residual, slack_voltage

# ---- Optimization -------------------------------------------------------- #
include("optimization/solver.jl")
include("optimization/radiality.jl")
include("optimization/screening.jl")
include("optimization/milp.jl")
include("optimization/warmstart.jl")

# ---- Search: the exhaustive greedy at q = 1 -------------------------------- #
include("search/search.jl")
include("search/cpu.jl")
include("search/gpu.jl")

export search_reduce, SearchSolution, ReductionSolution, gpu_available

# ---- End-to-end pipeline --------------------------------------------------- #
include("pipeline.jl")
export optikron, Reduction, load_case, case_directory, list_cases
export enforced_violation, held_out_violation
export export_matpower, export_reduced

export select_optimizer, gurobi_available, reset_solver_cache!
export solve_milp, MilpSolution, annulus_violation, bus_impedance
export ScreeningReport, screen!, deviation_bound, screen_voltage_spread!
export add_radiality_constraints!
export zero_injection_buses, zero_injection_warmstart

# Both optional dependencies are loaded at module init, never lazily inside a
# solve: loading a module and calling into it within one dynamic call chain
# raises a world-age MethodError, and the symptom is a silent fallback rather
# than a crash. CUDA's absence is normal and stays quiet.
function __init__()
    load_gurobi!()
    load_cuda!()
    return nothing
end

end # module
