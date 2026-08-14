# Papers

The two papers behind this repository. Both are open-access preprints; the
canonical version of each is the published one, linked first.

## Journal — single-phase MILP, radialization

> O. Mokhtari, S. Chevalier, and M. Almassalkhi, "Structure-preserving Optimal
> Kron-based Reduction of Radial Distribution Networks," *Electric Power Systems
> Research*, 2026.

- Published: [doi:10.1016/j.epsr.2026.113615](https://doi.org/10.1016/j.epsr.2026.113615)
- Preprint: [arXiv:2508.15006](https://arxiv.org/abs/2508.15006)

Introduces the MILP formulation, the linearized voltage-magnitude error
constraint, and the radialization step (Theorem 1) implemented in
[`src/core/radialization.jl`](../../src/core/radialization.jl).

## Conference — three-phase, exhaustive search, GPU

> O. Mokhtari, S. Chevalier, and M. Almassalkhi, "Optimal Kron-based Reduction of
> Networks (Opti-KRON) for Three-phase Distribution Feeders," *24th Power Systems
> Computation Conference (PSCC)*, 2026.

- Preprint: [arXiv:2510.19608](https://arxiv.org/abs/2510.19608)

Extends the method to unbalanced feeders: the phase-availability rule
(φ_j ⊆ φ_i) that [`admissible_pairs`](../../src/core/topology.jl) enforces, and
the exhaustive-search solver.

## Where each part of the papers lives in the code

Each piece is its own file, so a reader who knows the papers can go straight to
the one they care about.

| Paper concept | File | Entry point |
|---|---|---|
| **The optimization problem** — assignment matrix `A`, objective, error constraint | [`src/optimization/milp.jl`](../../src/optimization/milp.jl) | `solve_milp` |
| **Radialization (Theorem 1)** — reduce first, reinsert critical buses after | [`src/core/radialization.jl`](../../src/core/radialization.jl) | `radialize`, `critical_nodes` |
| **Radiality as a constraint** — never pick a meshing assignment at all | [`src/optimization/radiality.jl`](../../src/optimization/radiality.jl) | `add_radiality_constraints!` |
| **Constraint screening** — kill hopeless pairs, skip rows that cannot bind | [`src/optimization/screening.jl`](../../src/optimization/screening.jl) | `screen!` |
| **Zero-injection warm start** — the exact, cheap part of the reduction | [`src/optimization/warmstart.jl`](../../src/optimization/warmstart.jl) | `zero_injection_warmstart` |
| The method end to end | [`src/pipeline.jl`](../../src/pipeline.jl) | `optikron` |
| Kron reduction, `A ⊗ I₃`, assignment algebra | [`src/core/kron.jl`](../../src/core/kron.jl) | `kron_reduce`, `expand_assignment` |
| Admissible pairs, hop limit, phase availability (9h) | [`src/core/topology.jl`](../../src/core/topology.jl) | `admissible_pairs` |
| Operating point the model linearises around | [`src/core/powerflow.jl`](../../src/core/powerflow.jl) | `powerflow` |
| Ground-truth accuracy check (the nonconvex annulus) | [`src/optimization/milp.jl`](../../src/optimization/milp.jl) | `annulus_violation` |

The two radiality routes are alternatives, not stages. `radialize` is the
journal paper's: solve, then repair whatever meshed, which is minimal *for that
assignment* but chosen without knowing a repair was coming.
`add_radiality_constraints!` puts the condition inside the MILP, so the solver
can trade a merge it would have had to undo for one it can keep.

**`:in_model` is the default, and not only for reduction quality.** Repairing
afterwards can leave the error budget behind. Reinserting a bus takes it out of
its cluster, which changes `(A-I)C` and therefore `e`, so the certified solution
is not the one you end up with. Usually the change helps — on case69 at Ē=0.001
the worst violation moves from -8.1e-06 to -1.3e-04 — but on case533mt it moves
from -1.4e-06 to **+3.4e-05**, outside the budget the MILP had just proved.
Enforcing radiality inside the model has nothing to perturb afterwards.

Reduction quality is close either way, and neither dominates: at Ē=0.001,
in-model wins on case69 (56.5% vs 55.1%) and case85 (48.2% vs 43.5%), post-hoc
wins narrowly on case141 (74.5% vs 73.0%). Enforcing radiality requires
`direction=:downstream`, which itself removes candidate merges — that is what
pays for the guarantee.

## What this implementation does differently

The code is not a transcription of either paper. Where the two disagree, this
repository follows whichever is better, and in a few places neither:

| | Papers | Here |
|---|---|---|
| Phase handling | separate single- and three-phase methods | one phase-aware model; single-phase is the degenerate case |
| Solve structure | iterative, with cutting planes and a Ω/A_prev decomposition | single-shot: one model, one solve |
| Error linearization | Appendix-B big-M bounds | JuMP indicator constraints — no constant to size, no relaxation slack |
| Impedances | reconstructed per-branch as Π diag(z) Γ | slack-referenced inverse of the given Ybus directly |
| Objective | reduction traded against error via α | reduction alone; error is a hard constraint |

See the header of [`src/optimization/milp.jl`](../../src/optimization/milp.jl)
for the reasoning behind each.

## Citation

[`CITATION.cff`](../../CITATION.cff) records how to cite the software. If you use
the *methodology*, cite the papers above.
