using JuMP
using LinearAlgebra

const TEST3PH = joinpath(@__DIR__, "..", "data", "test3ph")

@testset "OptiKRON.milp" begin

    net4 = read_network_csv(TEST4)
    V4 = powerflow(net4)
    net3 = read_network_csv(TEST3PH)
    V3 = powerflow(net3)

    @testset "bus impedance" begin
        Z = bus_impedance(net4)
        nph = nphase_rows(net4)
        slack_rows = node_rows(net4)[net4.slack]
        other = setdiff(1:nph, slack_rows)

        # The slack neither moves nor propagates.
        @test all(iszero, Z[slack_rows, :])
        @test all(iszero, Z[:, slack_rows])
        @test Z[other, other] * Matrix(net4.Ybus[other, other]) ≈ I

        # test4 comes from branch data with no shunts, so the full Ybus is a
        # singular Laplacian -- the case that makes inv(Ybus) unusable and
        # bus_impedance necessary.
        @test abs(sum(Matrix(net4.Ybus)[2, :])) < 1e-9
        @test rank(Matrix(net4.Ybus)) < nph

        # Z reproduces the exact error model: hold the slack, solve the rest.
        A = identity_assignment(net4)
        A[1, 2] = 1.0                                  # bus 2 joins the slack
        A[2, 2] = 0.0
        C = conj.(net4.S ./ V4)
        shift = (expand_assignment(A, net4) - I) * C

        # Aggregation moves current between buses without creating any. This is
        # the property the whole error model leans on.
        @test maximum(abs.(sum(shift, dims=1))) < 1e-12

        exact = zeros(ComplexF64, nph, size(shift, 2))
        exact[other, :] = Matrix(net4.Ybus[other, other]) \ shift[other, :]
        @test Z * shift ≈ exact
    end

    @testset "a tight budget forbids every merge" begin
        for approach in (:R, :S), (net, V) in ((net4, V4), (net3, V3))
            sol = solve_milp(net, V, 1e-9; approach=approach, hops=3)
            @test sol.A == identity_assignment(net)
            @test length(sol.kept) == nnodes(net)
            @test sol.reduction == 0.0
            @test sol.status == MOI.OPTIMAL
        end
    end

    @testset "a loose budget collapses the feeder" begin
        sol = solve_milp(net4, V4, 1.0; hops=3)
        @test sol.kept == [net4.slack]
        @test sol.reduction ≈ 0.75
    end

    @testset "solutions satisfy the true nonconvex annulus" begin
        # The model only ever sees a linearisation; what has to hold is the
        # original constraint, re-derived independently from the sparse Ybus.
        for Ē in (0.002, 0.005, 0.02), (net, V) in ((net4, V4), (net3, V3))
            sol = solve_milp(net, V, Ē; hops=3)
            @test annulus_violation(net, sol.A, V, Ē) <= 0
        end
    end

    @testset "both linearisations agree here" begin
        for Ē in (0.002, 0.01), (net, V) in ((net4, V4), (net3, V3))
            r = solve_milp(net, V, Ē; approach=:R, hops=3)
            s = solve_milp(net, V, Ē; approach=:S, hops=3)
            @test length(r.kept) == length(s.kept)
        end
    end

    @testset "structural invariants" begin
        sol = solve_milp(net3, V3, 0.01; hops=3)
        A = sol.A
        B = nnodes(net3)

        @test all(sum(A[:, j]) ≈ 1 for j in 1:B)          # represented exactly once
        @test A[net3.slack, net3.slack] == 1              # slack always survives
        @test all(A[i, j] <= A[i, i] for i in 1:B, j in 1:B)
        @test all(x -> x in (0.0, 1.0), A)

        # Clusters are connected: everything between a bus and its super-node
        # belongs to the same cluster.
        tree = orient_radial(net3)
        for j in 1:B, i in 1:B
            (A[i, j] == 1 && i != j) || continue
            for k in interior_path(tree, i, j)
                @test A[i, k] == 1
            end
        end
    end

    @testset "phase availability is respected" begin
        # test3ph is N1(abc) - N2(abc) - N3(a). N3 may join either abc bus, but
        # neither abc bus may ever join the single-phase N3.
        sol = solve_milp(net3, V3, 1.0; hops=3)
        for j in 1:nnodes(net3), i in 1:nnodes(net3)
            sol.A[i, j] == 1 || continue
            for p in 1:3
                @test !(net3.phases[p, j] && !net3.phases[p, i])
            end
        end
        @test sol.A[3, 1] == 0 && sol.A[3, 2] == 0       # N3 cannot represent an abc bus
    end

    @testset "a looser budget never keeps more buses" begin
        counts = [length(solve_milp(net3, V3, E; hops=3).kept)
                  for E in (0.0005, 0.002, 0.01, 0.05)]
        @test issorted(counts, rev=true)
    end

    @testset "hops bounds the cluster radius" begin
        tree = orient_radial(net4)
        sol = solve_milp(net4, V4, 1.0; hops=1)
        for j in 1:nnodes(net4), i in 1:nnodes(net4)
            (sol.A[i, j] == 1 && i != j) || continue
            @test hop_distance(tree, i, j) <= 1
        end
        # At one hop the star's outer buses cannot all reach the slack.
        @test length(sol.kept) > 1
    end

    @testset "scenario selection" begin
        sol = solve_milp(net4, V4, 0.01; scenarios=[2], hops=3)
        @test sol.scenarios == [2]
        @test annulus_violation(net4, sol.A, V4, 0.01; scenarios=[2]) <= 0

        # A budget enforced on one scenario says nothing about the others; the
        # returned solution records which ones it was certified against.
        @test length(sol.scenarios) < nscenarios(net4)
    end

    @testset "reported metadata" begin
        sol = solve_milp(net4, V4, 0.01; hops=2, screen=false)
        @test sol.nbinaries > 0
        @test sol.nannulus == 2 * sol.nbinaries * length(sol.scenarios)
        @test sol.screening === nothing
        @test sol.objective ≈ length(sol.kept)
        @test sol.solver in (:gurobi, :highs)
        @test sol.reduction ≈ reduction_ratio(sol.A)
        @test sol.gap >= 0
    end

    @testset "a reduction radializes" begin
        # The end-to-end path the papers describe: reduce, then recover a radial
        # equivalent. Kron reduction meshes a radial feeder (Lemma 1), so this
        # is where the MILP's output meets radialization.
        net = net4
        sol = solve_milp(net, V4, 0.02; hops=3, max_reduction=0.5)
        A_radial, added = radialize(net, sol.A)

        @test is_radial(net, A_radial)
        @test all(A_radial[b, b] == 1 for b in added)
        @test reduction_ratio(A_radial) <= reduction_ratio(sol.A)

        # Reinserted buses take back their own injection, so nothing that stayed
        # a super-node changes hands.
        for b in super_nodes(sol.A)
            @test A_radial[b, b] == 1
        end
    end

    @testset "screening does not change the answer" begin
        # The screen removes only pairs that no assignment could have used and
        # rows that could never bind, so the optimum must survive it untouched.
        for Ē in (0.002, 0.005, 0.02), (net, V) in ((net4, V4), (net3, V3))
            off = solve_milp(net, V, Ē; hops=3, screen=false)
            on = solve_milp(net, V, Ē; hops=3, screen=true)
            @test length(on.kept) == length(off.kept)
            @test annulus_violation(net, on.A, V, Ē) <= 0
            @test on.nannulus <= off.nannulus
            @test on.nbinaries <= off.nbinaries
        end
    end

    @testset "screening reports what it removed" begin
        sol = solve_milp(net3, V3, 0.002; hops=3)
        report = sol.screening
        @test report isa ScreeningReport
        @test report.pairs_after == sol.nbinaries
        @test report.pairs_after <= report.pairs_before
        @test report.rows_kept + report.rows_skipped >= report.rows_kept
        @test report.passes >= 1
        @test size(report.deviation) == (nphase_rows(net3), length(sol.scenarios))
        @test all(report.deviation .>= 0)
    end

    @testset "the deviation bound really bounds the error" begin
        # dev must cover |e| for *every* assignment the mask allows, so a sample
        # of legal assignments must never exceed it.
        net, V, Ē = net3, V3, 0.05
        tree = orient_radial(net)
        admissible = admissible_pairs(net, tree; hops=3)
        C = conj.(net.S ./ V)
        Z = bus_impedance(net)
        dev = deviation_bound(Z, C, admissible, net, [1])

        for (i, j) in [(1, 2), (2, 3), (1, 3), (1, 2)]
            admissible[i, j] || continue
            A = identity_assignment(net)
            A[j, j] == 1 && A[i, i] == 1 && i != j && assign!(A, i, j)
            e = Z * ((expand_assignment(A, net) - I) * C)
            @test all(abs.(e[:, 1]) .<= dev[:, 1] .+ 1e-12)
        end
    end

    @testset "radiality can be enforced in the model" begin
        # The alternative to reduce-then-repair: the solver never picks an
        # assignment that would mesh the reduced network in the first place.
        sol = solve_milp(net4, V4, 0.02; hops=3, direction=:downstream,
            enforce_radiality=true, max_reduction=0.5)
        @test is_radial(net4, sol.A)
        @test radialize(net4, sol.A)[2] == Int[]     # nothing left to repair
        @test annulus_violation(net4, sol.A, V4, 0.02) <= 0
    end

    @testset "zero-injection warm start" begin
        net, V, Ē = net3, V3, 0.05
        tree = orient_radial(net)
        admissible = admissible_pairs(net, tree; hops=3)
        quiet = zero_injection_buses(net)

        @test !(net.slack in quiet)                     # the slack never counts
        for bus in quiet, s in 1:nscenarios(net)
            @test all(iszero, net.S[node_rows(net)[bus], s])
        end

        W = zero_injection_warmstart(net, V, Ē, admissible, tree, 1:nscenarios(net))
        B = nnodes(net)
        @test all(sum(W[:, j]) ≈ 1 for j in 1:B)        # a valid assignment
        @test W[net.slack, net.slack] == 1
        @test all(x -> x in (0.0, 1.0), W)

        # Only zero-injection buses may be absorbed, and merging them perturbs
        # nothing -- so the start is feasible for the real constraint.
        for j in 1:B, i in 1:B
            (W[i, j] == 1 && i != j) || continue
            @test j in quiet
        end
        @test annulus_violation(net, W, V, Ē) <= 0

        # e = Z(A-I)C is exactly zero for a zero-injection merge.
        shift = (expand_assignment(W, net) - I) * conj.(net.S ./ V)
        @test maximum(abs, bus_impedance(net) * shift) < 1e-12
    end

    @testset "warm start does not change the optimum" begin
        for Ē in (0.002, 0.02), (net, V) in ((net4, V4), (net3, V3))
            zi = solve_milp(net, V, Ē; hops=3, warm_start=:zero_injection)
            id = solve_milp(net, V, Ē; hops=3, warm_start=:identity)
            @test length(zi.kept) == length(id.kept)
            @test annulus_violation(net, zi.A, V, Ē) <= 0
        end
        @test_throws ErrorException solve_milp(net4, V4, 0.01; warm_start=:nonsense)
    end

    @testset "max_reduction floors the network size" begin
        sol = solve_milp(net4, V4, 1.0; hops=3, max_reduction=0.25)
        @test length(sol.kept) >= 3
    end

    @testset "the operating point is the caller's to judge" begin
        # `powerflow_residual` measures how far a V is from solving its own
        # network, which is what bounds how fine a budget that V can certify.
        # It is reported, not enforced: solve_milp accepts whatever it is given.
        net = read_network_ybus(joinpath(@__DIR__, "..", "data", "ieee34"))
        V = read_voltage(joinpath(@__DIR__, "..", "data", "ieee34"), net)
        @test powerflow_residual(net, V) < 1e-9

        skewed = V .* 1.02                       # plausible-looking, and wrong
        @test powerflow_residual(net, skewed) > 0.01
        @test solve_milp(net, skewed, 0.01; hops=3) isa MilpSolution
        @test solve_milp(net, skewed, 0.1; hops=3) isa MilpSolution

        # A well-conditioned inverse is not the same question: ieee123 inverts
        # poorly (cond ~7e14) yet its operating point is sound.
        poor = read_network_ybus(joinpath(@__DIR__, "..", "data", "ieee123"))
        Vp = read_voltage(joinpath(@__DIR__, "..", "data", "ieee123"), poor)
        @test OptiKRON.bus_impedance(poor) isa Matrix
        @test powerflow_residual(poor, Vp) < 1e-6
    end

    @testset "bad input is rejected" begin
        @test_throws ErrorException solve_milp(net4, V4, 0.01; approach=:Q)
        @test_throws ErrorException solve_milp(net4, V4, -0.01)
        @test_throws ErrorException solve_milp(net4, V4, 0.01; max_reduction=2.0)
        @test_throws ErrorException solve_milp(net4, V4[:, 1:1], 0.01)
        @test_throws ErrorException solve_milp(net4, V4, 0.01; scenarios=Int[])
        @test_throws ErrorException solve_milp(net4, V4, 0.01; scenarios=[99])
        @test_throws ErrorException solve_milp(net4, V4, 0.01; Z=zeros(ComplexF64, 2, 2))
    end
end

@testset "OptiKRON.pipeline" begin
    result = optikron(joinpath(@__DIR__, "..", "data", "test3ph"); Ē=0.02, hops=3)

    @test result isa Reduction
    @test result.Ē == 0.02
    @test result.scenarios == [1]
    @test size(result.V) == (nphase_rows(result.network), 1)
    @test is_radial(result.network, result.assignment)
    @test enforced_violation(result) <= 0
    @test isnan(held_out_violation(result))          # every scenario was enforced

    # The report is what a user reads; it must at least render.
    @test occursin("Reduction", sprint(show, MIME"text/plain"(), result))

    @test optikron(result.network; Ē=0.02, hops=3, radiality=:none).reinserted == Int[]
    @test_throws ErrorException optikron(result.network; radiality=:sideways)
    @test_throws ErrorException optikron(result.network; backend=:quantum)
    @test_throws ErrorException load_case("no-such-case")
end

@testset "OptiKRON.matpower export" begin
    net = read_network_csv(joinpath(@__DIR__, "..", "data", "case69"))
    V = powerflow(net)
    sol = solve_milp(net, V, 0.001; hops=5)
    A, _ = radialize(net, sol.A)
    kept = sort(super_nodes(A))
    path = joinpath(mktempdir(), "case69_reduced.m")

    write_matpower(path, net, A)
    @test occursin("Bus map", read(path, String))     # carries its own mapping back

    # Read it back with our own MATPOWER reader: the export is only useful if a
    # MATPOWER-format parser can consume it.
    mpc = read_matpower(path)
    @test length(mpc.bus) == length(kept)
    @test length(mpc.gen) == 1

    # Reassembling Ybus from the exported branches and shunts must reproduce the
    # Kron-reduced matrix the export came from -- that equality is what makes a
    # power flow on the exported case the same computation as on the reduction.
    Y_red, _ = kron_reduce(net, kept)
    position = Dict(Int(row[BUS_I]) => k for (k, row) in enumerate(mpc.bus))
    rebuilt = zeros(ComplexF64, length(kept), length(kept))
    for row in mpc.branch
        f, t = position[Int(row[F_BUS])], position[Int(row[T_BUS])]
        y = 1 / (row[BR_R] + im * row[BR_X])
        rebuilt[f, f] += y; rebuilt[t, t] += y
        rebuilt[f, t] -= y; rebuilt[t, f] -= y
    end
    for (k, row) in enumerate(mpc.bus)
        rebuilt[k, k] += (row[GS] + im * row[BS]) / mpc.base_mva
    end
    # Relative to the matrix scale: case69's admittances reach 3e4, so what is
    # left here is rounding at machine epsilon, not lost information.
    @test maximum(abs.(rebuilt .- Y_red)) < 1e-15 * maximum(abs, Y_red)

    # Loads: generation-positive here, demand-positive in MATPOWER.
    aggregated = aggregate_injections(A, net)[:, 1]
    @test sum(row[PD] for row in mpc.bus) / mpc.base_mva ≈ -sum(real, aggregated) atol = 1e-8
    @test sum(row[QD] for row in mpc.bus) / mpc.base_mva ≈ -sum(imag, aggregated) atol = 1e-8

    # Exactly one slack, and it is the feeder's.
    @test count(row -> Int(row[BUS_TYPE]) == 3, mpc.bus) == 1
    @test Int(mpc.bus[findfirst(r -> Int(r[BUS_TYPE]) == 3, mpc.bus)][BUS_I]) ==
          parse(Int, net.bus_ids[net.slack])

    @test_throws ErrorException write_matpower(path, net, A; scenario=99)
    three_phase = read_network_csv(joinpath(@__DIR__, "..", "data", "test3ph"))
    @test_throws ErrorException write_matpower(path, three_phase,
        identity_assignment(three_phase))
end

@testset "OptiKRON.reduced CSV export" begin
    # The three-phase feeders have no MATPOWER form, so their reduced networks go
    # out as CSV. The property that matters is that they come back: a reduced
    # network nobody can load is not an artifact, it is a receipt.
    directory = joinpath(@__DIR__, "..", "data", "ieee37")
    net = read_network_ybus(directory)
    V = read_voltage(directory, net)
    solution = solve_milp(net, V, 0.01; hops=5, enforce_radiality=true,
        direction=:downstream)

    out = mktempdir()
    files = write_reduced_csv(out, net, solution.A; V=V)

    @test length(files) == 5
    @test all(isfile, files)
    for name in ("bus.csv", "ybus.csv", "load.csv", "voltage.csv", "assignment.csv")
        @test isfile(joinpath(out, name))
    end

    # Round-trips through the reader the input format uses.
    reloaded = load_case(out)
    kept = sort(solution.kept)
    @test nnodes(reloaded) == length(kept)
    @test reloaded.bus_ids == net.bus_ids[kept]
    @test nscenarios(reloaded) == nscenarios(net)

    # And the admittance it carries is the Schur complement, to full precision.
    Y_reduced, kept_rows = kron_reduce(net, kept)
    @test Matrix(reloaded.Ybus) ≈ Y_reduced rtol = 1e-15
    @test nphase_rows(reloaded) == length(kept_rows)

    # Injections are the aggregated ones; eliminated buses contributed their load.
    aggregated = aggregate_injections(solution.A, net)[kept_rows, :]
    @test Matrix(reloaded.S) ≈ aggregated rtol = 1e-14
    @test sum(reloaded.S) ≈ sum(net.S) rtol = 1e-12

    # The slack survives as the slack.
    @test reloaded.bus_ids[reloaded.slack] == net.bus_ids[net.slack]

    # Every original bus appears in the map exactly once.
    map_lines = readlines(joinpath(out, "assignment.csv"))
    @test length(map_lines) == nnodes(net) + 1
    @test count(l -> endswith(l, ",1"), map_lines) == length(kept)

    @test_throws ErrorException write_reduced_csv(out, net, solution.A; scenarios=[99])
end
