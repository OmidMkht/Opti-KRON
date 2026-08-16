using JuMP

# Exported functions with no prior direct coverage -- each is exercised as a
# dependency of solve_milp/optikron elsewhere in the suite, but never called
# under its own name. This pins the contract each one documents.
@testset "OptiKRON.api" begin

    net4 = read_network_csv(TEST4)
    V4 = powerflow(net4)
    net3 = read_network_csv(TEST3PH)
    V3 = powerflow(net3)

    @testset "case_directory" begin
        dir = case_directory("test4")
        @test isdir(dir)
        @test isfile(joinpath(dir, "bus.csv"))
        @test_throws ErrorException case_directory("no-such-case")
    end

    @testset "list_cases" begin
        io = IOBuffer()
        list_cases(io; root=joinpath(@__DIR__, "..", "data"))
        out = String(take!(io))
        @test occursin("ieee37", out)
        @test !occursin("test4", out)      # unit-test fixtures are skipped
    end

    @testset "phases_of" begin
        # test3ph mixes abc/a -- N1 carries all three, N3 only phase a.
        @test phases_of(net3, 1) == [:a, :b, :c]
        @test phases_of(net3, 3) == [:a]
    end

    @testset "add_radiality_constraints!" begin
        tree = orient_radial(net4)
        factory, _ = select_optimizer(prefer=:highs)
        model = Model(factory)
        @variable(model, Akk[1:nnodes(net4)], Bin)
        n = add_radiality_constraints!(model, Akk, tree)
        # test4 is a star on bus 2 (1 <- 2 -> 3,4): exactly one branching node.
        @test n == 1
        @test n == count(j -> length(tree.children[j]) >= 2, eachindex(tree.children))
    end

    @testset "screen_voltage_spread! kills what the budget cannot reach" begin
        tree = orient_radial(net4)
        admissible = admissible_pairs(net4, tree; hops=3)
        tight = copy(admissible)
        @test screen_voltage_spread!(tight, net4, tree, abs.(V4), 1e-9, [1]) > 0
        loose = copy(admissible)
        @test screen_voltage_spread!(loose, net4, tree, abs.(V4), 1.0, [1]) == 0
    end

    @testset "screen! reports a consistent accounting" begin
        tree = orient_radial(net3)
        admissible = admissible_pairs(net3, tree; hops=3)
        Z = bus_impedance(net3)
        C = conj.(net3.S ./ V3)
        needed, report = screen!(copy(admissible), net3, tree, V3, 0.005, Z, C, [1])

        @test report.pairs_after <= report.pairs_before
        @test report.rows_kept == sum(length, values(needed); init=0)
        @test report.rows_kept + report.rows_skipped >= report.rows_kept
        @test all(>=(0), report.deviation)
        @test size(report.deviation) == (nphase_rows(net3), 1)
    end

    @testset "export_matpower and export_reduced dispatch on feeder kind" begin
        # Single-phase, branch-based: export_reduced must pick MATPOWER.
        net69 = read_network_csv(joinpath(@__DIR__, "..", "data", "case69"))
        V69 = powerflow(net69)
        sol = solve_milp(net69, V69, 0.001; hops=5)
        A, _ = radialize(net69, sol.A)
        r = Reduction(net69, V69, sol, A, Int[], 0.001, [1], :post, Device[], nothing)

        out1 = mktempdir()
        files1 = export_matpower(r, out1)
        @test length(files1) == 1
        @test endswith(only(files1), ".m")
        @test export_reduced(r, out1) == files1

        # Three-phase: export_reduced must fall back to CSV, since a MATPOWER
        # branch cannot carry phase a and c but not b.
        sol3 = solve_milp(net3, V3, 0.02; hops=3)
        r3 = Reduction(net3, V3, sol3, sol3.A, Int[], 0.02, [1], :none, Device[], nothing)
        out3 = mktempdir()
        files3 = export_reduced(r3, out3)
        @test any(endswith(f, "assignment.csv") for f in files3)
    end
end
