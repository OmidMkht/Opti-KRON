using Test
using LinearAlgebra, SparseArrays

include(joinpath(@__DIR__, "..", "src", "OptiKRON.jl"))
include(joinpath(@__DIR__, "..", "tools", "convert_matpower.jl"))
using .OptiKRON

const TEST4 = joinpath(@__DIR__, "..", "data", "test4")

@testset "OptiKRON.io" begin

    @testset "phase parsing" begin
        @test OptiKRON.parse_phases("abc") == [:a, :b, :c]
        @test OptiKRON.parse_phases("CA") == [:a, :c]      # normalised and sorted
        @test OptiKRON.parse_phases(" b ") == [:b]
        @test_throws ErrorException OptiKRON.parse_phases("abd")
        @test_throws ErrorException OptiKRON.parse_phases("aa")
        @test_throws ErrorException OptiKRON.parse_phases("")
    end

    @testset "branch-based CSV import" begin
        net = read_network_csv(TEST4; name="test4")

        @test nnodes(net) == 4
        @test nphase_rows(net) == 4
        @test nscenarios(net) == 2
        @test !is_three_phase(net)
        @test has_branch_data(net)
        @test net.slack == 1

        # 1 / (0.01 + 0.02im) == 20 - 40im
        @test net.Ybus[1, 1] ≈ 20.0 - 40.0im
        @test net.Ybus[1, 2] ≈ -20.0 + 40.0im
        # Every row of a grounded-at-slack ladder sums to zero except the slack.
        @test all(abs(sum(net.Ybus[i, :])) < 1e-9 for i in 2:4)

        # Radial star on bus 2, symmetric, hollow.
        @test Matrix(net.Lambda) == [0 1 0 0; 1 0 1 1; 0 1 0 0; 0 1 0 0]

        # Slack carries no injection; loads are negative (injection convention).
        @test net.S[1, 1] == 0
        @test net.S[2, 1] ≈ -0.05 - 0.02im
        @test net.bus_ids == ["N1", "N2", "N3", "N4"]
    end

    @testset "node_rows layout" begin
        # Two buses: abc then ac -> rows 1:3 and 4:5.
        phases = Bool[1 1; 1 0; 1 1]
        @test node_rows(phases) == [[1, 2, 3], [4, 5]]
    end

    @testset "matrices constructor" begin
        csv_net = read_network_csv(TEST4)
        net = network_from_matrices(csv_net.Ybus, csv_net.S; slack=1, name="from-matrices")

        @test nnodes(net) == 4
        @test !has_branch_data(net)                     # Ybus path drops line params
        @test Matrix(net.Lambda) == Matrix(csv_net.Lambda)
        @test net.Ybus ≈ csv_net.Ybus
    end

    @testset "validation catches bad input" begin
        csv_net = read_network_csv(TEST4)

        # An ungrounded Laplacian: every row sums to zero, so Ybus is singular.
        # Catch it at import, not inside a factorization.
        Y_float = Matrix(csv_net.Ybus)
        Y_float[1, :] .= 0.0 + 0.0im
        @test_throws ErrorException network_from_matrices(sparse(Y_float), csv_net.S; slack=1)

        # Slack index out of range.
        @test_throws ErrorException network_from_matrices(csv_net.Ybus, csv_net.S; slack=99)

        # A disconnected feeder cannot be reduced coherently.
        Y_split = Matrix(csv_net.Ybus)
        Y_split[2, 4] = Y_split[4, 2] = 0.0 + 0.0im
        Y_split[4, 4] = 1.0 - 1.0im          # keep the row nonzero so we hit the connectivity check
        @test_throws ErrorException network_from_matrices(sparse(Y_split), csv_net.S; slack=1)

        # Injection matrix whose row count disagrees with the phase mask.
        @test_throws ErrorException network_from_matrices(csv_net.Ybus, csv_net.S[1:3, :]; slack=1)
    end
end

include("test_core.jl")
include("test_powerflow.jl")
include("test_solver.jl")
include("test_milp.jl")
include("test_devices.jl")
include("test_contraction.jl")
include("test_search.jl")
include("test_api.jl")
