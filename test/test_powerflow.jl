@testset "OptiKRON.powerflow" begin

    net4 = read_network_csv(TEST4)
    net3 = read_network_csv(joinpath(@__DIR__, "..", "data", "test3ph"))

    @testset "balanced slack follows the phase, not the position" begin
        # A three-phase source is 1<0, 1<-120, 1<120.
        v = slack_voltage(net3)
        @test length(v) == 3
        @test all(abs.(v) .≈ 1.0)
        @test rad2deg(angle(v[1])) ≈ 0 atol = 1e-9
        @test rad2deg(angle(v[2])) ≈ -120 atol = 1e-9
        @test rad2deg(angle(v[3])) ≈ 120 atol = 1e-9
        @test sum(v) ≈ 0 atol = 1e-12          # balanced: the three sum to zero

        # A slack missing phase b is 1<0 and 1<120 -- not 1<0 and 1<-120, which
        # is what indexing by position rather than by phase would give.
        phases = Bool[1 1; 0 0; 1 1]
        partial = network_from_matrices(net3.Ybus[[1, 3, 4, 6], [1, 3, 4, 6]],
            net3.S[[1, 3, 4, 6], :]; phases=phases, slack=1)
        vp = slack_voltage(partial)
        @test length(vp) == 2
        @test rad2deg(angle(vp[2])) ≈ 120 atol = 1e-9

        @test all(abs.(slack_voltage(net3; magnitude=1.05)) .≈ 1.05)
    end

    @testset "solves the shipped cases" begin
        for net in (net4, net3)
            V = powerflow(net)
            @test size(V) == (nphase_rows(net), nscenarios(net))
            @test powerflow_residual(net, V) < 1e-9

            # The slack is held exactly where it was put.
            @test V[node_rows(net)[net.slack], :] ≈
                  repeat(slack_voltage(net), 1, nscenarios(net))
        end
    end

    @testset "the residual rejects a wrong operating point" begin
        net = net3
        V = powerflow(net)
        good = powerflow_residual(net, V)

        # A flat start ignores the 120-degree spacing; a rescaled V is in the
        # wrong per-unit base. Both must score far worse than the solution.
        @test powerflow_residual(net, fill(1.0 + 0.0im, size(V))) > 100 * good
        @test powerflow_residual(net, V .* sqrt(3)) > 100 * good

        # Rotating the feeder *away from* its slack is inconsistent. Rotating
        # everything including the slack is not -- power flow is invariant to the
        # angle reference, so that one must still pass.
        skewed = copy(V)
        other = setdiff(1:nphase_rows(net), node_rows(net)[net.slack])
        skewed[other, :] .*= cis(0.05)
        @test powerflow_residual(net, skewed) > 100 * good
        @test powerflow_residual(net, V .* cis(0.05)) ≈ good atol = 1e-9

        @test_throws ErrorException powerflow_residual(net4, powerflow(net4)[:, 1:1])
    end

    @testset "a raised slack lifts the whole feeder" begin
        V_nominal = powerflow(net3)
        V_raised = powerflow(net3; vslack=slack_voltage(net3; magnitude=1.05))
        @test all(abs.(V_raised) .> abs.(V_nominal))
        @test powerflow_residual(net3, V_raised) < 1e-9
    end

    @testset "bad input is rejected" begin
        # One value per slack phase, and no fewer.
        @test_throws ErrorException powerflow(net3; vslack=[1.0 + 0im])
        @test_throws ErrorException powerflow(net3; vslack=ComplexF64[])

        # An unreachable operating point must fail loudly rather than return
        # a half-converged answer that would invalidate every downstream bound.
        overloaded = network_from_matrices(net3.Ybus, net3.S .* 1e4; slack=net3.slack,
            phases=net3.phases)
        @test_throws ErrorException powerflow(overloaded)
    end
end
