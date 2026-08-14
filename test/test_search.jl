@testset "OptiKRON.search" begin
    directory = joinpath(@__DIR__, "..", "data", "ieee37")
    net = read_network_ybus(directory)
    V = read_voltage(directory, net)

    @testset "holds the budget it is given" begin
        for Ē in (0.002, 0.01)
            s = search_reduce(net, V, Ē)
            @test s isa SearchSolution
            @test s isa ReductionSolution
            @test s.objective <= Ē
            # The independent nonconvex check, not the one the search used.
            @test annulus_violation(net, s.A, V, Ē) <= 0
            @test s.reduction ≈ 1 - length(s.kept) / nnodes(net)
            @test net.slack in s.kept
        end
    end

    @testset "the assignment is well formed" begin
        s = search_reduce(net, V, 0.01)
        # Every bus represented exactly once, and only super-nodes represent.
        @test all(sum(s.A, dims=1) .≈ 1)
        @test sort(s.kept) == sort(super_nodes(s.A))
        for j in axes(s.A, 2), i in axes(s.A, 1)
            iszero(s.A[i, j]) && continue
            @test s.A[i, i] ≈ 1
            # φ_j ⊆ φ_i: a super-node carries every phase it stands for.
            @test all(p -> !net.phases[p, j] || net.phases[p, i], 1:3)
        end
    end

    @testset "a looser budget never reduces less" begin
        tight = search_reduce(net, V, 0.001)
        loose = search_reduce(net, V, 0.01)
        @test length(loose.kept) <= length(tight.kept)
    end

    @testset "it does not beat the optimum" begin
        # The greedy is a heuristic. On the same feasible set it can match the
        # MILP but must never do better -- if it does, the two are not measuring
        # the same error, which is a bug and not a result.
        Ē = 0.01
        s = search_reduce(net, V, Ē; radial=true, upstream=true)
        m = solve_milp(net, V, Ē; hops=nnodes(net), enforce_radiality=true,
            direction=:downstream)
        @test length(s.kept) >= length(m.kept)
    end

    @testset "pinned buses survive" begin
        pin = preserved_buses(read_devices(directory, net))
        s = search_reduce(net, V, 0.01; pin=pin)
        @test issubset(pin, s.kept)
        @test annulus_violation(net, s.A, V, 0.01) <= 0
        @test_throws ErrorException search_reduce(net, V, 0.01; pin=[0])
    end

    @testset "options are checked" begin
        @test_throws ErrorException search_reduce(net, V, 0.01; objective=:median)
        @test_throws ErrorException search_reduce(net, V, 0.01; backend=:tpu)
        @test_throws ErrorException search_reduce(net, V, -1.0)
        @test_throws ErrorException search_reduce(net, V, 0.01; scenarios=Int[])
        @test_throws ErrorException search_reduce(net, V, 0.01; scenarios=[99])
    end

    @testset "max_reduction is respected" begin
        s = search_reduce(net, V, 0.01; max_reduction=0.25)
        @test length(s.kept) >= ceil(Int, nnodes(net) * 0.75)
        @test s.status === :max_reduction
    end

    @testset "screening changes nothing but the work" begin
        # Screening is rigorous: it only drops pairs that provably cannot meet
        # the bound, so removing it must not change the answer.
        on = search_reduce(net, V, 0.01; screen=true)
        off = search_reduce(net, V, 0.01; screen=false)
        @test on.A == off.A
        @test on.screened > 0
        @test off.screened == 0
    end

    @testset "host and device agree" begin
        # Same formulas either way, so a disagreement here is a real defect
        # rather than a second implementation drifting.
        host = search_reduce(net, V, 0.01; backend=:zbus)
        @test host.A == search_reduce(net, V, 0.01; backend=:cpu).A
        if gpu_available()
            device = search_reduce(net, V, 0.01; backend=:gpu)
            @test device.A == host.A
            @test device.objective ≈ host.objective atol = 1e-12
            @test device.solver === :search_gpu
        else
            # No device is a normal state: warn and run on the CPU, never fail.
            fallback = @test_logs (:warn,) search_reduce(net, V, 0.01; backend=:gpu)
            @test fallback.A == host.A
            @test fallback.solver === :search_cpu
        end
    end

    @testset "through the pipeline" begin
        r = optikron("ieee37"; Ē=0.01, backend=:search_cpu, time_limit=120)
        @test r.solution isa SearchSolution
        @test enforced_violation(r) <= 0
        @test is_radial(r.network, r.assignment)
    end
end
