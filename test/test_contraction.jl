const C123 = joinpath(@__DIR__, "..", "data", "ieee123")
const C34 = joinpath(@__DIR__, "..", "data", "ieee34")

@testset "OptiKRON.contraction" begin

    net123 = read_network_ybus(C123)
    switches123 = read_devices(C123, net123; kinds=(:switch,))

    @testset "the index map is total and consistent" begin
        contracted, c = contract_switches(net123, switches123)

        @test length(c.parent) == nnodes(net123) == c.nbuses
        @test ncontracted(c) == nnodes(contracted)
        @test sort(vcat(c.members...)) == 1:nnodes(net123)   # a partition, exactly
        @test all(j in c.members[c.parent[j]] for j in 1:nnodes(net123))
        @test all(c.parent[j] == m for m in 1:ncontracted(c) for j in c.members[m])

        # Groups are ordered by their smallest member, so numbering is stable.
        @test issorted([first(g) for g in c.members])
        @test nnodes(contracted) == nnodes(net123) - length(switches123)
        @test count(g -> length(g) > 1, c.members) == length(switches123)

        # Both terminals of every switch land in one group.
        @test all(c.parent[d.from] == c.parent[d.to] for d in switches123)
    end

    @testset "contraction conserves the network" begin
        contracted, c = contract_switches(net123, switches123)

        # Injections are summed, never lost.
        @test sum(contracted.S) ≈ sum(net123.S)
        @test nscenarios(contracted) == nscenarios(net123)
        # A group's phases are the union of its members'.
        for m in 1:ncontracted(c), p in 1:3
            @test contracted.phases[p, m] == any(net123.phases[p, j] for j in c.members[m])
        end
        # The slack follows its own group, and the feeder stays radial.
        @test contracted.slack == c.parent[net123.slack]
        @test orient_radial(contracted) isa RadialTree
    end

    @testset "contraction is what fixes the conditioning" begin
        contracted, _ = contract_switches(net123, switches123)
        blocks(n) = setdiff(1:nphase_rows(n), node_rows(n)[n.slack])
        before = svdvals(Matrix{ComplexF64}(net123.Ybus[blocks(net123), blocks(net123)]))
        after = svdvals(Matrix{ComplexF64}(contracted.Ybus[blocks(contracted), blocks(contracted)]))

        # sigma_max is set by the switch jumpers; removing them drops it sharply,
        # while the smallest mode is physical and must survive untouched.
        @test after[1] < before[1] / 10
        @test after[end] ≈ before[end] rtol = 1e-6
        @test after[1] / after[end] < before[1] / before[end] / 10
    end

    @testset "a feeder with no switches is untouched" begin
        net34 = read_network_ybus(C34)
        contracted, c = contract_switches(net34, Device[])
        @test c.parent == 1:nnodes(net34)
        @test nnodes(contracted) == nnodes(net34)
        @test all(length(g) == 1 for g in c.members)
    end

    @testset "contract=true is safe wherever switch data is missing" begin
        # The option is driven by what the case ships, not by what the caller
        # promises. Asking for it where there is nothing to contract must be a
        # no-op, never an error -- three different ways of having no switches.
        for case in ("case69",          # MATPOWER: no device tables at all
            "R100",            # CSV: no device tables at all
            "ieee34")          # switch.csv present but empty
            r = optikron(case; Ē=0.05, hops=3, preserve=:none, contract=true,
                time_limit=300)
            @test r.contraction === nothing
            @test size(r.assignment) == (nnodes(r.network), nnodes(r.network))
            @test enforced_violation(r) <= 0
        end

        # A pre-loaded `Network` has no directory to read switches from, so the
        # option has nothing to act on and must still come back clean.
        loose = optikron(read_network_ybus(C34); Ē=0.05, hops=3, contract=true,
            time_limit=300)
        @test loose.contraction === nothing

        # And contracting nothing must give exactly the uncontracted answer.
        plain = optikron("ieee34"; Ē=0.05, hops=3, preserve=:none, time_limit=300)
        merged = optikron("ieee34"; Ē=0.05, hops=3, preserve=:none, contract=true,
            time_limit=300)
        @test merged.assignment == plain.assignment
    end

    @testset "uncontract restores every original bus" begin
        contracted, c = contract_switches(net123, switches123)
        tree = orient_radial(net123)
        A_c = identity_assignment(contracted)
        A = uncontract(A_c, c, tree)

        # Nothing reduced in, nothing reduced out.
        @test size(A) == (nnodes(net123), nnodes(net123))
        @test A ≈ identity_assignment(net123)

        @test_throws ErrorException uncontract(zeros(3, 3), c, tree)
    end

    @testset "an absorbed group sends each bus to its own side" begin
        contracted, c = contract_switches(net123, switches123)
        tree = orient_radial(net123)

        # Absorb one switch group into a neighbour and check where its members go.
        group = findfirst(g -> length(g) > 1, c.members)
        host = findfirst(m -> m != group && !iszero(contracted.Lambda[m, group]),
            1:ncontracted(c))
        @test host !== nothing

        A_c = identity_assignment(contracted)
        assign!(A_c, host, group)
        A = uncontract(A_c, c, tree)

        @test all(sum(A[:, j]) == 1 for j in 1:nnodes(net123))     # represented once
        for j in c.members[group]
            rep = findfirst(!iszero, view(A, :, j))
            @test rep in c.members[host]
            # It went to the nearest member of the host group -- the side it is on.
            @test hop_distance(tree, j, rep) ==
                  minimum(hop_distance(tree, j, g) for g in c.members[host])
        end
    end

    @testset "reducing contracted then expanding respects the original budget" begin
        Ē = 0.01
        plain = optikron("ieee123"; Ē=Ē, hops=5, preserve=:both, time_limit=600)
        merged = optikron("ieee123"; Ē=Ē, hops=5, preserve=:both, contract=true,
            time_limit=600)

        @test merged.contraction !== nothing
        @test plain.contraction === nothing
        # The expanded assignment is full size and honours the budget on the
        # ORIGINAL feeder, which is the only claim that matters.
        @test size(merged.assignment) == (nnodes(merged.network), nnodes(merged.network))
        @test enforced_violation(merged) <= 0
        @test enforced_violation(plain) <= 0

        # Preserving switches through contraction costs no reduction and fewer pins.
        @test length(super_nodes(merged.assignment)) <=
              length(super_nodes(plain.assignment))
    end

    @testset "a surviving switch comes through Kron reduction exactly" begin
        r = optikron("ieee123"; Ē=0.01, hops=5, preserve=:both, contract=true,
            time_limit=600)
        net = r.network
        kept = super_nodes(r.assignment)
        sorted = sort(kept)
        Y_reduced, _ = kron_reduce(net, kept)
        blocks = node_rows(net)
        lengths = [count(view(net.phases, :, b)) for b in sorted]

        checked = 0
        for d in read_devices(C123, net; kinds=(:switch,))
            ia = findfirst(==(d.from), sorted)
            ib = findfirst(==(d.to), sorted)
            # Contraction makes a switch all-or-nothing: both ends or neither.
            @test (ia === nothing) == (ib === nothing)
            (ia === nothing || ib === nothing) && continue
            checked += 1

            rows = (sum(lengths[1:ia-1])+1):sum(lengths[1:ia])
            cols = (sum(lengths[1:ib-1])+1):sum(lengths[1:ib])
            original = Matrix(net.Ybus[blocks[d.from], blocks[d.to]])
            @test maximum(abs, Y_reduced[rows, cols] .- original) <=
                  8 * eps(maximum(abs, original))
        end
        @test checked > 0
    end
end
