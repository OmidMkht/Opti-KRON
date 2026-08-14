const IEEE37 = joinpath(@__DIR__, "..", "data", "ieee37")
const IEEE123 = joinpath(@__DIR__, "..", "data", "ieee123")

"Rows of a Kron-reduced matrix belonging to a kept bus."
function reduced_block(net, kept_sorted, bus)
    position = findfirst(==(bus), kept_sorted)
    position === nothing && return nothing
    lengths = [count(view(net.phases, :, b)) for b in kept_sorted]
    offset = sum(lengths[1:position-1])
    return (offset+1):(offset+lengths[position])
end

@testset "OptiKRON.devices" begin

    @testset "reading device tables" begin
        net = read_network_ybus(IEEE37)
        devices = read_devices(IEEE37, net)

        # ieee37 ships two transformers, both explicit_ybus, no switches.
        @test length(devices) == 2
        @test all(d -> d.kind === :transformer, devices)
        @test all(d -> d.note == "explicit_ybus", devices)
        @test preserved_buses(devices) == sort(unique(vcat([d.from for d in devices],
            [d.to for d in devices])))
        @test length(preserved_buses(devices)) == 4

        # Terminals resolve to real buses, and a device is never a self-loop.
        for d in devices
            @test 1 <= d.from <= nnodes(net)
            @test 1 <= d.to <= nnodes(net)
            @test d.from != d.to
        end
    end

    @testset "switches are separable from transformers" begin
        net = read_network_ybus(IEEE123)
        all_devices = read_devices(IEEE123, net)
        transformers = read_devices(IEEE123, net; switches=false)

        @test length(transformers) == 1
        @test length(all_devices) == 7                       # 1 transformer + 6 closed switches
        @test count(d -> d.kind === :switch, all_devices) == 6
        @test isempty(read_devices(IEEE123, net; transformers=false, switches=false))

        # Preserving switches costs strictly more buses than not.
        @test length(preserved_buses(all_devices)) > length(preserved_buses(transformers))
    end

    @testset "collapsed devices are not preservable" begin
        # ieee8500 names 1178 transformers but all except one were collapsed into
        # their primary at export, so their secondary is not a bus here. Pinning
        # on the name alone would freeze 45% of the feeder for nothing.
        net = read_network_ybus(joinpath(@__DIR__, "..", "data", "ieee8500"))
        devices = read_devices(joinpath(@__DIR__, "..", "data", "ieee8500"), net;
            switches=false)

        @test length(devices) == 1
        @test devices[1].note == "explicit_ybus"
        @test length(preserved_buses(devices)) == 2
    end

    @testset "_apply_pins! clears only the pinned column" begin
        net = read_network_ybus(IEEE37)
        tree = orient_radial(net)
        admissible = admissible_pairs(net, tree; hops=5)
        before = copy(admissible)

        OptiKRON._apply_pins!(admissible, [7], nnodes(net))

        @test admissible[7, 7]                                # still represents itself
        @test all(!admissible[k, 7] for k in 1:nnodes(net) if k != 7)
        # Bus 7 may still represent others, and nothing else changed.
        @test all(admissible[7, j] == before[7, j] for j in 1:nnodes(net))
        @test all(admissible[i, j] == before[i, j]
                  for i in 1:nnodes(net), j in 1:nnodes(net) if j != 7)

        @test_throws ErrorException OptiKRON._apply_pins!(admissible, [0], nnodes(net))
        @test_throws ErrorException OptiKRON._apply_pins!(admissible, [nnodes(net) + 1],
            nnodes(net))
    end

    @testset "pinned buses survive the reduction" begin
        net = read_network_ybus(IEEE37)
        V = read_voltage(IEEE37, net)
        pin = preserved_buses(read_devices(IEEE37, net))

        solution = solve_milp(net, V, 0.01; hops=5, pin=pin, enforce_radiality=true,
            direction=:downstream)

        @test issubset(pin, solution.kept)
        @test all(solution.A[i, i] ≈ 1.0 for i in pin)
        # Pinning costs reduction but must not break the budget.
        @test annulus_violation(net, solution.A, V, 0.01) <= 0
    end

    @testset "preserved devices come through Kron reduction exactly" begin
        # The claim the whole feature rests on: with both terminals kept, the
        # Schur complement leaves the device's admittance block untouched --
        # not approximately, exactly. Y_rr is block-diagonal across the two
        # sides of the device edge, so the correction term at (i,j) is zero.
        for directory in (IEEE37, IEEE123)
            net = read_network_ybus(directory)
            V = read_voltage(directory, net)
            devices = read_devices(directory, net)

            solution = solve_milp(net, V, 0.01; hops=5, pin=preserved_buses(devices),
                enforce_radiality=true, direction=:downstream)
            kept = super_nodes(solution.A)
            Y_reduced, _ = kron_reduce(net, kept)
            kept_sorted = sort(kept)
            blocks = node_rows(net)

            @test !isempty(devices)
            for d in devices
                rows = reduced_block(net, kept_sorted, d.from)
                cols = reduced_block(net, kept_sorted, d.to)
                @test rows !== nothing && cols !== nothing

                original = Matrix(net.Ybus[blocks[d.from], blocks[d.to]])
                @test Y_reduced[rows, cols] == original       # bit-for-bit, not ≈
            end
        end
    end

    @testset "preserve option in the pipeline" begin
        loose = optikron("ieee123"; Ē=0.01, hops=5, preserve=:none, time_limit=600)
        transformers = optikron("ieee123"; Ē=0.01, hops=5, preserve=:transformers,
            time_limit=600)
        devices = optikron("ieee123"; Ē=0.01, hops=5, preserve=:devices, time_limit=600)

        @test isempty(loose.devices)
        @test length(transformers.devices) == 1
        @test length(devices.devices) == 7

        # More preserved equipment means more buses kept, never fewer.
        @test length(loose.solution.kept) <= length(transformers.solution.kept)
        @test length(transformers.solution.kept) <= length(devices.solution.kept)

        # And every one of them stays inside the budget.
        for r in (loose, transformers, devices)
            @test enforced_violation(r) <= 0
        end

        @test_throws ErrorException optikron("ieee123"; preserve=:everything)
    end

    @testset "a case with no device tables pins nothing" begin
        net = read_network_csv(TEST4)
        @test isempty(read_devices(TEST4, net))
        @test isempty(preserved_buses(read_devices(TEST4, net)))
    end
end
