const IEEE37 = joinpath(@__DIR__, "..", "data", "ieee37")
const IEEE34 = joinpath(@__DIR__, "..", "data", "ieee34")
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

    @testset "kinds come from the classifying tables" begin
        net = read_network_ybus(IEEE34)
        devices = read_devices(IEEE34, net)

        # ieee34 ships 8 transformer rows: 6 are named by regulator.csv, 1 by
        # phase_shift_equipment.csv, leaving 1 plain transformer.
        counts = Dict(k => count(d -> d.kind === k, devices) for k in OptiKRON.DEVICE_KINDS)
        @test counts[:regulator] == 6
        @test counts[:phase_shift] == 1
        @test counts[:transformer] == 1
        @test counts[:switch] == 0
        @test length(devices) == 8

        # The phase shifter is the delta-wye substation transformer.
        shifter = only(filter(d -> d.kind === :phase_shift, devices))
        @test shifter.id == "subxf"
        @test shifter.note == "delta-wye"

        for d in devices
            @test 1 <= d.from <= nnodes(net)
            @test 1 <= d.to <= nnodes(net)
            @test d.from != d.to
        end
    end

    @testset "kinds select independently" begin
        net = read_network_ybus(IEEE123)

        @test isempty(read_devices(IEEE123, net; kinds=()))
        @test all(d -> d.kind === :switch, read_devices(IEEE123, net; kinds=(:switch,)))
        @test length(read_devices(IEEE123, net; kinds=(:regulator,))) == 7

        # `:all` adds the shunt kind on top of the edge kinds the default reads.
        # ieee123 ships 4 capacitor banks.
        @test length(read_devices(IEEE123, net; kinds=:all)) ==
              length(read_devices(IEEE123, net)) + 4
        @test all(d -> d.from == d.to,
            read_devices(IEEE123, net; kinds=(:capacitor,)))

        # Selecting a subset never returns more than selecting everything.
        subset = read_devices(IEEE123, net; kinds=(:regulator, :switch))
        @test length(preserved_buses(subset)) <=
              length(preserved_buses(read_devices(IEEE123, net)))

        @test_throws ErrorException read_devices(IEEE123, net; kinds=(:inverter,))
    end

    @testset "open switches are not edges" begin
        net = read_network_ybus(IEEE123)
        switches = read_devices(IEEE123, net; kinds=(:switch,))

        # switch.csv lists 8, two of them open; an open switch is not in Ybus.
        @test length(switches) == 6
        @test all(d -> d.note == "closed", switches)
    end

    @testset "service transformers dominate a large feeder" begin
        # ieee8500 ships 2354 service transformers against 12 regulators and 38
        # switches. Pinning them all is a hard ceiling on reduction, which is why
        # `preserve` is per-kind rather than all-or-nothing.
        dir = joinpath(@__DIR__, "..", "data", "ieee8500")
        net = read_network_ybus(dir)

        everything = preserved_buses(read_devices(dir, net))
        stateful = preserved_buses(read_devices(dir, net;
            kinds=(:regulator, :phase_shift, :switch)))

        @test length(everything) / nnodes(net) > 0.4
        @test length(stateful) / nnodes(net) < 0.05
        @test issubset(stateful, everything)
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
        # The claim the whole feature rests on: with both terminals kept, Y_rr is
        # block-diagonal across the two sides of the device edge, so the Schur
        # correction at (i,j) is identically zero and the block is untouched.
        #
        # Asserted relative to the block's own scale rather than with `==`. The
        # correction is zero in exact arithmetic, but `pinv` reaches it through an
        # SVD, which does not preserve block-diagonality bit-for-bit; the residue
        # runs to 1.5e-14 absolute on a jumper with entries near 1.2e6, which is
        # under one ULP and as exact as Float64 gets.
        for directory in (IEEE37, IEEE34, IEEE123)
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
                scale = maximum(abs, original)
                @test maximum(abs, Y_reduced[rows, cols] .- original) <= 8 * eps(scale)
            end
        end
    end

    @testset "preserve option in the pipeline" begin
        bare = optikron("ieee123"; Ē=0.01, hops=5,
            preserve=(:center_tap, :phase_shift), time_limit=600)
        sw = optikron("ieee123"; Ē=0.01, hops=5, preserve=:switches, time_limit=600)
        reg = optikron("ieee123"; Ē=0.01, hops=5, preserve=:regulators, time_limit=600)
        required = optikron("ieee123"; Ē=0.01, hops=5, preserve=:required, time_limit=600)
        everything = optikron("ieee123"; Ē=0.01, hops=5, preserve=:all, time_limit=600)

        # Plain transformers and capacitor banks are opt-in, so only `:all` has them.
        for r in (bare, sw, reg, required)
            @test all(d -> d.kind !== :transformer, r.devices)
            @test all(d -> d.kind !== :capacitor, r.devices)
        end
        @test any(d -> d.kind === :transformer, everything.devices)
        @test any(d -> d.kind === :capacitor, everything.devices)

        @test all(d -> d.kind !== :regulator, sw.devices)
        @test all(d -> d.kind !== :switch, reg.devices)
        @test length(required.devices) == length(sw.devices) + length(reg.devices)

        # Each single-kind option preserves at least as much as the bare set,
        # `:required` at least as much as either, and `:all` most of all.
        for r in (sw, reg)
            @test length(bare.solution.kept) <= length(r.solution.kept)
            @test length(r.solution.kept) <= length(required.solution.kept)
        end
        @test length(required.solution.kept) <= length(everything.solution.kept)

        # And every one of them stays inside the budget.
        for r in (bare, sw, reg, required, everything)
            @test enforced_violation(r) <= 0
        end

        @test_throws ErrorException optikron("ieee123"; preserve=:everything)
    end

    @testset "_preserve_kinds resolves the shorthands" begin
        @test OptiKRON._preserve_kinds(:required) ==
              (:center_tap, :phase_shift, :regulator, :switch)
        @test OptiKRON._preserve_kinds(:all) == OptiKRON.DEVICE_KINDS
        @test OptiKRON._preserve_kinds(:switches) == (:center_tap, :phase_shift, :switch)
        @test OptiKRON._preserve_kinds(:regulators) == (:center_tap, :phase_shift, :regulator)
        @test OptiKRON._preserve_kinds((:switch,)) == (:switch,)

        # Neither a phase shift nor a center-tapped transformer can be recovered
        # from a Schur complement, so both stay pinned whatever the caller asks
        # for on the switch/regulator axis. Plain transformers and capacitor banks
        # cost reduction and are opt-in.
        for option in (:required, :switches, :regulators)
            kinds = OptiKRON._preserve_kinds(option)
            @test :phase_shift in kinds
            @test :center_tap in kinds
            @test :transformer ∉ kinds
            @test :capacitor ∉ kinds
        end
    end

    @testset "a case with no device tables pins nothing" begin
        net = read_network_csv(TEST4)
        @test isempty(read_devices(TEST4, net))
        @test isempty(preserved_buses(read_devices(TEST4, net)))
    end
end
