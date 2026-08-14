const TEST3PH = joinpath(@__DIR__, "..", "data", "test3ph")

@testset "OptiKRON.core" begin

    net = read_network_csv(TEST4)
    tree = orient_radial(net)

    @testset "radial orientation" begin
        # test4 is a star on N2: N1—N2, N2—N3, N2—N4.
        @test tree.parents == [0, 1, 2, 2]
        @test tree.depth == [0, 1, 2, 2]
        @test sort(tree.children[2]) == [3, 4]
        @test nedges(tree) == 3
        @test tree.slack == 1

        @test path_to_root(tree, 3) == [3, 2, 1]
        @test interior_path(tree, 1, 3) == [2]        # N2 sits between them
        @test interior_path(tree, 1, 2) == Int[]      # adjacent
        @test sort(interior_path(tree, 3, 4)) == [2]  # siblings meet at N2
        @test hop_distance(tree, 1, 3) == 2
        @test hop_distance(tree, 3, 4) == 2
        @test hop_distance(tree, 1, 1) == 0
    end

    @testset "meshed input is rejected" begin
        # Adding N3—N4 closes a cycle; the error bounds assume a radial feeder.
        meshed = Matrix(net.Lambda)
        meshed[3, 4] = meshed[4, 3] = 1
        @test_throws ErrorException orient_radial(meshed, 1)
    end

    @testset "admissible pairs" begin
        near = admissible_pairs(net, tree; hops=1, direction=:any)
        @test near[1, 2] && near[2, 1]                # adjacent, either way
        @test !near[1, 3]                             # two hops away
        @test all(near[i, i] for i in 1:4)            # a bus may always stay

        far = admissible_pairs(net, tree; hops=2, direction=:any)
        @test far[1, 3] && far[3, 1]
        @test far[3, 4] && far[4, 3]                  # siblings, two hops

        down = admissible_pairs(net, tree; hops=2, direction=:downstream)
        @test down[1, 3]                              # N3 absorbed upstream into N1
        @test !down[3, 1]                             # never the other way
        @test !down[3, 4]                             # siblings sit at equal depth
        @test_throws ErrorException admissible_pairs(net, tree; direction=:sideways)
    end

    @testset "assignment algebra" begin
        A = identity_assignment(net)
        @test reduction_ratio(A) == 0.0
        @test super_nodes(A) == [1, 2, 3, 4]

        assign!(A, 2, 3)                              # N3 -> N2
        @test A[2, 3] == 1 && A[3, 3] == 0
        @test super_nodes(A) == [1, 2, 4]
        @test reduction_ratio(A) == 0.25

        # Assigning N2 onward must carry N3 with it, not orphan it at a
        # bus that is no longer a super-node.
        assign!(A, 1, 2)
        @test A[1, 2] == 1 && A[1, 3] == 1
        @test A[2, 2] == 0 && A[2, 3] == 0
        @test super_nodes(A) == [1, 4]

        @test_throws ErrorException assign!(A, 1, 2)  # already reduced
        @test_throws ErrorException assign!(A, 2, 4)  # not a super-node
    end

    @testset "injections and voltage lifting" begin
        A = identity_assignment(net)
        assign!(A, 2, 3)

        # Single-phase: the expansion is the assignment matrix itself.
        @test Matrix(expand_assignment(A, net)) == A

        S_agg = aggregate_injections(A, net)
        @test S_agg[3, :] == [0, 0]                              # reduced bus emptied
        @test S_agg[2, 1] ≈ net.S[2, 1] + net.S[3, 1]            # load moved to N2
        @test sum(S_agg, dims=1) ≈ sum(net.S, dims=1)            # nothing lost

        V = ComplexF64[1.0; 0.99; 0.98; 0.97;;]
        V_lifted = lift_voltages(A, V, net)
        @test V_lifted[3] == V[2]                                # N3 inherits N2
        @test V_lifted[1] == V[1]
    end

    @testset "Kron reduction" begin
        Y_red, rows = kron_reduce(net, [1, 3, 4])
        @test rows == [1, 3, 4]
        @test size(Y_red) == (3, 3)

        # Eliminating N2 (degree 3) makes its neighbours a clique -- Lemma 1.
        adj = reduced_adjacency(Y_red, [1, 3, 4], net)
        @test adj == [0 1 1; 1 0 1; 1 1 0]

        @test_throws ErrorException kron_reduce(net, [3, 4])       # slack dropped
        @test_throws ErrorException kron_reduce(net, [1, 1, 3])    # duplicates
    end

    @testset "radialization" begin
        A = identity_assignment(net)
        assign!(A, 1, 2)                       # reduce N2 into N1 -> triangle
        @test !is_radial(net, A)

        # N2 has degree 3 in the sub-tree spanning {N1, N3, N4}, so it is the
        # one bus that must come back.
        @test critical_nodes(net, A) == [2]

        A_radial, critical = radialize(net, A)
        @test critical == [2]
        @test is_radial(net, A_radial)
        @test A_radial[2, 2] == 1
        @test reduction_ratio(A_radial) <= reduction_ratio(A)   # structure costs reduction

        # An already-radial reduction is left untouched.
        A_leaf = identity_assignment(net)
        assign!(A_leaf, 2, 3)                  # absorb a leaf; still a tree
        @test is_radial(net, A_leaf)
        @test critical_nodes(net, A_leaf) == Int[]
        @test radialize(net, A_leaf)[1] == A_leaf
    end

    @testset "fill-in threshold separates residue from coupling" begin
        Y_red, _ = kron_reduce(net, [1, 3, 4])

        # The Schur complement leaves rounding residue where the reduced network
        # has no branch. An absolute-zero threshold counts it as one, which is
        # what used to make radialized networks still report as meshed.
        Y_noisy = copy(Y_red)
        Y_noisy[1, 3] = Y_noisy[3, 1] = 1e-9 * maximum(abs, Y_red)
        @test reduced_adjacency(Y_noisy, [1, 3, 4], net) ==
              reduced_adjacency(Y_red, [1, 3, 4], net)
        @test reduced_adjacency(Y_noisy, [1, 3, 4], net; rtol=0.0)[1, 3] == 1

        # Transformers and regulators make Ybus legitimately asymmetric. Testing
        # each direction alone yields a non-symmetric adjacency, which Graph()
        # rejects -- so coupling is taken across both directions at once.
        Y_asym = copy(Y_red)
        Y_asym[1, 3] = 0
        adj = reduced_adjacency(Y_asym, [1, 3, 4], net)
        @test adj == transpose(adj)
        @test adj[1, 3] == 1                       # the [3,1] direction still couples
        @test OptiKRON.Graph(adj) isa OptiKRON.Graph   # no ArgumentError
    end

    @testset "sub-tree degree is local" begin
        # Degree inside the spanning sub-tree, not in the full feeder.
        members = spanning_subtree(tree, [1, 3, 4])
        @test members == [1, 2, 3, 4]
        degrees = subtree_degrees(tree, members)
        @test degrees[2] == 3                  # N2 joins all three
        @test degrees[1] == 1
    end

    @testset "three-phase network" begin
        net3 = read_network_csv(TEST3PH)
        @test is_three_phase(net3)
        @test nnodes(net3) == 3
        @test nphase_rows(net3) == 7           # abc + abc + a
        @test node_rows(net3) == [[1, 2, 3], [4, 5, 6], [7]]

        tree3 = orient_radial(net3)
        @test tree3.parents == [0, 1, 2]

        # N3 carries only phase a, so it may join either abc bus...
        ok = admissible_pairs(net3, tree3; hops=2, direction=:any)
        @test ok[2, 3] && ok[1, 3]
        # ...but an abc bus can never be absorbed into the single-phase N3.
        @test !ok[3, 2] && !ok[3, 1]

        # Expansion maps phase a of N3 onto phase a of N2 -- row 4, not 5 or 6.
        A3 = identity_assignment(net3)
        assign!(A3, 2, 3)
        E = expand_assignment(A3, net3)
        @test size(E) == (7, 7)
        @test E[4, 7] == 1
        @test iszero(E[5, 7]) && iszero(E[6, 7])

        S3 = aggregate_injections(A3, net3)
        @test S3[7, 1] == 0                              # N3 emptied
        @test S3[4, 1] ≈ net3.S[4, 1] + net3.S[7, 1]     # onto N2 phase a
        @test S3[5, 1] ≈ net3.S[5, 1]                    # phases b, c untouched
        @test sum(S3, dims=1) ≈ sum(net3.S, dims=1)
    end
end
