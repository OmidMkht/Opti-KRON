using JuMP

@testset "OptiKRON.solver" begin

    @testset "selection" begin
        _, auto_name = select_optimizer()
        @test auto_name in (:gurobi, :highs)
        @test auto_name === (gurobi_available() ? :gurobi : :highs)

        # HiGHS is a hard dependency, so forcing it always works.
        _, forced = select_optimizer(prefer=:highs)
        @test forced === :highs

        # :gurobi must never silently fall back -- a benchmark result that
        # quietly came from a different solver would be misleading.
        if gurobi_available()
            @test select_optimizer(prefer=:gurobi)[2] === :gurobi
        else
            @test_throws ErrorException select_optimizer(prefer=:gurobi)
        end

        @test_throws ErrorException select_optimizer(prefer=:cplex)
    end

    @testset "both solvers reach the same optimum" begin
        function solve_with(factory)
            m = Model(factory)
            @variable(m, 0 <= x[1:3] <= 1, Int)
            @constraint(m, sum(x) <= 2)
            @objective(m, Max, 3x[1] + 2x[2] + x[3])
            optimize!(m)
            return termination_status(m), objective_value(m)
        end

        status_h, obj_h = solve_with(select_optimizer(prefer=:highs)[1])
        @test status_h == MOI.OPTIMAL
        @test obj_h ≈ 5.0

        if gurobi_available()
            status_g, obj_g = solve_with(select_optimizer(prefer=:gurobi)[1])
            @test status_g == MOI.OPTIMAL
            @test obj_g ≈ obj_h
        end
    end

    @testset "attributes are solver-independent" begin
        # Same knobs, both backends -- these go through MOI rather than
        # solver-specific parameter names.
        for backend in (gurobi_available() ? (:highs, :gurobi) : (:highs,))
            factory, _ = select_optimizer(prefer=backend, time_limit=30, mip_gap=1e-3)
            m = Model(factory)
            @variable(m, 0 <= x <= 4, Int)
            @objective(m, Max, x)
            optimize!(m)
            @test termination_status(m) == MOI.OPTIMAL
            @test value(x) ≈ 4.0
            @test MOI.get(m, MOI.TimeLimitSec()) ≈ 30.0
        end
    end

    @testset "probe cache is resettable" begin
        first_probe = gurobi_available()
        reset_solver_cache!()
        @test gurobi_available() === first_probe
    end
end
