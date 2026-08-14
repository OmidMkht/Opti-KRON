# --------------------------------------------------------------------------- #
# AC power flow, and the consistency check that goes with it.
#
# Opti-KRON linearises around an operating point, so `V` is problem data on the
# same footing as `S`. The solver is a Z-bus fixed point: hold the slack, and
# repeatedly solve
#
#     Ybus[ns,ns] V[ns] = conj(S[ns] ./ V[ns]) - Ybus[ns,slack] V[slack]
#
# refactorising nothing. About seven iterations on the shipped feeders.
#
# It assumes constant-power wye injections and one slack -- deliberately not a
# general distribution power flow. No delta connections, ZIP loads or tap logic;
# a feeder needing those is solved elsewhere and its voltages handed over, which
# is why `V` is an input the package accepts rather than insists on.
#
# The start must respect phase angles. A three-phase slack sits at 1<0, 1<-120,
# 1<120, so starting every row at 1<0 is not a flat start but a zero-sequence one
# -- on R100 that left the iteration still moving after 100 passes, against 1e-13
# in seven from the correct angles.
# --------------------------------------------------------------------------- #

"Per-unit voltage of a balanced source on `phases`: a at 0 degrees, b at -120, c at +120."
balanced_source(phases::AbstractVector{Int}, magnitude::Number=1.0) =
    [ComplexF64(magnitude) * cis(-2pi * (p - 1) / 3) for p in phases]

"The phase indices (1=a, 2=b, 3=c) carried by bus `i`."
phase_indices(net::Network, i::Int) = [p for p in 1:3 if net.phases[p, i]]

"""
    slack_voltage(net; magnitude=1.0) -> Vector{ComplexF64}

Default slack voltage: a balanced source across whatever phases the slack bus
carries. Angles follow the *phase*, not the position -- a slack on `a` and `c` is
`1<0` and `1<120`, not `1<0` and `1<-120`.
"""
slack_voltage(net::Network; magnitude::Number=1.0) =
    balanced_source(phase_indices(net, net.slack), magnitude)

"""
    powerflow(net; vslack, maxiter, tol) -> Matrix{ComplexF64}

Solve the AC power flow for every scenario, laid out like `net.S`: one row per
node-phase, one column per scenario. `vslack` defaults to
[`slack_voltage`](@ref).

Errors rather than returning a half-converged answer, which would silently
invalidate every bound derived from it.

`tol` is 1e-9 pu because the iteration cannot beat the linear solve underneath:
near-zero-impedance jumpers push cond(Ybus) to ~1e9, and on R300 the step bottoms
out near 2e-10 -- still seven orders below the smallest budget anyone uses.
"""
function powerflow(net::Network;
    vslack::AbstractVector=slack_voltage(net),
    maxiter::Int=200,
    tol::Real=1e-9)

    nph = nphase_rows(net)
    slack_rows = node_rows(net)[net.slack]
    length(vslack) == length(slack_rows) ||
        error("vslack has $(length(vslack)) entries but the slack bus carries " *
              "$(length(slack_rows)) phase(s). A balanced three-phase slack needs three.")

    other = setdiff(1:nph, slack_rows)
    Y = Matrix{ComplexF64}(net.Ybus)
    factorization = lu(Y[other, other])
    coupling = Y[other, slack_rows]

    V = _balanced_start(net, vslack, slack_rows)
    V_slack = V[slack_rows, :]

    previous_step = Inf
    stalled = 0
    for iteration in 1:maxiter
        current = conj.(net.S ./ V)
        V_next = factorization \ (current[other, :] .- coupling * V_slack)
        step = maximum(abs.(V_next .- V[other, :]))
        V[other, :] .= V_next
        step < tol && return V

        # A step that stops shrinking means the accuracy floor of the linear
        # solve, not a bad feeder. Say so rather than spinning to `maxiter`.
        stalled = step < 0.99 * previous_step ? 0 : stalled + 1
        previous_step = step
        stalled >= 3 && error(
            "Power flow stalled at a step of $step per unit after $iteration " *
            "iterations, above the requested tol=$tol. This is the accuracy floor " *
            "of the admittance solve, not a modelling error -- a near-zero-impedance " *
            "branch makes Ybus ill-conditioned. Retry with tol above $step.")
    end

    error("Power flow did not converge in $maxiter iterations " *
          "(last step $previous_step). The feeder may be overloaded, or the " *
          "injections may not be reachable from this slack voltage.")
end

"Every row starts at its own phase's source voltage, which is what makes this converge."
function _balanced_start(net::Network, vslack, slack_rows)
    V = zeros(ComplexF64, nphase_rows(net), nscenarios(net))
    slack_phases = phase_indices(net, net.slack)
    reference = Dict(p => vslack[k] for (k, p) in enumerate(slack_phases))
    fallback = first(vslack)

    for i in 1:nnodes(net), (k, p) in enumerate(phase_indices(net, i))
        V[node_rows(net)[i][k], :] .= get(reference, p, fallback)
    end
    V[slack_rows, :] .= vslack
    return V
end

"""
    powerflow_residual(net, V) -> Float64

How far `V` is from solving the power flow, in per unit: the size of one Z-bus
step from `V`. A converged point gives ~1e-12, a flat start ~1.6.

Measured in *voltage* rather than as the mismatch current `Ybus V - conj(S./V)`,
which the smallest impedance dominates: on R100 a jumper of 4.7e6 pu admittance
turns a 5e-7 voltage difference into a mismatch of 2.5, so a correct `V` and one
wrongly scaled by sqrt(3) score 2.5 and 4.3 and cannot be told apart. Same reason
the MILP works with Z rather than Ybus.
"""
function powerflow_residual(net::Network, V::AbstractMatrix)
    size(V) == (nphase_rows(net), nscenarios(net)) ||
        error("V is $(size(V)) but the network implies " *
              "$((nphase_rows(net), nscenarios(net)))." )

    slack_rows = node_rows(net)[net.slack]
    other = setdiff(1:nphase_rows(net), slack_rows)
    Y = Matrix{ComplexF64}(net.Ybus)

    current = conj.(net.S ./ V)
    V_next = Y[other, other] \ (current[other, :] .- Y[other, slack_rows] * V[slack_rows, :])
    return maximum(abs.(V_next .- V[other, :]))
end
