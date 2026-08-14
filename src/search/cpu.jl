# --------------------------------------------------------------------------- #
# CPU backend for the q = 1 search.
#
# The arithmetic is the rank-1 Zbus update in gpu.jl; this only says where it
# runs. Host and device share one implementation on purpose -- two would drift,
# and a disagreement would then be ambiguous between hardware and code.
# --------------------------------------------------------------------------- #

_search_context(::Val{:cpu}, Z::Matrix{ComplexF64}, Vcur::Matrix{ComplexF64},
    absV::Matrix{Float64}, batch::Int) = _rank1_context(Z, Vcur, absV, batch, :cpu)
