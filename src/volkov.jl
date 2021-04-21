struct VolkovPhases{Time,T,U}
    t::Time
    ∫A::Vector{T}
    ∫A²::Vector{U}
end

function VolkovPhases(F::ElectricFields.AbstractField, t::AbstractVector)
    At = typeof(vector_potential(F, t[1]))
    Tt = eltype(t)
    nt = length(t)
    ∫Av = zeros(At, nt)
    ∫A²v = zeros(Tt, nt)

    tend = t[end]
    Aend = ElectricFields.vector_potential(F, tend)

    A = t -> ElectricFields.vector_potential(F, t) - Aend
    ∫A = AccumulatedIntegral(A, At, Tt; t=tend)
    ∫A² = AccumulatedIntegral(t -> sum(abs2, A(t)), Tt, Tt; t=tend)

    for i = nt-1:-1:1
        set!(∫A, t[i])
        set!(∫A², t[i])
        ∫Av[i] = ∫A.∫f
        ∫A²v[i] = ∫A².∫f
    end

    VolkovPhases(t .- tend, ∫Av, ∫A²v)
end

VolkovPhases(F, ndt::Integer) =
    VolkovPhases(F, timeaxis(F, ndt))

volkov_phase_k²(k, v, i) = norm(k)^2*v.t[i]

volkov_phase_2kA(k::Number, v::VolkovPhases{<:Any,<:Number,<:Number}, i) =
    2k*v.∫A[i]

volkov_phase_2kA(k::SVector{3}, v::VolkovPhases{<:Any,<:Number,<:Number}, i) =
    2k[3]*v.∫A[i]

volkov_phase_2kA(k::SVector{3}, v::VolkovPhases{<:Any,<:SVector{3},<:Number}, i) =
    2dot(k, v.∫A[i])

volkov_phase_A²(v, i) = v.∫A²[i]

volkov_phase(k, v::VolkovPhases, i) =
    -(volkov_phase_k²(k, v, i) +
      volkov_phase_2kA(k, v, i) +
      volkov_phase_A²(v, i))/2

volkov_phase(k, v::VolkovPhases, i, j) =
    -(volkov_phase(k, v, i) -
      volkov_phase(k, v, j))

function stationary_momentum(v::VolkovPhases, i, j)
    τ = v.t[i] - v.t[j]
    -1/τ*(v.∫A[i]-v.∫A[j])
end
