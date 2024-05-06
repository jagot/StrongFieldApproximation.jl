struct VolkovPhases{ReferenceTime,∫At,∫A²t}
    T::ReferenceTime
    ∫A::∫At
    ∫A²::∫A²t
end

function _VolkovPhases(F::ElectricFields.AbstractField, t, tend)
    Aend = ElectricFields.vector_potential(F, tend)

    A = t -> ElectricFields.vector_potential(F, t) - Aend
    ∫A = AutomaticIntegration(A, t, tend)
    square(z) = z^2
    ∫A² = AutomaticIntegration(t -> sum(square, A(t)), t, tend)

    VolkovPhases(tend, ∫A, ∫A²)
end

VolkovPhases(F::ElectricFields.AbstractField, t::AbstractVector) =
    _VolkovPhases(F, t, t[end])

VolkovPhases(F::ElectricFields.AbstractField, tre::AbstractVector, tim::AbstractVector) =
    _VolkovPhases(F, ComplexPlane(tre, tim), tre[end])

volkov_phase_k²(k, v::VolkovPhases, t) = norm(k)^2*t

kdotA(k::Number, A::Number) = k*A
kdotA(k::SVector{3}, A::Number) = k[3]*A
kdotA(k::SVector{3}, A::SVector{3}) = dot(k, A)
volkov_phase_2kA(k, v::VolkovPhases, t) = 2kdotA(k, v.∫A(t))

volkov_phase_A²(v::VolkovPhases, t) = v.∫A²(t)

function volkov_phase(k, v::VolkovPhases, t)
    # We only have to subtract the reference time from the integral
    # over k², since it has been properly accounted for when setting
    # up the integrals over A and A².
    -(volkov_phase_k²(k, v, t - v.T) +
      volkov_phase_2kA(k, v, t) +
      volkov_phase_A²(v, t))/2
end

volkov_phase(k, v::VolkovPhases, a, b) =
    -(volkov_phase(k, v, a) -
      volkov_phase(k, v, b))

function stationary_momentum(v::VolkovPhases, a, b)
    τ = a - b
    -1/τ*v.∫A(b, a)
end
