using Physical
using Calculus

a₀ = a_bohr.value
c_au = 1.0/alpha_fine_structure

function field_params_au(λ, I)
    ω = 2π*c_au/(λ/a₀)
    T = λ/(a₀*c_au)
    I_au = I/3.5094452e16
    E₀ = √(I_au)
    A₀ = E₀/ω
    ω,T,A₀
end

trapezoidal(flat,ramp) = t ->
((t .> 0) .== (t .<= ramp)) .* (t/ramp) +
((t .> ramp) .== (t .<= ramp+flat)) +
((t .> ramp+flat) .== (t .<= 2ramp+flat)) .* (2ramp + flat - t)/ramp

function field(λ, I;
               carrier = (ω,t) -> cos(ω*t),
               env = t -> 1)
    ω,T,A₀ = field_params_au(λ, I)
    A = t -> A₀*env(t/T).*carrier(ω,t)
    E = t -> -map(A', t)
    ω,T,A,E
end

export trapezoidal, field
