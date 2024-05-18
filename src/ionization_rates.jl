module IonizationRates
using SpecialFunctions
using HCubature
using ElectricFields

# * Common

Γ(x) = SpecialFunctions.gamma(x)

@doc raw"""
```math
\begin{equation}
\tag{ADK19}
\abs{C_{n^*\ell^*}}^2 =
\frac{2^{2n^*}}{n^*\Gamma(n^*+\ell^*+1)\Gamma(n^*-\ell^*)}
\end{equation}
```
"""
C2(n, ℓ) = 2.0^(2n)/(n*Γ(n+ℓ+1)*Γ(n-ℓ))

@doc raw"""
    Etilde(E0, Iₚ)
```math
\tilde{E} = \frac{E_0}{2(2I_p)^{3/2}} \equiv
\frac{\sqrt{I}}{2(2I_p)^{3/2}}.
```
"""
Etilde(E0, Iₚ) = E0/(2(2Iₚ)^1.5)

@doc raw"""
    f(ℓ, m)

Common prefactor appearing in e.g Eqs. (PPT54,59), (ADK12)

```math
f(\ell, m) =
\frac{(2\ell+1)(\ell+\abs{m})!}
{2^{\abs{m}}(\abs{m})!(\ell-\abs{m})!}
```
"""
function f(ℓ::Integer,m::Integer)
    ℓ >= 0 || error("Invalid value of ℓ ($(ℓ) < 0)")
    abs(m) <= ℓ || error("Invalid value of m ($(abs(m)) > $(ℓ))")
    (2ℓ+1)*(factorial(ℓ+abs(m)) /
            (2^abs(m)*factorial(abs(m))*factorial(ℓ-abs(m))))
end

Uₚ(I,ω) = I/4ω^2

keldysh(Iₚ, Uₚ) = √(Iₚ/2Uₚ)
keldysh(Iₚ, I, ω) = keldysh(Iₚ, Uₚ(I, ω))

@doc raw"""
    effective_n(Iₚ, Z)

Effective principal quantum number for ionization potential `Iₚ` and
residual nuclear charge `Z`.

```math
\begin{equation}
\tag{ADK2}
n^* = \frac{Z}{\sqrt{2I_p}}
\end{equation}
```
"""
effective_n(Iₚ, Z) =
    Z/√(2Iₚ)

# * PPT

@doc raw"""
    β(γ)

```math
\begin{equation}
\tag{PPT56ff}
\beta(\gamma) =
\frac{2\gamma}{\sqrt{1+\gamma^2}}
\end{equation}
```

Defined just below Eq. (PPT56).
"""
function β(γ::T) where T
    z = inv(γ)
    if abs(z) < eps(real(T))
        T(2) - z^2
    else
        2γ/√(1+γ^2)
    end
end

@doc raw"""
    α(γ)

```math
\begin{equation}
\tag{PPT34}
\alpha(\gamma) =
2\left(
\arcsinh\gamma -
\frac{\gamma}{\sqrt{1+\gamma^2}}
\right) \equiv
2\arcsinh\gamma -
\beta(\gamma)
\end{equation}
```
"""
α(γ) = 2asinh(γ) - β(γ)

@doc raw"""
    g(γ)

```math
\begin{equation}
\tag{PPT33}
g(\gamma) =
\frac{3}{2\gamma}
\left[
\left(
1+\frac{1}{2\gamma^2}
\right)
\arcsinh\gamma-
\frac{\sqrt{1+\gamma^2}}{2\gamma}
\right] \equiv
\frac{3}{2\gamma}
\left[
\left(
1+\frac{1}{2\gamma^2}
\right)
\arcsinh\gamma-
\frac{1}{\beta(\gamma)}
\right]
\end{equation}
```
"""
g(γ) = 3/(2γ)*((1+1/(2γ^2))*asinh(γ) - 1/β(γ))

@doc raw"""
    ν(Iₚ, I, ω)

Number of photons of energy ``\omega`` needed to reach above the
ionization potential ``I_p`` in the presence of a strong dressing
field of intensity ``I``.

```math
\begin{equation}
\tag{PPT26}
\nu = \frac{I_p}{\omega}
\left(1+\frac{1}{2\gamma^2}\right)
\end{equation}
```
"""
ν(Iₚ, I, ω) = Iₚ/ω * (1.0 + 1/2keldysh(Iₚ,I,ω)^2)

@doc raw"""
    w(m,x)

```math
\begin{equation}
\tag{PPT56}
w_m(x) =
\frac{x^{2\abs{m}+1}}{2}
\int_0^1\diff{t}
\frac{\exp(-x^2t)t^{\abs{m}}}{\sqrt{1-t}}
\end{equation}
```
"""
function w(m,x)
    val, = hquadrature(t -> exp(-x^2*t)*t^abs(m)/√(1-t), 0, 1)
    x^(2abs(m)+1)/2*val
end

@doc raw"""
    A(m, Iₚ, I, ω[; tol = √(eps()), maxterms=typemax(Int)])

```math
\begin{equation}
\tag{PPT55}
A_m(\omega,\gamma) =
\frac{4}{\sqrt{3\pi}}
\frac{1}{\abs{m}!}
\frac{\gamma^2}{1+\gamma^2}
\sum_{n\geq\nu}^\infty
\exp[-\alpha(n-\nu)]
w_m[\sqrt{\beta(n-\nu)}],
\end{equation}
```
where ``\nu`` is computed using [`ν`](@ref).
"""
function A(m, Iₚ, I, ω; tol = √(eps()), maxterms=typemax(Int))
    γ = keldysh(Iₚ,I,ω)
    nν = ν(Iₚ, I, ω)
    S = 0
    dS = Inf
    kappa = ceil(Int, nν)
    ii = 0
    α_γ = α(γ)
    β_γ = β(γ)
    while dS/S > tol && ii < maxterms
        dk = kappa - nν
        dS = exp(-α_γ*dk)*w(m,√(β_γ*dk))
        S += dS
        kappa += 1
        ii += 1
    end
    (4/√(3π))/factorial(abs(m)) * γ^2/(1+γ^2) * S
end

@doc raw"""
    PPT(Iₚ, I, ω, ℓ, m[, Z=1])

Compute the strong-field photoionization rate from an initial state
with ionization potential `Iₚ`, angular quantum numbers `ℓ` & `m`, and
a binding potential of charge `Z`, when subjected to monochromatic
light of intensity `I` and angular frequency `ω`, using the PPT
formalism (slightly generalized to effective principal quantum numbers
``n^*``):

```math
\begin{equation}
\tag{PPT54,ADK1}
w_{\ell m}(F, \omega) =
I_p
\abs{C_{n^*\ell^*}}^2
\sqrt{\frac{3}{2\pi}}
f(\ell, m)
\tilde{E}^{-(2n^* - \abs{m} - 3/2)}
(1+\gamma^2)^{-n^* + \abs{m}/2 + 3/4}
A_m(\omega,\gamma)
\exp\left[
-\frac{g(\gamma)}{3\tilde{E}}
\right],
\end{equation}
```
where ``\gamma`` is the Keldysh parameter, and ``n^*`` the [`effective_n`](@ref).
"""
function PPT(Iₚ, I, ω, ℓ, m, Z=1; kwargs...)
    γ = keldysh(Iₚ, I, ω)
    n⭑ = effective_n(Iₚ, Z)

    pre = √(3/2π)*C2(n⭑,n⭑-1)*f(ℓ,m)*Iₚ
    Ẽ = Etilde(√(I), Iₚ)

    Am = A(0, Iₚ, I, ω; kwargs...)

    w = pre*Ẽ^(-(2n⭑-abs(m)-3/2))*(1+γ^2)^(-n⭑+abs(m)/2+3/4)*Am*exp(-g(γ)/(3Ẽ))

    isnan(w) ? zero(w) : w
end

@doc raw"""
    Keldysh(Iₚ, I, ω)

Compute the strong-field photoionization rate from an initial state
with ionization potential `Iₚ`, when subjected to monochromatic light
of intensity `I` and angular frequency `ω`, using the Keldysh
formalism:

```math
\begin{equation}
\tag{Keldysh1}
w(F, \omega) =
\exp\left\{
-\frac{2I_p}{\omega}
\left[
  \left(1+\frac{1}{2\gamma^2}\right)\operatorname{arcsinh}\gamma -
  \frac{\sqrt{\gamma^2+1}}{2\gamma}
  \right]
\right\},
\end{equation}
```
where ``\gamma`` is the Keldysh parameter.
"""
function Keldysh(Iₚ, I, ω; η=1, kwargs...)
    γ = keldysh(Iₚ, I, ω)
    η*exp(-2*Iₚ/ω*((1 + 1/(2γ^2)) * asinh(γ) - 1/β(γ)))
end

function ionization_yield(F::ElectricFields.LinearField, (tmin,tmax)::Tuple{<:Number,<:Number}, Iₚ, args...; model=:ppt, kwargs...)
    Model = if model == :ppt
        PPT
    elseif model == :keldysh
        Keldysh
    else
        throw(ArgumentError("Unknown ionization model $(model)"))
    end

    f(t) = Model(Iₚ, abs2(field_amplitude(F, t)), ω, args...; kwargs...)

    s = span(F)
    a = max(tmin, s.left)
    b = min(tmax, s.right)

    ω = photon_energy(F)

    1 - exp(-first(hquadrature(f, a, b)))
end

function ionization_yield(F::ElectricFields.LinearField, Iₚ, args...; kwargs...)
    s = span(F)
    ionization_yield(F, (s.left, s.right), Iₚ, args...; kwargs...)
end

"""
    cumulative_integral!(v, f, t)

Compute the cumulative integral, such that at exit,

v[i] = ∫_{t[1]}^{t[i]} dt f(t).

We do this by recursively halving the interval until it is small
enough for brute force quadrature.
"""
function cumulative_integral!(v, f, t; min_length=3, kwargs...)
    @assert min_length > 2
    n = length(t)
    if n ≤ min_length
        # @info "Base case"
        # We may not write to v[1], since it is v[end] in the left
        # neighbouring interval, and the contribution from this
        # interval is zero anyway.
        for i = 2:n
            v[i] = f(t[1],t[i])
        end
    else
        s = n ÷ 2
        sel1 = 1:s
        sel2 = s:n
        # @info "Subdividing integral" n s sel1 sel2
        v1 =
        @sync begin
            Threads.@spawn cumulative_integral!(view(v, sel1), f, view(t, sel1); min_length=min_length, kwargs...)
            Threads.@spawn cumulative_integral!(view(v, sel2), f, view(t, sel2); min_length=min_length, kwargs...)
        end
        v[s+1:n] .+= v[s]
    end
    v
end

function ionization_yield(F::ElectricFields.LinearField, t::AbstractVector{T}, args...; kwargs...) where T
    f(a,b) = ionization_yield(F, a, b, args...; kwargs...)
    a,b = length(t) > 1 ? (t[1],t[2]) : (one(T),one(T)+eps(T))
    cumulative_integral!(zeros(typeof(f(a,b)), length(t)), f, t; kwargs...)
end

end
