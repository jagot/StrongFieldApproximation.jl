module IonizationRates
using SpecialFunctions
using HCubature

# * Common

Γ(x) = SpecialFunctions.gamma(x)

C2(n, ℓ) = 2.0^(2n)/(n*Γ(n+ℓ+1)*Γ(n-ℓ))

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
β(γ) = 2γ/√(1+γ^2)

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

function PPT(Iₚ, I, ω, ℓ, m, Z=1)
    γ = keldysh(Iₚ, I, ω)
    n⭑ = effective_n(Iₚ, Z)

    pre = √(3/2π)*C2(n⭑,n⭑-1)*f(ℓ,m)*Iₚ
    Ẽ = Etilde(√(I), Iₚ)

    Am = A(0, Iₚ, I, ω)

    w = pre*Ẽ^(-(2n⭑-abs(m)-3/2))*(1+γ^2)^(-n⭑+abs(m)/2+3/4)*Am*exp(-g(γ)/(3Ẽ))

    isnan(w) ? zero(w) : w
end

end
