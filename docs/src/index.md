```@meta
CurrentModule = StrongFieldApproximation
```

# StrongFieldApproximation.jl

Documentation for [StrongFieldApproximation.jl](https://github.com/jagot/StrongFieldApproximation.jl).

[Hartree atomic
units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used
throughout, unless otherwise specified.

At the moment, two kind of observables are supported: induced dipole
moment and photoelectron spectra. Time integrals are performed
numerically, with recursive time integrals limited by a `memory`,
i.e. how many time steps are considered (default is one cycle of the
fundamental). Integrals over intermediate photoelectron momenta are
performed using the saddle-point method, i.e. given two times, the
stationary momentum is given by
```math
\saddle{\vec{p}} = -\frac{1}{t_2-t_1}
\int_{t_1}^{t_2}\diff{\tau}
[\vec{A}(\tau)-\vec{A}(T)],
```
where ``\vec{A}(t)`` is the vector potential of the external field,
and ``T`` is the reference time, at which we usually have
``\vec{A}(T)=0``. The integral is thus approximated as
```math
\int\diff{\vec{p}}
f(\vec{p})
\ce^{-\im S(\vec{p},a,b)}
\approx
\left[
\frac{2\pi}{\im(a-b) + \epsilon}
\right]^{3/2}
f(\saddle{\vec{p}})
\ce^{-\im S(\saddle{\vec{p}},a,b)}
```
where
```math
S(\vec{p},a,b)\defd
-\frac{1}{2}
\int_a^b\diff{\tau}
[\vec{p}+\vec{A}(\tau)-\vec{A}(T)]^2,
```
and ``\epsilon`` is an infinitesimal regulator.

## High-order Harmonic Generation (HHG)

```jldoctest
julia> @field(F) do # Linear polarization
           λ = 800.0u"nm"
           I₀ = 3e14u"W/cm^2"
           τ = 6.2u"fs"
           toff = 10.0u"fs"
           tmax = 13.0u"fs"
           env = :trunc_gauss
       end;

julia> @field(F2) do # Slightly elliptical polarization
           λ = 800.0u"nm"
           I₀ = 3e14u"W/cm^2"
           τ = 6.2u"fs"
           toff = 10.0u"fs"
           tmax = 13.0u"fs"
           env = :trunc_gauss
           ξ = 0.2
       end;

julia> ndt = 300 # Steps per cycle
300

julia> Iₚ = 0.5 # Hydrogen
0.5

julia> # d will be a vector of scalars
       system,diagram,d = induced_dipole(Iₚ, F, ndt, memory=floor(Int, 0.65ndt));
┌ Info: Induced dipole calculation
│   system =
│    1-channel System:
│     1. IonizationChannel: Iₚ = 0.5 Ha = 13.6055 eV
│
│    Linearly polarized field with
│      - I₀ = 8.5484e-03 au = 3.0e14 W cm⁻² =>
│        - E₀ = 9.2457e-02 au = 47.5435 GV m⁻¹
│        - A₀ = 1.6234 au
│      – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV)
│      – and a Truncated Gaussian envelope of duration 256.3165 jiffies = 6.2000 fs (intensity FWHM; turn-off from 10.0000 fs to 13.0000 fs)
│      – Uₚ = 0.6588 Ha = 17.9276 eV => α = 28.5030 Bohr = 1.5083 nm
│   diagram =
│    Goldstone Diagram:
│       |0⟩
│      ╱   ╲⇜
│     1┃   │𝐤
└
Progress: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:00

julia> # d2 will be a vector of 3d vectors
       system2,diagram2,d2 = induced_dipole(Iₚ, F2, ndt, memory=floor(Int, 0.65ndt));
┌ Info: Induced dipole calculation
│   system =
│    1-channel System:
│     1. IonizationChannel: Iₚ = 0.5 Ha = 13.6055 eV
│
│    Transversely polarized field with
│      - I₀ = 8.5484e-03 au = 3.0e14 W cm⁻² =>
│        - E₀ = 9.2457e-02 au = 47.5435 GV m⁻¹
│        - A₀ = 1.6234 au
│      – a Elliptical carrier with ξ = 0.20 (right) @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV)
│      – and a Truncated Gaussian envelope of duration 256.3165 jiffies = 6.2000 fs (intensity FWHM; turn-off from 10.0000 fs to 13.0000 fs)
│      – Uₚ = 0.6588 Ha = 17.9276 eV => α = 28.5030 Bohr = 1.5083 nm
│   diagram =
│    Goldstone Diagram:
│       |0⟩
│      ╱   ╲⇜
│     1┃   │𝐤
└
Progress: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:00
```
![HHG example](figures/hhg_example.svg)

The dotted vertical line indicates the classical HHG cut-off:
```math
W_k^{\textrm{max}} = 3.17U_p + I_p.
```

## Photoelectron Spectra
