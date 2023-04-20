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
numerically, with recursive time integrals limited by a `window`,
i.e. how many time steps are considered (default is from the beginning
of the pulse). Integrals over intermediate photoelectron momenta are
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

julia> # d will be a vector of scalars; by limiting the "memory" of the
       # integrals, we can include only the short trajectory.
       d = induced_dipole(Iₚ, F, ndt, window=flat_window(floor(Int, 0.65ndt)));
┌ Info: Induced dipole calculation
│   system =
│    1-channel System:
│     1. IonizationChannel: Iₚ = 0.5 Ha = 13.6055 eV
│
│   diagram =
│    Goldstone Diagram:
│       |0⟩
│       ╱ ╲⇜
│     1┃   │𝐩
│       ╲ ╱⇝
│       |0⟩
└
Progress: 100%|████████████████████████████████████████████████████████████████| Time: 0:00:00

julia> # d_all includes all trajectories.
       d_all = induced_dipole(Iₚ, F, ndt);
┌ Info: Induced dipole calculation
│   system =
│    1-channel SFA System:
│     1. IonizationChannel: Iₚ = 0.5 Ha = 13.6055 eV
│
│   diagram =
│    Goldstone Diagram:
│       |0⟩
│       ╱ ╲⇜
│     1┃   │𝐩
│       ╲ ╱⇝
│       |0⟩
└
Progress: 100%|████████████████████████████████████████████████████████████████| Time: 0:00:08

julia> # d2 will be a vector of 3d vectors
       d2 = induced_dipole(Iₚ, F2, ndt, window=flat_window(floor(Int, 0.65ndt)));
┌ Info: Induced dipole calculation
│   system =
│    1-channel System:
│     1. IonizationChannel: Iₚ = 0.5 Ha = 13.6055 eV
│
│   diagram =
│    Goldstone Diagram:
│       |0⟩
│       ╱ ╲⇜
│     1┃   │𝐩
│       ╲ ╱⇝
│       |0⟩
└
Progress: 100%|████████████████████████████████████████████████████████████████| Time: 0:00:00
```
![HHG example](figures/hhg_example.svg)

The dotted vertical line indicates the classical HHG cut-off:
```math
W_k^{\textrm{max}} = 3.17U_p + I_p.
```

## Photoelectron Spectra

### Attosecond Streaking

Attosecond streaking is a method for measuring, among other things,
the vector potential of a ultrashort IR pulse. We can simulate this
experiment as follows:

```jldoctest
julia> @field(pump) do
           ħω = 80u"eV"
           I₀ = 1e8u"W/cm^2"
           τ = 200u"as"
           toff = 400u"as"
           tmax = 500u"as"
           env = :trunc_gauss
       end;

julia> @field(probe) do
           λ = 800.0u"nm"
           I₀ = 1e12u"W/cm^2"
           τ = 2.66u"fs"
           toff = 4.0u"fs"
           tmax = 5.0u"fs"
           env = :trunc_gauss
       end;

julia> F = pump + delay(probe, 1.0u"fs");

julia> ndt = 100
100

julia> nk = 100
100

julia> nθ = 2 # Forward and backward spectra only
2

julia> k,kmag,θ = momentum_grid(50u"eV", 80u"eV", nk, nθ,
                                spacing=:energy);

julia> Iₚ = 14u"eV" # "Krypton"
14 eV

julia> c = photoelectron_spectrum(k, Iₚ, F, ndt);
┌ Info: Photoelectrum spectrum calculation
│   system =
│    1-channel System:
│     1. IonizationChannel: Iₚ = 0.5144905104572604 Ha = 13.99980128005251 eV
│
│   diagram =
│    Goldstone Diagram:
│       |0⟩
│       ╱ ╲⇜
│     1┃   │𝐤
│
└   length(k) = 200
Progress: 100%|████████████████████████████████████████████████████████████████| Time: 0:00:01
```
![Streaking single delay](figures/streaking_single_delay.svg)

By sweeping the pump–probe delay, we can construct the following
streaking spectrogram:
![Streaking spectrogram](figures/streaking_spectrogram.svg)
