# * System
struct System{T,IonizationChannels<:AbstractVector{<:IonizationChannel{T}},
              Couplings<:AbstractVector{<:AbstractMatrix{<:AbstractCoupling}},
              Field<:ElectricFields.AbstractField,
              Times<:AbstractVector,
              Volkov<:VolkovPhases}
    ionization_channels::IonizationChannels

    couplings::Couplings

    F::Field
    ndt::Int
    t::Times
    dt::T

    volkov::Volkov
end

function System(ionization_channels::AbstractVector{<:IonizationChannel},
                couplings::AbstractVector,
                F::ElectricFields.AbstractField, ndt::Integer)
    t = timeaxis(F, ndt)
    volkov = VolkovPhases(F, ndt)

    System(ionization_channels, couplings, F, ndt, t, step(t), volkov)
end

System(ionization_channels::AbstractVector{<:IonizationChannel},
       F::ElectricFields.AbstractField, ndt::Integer;
       couplings=Matrix{AbstractCoupling}[], kwargs...) =
           System(ionization_channels, couplings, F, ndt)

System(ionization_channel::IonizationChannel, args...; kwargs...) =
    System([ionization_channel], args...; kwargs...)

System(Iₚ::Number, F, args...; kwargs...) =
    System(IonizationChannel(Iₚ, F), F, args...; kwargs...)

function Base.show(io::IO, system::System)
    n = length(system.ionization_channels)
    write(io, "$(n)-channel System")
end

function Base.show(io::IO, mime::MIME"text/plain", system::System)
    println(io, system, ":")
    for (i,c) in enumerate(system.ionization_channels)
        println(io, " ", i, ". ", c)
    end
    println(io)
    show(io, mime, system.F)
end

canonical_momentum_conservation(system::System, which::Integer) =
    iszero(which) ?
    CanonicalMomentumConservation() :
    canonical_momentum_conservation(first(system.couplings[which]))

# * Intermediate momenta

struct IntermediateMomentum{I₁<:Union{Nothing,Integer},I₂<:Union{Nothing,Integer},Times,Volkov}
    i₁::I₁ # Start time
    i₂::I₂ # End time = return time
    t::Times
    volkov::Volkov
end

IntermediateMomentum(i₁, i₂, system::System) =
    IntermediateMomentum(i₁, i₂, system.t, system.volkov)

Base.show(io::IO, 𝐩::IntermediateMomentum) =
    write(io, "IntermediateMomentum($(𝐩.i₁)..$(𝐩.i₂), ...)")


momentum(𝐩::IntermediateMomentum) = stationary_momentum(𝐩.volkov, 𝐩.i₂, 𝐩.i₁)
momentum(𝐤) = 𝐤

fix_momentum(𝐩::IntermediateMomentum{Nothing,<:Integer}, system::System, i) =
    IntermediateMomentum(i, 𝐩.i₂, system)

fix_momentum(𝐩::IntermediateMomentum{<:Integer,Nothing}, system::System, i) =
    IntermediateMomentum(𝐩.i₁, i, system)

fix_momentum(𝐤, args...) = 𝐤


excursion(𝐩::IntermediateMomentum) = 𝐩.t[𝐩.i₂] - 𝐩.t[𝐩.i₁]
prefactor(𝐩::IntermediateMomentum; ϵ=1e-2*√(eps(eltype(𝐩.t)))) where T = (2π/(im*excursion(𝐩) + ϵ))^(3/2)

prefactor(::Any) = true

kinematic_momentum(k::Number, A::Number) = k + A
kinematic_momentum(k::SVector{3}, A::SVector{3}) = k + A
kinematic_momentum(k::SVector{3}, A::Number) = SVector{3}(k[1], k[2], k[3]+A)

# * Diagrams

struct Diagram{Couplings}
    path::Vector{Tuple{Int,Int}}
    couplings::Couplings
    function Diagram(path, couplings::Couplings) where Couplings
        if !isempty(path)
            α,which = last(path)
            which == 0 ||
                throw(ArgumentError("Non-empty diagrams must have photoionization as the first interaction"))
        end
        new{Couplings}(path, couplings)
    end
end

Diagram(path::Tuple{Int,Int}, couplings) =
    Diagram([path], couplings)

function Diagram(path::AbstractVector, system::System)
    if isempty(path)
        length(system.ionization_channels) > 1 &&
            @warn "More than one ionization channel present in system, choosing the first"
        path = [(1,0)]
    end
    Diagram(path, [eltype(c) for c in system.couplings])
end

Diagram(system::System) = Diagram([], system)

for f in [:length, :first, :firstindex, :lastindex, :isempty]
    @eval Base.$f(d::Diagram) = $f(d.path)
end
Base.getindex(d::Diagram, i) = Diagram(d.path[i], d.couplings)

draw_exciton(io, ion, nd, electron="") =
    println(io, lpad(ion, nd)*"┃   │$(electron)")

draw_coulomb_interaction(io, nd) =
    println(io, repeat(" ",nd)*"┠┈┈┈┤")

function draw_dipole_interaction(io, nd)
    printstyled(io, lpad("⇝",nd), color=:light_red)
    println(io, "┃   │")
end

function Base.show(io::IO, ::MIME"text/plain", d::Diagram)
    println(io, "Goldstone Diagram:")
    if isempty(d)
        println(io, "|0⟩")
        return
    end
    nd = maximum(iw -> length(digits(iw[1])), d.path)+1
    println(io, lpad("|0⟩", nd+4))
    print(io, lpad("╱   ╲", nd+5))
    electrons = ["𝐩","𝐪"]
    cur_electron = 1
    for (i,(ion,which)) in Iterators.reverse(enumerate(d.path))
        which == 0 && printstyled(io, "⇜", color=:light_red)
        i == length(d.path) && println(io)
        if !iszero(which)
            c = d.couplings[which]
            if c <: CoulombCoupling
                draw_coulomb_interaction(io, nd)
            elseif c <: DipoleCoupling
                draw_dipole_interaction(io, nd)
            end
        end
        electron = if i == 1
            "𝐤"
        elseif !iszero(which) && d.couplings[which] <: DipoleCoupling
            ""
        else
            i = mod1(cur_electron,length(electrons))
            e = electrons[i]*repeat("′", fld1(cur_electron,length(electrons))-1)
            cur_electron += 1
            e
        end

        draw_exciton(io, ion, nd, electron)
    end
end

function get_interaction(system::System, diagram::Diagram)
    α, which = first(diagram)
    target_channel = system.ionization_channels[α]
    if iszero(which)
        @assert length(diagram) == 1
        return α, (𝐤, 𝐩, i) -> source_term(target_channel, system.t[i], 𝐤)
    end
end

# * Recursions

@doc raw"""
    recurse_common(fun, system::System, 𝐤, 𝐩, iref, irange, α, path; ϵ, memory)

At time ``t`` = `system.t[iref]`, we have an interaction, which
involves the photoelectrons and the ions:

```math
\begin{equation*}
\bra{\alpha}
\matrixel{\vec{k}}{V(t)}{\vec{p}}
\ket{\beta}
\end{equation*}
```

``\alpha`` and ``\vec{k}`` are given, i.e. this is what we want the
interaction to result in. If ``V(t)`` does not change the
photoelectron momentum, we have ``\vec{p}=\vec{k}``, but otherwise, we
have to determine the intermediate momentum ``\vec{p}`` from the
semiclassical action ending at time ``t`` and starting at an earlier
``t_{\textrm{prev}}``, which is the closes preceding event where the
photoelectron momentum changed, i.e. a previous scattering or the
initial photoionization.

"""
function recurse_common(::Type{Amp},
                        system::System{T}, 𝐤, 𝐩, iref, irange,
                        diagram::Diagram; weight::Function = 𝐪 -> true,
                        ϵ=1e-2*√(eps(T)), memory=system.ndt) where {Amp,T}
    α,interaction = get_interaction(system, diagram)

    Iₚ = system.ionization_channels[α].E
    tref = system.t[iref]
    amplitude = complex(zero(Amp))

    𝐀 = Base.Fix1(vector_potential, system.F)

    for i in irange
        tᵢ = system.t[i]
        τ = (tref-tᵢ) # Excursion time

        # At the time ti, the interaction takes us from ion
        # state β to α; we therefore need to accumulate phase in the α
        # channel in the interval ti..system.t[iref], and in
        # the β channel in the interval
        # system.t[irange[1]]..ti.
        Sᵢₒₙ = Iₚ*τ
        aᵢₒₙ = exp(-im*Sᵢₒₙ)

        # If 𝐩 is a determinate momentum, recurse will just return the
        # same value, i.e. this is the one we use to propagate the
        # Volkov waves. If however 𝐩 is an intermediate momentum, we
        # need to find its starting time, which we do by passing it
        # along recursively until we hit either an interaction that
        # changes the photoelectron momentum, or the initial time of
        # ionization.
        sub_amplitude,𝐩 = recurse(Amp, system, 𝐩, i, max(1,i-memory):i-1,
                                  diagram[2:end],
                                  ϵ=ϵ, memory=memory)

        # 𝐩 is now either a determinate momentum, or an indeterminate
        # momentum, but with two times: a starting time and a stopping
        # time, which we can use to find the saddle-point momentum.
        𝐩ₛₜ = momentum(𝐩)
        𝐤ᵢ = momentum(fix_momentum(𝐤, system, i))

        Sₑₗ = volkov_phase(𝐩ₛₜ, system.volkov, iref, i)
        𝐀ᵢ = 𝐀(tᵢ)
        a = weight(𝐩ₛₜ)*prefactor(𝐩)*aᵢₒₙ*exp(-im*Sₑₗ)*interaction(kinematic_momentum(𝐤ᵢ, 𝐀ᵢ), kinematic_momentum(𝐩ₛₜ, 𝐀ᵢ), i)

        amplitude += a*sub_amplitude
    end

    -im*system.dt*amplitude, 𝐤
end

# Canonical momentum is preserved, e.g. direct ionization in ATI or
# dipole interaction between ion states.
recurse(::CanonicalMomentumConservation,
        ::Type{Amp}, system::System, 𝐤, args...; kwargs...) where Amp =
    recurse_common(Amp, system, 𝐤, 𝐤, args...; kwargs...)

# Canonical momentum is not preserved in the interaction,
# e.g. scattering of the nucleus/remaining electrons. We therefore
# have to find the stationary momentum that leads back to the parent
# ion during the time interval i..iref.
function recurse(::NoCanonicalMomentumConservation,
                 ::Type{Amp}, system::System, 𝐤, iref, args...; kwargs...) where Amp
    𝐩 = IntermediateMomentum(nothing, iref, system)
    recurse_common(Amp, system, 𝐤, 𝐩, iref, args...; kwargs...)
end

function recurse(::Type{Amp}, system::System, 𝐤, iref, irange,
                 diagram::Diagram; kwargs...) where Amp
    # The base case corresponds to the ionization time, and hence we
    # fix the (possibly intermediate) momentum 𝐤 to have its starting
    # time at system.t[iref].
    isempty(diagram) && return true, fix_momentum(𝐤, system, iref)

    α, which = first(diagram)
    recurse(canonical_momentum_conservation(system, which),
            Amp, system, 𝐤, iref, irange,
            diagram; kwargs...)
end

# * High-level interface

function photoelectron_spectrum(k::AbstractArray{T},
                                system::System, diagram::Diagram;
                                verbosity=1, kwargs...) where T
    verbosity > 0 && @info "Photoelectrum spectrum calculation" system diagram length(k)

    cT = complex(eltype(T))
    c = similar(k, cT)
    iref = length(system.t)
    irange = 1:iref
    threaded_range_loop(eachindex(k)) do i
        c[i] = first(recurse(cT, system, k[i], iref, irange, diagram; kwargs...))
    end
    c
end

function photoelectron_spectrum(k, args...; kwargs...)
    system = System(args...; kwargs...)
    photoelectron_spectrum(k, system, Diagram(system); kwargs...)
end

"""
    induced_dipole(system, diagram[; kwargs...])

Compute the induced dipole moment as a function of time of the
[`System`](@ref), for a specific [`Diagram`](@ref).
"""
function induced_dipole(system::System, diagram::Diagram; verbosity = 1, kwargs...)
    if verbosity > 0
        display(system)
        println()
        display(diagram)
    end

    F = system.F
    t = system.t

    DT = typeof(field_amplitude(F, first(t)))
    𝐝 = zeros(DT, length(t))

    α = first(first(diagram))

    d = system.ionization_channels[α].st.d
    𝐀 = Base.Fix1(vector_potential, F)
    recombination = t -> (𝐩 -> conj(d(kinematic_momentum(𝐩, 𝐀(t)))))

    memory = get(kwargs, :memory, system.ndt)

    @showprogress for (i,t) in enumerate(t)
        𝐩 = IntermediateMomentum(nothing, i, system)
        𝐝̃,𝐩 = recurse(DT, system, 𝐩, i, max(1,i-memory):i-1,
                      diagram; weight=recombination(t), kwargs...)
        𝐝[i] = 2real(𝐝̃)
    end

    (system=system, diagram=diagram, dipole=𝐝)
end

function induced_dipole(args...; kwargs...)
    system = System(args...; kwargs...)
    induced_dipole(system, Diagram(system); kwargs...)
end

# * Exports

export Diagram, photoelectron_spectrum, induced_dipole
