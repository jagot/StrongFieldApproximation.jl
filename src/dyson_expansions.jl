# * System

"""
    System

Describes the combined system of atom/molecule (in terms of ionization
channels) and an external, time-dependent field. The ionization
channels may be coupled through various interactions, which may or may
not affect the photoelectron. `System` only describes the channels,
the external field, and the interactions possible, the actual process
is described by a [`Diagram`](@ref).
"""
struct System{T,IonizationChannels<:AbstractVector{<:IonizationChannel{T}},
              IonDipoleCouplings,
              Couplings<:AbstractVector{<:AbstractMatrix{<:AbstractCoupling}},
              VectorPotential,
              Times<:AbstractRange,
              Volkov<:VolkovPhases}
    ionization_channels::IonizationChannels
    𝐫ᵢₒₙ::IonDipoleCouplings

    couplings::Couplings

    𝐀::VectorPotential
    t::Times
    dt::T

    volkov::Volkov
end

IonDipoleCouplingsType = Union{AbstractMatrix{<:Number},UniformScaling,SVector{3},Nothing}
NoCouplings = Matrix{AbstractCoupling}[]

"""
    System(ionization_channels, 𝐫ᵢₒₙ, couplings, F, ndt)

Set up a [`System`](@ref) consisting of multiple
[`IonizationChannel`](@ref)s with ionic dipole moments `𝐫ᵢₒₙ` and
possible `couplings` between them, driven by an external field `F`,
sampled at a frequency of `fs`.
"""
function System(ionization_channels::AbstractVector{<:IonizationChannel},
                𝐫ᵢₒₙ::IonDipoleCouplingsType,
                couplings::AbstractVector,
                F::ElectricFields.AbstractField, fs::Number)
    t = timeaxis(F, fs)
    volkov = VolkovPhases(F, t)

    𝐀 = vector_potential.(F, t)
    System(ionization_channels, 𝐫ᵢₒₙ, couplings, 𝐀, t, step(t), volkov)
end

@doc raw"""
    System(ionization_channels, 𝐫ᵢₒₙ, couplings, Fv, Av, t)

Set up a [`System`](@ref) consisting of multiple
[`IonizationChannel`](@ref)s with with ionic dipole moments `𝐫ᵢₒₙ` and
possible `couplings` between them, driven by an external field `Fv`
with corresponding vector potential `Av`, both resolved on the times
given by `t`; it is up to the user to ensure that ``\vec{F} =
-\partial_t\vec{A}`` holds.
"""
function System(ionization_channels::AbstractVector{<:IonizationChannel},
                𝐫ᵢₒₙ::IonDipoleCouplingsType, couplings::AbstractVector,
                Fv::AbstractVector, Av::AbstractVector, t::AbstractRange)
    volkov = VolkovPhases(Av, t)

    System(ionization_channels, 𝐫ᵢₒₙ, couplings, Av, t, step(t), volkov)
end

System(ionization_channels::AbstractVector{<:IonizationChannel},
       F::ElectricFields.AbstractField, fs::Number) =
           System(ionization_channels, nothing, NoCouplings,
                  F, fs)

System(Iₚ::Number, args...; kwargs...) =
    System([IonizationChannel(Iₚ, args...)], args...; kwargs...)

num_channels(system::System) = length(system.ionization_channels)

Base.show(io::IO, system::System) =
    write(io, "$(num_channels(system))-channel SFA System")

function Base.show(io::IO, mime::MIME"text/plain", system::System)
    println(io, system, ":")
    for (i,c) in enumerate(system.ionization_channels)
        println(io, " ", i, ". ", c)
    end
    if !isnothing(system.𝐫ᵢₒₙ)
        nz(A) = count(!iszero, A)
        nz(A,i) = count(e -> !iszero(e[i]), A)
        nz(I::UniformScaling) = iszero(I) ? 0 : 1
        println(io, "Channel dipole couplings:")
        if system.𝐫ᵢₒₙ isa SVector{3}
            for (i,d) in enumerate(("x","y","z"))
                nzd = nz(system.𝐫ᵢₒₙ[i])
                println(io, "  - $d: $(nzd) non-zero couplings")
            end
        else
            nzz = nz(system.𝐫ᵢₒₙ)
            println(io, "  - z: $(nzz) non-zero couplings")
        end
    end
end

canonical_momentum_conservation(system::System, which::Integer) =
    iszero(which) ?
    CanonicalMomentumConservation() :
    canonical_momentum_conservation(first(system.couplings[which]))

kinematic_momentum(k::Number, A::Number) = k + A
kinematic_momentum(k::SVector{3}, A::SVector{3}) = k + A
kinematic_momentum(k::SVector{3}, A::Number) = SVector{3}(k[1], k[2], k[3]+A)

# * Diagrams

"""
    Diagram(path, couplings)

Goldstone diagram describing a strong-field process, i.e. after an
initial photoionization event, the ion and photoelectron are
propagated separately (at the chosen level of approximation, typically
Volkov waves are used for the photoelectrons). By interacting with the
external field, the ion state may change, and through electron
(re)scattering, both photoelectron and ion state may change.

The diagram is specified using a `path` of pairs (in antichronological
order), where each pair designates `(resultant ion channel,
interaction)`, where `interaction` is an index to one of the
`couplings`. The special value `interaction == 0` corresponds to the
initial ionization event, and should appear once at the end of `path`
(i.e. the first event); if it also appears as the first element of
`path` (i.e. the last event), it corresponds to recombination to the
initial state.

# Examples

```jldoctest
julia> couplings = (StrongFieldApproximation.DipoleCoupling,
                    StrongFieldApproximation.CoulombCoupling)
(StrongFieldApproximation.DipoleCoupling, StrongFieldApproximation.CoulombCoupling)

julia> Diagram([(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐤


julia> Diagram([(1,0),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ╲   ╱⇝
   |0⟩

julia> Diagram([(2,1),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
 ⇝┃   │
 2┃   │𝐤


julia> Diagram([(2,2),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ┠┈┈┈┤
 2┃   │𝐤


julia> Diagram([(1,2),(2,2),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ┠┈┈┈┤
 2┃   │𝐪
  ┠┈┈┈┤
 1┃   │𝐤


julia> Diagram([(3,1),(1,1),(2,2),(1,0)], couplings)
Goldstone Diagram:
   |0⟩
  ╱   ╲⇜
 1┃   │𝐩
  ┠┈┈┈┤
 2┃   │𝐪
 ⇝┃   │
 1┃   │
 ⇝┃   │
 3┃   │𝐤
```
"""
struct Diagram{Couplings}
    path::Vector{Tuple{Int,Int}}
    couplings::Couplings
    function Diagram(path, couplings::Couplings) where Couplings
        if !isempty(path)
            α,which = last(path)
            which == 0 ||
                throw(ArgumentError("Non-empty diagrams must have photoionization as the first interaction"))
        end
        ncoup = length(couplings)
        vcoup = 0:ncoup
        invalid_couplings = findall(∉(vcoup), last.(path))
        if !isempty(invalid_couplings)
            length(invalid_couplings) == 1 &&
                throw(ArgumentError("Coupling at vertex $(first(invalid_couplings)) not in valid range $(vcoup)"))
            throw(ArgumentError("Couplings at vertices $(invalid_couplings) not in valid range $(vcoup)"))
        end
        for i = 2:length(path)-1
            path[i][2] == 0 &&
                throw(ArgumentError("Only the first and last interaction in a diagram may have which = 0 (corresponding to ionization/recombination)"))
        end
        if length(path) > 1 && path[1][2] == 0
            path[1][1] == path[2][1] ||
                throw(ArgumentError("Must recombine from the last channel"))
        end
        new{Couplings}(path, couplings)
    end
end

Diagram(path::Tuple{Int,Int}, couplings) =
    Diagram([path], couplings)

"""
    Diagram(path::AbstractVector{Tuple{Int,Int}}, system)

Convenience constructor that sets up the [`Diagram`](@ref) given a
`path` and the couplings present in `system`. If `path` is empty, the
diagram corresponding to photoionization into the first ionization of
channel of `system` will automatically be generated.
"""
function Diagram(path::AbstractVector, system::System)
    nc = length(system.ionization_channels)
    if isempty(path)
        nc > 1 &&
            @warn "More than one ionization channel present in system, choosing the first"
        path = [(1,0)]
    end
    vc = 1:nc
    invalid_channels = findall(∉(vc), first.(path))
    if !isempty(invalid_channels)
        length(invalid_channels) == 1 &&
            throw(ArgumentError("Invalid channel at vertex $(first(invalid_channels))"))
        throw(ArgumentError("Invalid channels at vertices $(invalid_channels)"))
    end
    Diagram(path, [typeof(first(c)) for c in system.couplings])
end

Diagram(system::System) = Diagram([], system)

for f in [:length, :first, :firstindex, :lastindex, :isempty]
    @eval Base.$f(d::Diagram) = $f(d.path)
end
Base.getindex(d::Diagram, i) = Diagram(d.path[i], d.couplings)

function draw_ionization(io, nd)
    println(io, lpad("|0⟩", nd+4))
    print(io, lpad("╱ ╲", nd+4))
    printstyled(io, "⇜", color=:light_red)
    println(io)
end

function draw_recombination(io, nd)
    print(io, lpad("╲ ╱", nd+4))
    printstyled(io, "⇝", color=:light_red)
    println(io)
    println(io, lpad("|0⟩", nd+4))
end

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
    electrons = ["𝐩","𝐪"]
    cur_electron = 1
    for (i,(ion,which)) in Iterators.reverse(enumerate(d.path))
        if which == 0
            if i == length(d.path)
                draw_ionization(io, nd)
            elseif i == 1
                draw_recombination(io, nd)
            end
        else
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
            ei = mod1(cur_electron,length(electrons))
            e = electrons[ei]*repeat("′", fld1(cur_electron,length(electrons))-1)
            cur_electron += 1
            e
        end

        which == 0 && i == 1 && length(d) > 1 || draw_exciton(io, ion, nd, electron)
    end
end

# * Integrate diagrams

function analyze_diagram(system, diagram)
    α,which = first(diagram)
    ions = Int[]
    # For photoelectron spectra, time 1 is the reference time,
    # typically at which the laser pulse has ended; for direct
    # photoelectron diagrams, time 2 is the time of ionization. For
    # dipoles, time 1 is the time of recombination.
    unique_momenta = Tuple{Int,Int}[(1,2)]
    momenta = [1]
    indeterminate_momenta = Int[]
    ld = length(diagram)
    order = ld

    if ld > 1
        if which == 0
            push!(indeterminate_momenta, 1)
            # Recombination does not increase the order of the diagram, since
            # it only amounts to projecting the wavefunction on ⟨0|𝐫.
            order -= 1
        end
    else
        push!(ions, α)
    end

    i = 2
    for (α,which) in diagram.path[(iszero(which) ? 2 : 1):end]
        push!(ions, α)
        a,b = unique_momenta[end]
        unique_momenta[end] = (a,i)
        if !(iszero(which) || canonical_momentum_conservation(system, which) == CanonicalMomentumConservation())
            push!(unique_momenta, (i,i+1))
            push!(indeterminate_momenta, length(unique_momenta))
        end
        if !iszero(which)
            push!(momenta, length(unique_momenta))
        end

        i += 1
    end

    return ions, unique_momenta, momenta, indeterminate_momenta, order
end

# 𝐤 is nothing, we want a dipole amplitude
function diagram_amplitude(::Type{Amp}, system::System, diagram::Diagram, iref, i, ::Nothing) where Amp
    amplitude = complex(zero(Amp))

    amplitude
end

# 𝐤 is determinate, we want a photoelectron spectrum
function diagram_amplitude(::Type{Amp}, system::System, diagram::Diagram, iref, i, 𝐤) where Amp
    amplitude = complex(zero(Amp))

    amplitude
end

momentum_type(_, 𝐤) = typeof(𝐤)
momentum_type(system, ::Nothing) = eltype(system.volkov.∫A)

set_momentum!(𝐩s::AbstractVector{<:Number}, 𝐩ₛₜ::Number, i) =
    setindex!(𝐩s, 𝐩ₛₜ, i)
set_momentum!(𝐩s::AbstractVector{<:SVector{3}}, 𝐩ₛₜ::SVector{3}, i) =
    setindex!(𝐩s, 𝐩ₛₜ, i)
set_momentum!(𝐩s::AbstractVector{<:SVector{3}}, 𝐩ₛₜ::T, i) where {T<:Number} =
    setindex!(𝐩s, SVector{3,T}(zero(T), zero(T), 𝐩ₛₜ), i)

function evaluate_momenta!(𝐩s, prefactors, system, unique_momenta, indeterminate_momenta, i;
                           ϵ=1e-2*√(eps(eltype(system.t))))
    for idm in indeterminate_momenta
        uidm = unique_momenta[idm]
        a,b = i[uidm[1]],i[uidm[2]]
        set_momentum!(𝐩s, stationary_momentum(system.volkov, a, b), idm)
        τ = system.t[a]-system.t[b]
        ζ = (2π/(im*τ + ϵ))^(3/2)
        prefactors[idm] = ζ
    end
end

function ionization(system::System, diagram::Diagram, 𝐩, 𝐀, i)
    α,which = diagram.path[end]
    @assert which == 0
    source_term(system.ionization_channels[α],
                i,
                kinematic_momentum(𝐩, 𝐀[i]))
end

function recombination(system::System, diagram::Diagram, 𝐩, 𝐀, i)
    α,which = first(diagram)
    if which == 0 && length(diagram) > 1
        d = system.ionization_channels[α].st.d
        conj(d(kinematic_momentum(𝐩, 𝐀[i])))
    else
        true
    end
end

function integrate_diagram(::Type{Amp}, system::System, diagram::Diagram, iref, 𝐤=nothing; memory=typemax(Int), imin=1,
                           to=TimerOutput(), verbosity=1) where Amp
    ions, unique_momenta, momenta, indeterminate_momenta, order = @timeit to "Analyze diagram" analyze_diagram(system, diagram)

    @timeit to "Allocations" begin
        weight = (-im*system.dt)^order

        verbosity > 1 && @info "Integrating diagram up to" iref system diagram ions unique_momenta momenta indeterminate_momenta order weight 𝐤

        Eᵢₒₙₛ = [system.ionization_channels[ion].E for ion in ions]
        𝐩s = complex(zeros(momentum_type(system, 𝐤), length(unique_momenta)))
        prefactors = ones(complex(eltype(system.t)), length(unique_momenta))
        if !isnothing(𝐤)
            𝐩s[1] = 𝐤
        end

        𝐀 = system.𝐀
        amplitude = complex(zero(Amp))
        ctT = complex(eltype(system.t))

        is = vcat(iref, zeros(Int, order))
    end

    @timeit to "Time loop" begin
        for i in TelescopeIterator(max(1,imin):iref-1, order, memory)
            is[2:end] .= i
            # is = vcat(iref, i)
            # is = (iref,i...)
            # for i in 1:iref-1, j in 1:i-1
            #     is = (iref,i,j)
            # println(is)

            @timeit to "Evaluate momenta" evaluate_momenta!(𝐩s, prefactors, system, unique_momenta, indeterminate_momenta, is)
            verbosity > 10 && @show 𝐩s prefactors

            @timeit to "Evaluate propagators" begin
                Sᵢₒₙ = zero(ctT)
                Sₑₗ = zero(ctT)
                for j = 1:order
                    a,b = is[j], is[j+1]
                    τ = system.t[a] - system.t[b]
                    Sᵢₒₙ += Eᵢₒₙₛ[j]*τ
                    Sₑₗ += volkov_phase(𝐩s[j], system.volkov, a, b)
                end
                S = Sᵢₒₙ + Sₑₗ
                aₚᵣₒₚ = prod(prefactors)*exp(-im*S)
            end

            verbosity > 10 && @show is prefactors

            ∂a = @timeit to "Prefactor" (ionization(system, diagram, 𝐩s[end], 𝐀, is[end]) *
                                         aₚᵣₒₚ *
                                         recombination(system, diagram, 𝐩s[1], 𝐀, iref))

            @timeit to "Interior interactions" begin
                # Loop over "interior" interactions
                for j = (order>1 && first(diagram)[2]==0 ? 3 : 2):order
                    ion,which = diagram.path[j-1]
                    α = ions[j-1]
                    β = ions[j]
                    verbosity > 20 && @show j, ion, which α,β
                    interaction = system.couplings[which][α,β]

                    𝐤ᵢ = 𝐩s[momenta[j-1]]
                    𝐩ᵢ = 𝐩s[momenta[j]]
                    𝐀ᵢ = 𝐀[is[j]]

                    ∂a *= @timeit to "Interaction" interaction(kinematic_momentum(𝐤ᵢ, 𝐀ᵢ), kinematic_momentum(𝐩ᵢ, 𝐀ᵢ), is[j+1])
                end
            end

            amplitude += ∂a
        end
    end

    @timeit to "Weighting" weight*amplitude
end

# * High-level interface

function photoelectron_spectrum(k::AbstractArray{T},
                                system::System, diagram::Diagram;
                                iref=length(system.t),
                                verbosity=1, kwargs...) where T
    verbosity > 0 && @info "Photoelectrum spectrum calculation" system diagram length(k)

    cT = complex(eltype(T))
    c = similar(k, cT)
    to = TimerOutput()
    p = Progress(length(k))
    @timeit to "k loop" begin
        threaded_range_loop(eachindex(k)) do i
            tok = TimerOutput()
            c[i] = integrate_diagram(cT, system, diagram, iref, k[i]; to=tok, verbosity=verbosity-1, kwargs...)
            ProgressMeter.next!(p)
            merge!(to, tok, tree_point=["k loop"])
        end
    end
    TimerOutputs.complement!(to)
    verbosity > 0 && print_timer(to)
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
    verbosity > 0 && @info "Induced dipole calculation" system diagram

    t = system.t

    DT = eltype(system.𝐀)
    𝐝 = zeros(DT, length(t))

    memory = get(kwargs, :memory, typemax(Int))

    @showprogress for (i,t) in enumerate(t)
        𝐝̃ = integrate_diagram(DT, system, diagram, i; imin=i-memory, verbosity=verbosity-1, kwargs...)
        𝐝[i] = 2real(𝐝̃)
    end

    𝐝
end

function induced_dipole(args...; kwargs...)
    system = System(args...; kwargs...)
    diagram = Diagram([(1,0),(1,0)], system)
    induced_dipole(system, diagram; kwargs...)
end

# * Exports

export Diagram, photoelectron_spectrum, induced_dipole
