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
              Ions<:IonPropagator,
              Volkov<:VolkovPhases}
    ionization_channels::IonizationChannels
    𝐫ᵢₒₙ::IonDipoleCouplings

    couplings::Couplings

    𝐀::VectorPotential
    t::Times
    dt::T

    ions::Ions
    volkov::Volkov
end

IonDipoleCouplingsType = Union{AbstractMatrix{<:Number},UniformScaling,SVector{3},Nothing}
NoCouplings = Matrix{AbstractCoupling}[]

function System(ionization_channels, 𝐫ᵢₒₙ, couplings, 𝐅, 𝐀, t, volkov::VolkovPhases;
                Ions::Type{<:IonPropagator}=LaserDressedIons, kwargs...)
    ions = Ions(ionization_channels, 𝐫ᵢₒₙ, 𝐅, t)
    System(ionization_channels, 𝐫ᵢₒₙ, couplings, 𝐀, t, step(t), ions, volkov)
end

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
                F::ElectricFields.AbstractField, fs::Number; kwargs...)
    t = timeaxis(F, fs)
    volkov = VolkovPhases(F, t)

    𝐀 = vector_potential.(F, t)
    System(ionization_channels, 𝐫ᵢₒₙ, couplings, F, 𝐀, t, volkov; kwargs...)
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
                Fv::AbstractVector, Av::AbstractVector, t::AbstractRange;
                kwargs...)
    volkov = VolkovPhases(Av, t)

    System(ionization_channels, 𝐫ᵢₒₙ, couplings, Fv, Av, t, volkov; kwargs...)
end

System(ionization_channels::AbstractVector{<:IonizationChannel},
       F::ElectricFields.AbstractField, fs::Number; kwargs...) =
           System(ionization_channels, nothing, NoCouplings,
                  F, fs; kwargs...)

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
