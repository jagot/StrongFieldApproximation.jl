# * Source terms

abstract type AbstractSourceTerm end

# ** Dipole

struct DipoleSourceTerm{Dipole<:Function,Field<:ElectricFields.AbstractField} <: AbstractSourceTerm
    d::Dipole
    Fv::Field
end

dipole(F::Number, d::Number) = F*d
dipole(F::SVector{3}, d::SVector{3}) = dot(F, d)
dipole(F::Number, d::SVector{3}) = F*d[3]

source_term(d::DipoleSourceTerm, t, p) =
    dipole(field_amplitude(d.Fv, t), d.d(p))

# * Ionizations channels

"""
    IonizationChannel(E, st)

Represents an ionization channel with energy `E` above the neutral
(which is taken to have energy `0 Ha`).

"""
struct IonizationChannel{T,SourceTerm<:AbstractSourceTerm}
    E::T
    st::SourceTerm
end

source_term(ic::IonizationChannel, args...) =
    source_term(ic.st, args...)

IonizationChannel(E, F::ElectricFields.AbstractField) =
    IonizationChannel(austrip(E), DipoleSourceTerm(Base.Fix2(d_hyd, hyd_α(austrip(E))), F))

Base.show(io::IO, ic::IonizationChannel) =
    write(io, "IonizationChannel: Iₚ = $(ic.E) Ha = $(27.211ic.E) eV")

# * Couplings

abstract type AbstractCanonicalMomentumConservation end
struct CanonicalMomentumConservation <:  AbstractCanonicalMomentumConservation end
struct NoCanonicalMomentumConservation <:  AbstractCanonicalMomentumConservation end

abstract type AbstractCoupling end
Base.iszero(::AbstractCoupling) = false
canonical_momentum_conservation(::AbstractCoupling) = NoCanonicalMomentumConservation()

struct NoCoupling <: AbstractCoupling end
Base.iszero(::NoCoupling) = true
Base.zero(::AbstractCoupling) = NoCoupling()

# ** Dipole couplings

struct DipoleCoupling{DipoleMoment<:Union{Number,SVector{3}},Field<:ElectricFields.AbstractField} <: AbstractCoupling
    d::DipoleMoment
    F::Field
end

canonical_momentum_conservation(::DipoleCoupling) = CanonicalMomentumConservation()

# ** Coulomb couplings

momentum_transfer(k::Number, p::Number) = k-p
momentum_transfer(k::SVector{3}, p::SVector{3}) = k-p
momentum_transfer(k::SVector{3}, p::Number) = SVector{3}(k[1], k[2], k[3]-p)

struct CoulombCoupling{Coupling} <: AbstractCoupling
    coupling::Coupling
end
