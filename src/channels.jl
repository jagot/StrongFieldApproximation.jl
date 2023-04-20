# * Source terms

abstract type AbstractSourceTerm end

# ** Dipole

struct DipoleSourceTerm{Dipole<:Function,Field<:Union{ElectricFields.AbstractField,AbstractVector}} <: AbstractSourceTerm
    d::Dipole
    F::Field

    DipoleSourceTerm(d::Dipole, F::Field) where {Dipole<:Function, Field<:ElectricFields.AbstractField} =
        new{Dipole,Field}(d, F)

    function DipoleSourceTerm(d::Dipole, F::Field, fs::Number) where {Dipole<:Function, Field<:ElectricFields.AbstractField}
        t = timeaxis(F, fs)
        Fv = field_amplitude.(F, t)
        new{Dipole,typeof(Fv)}(d, Fv)
    end

    DipoleSourceTerm(d::Dipole, Fv::FieldAmplitude, Av::AbstractVector, t::AbstractRange) where {Dipole<:Function, FieldAmplitude<:AbstractVector} =
        new{Dipole,FieldAmplitude}(d, Fv)
end

Base.show(io::IO, d::DipoleSourceTerm) =
    write(io, "DipoleSourceTerm")

function Base.show(io::IO, mime::MIME"text/plain", d::DipoleSourceTerm)
    show(io, d)
    print(io, "\nDipole: ")
    show(io, mime, d.d)
    print(io, "\nField: ")
    show(io, d.F)
end

dipole(F::Number, d::Number) = F*d
dipole(F::SVector{3}, d::SVector{3}) = F[1]*d[1] + F[2]*d[2] + F[3]*d[3]
dipole(F::Number, d::SVector{3}) = F*d[3]

_field_amplitude(F::ElectricFields.AbstractField, t) =
    field_amplitude(F, t)
_field_amplitude(F::AbstractVector, i::Integer) = F[i]

source_term(d::DipoleSourceTerm, t, p) = dipole(_field_amplitude(d.F, t), d.d(p))

# * Ionizations channels

"""
    IonizationChannel(E, st)

Represents an ionization channel with energy `E` above the neutral
(which is taken to have energy `0 Ha`).

"""
struct IonizationChannel{T,SourceTerm<:AbstractSourceTerm}
    E::T
    st::SourceTerm
    IonizationChannel(E::T, st::SourceTerm) where {T,SourceTerm} =
        new{T,SourceTerm}(E, st)
end

source_term(ic::IonizationChannel, args...) =
    source_term(ic.st, args...)

IonizationChannel(E, args...) =
    IonizationChannel(austrip(E), DipoleSourceTerm(Base.Fix2(d_hyd, hyd_α(austrip(E))), args...))

Base.show(io::IO, ic::IonizationChannel) =
    write(io, "IonizationChannel: Iₚ = $(ic.E) Ha = $(27.211ic.E) eV")

function Base.show(io::IO, mime::MIME"text/plain", ic::IonizationChannel)
    show(io, ic)
    print(io, "\nSource term: ")
    show(io, mime, ic.st)
end

# * Couplings

abstract type AbstractCanonicalMomentumConservation end
struct CanonicalMomentumConservation <: AbstractCanonicalMomentumConservation end
struct NoCanonicalMomentumConservation <: AbstractCanonicalMomentumConservation end

abstract type AbstractCoupling end
Base.iszero(::AbstractCoupling) = false
canonical_momentum_conservation(::AbstractCoupling) = NoCanonicalMomentumConservation()
canonical_momentum_conservation(::Type{<:AbstractCoupling}) = NoCanonicalMomentumConservation()

struct NoCoupling <: AbstractCoupling end
Base.iszero(::NoCoupling) = true
Base.zero(::AbstractCoupling) = NoCoupling()
Base.zero(::Type{<:AbstractCoupling}) = NoCoupling()

Base.show(io::IO, ::NoCoupling) = write(io, "0")

# ** Dipole couplings

struct DipoleCoupling{DipoleMoment<:Union{Number,SVector{3}},Field<:ElectricFields.AbstractField} <: AbstractCoupling
    d::DipoleMoment
    F::Field
end

canonical_momentum_conservation(::DipoleCoupling) = CanonicalMomentumConservation()
canonical_momentum_conservation(::Type{<:DipoleCoupling}) = CanonicalMomentumConservation()

Base.show(io::IO, ::DipoleCoupling) = write(io, "𝐅⋅𝐫")

# ** Coulomb couplings

momentum_transfer(k::Number, p::Number) = k-p
momentum_transfer(k::SVector{3}, p::SVector{3}) = k-p
momentum_transfer(k::SVector{3}, p::Number) = SVector{3}(k[1], k[2], k[3]-p)

struct CoulombCoupling{Coupling} <: AbstractCoupling
    coupling::Coupling
end

(cc::CoulombCoupling)(𝐤, 𝐩, _) = cc.coupling(𝐤, 𝐩)

Base.show(io::IO, ::CoulombCoupling) = write(io, "ĝ")

function Base.show(io::IO, mime::MIME"text/plain", g::CoulombCoupling)
    write(io, "CoulombCoupling: ")
    show(io, mime, g.coupling)
end

# * Exports

export IonizationChannel
