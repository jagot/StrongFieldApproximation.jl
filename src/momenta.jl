kinematic_momentum(k::Number, A::Number) = k + A
kinematic_momentum(k::SVector{3}, A::SVector{3}) = k + A
kinematic_momentum(k::SVector{3}, A::Number) = SVector{3}(k[1], k[2], k[3]+A)

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
