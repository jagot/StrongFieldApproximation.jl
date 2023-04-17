function momentum_grid(kmin, kmax, nk; spacing=:momentum)
    if spacing == :momentum
        range(kmin, stop=kmax, length=nk)
    elseif spacing == :energy
        .√(2range(kmin^2/2, stop=kmax^2/2, length=nk))
    else
        throw(ArgumentError("Unknown spacing $(spacing)"))
    end
end

function momentum_grid(kmin, kmax, nk, nθ, nϕ=1; kwargs...)
    kmag = momentum_grid(kmin, kmax, nk; kwargs...)
    nθ == 1 && return (kmag,kmag,nothing)

    x = range(0, stop=1, length=nθ)
    k = Array{SVector{3,Float64}}(undef, length(kmag), nθ)
    for j = 1:nθ
        s,c = sincospi(x[j])
        k̂ = SVector{3}(s, 0, c)
        for (i,km) in enumerate(kmag)
            k[i,j] = km*k̂
        end
    end
    k,kmag,x*π
end

momentum_grid(Emin::Unitful.Energy, Emax::Unitful.Energy, args...; kwargs...) =
    momentum_grid(√(2austrip(Emin)), √(2austrip(Emax)), args...; kwargs...)

export momentum_grid
