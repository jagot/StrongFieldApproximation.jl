function momentum_grid(kmin, kmax, nk, nθ; spacing=:momentum)
    kmag = if spacing == :momentum
        range(kmin, stop=kmax, length=nk)
    elseif spacing == :energy
        .√(2range(kmin^2/2, stop=kmax^2/2, length=nk))
    else
        throw(ArgumentError("Unknown spacing $(spacing)"))
    end
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

export momentum_grid
