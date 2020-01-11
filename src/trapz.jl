module trapz

# Below routines borrowed from https://github.com/cmcbride/julia/blob/master/base/integrate.jl

trapzv(y::AbstractVector) =
    trapz(y, one(eltype(y)))

# integration for uniform points (in x)
function trapzv(y::AbstractVector, dx::Number)
    if length(y) == 0
        return 0
    end
    r = zero(zero(eltype(y)) * zero(dx))
    r += (y[1] + y[end])/2
    r += sum(y[2:end-1])
    r * dx
end

# integration for non-uniform points (in x)
function trapzv(y::AbstractVector,x::AbstractVector)
    n = length(y)
    if n != length(x)
        error("Input x,y must be of same length")
    end
    r = zero(zero(eltype(x))*zero(eltype(y)))
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0
end

# export Trapz, next!, reset!, trapzv
export trapzv
end
