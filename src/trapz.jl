module trapz
type Trapz
    f::Function
    x::AbstractVector
    i::Integer
    cur_integral
    cur_val
    I
end

Trapz(f,x) = Trapz(f,x,0,f(x[1]),0,0)

#=
Composite trapezoidal rule:
\[
\int\limits_a^b f(x) dx \approx
\frac{b-a}{2N} \left[f(x_1)+f(x_{N+1})+2\sum_{i=1}^N f(x_i)\right]
\]
=#

function next!(t::Trapz)
    if t.i == 0
        t.i = 1
        return 0
    elseif t.i < length(t.x)
        t.cur_integral += 2t.cur_val
        t.i += 1
        x = t.x[t.i]
        t.cur_val = t.f(x)
        t.I = ((x - t.x[1])/2(length(1:t.i)-1))*(t.cur_integral + t.cur_val)
    else
        t.I
    end
end

function reset!(t::Trapz, x::AbstractVector)
    t.x = x
    t.i = 0
    t.cur_integral = t.f(x[1])
    t.cur_val = 0
    t.I = 0
end

# Below routines borrowed from https://github.com/cmcbride/julia/blob/master/base/integrate.jl

function trapzv(y::AbstractVector)
    trapz(y, one(eltype(y)))
end

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

export Trapz, next!, reset!, trapzv
end
