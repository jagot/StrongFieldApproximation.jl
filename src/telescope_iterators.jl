"""
    TelescopeIterator(v, n[, memory=typemax])

Iterator that is formally equivalent to
```julia
for i = 1:length(v), j = max(1,i-memory):i-1, k = max(1,j-memory):j-1, ...
    e = v[i,j,k,...]
end
```
for `n` loop variables.

# Example

```jldoctest
julia> ti = StrongFieldApproximation.TelescopeIterator(1:6, 2, 3)
StrongFieldApproximation.TelescopeIterator{UnitRange{Int64}, Int64}(1:6, 2, 3)

julia> for e in ti
       println(e)
       end
[2, 1]
[3, 1]
[4, 1]
[3, 2]
[4, 2]
[5, 2]
[4, 3]
[5, 3]
[6, 3]
[5, 4]
[6, 4]
[6, 5]
```
"""
struct TelescopeIterator{V<:AbstractVector, I<:Integer}
    v::V
    n::I
    memory::I
    function TelescopeIterator(v::V, n::I, memory::I=typemax(I)) where {V<:AbstractVector, I<:Integer}
        (n < 0 || memory < 0) && throw(ArgumentError("Invalid values for n = $n or memory = $memory"))
        new{V,I}(v, n, memory)
    end
end

function Base.iterate(ti::TelescopeIterator, state=[ti.n+1-i for i=1:ti.n])
    (isempty(ti.v) || iszero(ti.n) || ti.n > length(ti.v) || ti.n > 1 && iszero(ti.memory)) && return nothing

    nv = length(ti.v)
    any(>(nv), state) && return nothing

    val = ti.v[state]

    if ti.n == 1 && state[1] < nv
        state[1] += 1
        return val, state
    end

    for j = 1:ti.n
        state[j] > nv || j == ti.n && state[j] == nv-(ti.n-2) && return nothing
        nsj = state[j] + 1
        j < ti.n && nsj - state[j+1] ≤ ti.memory && nsj ≤ nv-(j-1) || j == ti.n || continue
        state[j] = nsj
        for i = j-1:-1:1
            state[i] = state[i+1]+1
        end
        break
    end

    val, state
end
