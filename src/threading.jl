struct RangeIterator
    k::Int
    d::Int
    r::Int
end

"""
    RangeIterator(n::Int,k::Int)

Returns an iterator splitting the range `1:n` into `min(k,n)` parts of (almost) equal size.
"""
RangeIterator(n::Int, k::Int) = RangeIterator(min(n,k),divrem(n,k)...)
Base.length(it::RangeIterator) = it.k
endpos(it::RangeIterator, i::Int) = i*it.d+min(i,it.r)
Base.iterate(it::RangeIterator, i::Int=1) = i>it.k ? nothing : (endpos(it,i-1)+1:endpos(it,i), i+1)

function threaded_range_loop(fun::Function, n::Integer)
    @sync for r in RangeIterator(n, Threads.nthreads())
        Threads.@spawn for i in r
            fun(i)
        end
    end
end

threaded_range_loop(fun::Function, v::AbstractArray) =
    threaded_range_loop(i -> fun(v[i]), length(v))

threaded_enumerate_range_loop(fun::Function, v::AbstractArray) =
    threaded_range_loop(i -> fun(i,v[i]), length(v))
