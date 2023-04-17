# https://en.wikipedia.org/wiki/Breadth-first_search
function find_blocks(G::SparseMatrixCSC)
    m = size(G,1)
    visited = falses(m)

    rows = rowvals(G)
    vals = nonzeros(G)

    blocks = Vector{Vector{Int}}()
    while !all(visited)
        s = findfirst(!, visited)
        visited[s] = true
        nzs = nzrange(G, s)
        (isempty(nzs) || iszero(vals[first(nzs)])) && continue
        block = [s]

        q = Queue{Int}()
        enqueue!(q, s)
        while !isempty(q)
            v = dequeue!(q)
            for k in nzrange(G, v)
                w = rows[k]
                (visited[w] || iszero(vals[k])) && continue
                enqueue!(q, w)
                visited[w] = true
                push!(block, w)
            end
        end
        push!(blocks, sort(block))
    end

    blocks
end
