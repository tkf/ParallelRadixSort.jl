struct Pool{F,T}
    f::F
    capacity::Int
    buffers::Vector{Vector{T}}
end

Pool(f, capacity = 256) = Pool(f, capacity, empty!([[f()]]))

function initpool!(pool::Pool)
    buffers = pool.buffers
    resize!(buffers, Threads.nthreads())
    for i in eachindex(buffers)
        buffers[i] = empty!(eltype(buffers)(undef, pool.capacity))
    end
    return pool
end

function alloc!(pool::Pool, n::Integer)
    buffer = pool.buffers[Threads.threadid()]
    i = max(1, length(buffer) - n + 1)
    objects = buffer[i:end]
    resize!(buffer, i - 1)

    m = length(objects)
    if m < n
        resize!(objects, n)
        for k in m+1:n
            objects[k] = pool.f()
        end
    end

    return objects
end

function free!(pool::Pool{<:Any,T}, objects::AbstractVector{T}) where {T}
    buffer = pool.buffers[Threads.threadid()]
    append!(buffer, @view objects[1:min(end, pool.capacity - length(buffer))])
    return
end
