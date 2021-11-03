_default_basesize(n::Integer) = max(1, cld(n, Threads.nthreads()))
default_basesize(xs) = _default_basesize(length(xs))
const DEFAULT_SMALLSIZE = 64

Base.@kwdef struct ParallelStableMSDRadixSortAlg{Alg,BaseSize} <: RadixSortAlgorithm
    smallsort::Alg = Base.Sort.DEFAULT_STABLE
    smallsize::Int = DEFAULT_SMALLSIZE
    basesize::BaseSize = nothing  # lazily determined
end

Base.@kwdef struct StableMSDRadixSortAlg{Alg} <: RadixSortAlgorithm
    smallsort::Alg = Base.Sort.DEFAULT_STABLE
    smallsize::Int = DEFAULT_SMALLSIZE
end

const ParallelStableMSDRadixSort = ParallelStableMSDRadixSortAlg()
const StableMSDRadixSort = StableMSDRadixSortAlg()

const ParallelMSDRadixSort =
    ParallelStableMSDRadixSortAlg(smallsort = Base.Sort.DEFAULT_UNSTABLE)
const MSDRadixSort = StableMSDRadixSortAlg(smallsort = Base.Sort.DEFAULT_UNSTABLE)

_compose(::typeof(identity), ::typeof(identity)) = identity
_compose(::typeof(identity), f) = f
_compose(f, ::typeof(identity)) = f
_compose(::typeof(-), ::typeof(-)) = identity  # OK?
_compose(f, g) = f ∘ g

# TODO: reversing order should be done after the projection
_ord_by(o) = error("unsupported ordering: $o")
_ord_by(::ForwardOrdering) = identity
_ord_by(::ReverseOrdering{ForwardOrdering}) = -  # TODO; fix it for strings
_ord_by(order::Order.By) = _compose(_ord_by(order.order), order.by)

function Base.sort!(
    v::AbstractVector,
    lo::Integer,
    hi::Integer,
    a::ParallelStableMSDRadixSortAlg,
    o::Ordering,
)
    smallsort!(xs; by) = sort!(xs; alg = a.smallsort, by = by)
    xs = view(v, lo:hi)
    _stablemsd!(
        xs,
        copy(xs),
        _ord_by(o),
        static(1),
        something(a.basesize, _default_basesize(hi - lo + 1)),
        a.smallsize,
        smallsort!,
        Val(true),  # xs_mutable
    )
    return xs
end

function Base.sort!(
    v::AbstractVector,
    lo::Integer,
    hi::Integer,
    a::StableMSDRadixSortAlg,
    o::Ordering,
)
    smallsort!(xs; by) = sort!(xs; alg = a.smallsort, by = by)
    xs = view(v, lo:hi)
    _stablemsd_seq!(xs, copy(xs), _ord_by(o), static(1), a.smallsize, smallsort!, Val(true))
    return xs
end

# TODO: NaN ordering correctness
Base.Sort.Float.fpsort!(v::AbstractVector, a::RadixSortAlgorithm, o::Ordering) =
    sort!(v, firstindex(v), lastindex(v), a, o)

# https://cs.stackexchange.com/a/96420
# https://www.youtube.com/watch?v=zqs87a_7zxw

function _stablemsd!(
    ys,
    xs,
    by::BY,
    ibyte,
    basesize,
    smallsize,
    smallsort!,
    xs_mutable::Union{Val{true},Val{false}} = Val(false),
) where {BY}
    @assert firstindex(xs) == 1
    @assert firstindex(ys) == 1

    will_be_allpadded(eltype(xs), by, ibyte) && return ys

    if length(xs) <= basesize
        return _stablemsd_seq!(ys, xs, by, ibyte, smallsize, smallsort!, xs_mutable)
    end

    # Split the input into at most `nthreads` parts since there's no recursion
    # in `_countmsd!` to benefit from constructive cache sharing:
    chunksize = max(basesize, cld(length(xs), Threads.nthreads()))

    chunks = Iterators.partition(xs, chunksize)
    allcounts = [zeros(MVector{256,Int}) for _ in 1:length(chunks)]
    allpaddeds = Vector{Bool}(undef, length(chunks))
    let allcounts = allcounts
        @sync for (i, xs) in enumerate(chunks)
            @spawn (_, allpaddeds[i]) = _countmsd!(allcounts[i], xs, by, ibyte)
        end
    end
    alloffsets = allcounts
    @DBG allcounts = map(copy, alloffsets)
    all(allpaddeds) && return ys

    # TODO: parallel prefix
    # exclusive scan
    foldl(Iterators.product(alloffsets, 1:256); init = 0) do acc, (offsets, i)
        @inbounds (acc, offsets[i]) = (acc + offsets[i], acc)
        return acc
    end

    @sync for (idx, exclusive_offsets) in
              zip(Iterators.partition(1:length(xs), chunksize), alloffsets)
        @spawn _scatter!(ys, xs, idx, exclusive_offsets, by, ibyte)
    end
    @DBG @check vec(reduce(hcat, alloffsets)') == cumsum(vec(reduce(hcat, allcounts)'))

    _spawn_foreach_remaining_subrange(alloffsets[end]) do idx
        ys_chunk = view(ys, idx)
        if length(idx) <= smallsize
            smallsort!(ys_chunk, by = by)
        else
            xs_chunk = if xs_mutable === Val(true)
                view(xs, idx)
            else
                similar(ys_chunk)
            end
            _stablemsd!(
                ys_chunk,
                copyto!(xs_chunk, ys_chunk),
                by,
                static_if_leq(@stat(ibyte + 1), static(8)),
                basesize,
                smallsize,
                smallsort!,
                Val(true),
            )
        end
        return
    end

    return ys
end

function _stablemsd_seq!(
    ys,
    xs,
    by::BY,
    ibyte,
    smallsize,
    smallsort!::SORT,
    xs_mutable::Union{Val{true},Val{false}} = Val(false),
    bufs = [zeros(Int, 256)],
    ibuf = 1,
) where {BY,SORT}

    offsets, allpadded = _countmsd!(fill!(bufs[ibuf], 0), xs, by, ibyte)
    @DBG counts = copy(offsets)
    allpadded && return ys

    # exclusive scan
    (acc, offsets[1]) = (offsets[1], 0)
    for i in 2:length(offsets)
        @inbounds (acc, offsets[i]) = (acc + offsets[i], acc)
    end

    _scatter!(ys, xs, eachindex(xs), offsets, by, ibyte)
    @DBG @check offsets == cumsum(counts)

    # Skip next round if know we don't need it by the type information:
    will_be_allpadded(eltype(xs), by, @stat(ibyte + 1)) && return ys

    if ibuf + 1 > length(bufs)
        push!(bufs, zeros(Int, 256))
    end
    _foreach_remaining_subrange(offsets) do idx
        ys_chunk = view(ys, idx)
        if length(idx) <= smallsize
            smallsort!(ys_chunk, by = by)
        else
            if xs_mutable === Val(true)
                xs_chunk = copyto!(view(xs, idx), ys_chunk)
            else
                xs_chunk = copy(ys_chunk)
            end
            _stablemsd_seq!(
                ys_chunk,
                xs_chunk,
                by,
                static_if_leq(@stat(ibyte + 1), static(8)),
                smallsize,
                smallsort!,
                Val(true),
                bufs,
                ibuf + 1,
            )
        end
        return
    end

    return ys
end

@inline function will_be_allpadded(::Type{ElType}, by::F, ibyte) where {ElType,F}
    if ibyte == 1
        return will_be_allpadded(ElType, by, static(1))
    elseif ibyte == 2
        return will_be_allpadded(ElType, by, static(2))
    elseif ibyte == 3
        return will_be_allpadded(ElType, by, static(3))
    elseif ibyte == 4
        return will_be_allpadded(ElType, by, static(4))
    elseif ibyte == 5
        return will_be_allpadded(ElType, by, static(5))
    elseif ibyte == 6
        return will_be_allpadded(ElType, by, static(6))
    elseif ibyte == 7
        return will_be_allpadded(ElType, by, static(7))
    elseif ibyte == 8
        return will_be_allpadded(ElType, by, static(8))
    elseif ibyte == 9
        return will_be_allpadded(ElType, by, static(9))
    end
    return false
end
@inline function will_be_allpadded(
    ::Type{ElType},
    by::F,
    ibyte::StaticInteger,
) where {ElType,F}
    T = Core.Compiler.return_type(Tuple{ElType}) do x
        msd_at(by(x), ibyte)
    end
    return T === Nothing
end

"""
    _countmsd!(counts, xs, by, ibyte) -> (counts, allpadded::Bool)
"""
function _countmsd!(counts, xs, by, ibyte)
    allpadded = true
    i = firstindex(xs)
    n = lastindex(xs) - 8
    @inbounds while i <= n
        tryprefetch(xs[i+4], ibyte)
        tryprefetch(xs[i+5], ibyte)
        tryprefetch(xs[i+6], ibyte)
        tryprefetch(xs[i+7], ibyte)
        x1 = xs[i]
        x2 = xs[i+1]
        x3 = xs[i+2]
        x4 = xs[i+3]
        k1 = msd_at(by(x1), ibyte)
        k2 = msd_at(by(x2), ibyte)
        k3 = msd_at(by(x3), ibyte)
        k4 = msd_at(by(x4), ibyte)
        if k1 !== nothing
            allpadded = false
        end
        if k2 !== nothing
            allpadded = false
        end
        if k3 !== nothing
            allpadded = false
        end
        if k4 !== nothing
            allpadded = false
        end
        counts[something(k1, 1)] += 1
        counts[something(k2, 1)] += 1
        counts[something(k3, 1)] += 1
        counts[something(k4, 1)] += 1
        i += 4
    end
    while i <= lastindex(xs)
        k = msd_at(by(@inbounds xs[i]), ibyte)
        if k !== nothing
            allpadded = false
        end
        @inbounds counts[something(k, 1)] += 1
        i += 1
    end
    return counts, allpadded
end

function _scatter!(ys, xs, idx, exclusive_offsets, by, ibyte)
    offsets = exclusive_offsets  # exclusive scan of `counts`
    i = first(idx)
    n = last(idx) - 8
    @inbounds while i <= n
        tryprefetch(xs[i+4], ibyte)
        tryprefetch(xs[i+5], ibyte)
        tryprefetch(xs[i+6], ibyte)
        tryprefetch(xs[i+7], ibyte)
        x1 = xs[i]
        x2 = xs[i+1]
        x3 = xs[i+2]
        x4 = xs[i+3]
        k1 = something(msd_at(by(x1), ibyte), 1)
        k2 = something(msd_at(by(x2), ibyte), 1)
        k3 = something(msd_at(by(x3), ibyte), 1)
        k4 = something(msd_at(by(x4), ibyte), 1)
        j1 = offsets[k1] += 1
        j2 = offsets[k2] += 1
        j3 = offsets[k3] += 1
        j4 = offsets[k4] += 1
        ys[j1] = x1
        ys[j2] = x2
        ys[j3] = x3
        ys[j4] = x4
        i += 4
    end
    while i <= last(idx)
        j = offsets[something(msd_at(by(xs[i]), ibyte), 1)] += 1
        ys[j] = xs[i]
        i += 1
    end
    # `offsets` here is now inclusive scan of `counts`
end

function _foreach_remaining_subrange(f, inclusive_offsets)
    prev = 0
    for x in eachindex(inclusive_offsets)
        i = @inbounds inclusive_offsets[x]
        if !(prev == i || prev + 1 == i)
            f(prev+1:i)
        end
        prev = i
    end
end

function _spawn_foreach_remaining_subrange(
    f,
    inclusive_offsets,
    indices::UnitRange = firstindex(inclusive_offsets):lastindex(inclusive_offsets),
)
    @inline function maybe_nonempty_subrange(subindices::UnitRange)
        local prev = subindices[1] == 1 ? 0 : inclusive_offsets[subindices[1]-1]
        local curr = inclusive_offsets[subindices[end]]
        if prev == curr || prev + 1 == curr
            return nothing
        else
            return prev+1:curr
        end
    end
    if length(indices) == 0
    elseif length(indices) == 1
        subrange = maybe_nonempty_subrange(indices)
        if subrange !== nothing
            f(subrange)
        end
    else
        m = (last(indices) - first(indices) + 1) ÷ 2 + first(indices)
        left = first(indices):m-1
        right = m:last(indices)
        leftrange = maybe_nonempty_subrange(left)
        rightrange = maybe_nonempty_subrange(right)
        if leftrange === nothing
            if rightrange === nothing
            else
                _spawn_foreach_remaining_subrange(f, inclusive_offsets, right)
            end
        else
            if rightrange === nothing
                _spawn_foreach_remaining_subrange(f, inclusive_offsets, left)
            else
                @sync begin
                    @spawn _spawn_foreach_remaining_subrange(f, inclusive_offsets, right)
                    _spawn_foreach_remaining_subrange(f, inclusive_offsets, left)
                end
            end
        end
    end
    return
end

@inline msd_at(x, ibyte::Integer) = sizeof(x) < ibyte ? nothing : _msd_at(x, ibyte) + 1

@inline _msd_at(x::Unsigned, ibyte::Integer) =
    if ibyte == sizeof(x)
        Int(x % UInt8)
    else
        Int((x >> (8 * (sizeof(x) - ibyte))) % UInt8)
    end

# Taken from SortingAlgorithms.jl:
@inline _msd_at(x::Signed, ibyte::Integer) = _msd_at(unsigned(xor(x, typemin(x))), ibyte)
@inline function _msd_at(x::Float32, ibyte::Integer)
    y = reinterpret(Int32, x)
    return _msd_at(reinterpret(UInt32, ifelse(y < 0, ~y, xor(y, typemin(Int32)))), ibyte)
end
@inline function _msd_at(x::Float64, ibyte::Integer)
    y = reinterpret(Int64, x)
    return _msd_at(reinterpret(UInt64, ifelse(y < 0, ~y, xor(y, typemin(Int64)))), ibyte)
end
# TODO: normalize NaN?

@inline _msd_at(x::AbstractChar, ibyte::Integer) = _msd_at(UInt32(x), ibyte)

@inline function msd_at(x::AbstractString, ibyte::Integer)
    n, r = divrem(ibyte - 1, sizeof(eltype(x)))
    if ibyte isa StaticInteger
        c = _get(x, static(n + 1))
    else
        c = _get(x, n + 1)
    end
    c === nothing && return nothing
    return msd_at(c, r + 1)
end

@inline function _get(x::AbstractString, i::Integer)
    @DBG @check i > 0
    y = iterate(x)
    while true
        y === nothing && return nothing
        i -= 1
        i == 0 && return y[1]
        y = iterate(x, y[2])
    end
end

@inline function _get(x::AbstractString, i::StaticInteger)
    @DBG @check i > 1
    y0 = @inbounds iterate(x)
    y0 === nothing && return nothing
    i === static(1) && return y0[1]
    return _afoldl(y0[2], ntuple(identity, static(i - 1))...) do state, j
        j::StaticInteger
        y = @inbounds iterate(x, state)
        y === nothing && return Reduced(nothing)
        j === static(i - 1) && return Reduced(y[1])
        y[2]
    end
end

@inline function msd_at(x::AbstractArray, ibyte::Integer)
    n, r = divrem(ibyte - 1, sizeof(eltype(x)))
    checkbounds(Bool, x, n + 1) || return nothing
    return msd_at((@inbounds x[n+1]), r + 1)
end

msdsort!(
    xs;
    basesize = default_basesize(xs),
    by = identity,
    smallsize = DEFAULT_SMALLSIZE,
    smallsort! = sort!,
    _copy = copy(xs),
) = _stablemsd!(xs, _copy, by, static(1), basesize, smallsize, smallsort!, Val(true))

msdsort(
    xs;
    basesize = default_basesize(xs),
    by = identity,
    smallsize = DEFAULT_SMALLSIZE,
    smallsort! = sort!,
) = _stablemsd!(copy(xs), xs, by, static(1), basesize, smallsize, smallsort!)

msdsort_seq!(
    xs;
    by = identity,
    smallsize = DEFAULT_SMALLSIZE,
    smallsort! = sort!,
    _copy = copy(xs),
) = _stablemsd_seq!(xs, _copy, by, static(1), smallsize, smallsort!, Val(true))

msdsort_seq(xs; by = identity, smallsize = DEFAULT_SMALLSIZE, smallsort! = sort!) =
    _stablemsd_seq!(copy(xs), xs, by, static(1), smallsize, smallsort!)
