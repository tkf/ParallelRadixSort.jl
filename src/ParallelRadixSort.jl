module ParallelRadixSort

export MSDRadixSort, ParallelMSDRadixSort, ParallelStableMSDRadixSort, StableMSDRadixSort

# macro inbounds(ex)
#     esc(ex)
# end

using ArgCheck: @check
using Base.Order: Order, Ordering, ForwardOrdering, ReverseOrdering
using ConstructionBase: setproperties
using StaticArrays: MVector
using StaticNumbers: @stat, StaticInteger, maybe_static, static

if isdefined(Base.Experimental, :Tapir)
    using Base.Experimental.Tapir: @sync, @spawn
else
    using Base.Threads: @sync, @spawn
end

abstract type RadixSortAlgorithm <: Base.Sort.Algorithm end

(alg::RadixSortAlgorithm)(; kw...) = setproperties(alg; kw...)

include("utils.jl")
include("msd.jl")

end # module
