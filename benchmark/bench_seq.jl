module BenchSeq

import Random
using BenchmarkTools
using ParallelRadixSort: msdsort_seq

Random.seed!(1234)

n = 10_000
datasets = [
    # (label, setup)
    ("UInt8 (wide)", rand(UInt8, 100_000)),
    ("UInt64 (wide)", rand(UInt64, 100_000)),
    ("Float64 (wide)", randn(100_000)),
    # `Base` has `Base.Sort.sort_int_range!` for those cases:
    # ("UInt8 (narrow)", rand(UInt8(0):UInt8(9), 100_000)),
    # ("UInt64 (narrow)", rand(UInt64(0):UInt64(9), 100_000)),
    ("Float64 (narrow)", rand(0:0.1:0.9, 100_000)),
    ("String", [join(rand('a':'z', 10)) for _ in 1:10_000]),
    ("Vector{Int}", [rand(0:9, 10) for _ in 1:10_000]),
]
# Poor performance:
# * narrow
# * Vector{Int}

suite = BenchmarkGroup()

for (label, xs) in datasets
    s = suite[label] = BenchmarkGroup()
    s["radix"] = @benchmarkable(msdsort_seq($xs))
    s["base"] = @benchmarkable(sort($xs))
end

end  # module
BenchSeq.suite
