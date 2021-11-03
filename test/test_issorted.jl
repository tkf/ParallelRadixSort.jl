module TestIssorted

using ParallelRadixSort
using Test

testdata_integers = [
    ("$T[...]", () -> rand(T, 2^20)) for
    T in [UInt8, UInt16, UInt32, UInt64, Int8, Int16, Int32, Int64]
]

testdata_floats = [("$T[...]", () -> rand(T, 2^20)) for T in [Float32, Float64]]

raw_testdata = """
rand('a':'z', 100)
[rand('a':'z', 10) for _ in 1:100]
[join(rand('a':'z', 10)) for _ in 1:100]
[rand(1:10, 10) for _ in 1:100]
"""

# A vector of `(datalabel, function)`
testdata_misc = map(split(raw_testdata, "\n", keepempty = false)) do code
    fn = "() -> $code"
    @debug "Evaluating: $fn"
    (code, Base.include_string(@__MODULE__, fn))
end

testdata = vcat(testdata_integers, testdata_floats, testdata_misc)

@testset "$alglabel" for (alglabel, alg) in [
    (:MSDRadixSort, ParallelRadixSort.MSDRadixSort),
    (:StableMSDRadixSort, ParallelRadixSort.StableMSDRadixSort),
    (:ParallelMSDRadixSort, ParallelRadixSort.ParallelMSDRadixSort),
    (:ParallelStableMSDRadixSort, ParallelRadixSort.ParallelStableMSDRadixSort),
    (
        "ParallelMSDRadixSort(basesize = 2^10)",
        ParallelRadixSort.ParallelMSDRadixSort(basesize = 2^10),
    ),
]
    @testset "$datalabel" for (datalabel, fn) in testdata
        vector = fn()
        @test issorted(sort(vector; alg = alg))
        eltype(vector) isa Number || continue
        @test issorted(sort(vector; alg = alg, rev = true); rev = true)
        @test issorted(sort!(copy(vector); alg = alg, rev = true); rev = true)
    end
end

end  # module
