module TestIssorted

using ParallelRadixSort
using Test

testdata_integers = [
    ("$T[...]", rand(T(0):T(9), 100))
    for T in [UInt8, UInt16, UInt32, UInt64, Int8, Int16, Int32, Int64]
]

testdata_floats =
    [("$T[...]", rand(range(T(0), T(1), length = 11), 100)) for T in [Float32, Float64]]

raw_testdata = """
rand('a':'z', 100)
[rand('a':'z', 10) for _ in 1:100]
[join(rand('a':'z', 10)) for _ in 1:100]
[rand(1:10, 10) for _ in 1:100]
"""

# A vector of `(label, vector)`
testdata_misc = map(split(raw_testdata, "\n", keepempty = false)) do code
    @debug "Evaluating: $code"
    (code, Base.include_string(@__MODULE__, code))
end

testdata = vcat(testdata_integers, testdata_floats, testdata_misc)

@testset "$algname" for algname in [
    :MSDRadixSort,
    :StableMSDRadixSort,
    :ParallelMSDRadixSort,
    :ParallelStableMSDRadixSort,
]
    @testset "$label" for (label, vector) in testdata
        alg = getfield(ParallelRadixSort, algname)
        @test issorted(sort(vector; alg = alg))
        eltype(vector) isa Number || continue
        @test issorted(sort(vector; alg = alg, rev = true); rev = true)
        @test issorted(sort!(copy(vector); alg = alg, rev = true); rev = true)
    end
end

end  # module
