module TestParallelRadixSort
using Test

if lowercase(get(ENV, "CI", "false")) == "true" && VERSION ≥ v"1.7-"
    # Use https://github.com/perrutquist/StaticNumbers.jl/pull/12 as otherwise
    # the test would be extremely slow.
    Pkg = Base.require(Base.PkgId(Base.UUID(0x44cfe95a1eb252eab672e2afdf69b78f), "Pkg"))
    Pkg.add(
        Pkg.PackageSpec(;
            name = "StaticNumbers",
            url = "https://github.com/tkf/StaticNumbers.jl.git",
            rev = "kwargs.data",
        ),
    )
end

@testset "$file" for file in sort([
    file for file in readdir(@__DIR__) if match(r"^test_.*\.jl$", file) !== nothing
])
    include(file)
end

end  # module
