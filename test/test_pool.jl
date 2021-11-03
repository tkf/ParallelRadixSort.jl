module TestPool

using Test
using ParallelRadixSort: Pool, initpool!, alloc!, free!

zeros3() = zeros(3)

@testset begin
    capacity = 4
    pool = Pool(zeros3, capacity)
    initpool!(pool)
    xs1 = alloc!(pool, 2)
    @test length(xs1) == 2
    xs2 = alloc!(pool, 4)
    @test length(xs2) == 4
    free!(pool, xs1)
    free!(pool, xs2)
    @test length(pool.buffers[Threads.threadid()]) == capacity
end

end  # module
