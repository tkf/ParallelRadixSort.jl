"""
    @DBG expression_to_run_only_when_debugging
"""
macro DBG(ex)
    quote
        _debugging() && $(esc(ex))
        nothing
    end
end

_debugging() = false
enable_debug() = (@eval _debugging() = true; nothing)
disable_debug() = (@eval _debugging() = false; nothing)

struct Reduced{T}
    value::T
end

# ad-hoc terminatable foldl:
@inline _afoldl(op, x::Reduced, _...) = x.value
@inline _afoldl(op, acc, x, xs...) = _afoldl(op, op(acc, x), xs...)
@inline _afoldl(op, acc) = acc

"""
    unsafe_prefetch(
        address::Union{Ptr,Integer},
        ::Val{rw} = Val(:read),
        ::Val{locality} = Val(0),
        ::Val{cache_type} = Val(:data),
    )

# Arguments
- `rw`: `:read` (0) or `:write` (1)
- `locality`: no locality/NTA (0) -- extremely local/T0 (3)
- `cache_type`: `:data` (1) or `:instruction` (0)
"""
unsafe_prefetch

@generated function unsafe_prefetch(
    address::Union{Ptr,Integer},
    ::Val{rw} = Val(:read),
    ::Val{locality} = Val(0),
    ::Val{cache_type} = Val(:data),
) where {locality,rw,cache_type}

    rw = get(Dict(:read => 0, :write => 1), rw, rw)
    cache_type = get(Dict(:data => 1, :instruction => 0), cache_type, cache_type)

    @assert rw in (0, 1)
    @assert locality in 0:3
    @assert cache_type in (0, 1)

    declaration = "declare void @llvm.prefetch(i8*, i32, i32, i32)"

    typ = (Int === Int64 ? "i64" : "i32")
    instructions = """
    %addr = inttoptr $typ %0 to i8*
    call void @llvm.prefetch(i8* %addr, i32 $rw, i32 $locality, i32 $cache_type)
    ret void
    """
    if VERSION < v"1.6.0-DEV.674"
        IR = (declaration, instructions)
    else
        IR = (
            """
            $declaration

            define void @entry($typ) #0 {
            top:
                $instructions
            }

            attributes #0 = { alwaysinline }
            """,
            "entry",
        )
    end

    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($IR, Cvoid, Tuple{Ptr{Cvoid}}, address)
    end
end

# Since `tryprefetch` does not change the result of execution at all,
# using `return_type` like this should be OK?
@inline function tryprefetch(x::T, ibyte::Integer) where {T}
    R = Core.Compiler.return_type(pointer, Tuple{T})
    if R !== Union{} && R <: Ptr
        unsafe_prefetch(pointer(x) + ibyte)
    end
end


static_if_leq(x::Integer, ::Any) = x
static_if_leq(x::StaticInteger, y::StaticInteger) = x <= y ? x : Int(x)
