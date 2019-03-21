"""
    HalfInteger <: Real

Abstract supertype for half-integers.
"""
abstract type HalfInteger <: Real end

(T::Type{<:Rational})(x::HalfInteger)      = (t=twice(x); T(t, oftype(t,2)))
(T::Type{<:AbstractFloat})(x::HalfInteger) = T(float(x))
(T::Type{<:Integer})(x::HalfInteger)       = T(Integer(x))

Base.Bool(x::HalfInteger) = Bool(Integer(x))

function Base.Integer(x::HalfInteger)
    if isinteger(x)
        t = twice(x)
        t ÷ oftype(t,2)
    else
        throw(InexactError(:Integer, Integer, x))
    end
end

Base.:<(x::HalfInteger, y::HalfInteger)    = twice(x) < twice(y)
Base.:≤(x::HalfInteger, y::HalfInteger)    = twice(x) ≤ twice(y)
Base.:(==)(x::HalfInteger, y::HalfInteger) = twice(x) == twice(y)

# For correct hashing
Base.decompose(x::HalfInteger) = twice(x), 0, 2

Base.float(x::HalfInteger) = twice(x)/2

Base.isinteger(x::HalfInteger) = iseven(twice(x))

"""
    HalfInt <: HalfInteger

Number type for representing values `n/2` where `n` is an `Int`.

!!! warning

    Since the implementation of `HalfInt` is based on standard `Int` arithmetic, it is
    subject to integer overflow:

    ```jldoctest
    julia> x = typemax(HalfInt)
    $(typemax(Int))/2

    julia> x + 1
    -$(typemax(Int))/2

    julia> y = typemax(Int64)
    9223372036854775807

    julia> HalfInt(y) # 2*y == -2
    -1
    ```
"""
struct HalfInt <: HalfInteger
    twice::Int # the number 2x, where x is the value represented by the HalfInt object
    @doc """
        HalfInt(num, den)

    Create a `HalfInt` with the numerator `num` and denominator `den`. The denominator
    must be equal to `1` or `2`.

    # Examples

    ```jldoctest
    julia> HalfInt(5, 2)
    5/2

    julia> HalfInt(-3, 1)
    -3

    julia> HalfInt(1, 3)
    ERROR: DomainError with 3:
    denominator must be 1 or 2.
    Stacktrace:
    [...]
    ```
    """
    function HalfInt(num, den)
        den == 1 && return new(twice(num))
        den == 2 && return new(num)
        throw(DomainError(den, "denominator must be 1 or 2."))
    end
end

"""
    HalfInt(x)

Create a `HalfInt` equal to the number `x`.

# Examples

```jldoctest
julia> HalfInt(2)
2

julia> HalfInt(1.5)
3/2

julia> HalfInt(7//2)
7/2
```
"""
HalfInt(x)           = HalfInt(x, 1)
HalfInt(x::Rational) = HalfInt(numerator(x), denominator(x))
HalfInt(x::HalfInt)  = x

Base.:+(x::HalfInt, y::HalfInt) = HalfInt(twice(x)+twice(y), 2)
Base.:-(x::HalfInt, y::HalfInt) = HalfInt(twice(x)-twice(y), 2)
Base.:-(x::HalfInt)             = HalfInt(-twice(x), 2)

# Multiplication of two `HalfInt`s yields a Float64,
Base.:*(x::HalfInt, y::HalfInt) = twice(x)*twice(y)/4
# Multiplication of `HalfInt` and `Integer` yields a `HalfInt`
Base.:*(x::HalfInt, y::Integer) = HalfInt(twice(x)*y, 2)
Base.:*(x::Integer, y::HalfInt) = HalfInt(x*twice(y), 2)

Base.:/(x::HalfInt, y::HalfInt) = twice(x)/twice(y)

Base.://(x::HalfInt, y::HalfInt) = twice(x)//twice(y)
# Promotion doesn’t work for `Base.://`
Base.://(a::HalfInteger, b) = twice(a)//twice(b)
Base.://(a, b::HalfInteger) = twice(a)//twice(b)

Base.:^(x::Real, y::HalfInt) = x^float(y)

Base.gcd(x::HalfInt, y::HalfInt) = halfint_gcd(x,y)
Base.gcd(x::HalfInt, y::Integer) = halfint_gcd(x,y)
Base.gcd(x::Integer, y::HalfInt) = halfint_gcd(x,y)

Base.gcdx(x::HalfInt, y::HalfInt) = halfint_gcdx(x,y)
Base.gcdx(x::HalfInt, y::Integer) = halfint_gcdx(x,y)
Base.gcdx(x::Integer, y::HalfInt) = halfint_gcdx(x,y)

Base.lcm(x::HalfInt, y::HalfInt) = halfint_lcm(x,y)
Base.lcm(x::HalfInt, y::Integer) = halfint_lcm(x,y)
Base.lcm(x::Integer, y::HalfInt) = halfint_lcm(x,y)

# gcd where at least one of the arguments is a HalfInt
halfint_gcd(x, y) = HalfInt(gcd(twice(x), twice(y)), 2)
# lcm where at least one of the arguments is a HalfInt
halfint_lcm(x, y) = HalfInt(lcm(twice(x), twice(y)), 2)
# gcdx where at least one of the arguments is a HalfInt
function halfint_gcdx(x, y)
    d, u, v = gcdx(twice(x), twice(y))
    HalfInt(d, 2), u, v
end

Base.promote_rule(::Type{HalfInt}, T::Type)           = promote_type(Rational{Int}, T)
Base.promote_rule(::Type{HalfInt}, ::Type{<:Integer}) = HalfInt

Base.rem(x::HalfInt, y::HalfInt) = HalfInt(rem(twice(x),twice(y)), 2)

Base.round(x::HalfInt, ::RoundingMode{:Down}) = ifelse(isinteger(x), x, HalfInt(twice(x)-1, 2))

Base.show(io::IO, x::HalfInt) = isinteger(x) ? print(io, Integer(x)) : print(io, twice(x), "/2")

Base.typemax(::Type{HalfInt}) = HalfInt(typemax(Int), 2)
Base.typemin(::Type{HalfInt}) = HalfInt(typemin(Int), 2)

Base.intersect(r::StepRange{HalfInt}, s::StepRange{HalfInt})   = halfint_intersect(r, s)
Base.intersect(r::StepRange{HalfInt}, s::StepRange{<:Integer}) = halfint_intersect(r, s)
Base.intersect(r::StepRange{<:Integer}, s::StepRange{HalfInt}) = halfint_intersect(r, s)

Base.intersect(r::OrdinalRange{HalfInt}, s::OrdinalRange{HalfInt})   = halfint_intersect(r, s)
Base.intersect(r::OrdinalRange{HalfInt}, s::OrdinalRange{<:Integer}) = halfint_intersect(r, s)
Base.intersect(r::OrdinalRange{<:Integer}, s::OrdinalRange{HalfInt}) = halfint_intersect(r, s)

# Intersection of two ranges where the eltype of at least on of the ranges is HalfInt
function halfint_intersect(r::AbstractUnitRange, s::AbstractUnitRange)
    ifelse(isinteger(first(r)) ⊻ isinteger(first(s)),
           one(HalfInt):zero(HalfInt),
           max(first(r),first(s)):min(last(r),last(s)))
end

function halfint_intersect(r::OrdinalRange, s::OrdinalRange)
    twice_r = twice(first(r)):twice(step(r)):twice(last(r))
    twice_s = twice(first(s)):twice(step(s)):twice(last(s))
    twice_range = intersect(twice_r, twice_s)
    HalfInt(first(twice_range),2):HalfInt(step(twice_range),2):HalfInt(last(twice_range),2)
end

"""
    twice(x)

Return a value equal to `2*x`. For `HalfInteger`s, this returns an integer type. (not
exported)

# Examples

```jldoctest
julia> twice(3)
6

julia> twice(1.5)
3.0

julia> twice(HalfInt(1.5))
3
```
"""
twice(x)          = x + x
twice(x::HalfInt) = x.twice

# Conversion from and to AngularJ64, AngularM64

twice(x::Union{AngularJ64,AngularM64}) = ifelse(isone(x.den), twice(x.num), x.num)

AngularJ64(x::HalfInt) = isinteger(x) ? AngularJ64(Integer(x)) : AngularJ64(twice(x), 2)
AngularM64(x::HalfInt) = isinteger(x) ? AngularM64(Integer(x)) : AngularM64(twice(x), 2)
