using JAC: HalfInt
using JAC.Basics: twice

@testset "HalfInt" begin
    @testset "Construction" begin
        @test HalfInt(3, 1) === HalfInt(6, 2)
        @test HalfInt(3) === HalfInt(3, 1)
        @test HalfInt(3.0, 1//1) === HalfInt(3)
        @test HalfInt(5, 2.0) === HalfInt(5, 2)
        @test HalfInt(0.5) === HalfInt(1, 2)
        @test HalfInt(1.0) === HalfInt(1, 1)
        @test HalfInt(-1//2) === HalfInt(-1, 2)
        @test HalfInt(AngularJ64(-3, 2)) === HalfInt(-3, 2)
        @test HalfInt(AngularJ64(4)) === HalfInt(4)
        @test HalfInt(AngularM64(-3, 2)) === HalfInt(-3, 2)
        @test HalfInt(AngularM64(4)) === HalfInt(4)
        @test_throws DomainError HalfInt(1, 3)
        @test_throws DomainError HalfInt(5.0, 2.5)
        @test_throws DomainError HalfInt(1//3)
        @test_throws InexactError HalfInt(0.8)
    end

    @testset "Conversion" begin
        # To real types
        for T in (:AbstractFloat, :Float16, :Float32, :Float64, :Rational, :(Rational{Int}), :(Rational{Int8}))
            @eval @test @inferred($T(HalfInt(5//2))) === $T(5//2)
        end
        @test @inferred(Number(HalfInt(-2))) === HalfInt(-2)
        @test @inferred(Real(HalfInt(3/2))) === HalfInt(3/2)
        @test BigFloat(HalfInt(4)) isa BigFloat
        @test BigFloat(HalfInt(4)) == BigFloat(4)
        @test BigFloat(HalfInt(-5//2)) == BigFloat(-2.5)
        # To integer types
        for T in (:Integer, :Signed, :Unsigned,
                  :Int8, :Int16, :Int32, :Int64, :Int128,
                  :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
            @eval @test @inferred($T(HalfInt(3))) === $T(3)
            @eval @test_throws InexactError $T(HalfInt(3/2))
        end
        @test BigInt(HalfInt(3)) isa BigInt
        @test BigInt(HalfInt(3)) == BigInt(3)
        @test_throws InexactError BigInt(HalfInt(3/2))
        @test Bool(HalfInt(0)) === false
        @test Bool(HalfInt(1)) === true
        @test_throws InexactError Bool(HalfInt(3/2))
        # To angular momentum types
        @test AngularJ64(HalfInt(1)) === AngularJ64(1, 1)
        @test AngularJ64(HalfInt(-5/2)) === AngularJ64(-5, 2)
        @test AngularM64(HalfInt(1)) === AngularM64(1, 1)
        @test AngularM64(HalfInt(-5/2)) === AngularM64(-5, 2)
        # Other functions
        @test @inferred(complex(HalfInt(-13//2))) === HalfInt(-13/2) + HalfInt(0)*im
        @test @inferred(real(HalfInt(-13//2))) === HalfInt(-13/2)
        @test @inferred(imag(HalfInt(-13//2))) === HalfInt(0)
        @test @inferred(float(HalfInt(-3//2))) === -1.5
    end

    @testset "Comparison" begin
        for T in (:Int8, :Int16, :Int32, :Int64, :Int128, :BigInt,
                  :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
            @eval @test HalfInt(3)   == $T(3)
            @eval @test HalfInt(5/2) != $T(3)
            @eval @test HalfInt(5/2) ≤ $T(3)
            @eval @test HalfInt(3)   ≤ $T(3)
            @eval @test HalfInt(5/2) < $T(3)
            @eval @test HalfInt(5/2) ≥ $T(2)
            @eval @test HalfInt(2)   ≥ $T(2)
            @eval @test HalfInt(5/2) > $T(2)
        end
        for T in (:HalfInt, :Float16, :Float32, :Float64, :BigFloat, :Rational)
            @eval @test HalfInt(5/2) == $T(5//2)
            @eval @test HalfInt(5/2) != $T(3)
            @eval @test HalfInt(5/2) ≤ $T(5//2)
            @eval @test HalfInt(5/2) ≤ $T(3)
            @eval @test HalfInt(5/2) < $T(3)
            @eval @test HalfInt(5/2) ≥ $T(5//2)
            @eval @test HalfInt(5/2) ≥ $T(2)
            @eval @test HalfInt(5/2) > $T(2)
        end
    end

    @testset "Hashing" begin
        @test hash(HalfInt(1)) === hash(1)
        @test hash(HalfInt(5//2)) === hash(5//2)
    end

    @testset "Arithmetic" begin
        # Unary plus, minus
        @test @inferred(+HalfInt(3/2)) === HalfInt(3/2)
        @test @inferred(-HalfInt(3/2)) === HalfInt(-3/2)
        # Addition, subtraction
        @test HalfInt(3/2) + HalfInt(5/2) === HalfInt(4)
        @test HalfInt(3/2) - HalfInt(5/2) === HalfInt(-1)
        @test HalfInt(5/2) + 2 === HalfInt(9/2)
        @test HalfInt(5/2) - 2 === HalfInt(1/2)
        @test 2 + HalfInt(5/2) === HalfInt(9/2)
        @test 2 - HalfInt(5/2) === HalfInt(-1/2)
        @test HalfInt(2) + 1.5 === 3.5
        @test HalfInt(2) - 1.5 === 0.5
        @test 1.5 + HalfInt(2) === 3.5
        @test 1.5 - HalfInt(2) === -0.5
        @test HalfInt(-1/2) + 4//3 === 5//6
        @test HalfInt(-1/2) - 4//3 === -11//6
        @test 4//3 + HalfInt(-1/2) === 5//6
        @test 4//3 - HalfInt(-1/2) === 11//6
        # Multiplication
        @test HalfInt(3/2)*HalfInt(7/2) === 5.25
        @test HalfInt(3/2)*5 === HalfInt(15/2)
        @test HalfInt(3/2)*5.0 === 7.5
        @test HalfInt(3/2)*(5//2) === 15//4
        @test 3*HalfInt(-1/2) === HalfInt(-3/2)
        @test 3.0*HalfInt(-1/2) === -1.5
        @test (3//2)*HalfInt(-1/2) === -3//4
        # Division
        @test HalfInt(5/2)/HalfInt(3/2) === 5/3
        @test HalfInt(5/2)/2.5 === 1.0
        @test HalfInt(5/2)/2 === 5/4
        @test HalfInt(5/2)/(2//1) === 5//4
        @test 2/HalfInt(3/2) === 4/3
        @test 2.0/HalfInt(3/2) === 4/3
        @test (2//1)/HalfInt(3/2) === 4//3
        # ^
        @test HalfInt(5/2)^2 == 6.25
        @test HalfInt(2)^3 == 8.0
        @test HalfInt(5/2)^2.0 ≈ 6.25
        @test HalfInt(2)^3.0 ≈ 8.0
        @test (-1)^HalfInt(2) ≈ 1
        @test 2^HalfInt(5/2) ≈ √32
        @test (-1.0)^HalfInt(2) ≈ 1
        @test (2.0)^HalfInt(5/2) ≈ √32
        @test HalfInt(5/2)^HalfInt(2) ≈ 6.25
        @test HalfInt(2)^HalfInt(3) ≈ 8.0
        # //
        @test HalfInt(3/2)//2 === 3//4
        @test 1//HalfInt(-1/2) === -2//1
        @test HalfInt(3/2)//HalfInt(5/2) === 3//5
        # floor
        @test floor(HalfInt(3/2)) === HalfInt(1)
        @test floor(HalfInt(0)) === HalfInt(0)
        @test floor(HalfInt(-1)) === HalfInt(-1)
        @test floor(HalfInt(-3/2)) === HalfInt(-2)
        # isinteger
        @test isinteger(HalfInt(2))
        @test !isinteger(HalfInt(3/2))
        # lcm
        @test lcm(1, HalfInt(1/2)) === HalfInt(1)
        @test lcm(1, HalfInt(3/2)) === HalfInt(3)
        @test lcm(2, HalfInt(3/2)) === HalfInt(6)
        @test lcm(HalfInt(1), HalfInt(1/2)) === HalfInt(1)
        @test lcm(HalfInt(1), HalfInt(3/2)) === HalfInt(3)
        @test lcm(HalfInt(2), HalfInt(3/2)) === HalfInt(6)
        @test lcm(HalfInt(3/2), HalfInt(3/2)) === HalfInt(3/2)
        @test lcm(HalfInt(5/2), HalfInt(3/2)) === HalfInt(15/2)
        # gcd
        @test gcd(1, HalfInt(3/2)) === HalfInt(1/2)
        @test gcd(2, HalfInt(3/2)) === HalfInt(1/2)
        @test gcd(3, HalfInt(3/2)) === HalfInt(3/2)
        @test gcd(3, HalfInt(9/2)) === HalfInt(3/2)
        @test gcd(HalfInt(1), HalfInt(3/2)) === HalfInt(1/2)
        @test gcd(HalfInt(2), HalfInt(3/2)) === HalfInt(1/2)
        @test gcd(HalfInt(3), HalfInt(3/2)) === HalfInt(3/2)
        @test gcd(HalfInt(3), HalfInt(9/2)) === HalfInt(3/2)
        @test gcd(HalfInt(3/2), HalfInt(3/2)) === HalfInt(3/2)
        @test gcd(HalfInt(-5/2), HalfInt(3/2)) === HalfInt(1/2)
        # gcdx
        @test ((d,u,v) = gcdx(1, HalfInt(3/2)); u*1+HalfInt(3/2)*v == d === HalfInt(1/2))
        @test ((d,u,v) = gcdx(2, HalfInt(3/2)); u*2+HalfInt(3/2)*v == d === HalfInt(1/2))
        @test ((d,u,v) = gcdx(3, HalfInt(3/2)); u*3+HalfInt(3/2)*v == d === HalfInt(3/2))
        @test ((d,u,v) = gcdx(3, HalfInt(9/2)); u*3+HalfInt(9/2)*v == d === HalfInt(3/2))
        @test ((d,u,v) = gcdx(HalfInt(1), HalfInt(3/2)); u*1+HalfInt(3/2)*v == d === HalfInt(1/2))
        @test ((d,u,v) = gcdx(HalfInt(2), HalfInt(3/2)); u*2+HalfInt(3/2)*v == d === HalfInt(1/2))
        @test ((d,u,v) = gcdx(HalfInt(3), HalfInt(3/2)); u*3+HalfInt(3/2)*v == d === HalfInt(3/2))
        @test ((d,u,v) = gcdx(HalfInt(3), HalfInt(9/2)); u*3+HalfInt(9/2)*v == d === HalfInt(3/2))
        @test ((d,u,v) = gcdx(HalfInt(3/2), HalfInt(3/2)); u*HalfInt(3/2)+HalfInt(3/2)*v == d === HalfInt(3/2))
        @test ((d,u,v) = gcdx(HalfInt(-5/2), HalfInt(3/2)); u*HalfInt(-5/2)+HalfInt(3/2)*v == d === HalfInt(1/2))
        # rem
        @test rem(HalfInt(5/2), 2) === HalfInt(1/2)
        @test rem(HalfInt(5/2), HalfInt(1/2)) === HalfInt(0)
        @test rem(HalfInt(5/2), HalfInt(7/2)) === HalfInt(5/2)
        @test rem(4, HalfInt(3/2)) === HalfInt(1)
        # sign
        @test sign(HalfInt(5/2)) == 1
        @test sign(HalfInt(1)) == 1
        @test sign(HalfInt(0)) == 0
        @test sign(HalfInt(-1/2)) == -1
    end

    @testset "Ranges" begin
        # UnitRange
        @test HalfInt(1/2):HalfInt(5) isa UnitRange{HalfInt}
        @test HalfInt(1/2):5 isa UnitRange{HalfInt}
        @test HalfInt(1/2):HalfInt(5) === HalfInt(1/2):HalfInt(9/2)
        @test length(HalfInt(1/2):HalfInt(0)) == 0
        @test length(HalfInt(1/2):HalfInt(1/2)) == 1
        @test length(HalfInt(1/2):HalfInt(1)) == 1
        @test length(HalfInt(-1/2):HalfInt(1)) == 2
        @test 2.5 ∉ HalfInt(1):HalfInt(5)
        @test 2.5 ∈ HalfInt(1/2):HalfInt(5)
        @test 5//2 ∉ HalfInt(1):HalfInt(5)
        @test 5//2 ∈ HalfInt(1/2):HalfInt(5)
        @test 2 ∈ HalfInt(1):HalfInt(5)
        @test 2 ∉ HalfInt(1/2):HalfInt(5)
        # StepRange
        @test HalfInt(2):HalfInt(1/2):HalfInt(5) isa StepRange{HalfInt}
        @test HalfInt(2):HalfInt(2):HalfInt(9/2) isa StepRange{HalfInt}
        @test HalfInt(2):2:HalfInt(9/2) isa StepRange{HalfInt}
        @test 2:2:HalfInt(9/2) isa StepRange{HalfInt}
        @test HalfInt(2):HalfInt(2):HalfInt(9/2) === HalfInt(2):HalfInt(2):HalfInt(4)
        @test HalfInt(2):HalfInt(-1):HalfInt(1/2) === HalfInt(2):HalfInt(-1):HalfInt(1)
        @test length(HalfInt(2):HalfInt(1/2):HalfInt(5)) == 7
        @test length(HalfInt(2):HalfInt(2):HalfInt(9/2)) == 2
        @test 2.5 ∈ HalfInt(1/2):HalfInt(2):HalfInt(7)
        @test 2.5 ∉ HalfInt(1/2):HalfInt(5/2):HalfInt(7)
        @test 2.5 ∉ HalfInt(1/2):HalfInt(-1):HalfInt(-1/2)
        @test 5//2 ∈ HalfInt(1/2):HalfInt(2):HalfInt(7)
        @test 5//2 ∈ HalfInt(1/2):2:HalfInt(7)
        @test 5//2 ∈ HalfInt(1/2):2:7
        @test 5//2 ∉ HalfInt(1/2):HalfInt(5/2):HalfInt(7)
        @test 5//2 ∉ HalfInt(1/2):HalfInt(-1):HalfInt(-1/2)
        @test 5//2 ∉ HalfInt(1/2):-1:HalfInt(-1/2)
        @test 2 ∈ HalfInt(0):HalfInt(1):HalfInt(5/2)
        @test 2 ∉ HalfInt(1/2):HalfInt(1):HalfInt(5/2)
        # intersect
        @test @inferred((HalfInt(1):HalfInt(3)) ∩ (HalfInt(1/2):HalfInt(3))) isa UnitRange{HalfInt}
        @test @inferred((1:3) ∩ (HalfInt(1/2):HalfInt(3))) isa UnitRange{HalfInt}
        @test @inferred((HalfInt(1/2):HalfInt(3)) ∩ (1:3)) isa UnitRange{HalfInt}
        @test @inferred((HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (HalfInt(1):HalfInt(3))) isa StepRange{HalfInt}
        @test @inferred((HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (1:3)) isa StepRange{HalfInt}
        @test @inferred((1:3) ∩ (HalfInt(3/2):HalfInt(1):HalfInt(5))) isa StepRange{HalfInt}
        @test @inferred((1:3:7) ∩ (HalfInt(3/2):HalfInt(5))) isa StepRange{HalfInt}
        @test @inferred((HalfInt(3/2):HalfInt(5)) ∩ (1:3:7)) isa StepRange{HalfInt}
        @test @inferred((HalfInt(3/2):HalfInt(1/2):HalfInt(5)) ∩ (HalfInt(-1):HalfInt(3):HalfInt(5))) isa StepRange{HalfInt}
        @test @inferred((-1:3:5) ∩ (HalfInt(3/2):HalfInt(1/2):HalfInt(5))) isa StepRange{HalfInt}
        @test @inferred((HalfInt(3/2):HalfInt(1/2):HalfInt(5)) ∩ (-1:3:5)) isa StepRange{HalfInt}
        @test isempty((HalfInt(1):HalfInt(3)) ∩ (HalfInt(1/2):HalfInt(3)))
        @test (HalfInt(1/2):HalfInt(2)) ∩ (HalfInt(3/2):HalfInt(3)) == HalfInt(3/2):HalfInt(3/2)
        @test isempty((HalfInt(1/2):HalfInt(2)) ∩ (1:2))
        @test isempty((1:2) ∩ (HalfInt(1/2):HalfInt(2)))
        @test (HalfInt(2):HalfInt(5)) ∩ (1:3) == HalfInt(2):HalfInt(3)
        @test (1:3) ∩ (HalfInt(2):HalfInt(5)) == HalfInt(2):HalfInt(3)
        @test (HalfInt(2):HalfInt(1):HalfInt(5)) ∩ (1:3) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (1:3) ∩ (HalfInt(2):HalfInt(1):HalfInt(5)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test isempty((HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (1:3))
        @test isempty((1:3) ∩ (HalfInt(3/2):HalfInt(1):HalfInt(5)))
        @test (HalfInt(2):HalfInt(1/2):HalfInt(5)) ∩ (1:3) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (1:3) ∩ (HalfInt(2):HalfInt(1/2):HalfInt(5)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(3/2):HalfInt(1/2):HalfInt(5)) ∩ (1:3) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (1:3) ∩ (HalfInt(3/2):HalfInt(1/2):HalfInt(5)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(2):HalfInt(1):HalfInt(5)) ∩ (HalfInt(1):HalfInt(3)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(1):HalfInt(3)) ∩ (HalfInt(2):HalfInt(1):HalfInt(5)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test isempty((HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (HalfInt(1):HalfInt(3)))
        @test isempty((HalfInt(1):HalfInt(3)) ∩ (HalfInt(3/2):HalfInt(1):HalfInt(5)))
        @test (HalfInt(2):HalfInt(1/2):HalfInt(5)) ∩ (HalfInt(1):HalfInt(3)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(1):HalfInt(3)) ∩ (HalfInt(2):HalfInt(1/2):HalfInt(5)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(3/2):HalfInt(1/2):HalfInt(5)) ∩ (HalfInt(1):HalfInt(3)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(1):HalfInt(3)) ∩ (HalfInt(3/2):HalfInt(1/2):HalfInt(5)) == HalfInt(2):HalfInt(1):HalfInt(3)
        @test (HalfInt(5/2):HalfInt(3/2):HalfInt(7)) ∩ (HalfInt(1):HalfInt(8)) == HalfInt(4):HalfInt(3):HalfInt(7)
        @test (HalfInt(1):HalfInt(8)) ∩ (HalfInt(5/2):HalfInt(3/2):HalfInt(7)) == HalfInt(4):HalfInt(3):HalfInt(7)
        @test (1:3:7) ∩ (HalfInt(3):HalfInt(5)) == HalfInt(4):HalfInt(3):HalfInt(4)
        @test (HalfInt(3):HalfInt(5)) ∩ (1:3:7) == HalfInt(4):HalfInt(3):HalfInt(4)
        @test isempty((1:3:7) ∩ (HalfInt(3/2):HalfInt(5)))
        @test isempty((HalfInt(3/2):HalfInt(5)) ∩ (1:3:7))
        @test (HalfInt(3/2):HalfInt(1/2):HalfInt(5)) ∩ (-1:3:5) == HalfInt(2):HalfInt(3):HalfInt(5)
        @test (-1:3:5) ∩ (HalfInt(3/2):HalfInt(1/2):HalfInt(5)) == HalfInt(2):HalfInt(3):HalfInt(5)
        @test (HalfInt(2):HalfInt(1):HalfInt(5)) ∩ (-1:3:5) == HalfInt(2):HalfInt(3):HalfInt(5)
        @test (-1:3:5) ∩ (HalfInt(2):HalfInt(1):HalfInt(5)) == HalfInt(2):HalfInt(3):HalfInt(5)
        @test isempty((HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (-1:3:5))
        @test isempty((-1:3:5) ∩ (HalfInt(3/2):HalfInt(1):HalfInt(5)))
        @test isempty((HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (HalfInt(-1):HalfInt(3):HalfInt(5)))
        @test isempty((HalfInt(-1):HalfInt(3):HalfInt(5)) ∩ (HalfInt(3/2):HalfInt(1):HalfInt(5)))
        @test (HalfInt(3/2):HalfInt(1):HalfInt(5)) ∩ (HalfInt(-1):HalfInt(3/2):HalfInt(5)) == HalfInt(7/2):HalfInt(3):HalfInt(7/2)
    end

    @testset "string" begin
        @test string(HalfInt(5)) == "5"
        @test string(HalfInt(-1)) == "-1"
        @test string(HalfInt(1/2)) == "1/2"
        @test string(HalfInt(-7/2)) == "-7/2"
    end

    @testset "twice" begin
        @test @inferred(twice(HalfInt(3/2))) === 3
        @test @inferred(twice(-2)) === -4
        @test @inferred(twice(1.6)) === 3.2
        @test @inferred(twice(2//3)) === 4//3
    end
end
