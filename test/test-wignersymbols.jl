using JAC: HalfInt
using JAC.AngularMomentum: clebschgordan, wigner3j, wigner6j, wigner9j

@testset "Wigner symbols" begin
    @testset "clebschgordan" begin
        @test iszero(clebschgordan(1, 0, 1, 0, 1/2, 0))
        @test iszero(clebschgordan(1, 0, 1, 0, 1, 1))
        @test clebschgordan(5/2, 3/2, 0, 0, 5//2, 3//2) ≈ 1
        @test clebschgordan(AngularJ64(5//2), AngularM64(3//2),
                            AngularJ64(0),    AngularM64(0),
                            AngularJ64(5//2), AngularM64(3//2)) ≈ 1
        @test clebschgordan(1/2,  1/2, 1/2,  1/2, 1,  1) ≈ 1
        @test clebschgordan(1/2,  1/2, 1/2, -1/2, 1,  0) ≈ √(1/2)
        @test clebschgordan(1/2, -1/2, 1/2,  1/2, 1,  0) ≈ √(1/2)
        @test clebschgordan(1/2, -1/2, 1/2, -1/2, 1, -1) ≈ 1
        @test clebschgordan(1/2,  1/2, 1/2, -1/2, 0,  0) ≈ √(1/2)
        @test clebschgordan(1/2, -1/2, 1/2,  1/2, 0,  0) ≈ -√(1/2)
        @test clebschgordan(3/2, -1/2, 2, 1, 7/2,  1/2) ≈ √(12/35)
        @test clebschgordan(3/2, -1/2, 2, 1, 5/2,  1/2) ≈ -√(5/14)
        @test clebschgordan(3/2, -1/2, 2, 1, 3/2,  1/2) ≈ 0
        @test clebschgordan(3/2, -1/2, 2, 1, 1/2,  1/2) ≈ √(3/10)
        @test clebschgordan(2,  2, 3/2, -3/2, 7/2,  1/2) ≈ √(1/35)
        @test clebschgordan(2,  2, 3/2, -3/2, 5/2,  1/2) ≈ √(6/35)
        @test clebschgordan(2,  2, 3/2, -3/2, 3/2,  1/2) ≈ √(2/5)
        @test clebschgordan(2,  2, 3/2, -3/2, 1/2,  1/2) ≈ √(2/5)
        @test clebschgordan(2,  1, 3/2, -1/2, 7/2,  1/2) ≈ √(12/35)
        @test clebschgordan(2,  1, 3/2, -1/2, 5/2,  1/2) ≈ √(5/14)
        @test clebschgordan(2,  1, 3/2, -1/2, 3/2,  1/2) ≈ 0
        @test clebschgordan(2,  1, 3/2, -1/2, 1/2,  1/2) ≈ -√(3/10)
        @test clebschgordan(2,  0, 3/2,  1/2, 7/2,  1/2) ≈ √(18/35)
        @test clebschgordan(2,  0, 3/2,  1/2, 5/2,  1/2) ≈ -√(3/35)
        @test clebschgordan(2,  0, 3/2,  1/2, 3/2,  1/2) ≈ -√(1/5)
        @test clebschgordan(2,  0, 3/2,  1/2, 1/2,  1/2) ≈ √(1/5)
        @test clebschgordan(2, -1, 3/2,  3/2, 7/2,  1/2) ≈ √(4/35)
        @test clebschgordan(2, -1, 3/2,  3/2, 5/2,  1/2) ≈ -√(27/70)
        @test clebschgordan(2, -1, 3/2,  3/2, 3/2,  1/2) ≈ √(2/5)
        @test clebschgordan(2, -1, 3/2,  3/2, 1/2,  1/2) ≈ -√(1/10)
        @test clebschgordan(5/2,  3/2, 5/2, -3/2, 5,  0) ≈ √(25/252)
        @test clebschgordan(5/2,  3/2, 5/2, -3/2, 4,  0) ≈ √(9/28)
        @test clebschgordan(5/2,  3/2, 5/2, -3/2, 3,  0) ≈ √(49/180)
        @test clebschgordan(5/2,  3/2, 5/2, -3/2, 2,  0) ≈ √(1/84)
        @test clebschgordan(5/2,  3/2, 5/2, -3/2, 1,  0) ≈ -√(9/70)
        @test clebschgordan(5/2,  3/2, 5/2, -3/2, 0,  0) ≈ -√(1/6)
        @test clebschgordan(5/2, -3/2, 5/2,  3/2, 5,  0) ≈ √(25/252)
        @test clebschgordan(5/2, -3/2, 5/2,  3/2, 4,  0) ≈ -√(9/28)
        @test clebschgordan(5/2, -3/2, 5/2,  3/2, 3,  0) ≈ √(49/180)
        @test clebschgordan(5/2, -3/2, 5/2,  3/2, 2,  0) ≈ -√(1/84)
        @test clebschgordan(5/2, -3/2, 5/2,  3/2, 1,  0) ≈ -√(9/70)
        @test clebschgordan(5/2, -3/2, 5/2,  3/2, 0,  0) ≈ √(1/6)
        # Negative numbers as arguments for J
        @info "The following 20 warnings from GSL.jl are expected (testing error handling)"
        @test_throws DomainError clebschgordan(-5/2, 3/2, 0, 0, 5/2, 3/2)
        @test_throws DomainError clebschgordan(0, 0, -5/2, 3/2, 5/2, 3/2)
        # In the following case, GSL.jl does not print a warning (sqrt(-4.0) throws
        # DomainError before wigner3j is called)
        @test_throws DomainError clebschgordan(0, 0, 5/2, 3/2, -5/2, 3/2)
        # Arguments that are not half-integers
        @test_throws InexactError clebschgordan(5/3, 3/2, 0, 0, 5/2, 3/2)
        @test_throws InexactError clebschgordan(5/2, 4/3, 0, 0, 5/2, 3/2)
        @test_throws InexactError clebschgordan(0, 0, 5/3, 3/2, 5/2, 3/2)
        @test_throws InexactError clebschgordan(0, 0, 5/2, 4/3, 5/2, 3/2)
        @test_throws InexactError clebschgordan(5/3, 3/2, 0, 0, 5/3, 3/2)
        @test_throws InexactError clebschgordan(5/2, 4/3, 0, 0, 5/2, 4/3)
    end

    @testset "wigner3j" begin
        @test wigner3j(AngularJ64(2), AngularJ64(6), AngularJ64(4),
                       AngularM64(0), AngularM64(0), AngularM64(0)) ≈ sqrt(5/143)
        @test wigner3j(AngularJ64(2), AngularJ64(6), AngularJ64(4),
                       AngularM64(0), AngularM64(0), AngularM64(1)) == 0
        @test wigner3j(AngularJ64(2), AngularJ64(7), AngularJ64(4),
                       AngularM64(0), AngularM64(0), AngularM64(0)) == 0
        @test wigner3j(AngularJ64(3//2), AngularJ64(1), AngularJ64(5//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ -sqrt(1/10)
        @test wigner3j(AngularJ64(5//2), AngularJ64(1), AngularJ64(5//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ sqrt(1/210)
        @test wigner3j(AngularJ64(3//2), AngularJ64(0), AngularJ64(3//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ 1/2
        @test wigner3j(AngularJ64(3//2), AngularJ64(2), AngularJ64(3//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ -sqrt(1/20)
        @test wigner3j(AngularJ64(7//2), AngularJ64(2), AngularJ64(5//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ -sqrt(1/210)
        @test wigner3j(AngularJ64(7//2), AngularJ64(4), AngularJ64(5//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ sqrt(5/462)
        @test wigner3j(AngularJ64(7//2), AngularJ64(6), AngularJ64(5//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ -5*sqrt(1/858)
        @test wigner3j(AngularJ64(7//2), AngularJ64(0), AngularJ64(7//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ 1/2*sqrt(1/2)
        @test wigner3j(AngularJ64(7//2), AngularJ64(2), AngularJ64(7//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ -1/2*sqrt(5/42)
        @test wigner3j(AngularJ64(7//2), AngularJ64(4), AngularJ64(7//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ 3/2*sqrt(1/154)
        @test wigner3j(AngularJ64(7//2), AngularJ64(6), AngularJ64(7//2),
                       AngularM64(1//2), AngularM64(0), AngularM64(-1//2)) ≈ -5/2*sqrt(1/858)
        @test wigner3j(2, 6, 4, 0, 0, 0) ≈ sqrt(5/143)
        @test wigner3j(HalfInt(3/2), HalfInt(1), HalfInt(5/2), 1//2, 0, -1/2) ≈ -sqrt(1/10)
        @test wigner3j(3/2, 1, 5/2, 1/2, 0, -1/2) ≈ -sqrt(1/10)
        # Negative numbers as arguments for J
        @test_throws DomainError wigner3j(-3/2, 1, 5/2, 1//2, 0, -1/2)
        @test_throws DomainError wigner3j(3/2, -1, 5/2, 1//2, 0, -1/2)
        @test_throws DomainError wigner3j(3/2, 1, -5/2, 1//2, 0, -1/2)
        # Arguments that are not half-integers
        @test_throws InexactError wigner3j(4/3, 1, 5/2, 1//2, 0, -1/2)
        @test_throws InexactError wigner3j(3/2, 1.1, 5/2, 1//2, 0, -1/2)
        @test_throws InexactError wigner3j(3/2, 1, 5/3, 1//2, 0, -1/2)
        @test_throws InexactError wigner3j(3/2, 1, 5/2, 1//3, 0, -1/2)
        @test_throws InexactError wigner3j(3/2, 1, 5/2, 1//2, 0.1, -1/2)
        @test_throws InexactError wigner3j(3/2, 1, 5/2, 1//2, 0, -1/5)
    end

    @testset "wigner6j" begin
        @test wigner6j(AngularJ64(3), AngularJ64(3), AngularJ64(3),
                       AngularJ64(3), AngularJ64(3), AngularJ64(3)) ≈ -1/14
        @test wigner6j(AngularJ64(5), AngularJ64(5), AngularJ64(5),
                       AngularJ64(5), AngularJ64(5), AngularJ64(5)) ≈ 1/52
        @test wigner6j(AngularJ64(1), AngularJ64(3//2), AngularJ64(5//2),
                       AngularJ64(0), AngularJ64(5//2), AngularJ64(3//2)) ≈ -sqrt(1/24)
        @test wigner6j(1, 3/2, 5/2, 0, 5/2, 3/2) ≈ -sqrt(1/24)
        @test wigner6j(HalfInt(1), HalfInt(3/2), HalfInt(5/2), 0, 5//2, 3/2) ≈ -sqrt(1/24)
        # Negative numbers as arguments
        @test_throws DomainError wigner6j(-1, 3/2, 5/2, 0, 5//2, 3/2)
        @test_throws DomainError wigner6j(1, -3/2, 5/2, 0, 5//2, 3/2)
        @test_throws DomainError wigner6j(1, 3/2, -5/2, 0, 5//2, 3/2)
        @test_throws DomainError wigner6j(3/2, 5/2, 1, -5//2, 3/2, 0)
        @test_throws DomainError wigner6j(1, 3/2, 5/2, 0, -5//2, 3/2)
        @test_throws DomainError wigner6j(1, 3/2, 5/2, 0, 5//2, -3/2)
        # Arguments that are not half-integers
        @test_throws InexactError wigner6j(1.1, 3/2, 5/2, 0, 5//2, 3/2)
        @test_throws InexactError wigner6j(1, 4/3, 5/2, 0, 5//2, 3/2)
        @test_throws InexactError wigner6j(1, 3/2, 5/3, 0, 5//2, 3/2)
        @test_throws InexactError wigner6j(1, 3/2, 5/2, 0.1, 5//2, 3/2)
        @test_throws InexactError wigner6j(1, 3/2, 5/2, 0, 5//3, 3/2)
        @test_throws InexactError wigner6j(1, 3/2, 5/2, 0, 5//2, 4/3)
    end

    @testset "wigner9j" begin
        @test wigner9j(AngularJ64(1), AngularJ64(1), AngularJ64(1),
                       AngularJ64(1), AngularJ64(1), AngularJ64(1),
                       AngularJ64(1), AngularJ64(1), AngularJ64(0)) ≈ 1/18
        @test wigner9j(AngularJ64(1),    AngularJ64(3//2), AngularJ64(5//2),
                       AngularJ64(5//2), AngularJ64(0),    AngularJ64(5//2),
                       AngularJ64(3//2), AngularJ64(3//2), AngularJ64(0)) ≈ -1/24
        @test wigner9j(AngularJ64(1//2), AngularJ64(1//2), AngularJ64(1),
                       AngularJ64(3),    AngularJ64(3),    AngularJ64(2),
                       AngularJ64(5//2), AngularJ64(5//2), AngularJ64(1)) ≈ (1/21)*sqrt(2/5)
        @test wigner9j(1, 1, 1, 1, 1, 1, 1, 1, 0) ≈ 1/18
        @test wigner9j(1, AngularJ64(1), HalfInt(1), 1//1, 1.0, 1, 1, 1, 0) ≈ 1/18
        # Negative numbers as arguments
        @test_throws DomainError wigner9j(-1, 1, HalfInt(1), 1//1, 1.0, 1, 1, 1, 0)
        @test_throws DomainError wigner9j(1, -1, HalfInt(1), 1//1, 1.0, 1, 1, 1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(-1), 1//1, 1.0, 1, 1, 1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(1), -1//1, 1.0, 1, 1, 1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(1), 1//1, -1.0, 1, 1, 1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, -1, 1, 1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1, -1, 1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1, 1, -1, 0)
        @test_throws DomainError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1, 1, 1, -1)
        # Arguments that are not half-integers
        @test_throws InexactError wigner9j(1.1, 1, HalfInt(1), 1//1, 1.0, 1, 1, 1, 0)
        @test_throws InexactError wigner9j(1, 1.1, HalfInt(1), 1//1, 1.0, 1, 1, 1, 0)
        @test_throws InexactError wigner9j(1, 1, 1.1, 1//1, 1.0, 1, 1, 1, 0)
        @test_throws InexactError wigner9j(1, 1, HalfInt(1), 1//3, 1.0, 1, 1, 1, 0)
        @test_throws InexactError wigner9j(1, 1, HalfInt(1), 1//1, 1.6, 1, 1, 1, 0)
        @test_throws InexactError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1.1, 1, 1, 0)
        @test_throws InexactError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1, 1.2, 1, 0)
        @test_throws InexactError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1, 1, 1.1, 0)
        @test_throws InexactError wigner9j(1, 1, HalfInt(1), 1//1, 1.0, 1, 1, 1.3, 0.01)
    end
end
