using HalfIntegers: HalfInt
using JAC.AngularMomentum: clebschgordan, wigner3j, wigner6j, wigner9j

@testset "@racahsum" begin
    # Clebsch–Gordan coefficients
    @test @racahsum(clebschgordan(5/2, 3/2, 0, 0, 5/2, M), M) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(5/2, 3/2, 0, 0, 5/2,-M), M) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(5/2, 3/2, 0, 0, J, M), J, M) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(5/2, 3/2, 0, 0, J,-M), J, M) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(5/2, 3/2, 0, 0, J, M), M, J) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(5/2, 3/2, 0, 0, J,-M), M, J) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(j₁, m₁, 0, 0, 5/2, 3/2), m₁, j₁) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(j₁,-m₁, 0, 0, 5/2, 3/2), m₁, j₁) ≈ clebschgordan(5/2, 3/2, 0, 0, 5/2, 3/2)
    @test @racahsum(clebschgordan(5/2, 3/2, j₂, m₂, 5/2, 3/2), j₂, m₂) ≈
        sum(clebschgordan(5/2, 3/2, j₂, 0, 5/2, 3/2) for j₂=0:5)
    @test @racahsum(clebschgordan(5/2, 3/2, j₂,-m₂, 5/2, 3/2), j₂, m₂) ≈
        sum(clebschgordan(5/2, 3/2, j₂, 0, 5/2, 3/2) for j₂=0:5)
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J, 1/2), J) ≈
        sum(clebschgordan(2, 0, 3/2, 1/2, J, 1/2) for J = HalfInt(1/2):HalfInt(7/2))
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J, M), J, M) ≈
        sum(clebschgordan(2, 0, 3/2, 1/2, J, 1/2) for J = HalfInt(1/2):HalfInt(7/2))
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J,-M), J, M) ≈
        sum(clebschgordan(2, 0, 3/2, 1/2, J, 1/2) for J = HalfInt(1/2):HalfInt(7/2))
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J, M), J=HalfInt(3/2):HalfInt(9/2), M) ≈
        sum(clebschgordan(2, 0, 3/2, 1/2, J, 1/2) for J = HalfInt(3/2):HalfInt(7/2))
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J,-M), J=HalfInt(3/2):HalfInt(9/2), M) ≈
        sum(clebschgordan(2, 0, 3/2, 1/2, J, 1/2) for J = HalfInt(3/2):HalfInt(7/2))
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J, M), J=HalfInt(1):HalfInt(4), M) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(2, 0, 3/2, 1/2, J,-M), J=HalfInt(1):HalfInt(4), M) ≈ 0 atol=1e-15
    @test @racahsum((-1)^Integer(5/2-m) * clebschgordan(5/2, m, 5/2, -m, 0, 0), m) ≈ √6
    @test @racahsum((-1)^Integer(5/2-m) * clebschgordan(5/2, m, 5/2, -m, 2, 0), m) ≈ 0 atol=1e-15
    # Orthogonality relations for Clebsch–Gordan coefficients
    @test @racahsum(clebschgordan(2, 2, 3/2, -3/2, J, M)*clebschgordan(2, 2, 3/2, -3/2, J, M), J, M) ≈ 1
    @test @racahsum(clebschgordan(2, 1, 3/2, -1/2, J, M)*clebschgordan(2, 1, 3/2, -1/2, J, M), J, M) ≈ 1
    @test @racahsum(clebschgordan(2, 2, 3/2, -3/2, J, M)*clebschgordan(2, 1, 3/2, -1/2, J, M), J, M) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(2, 2, 3/2, -3/2, J,-M)*clebschgordan(2, 2, 3/2, -3/2, J,-M), J, M) ≈ 1
    @test @racahsum(clebschgordan(2, 1, 3/2, -1/2, J,-M)*clebschgordan(2, 1, 3/2, -1/2, J,-M), J, M) ≈ 1
    @test @racahsum(clebschgordan(2, 2, 3/2, -3/2, J,-M)*clebschgordan(2, 1, 3/2, -1/2, J,-M), J, M) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(5/2, m₁, 3/2, m₂, 2, 0)*clebschgordan(5/2, m₁, 3/2, m₂, 2, 0), m₁, m₂) ≈ 1
    @test @racahsum(clebschgordan(5/2, m₁, 3/2, m₂, 3, 0)*clebschgordan(5/2, m₁, 3/2, m₂, 3, 0), m₁, m₂) ≈ 1
    @test @racahsum(clebschgordan(5/2, m₁, 3/2, m₂, 3, 1)*clebschgordan(5/2, m₁, 3/2, m₂, 3, 1), m₁, m₂) ≈ 1
    @test @racahsum(clebschgordan(5/2, m₁, 3/2, m₂, 2, 0)*clebschgordan(5/2, m₁, 3/2, m₂, 3, 0), m₁, m₂) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(5/2, m₁, 3/2, m₂, 3, 0)*clebschgordan(5/2, m₁, 3/2, m₂, 3, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(5/2, m₁, 3/2, m₂, 2, 0)*clebschgordan(5/2, m₁, 3/2, m₂, 3, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(5/2,-m₁, 3/2, m₂, 2, 0)*clebschgordan(5/2,-m₁, 3/2, m₂, 2, 0), m₁, m₂) ≈ 1
    @test @racahsum(clebschgordan(5/2, m₁, 3/2,-m₂, 3, 0)*clebschgordan(5/2, m₁, 3/2,-m₂, 3, 0), m₁, m₂) ≈ 1
    @test @racahsum(clebschgordan(5/2,-m₁, 3/2,-m₂, 3, 1)*clebschgordan(5/2,-m₁, 3/2,-m₂, 3, 1), m₁, m₂) ≈ 1
    @test @racahsum(clebschgordan(5/2,-m₁, 3/2, m₂, 2, 0)*clebschgordan(5/2,-m₁, 3/2, m₂, 3, 0), m₁, m₂) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(5/2, m₁, 3/2,-m₂, 3, 0)*clebschgordan(5/2, m₁, 3/2,-m₂, 3, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test @racahsum(clebschgordan(5/2,-m₁, 3/2,-m₂, 2, 0)*clebschgordan(5/2,-m₁, 3/2,-m₂, 3, 1), m₁, m₂) ≈ 0 atol=1e-15
    # 3-j symbol
    @test @racahsum(wigner3j(2, j, 4, 0, 0, 0), j) ≈ sum(wigner3j(2, j, 4, 0, 0, 0) for j=2:2:6)
    @test @racahsum(wigner3j(2, j, 4, m, 0, m′), j, m, m′) ≈ sum(wigner3j(2, j, 4, m, 0, -m) for j=2:6, m=-2:2)
    @test @racahsum(wigner3j(2, j, 4,-m, 0, m′), j, m, m′) ≈ sum(wigner3j(2, j, 4, m, 0, -m) for j=2:6, m=-2:2)
    @test @racahsum(wigner3j(2, j, 4, m, 0,-m′), j, m, m′) ≈ sum(wigner3j(2, j, 4, m, 0, -m) for j=2:6, m=-2:2)
    @test @racahsum(wigner3j(2, j, 4,-m, 0,-m′), j, m, m′) ≈ sum(wigner3j(2, j, 4, m, 0, -m) for j=2:6, m=-2:2)
    @test @racahsum((-1)^Integer(3/2-m)*wigner3j(3/2, 3/2, 0,  m, -m, 0), m) ≈ 2
    @test @racahsum((-1)^Integer(3/2+m)*wigner3j(3/2, 3/2, 0, -m,  m, 0), m) ≈ 2
    @test @racahsum((-1)^Integer(3/2-m)*wigner3j(3/2, 3/2, 1,  m, -m, 0), m) ≈ 0 atol=1e-15
    @test @racahsum((-1)^Integer(3/2+m)*wigner3j(3/2, 3/2, 1, -m,  m, 0), m) ≈ 0 atol=1e-15
    @test @racahsum((-1)^Integer(4-m)*wigner3j(4, 4, 0,  m, -m, 0), m) ≈ 3
    @test @racahsum((-1)^Integer(4+m)*wigner3j(4, 4, 0, -m,  m, 0), m) ≈ 3
    @test @racahsum((-1)^Integer(4-m)*wigner3j(4, 4, 1,  m, -m, 0), m) ≈ 0 atol=1e-15
    @test @racahsum((-1)^Integer(4+m)*wigner3j(4, 4, 1, -m,  m, 0), m) ≈ 0 atol=1e-15
    # Orthogonality relations for 3-j symbols
    @test 9 * @racahsum(wigner3j(2, 6, 4,  m₁,  m₂, 0)*wigner3j(2, 6, 4,  m₁,  m₂, 0), m₁, m₂) ≈ 1
    @test 9 * @racahsum(wigner3j(2, 6, 4, -m₁,  m₂, 0)*wigner3j(2, 6, 4, -m₁,  m₂, 0), m₁, m₂) ≈ 1
    @test 9 * @racahsum(wigner3j(2, 6, 4,  m₁, -m₂, 0)*wigner3j(2, 6, 4,  m₁, -m₂, 0), m₁, m₂) ≈ 1
    @test 9 * @racahsum(wigner3j(2, 6, 4, -m₁, -m₂, 0)*wigner3j(2, 6, 4, -m₁, -m₂, 0), m₁, m₂) ≈ 1
    @test 9 * @racahsum(wigner3j(2, 6, 4,  m₁,  m₂, 0)*wigner3j(2, 6, 5,  m₁,  m₂, 0), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4, -m₁,  m₂, 0)*wigner3j(2, 6, 5, -m₁,  m₂, 0), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4,  m₁, -m₂, 0)*wigner3j(2, 6, 5,  m₁, -m₂, 0), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4, -m₁, -m₂, 0)*wigner3j(2, 6, 5, -m₁, -m₂, 0), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4,  m₁,  m₂, 0)*wigner3j(2, 6, 4,  m₁,  m₂, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4, -m₁,  m₂, 0)*wigner3j(2, 6, 4, -m₁,  m₂, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4,  m₁, -m₂, 0)*wigner3j(2, 6, 4,  m₁, -m₂, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test 9 * @racahsum(wigner3j(2, 6, 4, -m₁, -m₂, 0)*wigner3j(2, 6, 4, -m₁, -m₂, 1), m₁, m₂) ≈ 0 atol=1e-15
    @test 6 * @racahsum(wigner3j(3/2, 1, 5/2, m₁, m₂, -1/2)*wigner3j(3/2, 1, 5/2, m₁, m₂, -1/2), m₁, m₂) ≈ 1
    @test 4 * @racahsum(wigner3j(3/2, 1, 3/2, m₁, m₂, -1/2)*wigner3j(3/2, 1, 3/2, m₁, m₂, -1/2), m₁, m₂) ≈ 1
    @test 6 * @racahsum(wigner3j(3/2, 1, 5/2, m₁, m₂, -1/2)*wigner3j(3/2, 1, 3/2, m₁, m₂, -1/2), m₁, m₂) ≈ 0 atol=1e-15
    @test 6 * @racahsum(wigner3j(3/2, 1, 5/2, m₁, m₂, -1/2)*wigner3j(3/2, 1, 5/2, m₁, m₂,  1/2), m₁, m₂) ≈ 0 atol=1e-15
    @test 4 * @racahsum(wigner3j(3/2, 1, 3/2, m₁, m₂, -1/2)*wigner3j(3/2, 1, 3/2, m₁, m₂, -3/2), m₁, m₂) ≈ 0 atol=1e-15
    @test @racahsum((2j₃+1)*wigner3j(2, 6, j₃, 0, 0, m₃)*wigner3j(2, 6, j₃, 0, 0, m₃), j₃, m₃) ≈ 1
    @test @racahsum((2j₃+1)*wigner3j(2, 6, j₃, 0, 0, m₃)*wigner3j(2, 6, j₃, 0, 1, m₃), j₃, m₃) ≈ 0 atol=1e-15
    @test @racahsum((2j₃+1)*wigner3j(2, 6, j₃, 0, 0, m₃)*wigner3j(2, 6, j₃, 1, 0, m₃), j₃, m₃) ≈ 0 atol=1e-15
    @test @racahsum((2j₃+1)*wigner3j(3/2, 1, j₃, 1/2, 0, m₃)*wigner3j(3/2, 1, j₃,  1/2, 0, m₃), j₃, m₃) ≈ 1
    @test @racahsum((2j₃+1)*wigner3j(3/2, 1, j₃, 1/2, 0, m₃)*wigner3j(3/2, 1, j₃,  1/2, 1, m₃), j₃, m₃) ≈ 0 atol=1e-15
    @test @racahsum((2j₃+1)*wigner3j(3/2, 1, j₃, 1/2, 0, m₃)*wigner3j(3/2, 1, j₃, -1/2, 0, m₃), j₃, m₃) ≈ 0 atol=1e-15
    # 6-j symbol as sum over 3-j symbols
    @test @racahsum((-1)^Integer(18-m₁-m₂-m₃-m₄-m₅-m₆) * wigner3j(3, 3, 3, -m₁, -m₂, -m₃) *
                    wigner3j(3, 3, 3, m₁, -m₅, m₆) * wigner3j(3, 3, 3, m₄, m₂, -m₆) *
                    wigner3j(3, 3, 3, -m₄, m₅, m₃), m₁, m₂, m₃, m₄, m₅, m₆) ≈ wigner6j(3, 3, 3, 3, 3, 3)
    @test @racahsum((-1)^Integer(30-m₁-m₂-m₃-m₄-m₅-m₆) * wigner3j(5, 5, 5, -m₁, -m₂, -m₃) *
                    wigner3j(5, 5, 5, m₁, -m₅, m₆) * wigner3j(5, 5, 5, m₄, m₂, -m₆) *
                    wigner3j(5, 5, 5, -m₄, m₅, m₃), m₁, m₂, m₃, m₄, m₅, m₆) ≈ wigner6j(5, 5, 5, 5, 5, 5)
    @test @racahsum((-1)^Integer(9-m₁-m₂-m₃-m₄-m₅-m₆) * wigner3j(1, 3/2, 5/2, -m₁, -m₂, -m₃) *
                    wigner3j(1, 5/2, 3/2, m₁, -m₅, m₆) * wigner3j(0, 3/2, 3/2, m₄, m₂, -m₆) *
                    wigner3j(0, 5/2, 5/2, -m₄, m₅, m₃), m₁, m₂, m₃, m₄, m₅, m₆) ≈ wigner6j(1, 3/2, 5/2, 0, 5/2, 3/2)
    # 6-j symbol
    @test @racahsum(wigner6j(1,1,1,1,1,j), j) ≈ wigner6j(1,1,1,1,1,0) + wigner6j(1,1,1,1,1,1) + wigner6j(1,1,1,1,1,2)
    @test @racahsum(j*wigner6j(1,1,1,1,1,j), j) ≈ wigner6j(1,1,1,1,1,1) + 2wigner6j(1,1,1,1,1,2)
    @test @racahsum(wigner6j(1,3/2,j,0,5/2,3/2), j) ≈ wigner6j(1,3/2,5/2,0,5/2,3/2)
    # Orthogonality relation for 6-j symbols
    @test @racahsum((2j₃+1)*wigner6j(1,3/2,j₃,0,5/2,3/2)*wigner6j(1,3/2,j₃,0,5/2,3/2), j₃) ≈ 1/4
    @test @racahsum((2j₃+1)*wigner6j(1,3/2,j₃,0,5/2,3/2)*wigner6j(1,3/2,j₃,0,5/2,5/2), j₃) ≈ 0 atol=1e-15
    @test @racahsum((2j₃+1)*wigner6j(1,1,j₃,1,1,1)*wigner6j(1,1,j₃,1,1,0), j₃) ≈ 0 atol=1e-15
    @test @racahsum((2j₃+1)*wigner6j(1,1,j₃,1,1,1)*wigner6j(1,1,j₃,1,1,1), j₃) ≈ 1/3
    @test @racahsum((2j₃+1)*wigner6j(1,1,j₃,1,1,1)*wigner6j(1,1,j₃,1,1,2), j₃) ≈ 0 atol=1e-15
    # 9-j symbol as sum over 6-j symbols
    @test @racahsum((-1)^2j*(2j+1)*wigner6j(1,1,1,1,0,j)*wigner6j(1,1,1,1,j,1)*wigner6j(1,1,0,j,1,1), j) ≈
        wigner9j(1,1,1,1,1,1,1,1,0)
    @test @racahsum((-1)^2j*(2j+1)*wigner6j(1,5/2,3/2,3/2,0,j)*wigner6j(3/2,0,3/2,5/2,j,5/2)*
                    wigner6j(5/2,5/2,0,j,1,3/2), j) ≈ wigner9j(1,3/2,5/2,5/2,0,5/2,3/2,3/2,0)
    @test @racahsum((-1)^2j*(2j+1)*wigner6j(1,5/2,3/2,3/2,0,j)*wigner6j(3/2,0,3/2,5/2,j,5/2)*
                    wigner6j(5/2,5/2,0,j,1,3/2), j) ≈ wigner9j(1,3/2,5/2,5/2,0,5/2,3/2,3/2,0)
    @test @racahsum((-1)^2j*(2j+1)*wigner6j(1/2,3,5/2,5/2,1,j)*wigner6j(1/2,3,5/2,3,j,2)*
                    wigner6j(1,2,1,j,1/2,1/2), j) ≈ wigner9j(1/2,1/2,1,3,3,2,5/2,5/2,1)
    # 9-j symbol
    @test @racahsum(wigner9j(1,1,1,1,1,1,1,1,j), j) ≈ sum(wigner9j(1,1,1,1,1,1,1,1,j) for j=0:2)
    @test @racahsum(j^2*wigner9j(1,1,1,1,1,1,1,1,j), j) ≈ wigner9j(1,1,1,1,1,1,1,1,1) + 4wigner9j(1,1,1,1,1,1,1,1,2)
    @test @racahsum(wigner9j(1,3/2,5/2,5/2,j,5/2,3/2,3/2,0), j) ≈ sum(wigner9j(1,3/2,5/2,5/2,j,5/2,3/2,3/2,0) for j=0:3)
    @test @racahsum(wigner9j(1/2,1/2,1,j,3,j,5/2,5/2,1), j) ≈ wigner9j(1/2,1/2,1,2,3,2,5/2,5/2,1)
    # Orthogonality relation for 9-j symbols
    @test @racahsum((2j₇+1)*(2j₈+1)*wigner9j(1/2,1/2,1,3,3,2,j₇,j₈,1)*wigner9j(1/2,1/2,1,3,3,2,j₇,j₈,1), j₇, j₈) ≈ 1/15
    @test @racahsum((2j₇+1)*(2j₈+1)*wigner9j(1/2,1/2,0,3,3,1,j₇,j₈,1)*wigner9j(1/2,1/2,0,3,3,1,j₇,j₈,1), j₇, j₈) ≈ 1/3
    @test @racahsum((2j₇+1)*(2j₈+1)*wigner9j(1/2,1/2,0,3,3,1,j₇,j₈,1)*wigner9j(1/2,1/2,1,3,3,2,j₇,j₈,1), j₇, j₈) ≈ 0 atol=1e-15
    # If summation bounds cannot be determined and were not specified explicitly, sum over 0:1/2:30
    @test @racahsum(wigner6j(j,0,j,0,j,0), j) ≈ sum(wigner6j(j,0,j,0,j,0) for j=0:HalfInt(1/2):30)
    @test @racahsum(wigner6j(j,0,j,0,j,0), j=0:10) ≈ sum(wigner6j(j,0,j,0,j,0) for j=0:10)
end
