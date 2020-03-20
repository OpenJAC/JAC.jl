#
println("Fc) Symbolic simplification of recoupling coefficients by means of sum rules.")
#
using SymEngine
j1  = Basic(:j1);    j2  = Basic(:j2);    j3  = Basic(:j3);    j4  = Basic(:j4);    j5  = Basic(:j5)
j6  = Basic(:j6);    j7  = Basic(:j7);    j8  = Basic(:j8);    j9  = Basic(:j9);    j10 = Basic(:j10)
j11 = Basic(:j11);   j12 = Basic(:j12);   j13 = Basic(:j13);   j14 = Basic(:j14);   j15 = Basic(:j15)
j16 = Basic(:j16);   j17 = Basic(:j17);   j18 = Basic(:j18);   j19 = Basic(:j19);   j20 = Basic(:j20)
j21 = Basic(:j21);   j22 = Basic(:j22);   j23 = Basic(:j23);   j24 = Basic(:j24);   j25 = Basic(:j25)
J   = Basic(:J)

if  false
    #
    leftCsq  = RacahAlgebra.Csq( j1, j2, j3)
    wa       = RacahAlgebra.evaluate(leftCsq, leftCsq)    
elseif false
    #
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( j1, j2, j12), j3, J)
    rightCsq = RacahAlgebra.Csq( j1, RacahAlgebra.Csq( j2, j3, j23), J)
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)    
elseif true
    #
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j2,j5), RacahAlgebra.Csq(j3,j4,j6), j7 )
    rightCsq = RacahAlgebra.Csq( j1, RacahAlgebra.Csq( RacahAlgebra.Csq(j2,j3,j8), j4, j9), j7 )
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)
elseif false
    #
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j2,j9),j3,j10), 
                                 RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j4,j5,j11),j6,j12), RacahAlgebra.Csq(j7,j8,j13),j14), j15)
    rightCsq = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j4,j16), j7, j17), 
                                 RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j2,j5,j18), RacahAlgebra.Csq(j8,j3,j19), j20), j6,j21), j15)
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)
end
