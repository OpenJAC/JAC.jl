#
println("Fa) Symbolic evaluation by means of special values & recursion relations of the Wigner n-j symbols.")
#
using SymEngine

if  false
    # Test for special values of the Wigner 3-j symbols
    w3j = RacahAlgebra.selectW3j(17);                    println("w3j-original      = $w3j")
    rex = RacahAlgebra.equivalentForm(w3j);             println("rex-equivalent    = $rex")
    wa  = RacahAlgebra.evaluate(rex, special=true);     println("\nwa-special value  = $wa")
elseif  true
    # Test for special values of the Wigner 6-j symbols
    w6j = RacahAlgebra.selectW6j(17);                   println("w6j-original      = $w6j")
    rex = RacahAlgebra.equivalentForm(w6j);             println("rex-equivalent    = $rex")
    wa  = RacahAlgebra.evaluate(rex, special=true);     println("\nwa-special value  = $wa")
    ## ws  = [:ja => 2//1, :jb => 2, :jc => 1]
    ##wb  = RacahAlgebra.subs(w6j, ws);                   println("\nw6j-substitution  = $wb")
    ##wc  = RacahAlgebra.subs(wa, ws);                    println("\nrhs-substitution  = $wc")
    ##wd  = RacahAlgebra.evaluateNumerical(wb);           println("\nlhs-numerical     = $wd")
    ##we  = RacahAlgebra.evaluateNumerical(wc);           println("\nrhs-numerical     = $we")
elseif  false
    # Test for special values of the Wigner 9-j symbols
    w9j = RacahAlgebra.selectW9j(1);                    println("w9j-original      = $w9j")
    rex = RacahAlgebra.equivalentForm(w9j);             println("rex-equivalent    = $rex")
    wa  = RacahAlgebra.evaluate(rex, special=true);     println("\nwa-special value  = $wa")
elseif  false
    # Test for recursion relations of the Wigner 3-j symbols
    # RecursionW3jMagnetic(), RecursionW3jOneStep(), RecursionW3jHalfStep(), RecursionW3jLouck()
    w3j = RacahAlgebra.W9j(1,2,3, -1,-2, 3);                            println("w3j-original      = $w3j")
    wa  = RacahAlgebra.recursionW3j(w3j, RecursionW3jMagnetic());       println("wa-special value  = $wa")
    wb  = RacahAlgebra.evaluateNumerical(w3j);                          println("lhs-numerical     = $wb")
    wc  = RacahAlgebra.evaluateNumerical(wa);                           println("rhs-numerical     = $wc")  ## Array{RacahExpression,1}
end
