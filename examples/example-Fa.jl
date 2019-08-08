#
println("Fa) Symbolic evaluation by means of special values & recursion relations of the Wigner n-j symbols.")
#
using SymEngine

if  false
    # Test for special values of the Wigner 3-j symbols
    w3j = RacahAlgebra.selectW3j(9);                    println("\nw3j-original      = $w3j")
    ## rex = RacahAlgebra.equivalentForm(w3j);             println("\nrex-equivalent    = $rex")
    ## w3j = rex.w3js[1];                                  println("\nw3j-without phase = $w3j")
    wa  = RacahAlgebra.evaluate(w3j);                   println("\nwa-special value  = $wa")
    ws  = [:ja => 2//1, :ma => 1, :jb => 2, :jc => 1]
    wb  = RacahAlgebra.subs(w3j, ws);                   println("\nw3j-substitution  = $wb")
    wc  = RacahAlgebra.subs(wa, ws);                    println("\nrhs-substitution  = $wc")
    wd  = RacahAlgebra.evaluateNumerical(wb);           println("\nlhs-numerical     = $wd")
    we  = RacahAlgebra.evaluateNumerical(wc);           println("\nrhs-numerical     = $we")
elseif  true
    # Test for special values of the Wigner 6-j symbols
    w6j = RacahAlgebra.selectW6j(2);                    println("\nw6j-original      = $w6j")
    ## rex = RacahAlgebra.equivalentForm(w6j);             println("\nrex-equivalent    = $rex")
    ## w6j = rex.w6js[1];                                  println("\nw6j-without phase = $w6j")
    wa  = RacahAlgebra.evaluate(w6j);                   println("\nwa-special value  = $wa")
    ws  = [:ja => 2//1, :jb => 2, :jc => 1]
    wb  = RacahAlgebra.subs(w6j, ws);                   println("\nw6j-substitution  = $wb")
    wc  = RacahAlgebra.subs(wa, ws);                    println("\nrhs-substitution  = $wc")
    wd  = RacahAlgebra.evaluateNumerical(wb);           println("\nlhs-numerical     = $wd")
    we  = RacahAlgebra.evaluateNumerical(wc);           println("\nrhs-numerical     = $we")
elseif  false
    # Test for special values of the Wigner 9-j symbols
    w9j = RacahAlgebra.selectW9j(1);                    println("w9j-original      = $w9j")
    rex = RacahAlgebra.equivalentForm(w9j);             println("rex-equivalent    = $rex")
    w9j = rex.w9js[1];                                  println("w9j-without phase = $w9j")
    wa  = RacahAlgebra.evaluate(w9j);                   println("wa-special value  = $wa")
    ## wb  = SymEngine.subs(  ; w9j);                      println("w9j-substitution  = $wb")
    ## wc  = SymEngine.subs(  ; wa);                       println("rhs-numerical     = $wc")
    wd  = RacahAlgebra.evaluateNumerical(wb);           println("lhs-numerical     = $wd")
elseif  false
    # Test for recursion relations of the Wigner 3-j symbols
    # RecursionW3jMagnetic(), RecursionW3jOneStep(), RecursionW3jHalfStep(), RecursionW3jLouck()
    w3j = RacahAlgebra.W9j(1,2,3, -1,-2, 3);                            println("w3j-original      = $w3j")
    wa  = RacahAlgebra.recursionW3j(w3j, RecursionW3jMagnetic());       println("wa-special value  = $wa")
    wb  = RacahAlgebra.evaluateNumerical(w3j);                          println("lhs-numerical     = $wb")
    wc  = RacahAlgebra.evaluateNumerical(wa);                           println("rhs-numerical     = $wc")  ## Array{RacahExpression,1}
elseif  false 
    # Some earlier tests
    w3j = W3j(:a, :b, :c, :ma, :mb, :mc);    println("$w3j")
    w6j = W6j(:a, :b, :c, :d, :we, :f);      println("$w6j")

    wa  = RacahAlgebra.symmetricForms(w3j)

    xw3j = W3j(:ja, :ja, 0, :ma, -Basic(:ma), 0)
    yw3j = RacahAlgebra.equivalentForm(xw3j)
    wa   = RacahAlgebra.evaluate(yw3j.w3js[1])

    ma =  Basic(:ma)
    wa = [Basic(:ma), Basic(:mb), Basic(:mc), Basic(:md)]
end
