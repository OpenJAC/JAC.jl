#
println("Fb) Symbolic simplification of Racah expressions by means of sum rules.")
#
using SymEngine

if  true
    # Test for sum rules
    rex = RacahAlgebra.selectRacahExpression(41);        println("rex-original        = $rex")
    # rex = RacahAlgebra.equivalentForm(rex);              println("rex-equivalent      = $rex")
    rex = RacahAlgebra.evaluate(rex);                    println("wa-sum rule eval    = $rex")
    ## rex = RacahAlgebra.evaluate(rex);                    println("wb-sum rule eval    = $rex")
    ## ws  = [:j => 2//1, :m => 1, :jb => 2, :jc => 1]
    ## wb  = RacahAlgebra.subs(rex, ws);                   println("\nrex-substitution  = $wb")
    ## wc  = RacahAlgebra.subs(wa, ws);                    println("\nrhs-substitution  = $wc")
    ## wd  = RacahAlgebra.evaluateNumerical(wb);           println("lhs-numerical       = $wd")
    ## we  = RacahAlgebra.evaluateNumerical(wc);           println("rhs-numerical       = $we")
end
