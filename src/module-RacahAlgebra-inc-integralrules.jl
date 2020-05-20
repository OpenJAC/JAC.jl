
    
    """
    `RacahAlgebra.integralRulesForOneYlm(rex::RacahAlgebra.RacahExpression)`  
        ... attempts to find a simplification of the given Racah expression by using sum rules for one Wigner 3j symbol. 
            Once a simplification is found, no attempt is made to find another simplifcation for this set of rules.
            A (istrue, rex)::Tuple{Bool, RacahExpression} is returned but where rex has no meaning for !istrue.
    """
    function integralRulesForOneYlm(rex::RacahAlgebra.RacahExpression)
        
        error("Not yet implemented.")
        # Loop through all Wigner 3j symbols
        for  (iaW3j, aW3j) in enumerate(rex.w3js)
            aRexList = RacahAlgebra.symmetricForms(aW3j)
            for  xaRex in aRexList
                ww = xaRex.w3js[1]
                #
                #   Rule:                  -m   ( j  j J )         -j
                #   ----        Sum(m) (-1)     (        )  =  (-1)    (2j + 1)  delta(J,0)  delta(M,0)
                #                               ( m -m M )
                #
                specialW3j = W3j(ww.ja, ww.ja, ww.jc, ww.ma, -ww.ma, ww.mc)
                if  ww == specialW3j
                    newPhase  = rex.phase  + xaRex.phase;     testPhase = newPhase   + ww.ma  
                    newWeight = rex.weight * xaRex.weight;    testWeight = newWeight
                    newW3js   = W3j[];       for (ibW3j, bW3j) in enumerate(rex.w3js)  if  iaW3j != ibW3j   push!(newW3js, bW3j)   end   end
                    #
                    if  RacahAlgebra.hasIndexOld(ww.ma, rex.summations)          &&  RacahAlgebra.hasNoVars([ww.ma], testPhase)      &&  
                        RacahAlgebra.hasNoVars([ww.ma], testWeight)           &&  RacahAlgebra.hasNoVars([ww.ma], rex.deltas)     &&  
                        RacahAlgebra.hasNoVars([ww.ma], rex.triangles)        &&  RacahAlgebra.hasNoVars([ww.ma], newW3js)        &&  
                        RacahAlgebra.hasNoVars([ww.ma], rex.w6js)             &&  RacahAlgebra.hasNoVars([ww.ma], rex.w9js)
                        newSummations = RacahAlgebra.removeIndex(ww.ma, rex.summations)
                        newDeltas = rex.deltas;      push!( newDeltas, Kronecker(ww.jc, 0));     push!( newDeltas, Kronecker(ww.mc, 0))
                        wa = RacahExpression( newSummations, rex.integrals, RacahAlgebra.purifyPhase(newPhase + ww.ma - ww.ja), newWeight * (2*ww.ja+1), 
                                              newDeltas, rex.triangles, newW3js, rex.w6js, rex.w9js, rex.ylms, rex.djpqs )
                        println(">> Apply sum rule for one W3j -- Sum(m) (-1)^m ...")
                        return( (true, wa) )
                    end
                end
            end
        end
        
        return( (false, RacahExpression()) )
    end
