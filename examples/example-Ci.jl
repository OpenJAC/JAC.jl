
#
println("Ci) Test of the DecayYield module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-DecayYield.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Yields], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        yieldSettings=DecayYield.Settings(true, true, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


