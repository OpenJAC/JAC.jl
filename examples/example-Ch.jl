
#
println("Ch) Test of the AlphaVariation module with ASF from an internally generated multiplet.")
#
JAC.define("print summary: open", "zzz-AlphaVariation.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.AlphaX], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        alphaSettings=AlphaVariation.Settings(true, true, false, Int64[]) )

wb = perform(wa)
JAC.define("print summary: close", "")


