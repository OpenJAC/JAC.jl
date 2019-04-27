
#
println("Cj) Test of the GreenFunction module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-GreenFunction.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Green], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        greenSettings=GreenFunction.Settings(true, true, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


