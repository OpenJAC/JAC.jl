
#
println("Ck) Test of the MultipolePolarizibility module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-MultipolePolarizibility.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Polarity], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        polaritySettings=MultipolePolarizibility.Settings(EmMultipole[], 0, 0, Float64[], false, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


