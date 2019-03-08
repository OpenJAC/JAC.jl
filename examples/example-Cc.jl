
#
println("Cc) Test of the IsotopeShift module with ASF from an internally generated multiplet.")
#
JAC.define("print summary: open", "zzz-IsotopeShift.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Isotope], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        isotopeSettings=IsotopeShift.Settings(true, true, true, true, false, Int64[], "method-1") )

wb = perform(wa)
JAC.define("print summary: close", "")


