
#
println("Cd) Test of the FormFactor module with ASF from an internally generated multiplet.")
#
JAC.define("print summary: open", "zzz-FormFactor.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.FormF], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        formSettings=FormFactor.Settings([0., 0.1, 1.0, 10., 100., 1000.], true, false, Int64[]) )

wb = perform(wa)
JAC.define("print summary: close", "")


