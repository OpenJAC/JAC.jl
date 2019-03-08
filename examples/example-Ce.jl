
#
println("Ce) Test of the PlasmaShift module with ASF from an internally generated multiplet.")
#
JAC.define("print summary: open", "zzz-PlasmaShift.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Plasma], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        plasmaSettings=PlasmaShift.Settings() )

wb = perform(wa)
JAC.define("print summary: close", "")


