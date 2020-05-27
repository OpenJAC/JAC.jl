
#
println("Ch) Test of the AlphaVariation module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-AlphaVariation.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                        properties=[AlphaX()], 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        alphaSettings=AlphaVariation.Settings(true, true, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


