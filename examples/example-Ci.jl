
#
println("Ci) Test of the DecayYield module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-DecayYield.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(false), nuclearModel=Nuclear.Model(12.), 
                        properties=[JAC.Yields],
                        configs=[Configuration("1s 2s^2 2p^6")],
                        yieldSettings=DecayYield.Settings("SCA", true, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


