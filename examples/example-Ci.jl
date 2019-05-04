
#
println("Ci) Test of the DecayYield module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-DecayYield.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(12.); properties=[JAC.Yields],
                        grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600),
                        configs=[Configuration("1s 2s^2 2p^6")],
                        yieldSettings=DecayYield.Settings("SCA", true, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


