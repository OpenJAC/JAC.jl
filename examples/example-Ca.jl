#
println("Ca) Test of the Einstein module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-Einstein.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.EinsteinX], 
                        configs=[Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d"), Configuration("[Ne] 3s^2 3p^5")],
                        einsteinSettings=Einstein.Settings([E1, M1, E2], true, false, Tuple{Int64,Int64}[], 0., 0., 10000. ) )

wb = perform(wa; output=true)
setDefaults("print summary: close", "")
