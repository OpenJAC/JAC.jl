#
println("Ca) Test of the Einstein module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-Einstein.sum")
grid=JAC.Radial.Grid(true)
## grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 1.0e-2, hp = 1.0e-2, NoPoints = 2000)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), 
                        properties=[JAC.EinsteinX], 
                        configs=[Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d"), Configuration("[Ne] 3s^2 3p^5")],
                        einsteinSettings=Einstein.Settings([E1, M1, E2], true, false, Tuple{Int64,Int64}[], 0., 0., 10000. ) )

wb = perform(wa; output=true)
setDefaults("print summary: close", "")
