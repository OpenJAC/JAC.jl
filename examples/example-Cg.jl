
#
println("Cg) Test of the reduced 1- and 2-particle density matrices & natural orbitals.")
#
setDefaults("print summary: open", "zzz-ReducedDensityMatrix.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                        configs=[Configuration("[Ne] 3s^2 3p^5")],
                        propertySettings=[ ReducedDensityMatrix.Settings() ] )

wb = perform(wa)
setDefaults("print summary: close", "")


