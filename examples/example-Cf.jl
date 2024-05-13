
println("Cf) Apply & test the reduced 1- and 2-particle density matrices & natural orbitals.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-ReducedDensityMatrix.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("[Ne] 3s^2 3p^5")],
                            propertySettings=[ ReducedDensityMatrix.Settings(true, true, true, true, true, LevelSelection(true, indices=[(2)])) ] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


