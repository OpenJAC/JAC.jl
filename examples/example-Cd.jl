
println("Cd) Apply & test the LandeZeeman module with ASF from an internally generated multiplet.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-LandeZeeman.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), 
                            nuclearModel=Nuclear.Model(26., "Fermi", 58., 3.75, AngularJ64(5//2), 1.0, 2.0), 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings=[ LandeZeeman.Settings(true, true, true, true, 0., true, LevelSelection()) ] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


