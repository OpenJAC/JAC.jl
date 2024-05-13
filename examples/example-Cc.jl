
println("Cc) Apply & test the IsotopeShift module with ASF from an internally generated multiplet.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-IsotopeShift.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), 
                            nuclearModel=Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0),
                            configs=[Configuration("[Ne] 3s^2"), Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3p^2")],
                            propertySettings=[ IsotopeShift.Settings(true, true, true, true, false, 0., LevelSelection() )] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end
