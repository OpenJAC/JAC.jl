
println("Ce) Test of the FormFactor module with ASF from an internally generated multiplet.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-FormFactor.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings=[ FormFactor.Settings([0., 0.1, 1.0, 10., 100., 1000.], true, LevelSelection())] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


