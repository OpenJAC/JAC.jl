
println("Cg) Apply & test the PlasmaShift module with ASF from an internally generated multiplet.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-PlasmaShift.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings=[ PlasmaShift.Settings(PlasmaShift.DebyeHueckel(), 0.25, 0., 0) ] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


