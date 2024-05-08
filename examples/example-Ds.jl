#
println("Ds) Test of the CoulombExcitation module with ASF from an internally generated initial- and final-state multiplet.")


if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-CoulombExcitation.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(36.), 
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                            processSettings=CoulombExcitation.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [3000., 4000.], false, true, 
                                            LineSelection()) )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end

