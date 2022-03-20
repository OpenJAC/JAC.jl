
println("Db) Test of the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoExcitation.sum")
if false
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(3.), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            finalConfigs  =[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                            processSettings=PhotoExcitation.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, true, true, true, 
                                                                     LineSelection(), 0., 0., 1.0e6, ExpStokes(0., 0., 0.) ) )

    wb = perform(wa)
    setDefaults("print summary: close", "")

elseif true
    # Test of absorption f-values for 1s^2 --> 1s2p ^1P_1 resonance  in the helium isoelectronic sequence; cf. Table 8.1 in Section 8.1.a
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(20.), 
                            initialConfigs=[Configuration("1s^2")],
                            finalConfigs  =[Configuration("1s 2p")], 
                            processSettings=PhotoExcitation.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, true, true, true, 
                                                                     LineSelection(), 0., 0., 1.0e6, ExpStokes(0., 0., 0.) ) )

    wb = perform(wa)
    setDefaults("print summary: close", "")
end
