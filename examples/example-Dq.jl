
#
println("Dq)  Test of the TwoElectronOnePhoton module with ASF from an internally generated initial and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-TwoElectronOnePhoton.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")


grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

if  true
    # Two-electron-one-photon transition of Li-like Ar: 1s 2p 3d  -->  1s^2 2p: Generation of GreenFunction
    rydbergShells = Basics.generateShellList(6,  6, "p") 
    teopSettings  = TwoElectronOnePhoton.Settings(5.0, 0.1, LineSelection(), CoulombInteraction())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(18.), 
                            initialConfigs = [Configuration("1s 2p 3d")],
                            finalConfigs   = [Configuration("1s^2 2s")],  
                            processSettings= teopSettings )

    wb = perform(wa)
elseif  true
    # Two-electron-one-photon transition of Li-like Ar: 1s 2p 3d  -->  1s^2 2p: Computation of rates
    rydbergShells = Basics.generateShellList(6,  6, "p") 
    teopSettings  = TwoElectronOnePhoton.Settings(5.0, 0.1, LineSelection(), CoulombInteraction())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(18.), 
                            initialConfigs = [Configuration("1s 2p 3d")],
                            finalConfigs   = [Configuration("1s^2 2s")],  
                            processSettings= teopSettings )

    wb = perform(wa)
end
    
