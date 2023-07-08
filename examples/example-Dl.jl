#
println("Dl) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoDoubleIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if true
    # Single-photon 1s^2 photoionization of helium
    activeShells   = Basics.generateShellList(4, 4, 1)
    doubleSettings = PhotoDoubleIonization.Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], 
                                                    [3., 4.],   ## , 4.5, 5., 6., 7., 8., 9.], 
                                                    activeShells, 1.0, 2, true, 2, LineSelection(true, indexPairs=[(1,0)]) )
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(2.), 
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s^0")], 
                            processSettings = doubleSettings )

    wb = perform(wa)
    #
elseif true
    # Single-photon 2p^2 photoionization of neon
    doubleSettings   = PhotoDoubleIonization.Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], [10., 20., 30.], doubleGreen, 2, 
                                                      false, true, 2, LineSelection(true, indexPairs=[(1,2)]))
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(5.), 
                            initialConfigs  =[Configuration("1s^2 2s^2")],
                            finalConfigs    =[Configuration("1s 2s")], 
                            processSettings = doubleSettings )

    wb = perform(wa)
    
end
setDefaults("print summary: close", "")


