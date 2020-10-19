#
println("Dl) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoDoubleIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  false
    # Green function for single-photon 2p^2 photoionization of neon
    name             = "Double Auger decay of 1s photoionized neon"
    refConfigs       = [Configuration("1s^2 2s^2 2p^6")]
    levelSymmetries  = [LevelSymmetry(0, Basics.minus), LevelSymmetry(1, Basics.minus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())
    greenRep         = Representation(name, Nuclear.Model(10.), Radial.Grid(true), refConfigs, 
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 3, greenSettings) ) 
    greenOut         = generate(greenRep, output=true)
    doubleGreen      = greenOut["Green channels"]
    
elseif true
    # Single-photon 2p^2 photoionization of neon
    doubleSettings   = PhotoDoubleIonization.Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], [50.0], doubleGreen, 2, 
                                                      false, true, LineSelection(true, indexPairs=[(1,0)]))
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  =[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^4")], 
                            process = PhotoDouble(),  processSettings = doubleSettings )

    wb = perform(wa)
    
    
end
setDefaults("print summary: close", "")


