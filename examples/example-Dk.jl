#
println("Dk) Test of the DoubleAutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-DoubleAutoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  false
    # Green function for double Auger decay of 1s photoionized neon
    name             = "Double Auger decay of 1s photoionized neon"
    refConfigs       = [Configuration("1s 2s^2 2p^6")]
    levelSymmetries  = [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.minus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())
    greenRep         = Representation(name, Nuclear.Model(10.), Radial.Grid(true), refConfigs, 
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 3, greenSettings) ) 
    greenOut         = generate(greenRep, output=true)
    daGreen          = greenOut["Green channels"]
    
elseif true
    # Double Auger decay of 1s photoionized neon
    daSettings = DoubleAutoIonization.Settings(daGreen, 2, true, 1.0, 100., 2, CoulombInteraction(), LineSelection(true, indexPairs=[(1,0)]))
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  = [Configuration("1s 2s^2 2p^6")],
                            finalConfigs    = [Configuration("1s^2 2s^2 2p^3"), Configuration("1s^2 2s 2p^4"), Configuration("1s^2 2p^5")], 
                            process = DoubleAuger(),  processSettings = daSettings )

    wb = perform(wa)
    
    
end
setDefaults("print summary: close", "")


