#
println("De) Test of the AutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-AutoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  true
    # K-LL Auger spectrum of neon: Comparison with PhD and related work
    augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(1,0)]), 0., 1.0e6, 4, CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "a.u.")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  =[Configuration("1s 2s^2 2p^6")],
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                            processSettings = augerSettings )

    wb = perform(wa)
    
elseif  false
    # K-LL Auger spectrum of Be-like neon: Comparison with Bruch (PRA, 1991)
    ## augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(1,0)]), 0., 1.0e6, 4, CoulombInteraction())
    augerSettings = AutoIonization.Settings(false, true, LineSelection(), 0., 1.0e6, 4, CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  =[Configuration("1s 2s^2 2p")],
                            finalConfigs    =[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
elseif  false
    # K-LL Auger spectrum of Ne-like aluminium: Comparison with Fan et al. (PRA, 2018)
    ## augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(1,0)]), 0., 1.0e6, 4, CoulombInteraction())
    augerSettings = AutoIonization.Settings(false, true, LineSelection(), 0., 1.0e6, 4, CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(13.), 
                            ## initialConfigs  =[Configuration("1s 2s^2 2p^6")],
                            ## finalConfigs    =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                            initialConfigs  = [Configuration("1s 2s 2p^6"), Configuration("2s^2 2p^6")],
                            finalConfigs    = [Configuration("1s 2s^2 2p^4"), Configuration("1s 2s 2p^5"), Configuration("1s 2p^6")], 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
elseif  true
    # K-LL Auger spectrum of Be-like aluminium: Comparison with Fan et al. (PRA, 2018)
    ## augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(1,0)]), 0., 1.0e6, 4, CoulombInteraction())
    augerSettings = AutoIonization.Settings(false, true, LineSelection(), 0., 1.0e6, 4, CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(13.), 
                            initialConfigs  = [Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs    = [Configuration("1s^2")], 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
elseif  false
    # Resonant L_23 - M_23 M_23 Auger spectrum of argon: Comparison with Chen (PRA, 1991)
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())
    augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(2,0), (4,0)]), 0., 1.0e6, 4, CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(18.), 
                            initialConfigs  =[Configuration("1s^2 2s^2 2p^5 3s^2 3p^6 4s")], initialAsfSettings=asfSettings,
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^4 4s")], finalAsfSettings=asfSettings, 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
elseif  false
    # Resonant M_45 - N_23 N_23 Auger spectrum of krypton: Comparison with Chen (PRA, 1991)
    augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(4,0), (9,0), (11,0)]), 0., 1.0e6, 4, CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.2e-2, rbox = 15.0)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(36.), 
                            initialConfigs  =[Configuration("[Ar] 3d^9 4s^2 4p^6 5p")],
                            finalConfigs    =[Configuration("[Ar] 3d^10 4s^2 4p^4 5p")], 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
end
setDefaults("print summary: close", "")


