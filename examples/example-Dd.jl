
#
println("Dd) Test of the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoRecombination.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  false
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.0e-2, NoPoints = 900)

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(12.), 
                            initialConfigs=[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^6")], 
                            process = Rec(), 
                            processSettings=PhotoRecombination.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 30., 50.], [0.], 
                                                                        false, true, true, true, true, Tuple{Int64,Int64}[(1,1)]) )

    wb = perform(wa)
    
elseif  false
    # REC calculations for the capture into 1s, 2s, 2p, 3s, 3p, 3d and Z=5: Comparison with Eichler (ADNDT, 2000)
    setDefaults("unit: energy", "eV") # then, give photonEnergies in eV
    setDefaults("unit: cross section", "barn")   
    
    recSettings = PhotoRecombination.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [1., 10., 100., 1000., 10000.], 
                                              Float64[], false, false, false, true, LineSelection())
                                                                        
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(9.), 
                            ## initialConfigs = [Configuration("1s^2")],
                            ## finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p"), 
                            ##                   Configuration("1s^2 3s"), Configuration("1s^2 3p"), Configuration("1s^2 3d")], 
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("1s^2")], 
                            process = Rec(), processSettings=recSettings )

    wb = perform(wa)

elseif  true
    # REC calculations for the capture into 2s at given ion energies and Z=92: Comparison with Fritzsche (PRA, 2005)
    setDefaults("unit: energy", "eV") # then, give photonEnergies in eV
    setDefaults("unit: cross section", "barn")   
    
    recSettings = PhotoRecombination.Settings([E1, M1, E2, M2, E3, M3], [JAC.UseBabushkin], Float64[], 
                                              [2.18, 21.8, 218.0], true, false, false, true, LineSelection())
                                                                        
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.3e-2, rbox = 10.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(92.), 
                            ## initialConfigs = [Configuration("1s^2")],
                            ## finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p"), 
                            ##                   Configuration("1s^2 3s"), Configuration("1s^2 3p"), Configuration("1s^2 3d")], 
                            initialConfigs = [Configuration("1s^2")],
                            finalConfigs   = [Configuration("1s^2 2s")], 
                            process = Rec(), processSettings=recSettings )

    wb = perform(wa)

end
setDefaults("print summary: close", "")

