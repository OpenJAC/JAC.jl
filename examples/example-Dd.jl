
#
println("Dd) Test of the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoRecombination.sum")
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, Galerkin") setDefaults("method: continuum, asymptotic Coulomb") 
setDefaults("method: normalization, pure Coulomb")    ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

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
    grid = Radial.Grid(Radial.Grid(true), rnt = 1.0e-6,h = 8.0e-3, hp = 4.0e-3, NoPoints = 2600)

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(92.), 
                            initialConfigs=[Configuration("1s")],
                            finalConfigs  =[Configuration("1s^2")], 
                            process = Rec(), 
                            processSettings=PhotoRecombination.Settings([E1, M1, E2, M2], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 30., 50.], 
                                                                        [2.18, 21.8, 218.0], true, true, true, true, false, Tuple{Int64,Int64}[(1,1)]) )

    wb = perform(wa)

    
elseif  true
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.0e-2, NoPoints = 900)

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(12.), 
                            initialConfigs=[Configuration("1s^2")],
                            finalConfigs  =[Configuration("1s^2 2s"), Configuration("1s^2 3s"), Configuration("1s^2 3p"), Configuration("1s^2 3d")], 
                            process = Rec(), 
                            processSettings=PhotoRecombination.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10.], [2.18, 21.8, 218.0], 
                                                                        false, false, false, true, LineSelection(false, indexPairs=[(1,0)]) ) )

    wb = perform(wa)

end
setDefaults("print summary: close", "")

