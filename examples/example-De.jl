#
println("De) Test of the AutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-AutoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
setDefaults("method: normalization, pure Coulomb")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  false
    augerSettings = AutoIonization.Settings(true, true, true, Tuple{Int64,Int64}[(1,0)], 0., 1.0e6, 2, "Coulomb")
    grid          = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.3e-2, NoPoints = 800)
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(12.), 
                            initialConfigs  =[Configuration("1s 2s^2 2p^6")],
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
else
    augerSettings = AutoIonization.Settings(true, true, true, Tuple{Int64,Int64}[(2,0), (4,0)], 0., 1.0e6, 4, "Coulomb")
    grid          = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.1e-2, NoPoints = 1100)
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(18.), 
                            initialConfigs  =[Configuration("1s^2 2s^2 2p^5 3s^2 3p^6 4s")],
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^4 4s")], 
                            process = Auger(),  processSettings = augerSettings )

    wb = perform(wa)
    
end
setDefaults("print summary: close", "")


