
#
#=
Example: Test of the Dielectronic module      
=#
println("Df) Test of the Dielectronic module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
setDefaults("print summary: open", "zzz-Dielectronic.sum")
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, Galerkin") setDefaults("method: continuum, asymptotic Coulomb") 
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  false
grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.0e-2, NoPoints = 900)
## grid = Radial.Grid(false)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(18.), 
                        initialConfigs=[Configuration("1s^2 2s")],
                        intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],  
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ],  
                        process = JAC.Dierec, 
                        processSettings=Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                              false, Tuple{Int64,Int64,Int64}[(1,5,0), (1,6,0)], 0., 0., 0., "Coulomb")  )

wb = perform(wa)

elseif  false

augerSettings = AutoIonization.Settings(true, true, false, Tuple{Int64,Int64}[(1,0)], 0., 1.0e6, 2, "Coulomb")
grid          = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.3e-2, NoPoints = 800)
grid          = Radial.Grid(false)
    
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), 
                        initialConfigs  =[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p 3p")],
                        finalConfigs    =[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                        process = JAC.Auger,  processSettings = augerSettings )

wb = perform(wa)

elseif  false

grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 1500)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(14.0), 
                        initialConfigs=[Configuration("[Ne] 3s")],
                        intermediateConfigs=[Configuration("[Ne] 3s^2"), Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3s 3d"),
                                             Configuration("[Ne] 3p 4d"), Configuration("[Ne] 3p 4f") ],  
                        finalConfigs  =[Configuration("[Ne] 3s^2"),
                                        Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3s 3d"), Configuration("[Ne] 3s 4d"), Configuration("[Ne] 3s 4f") ],
                        process = JAC.Dierec, 
                        processSettings=Dielectronic.Settings([E1, M1, M2, E2], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                              false, Tuple{Int64,Int64,Int64}[(1,10,1), (1,10,7)], 0., 0., 0., "Coulomb")  )

wb = perform(wa)

elseif  false

refConfigs = [Configuration("[Ne] 3s")]
fromShells = [Shell("3s")]
toShells   = [Shell("3s"), Shell("3p"), Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f"), Shell("5s")]
intermediateConfigs = Basics.generateConfigurations(refConfigs, fromShells, toShells)
grid=JAC.Radial.Grid(true)
#== grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 2000)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.0, "Fermi"), 
                        properties=JAC.AtomicLevelProperty[],
                        configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p^2"), Configuration("1s^2 2s 3s"), 
                                 Configuration("1s^2 2s 3p"), Configuration("1s^2 2s 3d"), Configuration("1s^2 2p 3s"), Configuration("1s^2 2p 3p"), 
                                 Configuration("1s^2 2p 3d")],
                        ## configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), 
                        ##          Configuration("1s^2 7s"), Configuration("1s^2 7p"), Configuration("1s^2 7d"), Configuration("1s^2 7f")], 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, true, QedPetersburg(), "yyy", LSjjSettings(false),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  ) ==#

grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 1500)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.0), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), 
                                        Configuration("1s^2 7s"), Configuration("1s^2 7p"), Configuration("1s^2 7d"), Configuration("1s^2 7f")],
                        intermediateConfigs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p 7s"), Configuration("1s^2 2p 7p"),  
                                             Configuration("1s^2 2p 7d"), Configuration("1s^2 2p 7f")],  
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2s 7s"), Configuration("1s^2 2s 7p"),  
                                        Configuration("1s^2 2s 7d"), Configuration("1s^2 2s 7f")],
                        process = JAC.Dierec, 
                        processSettings=Dielectronic.Settings([E1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                              false, Tuple{Int64,Int64,Int64}[(1,10,1), (1,10,7)], 0., 0., 0., "Coulomb")  )

wb = perform(wa)

@show wa
wb = perform(wa)

end
setDefaults("print summary: close", "")


