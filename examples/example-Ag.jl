#
println("Ag) Test of the QED model corrections to the level structure of atoms and ions.")
#
wa = Atomic.Computation(Atomic.Computation(), name="QED estimates for carbon-like Xe", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(54.), 
                        configs=[Configuration("1s^2 2s^2 2p^6")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", LSjjSettings(true),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

wb = perform(wa)



