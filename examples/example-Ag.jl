#
println("Ag) Test of the QED model corrections to the level structure of atoms and ions.")
#
if   true
wa = Atomic.Computation(Atomic.Computation(), name="QED estimates for xenon-like Sn", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(50.), 
                        configs=[Configuration("1s^2 2s^2 2p^6")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, QedPetersburg(), "yyy", LSjjSettings(false),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )
                                                
elseif  false
wa = Atomic.Computation(Atomic.Computation(), name="QED estimates for hydrogen-like Sn", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(50.), 
                        configs=[Configuration("1s"), Configuration("2s"), Configuration("2p")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    0, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, QedPetersburg(), "yyy", LSjjSettings(false),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )
                                                

end

wb = perform(wa)


# NoneQed()
#  
