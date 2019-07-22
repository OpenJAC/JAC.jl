#
println("Af) Test of the QED model corrections to the level structure of atoms and ions.")
#
wa = Atomic.Computation("QED estimates for carbon-like Xe",  Nuclear.Model(54.); 
                        configs=[Configuration("1s^2 2s^2 2p^6")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", "xxx", [1],     0, 1.0e-6, JAC.Subshell[], 
                                                true, false, QedPetersburg(), "yyy", LSjjSettings(false), false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )

wb = perform(wa)



