#
println("Af) Test of the QED model corrections to the level structure of atoms and ions.")
#
setDefaults("QED model: Petersburg")  ## Petersburg, Sydney
wa = Atomic.Computation("xx",  Nuclear.Model(19.); 
                        configs=[Configuration("[Ar] 4s"), Configuration("[Ar] 4p")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", "xxx", [1],     0, 1.0e-6, JAC.Subshell[], 
                                                true, false, true, "yyy", false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )

wb = perform(wa)



