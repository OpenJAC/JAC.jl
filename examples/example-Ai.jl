#
println("Ai) Test of the jj-LS transformation of levels from a given multiplet.")
#
wa = Atomic.Computation("jj-LS level transformation",  Nuclear.Model(26.); 
                        configs=[Configuration("[Ne] 3p 3d^2")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],     0, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", LSjjSettings(true),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )

wb = perform(wa)



##                        configs=[Configuration("[Ne] 3s^2 3p^4"), Configuration("[Ne] 3s^2 3d^4")],  , Configuration("1s^2 2s 2p^6")
##                      , Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^5")
