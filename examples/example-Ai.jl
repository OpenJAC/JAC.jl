#
println("Ai) Test of the jj-LS transformation of levels from a given multiplet.")
#
wa = Atomic.Computation("jj-LS level transformation",  Nuclear.Model(26.); 
                        configs=[Configuration("[Ne] 3s^2 3p^4 3d")],  
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", "xxx", [1],     0, 1.0e-6, JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", ManyElectron.LSjjSettings(true),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )

wb = perform(wa)



##                        configs=[Configuration("[Ne] 3s^2 3p^4"), Configuration("[Ne] 3s^2 3d^4")],  , Configuration("1s^2 2s 2p^6")
