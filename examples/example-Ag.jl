#
println("Ag) Apply & test the jj-LS transformation of levels from a given multiplet.")

if  true
    # Last successful:  unknown 
    # LS-jj transformation for ... 
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS level transformation", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("[Ne] 3s 3p 3d")],  
                            asfSettings=AsfSettings(true, false, Basics.DFSField(), "hydrogenic", Dict{Subshell, Orbital}(), [1],     0, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                    true, false, NoneQed(), LSjjSettings(true),
                                                    false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )
    wb = perform(wa)
    ##                        configs=[Configuration("[Ne] 3s^2 3p^4"), Configuration("[Ne] 3s^2 3d^4")],  , Configuration("1s^2 2s 2p^6")
    ##                      , Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^5")
end
