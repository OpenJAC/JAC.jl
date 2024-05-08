#
println("Ad) Test of the frequency-independent Breit interaction for an internally generated chlorine-like multiplet.")

if  true
    # Last successful:  unknown ... need to be adapted
    # Compute the level structure of Cl-like Fe^10+ ions including the Breit interaction
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(54.), 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],  
                            asfSettings=AsfSettings(true, false, Basics.DFSField(), "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, 
                                                    JAC.Subshell[], JAC.Subshell[], true, false, NoneQed(), LSjjSettings(true),
                                                    false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

    wb = perform(wa)
    #
end



