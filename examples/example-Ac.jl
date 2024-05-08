#
println("Ac) Test of the CI part for an internally generated neon multiplet without Breit interaction.")

if  true
    # Last successful:  unknown ... need to be adapted
    # Compute different direct potentials for the charge density of a given level 
    levelSelection = LevelSelection(true, indices = [1,2, 4,5, 7,8]) 
                                    ## symmetries = [ LevelSymmetry(1//2,Basics.plus), LevelSymmetry(1//2,Basics.minus),  LevelSymmetry(5//2,Basics.plus)])
        
    settings1 = AsfSettings()
    settings2 = AsfSettings(AsfSettings(), levelSelectionCI = levelSelection, eeInteractionCI=CoulombInteraction(), jjLS = LSjjSettings(false) )

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(36.), 
                            configs=[Configuration("[Ar] 4f^7")], 
                            ## configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d"), 
                            ##         Configuration("[Ne] 3s^2 3p^3 3d^4"), Configuration("[Ne] 3s^2 3p^2 3d^3"), ],  #, 
                            ##         # Configuration("[Ne] 3s^2 3p^2 3d^3"), Configuration("[Ne] 3s 3p^2 3d^4"), Configuration("[Ne] 3p^2 3d^5")], 
                            asfSettings=settings1 )

    @time wb = perform(wa)
end



