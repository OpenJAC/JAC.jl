#
println("Ae) Test of the CI part for an internally generated neon multiplet without Breit interaction.")
#
settings = AsfSettings(AsfSettings(),selectLevelsCI = true, selectedLevelsCI = [1,2, 4,5, 7,8], breitCI=false,
                                     selectSymmetriesCI = true, selectedSymmetriesCI = [ LevelSymmetry(1//2,Basics.plus),  
                                                          LevelSymmetry(1//2,Basics.minus),  LevelSymmetry(5//2,Basics.plus)], jjLS = LSjjSettings(false) )

wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],  #, 
                                 # Configuration("[Ne] 3s^2 3p^2 3d^3"), Configuration("[Ne] 3s 3p^2 3d^4"), Configuration("[Ne] 3p^2 3d^5")], 
                        asfSettings=settings  )

@time wb = perform(wa)



