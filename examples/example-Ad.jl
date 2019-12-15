#
println("Ad) Tests towards a SCF field: B-spline primitives and one-particle spectra in a local potential.")
#
wa = Atomic.Computation("xx",  Nuclear.Model(18., "point"); grid=JAC.Radial.Grid(true), 
                        properties=JAC.AtomicLevelProperty[],
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],  ## 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", LSjjSettings(true),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

wb = perform(wa)
