#
println("Ad) Tests towards a SCF field: B-spline primitives and one-particle spectra in a local potential.")
#
wa = Atomic.Computation("xx",  Nuclear.Model(18., "point"); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, NoPoints=390), 
                        properties=JAC.AtomicLevelProperty[],
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],  ## 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", "xxx", [1], 40, 1.0e-6, JAC.Subshell[], 
                                                true, false, false, "yyy", false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )
                        ##x scfSettings=ScfSettings(false, "hydrogenic", " ", [1], false, 24, 1.0e-8, JAC.Subshell[]) )

wb = perform(wa)
error()
