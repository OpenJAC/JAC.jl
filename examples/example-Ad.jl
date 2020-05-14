#
println("Ad) Tests towards a SCF field: B-spline primitives and one-particle spectra in a local potential.")
#

if  false
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(18., "point"), 
                        properties=JAC.AtomicLevelProperty[],
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],  ## 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", LSjjSettings(true),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

wb = perform(wa)
elseif true
#
# Test for Björn, 7. Mai 2020
## grid=JAC.Radial.Grid(true)
grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 2000)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(8.0, "Fermi"), 
                        properties=JAC.AtomicLevelProperty[],
                        configs=[Configuration("1s^2 2s"), Configuration("1s 2s^2"), Configuration("1s 2s 3p"), Configuration("1s 2s 4p")], 
                        ## configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s 2s 3s"), Configuration("1s 2p 3s"), 
                        ##          Configuration("1s 2s 2p")],  ## 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", LSjjSettings(false),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

wb = perform(wa)

end
