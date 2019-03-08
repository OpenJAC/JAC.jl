
println("Ag) Generate and normalize continuum orbitals in a local potential.")
#
wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 8.0e-2, NoPoints = 500),
                        configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", "xxx", [1], 40, 1.0e-6, JAC.Subshell[], 
                                                true, false, false, "no-file", false, [1,2,3,4], false, JAC.LevelSymmetry[] ) )

wb = perform(wa, output=true)
#
multiplet = wb["multiplet:"]
grid      = wb["grid:"]
level     = multiplet.levels[1]
basis     = multiplet.levels[1].basis  
#
println("\nCompute the nuclear and Hartree potentials.")
w1 = JAC.Nuclear.pointNucleus(2.0, grid)
w2 = compute("radial potential: Dirac-Fock-Slater", grid, level)

settings = Continuum.Settings(false, grid.nr-30)
JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalLocalPotential(1.0, Subshell("3s_3/2"), w1, settings)
JAC.define("method: normalization, pure Coulomb")
wc       = JAC.Continuum.generateOrbitalLocalPotential(1.0, Subshell("3s_3/2"), w1, settings)

error("stop here")

settings = Continuum.Settings(false, grid.nr)
JAC.define("method: continuum, pure sine")
JAC.define("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalLocalPotential(100.0, Subshell("3d_3/2"), w1, settings)

JAC.define("method: continuum, asymptotic Coulomb")
JAC.define("method: normalization, pure Coulomb")
wc       = JAC.Continuum.generateOrbitalLocalPotential(100.0, Subshell("3d_3/2"), w1, settings)

settings = Continuum.Settings(false, grid.nr)
JAC.define("method: continuum, nonrelativistic Coulomb")
JAC.define("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalNonrelativisticCoulomb(100.0, Subshell("3s_1/2"), 1.0, w1.grid, settings)

settings = Continuum.Settings(false, grid.nr)
JAC.define("method: continuum, nonrelativistic Coulomb")
JAC.define("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalNonrelativisticCoulomb(100.0, Subshell("3s_1/2"), w1, settings)
