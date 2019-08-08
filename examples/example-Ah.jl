using Plots
pyplot()

println("Ah) Generate and normalize continuum orbitals in a local potential.")

@warn("\n\n !!! This example does not work properly at present !!! \n\n")
#
wa = Atomic.Computation("xx",  Nuclear.Model(3.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 800),
                        configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                        asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", "xxx", [1],    40, 1.0e-6, JAC.Subshell[], 
                                                true, false, NoneQed(), "yyy", LSjjSettings(true),
                                                false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

wb = perform(wa, output=true)
#
multiplet = wb["multiplet:"]
grid      = wb["grid:"]
level     = multiplet.levels[1]
basis     = multiplet.levels[1].basis  
#
println("\nCompute the nuclear and Hartree potentials.")
w1 = JAC.Nuclear.pointNucleus(5.0, grid)
## w2 = compute("radial potential: Dirac-Fock-Slater", grid, level)
## w2 = JAC.add(w1,w2)
w2 = JAC.Nuclear.nuclearPotentialDH(Nuclear.Model(5.), grid, 0.9)
w3 = JAC.Nuclear.nuclearPotentialDH(Nuclear.Model(5.), grid, 2.5)
plot("radial potentials", [w1, w2, w3], grid; N=500)

settings = Continuum.Settings(false, grid.nr-100)
setDefaults("method: continuum, Galerkin")
## setDefaults("method: normalization, pure sine")
setDefaults("method: normalization, pure Coulomb")
energy   = 30.;   sh = Subshell("3f_5/2")
wc       = JAC.Continuum.generateOrbitalLocalPotential(energy, sh, w1, settings)
wd       = JAC.Continuum.generateOrbitalLocalPotential(energy, sh, w2, settings)
we       = JAC.Continuum.generateOrbitalLocalPotential(energy, sh, w3, settings)
wavenb   = sqrt( 2energy + energy * JAC.give("alpha")^2 )
wavelgth = 2pi / wavenb
println("wavelength = $wavelgth")
plot("radial orbitals: large", [wc[1], wd[1], we[1]], grid; N=900)

error("stop here")

settings = Continuum.Settings(false, grid.nr)
setDefaults("method: continuum, pure sine")
setDefaults("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalLocalPotential(100.0, Subshell("3d_3/2"), w1, settings)

setDefaults("method: continuum, asymptotic Coulomb")
setDefaults("method: normalization, pure Coulomb")
wc       = JAC.Continuum.generateOrbitalLocalPotential(100.0, Subshell("3d_3/2"), w1, settings)

settings = Continuum.Settings(false, grid.nr)
setDefaults("method: continuum, nonrelativistic Coulomb")
setDefaults("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalNonrelativisticCoulomb(100.0, Subshell("3s_1/2"), 1.0, w1.grid, settings)

settings = Continuum.Settings(false, grid.nr)
setDefaults("method: continuum, nonrelativistic Coulomb")
setDefaults("method: normalization, pure sine")
wc       = JAC.Continuum.generateOrbitalNonrelativisticCoulomb(100.0, Subshell("3s_1/2"), w1, settings)
