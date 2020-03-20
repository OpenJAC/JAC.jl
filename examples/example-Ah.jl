using PyPlot: plot
## pyplot()

println("Ah) Generate and normalize continuum orbitals in a local potential.")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, NoPoints = 1100)
wa   = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(3.), 
                          configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                          asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                  true, false, NoneQed(), "yyy", LSjjSettings(false), false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

wb = perform(wa, output=true)

# Extract the ground level to generate a proper potential
multiplet = wb["multiplet:"]
grid      = wb["grid:"]
level     = multiplet.levels[1]
basis     = multiplet.levels[1].basis  

#
println("\n\n>> Compute the nuclear and Hartree potentials.")
w1  = JAC.Nuclear.pointNucleus(3.0, grid)
w2  = compute("radial potential: Dirac-Fock-Slater", grid, level)
# w3  = compute("radial potential: Hartree",           grid, level)
w3  = Radial.Potential(" ", 0.95*w2.Zr, w2.grid)
ww2 = JAC.add(w1,w2)
ww3 = JAC.add(w1,w3)
## plot(grid.r, w1.Zr)
## plot(grid.r, ww2.Zr)
## plot(grid.r, ww3.Zr)

## error("stop aa")
# Basics.plot("radial potentials", [w1, w2, w3], grid; N=900)
# Basics.plot("radial orbitals: large", [wc[1], wd[1], we[1]], grid; N=900)

## setDefaults("method: continuum, asymptotic Coulomb")
## setDefaults("method: continuum, pure sine")
## setDefaults("method: continuum, spherical Bessel")
setDefaults("method: continuum, Galerkin")
## setDefaults("method: normalization, pure sine")
setDefaults("method: normalization, pure Coulomb")

settings = Continuum.Settings(false, grid.nr)
#
energy   = 1.0;    sh = Subshell("1000p_3/2");      setDefaults("method: normalization, pure sine")
wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
wc2      = Continuum.generateOrbitalLocalPotential(energy, sh, ww2, settings)
#
energy   = 1.0;   sh = Subshell("1000p_3/2");       setDefaults("method: normalization, pure Coulomb")
wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
wc3      = Continuum.generateOrbitalLocalPotential(energy, sh, ww2, settings)
#
## energy   = 1.0;   sh = Subshell("1000p_3/2");       setDefaults("method: continuum, Galerkin")
## wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
## println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
## wc4      = Continuum.generateOrbitalLocalPotential(energy, sh, ww2, settings)
#
plot(grid.r[1:1100], wc2[1].P[1:1100])
plot(grid.r[1:1100], wc3[1].P[1:1100])
## plot(grid.r[1:1100], wc4[1].P[1:1100])

error("stop bb")

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
