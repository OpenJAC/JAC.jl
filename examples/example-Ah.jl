using PyPlot: plot
## pyplot()

println("Ah) Generate and normalize continuum orbitals in a local potential.")

nuclearZ = 5.0;     noPoints = 1100
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, NoPoints = noPoints)

if  false
    # Define a full nuclear + Dirac-Fock-Slater potential
    wa   = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(nuclearZ), 
                              configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                              asfSettings=AsfSettings(true, false, "meanDFS", "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                                                      true, false, NoneQed(), "yyy", LSjjSettings(false), false, [1,2,3,4], false, JAC.LevelSymmetry[] )  )

    wb = perform(wa, output=true)

    # Extract the ground level to generate a proper potential
    multiplet = wb["multiplet:"]
    ## grid      = wb["grid:"]
    level     = multiplet.levels[1]
    basis     = multiplet.levels[1].basis  

    #
    println("\n\n>> Compute the nuclear and Hartree potentials.")
    w1  = JAC.Nuclear.pointNucleus(nuclearZ, grid)
    w2  = compute("radial potential: Dirac-Fock-Slater", grid, level)     ## compute("radial potential: Hartree",grid, level)
    ww = JAC.add(w1,w2)                                                   ## ww = Radial.Potential(" ", 0.95*ww.Zr, w2.grid)
else
    # Define a Coulomb potential
    ww  = JAC.Nuclear.pointNucleus(nuclearZ, grid)                         ## compute("radial potential: Hartree",grid, level)
    plot(grid.r, ww.Zr)
end
## plot(grid.r, ww.Zr);   error("stop aa")


if  true
    # Compute and normalize continuum orbitals for a range of energies
    setDefaults("method: continuum, Galerkin")
    ## setDefaults("method: continuum, asymptotic Coulomb")   "method: continuum, pure sine"   "method: continuum, spherical Bessel"
    settings = Continuum.Settings(false, grid.nr)
    sh       = Subshell("1000s_1/2")
    energies = [x  for  x = 10.0:0.15:13.0];   
    phases1  = Float64[];   phases2 = Float64[];   phases3 = Float64[];   norms1 = Float64[];   norms2 = Float64[];   norms3 = Float64[]
    #
    for  energy in energies   
        wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
        println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
        setDefaults("method: normalization, pure sine")
        wc1      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        push!(phases1, wc1[2]);   push!(norms1, abs(wc1[3]))
        #
        setDefaults("method: normalization, pure Coulomb")
        wc2      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        push!(phases2, wc2[2]);   push!(norms2, abs(wc2[3]))
        #
        setDefaults("method: normalization, Ong-Russek")
        wc3      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        push!(phases3, wc3[2]);   push!(norms3, abs(wc3[3]))
    end
    #
    plot(energies, phases3);   ## plot(energies, phases3);   ## plot(energies, phases3)
    plot(energies, norms1);    plot(energies, norms3) ## ;    plot(energies, norms3)
    
else
    # Compute and display single continuum orbital
    setDefaults("method: continuum, Galerkin")
    ## setDefaults("method: continuum, asymptotic Coulomb")   "method: continuum, pure sine"   "method: continuum, spherical Bessel"
    settings = Continuum.Settings(false, grid.nr)
    sh       = Subshell("1000s_1/2")
    energies = [10.0];   
    for  energy in energies   
        wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
        println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
        setDefaults("method: normalization, pure sine")
        wc1      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        plot(grid.r[1:noPoints], wc1[1].P[1:noPoints])
    end
end
