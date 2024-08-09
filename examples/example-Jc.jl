#
println("Jc) Apply & test the Saha-Boltzmann computations for an ionic mixture in LTE.")

setDefaults("print summary: open", "zzz-saha-boltzmann.sum")


if  false
    # Last successful:  20Jul2024
    # Compute the Saha-Boltzmann equilibrium densities for a mixture of Carbon ions
    grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    nm          = Nuclear.Model(6.0)
    rho         = 2.0e-5      # number density
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 200.)
    settings    = Plasma.Settings(temp_au, rho, true)
    ionMixture  = [IsotopicFraction(6., 12.011, 1.0)]  ## , IsotopicFraction(9., 20.2, 0.4)
    scheme      = Plasma.SahaBoltzmannScheme(true, true, 4:6, 10000, 1, 4, ionMixture, ["IonicLevelDataZ6A12RobinA.jld"])
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)

elseif  true
    # Last successful:  20Jul2024
    # Compute the Saha-Boltzmann equilibrium densities for a mixture of Carbon ions including all charge states,
    # NoExcitations = 2, upperShellNo = 8
    grid        = Radial.Grid(Radial.Grid(true), rnt = 1.0e-4, h = 5.0e-2, hp = 0., rbox = 100.0)
    nm          = Nuclear.Model(6.0)
    rho         = 2.0e-5      # number density
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 200.)
    settings    = Plasma.Settings(temp_au, rho, true)
    ionMixture  = [IsotopicFraction(6., 12.2, 1.0)]
    ## scheme   = Plasma.SahaBoltzmannScheme(Basics.NoPlasmaModel(), true, false, 10000, 7, 2, 3, ionMixture, String["IonicLevelDataZ6A12-allq-n3.jld"])
    scheme      = Plasma.SahaBoltzmannScheme(Basics.NoPlasmaModel(), true, false, 2:6, 10000, 2, 8, ionMixture, String["IonicLevelDataZ6A12-allq-n8.jld"])
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)
    
end
#
setDefaults("print summary: close", "")

# Filenames generated
# "IonicLevelDataZ6A12.jld", "IonicLevelDataZ9A20.jld"
