#
println("Jc) Apply & test the Saha-Boltzmann computations for an ionic mixture in LTE.")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

if  true
    # Last successful:  15Jun2024
    # Compute the ...
    nm          = Nuclear.Model(6.0)
    rho         = 2.0e-5      # number density
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 200.)
    settings    = Plasma.Settings(temp_au, rho, true)
    ionMixture  = [IsotopicFraction(6., 12.2, 1.0)]
    scheme      = Plasma.SahaBoltzmannScheme(true, 10, 3, 1, 4, ionMixture, String[])
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)
    
end
