#
println("Ia) Tests of average-atom computations.")

if  true
    # Last successful:  unknown ...
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 10.0 * 1)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm) 
    ## gridx    = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 3.0)
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = radiusWS +1.0)
    settings    = Plasma.Settings(temperature, rho)
    qValues     = [ q for q in 1:10 ]
    scheme      = Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField(), false, true, false, Subshell[], Float64[], qValues )
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, 
                                     nuclearModel=nm, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)
    
    ## chemMu      = AverageAtom.determineChemicalPotential(orbitals, temp_au, nm::Nuclear.Model, grid)
end
