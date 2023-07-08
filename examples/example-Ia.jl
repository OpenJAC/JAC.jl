#
println("Ia) Tests of average-atom computations.")
#

if  true
    nm          = Nuclear.Model(18.0)
    rho         = 100.0      # [g/cm^3]
    temperature = 1000000.  # [K]
    temp_au     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperature);   @show temp_au
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm) 
    grid        = Radial.generateGrid(Radial.Grid(true), boxSize=3.0)
    settings    = Plasma.Settings(temperature, rho)
    
    wa          = Plasma.Computation(Plasma.Computation(), 
                                     scheme=Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField() ), nuclearModel=nm, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa)
    
    ## chemMu      = AverageAtom.determineChemicalPotential(orbitals, temp_au, nm::Nuclear.Model)
end
