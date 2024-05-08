
println("Dm)  Test of the InternalRecombination module with ASF from an internally generated initial and final-state multiplet.")

setDefaults("print summary: open", "zzz-InternalRecombination.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


if  true
    # Last successful:  unknown ...
    # Internal recombination of Xe: [Ar] 3d (10l + 11l)  -->  [Ar] 4p (4l + 5l)
    rydbergShells = Basics.generateShellList(6,  6, "p") 
    irSettings    = InternalRecombination.Settings(rydbergShells, true, 5.0, 0.1, LineSelection(), CoulombInteraction())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(54.), 
                            initialConfigs = [Configuration("[Ne] 3d")],
                            finalConfigs   = [Configuration("[Ne] 4p 5s")],  ## , Configuration("[Ar] 4p 5p"), Configuration("[Ar] 4p 5d")],
                            processSettings= irSettings )

    wb = perform(wa)
    #
end
    
