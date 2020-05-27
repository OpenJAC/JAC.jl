#
println("Hb) Tests of the StrongField module to calculate energy and momentum distributions of photoelectron in SFA.")
#
## setDefaults("print summary: open", "zzz-HHG.sum")

if  true
    asfSettings  = AsfSettings(true, false, Basics.DFSField(), "hydrogenic", Dict{Subshell, Orbital}(), [1],    40, 1.0e-6, JAC.Subshell[], JAC.Subshell[], 
                               true, false, NoneQed(), "yyy", LSjjSettings(true), false, [1,2,3,4], false, JAC.LevelSymmetry[] )
    asfSettings  = AsfSettings()                           
    grid         = Radial.Grid(false)
    nuclearModel = Nuclear.Model(3.)
    # Compute the initial and final levels for the SFA computation
    wa           = Atomic.Computation(Atomic.Computation(), name="Li ground-state levels", grid=grid, nuclearModel=nuclearModel, 
                                     configs=[Configuration("[He] 2s")],  asfSettings=asfSettings )
    wb           = perform(wa, output=true)
    #
    observable   = StrongField.SfaEnergyDistribution(0., 0., [0.1, 0.3, 0.5])
    beam         = Pulse.PlaneWaveBeam(1.0, 0.3, 0.)  # (A0, omega, cep)
    envelope     = Pulse.InfiniteEnvelope()
    polarization = Basics.RightCircular()
    sfaSettings  = StrongField.Settings()
    
    wa = StrongField.Computation(observable, nuclearModel, grid, initialLevel, finalLevel, beam, envelope, polarization, sfaSettings)
    wb = StrongField.compute(wa)
end
