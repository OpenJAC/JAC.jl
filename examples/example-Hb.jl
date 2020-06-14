#
println("Hb) Tests of the StrongField module to calculate energy and momentum distributions of photoelectron in SFA.")
#
## setDefaults("print summary: open", "zzz-HHG.sum")

if  false
    asfSettings  = AsfSettings(AsfSettings(), generateScf=true)                           
    grid         = Radial.Grid(Radial.Grid(false), hp=5.0e-2)
    nuclearModel = Nuclear.Model(3.)
    # Compute the initial and final levels for the SFA computation
    wa           = Atomic.Computation(Atomic.Computation(), name="Li ground-state levels", grid=grid, nuclearModel=nuclearModel, 
                                      configs=[Configuration("[He] 2s")],  asfSettings=asfSettings )
    wb           = perform(wa, output=true)
    initialLevel = wb["multiplet:"].levels[1]
    #
    wa           = Atomic.Computation(Atomic.Computation(), name="He ground-state level", grid=grid, nuclearModel=nuclearModel, 
                                      configs=[Configuration("[He]")],  asfSettings=asfSettings )
    wb           = perform(wa, output=true)
    finalLevel   = wb["multiplet:"].levels[1]
    #
    omega        = convertUnits("energy: from wavelength [nm] to atomic", 800.)
    intensity    = convertUnits("intensity: from W/cm^2 to atomic", 1.0e14)
    A0           = Pulse.computeFieldAmplitude(intensity, omega)
    observable   = StrongField.SfaEnergyDistribution(0., pi/2., 10, 2*omega)  
    beam         = Pulse.PlaneWaveBeam(A0, omega, 0.)  # (A0, omega, cep)
    envelope     = Pulse.InfiniteEnvelope()
    polarization = Basics.RightCircular()
    volkov       = StrongField.FreeVolkov()  ## CoulombVolkov(), DistortedVolkov()
    sfaSettings  = StrongField.Settings()
    
elseif true
    
    wa = StrongField.Computation(observable::StrongField.AbstractSFAObservable, nuclearModel, grid, initialLevel, finalLevel, 
                                 beam, envelope, polarization, volkov, sfaSettings)
    wb = StrongField.perform(wa, output=true)
    
    distribution = wb["energy disribution"]
    @show   distribution
    nothing
end
