#
println("Hb) Tests of the StrongField module to calculate energy and momentum distributions of photoelectron in SFA.")
#
## setDefaults("print summary: open", "zzz-HHG.sum")

#Pulse shape integrals: Properties of pulses

using DelimitedFiles

if  true
    asfSettings  = AsfSettings(AsfSettings(), generateScf=true)                           
    grid         = Radial.Grid(Radial.Grid(false), hp=5.0e-2)
    nuclearModel = Nuclear.Model(3.)
    # Compute the initial and final levels for the SFA computation
    wa           = Atomic.Computation(Atomic.Computation(), name="Li ground-state levels", grid=grid, nuclearModel=nuclearModel, configs=[Configuration("[He] 2s")],  asfSettings=asfSettings )
    wb           = perform(wa, output=true)
    initialLevel = wb["multiplet:"].levels[1]
    #
    wa           = Atomic.Computation(Atomic.Computation(), name="He ground-state level", grid=grid, nuclearModel=nuclearModel, configs=[Configuration("[He]")],  asfSettings=asfSettings )
    wb           = perform(wa, output=true)
    finalLevel   = wb["multiplet:"].levels[1]
    #
    #omega        = convertUnits("energy: from wavelength [nm] to atomic", 800.)
    omega        = 0.056939228883964354555047293615643866360187530517578125 #Since we use different conversion factors (precision!) in Mathematica and JAC
    intensity    = convertUnits("intensity: from W/cm^2 to atomic", 1.0e14)
    #A0           = Pulse.computeFieldAmplitude(intensity, omega)
    A0           = 0.4015398801828808927893987856805324554443359375 #Since we use different conversion factors (precision!) in Mathematica and JAC
    
    observable   = StrongField.SfaEnergyDistribution(pi/2, pi/2, 10, 6*omega) #Arguments are: theta, phi, number of energies, maximum energy
    #observable    = StrongField.SfaMomentumDistribution(pi/2, 360, 5, 10*omega) #Arguments are: theta, number of azimuthal angles, number of energies, masimum energy
    
    beam         = Pulse.PlaneWaveBeam(A0, omega, 0.)  # (A0, omega, cep)
    #envelope     = Pulse.InfiniteEnvelope()
    envelope     = Pulse.SinSquaredEnvelope(14)
    #polarization = Basics.RightCircular()
    polarization = Basics.RightElliptical(0.5)
    volkov       = StrongField.CoulombVolkov(1)  ## FreeVolkov(), CoulombVolkov(Z) with the charge of the ion Z, DistortedVolkov()

    sfaSettings  = StrongField.Settings()
end
#elseif true
 
if true
    wa = StrongField.Computation(observable::StrongField.AbstractSFAObservable, nuclearModel, grid, initialLevel, finalLevel, beam, envelope, polarization, volkov, sfaSettings)
    wb = StrongField.perform(wa, output=true)
    
    distribution = wb["energy distribution"]
    
    writedlm("Probabilities.csv",distribution.probabilities)
    
    nothing
end

if true
    StrongField.plot(wa,wb,"omega","linear")
end

#Convergence test for pulse shape integrals
if false
    Tp = 2;
    envelope = Pulse.GaussianEnvelope(Tp * 110.3489708296542008205005779602025351176903820678791359884044538)
    FPlus = ComplexF64[]
    
    thetap = pi/2
    phip = 0.
    
    Ep = 0.1138784578
    initialEn = initialLevel.energy
    
    orderGH = Int[]
    
    for j = 1:500
        push!(orderGH, j*100)
        push!(FPlus, Pulse.envelopeVolkovIntegral( true, envelope, beam, polarization, thetap, phip, 0.1138784578, initialEn, j*100  ) )
    end
    
    writedlm("ReFPlus.csv",real(FPlus))
    writedlm("ImFPlus.csv",imag(FPlus))
    
end
