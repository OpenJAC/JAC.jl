
"""
`module  JAC.StrongField
    ... a submodel of JAC that contains all methods to set-up and perform strong-field computations. 
"""
module StrongField

    using    GSL, HypergeometricFunctions, Plots, Printf, SpecialFunctions, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..Radial, ..ManyElectron, 
             ..Nuclear, ..Pulse, ..TableStrings

    export   aaaa
    
    """
    `abstract type StrongField.AbstractSFAObservable` 
        ... defines an abstract and a number of types for the observables that can be computed with SFA amplitudes.

        + struct SfaNoObservable            ... an empty instance of an observable that does not help compute anything.
        + struct SfaEnergyDistribution      ... to compute the energy distribution of the photoelectrons.
        + struct SfaMomentumDistribution    ... to compute the momentum distribution of the photoelectrons.
        + struct SfaAzimuthalAngularDistribution     ... to compute the angular distribution (azimuthal angle phi) of photoelectrons with fixed energy and fixed polar angle.
        + struct SfaPolarAngularDistribution     ... to compute the angular distribution (polar angle theta) of photoelectrons with fixed energy and fixed azimuthal angle.
    """
    abstract type  AbstractSFAObservable                                 end
    struct   SfaNoObservable     <: StrongField.AbstractSFAObservable    end


    """
    `struct  StrongField.SfaEnergyDistribution  <: StrongField.AbstractSFAObservable`   
        ... to compute in SFA the energy distribution of photoelectrons at given energies and angles.

        + theta         ::Float64             ... polar angle of the energy distribution
        + phi           ::Float64             ... azimuthal angle of the energy distribution
        + energies      ::Array{Float64,1}    ... specifies the photoelectron energy in the current units. 
    """
    struct   SfaEnergyDistribution  <: StrongField.AbstractSFAObservable
        theta           ::Float64
        phi             ::Float64
        energies        ::Array{Float64,1}
    end

    """
    `StrongField.SfaEnergyDistribution(theta::Float64, phi::Float64, NoEnergies::Int64, maxEnergy::Float64)`  
        ... defines an energy distribution for given (theta,phi) and for NoEnergies between 0 < energy <= maxEnergy;
            a dist::SfaEnergyDistribution is returned.
    """
    function SfaEnergyDistribution(theta::Float64, phi::Float64, NoEnergies::Int64, maxEnergy::Float64)
        energies = Float64[];      for  i=1:NoEnergies   push!(energies, i*maxEnergy / NoEnergies)   end
        SfaEnergyDistribution(theta, phi, energies)
    end

    function Base.string(obs::SfaEnergyDistribution)
        sa = "Compute photo-electron energy distribution at angles (theta,phi) = ($(obs.theta), $(obs.phi)) for the energies [a.u.]:" 
        return( sa )
    end

    function Base.show(io::IO, obs::SfaEnergyDistribution)
        sa = string(obs);       print(io, sa, "\n");     print(io, "   ", obs.energies)
    end
    
    """
    `struct  StrongField.SfaMomentumDistribution  <: StrongField.AbstractSFAObservable`   
        ... to compute in SFA the momentum distribution of photoelectrons for given polar angle theta and at given energies and azimuthal angles.

        + theta         ::Float64             ... polar angle of the energy distribution
        + phi           ::Array{Float64,1}    ... specifies the azimuthal angles
        + energies      ::Array{Float64,1}    ... specifies the photoelectron energy in the current units. 
    """
    struct   SfaMomentumDistribution  <: StrongField.AbstractSFAObservable
        theta           ::Float64
        phis            ::Array{Float64,1}
        energies        ::Array{Float64,1}
    end
    
    """
    `StrongField.SfaMomentumDistribution(theta::Float64, NoPhi::Int64, NoEnergies::Int64, maxEnergy::Float64)`  
        ... defines a momentum distribution for given theta and for NoEnergies between 0 < energy <= maxEnergy and NoPhi between 0 <= phi < 2pi;
            a dist::SfaMomentumDistribution is returned.
    """
    function SfaMomentumDistribution(theta::Float64, NoPhi::Int64, NoEnergies::Int64, maxEnergy::Float64)
        energies = Float64[];      for  i=1:NoEnergies   push!(energies, i*maxEnergy / NoEnergies)   end
        phis = Float64[];          for  i=1:NoPhi        push!(phis, i*2pi / NoPhi)                  end
        SfaMomentumDistribution(theta, phis, energies)
    end

    function Base.string(obs::SfaMomentumDistribution)
        sa = "Compute photo-electron momentum distribution at polar angle theta = $(obs.theta) for (energy [a.u.], phi):" 
        return( sa )
    end

    function Base.show(io::IO, obs::SfaMomentumDistribution)
        sa = string(obs);       print(io, sa, "\n");     print(io, "   ", obs.energies);    print(io, "   ", obs.phi)
    end

    
    """
    `struct  StrongField.SfaAzimuthalAngularDistribution  <: StrongField.AbstractSFAObservable`   
        ... to compute in SFA the angular distribution of photoelectrons for given polar angle theta and at given azimuthal angles for a fixed energy

        + theta         ::Float64             ... polar angle of the energy distribution
        + phi           ::Array{Float64,1}    ... specifies the azimuthal angles
        + energy        ::Float64             ... specifies the photoelectron energy in the current units. 
    """
    struct   SfaAzimuthalAngularDistribution  <: StrongField.AbstractSFAObservable
        theta           ::Float64
        phis            ::Array{Float64,1}
        energy          ::Float64
    end
    
    """
    `StrongField.SfaAzimuthalAngularDistribution(theta::Float64, NoPhi::Int64, energy::Float64)`  
        ... defines a angular distribution for given theta and for NoPhi between 0 <= phi < 2pi and a fixed energy
            a dist::SfaAzimuthalAngularDistribution is returned.
    """
    function SfaAzimuthalAngularDistribution(theta::Float64, NoPhi::Int64, energy::Float64)
        phis = Float64[];          for  i=1:NoPhi        push!(phis, i*2pi / NoPhi)                  end
        SfaAzimuthalAngularDistribution(theta, phis, energy)
    end

    function Base.string(obs::SfaAzimuthalAngularDistribution)
        sa = "Compute photoelectron angular distribution at polar angle theta = $(obs.theta) and energy $(obs.energy) for phi:" 
        return( sa )
    end

    function Base.show(io::IO, obs::SfaAzimuthalAngularDistribution)
        sa = string(obs);       print(io, sa, "\n");    print(io, "   ", obs.phi)
    end
    
    
    """
    `struct  StrongField.SfaPolarAngularDistribution  <: StrongField.AbstractSFAObservable`   
        ... to compute in SFA the angular distribution of photoelectrons for given azimuthal angle phi and at given polar angles for a fixed energy

        + phi         ::Float64             ... phi angle of the energy distribution
        + theta       ::Array{Float64,1}    ... specifies the polar angles
        + energy      ::Float64             ... specifies the photoelectron energy in the current units. 
    """
    struct   SfaPolarAngularDistribution  <: StrongField.AbstractSFAObservable
        phi           ::Float64
        thetas        ::Array{Float64,1}
        energy        ::Float64
    end
    
    """
    `StrongField.SfaPolarAngularDistribution(phi::Float64, NoTheta::Int64, energy::Float64)`  
        ... defines a angular distribution for given phi and for NoPhi between 0 <= phi < 2pi and a fixed energy
            a dist::SfaPolarAngularDistribution is returned.
    """
    function SfaPolarAngularDistribution(phi::Float64, NoTheta::Int64, energy::Float64)
        thetas = Float64[];          for  i=1:NoTheta        push!(thetas, i*2pi / NoTheta)                  end
        SfaPolarAngularDistribution(phi, thetas, energy)
    end

    function Base.string(obs::SfaPolarAngularDistribution)
        sa = "Compute photoelectron angular distribution at azimuthal angle phi = $(obs.phi) and energy $(obs.energy) for theta:" 
        return( sa )
    end

    function Base.show(io::IO, obs::SfaPolarAngularDistribution)
        sa = string(obs);       print(io, sa, "\n");    print(io, "   ", obs.theta)
    end
    
    
    """
    `abstract type StrongField.AbstractVolkovState` 
        ... defines an abstract and a number of types for dealing with the Volkov states in the computation of the SFA amplitudes.

        + struct FreeVolkov       ... to apply the free-Volkov states.
        + struct CoulombVolkov    ... apply the Coulomb-Volkov states.
        + struct DistortedVolkov  ... apply the Distorted-Volkov states.
    """
    abstract type  AbstractVolkovState                         end
    struct         FreeVolkov        <:  AbstractVolkovState   end
    
    
    """
    `struct StrongField.CoulombVolkov     <:  StrongField.AbstractVolkovState`
    
        + Z ::Float64   ... Charge that generates the Coulomb potential
    """
    struct  CoulombVolkov   <: StrongField.AbstractVolkovState
        Z   ::Float64
    end
    
    struct         DistortedVolkov   <:  AbstractVolkovState   end
    
    
    
    function Base.string(volkov::FreeVolkov)         return( "Free-Volkov states" )         end
    function Base.string(volkov::CoulombVolkov)      return( "Coulomb-Volkov states in potential of charge Z = $(volkov.Z)" )      end
    function Base.string(volkov::DistortedVolkov)    return( "Distorted-Volkov states" )    end
    
    
    """
    `struct  StrongField.Settings`  
        ... defines a type for the details and parameters of computing SFA amplitudes.

        + multipoles            ::Array{EmMultipoles}  ... Multipoles of the radiation field that are to be included.
        + gauges                ::Array{UseGauge}      ... Specifies the gauges to be included into the computations.
        + printAmplitudes       ::Bool                 ... True: the amplitudes are printed during computation
        + coupledBasis          ::Bool                 ... True: compute in |j l mj> angular momentum basis; False: compute in |l ml> angular momentum basis
        + hydrogenic            ::Bool                 ... True: Use hydrogenic wave functions for the initial state
        + hydrogenic1s          ::Bool                 ... True: Use hydrogenic wave functions for the initial state
        + mAverage              ::Bool                 ... True: Average over initialState projections mj or ml and sum over continuum spin projection msp
    """
    struct Settings 
        multipoles              ::Array{EmMultipole,1}
        gauges                  ::Array{UseGauge}
        printAmplitudes         ::Bool
        coupledBasis            ::Bool
        hydrogenic              ::Bool
        hydrogenic1s            ::Bool
        mAverage                ::Bool
    end 


    """
    `StrongField.Settings()`  
        ... constructor for the default values of SFA computations.
    """
    function Settings()
        Settings([E1], [UseCoulomb], true, true, false, false, true)
    end


    # `Base.show(io::IO, settings::StrongField.Settings)`  ... prepares a proper printout of the variable settings::StrongField.Settings.
    function Base.show(io::IO, settings::StrongField.Settings) 
        println(io, "multipoles:                 $(settings.multipoles)  ")
        println(io, "use-gauges:                 $(settings.gauges)  ")
        println(io, "printBefore:                $(settings.printBefore)  ")
        println(io, "printAmplitudes:            $(settings.printAmplitudes)  ")
        println(io, "coupledBasis:               $(settings.coupledBasis)  ")
        println(io, "hydrogenic:                 $(settings.hydrogenic)  ")
        println(io, "mAverage:                   $(settings.mAverage)  ")
    end


    """
    `struct  StrongField.Computation`  
        ... defines a type for the computation of strong-field amplitudes and observables (properties).
        
        + observable         ::AbstractSFAObservable        ... The obserable to be calculated in SFA.
        + nuclearModel       ::Nuclear.Model                ... Model, charge and parameters of the nucleus.
        + grid               ::Radial.Grid                  ... The radial grid to be used for the computation.
        + initialLevel       ::Level                        ... Initial level of the atom
        + finalLevel         ::Level                        ... Final level of the atom
        + beam               ::Pulse.AbstractBeam           ... Beam character and properties of the incident light pulse.
        + envelope           ::Pulse.AbstractEnvelope       ... Envelope of the incident light pulse.   
        + polarization       ::Basics.AbstractPolarization  ... Envelope of the incident light pulse. 
        + volkov             ::AbstractVolkovState          ... Specify the treatment/approach of the Volkov states.
        + settings           ::StrongField.Settings         ... Settings for the controlling the SFA amplitudes.
    """
    struct  Computation
        observable           ::AbstractSFAObservable
        nuclearModel         ::Nuclear.Model
        grid                 ::Radial.Grid 
        initialLevel         ::Level
        finalLevel           ::Level
        beam                 ::Pulse.AbstractBeam 
        envelope             ::Pulse.AbstractEnvelope
        polarization         ::Basics.AbstractPolarization 
        volkov               ::AbstractVolkovState
        settings             ::StrongField.Settings
    end 

    
    """
    `StrongField.Computation()`  ... constructor for an `empty` instance of StrongField.Computation().
    """
    function Computation()
        Computation( SfaNoObservable(), Nuclear.Model(1.0), Radial.Grid(true), Orbital("1s_1/2",0.5), Orbital("1s_1/2",0.5), 
                     Pulse.PlaneWaveBeam(), Pulse.InfiniteEnvelope(), Basics.RightCircular(), FreeVolkov(), Settings()  )
    end


    # `Base.string(comp::StrongField.Computation)`  ... provides a String notation for the variable comp::StrongField.Computation.
    function Base.string(comp::StrongField.Computation)
        sa = "Strong-field computation:  " * string(comp.observable) * " for Z = $(comp.nuclearModel.Z)  with " * string(comp.volkov) * "\n"
        sa = sa * " initial level with (J, P, energy) = ($(comp.initialLevel.J), $(comp.initialLevel.parity), $(comp.initialLevel.energy)) \n" 
        sa = sa * " final level with   (J, P, energy) = ($(comp.finalLevel.J), $(comp.finalLevel.parity), $(comp.finalLevel.energy))   \n\n" 
        sa = sa * "The incident laser pulse is described by the: \n "
        sa = sa * string(comp.beam) * "\n " * string(comp.polarization) * "\n " * string(comp.envelope)
        return( sa )
    end


    # `Base.show(io::IO, comp::StrongField.Computation)`  ... prepares a proper printout of the variable comp::StrongField.Computation.
    function Base.show(io::IO, comp::StrongField.Computation) 
        sa = Base.string(comp);             print(io, sa, "\n\n")
        println(io, "Settings:              \n$(comp.settings)    ")  
        println(io, "nuclearModel:          $(comp.nuclearModel)  ")
        println(io, "grid:                  $(comp.grid)  ")
    end


    """
    `struct  StrongField.SphericalAmplitude`   
        ... to keep the amplitude at a given  energy-angular point (energy, theta, phi) in momentum space.

        + energy        ::Float64      ... Kinetic energy of the (outgoing) photoelectron.
        + theta         ::Float64      ... Polar angle of the (outgoing) photoelectron.
        + phi           ::Float64      ... Azimuthal angle of the (outgoing) photoelectron.
        + value         ::ComplexF64   ... (Total) Amplitude.
    """
    struct   SphericalAmplitude
        energy          ::Float64
        theta           ::Float64
        phi             ::Float64 
        value           ::ComplexF64
    end


    function Base.string(amp::StrongField.SphericalAmplitude)
        sa = "Total SFA amplitude at (energy, theta, phi) = ($(amp.energy),$(amp.theta),$(amp.phi)) is $(amp.value)." 
        return( sa )
    end

    function Base.show(io::IO, amp::StrongField.SphericalAmplitude)
        sa = string(amp);       print(io, sa)
    end


    """
    `struct  StrongField.OutcomeEnergyDistribution`   
        ... to comprise the energy distribution of photoelectrons at given energies and angles, for instance,
            for later graphical representation.

        + theta         ::Float64             ... polar angle of the energy distribution
        + phi           ::Float64             ... azimuthal angle of the energy distribution
        + energies      ::Array{Float64,1}    ... selected energies of the distribution.
        + probabilities ::Array{Float64,1}    ... calculated probabilities of the energy distribution.
    """
    struct   OutcomeEnergyDistribution
        theta           ::Float64
        phi             ::Float64
        energies        ::Array{Float64,1}
        probabilities   ::Array{Float64,1}
    end


    # `Base.show(io::IO, outcome::StrongField.OutcomeEnergyDistribution)`  ... prepares a proper printout of the variable outcome::StrongField.OutcomeEnergyDistribution.
    function Base.show(io::IO, outcome::StrongField.OutcomeEnergyDistribution) 
        println(io, "theta:             $(outcome.theta)  ")
        println(io, "phi:               $(outcome.phi)  ")
        println(io, "energies:          $(outcome.energies)  ")
        println(io, "probabilities:     $(outcome.probabilities)  ")
    end
    
    
    """
    `struct  StrongField.OutcomeMomentumDistribution`   
        ... to comprise the momentum distribution of photoelectrons at given energies and angles, for instance,
            for later graphical representation.

        + theta         ::Float64             ... polar angle of the momentum distribution
        + momenta       ::Array{Float64}      ... pairs [px py] of momenta of the momentum distribution
        + probabilities ::Array{Float64,1}    ... calculated probabilities of the momentum distribution.
    """
    struct   OutcomeMomentumDistribution
        theta           ::Float64
        momenta         ::Array{Float64}
        probabilities   ::Array{Float64,1}
    end


    # `Base.show(io::IO, outcome::StrongField.OutcomeMomentumDistribution)`  ... prepares a proper printout of the variable outcome::StrongField.OutcomeMomentumDistribution.
    function Base.show(io::IO, outcome::StrongField.OutcomeMomentumDistribution) 
        println(io, "theta:             $(outcome.theta)  ")
        println(io, "momenta:           $(outcome.momenta)  ")
        println(io, "probabilities:     $(outcome.probabilities)  ")
    end
    
    
    """
    `struct  StrongField.OutcomeAzimuthalAngularDistribution`   
        ... to comprise the angular distribution of photoelectrons at given energy and angles,
            for later graphical representation.

        + theta         ::Float64             ... polar angle of the angular distribution
        + energy        ::Float64             ... photoelectron energy 
        + phis          ::Array{Float64,1}      ... azimuthal angles phi
        + probabilities ::Array{Float64,1}    ... calculated probabilities of the angular distribution.
    """
    struct   OutcomeAzimuthalAngularDistribution
        theta           ::Float64
        energy          ::Float64
        phis            ::Array{Float64,1}
        probabilities   ::Array{Float64,1}
    end


    # `Base.show(io::IO, outcome::StrongField.OutcomeAzimuthalAngularDistribution)`  ... prepares a proper printout of the variable outcome::StrongField.OutcomePolarAngularDistribution.
    function Base.show(io::IO, outcome::StrongField.OutcomeAzimuthalAngularDistribution) 
        println(io, "theta:             $(outcome.theta)  ")
        println(io, "energy:            $(outcome.energy)  ")
        println(io, "phis:              $(outcome.phis)  ")
        println(io, "probabilities:     $(outcome.probabilities)  ")
    end
    
    """
    `struct  StrongField.OutcomePolarAngularDistribution`   
        ... to comprise the angular distribution of photoelectrons at given energy and angles,
            for later graphical representation.

        + phi           ::Float64             ... azimuthal angle of the angular distribution
        + energy        ::Float64             ... photoelectron energy 
        + thetas        ::Array{Float64,1}      ... polar angles phi
        + probabilities ::Array{Float64,1}    ... calculated probabilities of the angular distribution.
    """
    struct   OutcomePolarAngularDistribution
        phi             ::Float64
        energy          ::Float64
        thetas          ::Array{Float64,1}
        probabilities   ::Array{Float64,1}
    end


    # `Base.show(io::IO, outcome::StrongField.OutcomePolarAngularDistribution)`  ... prepares a proper printout of the variable outcome::StrongField.OutcomePolarAngularDistribution.
    function Base.show(io::IO, outcome::StrongField.OutcomePolarAngularDistribution) 
        println(io, "phi:             $(outcome.phi)  ")
        println(io, "energy:            $(outcome.energy)  ")
        println(io, "thetas:              $(outcome.thetas)  ")
        println(io, "probabilities:     $(outcome.probabilities)  ")
    end
    
    ##################################################################################################################################################
    ##################################################################################################################################################
    ##################################################################################################################################################
    
    """
    `StrongField.VolkovP(epsilonp::Float64, lp::Int, rGrid::Array{Float64,1})`  
        ... returns the (plane-wave) Volkov radial wave function
                at the radial grid points rGrid
    """
    function  VolkovP(epsilonp::Float64, lp::Int, rGrid::Array{Float64,1})
        P = Float64[]
        
        for j = 1:length(rGrid)
            r = rGrid[j]
            p = r * GSL.sf_bessel_jl( lp, sqrt(2*epsilonp)*r )
            push!( P, p )
        end
        
        return( P )
    end
    
    """
    `StrongField.CoulombVolkovP(epsilonp::Float64, lp::Int, ZContinuum::Float64, rGrid::Array{Float64,1})`  
        ... returns the Coulomb-Volkov radial wave function
                at the radial grid points rGrid
    """
    function  CoulombVolkovP(epsilonp::Float64, lp::Int, ZContinuum::Float64, rGrid::Array{Float64,1})
        sqrtTwoEps = sqrt(2*epsilonp)
    
        P = ComplexF64[]
        
        for j = 1:length(rGrid)
            r = rGrid[j]
            rho = sqrtTwoEps*r
            eta = -ZContinuum/sqrtTwoEps
            GammaValue = gamma( lp + 1 + eta * im )
            sigma = angle( GammaValue )
            Fl = 2.0^lp * exp(-0.5*pi*eta) * abs( GammaValue ) / factorial(2*lp+1) * rho^(lp+1) * exp( rho * im ) * HypergeometricFunctions.drummond1F1(lp+1+eta*im,2*lp+2,-2*rho*im)
        
        
            p = Fl * exp( sigma * im ) / sqrtTwoEps
            
            push!( P, p )
        end
        
        return( P )
    end
    
    ####################################################################################################################
    ####################################################################################################################

    """
    `StrongField.computeEnvelopeIntegrals(envelope::Pulse.AbstractEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                          thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64)`  
        ... computes the pulse-envelope integrals for the given laser pulse; this pulse is completely specified by its beam character
            and parameters, the pulse-envelope parameters and the polarization of the light. A tuple of the two integrals
            (fVolkov, fVolkovSquared)::Tuple{Complex{Float64},Complex{Float64}} is returned.
    """
    function computeEnvelopeIntegrals(envelope::Pulse.AbstractEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                      thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64)
        fVolkovPlus    = Pulse.envelopeVolkovIntegral(true,  envelope, beam, polarization, thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64)
        fVolkovMinus   = Pulse.envelopeVolkovIntegral(false, envelope, beam, polarization, thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64)
        fVolkovSquared = Pulse.envelopeQuadVolkovIntegral(envelope, beam, polarization,    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64)
        
        #println("> Envelope integrals for a $(string(envelope)) with A0=$(beam.A0), omega=$(beam.omega) [a.u.], cep=$(beam.cep)" * 
        #        " and for $(string(polarization)) light are:" )
        #println("    F^(Volkov) [+; omega; f^(env); A]   = $fVolkovPlus" )
        #println("    F^(Volkov) [-; omega; f^(env); A]   = $fVolkovMinus" )
        #println("    F^(quad, Volkov) [f^(env); A]       = $fVolkovSquared" )
        
        return( (fVolkovPlus, fVolkovMinus, fVolkovSquared) )
    end


    """
    `StrongField.computeOutcome(observable::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,3})`  
        ... computes the requested photoelectron energy distribution and returns the results in the variable
            outcome::StrongField.OutcomeEnergyDistribution
    """
    function  computeOutcome(obs::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,3})
        probabilities = Float64[]
        dimensions = size(amplitudes)
        for  jAmp  = 1:dimensions[3]
            pp = sqrt(2*amplitudes[1,1,jAmp].energy)
            amp = 0.
            for jmi = 1:dimensions[1]
                for jmf = 1:dimensions[2]
                    amp += amplitudes[jmi,jmf,jAmp].value * conj(amplitudes[jmi,jmf,jAmp].value) #sum of probabilities
                end
            end
            push!(probabilities, pp * amp )
        end
        outcome = OutcomeEnergyDistribution(obs.theta, obs.phi, obs.energies, probabilities)
        return( outcome )
    end

    """
    `StrongField.computeOutcome(observable::StrongField.SfaMomentumDistribution, amplitudes::Array{SphericalAmplitude,3})`  
        ... computes the requested photoelectron momentum distribution and returns the results in the variable
            outcome::StrongField.OutcomeMomentumDistribution
    """
    function  computeOutcome(obs::StrongField.SfaMomentumDistribution, amplitudes::Array{SphericalAmplitude,3})
        probabilities   = Float64[]
        dimensions      = size(amplitudes)
        momenta         = Array{Float64}(undef,0,2)
        
        for  jAmp  = 1:dimensions[3]
            pp = sqrt(2*amplitudes[1,1,jAmp].energy)
            px = pp*sin(amplitudes[1,1,jAmp].theta)*cos(amplitudes[1,1,jAmp].phi)
            py = pp*sin(amplitudes[1,1,jAmp].theta)*sin(amplitudes[1,1,jAmp].phi)
            momenta = [momenta ; [px py]]
            amp = 0.
            for jmi = 1:dimensions[1]
                for jmf = 1:dimensions[2]
                    amp += amplitudes[jmi,jmf,jAmp].value * conj(amplitudes[jmi,jmf,jAmp].value) #sum of probabilities
                end
            end
            push!(probabilities, pp * amp )
        end
        outcome = OutcomeMomentumDistribution(obs.theta, momenta, probabilities)
        return( outcome )
    end
    
    """
    `StrongField.computeOutcome(observable::StrongField.SfaAzimuthalAngularDistribution, amplitudes::Array{SphericalAmplitude,3})`  
        ... computes the requested photoelectron angular distribution and returns the results in the variable
            outcome::StrongField.OutcomeAzimuthalAngularDistribution
    """
    function  computeOutcome(obs::StrongField.SfaAzimuthalAngularDistribution, amplitudes::Array{SphericalAmplitude,3})
        probabilities = Float64[]
        dimensions = size(amplitudes)
        for  jAmp  = 1:dimensions[3]
            pp = sqrt(2*amplitudes[1,1,jAmp].energy)
            amp = 0.
            for jmi = 1:dimensions[1]
                for jmf = 1:dimensions[2]
                    amp += amplitudes[jmi,jmf,jAmp].value * conj(amplitudes[jmi,jmf,jAmp].value) #sum of probabilities
                end
            end
            push!(probabilities, pp * amp )
        end
        outcome = OutcomeAzimuthalAngularDistribution(obs.theta, obs.energy, obs.phis, probabilities)
        return( outcome )
    end
    
    """
    `StrongField.computeOutcome(observable::StrongField.SfaPolarAngularDistribution, amplitudes::Array{SphericalAmplitude,3})`  
        ... computes the requested photoelectron angular distribution and returns the results in the variable
            outcome::StrongField.OutcomePolarAngularDistribution
    """
    function  computeOutcome(obs::StrongField.SfaPolarAngularDistribution, amplitudes::Array{SphericalAmplitude,3})
        probabilities = Float64[]
        dimensions = size(amplitudes)
        for  jAmp  = 1:dimensions[3]
            pp = sqrt(2*amplitudes[1,1,jAmp].energy)
            amp = 0.
            for jmi = 1:dimensions[1]
                for jmf = 1:dimensions[2]
                    amp += amplitudes[jmi,jmf,jAmp].value * conj(amplitudes[jmi,jmf,jAmp].value) #sum of probabilities
                end
            end
            push!(probabilities, pp * amp )
        end
        outcome = OutcomePolarAngularDistribution(obs.phi, obs.energy, obs.thetas, probabilities)
        return( outcome )
    end

    
    #---------------------------------------------------------------------------------------------------
    #--------------------------------------FORMULATION in j-l-mj-BASIS:---------------------------------
    #---------------------------------------------------------------------------------------------------
    
    """
    `StrongField.pReducedME(Pepsplp::Array{ComplexF64,1}, lp::Int, jp::Float64, n::Int, l::Int, j::Float64, initialOrbital::Orbital, grid::Radial.Grid)`  
        ... computes the reduced matrix elements of the momentum operator <epsilonp lp jp||p||n l j> in the one-particle picture
            where the radial wave function Pepsplp of |epsilonp lp jp> needs to be provided as an argument on the grid rGrid
    """
    function  pReducedME(Pepsplp::Array{ComplexF64,1}, lp::Int, jp::Float64, l::Int, j::Float64, initialOrbital::Orbital, grid::Radial.Grid)
    
        #---- Compute the factor in <epsilonp lp jp||p||n l j> = fac * <epsilonp lp||p||n l>-----
        fac = 0.
        for    kms = -1:0 #ms = 1/2 + kms
            ms = 1/2 + kms
            for    ml = -l:l
                mj = ml + ms
                for    mlp = -lp:lp
                    mjp = mlp + ms
                    fac = fac + AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(lp)),   AngularM64(Rational(mlp)) , 
                                                                    AngularJ64(1//2),           AngularM64(Rational(ms)), 
                                                                    AngularJ64(Rational(jp)),   AngularM64(Rational(mjp))       )       * 
                                AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(l)),    AngularM64(Rational(ml)),  
                                                                    AngularJ64(1//2),           AngularM64(Rational(ms)), 
                                                                    AngularJ64(Rational(j)),    AngularM64(Rational(mj))       )
                end
            end
        end
        #----------------------------------------------------------------------------------------
        
    
        #--------------------------------Compute <epsilonp lp||p||n l>---------------------------        
        rvalues = grid.r
        rweights = grid.wr
        orderGL = min(size(Pepsplp, 1), size(initialOrbital.P, 1))
        
        integral = 0. * im
        
        #Sum over grid and compute Gauss-Legendre sum
        for    k = 1:orderGL
            r = rvalues[k]
            integrand = conj( Pepsplp[k] )/r * ( r*initialOrbital.Pprime[k] - ((lp-l)*(lp+l+1))/2 * initialOrbital.P[k] )

            #Gauss-Legendre sum
            integral = integral + rweights[k] * integrand
        end

        #Note that GSL.sf_coupling_3j takes takes the input (2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
        integral = integral * (-im)^(lp+1) * (-1)^lp * GSL.sf_coupling_3j( 2*lp, 2*1, 2*l, 0, 0, 0 )
        #----------------------------------------------------------------------------------------
        
        return( fac * integral )
    end
    
    
    """
    `StrongField.scalarProdBoundCont(Pepsplp::Array{ComplexF64,1}, l::Int, initialOrbital::Orbital, grid::Radial.Grid)`  
        ... computes the scalar product of the bound and continuum states in the one-particle picture
            where the radial wave function Pepsplp of |epsilonp lp jp> needs to be provided as an argument on the grid rGrid
    """
    function  scalarProdBoundCont(Pepsplp::Array{ComplexF64,1}, l::Int, initialOrbital::Orbital, grid::Radial.Grid)     
        rvalues = grid.r
        rweights = grid.wr
        
        orderGL = min(size(Pepsplp, 1), size(initialOrbital.P, 1))
        
        integral = 0. * im
        
        #Sum over grid and compute Gauss-Legendre sum
        for    k = 1:orderGL
            r = rvalues[k]
            integrand = conj( Pepsplp[k] ) * initialOrbital.P[k]

            #Gauss-Legendre sum
            integral = integral + rweights[k] * integrand
        end
        
        return((-im)^l * integral)
    end
    
    
    """
    `StrongField.computeSphericalAmplitudes(comp::StrongField.Computation)`  
        ... computes all necessary spherical amplitudes for the given initial and final levels, the Volkov states
            and polarization and by using the pulse envelope integrals - using a coupled angular-momentum basis for the active electron,
            A sfaAmplitudes::Array{SphericalAmplitude,3} is returned.
    """
    function  computeSphericalAmplitudes(comp::StrongField.Computation)
        
        #Extract the initial orbital of the active electron from the many-electron comp.initialLevel and set quantum numbers
        initialOrbitals = comp.initialLevel.basis.orbitals
        
        #Find highest lying orbital (smallest ionization potential)
        defaultSubshell = [sh for (sh,or) in initialOrbitals][1] #This is not nice; must be a better way to simply get a default element from a Dict
        initialOrbital = initialOrbitals[defaultSubshell]
        minIonPotential = abs(initialOrbital.energy)
        for (subshell,orbital) in initialOrbitals
            if   abs(orbital.energy) < minIonPotential
                initialOrbital = orbital
                minIonPotential = abs(orbital.energy)
            end
        end
        
        initialEneV = convertUnits("energy: from atomic to eV", initialOrbital.energy)
        println("")
        println("Selected orbital $(initialOrbital.subshell) as active electron for the strong-field ionization process. Binding energy: $(initialEneV) eV.")
        println("")
        
        ls = LevelSymmetry(initialOrbital.subshell)
        n = initialOrbital.subshell.n;      l = (ls.J.num+1)/2;    j = ls.J.num/2;     initialEn = initialOrbital.energy
        
        if  (sign((-1)^l) == -1 && ls.parity == plus::Parity) || (sign((-1)^l) == 1 && ls.parity == minus::Parity)
            l = l - 1
        end
        l = floor(Int,l)
        
        if comp.settings.hydrogenic && comp.settings.hydrogenic1s
            n = 1
            l = 0
            j = 1/2
        end
        
        lpMin = abs(l-1)
        lpMax = l+1
           
        #Determine which spherical SFA amplitudes need to be computed before the actual computation starts
        mjNum = 1; mspNum = 1    #If comp.Average == false
        mjAverageFactor = 1.     #If comp.Average == false
        if( comp.settings.mAverage )
            mjNum = Int(2*j + 1)
            mspNum = 2
            mjAverageFactor = 1.0/sqrt( mjNum )
        end
        
        sfaAmplitudes = StrongField.determineSphericalAmplitudes(comp.observable,mjNum,mspNum)
        
        #-----Generate continuum wave functions for all energies and values of lp that are needed below-----
        
        println("")
        println("Computing continuum wave functions...")
        
        energies = Array{Float64}(undef,length(sfaAmplitudes[1,1,:]))
        for  jAmp = 1:length(sfaAmplitudes[1,1,:])  #Get all different energy values
            energies[jAmp] = sfaAmplitudes[1,1,jAmp].energy
        end
        energies = unique(energies)
        lpArray = unique( [lpMin l lpMax] )
        
        Pepsplp = Array{ComplexF64}(undef,size(energies)[1],size(lpArray)[1],size(comp.grid.r)[1])
        
        for jlp = 1:size(lpArray)[1]

            lp = lpArray[jlp]
            
            for jE = 1:size(energies)[1]
            
                epsilonp = energies[jE]
        
                if  typeof(comp.volkov) == FreeVolkov
                
                    Pepsplp[jE,jlp,:] = VolkovP( epsilonp, lp, comp.grid.r )
                    
                elseif  typeof(comp.volkov) == CoulombVolkov
                
                    Pepsplp[jE,jlp,:] = CoulombVolkovP( epsilonp, lp, comp.volkov.Z, comp.grid.r )
                    
                elseif  typeof(comp.volkov) == DistortedVolkov
                
                    eta = (-1)^(lp+1)
                    kappa = eta*(lp+1/2)-1/2
                
                    nrContinuum = Continuum.gridConsistency(epsilonp, comp.grid)
                    contSettings = Continuum.Settings(false, comp.grid.NoPoints)

                    newiLevel = Basics.generateLevelWithSymmetryReducedBasis(comp.initialLevel, comp.initialLevel.basis.subshells)
                    newfLevel = Basics.generateLevelWithSymmetryReducedBasis(comp.finalLevel, newiLevel.basis.subshells)
                    newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, kappa), newiLevel)
                    cOrbital, phase  = Continuum.generateOrbitalForLevel(epsilonp, Subshell(101, kappa), newfLevel, comp.nuclearModel, comp.grid, contSettings)
                    
                    jr = size(cOrbital.P)[1]
                    while jr < size(comp.grid.r)[1]
                        push!(cOrbital.P,0)
                        jr += 1
                    end
                    Pepsplp[jE,jlp,:] = cOrbital.P * exp(im*phase)
                    
                end
                
            end #end for energy
            
        end #end for lp
        
        #---------------------------------------------------------------------------------------------------
        
        println("")
        println("Computing SFA amplitudes...")
        println("")
        
        # Compute the requested amplitudes
        for jmj = 1:mjNum
            mj = -j + (jmj-1) 
            for jmsp = 1:mspNum
                msp = -1/2 + (jmsp-1)
                for  jAmp = 1:length(sfaAmplitudes[jmj,jmsp,:]) #amp in  sfaAmplitudes[jmj,jmsp]
                    amp = sfaAmplitudes[jmj,jmsp,jAmp]
                    thetap        = amp.theta;   phip = amp.phi;     
                    energyp = amp.energy; jE = findall(x->x==energyp, vec(energies))[1]
                    envIntegrals  = StrongField.computeEnvelopeIntegrals(comp.envelope, comp.beam, comp.polarization, thetap, phip, energyp, initialEn)
                    fVolkovPlus, fVolkovMinus, fVolkovSquared = envIntegrals
                    #
                    
                    wminus = 0.0im;    wplus = 0.0im
                    # Collect contributions from all l_p, j_p, q terms
                    for  jlp = 1:size(lpArray)[1]
                        lp = lpArray[jlp]
                        
                        kjpMin = floor(Int, ( 2 * abs(lp-1/2) - 1 ) / 2) # jp = (2*kjp + 1)/2, so that the sum below only runs over the correct values. THERE IS A BETTER WAY
                        kjpMax = floor(Int, ( 2 * abs(lp+1/2) - 1 ) / 2)
                        
                        #Compute the sum over total angular momentum jp
                        for kjp = kjpMin:kjpMax 
                            jp = 1/2 * (2 * kjp + 1)
                            
                            KWignerEckart = 1.0/sqrt(2*jp+1) #Factor in the Wigner-Eckart theorem (CONVENTION that needs to fit the reduced ME!)
                            
                            if comp.settings.hydrogenic
                                reducedME = pReducedMEHydrogenic(Pepsplp[jE,jlp,:], lp, jp, initialEn, n, l, j, comp.grid)
                            else
                                reducedME = pReducedME(Pepsplp[jE,jlp,:], lp, jp, l, j, initialOrbital, comp.grid)
                            end
                            
                            
                            #Compute the sum over spherical basis components
                            #Note: sphericalYlm needs an integer in the second argument. mj-msp-q is always a whole number, so we can convert to Int using floor without problems
                            for  q in [-1, 0, 1]
                                wplus = wplus +  KWignerEckart * (-1)^q * Basics.determinePolarizationVector(q, comp.polarization, star=false)   * 
                                                    AngularMomentum.sphericalYlm(lp, floor(Int,mj-msp-q), amp.theta, amp.phi)                           *
                                                    AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(lp)),   AngularM64(Rational(mj-msp-q)), 
                                                                                        AngularJ64(1//2),           AngularM64(Rational(msp)), 
                                                                                        AngularJ64(Rational(jp)),   AngularM64(Rational(mj-q)) )        *
                                                    AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(j)),    AngularM64(Rational(mj)), 
                                                                                        AngularJ64(1),              AngularM64(Rational(-q)), 
                                                                                        AngularJ64(Rational(jp)),   AngularM64(Rational(mj-q)) )         * 
                                                    reducedME
                                        
                                wminus = wminus +  KWignerEckart * Basics.determinePolarizationVector(q, comp.polarization, star=true)           * 
                                                    AngularMomentum.sphericalYlm(lp, floor(Int,mj-msp+q), amp.theta, amp.phi)               *
                                                    AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(lp)),   AngularM64(Rational(mj-msp+q)), 
                                                                                        AngularJ64(1//2),           AngularM64(Rational(msp)), 
                                                                                        AngularJ64(Rational(jp)),   AngularM64(Rational(mj+q)) )        *
                                                    AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(j)),    AngularM64(Rational(mj)), 
                                                                                        AngularJ64(1),              AngularM64(Rational(q)), 
                                                                                        AngularJ64(Rational(jp)),   AngularM64(Rational(mj+q)) )         * 
                                                    reducedME
                            end
                            
                        end # for jp
                        
                    end # for lp
                    
                    #Compute the scalar product in the term proportional to F_2 in the spherical amplitude
                    if comp.settings.hydrogenic
                        jlp = findall(x->x==l, vec(lpArray))[1]
                        scalarProd = scalarProdBoundContHydrogenic(Pepsplp[jE,jlp,:], initialEn, n, l, comp.grid)
                    else
                        jlp = findall(x->x==l, vec(lpArray))[1]
                        scalarProd = scalarProdBoundCont(Pepsplp[jE,jlp,:], l, initialOrbital, comp.grid)
                    end
                    
                    #Compute the total amplitude
                    wa = (fVolkovPlus * wminus + fVolkovMinus * wplus) + fVolkovSquared * AngularMomentum.sphericalYlm(l, floor(Int,mj-msp), amp.theta, amp.phi) *
                                                                            AngularMomentum.ClebschGordan_old(    AngularJ64(Rational(l)),    AngularM64(Rational(mj-msp)), 
                                                                                                                    AngularJ64(1//2),           AngularM64( Rational(msp) ), 
                                                                                                                    AngularJ64(Rational(j)),    AngularM64(Rational(mj))       )  *
                                                                            scalarProd
                    
                    wa = - im/sqrt(2*pi) * wa
                    wa = mjAverageFactor * wa
                    
                    sfaAmplitudes[jmj,jmsp,jAmp] = SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa)
                
                    if( comp.settings.printAmplitudes )
                        println(">> $(SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa))")
                    end
                end #for amp
            end #for jmsp
        end #for jmj
                
        return( sfaAmplitudes )
    end
    
    #---------------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------
    

    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaEnergyDistribution, miNum::Int64, mfNum::Int64)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            A list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaEnergyDistribution, miNum::Int64, mfNum::Int64 )
        ENum = length(observable.energies)
        amplitudes = Array{StrongField.SphericalAmplitude}(undef, miNum, mfNum, ENum)
        
        for jmi = 1:miNum
            for jmf = 1:mfNum
                for jE = 1:ENum 
                    energy = observable.energies[jE]
                    amplitudes[jmi,jmf,jE] = SphericalAmplitude(observable.energies[jE], observable.theta, observable.phi, 0.)    
                end
            end
        end
        
        println("> A total of $(length(amplitudes)) spherical amplitudes need to be calculated.")
        return( amplitudes )
    end
    
    
    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaMomentumDistribution, miNum::Int64, mfNum::Int64)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            A list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaMomentumDistribution, miNum::Int64, mfNum::Int64)
        ENum = length(observable.energies)
        phiNum = length(observable.phis)
    
        amplitudes = Array{StrongField.SphericalAmplitude}(undef, miNum, mfNum, ENum*phiNum)
        
        for jmi = 1:miNum
            for jmf = 1:mfNum
                ctr = 1
                for  jE = 1:ENum
                    for  jphi = 1:phiNum
                        amplitudes[jmi,jmf,ctr] = SphericalAmplitude(observable.energies[jE], observable.theta, observable.phis[jphi], 0.)
                        ctr = ctr + 1
                    end
                end
            end
        end
        
        println("> A total of $(length(amplitudes)) spherical amplitudes need to be calculated.")
        return( amplitudes )
    end
    
    
    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaAzimuthalAngularDistribution, miNum::Int64, mfNum::Int64)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            A list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaAzimuthalAngularDistribution, miNum::Int64, mfNum::Int64)
        phiNum = length(observable.phis)
        amplitudes = Array{StrongField.SphericalAmplitude}(undef, miNum, mfNum, phiNum)
        
        for jmi = 1:miNum
            for jmf = 1:mfNum
                for  jphi = 1:phiNum
                    amplitudes[jmi,jmf,jphi] = SphericalAmplitude(observable.energy, observable.theta, observable.phis[jphi], 0.)
                end
            end
        end
        
        println("> A total of $(length(amplitudes)) spherical amplitudes need to be calculated.")
        return( amplitudes )
    end
    
    
    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaPolarAngularDistribution, miNum::Int64, mfNum::Int64)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            A list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaPolarAngularDistribution, miNum::Int64, mfNum::Int64)
        thetaNum = length(observable.thetas)
        amplitudes = Array{StrongField.SphericalAmplitude}(undef, miNum, mfNum, thetaNum)
        
        for jmi = 1:miNum
            for jmf = 1:mfNum
                for  jtheta = 1:thetaNum
                    amplitudes[jmi,jmf,thetaNum] = SphericalAmplitude(observable.energy, observable.thetas[jtheta], observable.phi, 0.)
                end
            end
        end
        
        println("> A total of $(length(amplitudes)) spherical amplitudes need to be calculated.")
        return( amplitudes )
    end

    

    """
    `StrongField.perform(comp::StrongField.Computation; output::Bool=false)` 
        ... to perform a computation of (one) selected observable that is related to SFA, such as a photoelectron energy 
            distribution, momentum distribution or the computation of sidebands.
    """
    function perform(comp::StrongField.Computation; output::Bool=false)
        if  output    results = Dict{String, Any}()    else    results = nothing    end
        nModel = comp.nuclearModel
        
        if       typeof(comp.observable) == StrongField.SfaEnergyDistribution
            if  comp.settings.coupledBasis
                sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp)
            else
                sfaAmplitudes = StrongField.computeSphericalAmplitudesUncoupled(comp)
            end
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "energy distribution" => sfaOutcome) )  end
        elseif   typeof(comp.observable) == StrongField.SfaMomentumDistribution
            if  comp.settings.coupledBasis
                sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp)
            else
                sfaAmplitudes = StrongField.computeSphericalAmplitudesUncoupled(comp)
            end
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "momentum distribution" => sfaOutcome) )  end
        elseif   typeof(comp.observable) == StrongField.SfaAzimuthalAngularDistribution
            if  comp.settings.coupledBasis
                sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp)
            else
                sfaAmplitudes = StrongField.computeSphericalAmplitudesUncoupled(comp)
            end
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "angular distribution" => sfaOutcome) )  end
        elseif   typeof(comp.observable) == StrongField.SfaPolarAngularDistribution
            if  comp.settings.coupledBasis
                sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp)
            else
                sfaAmplitudes = StrongField.computeSphericalAmplitudesUncoupled(comp)
            end
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "angular distribution" => sfaOutcome) )  end
        else     error("Undefined observable for strong-field computations.")
        end
        
        if  output   return( results )   else    return( nothing )   end
    end
    
    #######################################################################################################################################
    #######################################################################################################################################
    
    include("module-StrongField-inc-hydrogenic.jl") #StrongField routines for hydrogenic initial states
    include("module-StrongField-inc-uncoupled.jl")  #StrongField routines for computations in the basis of "uncoupled" angular momentum l,ml instead of j,l,mj
    ## include("module-StrongField-inc-postProcessing.jl")       #StrongField.plot(...) routine to plot results of StrongField computations

end # module

