
"""
`module  JAC.StrongField
    ... a submodel of JAC that contains all methods to set-up and perform strong-field computations. 
"""
module StrongField

    using    Printf, ..Basics, ..Defaults, ..Radial, ..ManyElectron, ..Nuclear, ..Pulse, ..TableStrings

    export   aaaa

    
    """
    `abstract type StrongField.AbstractSFAObservable` 
        ... defines an abstract and a number of types for the observables that can be computed with SFA amplitudes.

        + struct SfaEnergyDistribution      ... to compute the energy distribution of the photoelectrons.
        + struct SfaMomentumDistribution    ... to compute the momentum distribution of the photoelectrons (not yet).
    """
    abstract type  AbstractSFAObservable    end


    """
    `struct  StrongField.SfaEnergyDistribution  <: StrongField.AbstractSFAObservable`   
        ... to compute in SFA the energy distribution of photoelectrons at given energies and angles.

        + theta         ::Float64             ... polar angles of the energy distribution
        + phi           ::Float64
        + energies      ::Array{Float64,1}    ... specifies the photon energy in the current units. 
    """
    struct   SfaEnergyDistribution  <: StrongField.AbstractSFAObservable
        theta           ::Float64
        phi             ::Float64
        energies        ::Array{Float64,1}
    end


    function Base.string(obs::SfaEnergyDistribution)
        sa = "Compute photo-electron energy distribution at angles (theta,phi) = ($(obs.theta), $(obs.phi)) for the energies [a.u.]:" 
        return( sa )
    end

    function Base.show(io::IO, obs::SfaEnergyDistribution)
        sa = string(obs);       print(io, sa, "\n");     print(io, "   ", obs.energies)
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
    struct         CoulombVolkov     <:  AbstractVolkovState   end
    struct         DistortedVolkov   <:  AbstractVolkovState   end

    function Base.string(volkov::FreeVolkov)         return( "apply free-Volkov states" )         end
    function Base.string(volkov::CoulombVolkov)      return( "apply Coulomb-Volkov states" )      end
    function Base.string(volkov::DistortedVolkov)    return( "apply distorted-Volkov states" )    end

    
    """
    `struct  StrongField.Settings`  
        ... defines a type for the details and parameters of computing SFA amplitudes.

        + multipoles            ::Array{EmMultipoles}  ... Multipoles of the radiation field that are to be included.
        + gauges                ::Array{UseGauge}      ... Specifies the gauges to be included into the computations.
        + printBefore           ::Bool                 ... True, if all the requested amplitudes are printed before the computation start.
        + printAmplitudes       ::Bool                 ... True, if the amplitudes are to be printed.
    """
    struct Settings 
        multipoles              ::Array{EmMultipole,1}
        gauges                  ::Array{UseGauge}
        printBefore             ::Bool
        printAmplitudes         ::Bool
    end 


    """
    `StrongField.Settings()`  
        ... constructor for the default values of SFA computations.
    """
    function Settings()
        Settings([E1], UseGauge[], false, false)
    end


    # `Base.show(io::IO, settings::StrongField.Settings)`  ... prepares a proper printout of the variable settings::StrongField.Settings.
    function Base.show(io::IO, settings::StrongField.Settings) 
        println(io, "multipoles:                 $(settings.multipoles)  ")
        println(io, "use-gauges:                 $(settings.gauges)  ")
        println(io, "printBefore:                $(settings.printBefore)  ")
        println(io, "printAmplitudes:            $(settings.printAmplitudes)  ")
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


    # `Base.show(io::IO, comp::StrongField.Computation)`  ... prepares a proper printout of the variable comp::StrongField.Computation.
    function Base.show(io::IO, comp::StrongField.Computation) 
        println(io, "observable:         $(comp.observable)  ")
        println(io, "nuclearModel:       $(comp.nuclearModel)  ")
        println(io, "grid:               $(comp.grid)  ")
        println(io, "initialLevel:       $(comp.initialLevel)  ")
        println(io, "finalLevel:         $(comp.finalLevel)  ")
        println(io, "beam:               $(comp.beam)  ")
        println(io, "envelope:           $(comp.envelope)  ")
        println(io, "polarization:       $(comp.polarization)  ")
        println(io, "volkov:             $(comp.volkov)  ")
        println(io, "settings:           $(comp.settings)  ")
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
        sa = string(op);       print(io, sa, "\n")
    end


    """
    `struct  StrongField.OutcomeEnergyDistribution`   
        ... to comprise the energy distribution of photoelectrons at given energies and angles, for instance,
            for later graphical representation.

        + theta         ::Float64             ... polar angles of the energy distribution
        + phi           ::Float64
        + probabilities ::Array{Float64,1}    ... calculate probabilities of the energy distribution.
    """
    struct   OutcomeEnergyDistribution
        theta           ::Float64
        phi             ::Float64
        probabilities   ::Array{Float64,1}
    end


    # `Base.show(io::IO, outcome::StrongField.OutcomeEnergyDistribution)`  ... prepares a proper printout of the variable outcome::StrongField.OutcomeEnergyDistribution.
    function Base.show(io::IO, outcome::StrongField.OutcomeEnergyDistribution) 
        println(io, "theta:             $(outcome.theta)  ")
        println(io, "phi:               $(outcome.phi)  ")
        println(io, "probabilities:     $(outcome.probabilities)  ")
    end
    
    
    
    ##################################################################################################################################################
    ##################################################################################################################################################
    ##################################################################################################################################################


    """
    `StrongField.computeEnvelopeIntegrals(observable::AbstractSFAObservable, beam::AbstractBeam, envelope::Pulse.AbstractEnvelope, 
                                          polarization::Basics.AbstractPolarization)`  
        ... computes the pulse-envelope integrals for the given laser pulse; this pulse is completely specified by its beam character
            and parameters, the pulse-envelope parameters and the polarization of the light. A tuple of the two integrals
            (fVolkov, fVolkovSquared)::Tuple{Complex{Float64},Complex{Float64}} is returned.
    """
    function computeEnvelopeIntegrals(observable::AbstractSFAObservable, beam::AbstractBeam, envelope::Pulse.AbstractEnvelope, 
                                      polarization::Basics.AbstractPolarization)
        fVolkov        = 1.0
        fVolkovSquared = 1.0 
        return( (fVolkov, fVolkovSquared) )
    end


    """
    `StrongField.computeOutcome(observable::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,1})`  
        ... computes the requested photoelectron energy distribution and returns the results in the variable
            outcome::StrongField.OutcomeEnergyDistribution
    """
    function  computeOutcome(obs::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,1})
        probabilities = Float64[]
        outcome = OutcomeEnergyDistribution(obs.theta, obs.phi, obs.energies, probabilities)
        return( outcome )
    end


    """
    `StrongField.computeSphericalAmplitudes(comp::StrongField.Computation, envIntegrals::Tuple{Complex{Float64},Complex{Float64}})`  
        ... computes all necessary spherical amplitudes for the given initial and final levels, the Volkov states
            and polarization and by using the pulse envelope integrals. A newAmplitudes::Array{SphericalAmplitude,1} is returned.
    """
    function  computeSphericalAmplitudes(comp::StrongField.Computation, envIntegrals::Tuple{Complex{Float64},Complex{Float64}})
        # First determine which spherical SFA amplitudes need to be computed before the actual computation starts
        sfaAmplitudes = StrongField.determineSphericalAmplitudes(comp.observable)
        newAmplitudes = SphericalAmplitude[]
        return( newAmplitudes )
    end


    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaEnergyDistribution)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            An list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaEnergyDistribution)
        amplitudes = StrongField.SphericalAmplitude[]
        for  energy in observable.energies      push!(amplitudes, SphericalAmplitude(energy, theta, phi, 0.))     end
        
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
            envIntegrals  = StrongField.computeEnvelopeIntegrals(comp.observable, comp.beam, comp.envelope, comp.polarization)
            sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp, envIntegrals)
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "energy disribution" => sfaOutcome) )  end
        elseif   typeof(comp.observable) == StrongField.SfaMomentumDistribution
                error("not yet implemented.")
        else     error("Undefined observable for strong-field computations.")
        end
        
        return( 0. )
    end

end # module

