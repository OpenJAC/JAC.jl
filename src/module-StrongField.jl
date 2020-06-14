
"""
`module  JAC.StrongField
    ... a submodel of JAC that contains all methods to set-up and perform strong-field computations. 
"""
module StrongField

    using    Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..Radial, ..ManyElectron, 
             ..Nuclear, ..Pulse, ..TableStrings

    export   aaaa

    
    """
    `abstract type StrongField.AbstractSFAObservable` 
        ... defines an abstract and a number of types for the observables that can be computed with SFA amplitudes.

        + struct SfaNoObservable            ... an empty instance of an observable that does not help compute anything.
        + struct SfaEnergyDistribution      ... to compute the energy distribution of the photoelectrons.
        + struct SfaMomentumDistribution    ... to compute the momentum distribution of the photoelectrons (not yet).
    """
    abstract type  AbstractSFAObservable                                 end
    struct   SfaNoObservable     <: StrongField.AbstractSFAObservable    end


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

    function Base.string(volkov::FreeVolkov)         return( "free-Volkov states" )         end
    function Base.string(volkov::CoulombVolkov)      return( "Coulomb-Volkov states" )      end
    function Base.string(volkov::DistortedVolkov)    return( "distorted-Volkov states" )    end

    
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
        Settings([E1], [UseCoulomb], false, false)
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

    
    """
    `StrongField.Computation()`  ... constructor for an `empty` instance of StrongField.Computation().
    """
    function Computation()
        Computation( SfaNoObservable(), Nuclear.Model(1.0), Radial.Grid(true), Level(), Level(), 
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

        + theta         ::Float64             ... polar angles of the energy distribution
        + phi           ::Float64
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
    
    
    
    ##################################################################################################################################################
    ##################################################################################################################################################
    ##################################################################################################################################################


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
        
        println("> Envelope integrals for an $(string(envelope)) with A0=$(beam.A0), omega=$(beam.omega) [a.u.], cep=$(beam.cep)" * 
                " and for $(string(polarization)) light are:" )
        println("    F^(Volkov) [+; omega; f^(env); A]   = $fVolkovPlus" )
        println("    F^(Volkov) [-; omega; f^(env); A]   = $fVolkovMinus" )
        println("    F^(quad, Volkov) [f^(env); A]       = $fVolkovSquared" )
        
        return( (fVolkovPlus, fVolkovMinus, fVolkovSquared) )
    end


    """
    `StrongField.computeOutcome(observable::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,1})`  
        ... computes the requested photoelectron energy distribution and returns the results in the variable
            outcome::StrongField.OutcomeEnergyDistribution
    """
    function  computeOutcome(obs::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,1})
        probabilities = Float64[]
        for  amp  in  amplitudes   pp = sqrt(2*amp.energy);  push!(probabilities, pp * amp.value * conj(amp.value) )      end
        outcome = OutcomeEnergyDistribution(obs.theta, obs.phi, obs.energies, probabilities)
        return( outcome )
    end


    """
    `StrongField.computeSphericalAmplitudes(comp::StrongField.Computation)`  
        ... computes all necessary spherical amplitudes for the given initial and final levels, the Volkov states
            and polarization and by using the pulse envelope integrals. A newAmplitudes::Array{SphericalAmplitude,1} is returned.
    """
    function  computeSphericalAmplitudes(comp::StrongField.Computation)
        newAmplitudes = SphericalAmplitude[]
        #
        # First determine which spherical SFA amplitudes need to be computed before the actual computation starts
        sfaAmplitudes = StrongField.determineSphericalAmplitudes(comp.observable)
        # Determine quantum numbers of the initial and final state
        nqn = 1;    lqn = 0;    m = 0;    initialEn = comp.initialLevel.energy
        #
        # Compute the requested amplitudes
        for  amp in  sfaAmplitudes
            thetap        = amp.theta;   phip = amp.phi;     energyp = amp.energy
            envIntegrals  = StrongField.computeEnvelopeIntegrals(comp.envelope, comp.beam, comp.polarization, thetap, phip, energyp, initialEn)
            fVolkovPlus, fVolkovMinus, fVolkovSquared = envIntegrals
            #
            wminus = 0.0im;    wplus = 0.0im
            # Collect contributions from all l_p, q terms; this summation will change in a (nljm) representation
            for  lp = 0:2
                for  q in [-1, 0, 1]
                    reducedME = 1.0  ##  InteractionStrength.nrMomentum(amp.energy, lp, nqn, lqn)
                    wminus = wminus + AngularMomentum.sphericalYlm(lp, m-q, amp.theta, amp.phi) * (-1)^q  * 
                             Basics.determinePolarizationVector(q, comp.polarization)                     *
                             AngularMomentum.ClebschGordan(lp, m, 1, -q, lp, m-q) * reducedME
                    wminus = wminus + AngularMomentum.sphericalYlm(lqn, m+q, amp.theta, amp.phi)          * 
                             Basics.determinePolarizationVector(q, comp.polarization, star=true)          *
                             AngularMomentum.ClebschGordan(lqn, m, 1, q, lp, m+q) * reducedME 
                end
            end
            wa = -im * sqrt(2/pi) * (fVolkovMinus * wminus + fVolkovPlus * wplus) - im / sqrt(2*pi) * fVolkovSquared 
            push!(newAmplitudes, SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa))
            
            println(">> $(SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa))")
        end
        
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
        for  energy in observable.energies      push!(amplitudes, SphericalAmplitude(energy, observable.theta, observable.phi, 0.))     end
        
        println("> A totel of $(length(amplitudes)) spherical amplitudes need to be calculated.")
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
            sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp)
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "energy disribution" => sfaOutcome) )  end
        elseif   typeof(comp.observable) == StrongField.SfaMomentumDistribution
                 error("not yet implemented.")
        else     error("Undefined observable for strong-field computations.")
        end
        
        if  output   return( results)   else    return( nothing )   end
    end

end # module

