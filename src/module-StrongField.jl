
"""
`module  JAC.StrongField
    ... a submodel of JAC that contains all methods to set-up and perform strong-field computations. 
"""
module StrongField

    using    GSL, ## HypergeometricFunctions, Plots, 
             Printf, SpecialFunctions, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..Radial, ..ManyElectron, 
             ..Nuclear, ..Pulse, ..TableStrings

    export   aaaa

    
    """
    `abstract type StrongField.AbstractSFAObservable` 
        ... defines an abstract and a number of types for the observables that can be computed with SFA amplitudes.

        + struct SfaNoObservable            ... an empty instance of an observable that does not help compute anything.
        + struct SfaEnergyDistribution      ... to compute the energy distribution of the photoelectrons.
        + struct SfaMomentumDistribution    ... to compute the momentum distribution of the photoelectrons.
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
    
    
    
    function Base.string(volkov::FreeVolkov)         return( "free-Volkov states" )         end
    function Base.string(volkov::CoulombVolkov)      return( "Coulomb-Volkov states in potential of charge Z = $(volkov.Z)" )      end
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
        + probabilities ::Array{Float64,1}    ... calculated probabilities of the energy distribution.
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
    
    
    ##################################################################################################################################################
    ##################################################################################################################################################
    ##################################################################################################################################################

    
    ####################################################################################################################
    ######The following four functions will be replaced later and serve only for comparison with previous results#######
    ####################################################################################################################
 
    """
    `StrongField.HydrogenPnl(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})`  
        ... returns the non-relativistic hydrogen-like radial wave function with modified binding energy epsiloni
                at the radial grid points rGrid
    """
    function  HydrogenPnl(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})
        Pnl = Float64[]
        Z = n*sqrt(-2*epsiloni)
        
        for j = 1:length(rGrid)
            r = rGrid[j]
            p = r * sqrt( (2*Z/n)^3 * factorial(n-l-1)/(2*n*factorial(n+l)) ) * exp(-Z*r/n) * (2*Z*r/n)^l * GSL.sf_laguerre_n(n-l-1,2*l+1,2*Z*r/n)
            push!( Pnl, p )
        end
        
        return( Pnl )
    end
    
    """
    `StrongField.HydrogenDPnlDr(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})`  
        ... returns the r-derivative of the non-relativistic hydrogen-like radial wave function with modified binding energy epsiloni
                at the radial grid points rGrid
    """
    function  HydrogenDPnlDr(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})
        DPnl = Float64[]
        Z = n*sqrt(-2*epsiloni)
        
        for j = 1:length(rGrid)
            r = rGrid[j]
            if  n-l-2 >= 0
                p = 1/n * 2^(l+1) * exp(-r*Z/n) * (r*Z/n)^l * sqrt( Z^3 * factorial(n-l-1) / (n^4*factorial(l+n)) ) * ( (n+l*n-r*Z) * GSL.sf_laguerre_n(n-l-1,2*l+1,2*r*Z/n) - 2*r*Z * GSL.sf_laguerre_n(n-l-2,2+2*l,2*r*Z/n) )
            elseif  n == 1 && l == 0
                p = -2 * exp(-r*Z) * sqrt(Z^3) * (r*Z - 1)
            else #not implemented
                p = 0
            end
            push!( DPnl, p )
        end
        
        return( DPnl )
    end
    
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
    
    
    """
    `StrongField.pReducedME(epsilonp::Float64, lp::Int, n::Int, l::Int, epsiloni::Float64, volkov::AbstractVolkovState)`  
        ... computes the reduced matrix elements of the momentum operator <epsilonp lp ||p||n l> in the one-particle picture
    """
    function  pReducedME(epsilonp::Float64, lp::Int, n::Int, l::Int, epsiloni::Float64, volkov::AbstractVolkovState)
        rmax = 100.
        orderGL = 1000
        gaussLegendre = Radial.GridGL("Finite",0.0,rmax,orderGL)
        rgrid = gaussLegendre.t
        weights = gaussLegendre.wt
        
        Pnl = HydrogenPnl( epsiloni, n, l, rgrid )
        if  typeof(volkov) == FreeVolkov
            Pepsplp = VolkovP( epsilonp, lp, rgrid )
        elseif  typeof(volkov) == CoulombVolkov
            Pepsplp = CoulombVolkovP( epsilonp, lp, volkov.Z, rgrid )
        end
        
        DPnl = HydrogenDPnlDr( epsiloni, n, l, rgrid )
        
        integral = 0. * im
        
        #Sum over grid and compute Gauss-Legendre sum
        for    j = 1:orderGL
            r = rgrid[j]
            integrand = conj( Pepsplp[j] )/r * ( r*DPnl[j] - ((lp-l)*(lp+l+1))/2 * Pnl[j] )

            #Gauss-Legendre sum
            integral = integral + weights[j] * integrand
        end

        #Note that GSL.sf_coupling_3j takes takes the input (2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
        integral = integral * (-im)^(lp+1) * (-1)^lp * GSL.sf_coupling_3j( 2*lp, 2*1, 2*l, 0, 0, 0 )
        
        return( integral )
    end
    
    """
    `StrongField.scalarProdBoundCont(epsilonp::Float64, lp::Int, n::Int, l::Int, epsiloni::Float64, volkov::AbstractVolkovState)`  
        ... computes the scalar product of the bound (hydrogenic) and continuum (Volkov) states in the one-particle picture
    """
    function  scalarProdBoundCont(epsilonp::Float64, lp::Int, n::Int, l::Int, m::Int, epsiloni::Float64, volkov::AbstractVolkovState)
        rmax = 1000.
        orderGL = 1000
        gaussLegendre = Radial.GridGL("Finite",0.0,rmax,orderGL)
        rgrid = gaussLegendre.t
        weights = gaussLegendre.wt
        
        Pnl = HydrogenPnl( epsiloni, n, l, rgrid )
        
        if  typeof(volkov) == FreeVolkov
            Pepsplp = VolkovP( epsilonp, lp, rgrid )
        elseif  typeof(volkov) == CoulombVolkov
            Pepsplp = CoulombVolkovP( epsilonp, lp, volkov.Z, rgrid )
        end
        
        DPnl = HydrogenDPnlDr( epsiloni, n, l, rgrid )
        
        integral = 0. * im
        
        #Sum over grid and compute Gauss-Legendre sum
        for    j = 1:orderGL
            r = rgrid[j]
            integrand = conj( Pepsplp[j] ) * Pnl[j]

            #Gauss-Legendre sum
            integral = integral + weights[j] * integrand
        end
        
        return((-im)^l * integral)
    end
    
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
        
        #println("> Envelope integrals for an $(string(envelope)) with A0=$(beam.A0), omega=$(beam.omega) [a.u.], cep=$(beam.cep)" * 
        #        " and for $(string(polarization)) light are:" )
        #println("    F^(Volkov) [+; omega; f^(env); A]   = $fVolkovPlus" )
        #println("    F^(Volkov) [-; omega; f^(env); A]   = $fVolkovMinus" )
        #println("    F^(quad, Volkov) [f^(env); A]       = $fVolkovSquared" )
        
        return( (fVolkovPlus, fVolkovMinus, fVolkovSquared) )
    end


    """
    `StrongField.computeOutcome(observable::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,1})`  
        ... computes the requested photoelectron energy distribution and returns the results in the variable
            outcome::StrongField.OutcomeEnergyDistribution
    """
    function  computeOutcome(obs::StrongField.SfaEnergyDistribution, amplitudes::Array{SphericalAmplitude,1})
        probabilities = Float64[]
        for  amp  in  amplitudes   pp = sqrt(2*amp.energy);  push!(probabilities, pp * (amp.value * conj(amp.value)) )      end
        outcome = OutcomeEnergyDistribution(obs.theta, obs.phi, obs.energies, probabilities)
        return( outcome )
    end

    """
    `StrongField.computeOutcome(observable::StrongField.SfaMomentumDistribution, amplitudes::Array{SphericalAmplitude,1})`  
        ... computes the requested photoelectron momentum distribution and returns the results in the variable
            outcome::StrongField.OutcomeMomentumDistribution
    """
    function  computeOutcome(obs::StrongField.SfaMomentumDistribution, amplitudes::Array{SphericalAmplitude,1})
        probabilities = Float64[]
        momenta       = Array{Float64}(undef,0,2)
        for  amp  in  amplitudes   
            pp = sqrt(2*amp.energy)
            px = pp*sin(amp.theta)*cos(amp.phi)
            py = pp*sin(amp.theta)*sin(amp.phi)
            momenta = [momenta ; [px py]]
            push!( probabilities, pp * (amp.value * conj(amp.value)) )
        end
        outcome = OutcomeMomentumDistribution(obs.theta, momenta, probabilities)
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
        nqn = 1;    lqn = 0;    n = 1;      l = 0;      m = 0;    initialEn = convertUnits("energy: from eV to atomic", comp.initialLevel.energy)
        #
        
        lpMin = abs(l-1)
        lpMax = l+1
        
        # Compute the requested amplitudes
        for  amp in  sfaAmplitudes
            thetap        = amp.theta;   phip = amp.phi;     energyp = amp.energy
            envIntegrals  = StrongField.computeEnvelopeIntegrals(comp.envelope, comp.beam, comp.polarization, thetap, phip, energyp, initialEn)
            fVolkovPlus, fVolkovMinus, fVolkovSquared = envIntegrals
            #
            reducedMEArray = ComplexF64[]
            for lp = lpMin:lpMax
                push!( reducedMEArray, pReducedME(energyp, lp, n, l, initialEn, comp.volkov) )
            end
            #
                
            wminus = 0.0im;    wplus = 0.0im
            # Collect contributions from all l_p, q terms; this summation will change in a (nljm) representation
            for  lp = lpMin:lpMax
                reducedME = reducedMEArray[lp-lpMin+1]
                for  q in [-1, 0, 1]
                    wplus = wplus + AngularMomentum.sphericalYlm(lp, m-q, amp.theta, amp.phi) * (-1)^q  * 
                             Basics.determinePolarizationVector(q, comp.polarization, star=false)         *
                             AngularMomentum.ClebschGordan(l, m, 1, -q, lp, m-q) * reducedME
                    wminus = wminus + AngularMomentum.sphericalYlm(lp, m+q, amp.theta, amp.phi)           * 
                             Basics.determinePolarizationVector(q, comp.polarization, star=true)          *
                             AngularMomentum.ClebschGordan(l, m, 1, q, lp, m+q) * reducedME
                end
            end
            scalarProd = scalarProdBoundCont(energyp, l, n, l, m, initialEn, comp.volkov)
            wa = -im * sqrt(2/pi) * (fVolkovPlus * wminus + fVolkovMinus * wplus) - im / sqrt(2*pi) * fVolkovSquared * AngularMomentum.sphericalYlm(l, m, amp.theta, amp.phi) * scalarProd
            push!(newAmplitudes, SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa))
        
            println(">> $(SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa))")
        end
                
        return( newAmplitudes )
    end


    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaEnergyDistribution)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            A list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaEnergyDistribution)
        amplitudes = StrongField.SphericalAmplitude[]
        for  energy in observable.energies      push!(amplitudes, SphericalAmplitude(energy, observable.theta, observable.phi, 0.))     end
        
        println("> A total of $(length(amplitudes)) spherical amplitudes need to be calculated.")
        return( amplitudes )
    end
    
    
    """
    `StrongField.determineSphericalAmplitudes(observable::StrongField.SfaMomentumDistribution)`  
        ... determines which direct (and other) SFA amplitudes need to be computed; these amplitudes are not yet computed here but can be arranged
            so that only a minimum number of Volkov states and/or reduced many-electron matrix elements need to be computed.
            A list of amplitudes::Array{StrongField.SphericalAmplitude,1} is returned.
    """
    function  determineSphericalAmplitudes(observable::StrongField.SfaMomentumDistribution)
        amplitudes = StrongField.SphericalAmplitude[]
        for  energy in observable.energies
            for  phi in observable.phis
                push!(amplitudes, SphericalAmplitude(energy, observable.theta, phi, 0.))
            end
        end
        
        println("> A total of $(length(amplitudes)) spherical amplitudes need to be calculated.")
        return( amplitudes )
    end
    
    
    """
    `StrongField.plot( comp::StrongField.Computation, results::Dict{String,Any}, energyScale::String = "atomic", probabilityScaling::String = "linear" )`  
        ... generates a graphical representation of the observable (either StrongField.SfaEnergyDistribution or StrongField.SfaMomentumDistribution).
        energyScale determines the scaling of the energy axis either in atomic units (energyScale = "atomic") or in units of hbar*omega (energyScale = "omega").
        The y-axis is scaled either linearly (probabilityScaling = "linear") or logarithmically (probabilityScaling = "log")
    """
    function plot( comp::StrongField.Computation, results::Dict{String,Any}, energyScale::String = "atomic", probabilityScaling::String = "linear" )
    
        #---Photoelectron energy spectrum---
        if  typeof(comp.observable) == StrongField.SfaEnergyDistribution
                #prepare data and rescale the x-axis if neccessary
                if  energyScale == "atomic"
                    energyLabel = "E (a.u.)"
                    energies = results["energy distribution"].energies
                elseif  energyScale == "omega"
                    energyLabel = "E/omega"
                    energies = results["energy distribution"].energies / comp.beam.omega
                end
                probabilities = results["energy distribution"].probabilities
                gr() #sets the plotting backend to the package "GR"
                
                #set scaling of the y-axis
                scaling = :identity
                if  probabilityScaling == "log"
                            scaling = :log10
                end
                
                #generate the plot
                p = Plots.plot(energies, probabilities, 
                                title = "Photoelectron energy spectrum",
                                xlabel = energyLabel,
                                ylabel = "P(E)",
                                xscale = :identity,
                                yscale = scaling,
                                framestyle = :box,
                                legend = :none,
                                #markershape = :circle,
                                #markershape = :none,
                                line = 2,
                                linecolor = :black,
                                gridlinewidth = 2,
                                tickfontsize = 10,
                                labelfontsize = 10,
                                labelfontfamily = "Latin Modern Roman",
                                titlefontfamily = "Latin Modern Roman"
                                
                              )
                              
                #export the plot as png-file
                png("energy_spectrum")
        #-----------------------------------
        
        #---Photoelectron momentum distribution---
        elseif  typeof(comp.observable) == StrongField.SfaMomentumDistribution
                #prepare data
                pList = results["momentum distribution"].momenta
                probList = results["momentum distribution"].probabilities

                #generate the plot
                #pyplot()
                #r = range(0,stop=10,length=11)
                #theta = range(0,stop=360,length=361)
                #f(r,theta) = r^2
                #println("$(f.(r,theta'))")
                #Plots.plot( heatmap( f.(r,theta'), proj=:polar ) )
                #println("$(transpose(probList))")
                #Plots.plot( heatmap( pList, transpose(probList), proj=:polar ) )
                #imshow(pList,transpose(probList))
                #probList = probList[2:-1]
                #matplotlib.pyplot.pcolormesh(pList[:,1], pList[:,2], transpose(probList))
                contourf( pList[:,1], pList[:,2], probList )
                #matplotlib.pyplot.imshow(pList)
                #heatmap( pList[:,1], pList[:,2], probList, proj=:polar, legend=true
                                #heatmap( probList
                #                )

                #p = Plots.plot( #heatmap( pList[:,1], pList[:,2], probList
                #               Plots.GR.polarheatmap( probList 
                            ##heatmap( probList
                #                 )
                #              )
                
                
                #export the plot as png-file
                png("momentum_distribution")
        #-----------------------------------------
                
        #---Not a valid obserable---
        else     
                error("Undefined observable for strong-field computations.")
        end
        #---------------------------
    
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
            if output    results = Base.merge( results, Dict("computation" => comp, "energy distribution" => sfaOutcome) )  end
        elseif   typeof(comp.observable) == StrongField.SfaMomentumDistribution
            sfaAmplitudes = StrongField.computeSphericalAmplitudes(comp)
            sfaOutcome    = StrongField.computeOutcome(comp.observable, sfaAmplitudes)
            if output    results = Base.merge( results, Dict("computation" => comp, "momentum distribution" => sfaOutcome) )  end
        else     error("Undefined observable for strong-field computations.")
        end
        
        if  output   return( results )   else    return( nothing )   end
    end

end # module

