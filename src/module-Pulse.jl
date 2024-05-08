
"""
`module  JAC.Pulse`  
... a submodel of JAC that contains all structs and methods to deal with time-dependent pulses of the em field.
"""
module Pulse


using   ..Basics, ..Defaults, ..Radial, GSL, SpecialFunctions

export  AbstractEnvelope, AbstractBeam


"""
`abstract type Pulse.AbstractEnvelope` 
    ... defines an abstract type to comprise various envelopes of (possible) laser pulses in terms of their shape,  
        pulse duration or number of cycles, etc.

    + InfiniteEnvelope         ... to represent an infinte (plane-wave) pulse.
    + RectangularEnvelope      ... to represent a finite rectangular pulse.
    + SinSquaredEnvelope       ... to represent a finite sin^2 pulse.
    + GaussianEnvelope         ... to represent a Gaussian light pulse.
"""
abstract type  AbstractEnvelope  end


"""
`struct  Pulse.InfiniteEnvelope  <: Pulse.AbstractEnvelope`   ... to represent an infinte (plane-wave) pulse.
"""
struct   InfiniteEnvelope        <: Pulse.AbstractEnvelope    end

function Base.string(env::InfiniteEnvelope)
    sa = "infinite pulse"
    return( sa )
end

function Base.show(io::IO, env::InfiniteEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`struct  Pulse.RectangularEnvelope  <: Pulse.AbstractEnvelope`   ... to represent a finte rectangular pulse.

    + cycles      ::Int64     ... Number of cycles of the pulse.
"""
struct   RectangularEnvelope       <: Pulse.AbstractEnvelope
    cycles        ::Int64
end

function Base.string(env::RectangularEnvelope)
    sa = "rectangular pulse of $(env.cycles) cycles"
    return( sa )
end

function Base.show(io::IO, env::RectangularEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`struct  Pulse.SinSquaredEnvelope  <: Pulse.AbstractEnvelope`   ... to represent a finite sin^2 pulse.

    + cycles      ::Int64     ... Number of cycles of the pulse.
"""
struct   SinSquaredEnvelope        <: Pulse.AbstractEnvelope
    cycles      ::Int64
end


function Base.string(env::SinSquaredEnvelope)
    sa = "sin^2 pulse of $(env.cycles) cycles"
    return( sa )
end

function Base.show(io::IO, env::SinSquaredEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`struct  Pulse.GaussianEnvelope  <: Pulse.AbstractEnvelope`   ... to represent a (infinite) Gaussian pulse.

    + fwhm        ::Float64     ... FWHM which is often taken as pulse duration
"""
struct   GaussianEnvelope        <: Pulse.AbstractEnvelope
    fwhm          ::Float64
end


function Base.string(env::GaussianEnvelope)
    sa = "Gaussian pulse with pulse duration or FWHM = $(env.fwhm)"
    return( sa )
end

function Base.show(io::IO, env::GaussianEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`abstract type Pulse.AbstractBeam` 
    ... defines an abstract type to comprise various basic laser pulses that are characterized in terms of their amplitude, 
        frequency, carrier-envelope phase, etc. In general, the basic beam properties are independent of the (pulse) envelope 
        and the polarization properties which are handled and communicated separately (to and within the program).

    + PlaneWaveBeam            ... to represent a plane-wave beam.
    + BesselBeam               ... to represent a Bessel beam (not yet).
"""
abstract type  AbstractBeam  end


"""
`struct  Pulse.PlaneWaveBeam  <: AbstractBeam`   
    ... to represent a plane-wave beam with given (real) amplitude, frequency and carrier-envelope phase.

    + A0            ::Float64            ... (Constant) Amplitude of the light pulse. 
    + omega         ::Float64            ... Central frequency. 
    + cep           ::Float64            ... Carrier-envelope phase. 
"""
struct   PlaneWaveBeam  <: AbstractBeam
    A0              ::Float64
    omega           ::Float64
    cep             ::Float64
end


"""
`Pulse.PlaneWaveBeam()`  ... constructor for an `empty` instance of Pulse.PlaneWaveBeam().
"""
function PlaneWaveBeam()
    PlaneWaveBeam( 0., 0., 0. )
end


function Base.string(beam::PlaneWaveBeam)
    sa = "Plane-wave beam/pulse with amplitude A0=$(beam.A0), frequency omega=$(beam.omega) a.u. and carrier-envelope phase cep=$(beam.cep)"
    return( sa )
end

function Base.show(io::IO, beam::PlaneWaveBeam)
    sa = string(beam);       print(io, sa, "\n")
end


#====================================================================================== 
#  These data types need to be worked through

"""
`struct  Pulse.Envelope`  ... defines a type for the envelope (function) of an em pulse.

    + timeDelay          ::Float64             ... Time delay with regard to an arbitrary time t=0, often used for the first pulse.
    + f0                 ::Float64             ... Amplitude of the envelope function f(t) = f0 * f_shape(t)
    + shapeFunction      ::Array{Float64,1}    ... (Normalized) shape function of the pulse.
"""
struct Envelope
    timeDelay            ::Float64
    f0                   ::Float64
    shapeFunction        ::Array{Float64,1}
end 


"""
`Pulse.Envelope()`  ... constructor for an `empty` instance of Pulse.Envelope().
"""
function Envelope()
    Envelope( 0., 0., Float64[] )
end


# `Base.show(io::IO, envelope::Pulse.Envelope)`  ... prepares a proper printout of the variable envelope::Pulse.Envelope.
function Base.show(io::IO, envelope::Pulse.Envelope) 
    println(io, "timeDelay:             $(envelope.timeDelay)  ")
    println(io, "f0:                    $(envelope.f0)  ")
    println(io, "shapeFunction:         $(envelope.shapeFunction)  ")
end


"""
`@enum   PolarizationType`  ... defines a type of polarization of an experimentally described light pulse

+ NoType, Linear, LeftCircular, RightCircular, Elliptical      ... with obvious meaning.
"""
@enum   PolarizationType    NoType    Linear    LeftCircular    RightCircular    Elliptical


"""
`Pulse.PolarizationType(sa::String)`  ... constructor for a given String.
"""
function PolarizationType(sa::String)
    if       sa == "linear"                  wa = Linear
    elseif   sa == "left-circular"           wa = LeftCircular
    elseif   sa == "right-circular"          wa = RightCircular
    elseif   sa == "elliptical"              wa = Elliptical
    else
        error("stop a")
    end
    PolarizationType(wa)
end


# `Base.show(io::IO, ptyp::Pulse.PolarizationType)`  ... prepares a proper printout of the variable ptyp::Pulse.PolarizationType.
function Base.show(io::IO, ptyp::Pulse.PolarizationType) 
    print(io, string(ptyp) )
end


"""
`Base.string(ptyp::Pulse.PolarizationType)`  ... provides a proper printout of the variable ptyp::Pulse.PolarizationType.
"""
function Base.string(ptyp::Pulse.PolarizationType) 
    if      ptyp == NoType              return( "undefined polarization" )
    elseif  ptyp == Linear              return( "linear polarization" )
    elseif  ptyp == LeftCircular        return( "left-circular polarization" )
    elseif  ptyp == RightCircular       return( "right-circular polarization" )
    elseif  ptyp == Elliptical          return( "elliptical polarization" )
    elseif  ptyp == Linear              return( "linear polarization" )
    else    error("stop a")
    end
end


"""
`struct  Pulse.Polarization`  ... defines a type for the polarization of an em pulse.

    + type               ::Pulse.PolarizationType    ... General type of polarization.
    + circularDegree     ::Float64                   ... Degree of circular polarization.
    + linearDegree       ::Float64                   ... Degree of linear polarization.
    + linearPolarization ::SolidAngle                ... Direction of the (linear) polarization vector as described by the
                                                            unit vector u = u(Omega') with regard to the x'-z' plane of the incident pulse.
    + has_gCoefficients  ::Bool                      ... True, of the complex gPlus and gMinus coefficients are properly defined.
    + gPlus              ::Complex{Float64}          ... g_+1 coefficient.    
    + gMinus             ::Complex{Float64}          ... g_-1 coefficient.    
"""
struct Polarization
    typ                  ::Pulse.PolarizationType 
    circularDegree       ::Float64
    linearDegree         ::Float64
    linearPolarization   ::SolidAngle
    has_gCoefficients    ::Bool
    gPlus                ::Complex{Float64} 
    gMinus               ::Complex{Float64}
end 


"""
`Pulse.Polarization()`  ... constructor for an `empty` instance of Pulse.Polarization().
"""
function Polarization()
    Polarization( NoType, 0., 0., SolidAngle(0., 0.), false, Complex(0.), Complex(0.) )
end


# `Base.show(io::IO, polarization::Pulse.Polarization)`  ... prepares a proper printout of the variable polarization::Pulse.Polarization.
function Base.show(io::IO, polarization::Pulse.Polarization) 
    println(io, "type:                    $(polarization.typ)  ")
    println(io, "circularDegree:          $(polarization.circularDegree)  ")
    println(io, "linearDegree:            $(polarization.linearDegree)  ")
    println(io, "linearPolarization:      $(polarization.linearPolarization)  ")
    println(io, "has_gCoefficients:       $(polarization.has_gCoefficients)  ")
    println(io, "gPlus:                   $(polarization.gPlus)  ")
    println(io, "gMinus:                  $(polarization.gMinus)  ")
end


"""
`struct  Pulse.ExperimentalCharacterization`  ... defines a type for characterizing an experimental or physically described em pulse.

    + kind             ::Pulse.Shape            ... General shape of a given pulse.
    + propagation      ::SolidAngle             ... Propagation direction of the pulse as described by the unit vector u = u(Omega).
    + omega            ::Float64                ... Central frequency of the em pulse.
    + maxIntensity     ::Float64                ... maximum intensity in units ...
    + fwhm             ::Float64                ... approximate full-width-half maximum (fwhm) for rather long pulses; the nearest widths
                                                    with an integer No. of cycles is taken.
    + timeDelay        ::Float64                ... time delay with regard to some arbitrarely chosen time t = 0 
                                                    (often taken for the first pulse).
    + NoCycles         ::Int64                  ... No. of cycles which always determines the lengths of the pulse.
    + multipoles       ::Array{EmMultipole,1}   ... Multipoles of the em field to be included in the description of the light field.
    + polarization     ::Pulse.Polarization     ... Describes the polarization properties of the em pulse.
"""
struct ExperimentalCharacterization
    kind               ::Pulse.Shape
    propagation        ::SolidAngle
    omega              ::Float64
    maxIntensity       ::Float64
    fwhm               ::Float64
    timeDelay          ::Float64
    NoCycles           ::Int64 
    multipoles         ::Array{EmMultipole,1}
    polarization       ::Pulse.Polarization
end 


"""
`Pulse.ExperimentalCharacterization()`  ... constructor for an `empty` instance of ExperimentalCharacterization().
"""
function ExperimentalCharacterization()
    ExperimentalCharacterization( NoShape, SolidAngle(0., 0.), 0., 0., 0., 0., 0, EmMultipole[], Polarization() )
end


# `Base.show(io::IO, chz::Pulse.ExperimentalCharacterization)`  
#		... prepares a proper printout of the variable  chz::Pulse.ExperimentalCharacterization.
function Base.show(io::IO, characterization::Pulse.ExperimentalCharacterization) 
    println(io, "kind:              $(chz.kind)  ")
    println(io, "propagation:       $(chz.propagation)  ")
    println(io, "omega:             $(chz.omega)  ")
    println(io, "maxIntensity:      $(chz.maxIntensity)  ")
    println(io, "fwhm:              $(chz.fwhm)  ")
    println(io, "timeDelay:         $(chz.timeDelay)  ")
    println(io, "NoCycles:          $(chz.NoCycles)  ")
    println(io, "multipoles:        $(chz.multipoles)  ")
    println(io, "polarization:      $(chz.polarization)  ")
end


"""
`struct  Pulse.Gaussian`  ... defines a type for a Gaussian light pulse that is used for evaluating time-dependent statistical tensors.

    + propagation      ::SolidAngle             ... Propagation direction of the pulse as described by the unit vector u = u(Omega).
    + omega            ::Float64                ... Central frequency of the em pulse.
    + multipoles       ::Array{EmMultipole,1}   ... Multipoles of the em field to be included in the description of the light field.
    + envelope         ::Pulse.Envelope         ... Envelope (function) of the light pulse.
    + polarization     ::Pulse.Polarization     ... Polarization of the light pulse with typically well-defined g_+1 and g_-1 coefficients.
"""
struct Gaussian
    propagation        ::SolidAngle 
    omega              ::Float64 
    multipoles         ::Array{EmMultipole,1}
    envelope           ::Pulse.Envelope 
    polarization       ::Pulse.Polarization
end 


"""
`Pulse.Gaussian()`  ... constructor for an `empty` instance of Pulse.Gaussian().
"""
function Gaussian()
    Gaussian( SolidAngle(0., 0.), 0., EmMultipole[], Pulse.Envelope(), Pulse.Polarization()  )
end


# `Base.show(io::IO, gaussian::Pulse.Gaussian)`  ... prepares a proper printout of the variable gaussian::Pulse.Gaussian.
function Base.show(io::IO, gaussian::Pulse.Gaussian) 
    println(io, "propagation:        $(gaussian.propagation)  ")
    println(io, "omega:              $(gaussian.omega)  ")
    println(io, "multipoles:         $(gaussian.multipoles)  ")
    println(io, "envelope:           $(gaussian.envelope)  ")
    println(io, "polarization:       $(gaussian.polarization)  ")
end
============================================================================================#

###################################################################################################################
###################################################################################################################
###################################################################################################################



"""
`Pulse.computeFieldAmplitude(intensity::Float64, omega::Float64)`  
    ... compute the field amplitude from the given (maximum) intensity [in a.u.] and the central frequency [in a.u.]
        of the light field: A0 = sqrt( 8pi * alpha * intensity) / omega.
"""
function computeFieldAmplitude(intensity::Float64, omega::Float64)
    wa = sqrt(8 * pi * Defaults.getDefaults("alpha") * intensity) / omega
    return( wa )
end


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(infinite); A; angles & energies] for an infinite pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;                            
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) / beam.omega;   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    #
    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -10:10
        if  plus   wb = Basics.diracDelta((s-1)*beam.omega + energyp - initialEn + Up, 1.0e-3)
        else       wb = Basics.diracDelta((s+1)*beam.omega + energyp - initialEn + Up, 1.0e-3)       
        end
        #
        if  wb != 0.
            wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * (phiCep - lambda*phip)) * wb
        end
    end
    if  plus   wa = wa * 2pi * beam.A0 * exp(-im*phiCep)
    else       wa = wa * 2pi * beam.A0 * exp(im*phiCep)
    end
end
    
    return( wa )
end


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for an infinite pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   Tp = 2pi * envelope.cycles / beam.omega
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / ( sqrt(2) * beam.omega );   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    #
    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -20:20
        if  plus   wb = (s-1)*beam.omega + energyp - initialEn + Up
        else       wb = (s+1)*beam.omega + energyp - initialEn + Up       
        end
        #
        if  wb != 0.
            wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * phaseConstant) / wb * ( exp(im*wb*Tp) - 1 )
        end
    end
    if  plus   wa = wa * (-im) * exp(-im * a * sin(phaseConstant)) * beam.A0 * exp(-im*phiCep)
    else       wa = wa * (-im) * exp(-im * a * sin(phaseConstant)) * beam.A0 * exp(im*phiCep)
    end
    end
    
    return( wa )
end


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for a sine-squared pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   np = envelope.cycles;   Tp = 2pi * np / beam.omega
    omega = beam.omega
    phiCep = beam.cep;
    sinSqrArg = 0.5 * omega / np;
    lambda = Basics.determinePolarizationLambda(polarization)
    
    p = sqrt(2.0*energyp)
    px = p*sin(thetap)*cos(phip)
    py = p*sin(thetap)*sin(phip)
    
    epsilon = 1.0
    if polarization != Basics.RightCircular() && polarization != Basics.LeftCircular()
        epsilon = polarization.ellipticity
    end
    
    A0eps = beam.A0/sqrt(1.0 + epsilon^2)
    
    #Define Gauss-Legendre grid, convergence is typically good for orderGL = 100 * np (time consuming for np > 10); tested up to np = 20
    if  np <= 10     orderGL = 100*np
    else             orderGL = 1000
    end
    gaussLegendre = Radial.GridGL("Finite",0.0,Tp,orderGL)
    tgrid = gaussLegendre.t
    weights = gaussLegendre.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        t = tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        cosIntegral = 0.25 / (omega * (np^2-1)) * (  2*sin(phiCep) + 2 * (np^2-1) * sin(phiCep + omega*t) - np * ( (1+np)*sin( phiCep + (np-1)/np * omega*t ) + (np-1) * sin( phiCep + (np+1)/np * omega*t )  ) )
        
        sinIntegral = 0.25 / (omega * (np^2-1)) * ( -2*cos(phiCep) - 2 * (np^2-1) * cos(phiCep + omega*t) + np * ( (1+np)*cos( phiCep + (np-1)/np * omega*t ) + (np-1) * cos( phiCep + (np+1)/np * omega*t )  ) )
        
        cos2Integral = sin(2*phiCep)/omega * ( -6 - np/(np-1) - np/(np+1) + 8*np/(2*np-1) + 8*np/(2*np+1) )
                        + 12*t + 6/omega * cos(2*omega*t) * sin(2*phiCep) + 6/omega * cos(2*phiCep) * sin(2*omega*t) - 16/omega*np*sin(omega*t/np) + 2/omega*np*sin(2*omega*t/np)
                        - 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t) + np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t))
                        + np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np * omega*t)) - 8*np/(omega*(2*np-1))*sin(2*phiCep + (2*np-1)/np * omega*t )
        cos2Integral = cos2Integral / 64
        
        sin2Integral = 12*t + 6/omega * sin(2*phiCep) * ( 1/(1-5*np^2+4*np^4) - cos(2*omega*t) ) - 6/omega * cos(2*phiCep)*sin(2*omega*t)
                        - 16/omega * np * sin(omega*t/np) + 2/omega * np * sin(2*omega*t/np) + 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t)
                        - np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t )) - np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np*omega*t)) + 8/omega * np/(2*np-1) * sin(2*phiCep + (2*np-1)/np*omega*t)
        sin2Integral = sin2Integral / 64
        
        SVolkov = energyp*t + A0eps*px*cosIntegral + A0eps*lambda*epsilon*py*sinIntegral + 0.5 * A0eps^2 * ( cos2Integral + epsilon^2 * sin2Integral )
        
        #Compute integrand at gridpoint t
        if  plus    integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( initialEn + beam.omega ) * t - SVolkov ) )
        else        integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( initialEn - beam.omega ) * t - SVolkov ) )
        end
        
        #Gauss-Legendre sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    if  plus    wa = beam.A0 * exp(-im * phiCep) * wa
    else        wa = beam.A0 * exp(im * phiCep) * wa
    end
    end
    
    return( wa )
end    


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for a gaussian pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int )
if orderSFA == 0
    wa = 0. * im;  Tp = envelope.fwhm
    phiCep = beam.cep;   Up = beam.A0^2 / 4
    a1 = 0.25 * Up * sqrt(pi/log(4)) * Tp
    a2 = 0.125 * Tp * sqrt(pi/log(2)) * beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) * exp( -beam.omega^2 * Tp^2 / log( 65536 ) )
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    
    transformFac = Tp/(2*sqrt(log(2))) #factor from (linear) variable transformation tau -> t in order to allow Gauss-Hermite integration ( weight function e^(-t^2) )
    
    #Define Gauss-Hermite grid, convergence is typically good for orderGH = 10000
    orderGH = 10000
    gaussHermite = Radial.GridGH(orderGH)
    tgrid = gaussHermite.t
    weights = gaussHermite.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGH
        t = transformFac * tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        SVolkov = energyp * t + a1 * erf(2 * sqrt(log(4)) * t/Tp) 
                    + a2 * ( exp( -im*phaseConstant ) * erf( (im*Tp^2*beam.omega + 8*log(2)*t)/(4*Tp*sqrt(log(2))) )
                                - exp( im*phaseConstant ) * erf( (im*Tp^2*beam.omega - 8*log(2)*t)/(4*Tp*sqrt(log(2))) ) )
        
        #Compute integrand at gridpoint t
        if  plus    integrand = exp( -im * ( ( initialEn + beam.omega ) * t - SVolkov ) )
        else        integrand = exp( -im * ( ( initialEn - beam.omega ) * t - SVolkov ) )
        end
        
        #Gauss-Hermite sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    if  plus    wa = transformFac * beam.A0 * exp(-im * phiCep) * wa
    else        wa = transformFac * beam.A0 * exp(im * phiCep) * wa
    end
end
    
    return( wa )
end  


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderGH::Int64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for a gaussian pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderGH::Int64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;  Tp = envelope.fwhm
    phiCep = beam.cep;   Up = beam.A0^2 / 4
    a1 = 0.25 * Up * sqrt(pi/log(4)) * Tp
    a2 = 0.125 * Tp * sqrt(pi/log(2)) * beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) * exp( -beam.omega^2 * Tp^2 / log( 65536 ) )
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    
    transformFac = Tp/(2*sqrt(log(2))) #factor from (linear) variable transformation tau -> t in order to allow Gauss-Hermite integration ( weight function e^(-t^2) )
    
    gaussHermite = Radial.GridGH(orderGH)
    tgrid = gaussHermite.t
    weights = gaussHermite.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGH
        t = transformFac * tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        SVolkov = energyp * t + a1 * erf(2 * sqrt(log(4)) * t/Tp) 
                    + a2 * ( exp( -im*phaseConstant ) * erf( (im*Tp^2*beam.omega + 8*log(2)*t)/(4*Tp*sqrt(log(2))) )
                                - exp( im*phaseConstant ) * erf( (im*Tp^2*beam.omega - 8*log(2)*t)/(4*Tp*sqrt(log(2))) ) )
        
        #Compute integrand at gridpoint t
        if  plus    integrand = exp( -im * ( ( initialEn + beam.omega ) * t - SVolkov ) )
        else        integrand = exp( -im * ( ( initialEn - beam.omega ) * t - SVolkov ) )
        end
        
        #Gauss-Hermite sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    if  plus    wa = transformFac * beam.A0 * exp(-im * phiCep) * wa
    else        wa = transformFac * beam.A0 * exp(im * phiCep) * wa
    end
end
    
    return( wa )
end  


"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies] for an infinite pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
if orderSFA == 0
    wa = 0. * im;                            
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) / beam.omega;   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    #
    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -10:10
        wb = Basics.diracDelta(s*beam.omega + energyp - initialEn + Up, 1.0e-3)
        #
        if  wb != 0.    wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * (phiCep - lambda*phip)) * wb     
        end
    end
    wa = wa * 4pi * Up
end
return( wa )
end


"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies]  for a rectangular pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   Tp = 2pi * envelope.cycles / beam.omega
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) / beam.omega;   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip

    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -20:20
        wb = s*beam.omega + energyp - initialEn + Up
        #
        if  wb != 0.    wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * phaseConstant) / wb * ( exp(im*wb*Tp) - 1 )   
        end
    end
    wa = wa * (-im) * 2 * Up * exp(-im * a * sin(phaseConstant))
    end
    
    return( wa )
end


"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies] for a sine-squared pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   np = envelope.cycles;   Tp = 2pi * np / beam.omega
    omega = beam.omega
    phiCep = beam.cep
    sinSqrArg = 0.5 * omega / np
    lambda = Basics.determinePolarizationLambda(polarization)
    
    p = sqrt(2.0*energyp)
    px = p*sin(thetap)*cos(phip)
    py = p*sin(thetap)*sin(phip)
    
    epsilon = 1.0
    if polarization != Basics.RightCircular() && polarization != Basics.LeftCircular()
        epsilon = polarization.ellipticity
    end
    
    A0eps = beam.A0/sqrt(1.0 + epsilon^2)
    
    #Define Gauss-Legendre grid, convergence is typically good for orderGL = 100 * np (time consuming for np > 10); tested up to np = 20
    if  np <= 10     orderGL = 100*np
    else             orderGL = 1000
    end
    gaussLegendre = Radial.GridGL("Finite",0.0,Tp,orderGL)
    tgrid = gaussLegendre.t
    weights = gaussLegendre.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        t = tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        cosIntegral = 0.25 / (omega * (np^2-1)) * (  2*sin(phiCep) + 2 * (np^2-1) * sin(phiCep + omega*t) - np * ( (1+np)*sin( phiCep + (np-1)/np * omega*t ) + (np-1) * sin( phiCep + (np+1)/np * omega*t )  ) )
        
        sinIntegral = 0.25 / (omega * (np^2-1)) * ( -2*cos(phiCep) - 2 * (np^2-1) * cos(phiCep + omega*t) + np * ( (1+np)*cos( phiCep + (np-1)/np * omega*t ) + (np-1) * cos( phiCep + (np+1)/np * omega*t )  ) )
        
        cos2Integral = sin(2*phiCep)/omega * ( -6 - np/(np-1) - np/(np+1) + 8*np/(2*np-1) + 8*np/(2*np+1) )
                        + 12*t + 6/omega * cos(2*omega*t) * sin(2*phiCep) + 6/omega * cos(2*phiCep) * sin(2*omega*t) - 16/omega*np*sin(omega*t/np) + 2/omega*np*sin(2*omega*t/np)
                        - 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t) + np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t))
                        + np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np * omega*t)) - 8*np/(omega*(2*np-1))*sin(2*phiCep + (2*np-1)/np * omega*t )
        cos2Integral = cos2Integral / 64
        
        sin2Integral = 12*t + 6/omega * sin(2*phiCep) * ( 1/(1-5*np^2+4*np^4) - cos(2*omega*t) ) - 6/omega * cos(2*phiCep)*sin(2*omega*t)
                        - 16/omega * np * sin(omega*t/np) + 2/omega * np * sin(2*omega*t/np) + 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t)
                        - np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t )) - np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np*omega*t)) + 8/omega * np/(2*np-1) * sin(2*phiCep + (2*np-1)/np*omega*t)
        sin2Integral = sin2Integral / 64
        
        SVolkov = energyp*t + A0eps*px*cosIntegral + A0eps*lambda*epsilon*py*sinIntegral + 0.5 * A0eps^2 * ( cos2Integral + epsilon^2 * sin2Integral )
        
        #Compute integrand at gridpoint t
        integrand = sin( sinSqrArg * t )^4 * ( cos(omega*t+phiCep)^2 + epsilon^2 * sin(omega*t+phiCep)^2 ) * exp( -im * ( initialEn * t - SVolkov ) )
        
        #Gauss-Legendre sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    wa = A0eps^2 * wa
    end
    
    return( wa )
end



"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)` 
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies] for a gaussian pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;  Tp = envelope.fwhm
    phiCep = beam.cep;   Up = beam.A0^2 / 4
    a1 = 0.25 * Up * sqrt(pi/log(4)) * Tp
    a2 = 0.125 * Tp * sqrt(pi/log(2)) * beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) * exp( -beam.omega^2 * Tp^2 / log( 65536 ) )
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    
    transformFac = Tp/(2*sqrt(2*log(2))) #factor from (linear) variable transformation tau -> t in order to allow Gauss-Hermite integration ( weight function e^(-t^2) )
    
    #Define Gauss-Hermite grid, convergence is typically good for ......
    orderGH = 10000
    gaussHermite = Radial.GridGH(orderGH)
    tgrid = gaussHermite.t
    weights = gaussHermite.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGH
        t = transformFac * tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        SVolkov = energyp * t + a1 * erf(2 * sqrt(log(4)) * t/Tp) 
                    + a2 * ( exp( -im*phaseConstant ) * erf( (im*Tp^2*beam.omega + 8*log(2)*t)/(4*Tp*sqrt(log(2))) )
                                - exp( im*phaseConstant ) * erf( (im*Tp^2*beam.omega - 8*log(2)*t)/(4*Tp*sqrt(log(2))) ) )
        
        #Compute integrand at gridpoint t
        integrand = exp( -im * ( initialEn * t - SVolkov ) )
        
        #Gauss-Hermite sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    wa = 0.5 * transformFac * beam.A0^2 * wa
    end
    
    return( wa )
end  

"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.AbstractEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies]  for all pulses for which no analytical
        expression is so easily available.an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.AbstractEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   
    # Collect parameters that are specific to a given pulse envelope
    if       typeof(envelope) == SinSquaredEnvelope
    elseif   typeof(envelope) == GaussianEnvelope
    end
    #
    # Determine first the integrant; the timeGrid must still be adapted to thos
    timeGrid = [0.1i for i = 1:10]
    for   t  in  timeGrid 
        phase = Pulse.volkovPhase(t, envelope) ## , ...)
    end
    #
    error("Not yet implemented.")
end
return( wa )
end

end # module
