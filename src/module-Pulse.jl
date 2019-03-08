
"""
`module  JAC.Pulse`  ... a submodel of JAC that contains all structs and methods to deal with time-dependent pulses of the em field; 
                         it is using JAC.
"""
module Pulse

    using JAC


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
    `JAC.Pulse.Envelope()`  ... constructor for an `empty` instance of Pulse.Envelope().
    """
    function Envelope()
        Envelope( 0., 0., Float64[] )
    end


    """
    `Base.show(io::IO, envelope::Pulse.Envelope)`  ... prepares a proper printout of the variable envelope::Pulse.Envelope.
    """
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
    `JAC.Pulse.PolarizationType(sa::String)`  ... constructor for a given String.
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


    """
    `Base.show(io::IO, ptyp::Pulse.PolarizationType)`  ... prepares a proper printout of the variable ptyp::Pulse.PolarizationType.
    """
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
    `@enum   Shape`  ... defines a shape of a general em pulse

    + Gaussian               ... with obvious meaning
    """
    @enum   Shape    NoShape    GaussianX   SinSquared    SinSquaredCycles


    """
    `JAC.Pulse.Shape(sa::String)`  ... constructor for a given String.
    """
    function Shape(sa::String)
        if       sa == "Gaussian"                wa = GaussianX
        elseif   sa in ["sin^2 (fwhm)"]          wa = SinSquared
        elseif   sa in ["sin^2 (No. cycles)"]    wa = SinSquaredCycles
        else
            error("stop a")
        end
        Shape(wa)
    end


    """
    `Base.show(io::IO, shape:Pulse.Shape)`  ... prepares a proper printout of the variable shape::Pulse.Shape.
    """
    function Base.show(io::IO, shape::Pulse.Shape) 
        print(io, string(shape) )
    end


    """
    `Base.string(shape:Pulse.Shape)`  ... provides a proper printout of the variable shape::Pulse.Shape.
    """
    function Base.string(shape::Pulse.Shape) 
        if      shape == Shape              return( "undefined pulse" )
        elseif  shape == GaussianX          return( "Gaussian pulse" )
        elseif  shape == SinSquared         return( "sin^2 (fwhm) pulse" )
        elseif  shape == SinSquaredCycles   return( "sin^2 (No. cycles) pulse" )
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
    `JAC.Pulse.Polarization()`  ... constructor for an `empty` instance of Pulse.Polarization().
    """
    function Polarization()
        Polarization( NoType, 0., 0., SolidAngle(0., 0.), false, Complex(0.), Complex(0.) )
    end


    """
    `Base.show(io::IO, polarization::Pulse.Polarization)`  ... prepares a proper printout of the variable polarization::Pulse.Polarization.
    """
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
    `JAC.Pulse.ExperimentalCharacterization()`  ... constructor for an `empty` instance of ExperimentalCharacterization().
    """
    function ExperimentalCharacterization()
        ExperimentalCharacterization( NoShape, SolidAngle(0., 0.), 0., 0., 0., 0., 0, EmMultipole[], Polarization() )
    end


    """
    `Base.show(io::IO, chz::Pulse.ExperimentalCharacterization)`  ... prepares a proper printout of the variable
                                                                      chz::Pulse.ExperimentalCharacterization.
    """
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
    `JAC.Pulse.Gaussian()`  ... constructor for an `empty` instance of Pulse.Gaussian().
    """
    function Gaussian()
        Gaussian( SolidAngle(0., 0.), 0., EmMultipole[], Pulse.Envelope(), Pulse.Polarization()  )
    end


    """
    `Base.show(io::IO, gaussian::Pulse.Gaussian)`  ... prepares a proper printout of the variable gaussian::Pulse.Gaussian.
    """
    function Base.show(io::IO, gaussian::Pulse.Gaussian) 
        println(io, "propagation:        $(gaussian.propagation)  ")
        println(io, "omega:              $(gaussian.omega)  ")
        println(io, "multipoles:         $(gaussian.multipoles)  ")
        println(io, "envelope:           $(gaussian.envelope)  ")
        println(io, "polarization:       $(gaussian.polarization)  ")
    end

end # module
