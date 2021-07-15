
"""
`module  JAC.ReducedDensityMatrix`  
    ... a submodel of JAC that contains all methods for computing reduced density matrices, natural orbitals, radial density
        distributions and related information for some level(s).
"""
module ReducedDensityMatrix

    using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
                  ..SpinAngular, ..TableStrings

    """
    `struct  ReducedDensityMatrix.Settings  <:  AbstractPropertySettings`  
        ... defines a type for providing and parameters for reduced density matrices and other parameters.

        + calc1pRDM                ::Bool             ... True if 1-particle RDM need to be calculated, and false otherwise.
        + calc2pRDM                ::Bool             ... True if 2-particle RDM need to be calculated, and false otherwise.
        + calcCumulant             ::Bool             ... True if the cumulant need to be calculated, and false otherwise.
        + calcNatural              ::Bool             ... True if natural orbitals need to be calculated, and false otherwise.
        + calcDensity              ::Bool             ... True if the radial density need to be calculated, and false otherwise.
        + calcIpq                  ::Bool             ... True if orbital interaction I_pq need to be calculated, and false otherwise.
        + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
        + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractPropertySettings 
        calc1pRDM                  ::Bool 
        calc2pRDM                  ::Bool 
        calcCumulant               ::Bool
        calcNatural                ::Bool 
        calcDensity                ::Bool 
        calcIpq                    ::Bool 
        printBefore                ::Bool 
        levelSelection             ::LevelSelection
    end 


    """
    `ReducedDensityMatrix.Settings(; calc1pRDM::Bool=true,` calc2pRDM::Bool=false, calcCumulant::Bool=false, calcNatural::Bool=false, 
                                     calcDensity::Bool=false, calcIpq::Bool=false, printBefore::Bool=true, 
                                     levelSelection::LevelSelection=LevelSelection()) 
        ... keyword constructor to overwrite selected value of reduced-density matrix computations.
    """
    function Settings(; calc1pRDM::Bool=true, calc2pRDM::Bool=false, calcCumulant::Bool=false, calcNatural::Bool=false, 
                        calcDensity::Bool=false, calcIpq::Bool=false, printBefore::Bool=true, 
                        levelSelection::LevelSelection=LevelSelection())
        Settings(calc1pRDM, calc2pRDM, calcCumulant, calcNatural, calcDensity, calcIpq, printBefore, levelSelection)
    end

    # `Base.show(io::IO, settings::ReducedDensityMatrix.Settings)`  
    #       ... prepares a proper printout of the variable settings::ReducedDensityMatrix.Settings.
    function Base.show(io::IO, settings::ReducedDensityMatrix.Settings) 
        println(io, "calc1pRDM:                $(settings.calc1pRDM)  ")
        println(io, "calc2pRDM:                $(settings.calc2pRDM)  ")
        println(io, "calcCumulant:             $(settings.calcCumulant)  ")
        println(io, "calcNatural:              $(settings.calcNatural)  ")
        println(io, "calcDensity:              $(settings.calcDensity)  ")
        println(io, "calcIpq:                  $(settings.calcIpq)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "levelSelection:           $(settings.levelSelection)  ")
    end


    """
    `struct  ReducedDensityMatrix.Outcome`  
        ... defines a type to keep the outcome of a reduced density-matrix and natural-orbital computation.

        + level                     ::Level              ... Atomic level to which the outcome refers to.
        + Knms                      ::Float64            ... K_nms parameter
        + Ksms                      ::Float64            ... K_sms parameter
        + Fme                       ::Float64            ... F parameter from matrix element.
        + Fdensity                  ::Float64            ... F parameter from density.
        + Xboson                    ::Float64            ... X boson-field shift constant.
        + amplitudeKnms             ::Complex{Float64}   ... K_nms amplitude
        + amplitudeKsmsA            ::Complex{Float64}   ... K_sms,A amplitude
        + amplitudeKsmsB            ::Complex{Float64}   ... K_sms,B amplitude
        + amplitudeKsmsC            ::Complex{Float64}   ... K_sms,C amplitude
    """
    struct Outcome 
        level                       ::Level 
        Knms                        ::Float64
        Ksms                        ::Float64
        Fme                         ::Float64
        Fdensity                    ::Float64
        Xboson                      ::Float64
        amplitudeKnms               ::Complex{Float64}
        amplitudeKsmsA              ::Complex{Float64}
        amplitudeKsmsB              ::Complex{Float64}
        amplitudeKsmsC              ::Complex{Float64}
    end 


    """
    `ReducedDensityMatrix.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
    """
    function Outcome()
        Outcome(Level(), 0., 0., 0., 0., 0.,   0., 0., 0., 0.)
    end


    # `Base.show(io::IO, outcome::IsotopeShift.Outcome)`  ... prepares a proper printout of the variable outcome::IsotopeShift.Outcome.
    function Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome) 
        println(io, "level:                   $(outcome.level)  ")
        println(io, "Knms:                    $(outcome.Knms)  ")
        println(io, "Ksms:                    $(outcome.Ksms)  ")
        println(io, "Fme:                     $(outcome.Fme)  ")
        println(io, "Fdensity:                $(outcome.Fdensity)  ")
        println(io, "Xboson:                  $(outcome.Xboson)  ")
        println(io, "amplitudeKnms:           $(outcome.amplitudeKnms)  ")
        println(io, "amplitudeKsmsA:          $(outcome.amplitudeKsmsA)  ")
        println(io, "amplitudeKsmsB:          $(outcome.amplitudeKsmsB)  ")
        println(io, "amplitudeKsmsC:          $(outcome.amplitudeKsmsC)  ")
    end

end # module
