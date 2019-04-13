
"""
`module  JAC.Semiempirical`  ... a submodel of JAC that contains all methods to set-up and process (simple)  semiempirical estimations of 
                                 atomic properties that cannot be calculated so easily in a rigorious manner; it is using JAC, 
                                 JAC.ManyElectron, JAC.ImpactIonization.
"""
module Semiempirical

    using JAC, JAC.ManyElectron, JAC.ImpactIonization


    """
    `@enum   AtomicCrossSection`  ... defines a enumeration for (supported) cross section estimation.
    """
    @enum   AtomicCrossSection    NoCrossSection    ImpactIonizationCS


    """
    `JAC.AtomicCrossSection(sa::String)`  ... constructor for a given String.
    """
    function AtomicCrossSection(sa::String)
        if       sa == "No cross section"      wa = NoCrossSection
        elseif   sa == "impact-ionization"     wa = ImpactIonizationCS
        else
            error("stop a")
        end

        return( wa )
    end  


    # `Base.show(io::IO, cs::AtomicCrossSection)`  ... prepares a proper printout of the variable cs::AtomicCrossSection.
    function Base.show(io::IO, cs::AtomicCrossSection) 
        print(io, string(cs) )
    end


    # `Base.string(cs::AtomicCrossSection)`  ... provides a proper printout of the variable cs::AtomicCrossSection.
    function Base.string(cs::AtomicCrossSection) 
        if       cs == NoCrossSection       return("No cross section")
        elseif   cs == ImpactIonizationCS   return("impact-ionization cross section")
        else     error("stop a")
        end
    end


    """
    `struct  Estimation`  ... defines a type for defining  (a framework for) semiempirical estimates of different properties that cannot
                              be calculated rigoriously in JAC.

        + calcCrossSections    ::Bool                        ... True, if atomic cross sections of some kind are to be calculated
        + crossSection         ::AtomicCrossSection          ... Selected cross sections type.
        + config               ::Configuration               ... A non-relativistic configurations for which computations are to be carried out.
        + crossSectionSettings ::Union{JAC.ImpactIonization.Settings} ... Provides the settings for the selected cross section estimates.
    """
    struct  Estimation
        calcCrossSections      ::Bool  
        crossSection           ::AtomicCrossSection
        config                 ::Configuration
        crossSectionSettings   ::Union{JAC.ImpactIonization.Settings}
    end 


    """
    `JAC.Semiempirical.Estimation()`  ... constructor for an 'empty' instance::Semiempirical.Estimation.
    """
    function Estimation()
        Estimation(false, NoCrossSection, Configuration(), JAC.ImpactIonization.Settings() )
    end


    # `Base.show(io::IO, estimation::Semiempirical.Estimation)`  ... prepares a proper printout of the variable estimation::Semiempirical.Estimation.
    function Base.show(io::IO, estimation::Semiempirical.Estimation) 
        println(io, "calcCrossSections:        $(estimation.calcCrossSections)  ")
        println(io, "crossSection:             $(estimation.crossSection)  ")
        println(io, "config:                   $(estimation.config)  ")
        println(io, "crossSectionSettings:     $(estimation.crossSectionSettings)  ")
    end


    """
    `JAC.Semiempirical.Estimation("interactive")`  ... constructor to generate a new instance of Semiempirical.Estimation interactively by 
         replying to some detailed dialog.  **Not yet fully implemented !**
    """
    function Semiempirical.Estimation(sa::String)
        sa != "interactive"   &&    error("Unsupported keystring = $sa.")
    end

end # module
