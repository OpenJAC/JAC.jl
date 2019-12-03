
"""
`module  JAC.Semiempirical`  ... a submodel of JAC that contains all methods to set-up and process (simple)  semiempirical estimations of 
                                 atomic properties that cannot be calculated so easily in a rigorious manner.
"""
module Semiempirical

    using ..Basics, ..Defaults, ..ManyElectron, ..PeriodicTable    ## ..ImpactIonization


    """
    `@enum   AtomicCrossSection`  ... defines a enumeration for (supported) cross section estimation.
    """
    @enum   AtomicCrossSection    NoCrossSection    ImpactIonizationCS


    """
    `Semiempirical.AtomicCrossSection(sa::String)`  ... constructor for a given String.
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
    ## function Base.show(io::IO, cs::AtomicCrossSection) 
    ## print(io, string(cs) )
    ## end


    # `Base.string(cs::AtomicCrossSection)`  ... provides a proper printout of the variable cs::AtomicCrossSection.
    function Base.string(cs::AtomicCrossSection) 
        if       cs == NoCrossSection       return("No cross section")
        elseif   cs == ImpactIonizationCS   return("impact-ionization cross section")
        else     error("stop a")
        end
    end

    #=========================================
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
    =====================================================#
    

    """
    `Semiempirical.estimate("ionization potential: inner-shell", shell::Shell, Z::Int64)`  
        ... to estimate the ionization potential = mean binding energy of an electron in the given shell; an energy value::Float64 
            is returned. These ionization potentials are taken from the tabulation ....
    """
    function  estimate(sa::String, shell::Shell, Z::Int64)
        if      sa == "ionization potential: inner-shell"
            shellIndex = Dict( Subshell("1s_1/2") => 1, Subshell("2s_1/2") => 2, Subshell("2p_1/2") => 3, Subshell("2p_3/2") => 4, 
                            Subshell("3s_1/2") => 5, Subshell("3p_1/2") => 6, Subshell("3p_3/2") => 7, 
                            Subshell("3d_3/2") => 8, Subshell("3d_3/2") => 9  )
            idx = shellIndex[shell]
            return( PeriodicTable.store("binding energies: Williams (2000)", Z)[idx] )
        else
            error("Unsupported keystring = $sa")
        end
    end



    """
    `Semiempirical.estimate("binding energy", Z::Float64, sh::Subshell)`   
        ... to provide the binding energy of a subshell electron, taken from a semi-empirical tabulations by Williams et al., 
            https://userweb.jlab.org/~gwyn/ebindene.html. A energy::Float64 in  Hartree is returned.
    """
    function estimate(sa::String, Z::Int64, sh::Subshell)
        if     sa == "binding energy"
            wa = PeriodicTable.store_Williams2000(Z)
            if      sh == Subshell("1s_1/2")    wb = wa[1]
            elseif  sh == Subshell("2s_1/2")    wb = wa[2]
            elseif  sh == Subshell("2p_1/2")    wb = wa[3]
            elseif  sh == Subshell("2p_3/2")    wb = wa[4]
            elseif  sh == Subshell("3s_1/2")    wb = wa[5]
            elseif  sh == Subshell("3p_1/2")    wb = wa[6]
            elseif  sh == Subshell("3p_3/2")    wb = wa[7]
            elseif  sh == Subshell("3d_3/2")    wb = wa[8]
            elseif  sh == Subshell("3d_5/2")    wb = wa[9]
            elseif  sh == Subshell("4s_1/2")    wb = wa[10]
            elseif  sh == Subshell("4p_1/2")    wb = wa[11]
            elseif  sh == Subshell("4p_3/2")    wb = wa[12]
            else    error("No binding energy available for Z = $Z and subshell $sh ")
            end
        else   error("Unsupported keystring")
        end
        #
        if     wb == -1.   error("No binding energy available for Z = $Z and subshell $sh ")
        else   wb = Defaults.convertUnits("energy: from eV to atomic", wb)
        end
        # 
        return( wb )
    end



    """
    `Semiempirical.estimate("binding energy", Z::Float64, conf::Configuration)`  
        ... to provide an approximate binding energy of a given electron configuration. This estimate adds the binding 
            energies of all subshell, taken frogm a semi-empirical tabulations by Williams et al., 
            https://userweb.jlab.org/~gwyn/ebindene.html. If no binding energy is provided by this table, it simply scales the 
            binding energies (e/ 12.) from the next lower shell of the same symmetry
            No relaxation effects are included if several hole states
            occur with regard to the neutral atom. An energy::Float64 in  Hartree is returned.
    """
    function estimate(sa::String, Z::Int64, conf::Configuration)
        if     sa == "binding energy"
            wa = PeriodicTable.bindingEnergies_Williams2000(Z);    wb = 0.
            for (sh,v) in  conf.shells
                if      sh == Shell("1s")    if wa[1]  == -1.   error("stop aa")          else    wb = wb + v * wa[1]    end
                elseif  sh == Shell("2s")    if wa[2]  == -1.   error("stop ab")          else    wb = wb + v * wa[2]    end
                elseif  sh == Shell("2p")    if wa[3]  == -1.   error("stop ac")          else    wb = wb + v * wa[3]    end
                elseif  sh == Shell("3s")    if wa[5]  == -1.   wb = wb + v * wa[2]/12.   else    wb = wb + v * wa[5]    end
                                             if wa[2]  == -1.   error("stop ad")   end
                elseif  sh == Shell("3p")    if wa[6]  == -1.   wb = wb + v * wa[3]/12.   else    wb = wb + v * wa[6]    end
                                             if wa[3]  == -1.   error("stop af")   end
                elseif  sh == Shell("3d")    if wa[8]  == -1.   error("stop ah")          else    wb = wb + v * wa[8]    end
                elseif  sh == Shell("4s")    if wa[10] == -1.   error("stop aj")          else    wb = wb + v * wa[10]   end
                elseif  sh == Shell("4p")    if wa[11] == -1.   error("stop ak")          else    wb = wb + v * wa[11]   end
                else    error("No binding energy available for Z = $Z and subshell $sh ")
                end
            end
        else   error("Unsupported keystring")
        end
        #
        wb = Defaults.convertUnits("energy: from eV to atomic", wb)
        # 
        return( wb )
    end

    
end # module
