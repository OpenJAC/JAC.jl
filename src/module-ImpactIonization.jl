
"""
`module  JAC.ImpactIonization`  
    ... a submodel of JAC that contains all methods for computing electron impact ionization cross sections.
"""
module ImpactIonization

    using ..Basics, ..ManyElectron, ..Nuclear
    

    """
    `abstract type ImpactIonization.AbstractApproximation` 
        ... defines an abstract and a number of singleton types for the computation of empirical EII cross sections in different 
            approximations.

      + struct BEBapproximation        
        ... to calculate the EII cross sections in the binary-encounter-Bethe (BEB) approximation.

      + struct BEDapproximation        
        ... to calculate the EII cross sections in the binary-encounter-... (BED) approximation.

      + struct MottApproximation        
        ... to calculate the EII cross sections in the Mott approximation.
    """
    abstract type  AbstractApproximation                                      end
    struct         BEBapproximation               <:  AbstractApproximation   end
    struct         BEDapproximation               <:  AbstractApproximation   end
    struct         RelativisticBEBapproximation   <:  AbstractApproximation   end
    struct         RelativisticBEDapproximation   <:  AbstractApproximation   end
    
    export BEBapproximation, BEDapproximation, RelativisticBEBapproximation, RelativisticBEDapproximation


    """
    `struct  ImpactIonization.Settings  <:  AbstractEmpiricalSettings`  
        ... defines a type for the settings of empirical electron-impact ionization cross sections.

        + approximation          ::ImpactIonization.AbstractApproximation ... approximation for computing empirical EII cross sections.
        + impactEnergies         ::Array{Float64,1}                       ... List of impact-energies of the incoming elecgtrons.
        + finalElectronEnergies  ::Array{Float64,1}                       ... List of final-electron energies.
        + calcDifferentialCs     ::Bool               ... True if energy-differential cross sections are to be calculated, and false otherwise.
        + calcPartialCs          ::Bool               ... True if partial (shell-dependent) cross sections are to be calculated, and false otherwise.
        + calcTotalCs            ::Bool               ... True if total cross sections are to be calculated, and false otherwise.
        + shellSelection         ::ShellSelection    ... Describe the selected shells, if any.
    """
    struct Settings   <:  AbstractEmpiricalSettings
        approximation            ::ImpactIonization.AbstractApproximation
        impactEnergies           ::Array{Float64,1}
        finalElectronEnergies    ::Array{Float64,1} 
        calcDifferentialCs       ::Bool 
        calcPartialCs            ::Bool 
        calcTotalCs              ::Bool 
        shellSelection           ::ShellSelection
    end 


    """
    `JAC.ImpactIonization.Settings()`  ... constructor for the default values of empirical electron-impact ionization cross sections.
    """
    function Settings()
       Settings( BEBapproximation(), Float64[], Float64[], false, false, false, ShellSelection() )
    end


    ## `Base.show(io::IO, settings::ImpactIonization.Settings)`  ... prepares a proper printout of the variable settings::ImpactIonization.Settings.
    function Base.show(io::IO, settings::ImpactIonization.Settings) 
        println(io, "approximation:           $(settings.approximation)  ")
        println(io, "impactEnergies:          $(settings.impactEnergies)  ")
        println(io, "finalElectronEnergies:   $(settings.finalElectronEnergies)  ")
        println(io, "calcDifferentialCs:      $(settings.calcDifferentialCs)  ")
        println(io, "calcPartialCs:           $(settings.calcPartialCs)  ")
        println(io, "calcTotalCs:             $(settings.calcTotalCs)  ")
        println(io, "shellSelection:         $(settings.shellSelection)  ")
    end


    """
    `struct  ImpactIonization.CrossSection`  
        ... defines a type for dealing with (shell-dependent) partial EII cross sections that are based on empirical computations.

        + shell             ::Shell             ... shell for which the (differential and partial) cross sections are calculated.
        + electronEnergies  ::Array{Float64,1}  ... energies of the ionized electron.
        + differentialCS    ::Array{Float64,1}  ... differential cross sections for this shell.
        + partialCS         ::Float64           ... differential cross sections for this shell.
    """
    struct  CrossSection
        shell               ::Shell
        electronEnergies    ::Array{Float64,1}
        differentialCS      ::Array{Float64,1}
        partialCS           ::Float64
    end 


    # `Base.show(io::IO, line::ImpactIonization.CrossSection)`  ... prepares a proper printout of the variable line::ImpactIonization.CrossSection.
    function Base.show(io::IO, cs::ImpactIonization.CrossSection) 
        println(io, "shell:              $(cs.shell)  ")
        println(io, "electronEnergies:   $(cs.electronEnergies)  ")
        println(io, "differentialCS:     $(cs.differentialCS)  ")
        println(io, "partialCS:          $(cs.partialCS)  ")
    end


    """
    `ImpactIonization.computeCrossSections(basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)`  
        ... to compute the EII cross sections for all requested shells of the given multiplett; 
            a list of crossSections::Array{ImpactIonization.CrossSection,1} is returned.
    """
    function  computeCrossSections(basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)
        # ... follow the PhotoIonization module
        # determine cross sections
        # compute cross sections
        # display 
        # return results
        
        println("\n\nEnter computeCrossSections()")
        
        ImpactIonization.determineCrossSections(basis, settings)
        # ImpactIonization.displayCrossSections()
    
        println("\n\nNow do it, please !!")
        
        return( false )
    end


    """
    `ImpactIonization.determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)`  
        ... to determine an array of EII cross sections for all requested shells of the given multiplett that need to be calculated in these
            computations; an Array{ImpactIonization.CrossSection,1} is returned.
    """
    function  determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)

        # ... follow the PhotoIonization module ... make some nice tables about the partial and total cross sections
        # ... use shellSelection in settings
       
        println("Enter determineCrossSections()")
        
        return( false )
    end


    """
    `ImpactIonization.displayCrossSections(cs::Array{ImpactIonization.CrossSection,1}, settings::ImpactIonization.Settings)`  
        ... displays the EII cross sections in a neat tabular form; nothing is returned.
    """
    function  displayCrossSections(cs::Array{ImpactIonization.CrossSection,1}, settings::ImpactIonization.Settings)
    
        # ... follow the PhotoIonization module ... make some nice tables about the partial and total cross sections
       
        println("Enter displayCrossSections()")
    
        return( nothing )
    end

end # module
