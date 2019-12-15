
"""
`module  JAC.PlasmaShift`  
    ... a submodel of JAC that contains all methods for computing plasma shifts and related properties between for some level or 
        between some initial and final-state multiplets.
"""
module PlasmaShift

    using Printf, ..Basics, ..Defaults, ..Nuclear, ..ManyElectron, ..Radial, ..TableStrings

    #============================
    """
    `@enum   PlasmaModel`  ... defines a enumeration for the (allowed) plasma models.

        + NoPlasmaModel      ... No plasma model defined.
        + DebyeHueckel       ... Debye-Hueckel plasma model.
        + IonSphere          ... Ion-sphere (not yet supported).
    """
    @enum   PlasmaModel    NoPlasmaModel    DebyeHueckel    DebeyBox    IonSphere


    """
    `PlasmaModel(sa::String)`  ... constructor for a given String.
    """
    function PlasmaModel(sa::String)
        if       sa == "No Model"              wa = NoPlasmaModel
        elseif   sa == "Debye-Hueckel"         wa = DebyeHueckel
        elseif   sa == "DebeyBox"              wa = DebeyBox
        elseif   sa == "IonSphere"             wa = IonSphere
        else     error("stop a")
        end

        return( UseGauge(wa) )
    end  


    # `Base.show(io::IO, model::PlasmaModel)`  ... prepares a proper printout of the variable model::PlasmaModel.
    function Base.show(io::IO, model::PlasmaModel) 
        print(io, string(model) )
    end


    # `Base.string(model::PlasmaModel)`  ... provides a proper printout of the variable model::PlasmaModel.
    function Base.string(model::PlasmaModel) 
        if       model == NoPlasmaModel     return("No plasma model")
        elseif   model == DebyeHueckel      return("Debye-Hueckel model")
        elseif   model == DebeyBox          return("Debey-box model")
        elseif   model == IonSphere         return("Ion-sphere model")
        else     error("stop a")
        end
    end          ==========================================#

    
    """
    `abstract type PlasmaShift.AbstractPlasmaModel` 
        ... defines an abstract and a number of singleton types for the the (allowed) plasma models.

        + NoPlasmaModel      ... No plasma model defined.
        + DebyeHueckel       ... Debye-Hueckel plasma model.
        + IonSphere          ... Ion-sphere (not yet supported).
    """
    abstract type  AbstractPlasmaModel                      end
    struct         NoPlasmaModel  <:  AbstractPlasmaModel   end
    struct         DebyeHueckel   <:  AbstractPlasmaModel   end
    struct         IonSphere      <:  AbstractPlasmaModel   end

    
    # `Base.show(io::IO, model::AbstractPlasmaModel)`  ... prepares a proper printout of the variable model::AbstractPlasmaModel.
    function Base.show(io::IO, model::AbstractPlasmaModel) 
        print(io, string(model) )
    end

    
    # `Base.string(model::AbstractPlasmaModel)`  ... provides a proper printout of the variable model::AbstractPlasmaModel.
    function Base.string(model::AbstractPlasmaModel) 
        if       model == NoPlasmaModel()     return("No plasma model")
        elseif   model == DebyeHueckel()      return("Debye-Hueckel model")
        elseif   model == DebeyBox()          return("Debey-box model")
        elseif   model == IonSphere()         return("Ion-sphere model")
        else     error("stop a")
        end
    end 

    
    """
    `struct  PlasmaShift.Settings`  ... defines a type for the details and parameters of computing level energies with plasma interactions.

        + plasmaModel      ::AbstractPlasmaModel        ... Specify a particular plasma model, e.g. ion-sphere, debye.
        + lambdaDebye      ::Float64                    ... The lambda parameter of different plasma models.
        + ionSphereR0      ::Float64                    ... The effective radius of the ion-sphere model.
        + NoBoundElectrons ::Int64                      ... Effective number of bound electrons.
    """
    struct Settings 
        plasmaModel        ::AbstractPlasmaModel
        lambdaDebye        ::Float64 
        ionSphereR0        ::Float64
        NoBoundElectrons   ::Int64
    end 


    """
    `PlasmaShift.Settings()`  ... constructor for a standard instance of PlasmaShift.Settings.
    """
    function Settings()
        Settings(DebyeHueckel(), 0.25, 0., 0)
    end


    # `Base.show(io::IO, settings::PlasmaShift.Settings)`  ... prepares a proper printout of the settings::PlasmaShift.Settings.
    function Base.show(io::IO, settings::PlasmaShift.Settings)
        println(io, "plasmaModel:            $(settings.plasmaModel)  ")
        println(io, "lambdaDebye:            $(settings.lambdaDebye)  ")
        println(io, "ionSphereR0:            $(settings.ionSphereR0)  ")
        println(io, "NoBoundElectrons:       $(settings.NoBoundElectrons)  ")
    end


    """
    `struct  PlasmaShift.AugerSettings`  
        ... defines a type for the details and parameters of computing Auger rates with plasma interactions.

        + plasmaModel            ::AbstractPlasmaModel           ... Specify a particular plasma model, e.g. ion-sphere, debye.
        + lambdaDebye            ::Float64                       ... The lambda parameter of different plasma models.
        + ionSphereR0            ::Float64                       ... The effective radius of the ion-sphere model.
        + NoBoundElectrons       ::Int64                         ... Effective number of bound electrons.
        + printBefore ::Bool                          ... True, if all energies and lines are printed before their evaluation.
        + selectLines            ::Bool                          ... True, if lines are selected individually for the computations.
        + selectedLines          ::Array{Tuple{Int64,Int64},1}   ... List of lines, given by tupels (inital-level, final-level).
    """
    struct AugerSettings 
        plasmaModel              ::AbstractPlasmaModel
        lambdaDebye              ::Float64 
        ionSphereR0              ::Float64
        NoBoundElectrons         ::Int64
        printBefore   ::Bool 
        selectLines              ::Bool
        selectedLines            ::Array{Tuple{Int64,Int64},1}
    end 


    """
    `PlasmaShift.AugerSettings()`  ... constructor for a standard instance of PlasmaShift.AugerSettings.
    """
    function AugerSettings()
        AugerSettings(DebyeHueckel(), 0.25, 0., 0, true, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::PlasmaShift.AugerSettings)`  ... prepares a proper printout of the settings::PlasmaShift.AugerSettings.
    function Base.show(io::IO, settings::PlasmaShift.AugerSettings)
        println(io, "plasmaModel:             $(settings.plasmaModel)  ")
        println(io, "lambdaDebye:             $(settings.lambdaDebye)  ")
        println(io, "ionSphereR0:             $(settings.ionSphereR0)  ")
        println(io, "NoBoundElectrons:        $(settings.NoBoundElectrons)  ")
        println(io, "printBefore:  $(settings.printBefore)  ")
        println(io, "selectLines:             $(settings.selectLines)  ")
        println(io, "selectedLines:           $(settings.selectedLines)  ")
    end


    """
    `struct  PlasmaShift.PhotoSettings`  
        ... defines a type for the details and parameters of computing photoionization rates with plasma interactions.

        + plasmaModel            ::AbstractPlasmaModel           ... Specify a particular plasma model, e.g. ion-sphere, debye.
        + lambdaDebye            ::Float64                       ... The lambda parameter of different plasma models.
        + ionSphereR0            ::Float64                       ... The effective radius of the ion-sphere model.
        + NoBoundElectrons       ::Int64                         ... Effective number of bound electrons.
        + multipoles             ::Array{Basics.EmMultipole}     ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                 ::Array{Basics.UseGauge}        ... Specifies the gauges to be included into the computations.
        + photonEnergies         ::Array{Float64,1}              ... List of photon energies.  
        + printBefore ::Bool                          ... True, if all energies and lines are printed before their evaluation.
        + selectLines            ::Bool                          ... True, if lines are selected individually for the computations.
        + selectedLines          ::Array{Tuple{Int64,Int64},1}   ... List of lines, given by tupels (inital-level, final-level).
    """
    struct PhotoSettings 
        plasmaModel              ::AbstractPlasmaModel
        lambdaDebye              ::Float64 
        ionSphereR0              ::Float64
        NoBoundElectrons         ::Int64
        multipoles               ::Array{Basics.EmMultipole}
        gauges                   ::Array{Basics.UseGauge}  
        photonEnergies           ::Array{Float64,1} 
        printBefore   ::Bool 
        selectLines              ::Bool
        selectedLines            ::Array{Tuple{Int64,Int64},1}
    end 


    """
    `PlasmaShift.PhotoSettings()`  ... constructor for a standard instance of PlasmaShift.PhotoSettings.
    """
    function PhotoSettings()
        PhotoSettings(DebyeHueckel(), 0.25, 0., 0, [E1], [Basics.UseCoulomb], Float64[], true, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::PlasmaShift.PhotoSettings)`  ... prepares a proper printout of the settings::PlasmaShift.PhotoSettings.
    function Base.show(io::IO, settings::PlasmaShift.PhotoSettings)
        println(io, "plasmaModel:             $(settings.plasmaModel)  ")
        println(io, "lambdaDebye:             $(settings.lambdaDebye)  ")
        println(io, "ionSphereR0:             $(settings.ionSphereR0)  ")
        println(io, "NoBoundElectrons:        $(settings.NoBoundElectrons)  ")
        println(io, "multipoles:              $(settings.multipoles)  ")
        println(io, "gauges:                  $(settings.gauges)  ")
        println(io, "photonEnergies:          $(settings.photonEnergies)  ")
        println(io, "printBefore:  $(settings.printBefore)  ")
        println(io, "selectLines:             $(settings.selectLines)  ")
        println(io, "selectedLines:           $(settings.selectedLines)  ")
    end


    """
    `PlasmaShift.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                     asfSettings::AsfSettings, settings::PlasmaShift.Settings; output=true)` 
        ... to compute a new multiplet from the plasma-modified Hamiltonian matrix as specified by the given 
            settings. The diagonalization of the (plasma-modified) Hamiltonian matrix follows the asfSettings that
            were applied before for the computation of the multiplet; warnings are issued if features were selected
            that are not supported for plasma calculations, such as the Breit interaction, QED shifts and others.
            The energies and plasma shifts (with regard to the unperturbed multiplet) are printed to screen 
            and into the summary file. The generated (plasma) multiplet is returned if output=true and nothing otherwise,
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                             asfSettings::AsfSettings, settings::PlasmaShift.Settings; output=true)
        println("")
        printstyled("PlasmaShift.computeOutcomes(): The computation of plasma-shifts starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------- \n", color=:light_green)
        #
        basis        = multiplet.levels[1].basis
        newMultiplet = Basics.perform("computation: CI for plasma", basis, nm, grid, asfSettings, settings)
        # Print all results to screen
        PlasmaShift.displayResults(stdout, multiplet, newMultiplet, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    PlasmaShift.displayResults(iostream, multiplet, newMultiplet, settings)   end
        #
        if    output    return( newMultiplet )
        else            return( nothing )
        end
    end


    """
    `PlasmaShift.displayResults(stream::IO, multiplet::Multiplet, newMultiplet::Multiplet, settings::PlasmaShift.Settings)`  
        ... to display the energies, M_ms and F-parameters, etc. for the  selected levels. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, multiplet::Multiplet, pMultiplet::Multiplet, settings::PlasmaShift.Settings)
        println(stream, " ")
        println(stream, " ")
        println(stream, "  Plasma shifts for $(settings.plasmaModel):")
        
        # Write out the relevant plasma parameters
        if settings.plasmaModel == PlasmaShift.DebyeHueckel()
            println(stream, "\n     + lambda = $(settings.lambdaDebye)")
            println(stream,   "     + Plasma screening included perturbatively in CI matrix but not in SCF field.")
        else  error("Unsupported plasma model = $(plasmaSettings.plasmaModel)")
        end
        
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(64))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(18, "Energy w/o plasma"; na=4)              
        sb = sb * TableStrings.center(18, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "Delta E";     na=4)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(64)) 
        #  
        for  i  in  1:length(multiplet.levels)
            sa  = "  ";    sym = LevelSymmetry( multiplet.levels[i].J, multiplet.levels[i].parity);    
                        newsym = LevelSymmetry( pMultiplet.levels[i].J, pMultiplet.levels[i].parity)
            sa = sa * TableStrings.center(10, TableStrings.level(multiplet.levels[i].index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=6)
            energy = multiplet.levels[i].energy
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))              * "     "
            if  sym == newsym   deltaE = pMultiplet.levels[i].energy - energy;    
                                sa     = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", deltaE))   
            else                sa     = sa * "Level crossing."
            end
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(64), "\n\n")
        #
        return( nothing )
    end

end # module
