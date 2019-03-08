
"""
`module  JAC.PlasmaShift`  ... a submodel of JAC that contains all methods for computing plasma shifts and related properties between for 
                               some level or between some initial and final-state multiplets; it is using JAC, JAC.ManyElectron.
"""
module PlasmaShift

    using Printf, JAC, JAC.ManyElectron


    """
    `@enum   PlasmaModel`  ... defines a enumeration for the (allowed) plasma models.
    """
    @enum   PlasmaModel    NoPlasmaModel    DebyeHueckel    DebeyBox    IonSphere


    """
    `JAC.PlasmaModel(sa::String)`  ... constructor for a given String.
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


    """
    `Base.show(io::IO, model::PlasmaModel)`  ... prepares a proper printout of the variable model::PlasmaModel.
    """
    function Base.show(io::IO, model::PlasmaModel) 
        print(io, string(model) )
    end


    """
    `Base.string(model::PlasmaModel)`  ... provides a proper printout of the variable model::PlasmaModel.
    """
    function Base.string(model::PlasmaModel) 
        if       model == NoPlasmaModel     return("No plasma model")
        elseif   model == DebyeHueckel      return("Debye-Hueckel model")
        elseif   model == DebeyBox          return("Debey-box model")
        elseif   model == IonSphere         return("Ion-sphere model")
        else     error("stop a")
        end
    end


    """
    `struct  PlasmaShift.Outcome`  ... defines a type to keep the results of a plasma-shift computation for a fine-structure level

        + level                 ::Level            ... Fine-structure levels to which the results refer to.
        + K                         ::Float64            ... K enhancement parameter
    """
    struct Outcome 
        level                   ::Level
        K                       ::Float64      
    end


    """
    `Base.show(io::IO, outcome::PlasmaShift.Outcome)`  ... prepares a proper printout of the variable outcome::PlasmaShift.Outcome.
    """
    function Base.show(io::IO, outcome::PlasmaShift.Outcome) 
        println(io, "level:           $(outcome.level)  ")
    end


    """
    `struct  PlasmaShift.Settings`  ... defines a type for the details and parameters of computing level energies with plasma interactions.

        + plasmaModel            ::PlasmaModel             ... Specified the particular type of plasma model; from {ion-sphere, debye}
        + lambdaDebye            ::Float64                 ... The lambda parameter of different plasma models.
        + ionSphereR0            ::Float64                 ... The effective radius of the ion-sphere model.
        + NoBoundElectrons       ::Int64                   ... Effective number of bound electrons.
        + diagonalizationMethod  ::String                  ... method for diagonalizing the matrix.
        + printBeforeComputation ::Bool                    ... True if a list of selected levels is printed before the actual computations start. 
        + selectLevels           ::Bool                    ... true, if specific level (number)s have been selected for computation.
        + selectedLevels         ::Array{Int64,1}          ... Level number that have been selected.
        + selectSymmetries       ::Bool                    ... true, if specific level symmetries have been selected for computation.
        + selectededSymmetries   ::Array{LevelSymmetry,1}  ... Level symmetries that have been selected.
    """
    struct Settings 
        plasmaModel              ::PlasmaModel
        lambdaDebye              ::Float64 
        ionSphereR0              ::Float64
        NoBoundElectrons         ::Int64
        diagonalizationMethod    ::String
        printBeforeComputation   ::Bool  
        selectLevels             ::Bool
        selectedLevels           ::Array{Int64,1}
    end 


    """
    `JAC.PlasmaShift.Settings()`  ... constructor for a standard instance of PlasmaShift.Settings.
    """
    function Settings()
        Settings(NoPlasmaModel, 0., 0., 0, "", false, false, Int64[])
    end


    """
    `Base.show(io::IO, settings::PlasmaShift.Settings)`  ... prepares a proper printout of the settings::PlasmaShift.Settings.
    """
    function Base.show(io::IO, settings::PlasmaShift.Settings)
        println(io, "plasmaModel:            $(settings.plasmaModel)  ")
        println(io, "lambda:                 $(settings.lambda)  ")
        println(io, "diagonalizationMethod:  $(settings.diagonalizationMethod)  ")
        println(io, "printBeforeComputation: $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:           $(settings.selectLevels)  ")
        println(io, "selectedLevels:         $(settings.selectedLevels)  ")
    end


    """
    `JAC.PlasmaShift.computeMatrixIonSphere(basis::Basis, settings::PlasmaShift.Settings)`  ... to compute the CI matrix for an ion-sphere
                     interaction. A quadratic matrix of dimension length(basis.csfs) x length(basis.csfs) is returned. 
    """
    function computeMatrixIonSphere(basis::Basis, settings::PlasmaShift.Settings)   
        error("Not yet implemented.") 
    end


    """
    `JAC.PlasmaShift.computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                                 settings::PlasmaShift.Settings; output=true)` 
         ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by the given settings. 
         The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
             settings::PlasmaShift.Settings; output=true)
        println("")
        printstyled("JAC.PlasmaShift.computeOutcomes(): The computation of plasma-shifts starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = JAC.PlasmaShift.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    JAC.PlasmaShift.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = PlasmaShift.Outcome[]
        for  outcome in outcomes
            ## newOutcome = JAC.PlasmaShift.computeAmplitudesProperties(outcome, nm, grid, settings) 
            newOutcome = JAC.PlasmaShift.Outcome(outcome.level, 3.0 )
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        JAC.PlasmaShift.displayResults(stdout, newOutcomes)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    JAC.PlasmaShift.displayResults(iostream, newOutcomes)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `JAC.PlasmaShift.determineOutcomes(multiplet::Multiplet, settings::PlasmaShift.Settings)`  
         ... to determine a list of Outcomes's for the computation of the multipole polarizibilities for the given multiplet. It takes 
         into account the particular selections and settings. An Array{PlasmaShift.Outcome,1} is returned. Apart from the 
         level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::PlasmaShift.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        outcomes = PlasmaShift.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, PlasmaShift.Outcome(multiplet.levels[i], 3. ) )
        end
        return( outcomes )
    end


    """
    `JAC.PlasmaShift.displayOutcomes(outcomes::Array{PlasmaShift.Outcome,1})`  ... to display a list of levels that have been selected 
         for the computations. A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{PlasmaShift.Outcome,1})
        println(" ")
        println("  Selected PlasmaShift levels:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(43))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(43)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.level.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", outcome.level.energy)) * "    "
            println( sa )
        end
        println("  ", JAC.TableStrings.hLine(43))
        #
        return( nothing )
    end


    """
    `JAC.PlasmaShift.displayResults(stream::IO, outcomes::Array{PlasmaShift.Outcome,1})`  ... to display the energies, M_ms and F-parameters, etc. for the 
         selected levels. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{PlasmaShift.Outcome,1})
        println(stream, " ")
        println(stream, "  Plasma shift parameters:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(64))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4)              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(14, "xxx";     na=4)              
        sb = sb * JAC.TableStrings.center(14, "    " ; na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(64)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.level.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            energy = 1.0
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", energy))              * "    "
            sa = sa * @sprintf("%.8e", outcome.K)                                               * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(64), "\n\n")
        #
        return( nothing )
    end

end # module
