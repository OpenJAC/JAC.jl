
"""
`module  JAC.MultipolePolarizibility`  ... a submodel of JAC that contains all methods for computing dynamic multipole polarizibilities, including 
                                         (dynanmic) electric- & magnetic-dipole polarizibilities and others; it is using JAC, JAC.ManyElectron, 
                                         JAC.Radial.
"""
module MultipolePolarizibility

    using Printf,  JAC, JAC.BasicTypes, JAC.ManyElectron, JAC.Radial, JAC.Nuclear
    global JAC_counter = 0
    

    """
    `struct  MultipolePolarizibility.Amplitude`  ... defines a type for an individual multipole term as defined by the multipole and frequency.

        + multipole  ::EmMultipole         ... Multipole of the dynamic polarizibility contribution.
        + omega      ::Float64             ... (dynamic) frequency of the polarizibility.
        + value      ::Complex{Float64}    ... Value of this multipole amplitude.

    """
    struct  Amplitude
        multipole    ::EmMultipole
        omega        ::Float64
        value        ::Complex{Float64}
    end 


    """
    `JAC.MultipolePolarizibility.Amplitude()`  ... constructor for an `empty` instance of MultipolePolarizibility.Amplitude.
    """
    function Amplitude()
        Amplitude(JAC.E1, 0., 0.)
    end


    # `Base.show(io::IO, amplitude::MultipolePolarizibility.Amplitude)`  ... prepares a proper printout of the variable amplitude.
    function Base.show(io::IO, amplitude::MultipolePolarizibility.Amplitude) 
        println(io, "multipole :       $(amplitude.multipole)  ")
        println(io, "omega:            $(amplitude.omega)  ")
        println(io, "value:            $(amplitude.value)  ")
    end


    """
    `struct  MultipolePolarizibility.Outcome`  ... defines a type to keep the outcome of a dynamic-polarizibility computation.

        + level           ::Level                                     ... Atomic level to which the outcome refers to.
        + amplitudes      ::Array{MultipolePolarizibility.Amplitude,1}  ... multipole (and frequency-dependent) contribution to the polarizibility.
    """
    struct Outcome 
        level             ::Level  
        amplitudes        ::Array{MultipolePolarizibility.Amplitude,1}
    end 


    """
    `JAC.MultipolePolarizibility.Outcome()`  ... constructor for an `empty` instance of MultipolePolarizibility.Outcome for the computation of 
                                               dynamic polarizibilities.
    """
    function Outcome()
        Outcome(Level(), MultipolePolarizibility.Amplitude[])
    end


    # `Base.show(io::IO, outcome::MultipolePolarizibility.Outcome)` 
    #		 ... prepares a proper printout of the variable outcome::MultipolePolarizibility.Outcome.
    function Base.show(io::IO, outcome::MultipolePolarizibility.Outcome) 
        println(io, "level:                     $(outcome.level)  ")
        println(io, "amplitudes:                $(outcome.amplitudes)  ")
    end


    """
    `struct  MultipolePolarizibility.Settings`  ... defines a type for the details and parameters of computing multipolar polarizibilities.

        + multipoles                ::Array{EmMultipole,1}  ... Multipoles to be considered for the polarity.
        + nLower                    ::Int64                 ... Lower and upper indices in the configurations for the summation 
        + nUpper                    ::Int64                     over the intermediate levels.
        + omegas                    ::Array{Float64,1}      ... List of omegas (energies) of the dynamic polarizibility.
        + printBeforeComputation    ::Bool                  ... True if a list of selected levels is printed before the actual computations start. 
        + selectLevels              ::Bool                  ... True if individual levels are selected for the computation.
        + selectedLevels            ::Array{Int64,1}        ... List of selected levels.
    """
    struct Settings
        multipoles                  ::Array{EmMultipole,1}
        nLower                      ::Int64  
        nUpper                      ::Int64
        omegas                      ::Array{Float64,1}
        printBeforeComputation      ::Bool 
        selectLevels                ::Bool
        selectedLevels              ::Array{Int64,1}
    end 


    """
    `JAC.MultipolePolarizibility.Settings()`  ... constructor for an `empty` instance of MultipolePolarizibility.Settings for the computation 
                                                dynamic polarizibilities.
    """
    function Settings()
        Settings(EmMultipole[], 0, 0, Float64[], false, false, Int64[])
    end


    # `Base.show(io::IO, settings:MultipolePolarizibility.Settings)`  
    #		... prepares a proper printout of the variable settings::MultipolePolarizibility.Settings.
    function Base.show(io::IO, settings::MultipolePolarizibility.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "nLower:                   $(settings.nLower)  ")
        println(io, "nUpper:                   $(settings.nUpper)  ")
        println(io, "omegas:                   $(settings.omegas)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    """
    `JAC.MultipolePolarizibility.computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, grid::Radial.Grid, 
                               settings::MultipolePolarizibility.Settings) ... to compute all amplitudes and properties for a given level; 
         an outcome::MultipolePolarizibility.Outcome is returned for which the amplitudes and properties are now evaluated explicitly.
    """
    function  computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, grid::Radial.Grid, settings::MultipolePolarizibility.Settings)
        ## global JAC_counter
        ## return( newOutcome )
    end


    """
    `JAC.MultipolePolarizibility.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                                 settings::MultipolePolarizibility.Settings; output=true)` 
         ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by the given settings. 
         The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
             settings::MultipolePolarizibility.Settings; output=true)
        println("")
        printstyled("JAC.MultipolePolarizibility.computeOutcomes(): The computation of multipole polarizibilities starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = JAC.MultipolePolarizibility.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    JAC.MultipolePolarizibility.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = MultipolePolarizibility.Outcome[]
        for  outcome in outcomes
            ## newOutcome = JAC.MultipolePolarizibility.computeAmplitudesProperties(outcome, nm, grid, settings) 
            newOutcome = JAC.MultipolePolarizibility.Outcome(outcome.level, [ MultipolePolarizibility.Amplitude(JAC.E1, 1., 2.)] )
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        JAC.MultipolePolarizibility.displayResults(stdout, newOutcomes)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    JAC.MultipolePolarizibility.displayResults(iostream, newOutcomes)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `JAC.MultipolePolarizibility.determineOutcomes(multiplet::Multiplet, settings::MultipolePolarizibility.Settings)`  
         ... to determine a list of Outcomes's for the computation of the multipole polarizibilities for the given multiplet. It takes 
         into account the particular selections and settings. An Array{MultipolePolarizibility.Outcome,1} is returned. Apart from the 
         level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::MultipolePolarizibility.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        outcomes = MultipolePolarizibility.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, MultipolePolarizibility.Outcome(multiplet.levels[i], MultipolePolarizibility.Amplitude[]) )
        end
        return( outcomes )
    end


    """
    `JAC.MultipolePolarizibility.displayOutcomes(outcomes::Array{MultipolePolarizibility.Outcome,1})`  ... to display a list of levels that have been selected 
         for the computations. A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{MultipolePolarizibility.Outcome,1})
        println(" ")
        println("  Selected MultipolePolarizibility levels:")
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
    `JAC.MultipolePolarizibility.displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})`  ... to display the energies, M_ms and F-parameters, etc. for the 
         selected levels. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})
        println(stream, " ")
        println(stream, "  Multipole polarizibilities and amplitudes:")
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
            ## sa = sa * @sprintf("%.8e", outcome.K)                                               * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(64), "\n\n")
        #
        return( nothing )
    end

end # module


