
"""
`module  JAC.StarkShift`  
... a submodel of JAC that contains all methods for computing Lande factors and Zeeman properties for some level(s).
"""
module StarkShift


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..MultipoleMoment, 
              ..Nuclear, ..Radial, ..SpinAngular, ..TableStrings


"""
`struct  StarkShift.SublevelJ`  ... defines a type to specify a Stark-shifted sublevel with well-defined J.

    + M                      ::AngularM64        ... M_J-value
    + energy                 ::Float64           ... energy of this sublevel
"""
struct SublevelJ 
    M                        ::AngularM64
    energy                   ::Float64
end 


# `Base.show(io::IO, Jsublevel::StarkShift.SublevelJ)`  ... prepares a proper printout of the variable Jsublevel::StarkShift.SublevelJ.
function Base.show(io::IO, Jsublevel::StarkShift.SublevelJ) 
    println(io, "Sublevel [M=$(Jsublevel.M); energy = $(Jsublevel.energy)]")
end


"""
`struct  StarkShift.Outcome`  
    ... defines a type to keep the Stark shift and electric-quadrupole-moment parameters of a fine-structure level.

    + Jlevel                 ::Level             ... Fine-structure levels to which the results refer to.
    + Theta1                 ::Float64           ... Electric-dipole moment of the atom in this level. 
    + Theta2                 ::Float64           ... Electric-quadrupole moment of the atom in this level. 
    + Jsublevels             ::Array{StarkShift.SublevelJ,1}   
        ... List of the Stark-shifted fine-structure sublevels and data with well-defined J-value 
"""
struct Outcome 
    Jlevel                   ::Level
    Theta1                   ::Float64
    Theta2                   ::Float64
    Jsublevels               ::Array{StarkShift.SublevelJ,1}   
end 


"""
`StarkShift.Outcome()`  ... constructor for an `empty` instance of LStarkShift.Outcome.
"""
function Outcome()
    Outcome(Level(), 0., 0., StarkShiftJsublevels[])
end


# `Base.show(io::IO, outcome::StarkShift.Outcome)`  ... prepares a proper printout of the variable StarkShift.Outcome.
function Base.show(io::IO, outcome::StarkShift.Outcome) 
    println(io, "Jlevel:           $(outcome.Jlevel)  ")
    println(io, "Theta1:           $(outcome.Theta1)  ")
    println(io, "Theta2:           $(outcome.Theta2)  ")
    println(io, "Jsublevels:       $(outcome.Jsublevels)  ")
end


"""
`struct  StarkShift.Settings  <:  AbstractPropertySettings`  
    ... defines a type for the details and parameters of computing the Stark shift and electric-quadrupole-moments
        of fine-structure levels.

    + calcEDM                ::Bool              ... True if the EDM of selected levels is to be calculated.
    + calcEQM                ::Bool              ... True if the EQM of selected levels is to be calculated.
    + calcStarkshifts        ::Bool              ... True if the energy shifts are to be printed (not yet).
    + printBefore            ::Bool              ... True if a list of selected levels is printed before the actual computations start. 
    + EField                 ::Float64           ... Strength of the electric field in [V/cm]
    + levelSelection         ::LevelSelection    ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings 
    calcEDM                  ::Bool
    calcEQM                  ::Bool
    calcStarkshifts          ::Bool
    printBefore              ::Bool 
    EField                   ::Float64 
    levelSelection           ::LevelSelection
end 


"""
`StarkShift.Settings()`  
    ... constructor for an `empty` instance of StarkShift.Settings for the computation of ... .
"""
function Settings()
        Settings(false, false, false, false, 0., LevelSelection() )
end


# `Base.show(io::IO, settings::StarkShift.Settings)`  ... prepares a proper printout of the variable settings::StarkShift.Settings.
function Base.show(io::IO, settings::StarkShift.Settings) 
    println(io, "calcEDM:               $(settings.calcEDM)  ")
    println(io, "calcEQM:               $(settings.calcEQM)  ")
    println(io, "calcStarkshifts:       $(settings.calcStarkshifts)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "EField:                $(settings.EField)  ")
    println(io, "levelSelection:        $(settings.levelSelection)  ")
end


"""
`StarkShift.computeAmplitudesProperties(outcome::StarkShift.Outcome, grid::Radial.Grid, settings::StarkShift.Settings)`  
    ... to compute all amplitudes and properties of for a given level. An outcome::StarkShift.Outcome is returned for 
        which the amplitudes and all requested properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::StarkShift.Outcome, grid::Radial.Grid, settings::StarkShift.Settings)
    edm = eqm = 0.
    if  settings.calcEDM    edm = StarkShift.computeTheta1(outcome.Jlevel, grid)    end
    if  settings.calcEQM    eqm = StarkShift.computeTheta2(outcome.Jlevel, grid)    end
        
    if  settings.calcStarkshifts
        println(">>>> Warning: Stark shifts of level energies not yet implemented.")
    end
    
    newOutcome = StarkShift.Outcome( outcome.Jlevel, edm, eqm, outcome.Jsublevels)
    return( newOutcome )
end


"""
`StarkShift.computeTheta1(level::Level, grid::Radial.Grid)`  
    ... to compute the electric-dipole-moment [Theta^(1)] for the given level; an edm::Float64 is returned.
"""
function  computeTheta1(level::Level, grid::Radial.Grid)
    print(">>> Calculate the Theta^(1) (shift) coefficient for level (index=$(level.index), J=$(level.J)) ...")
    M = AngularM64(level.J.num, level.J.den)
    Theta1 = MultipoleMoment.emmStaticAmplitude(1, level, level, grid; display=false)
             ## AngularMomentum.ClebschGordan(level.J, M, AngularJ64(1), AngularM64(0), level.J, M) / sqrt(Basics.twice(level.J) + 1)
    J      = level.J.num/level.J.den
    w3j    = AngularMomentum.Wigner_3j(J, 2, J, J, 0, -J)
    Theta1 = Theta1 * (-1)^(2*J - 1 + 1) * w3j * sqrt(2*J+1)
         
    Theta1JG = MultipoleMoment.amplitude(1, level, grid, display=false)
    println("   Theta^(1) (level=$(level.index) [J=$(level.J)$(string(level.parity))]) = $(Theta1)  $(Theta1JG)  e a_o")
    return( Theta1 )
end



"""
`StarkShift.computeTheta2(level::Level, grid::Radial.Grid)`  
    ... to compute the electric-quadrupole-moment [Theta^(2)] for the given level; an eqm::Float64 is returned.
        Until July 2024, it remains unclear why the Clebsch-Gordan coefficient does not give a proper result
"""
function  computeTheta2(level::Level, grid::Radial.Grid)
    print(">>> Calculate the Theta^(2) (shift) coefficient for level (index=$(level.index), J=$(level.J)) ...")
    M = AngularM64(level.J.num, level.J.den)
    Theta2 = MultipoleMoment.emmStaticAmplitude(2, level, level, grid; display=false) 
             ## AngularMomentum.ClebschGordan(level.J, M, AngularJ64(2), AngularM64(0), level.J, M) / sqrt(Basics.twice(level.J) + 1)
    # Use the conversion factor by Jan Gilles (2024)
    J      = level.J.num/level.J.den
    w3j    = AngularMomentum.Wigner_3j(J, 2, J, J, 0, -J)
    Theta2 = Theta2 * (-1)^(2*J - 2 + 1) * w3j * sqrt(2*J+1)

    Theta2JG = MultipoleMoment.amplitude(2, level, grid, display=false)
    println("   Theta^(2) (level=$(level.index) [J=$(level.J)$(string(level.parity))]) = $(Theta2)  $(Theta2JG)  e a_o^2")
    return( Theta2 )
end


"""
`StarkShift.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkShift.Settings; output=true)`  
    ... to compute (as selected) the EDM and EQM factors for the levels of the given multiplet and as specified by the given 
        settings. The results are returned in table (if required) or nothing is returned otherwise.
        The nuclear model is not used at present.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkShift.Settings; output=true)
    println("")
    printstyled("StarkShift.computeOutcomes(): The computation of the electric multipole moments starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    outcomes = StarkShift.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    StarkShift.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = StarkShift.Outcome[]
    for  outcome in outcomes
        newOutcome = StarkShift.computeAmplitudesProperties(outcome, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    StarkShift.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    StarkShift.displayResults(iostream, newOutcomes, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`StarkShift.determineOutcomes(multiplet::Multiplet, settings::StarkShift.Settings)`  
    ... to determine a list of Outcomes's for the computation of electric-multipole-moments and energy shifts for 
        levels from the given multiplet. It takes into account the particular selections and settings. 
        An Array{StarkShift.Outcome,1} is returned. Apart from the level specification, all physical properties are 
        still set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::StarkShift.Settings) 
    # Define values that depend on the requested computations
    outcomes   = StarkShift.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            Jsublevels = StarkShift.SublevelJ[];   Mvalues = AngularMomentum.m_values(level.J) 
            for  M in Mvalues   push!(Jsublevels, StarkShift.SublevelJ(M, 0.) )    end
            push!( outcomes, StarkShift.Outcome(level, 0., 0., Jsublevels) )
        end
    end
    return( outcomes )
end


"""
`StarkShift.displayOutcomes(outcomes::Array{StarkShift.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected 
        levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{StarkShift.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected levels for EMM and Stark-shift computations:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.Jlevel.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`StarkShift.displayResults(stream::IO, outcomes::Array{StarkShift.Outcome,1}, settings::StarkShift.Settings)`  
    ... to display the energies, Lande factors, Zeeman amplitudes etc. for the selected levels. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{StarkShift.Outcome,1}, settings::StarkShift.Settings)
    #
    if  settings.calcEDM  ||  settings.calcEQM
        nx = 82
        println(stream, " ")
        println(stream, "  Atomic electric-multipole-moment (EMM) shift coefficients:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(14, "Theta^(1)/EDM"; na=5)                          
        sb = sb * TableStrings.center(14, "[e a_o]"; na=5)        
        sa = sa * TableStrings.center(14, "Theta^(2)/EQM"; na=5)                          
        sb = sb * TableStrings.center(14, "[e a_o^2]"; na=5)        
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = outcome.Jlevel.energy
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.Theta1) )              * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.Theta2) )              * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    if  settings.calcStarkshifts
        nx = 80
        println(stream, " ")
        println(stream, "  Stark-shifts (energy-shifts) of the atomic sublevels:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #
        sa = "      *** Warning:  Not yet implemented.  ***"
        println(stream, sa);    println(stream, sa);    
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end

end # module
