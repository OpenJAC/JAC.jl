
"""
`module  JAC.TwoElectronOnePhoton`  
    ... a submodel of JAC that contains all methods for computing two-electron-one-photon (TEOP) transition rates between 
        some initial and final-state multiplets. 
"""
module TwoElectronOnePhoton

    using Printf, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial, 
                  ..SpinAngular, ..TableStrings


    """
    `struct  TwoElectronOnePhoton.Settings  <:  AbstractProcessSettings`  
        ... defines a type for the details and parameters of computing radiative lines.

        + multipoles              ::Array{EmMultipoles}     ... Specifies the (radiat. field) multipoles to be included.
        + gauges                  ::Array{UseGauge}         ... Gauges to be included into the computations.
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before comput.
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
        + photonEnergyShift       ::Float64                 ... An overall energy shift for all photon energies.
        + green                   ::Array{AtomicState.GreenChannel,1} ... Green channel of different symmetry.
    """
    struct Settings  <:  AbstractProcessSettings 
        multipoles                ::Array{EmMultipole,1}
        gauges                    ::Array{UseGauge}
        printBefore               ::Bool 
        lineSelection             ::LineSelection 
        photonEnergyShift         ::Float64
        green                     ::Array{AtomicState.GreenChannel,1}
    end 


    """
    `TwoElectronOnePhoton.Settings()`  ... constructor for the default values of TEOP line computations
    """
    function Settings()
        Settings(EmMultipole[E1], UseGauge[Basics.UseCoulomb], false, LineSelection(), 0., AtomicState.GreenChannel[])
    end


    """
    `TwoElectronOnePhoton.Settings(set::TwoElectronOnePhoton.Settings;`
    
            multipoles::=..,        gauges=..,                printBefore=..,         lineSelection=..,         
            photonEnergyShift=..,   green=..) 
                        
        ... constructor for modifying the given TwoElectronOnePhoton.Settings by 'overwriting' the previously selected parameters.
    """
    function Settings(set::TwoElectronOnePhoton.Settings;    
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,    gauges::Union{Nothing,Array{UseGauge}}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,                   lineSelection::Union{Nothing,LineSelection}=nothing, 
        photonEnergyShift::Union{Nothing,Float64}=nothing,          green::Union{Nothing,Array{AtomicState.GreenChannel,1}}=nothing)
        
        if  multipoles          == nothing   multipolesx          = set.multipoles          else  multipolesx          = multipoles         end 
        if  gauges              == nothing   gaugesx              = set.gauges              else  gaugesx              = gauges             end 
        if  printBefore         == nothing   printBeforex         = set.printBefore         else  printBeforex         = printBefore        end 
        if  lineSelection       == nothing   lineSelectionx       = set.lineSelection       else  lineSelectionx       = lineSelection      end 
        if  photonEnergyShift   == nothing   photonEnergyShiftx   = set.photonEnergyShift   else  photonEnergyShiftx   = photonEnergyShift  end 
        if  green               == nothing   greenx               = set.green               else  greenx               = green              end 
         
        Settings( multipolesx, gaugesx, printBeforex, lineSelectionx, photonEnergyShiftx, greenx)
    end


    # `Base.show(io::IO, settings::TwoElectronOnePhoton.Settings)`  
    #    ... prepares a proper printout of the variable settings::TwoElectronOnePhoton.Settings.
    function Base.show(io::IO, settings::TwoElectronOnePhoton.Settings) 
        println(io, "multipoles:             $(settings.multipoles)  ")
        println(io, "gauges:                 $(settings.gauges)  ")
        println(io, "printBefore:            $(settings.printBefore)  ")
        println(io, "lineSelection:          $(settings.lineSelection)  ")
        println(io, "photonEnergyShift:      $(settings.photonEnergyShift)  ")
        println(io, "green:                  $(settings.green)  ")
    end


    """
    `struct  TwoElectronOnePhoton.Channel`  
        ... defines a type for a single radiative emission/absorption channel that specifies the multipole, gauge and amplitude.

        + multipole         ::EmMultipole        ... Multipole of the photon emission/absorption.
        + gauge             ::EmGauge            ... Gauge for dealing with the (coupled) radiation field.
        + symmetry          ::LevelSymmetry      ... total angular momentum and parity of the intermediate Green function levels.
        + amplitude         ::Complex{Float64}   ... Amplitude of this multiple channel.
    """
    struct Channel 
        multipole           ::EmMultipole
        gauge               ::EmGauge
        symmetry            ::LevelSymmetry
        amplitude           ::Complex{Float64}
    end 


    # `Base.show(io::IO, channel::TwoElectronOnePhoton.Channel)`  
    #    ... prepares a proper printout of the variable channel::TwoElectronOnePhoton.Channel.
    function Base.show(io::IO, channel::TwoElectronOnePhoton.Channel) 
        print(io, "TwoElectronOnePhoton.Channel($(channel.multipole), $(channel.gauge), amp = $(channel.amplitude)) ") 
    end


    """
    `struct  TwoElectronOnePhoton.Line`  
        ... defines a type for a TEOP line that may include the definition of sublines and their corresponding amplitudes.

        + initialLevel   ::Level               ... initial-(state) level
        + finalLevel     ::Level               ... final-(state) level
        + omega          ::Float64             ... Transition frequency of this line; can be shifted w.r.t. the level energies.
        + teopRate       ::EmProperty          ... Total TEOP rate of this line.
        + channels       ::Array{TwoElectronOnePhoton.Channel,1}  ... List of TEOP (photon) channels
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        omega            ::Float64
        teopRate         ::EmProperty
        channels         ::Array{TwoElectronOnePhoton.Channel,1}
    end 


    """
    `TwoElectronOnePhoton.Line()`  
        ... constructor an empty TEOP line.
    """
    function Line()
       Line(Level(), Level(), 0., EmProperty(0., 0.), TwoElectronOnePhoton.Channel[])
    end


    # `Base.show(io::IO, line::TwoElectronOnePhoton.Line)`  ... prepares a proper printout of the variable line::TwoElectronOnePhoton.Line.
    function Base.show(io::IO, line::TwoElectronOnePhoton.Line) 
        println(io, "initialLevel:         $(line.initialLevel)  ")
        println(io, "finalLevel:           $(line.finalLevel)  ")
        println(io, "omega:                $(line.omega)  ")
        println(io, "teopRate:             $(line.teopRate)  ")
        println(io, "channels:             $(line.channels)  ")
    end


    """
    `TwoElectronOnePhoton.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                                       settings::TwoElectronOnePhoton.Settings; output=true)`  
        ... to compute the TEOP transition amplitudes and rates as requested by the given settings. A list of 
            lines::Array{TwoElectronOnePhoton.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::TwoElectronOnePhoton.Settings; output=true) 
        # Define a common subshell list for both multiplets
        subshellList = Basics.generate("subshells: ordered list for two bases", finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
        Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
        error("Code need to be adapted.")
        #
        #
        println("")
        printstyled("TwoElectronOnePhoton.computeLines(): The computation of the TEOP amplitudes and rates starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        lines = TwoElectronOnePhoton.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    TwoElectronOnePhoton.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = TwoElectronOnePhoton.Line[]
        for  line in lines
            newLine = TwoElectronOnePhoton.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        TwoElectronOnePhoton.displayRates(stdout, newLines, settings)
        if  settings.calcAnisotropy    TwoElectronOnePhoton.displayAnisotropies(stdout, newLines, settings)    end
        #
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   TwoElectronOnePhoton.displayRates(iostream, newLines, settings)
        end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end


    """
    `TwoElectronOnePhoton.determineChannels(finalLevel::Level, initialLevel::Level, settings::TwoElectronOnePhoton.Settings)`  
        ... to determine a list of TwoElectronOnePhoton.Channel for a transitions from the initial to final level and by taking into 
            account the particular settings of for this computation; an Array{TwoElectronOnePhoton.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::TwoElectronOnePhoton.Settings)
        channels = TwoElectronOnePhoton.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        error("Code need to be adapted.")
        #
        #
        for  mp in settings.multipoles
            if   AngularMomentum.isAllowedMultipole(symi, mp, symf)
                hasMagnetic = false
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      push!(channels, TwoElectronOnePhoton.Channel(mp, Basics.Coulomb,   0.) )
                    elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    push!(channels, TwoElectronOnePhoton.Channel(mp, Basics.Babushkin, 0.) )  
                    elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)                  push!(channels, TwoElectronOnePhoton.Channel(mp, Basics.Magnetic,  0.) );
                                                        hasMagnetic = true; 
                    end 
                end
            end
        end
        return( channels )  
    end


    """
    `TwoElectronOnePhoton.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::TwoElectronOnePhoton.Settings)`  
        ... to determine a list of TwoElectronOnePhoton Line's for transitions between the levels from the given initial- and 
            final-state multiplets and by taking into account the particular selections and settings for this computation; 
            an Array{TwoElectronOnePhoton.Line,1} is returned. Apart from the level specification, all physical properties are set 
            to zero during the initialization process.  
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::TwoElectronOnePhoton.Settings)
        lines = TwoElectronOnePhoton.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    omega = iLevel.energy - fLevel.energy   + settings.photonEnergyShift

                    channels = TwoElectronOnePhoton.determineChannels(fLevel, iLevel, settings) 
                    if   length(channels) == 0   continue   end
                    push!( lines, TwoElectronOnePhoton.Line(iLevel, fLevel, omega, EmProperty(0., 0.), EmProperty(0., 0.), true, channels) )
                end
            end
        end
        return( lines )
    end

end # module

