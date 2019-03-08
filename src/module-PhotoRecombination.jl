
"""
`module  JAC.PhotoRecombination`  ... a submodel of JAC that contains all methods for computing photo-recomnbination, i.e. radiative 
                                      recombination and radiative electron capture properties between some initial and final-state multiplets; 
                                      it is using JAC, JAC.ManyElectron, JAC.Radial.
"""
module PhotoRecombination

    using Printf, JAC, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  Settings`  ... defines a type for the details and parameters of computing photo recombination lines.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + electronEnergies        ::Array{Float64,1}             ... List of electron energies.
        + ionEnergies             ::Array{Float64,1}             ... List of ion energies.
        + useIonEnergies          ::Bool                         ... Make use of ion energies in [MeV/u] to obtain the electron energies.
        + calcAnisotropy          ::Bool                         ... True, if the overall anisotropy is to be calculated.
        + calcTensors             ::Bool                         ... True, if the statistical tensors are to be calculated and false otherwise.
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
    """
    struct Settings
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        electronEnergies          ::Array{Float64,1} 
        ionEnergies               ::Array{Float64,1}
        useIonEnergies            ::Bool
        calcAnisotropy            ::Bool
        calcTensors               ::Bool 
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
    end 


    """
    `JAC.PhotoRecombination.Settings()`  ... constructor for the default values of photo recombination line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], Float64[], false, false, false, false, false, Tuple{Int64,Int64}[])
    end


    """
    `Base.show(io::IO, settings::PhotoRecombination.Settings)`  ... prepares a proper printout of the variable 
                                                                    settings::PhotoRecombination.Settings.
    """
    function Base.show(io::IO, settings::PhotoRecombination.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "electronEnergies:         $(settings.electronEnergies)  ")
        println(io, "ionEnergies:              $(settings.ionEnergies)  ")
        println(io, "useIonEnergies:           $(settings.useIonEnergies)  ")
        println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
        println(io, "calcTensors:              $(settings.calcTensors)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end


    """
    `struct  PhotoRecombination.Channel`  ... defines a type for a photorecombination channel to help characterize a single multipole 
                                              and scattering (continuum) state of many electron-states with a single free electron.

        + multipole      ::EmMultipole          ... Multipole of the photon emission/absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Rec amplitude associated with the given channel.
    """
    struct  Channel
        multipole        ::EmMultipole
        gauge            ::EmGauge
        kappa            ::Int64
        symmetry         ::LevelSymmetry
        phase            ::Float64
        amplitude        ::Complex{Float64}
    end


    """
    `struct  PhotoRecombination.Line`  ... defines a type for a Photorecombination line that may include the definition of channels.

        + initialLevel   ::Level                  ... initial-(state) level
        + finalLevel     ::Level                  ... final-(state) level
        + electronEnergy ::Float64                ... Energy of the (incoming free) electron.
        + photonEnergy   ::Float64                ... Energy of the emitted photon.
        + crossSection   ::EmProperty             ... Cross section for this electron capture.
        + hasChannels    ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      free-electron energy, kappa, multipole, etc., or not.
        + channels       ::Array{PhotoRecombination.Channel,1}    ... List of photorecombination channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        electronEnergy   ::Float64
        photonEnergy     ::Float64
        crossSection     ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{PhotoRecombination.Channel,1}
    end 


    """
    `JAC.PhotoRecombination.Line()`  ... constructor for an `empty instance` of a photorecombination line between a specified initial 
                                         and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level)
        Line(initialLevel::Level, finalLevel::Level, 0., 0., EmProperty(0., 0.), false, Channel[])
    end


    """
    `Base.show(io::IO, line::PhotoRecombination.Line)`  ... prepares a proper printout of the variable line::PhotoRecombination.Line.
    """
    function Base.show(io::IO, line::PhotoRecombination.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "electronEnergy:    $(line.electronEnergy)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end


    """
    `JAC.PhotoRecombination.amplitude(kind::String, channel::PhotoRecombination.Channel, energy::Float64, finalLevel::Level, 
                                      continuumLevel::Level, grid::Radial.Grid)`  
         ... to compute the kind = (photorecombination) amplitude  < alpha_f J_f || O^(photorecombination) || (alpha_i J_i, epsilon kappa) J_t>  
         due to the electron-photon interaction for the given final and continuum level, the partial wave of the outgoing electron as well as 
         the given multipole and gauge. A value::ComplexF64 is returned.
    """
    function amplitude(kind::String, channel::PhotoRecombination.Channel, energy::Float64, finalLevel::Level, continuumLevel::Level, grid::Radial.Grid)
        if      kind in [ "photorecombination"]
        #-----------------------------------
            amplitude = JAC.Radiative.amplitude("emission", channel.multipole, channel.gauge, energy, finalLevel, continuumLevel, grid)
            amplitude = im^JAC.subshell_l(Subshell(101, channel.kappa)) * exp( im*channel.phase ) * amplitude
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `JAC.PhotoRecombination.computeAmplitudesProperties(line::PhotoRecombination.Line, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::PhotoRecombination.Settings)`  
         ... to compute all amplitudes and properties of the given line; a line::PhotoRecombination.Line is returned for which the amplitudes 
         and properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoRecombination.Line, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::PhotoRecombination.Settings)
        newChannels = PhotoRecombination.Channel[];;   contSettings = JAC.Continuum.Settings(false, grid.nr-50);    csC = 0.;    csB = 0.
        for channel in line.channels
            newfLevel = JAC.generateLevelWithSymmetryReducedBasis(line.finalLevel)
            newfLevel = JAC.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newfLevel)
            newiLevel = JAC.generateLevelWithSymmetryReducedBasis(line.initialLevel)
            cOrbital, phase  = JAC.Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newiLevel, nm, grid, contSettings)
            newcLevel  = JAC.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newiLevel)
            newChannel = PhotoRecombination.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
            amplitude  = JAC.PhotoRecombination.amplitude("photorecombination", channel, line.photonEnergy, newfLevel, newcLevel, grid)
            push!( newChannels, PhotoRecombination.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry, 
                                                           newChannel.phase, amplitude) )
            if       channel.gauge == JAC.Coulomb     csC = csC + abs(amplitude)^2
            elseif   channel.gauge == JAC.Babushkin   csB = csB + abs(amplitude)^2
            elseif   channel.gauge == JAC.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
            end
        end
        Ji2 = JAC.AngularMomentum.twoJ(line.initialLevel.J)
        csFactor     = 8 * pi^3 * JAC.give("alpha")^3 * line.photonEnergy / (Ji2 + 1)
        crossSection = EmProperty(csFactor * csC, csFactor * csB)
        newline = PhotoRecombination.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, 
                                           crossSection, true, newChannels)
        return( line )
    end


    """
    `JAC.PhotoRecombination.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                         settings::PhotoRecombination.Settings; output::Bool=true)` ... to compute the photo recombination 
         transition amplitudes and all properties as requested by the given settings. A list of lines::Array{PhotoRecombination.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::PhotoRecombination.Settings; 
                           output::Bool=true)
        println("")
        printstyled("JACPhotoRecombination.computeLines(): The computation of photo-recombination properties starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------------ \n", color=:light_green)
        println("")
        lines = JAC.PhotoRecombination.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.PhotoRecombination.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoRecombination.Line[]
        for  line in lines
            newLine = JAC.PhotoRecombination.computeAmplitudesProperties(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.PhotoRecombination.displayResults(stdout, lines, settings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.PhotoRecombination.displayResults(iostream, lines, settings)    end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end

    

    """
    `JAC.PhotoRecombination.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombination.Settings)`  ... to 
         determine a list of RecChannel for a transitions from the initial to final level and by taking into account the particular settings 
         of for this computation; an Array{PhotoRecombination.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombination.Settings)
        channels = PhotoRecombination.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp in settings.multipoles
            for  gauge in settings.gauges
                symList = JAC.AngularMomentum.allowedMultipoleSymmetries(symf, mp)
                for  symt in symList
                    kappaList = JAC.AngularMomentum.allowedKappaSymmetries(symi, symt)
                    for  kappa in kappaList
                        # Include further restrictions if appropriate
                        if     string(mp)[1] == 'E'  &&   gauge == JAC.UseCoulomb      
                            push!(channels, PhotoRecombination.Channel(mp, JAC.Coulomb,   kappa, symt, 0., Complex(0.)) )
                        elseif string(mp)[1] == 'E'  &&   gauge == JAC.UseBabushkin    
                            push!(channels, PhotoRecombination.Channel(mp, JAC.Babushkin, kappa, symt, 0., Complex(0.)) )  
                        elseif string(mp)[1] == 'M'                                
                            push!(channels, PhotoRecombination.Channel(mp, JAC.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                        end 
                    end
                end
            end
        end
        return( channels )  
    end


    """
    `JAC.PhotoRecombination.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoRecombination.Settings)`  
         ... to determine a list of PhotoRecombination.Line's for transitions between levels from the initial- and final-state multiplets, and  
         by taking into account the particular selections and settings for this computation; an Array{PhotoRecombination.Line,1} is returned. 
         Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoRecombination.Settings)
        if    settings.selectLines    selectLines   = true;   
                       selectedLines = JAC.determineSelectedLines(settings.selectedLines, initialMultiplet, finalMultiplet)
        else                          selectLines   = false
        end
    
        lines = PhotoRecombination.Line[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !((i,f) in selectedLines )    continue   end
                ##x println("PhotoRecombination.determineLines-aa:  i = $i, f = $f")
                #
                if  settings.useIonEnergies
                    error("Not yet implemented: Given ion energies.")
                else
                    for  en in settings.electronEnergies
                        if  en < 0    continue   end 
                        # Electron energies are still in 'pre-defined' units; convert to Hartree
                        en_au = JAC.convert("energy: to atomic", en)
                        omega = en_au + initialMultiplet.levels[i].energy - finalMultiplet.levels[f].energy

                        channels = JAC.PhotoRecombination.determineChannels(finalMultiplet.levels[f], initialMultiplet.levels[i], settings) 
                        push!( lines, PhotoRecombination.Line(initialMultiplet.levels[i], finalMultiplet.levels[f], en_au, omega, 
                                                              EmProperty(0., 0.), true, channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `JAC.PhotoRecombination.displayLines(lines::Array{PhotoRecombination.Line,1})`  ... to display a list of lines and channels that have been 
         selected due to the prior settings. A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoRecombination.Line,1})
        println(" ")
        println("  Selected photorecombination lines:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(175))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(10, "Energy_if"; na=2);              
        sb = sb * JAC.TableStrings.center(10, JAC.TableStrings.inUnits("energy"); na=2)
        sa = sa * JAC.TableStrings.center(12, "Energy e_r"; na=1);              
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=2)
        sa = sa * JAC.TableStrings.center(10, "omega"; na=5);              
        sb = sb * JAC.TableStrings.center(10, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * JAC.TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(175)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", energy))              * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", line.photonEnergy))   * "    "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaMultipoleSymmetryList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge, 
                                                    line.channels[i].symmetry) )
            end
            wa = JAC.TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
            sb = sa * wa[1];    println( sb )  
            for  i = 2:length(wa)
                sb = JAC.TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
        end
        println("  ", JAC.TableStrings.hLine(175), "\n")
        #
        return( nothing )
    end


    """
    `JAC.PhotoRecombination.displayResults(stream::IO, lines::Array{PhotoRecombination.Line,1}, settings::PhotoRecombination.Settings)`  
         ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
         returned otherwise.
    """
    function  displayResults(stream::IO, lines::Array{PhotoRecombination.Line,1}, settings::PhotoRecombination.Settings)
        println(stream, " ")
        println(stream, "  Photorecombination cross sections for beta^2 * gamma^2 = 1:")
        println(stream, "  ")
        println(stream, "  ", JAC.TableStrings.hLine(133))
        sa = " ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"   ; na=1);                       sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"   ; na=3);                       sb = sb * JAC.TableStrings.hBlank(21)
        sa = sa * JAC.TableStrings.center(12, "i--Energy--f"; na=4)               
        sb = sb * JAC.TableStrings.center(12,JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * JAC.TableStrings.hBlank(13)
        sa = sa * JAC.TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)      
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("cross section") * "          " * 
                                              JAC.TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(133)) 
        #   
        for  line in lines
            sa  = " ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                          fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
            sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=2)
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.crossSection.Babushkin))   * "    "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(133))
        #
        #
        if  settings.calcAnisotropy  
            println(stream, " ")
            println(stream, "  Anisotropy angular parameters beta_k for  ...  ")
            println(stream, " ")
            println(stream, "  Not yet implemented !!")
        end
        #
        #
        if  settings.calcTensors   
            println(stream, " ")
            println(stream, "  Reduced statistical tensors of the recombined ion ...  ")
            println(stream, " ")
            println(stream, "  Not yet implemented !!")
        end
        #
        return( nothing )
    end

end # module
