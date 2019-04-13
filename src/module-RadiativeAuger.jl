
"""
`module  JAC.RadiativeAuger`  ... a submodel of JAC that contains all methods for computing radiative Auger and autoionization amplitudes 
                                  and rates; it is using JAC, JAC.Radial, JAC.ImpactExcitation, JAC.Auger.
"""
module RadiativeAuger

    using Printf, JAC, JAC.ManyElectron, JAC.Radial, JAC.ImpactExcitation, JAC.AutoIonization
    global JAC_counter = 0


    """
    `struct  RadiativeAuger.Settings`  ... defines a type for the settings in estimating radiative-Auger and autoionization rates.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + NoEnergySharings        ::Int64                        ... Number of energy sharings that are used in the computations for each line.
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
        + minAugerEnergy          ::Float64                      ... Minimum energy of free (Auger) electrons to be included.
        + maxAugerEnergy          ::Float64                      ... Maximum energy of free (Auger) electrons to be included.
        + maxKappa                ::Int6464                      ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        NoEnergySharings          ::Int64
        printBeforeComputation    ::Bool
        selectLines               ::Bool  
        selectedLines             ::Array{Tuple{Int64,Int64},1}
        minAugerEnergy            ::Float64
        maxAugerEnergy            ::Float64
        maxKappa                  ::Int64
    end 


    """
    `JAC.RadiativeAuger.Settings()`  ... constructor for the default values of RadiativeAuger line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], 0, false, false, Array{Tuple{Int64,Int64},1}[], 0., 10e5, 100)
    end


    # `Base.show(io::IO, settings::RadiativeAuger.Settings)`  ... prepares a proper printout of the variable settings::RadiativeAuger.Settings.
    function Base.show(io::IO, settings::RadiativeAuger.Settings) 
        println(io, "multipoles:                   $(settings.multipoles)  ")
        println(io, "gauges:                       $(settings.gauges)  ")
        println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
        println(io, "printBeforeComputation:       $(settings.printBeforeComputation)  ")
        println(io, "selectLines:                  $(settings.selectLines)  ")
        println(io, "selectedLines:                $(settings.selectedLines)  ")
        println(io, "minAugerEnergy:               $(settings.minAugerEnergy)  ")
        println(io, "maxAugerEnergy:               $(settings.maxAugerEnergy)  ")
        println(io, "maxKappa:                     $(settings.maxKappa)  ")
    end


    """
    `struct  Channel`  ... defines a type for a RadiativeAuger channel to help characterize a scattering (continuum) state of many 
                           electron-states with a single free electron.

        + multipole      ::EmMultipole          ... Multipole of the photon absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
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
    `struct  Sharing`  ... defines a type for a RadiativeAuger sharing to help characterize energy sharing between the emitted photon and
             the scattering (continuum) state with a single free electron.

        + photonEnergy   ::Float64         ... Energy of the emitted photon.
        + electronEnergy ::Float64         ... Energy of the (outgoing free) electron.
        + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
        + hasChannels    ::Bool            ... Determines whether the individual scattering (sub-) channels are defined in terms of their free-
                                               electron energies, kappas and the total angular momentum/parity as well as the amplitude, or not.
        + channels       ::Array{RadiativeAuger.Channel,1}  ... List of RadiativeAuger channels of this line.
    """
    struct  Sharing
        photonEnergy     ::Float64
        electronEnergy   ::Float64
        differentialCs   ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{RadiativeAuger.Channel,1}
    end




    # RadiativeAuger line between an initial and final (bound-state) level
    """
    `struct  Line`  ... defines a type for a radiative Auger line that may include the definition of sublines and their corresponding amplitudes.

        + initialLevel   ::Level          ... initial-(state) level
        + finalLevel     ::Level          ... final-(state) level
        + totalRate      ::EmProperty     ... Total rate of this line.
        + hasSharings    ::Bool           ... Determines whether the individual energy sharings are defined in terms of their photon
                                              and electron energies as well as their channels , or not.
        + sharings       ::Array{RadiativeAuger.Sharing,1}  ... List of RadiativeAuger sharings of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        totalRate        ::EmProperty
        hasSharings      ::Bool
        sharings         ::Array{RadiativeAuger.Sharing,1}
    end 


    # `Base.show(io::IO, line::RadiativeAuger.Line)`  ... prepares a proper printout of the variable line::RadiativeAuger.Line.
     function Base.show(io::IO, line::RadiativeAuger.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "totalRate:              $(line.totalRate)  ")
        println(io, "hasSharings:            $(line.hasSharings)  ")
        println(io, "sharings:               $(line.sharings)  ")
    end


    """
    `JAC.RadiativeAuger.computeAmplitudesProperties(line::RadiativeAuger.Line, grid::Radial.Grid, settings::RadiativeAuger.Settings)`  
         ... to compute all amplitudes and properties of the given line; a line::RadiativeAuger.Line is returned for which the amplitudes 
         and properties have now been evaluated.
    """
    function  computeAmplitudesProperties(line::RadiativeAuger.Line, grid::Radial.Grid, settings::RadiativeAuger.Settings)
        global JAC_counter
        newSharings = RadiativeAuger.Sharing[]
        for sharing in line.sharings
            newChannels = RadiativeAuger.Channel[]
            for channel in sharing.channels
                # Generate a continuum orbital
                JAC_counter = JAC_counter + 1
                if   JAC_counter < 20   println("RadiativeAuger.computeAmplitudesProperties-aa: warning ... no coninuum orbital is generated.")  end
                phase  = 0.
                # Define a proper continuum basis from the finalLevel.basis and the continuum orbital
                JAC_counter = JAC_counter + 1
                if   JAC_counter < 20   println("RadiativeAuger.computeAmplitudesProperties-ab: warning ... no coninuum basis is generated.") end
                # Compute the transition matrix for the continuum and the initial-state basis
                JAC_counter = JAC_counter + 1
                if   JAC_counter < 20   println("RadiativeAuger.computeAmplitudesProperties-ac: warning ... no transition matrix is computed.") end
                # matrix    = JAC.RadiativeAuger.computeMatrix(channel.multipole, channel.gauge, line.omega, line.finalLevel.basis, 
                #                                              line.initialLevel.basis, grid, settings)
                # amplitude = line.finalLevel.mc * matrix * line.initialLevel.mc 
                amplitude = 1.0 
                push!( newChannels, RadiativeAuger.Channel( channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, amplitude) )
            end
            push!( newSharings, RadiativeAuger.Sharing( sharing.photonEnergy, sharing.electronEnergy, EmProperty(-1., -1.), true, newChannels) )
        end
        # Calculate the totalRate 
        JAC_counter = JAC_counter + 1
        if   JAC_counter < 20   println("RadiativeAuger.computeAmplitudesProperties-ba: warning ... totalRate set to -1.") end
        totalRate = EmProperty(-1., -1.)
        line = RadiativeAuger.Line( line.initialLevel, line.finalLevel, totalRate, true, newSharings)
        return( line )
    end


    """
    `JAC.RadiativeAuger.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                     settings::RadiativeAuger.Settings; output=true)`  ... to compute the radiative Auger transition amplitudes and 
         all properties as requested by the given settings. A list of lines::Array{RadiativeAuger.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::RadiativeAuger.Settings; output=true)
        println("")
        printstyled("JAC.RadiativeAuger.computeLines(): The computation of radiative Auger rates starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------ \n", color=:light_green)
        println("")
        #
        lines = JAC.RadiativeAuger.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.RadiativeAuger.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = RadiativeAuger.Line[]
        for  line in lines
            newLine = JAC.RadiativeAuger.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.RadiativeAuger.displayResults(lines)
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `JAC.RadiativeAuger.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RadiativeAuger.Settings)` ... to determine 
         a list of RadiativeAuger.Line's for transitions between levels from the initial- and final-state multiplets, and  by taking into account 
         the particular selections and settings for this computation; an Array{RadiativeAuger.Line,1} is returned. Apart from the level 
         specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RadiativeAuger.Settings)
        if    settings.selectLines    selectLines   = true;   selectedLines = JAC.determine("selected lines", settings.selectedLines)
        else                          selectLines   = false
        end
    
        lines = RadiativeAuger.Line[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !(haskey(selectedLines, (i,f)) )    continue   end
                energy    = initialMultiplet.levels[i].energy - finalMultiplet.levels[f].energy
                if  energy < settings.minAugerEnergy  ||  energy > settings.maxAugerEnergy    continue   end  

                channels = JAC.RadiativeAuger.determineSharingsAndChannels(finalMultiplet.levels[f], initialMultiplet.levels[i], energy, settings) 
                push!( lines, RadiativeAuger.Line(initialMultiplet.levels[i], finalMultiplet.levels[f], EmProperty(0., 0.,), true, channels) )
            end
        end
        return( lines )
    end


    """
    `JAC.RadiativeAuger.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::RadiativeAuger.Settings)`  
         ... to determine a list of RadiativeAuger Sharing's and Channel's for a transitions from the initial to final level and by taking into 
         account the particular settings of for this computation; an Array{RadiativeAuger.Sharing,1} is returned.
    """
    function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::RadiativeAuger.Settings)
        sharings  = RadiativeAuger.Sharing[];    eSharings = JAC.determineEnergySharings(energy, settings.NoEnergySharings) 
        for  es in eSharings
            pEnergy   = es[1];    eEnergy = es[2] 
            channels  = RadiativeAuger.Channel[];   
            symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
            for  mp in settings.multipoles
                for  gauge in settings.gauges
                    symList = JAC.AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                    for  symt in symList
                        kappaList = JAC.AngularMomentum.allowedKappaSymmetries(symt, symf)
                        for  kappa in kappaList
                            # Include further restrictions if appropriate
                            if     string(mp)[1] == 'E'  &&   gauge == JAC.UseCoulomb      
                                push!(channels, RadiativeAuger.Channel(mp, JAC.Coulomb,   kappa, symt, 0., Complex(0.)) )
                            elseif string(mp)[1] == 'E'  &&   gauge == JAC.UseBabushkin    
                                push!(channels, RadiativeAuger.Channel(mp, JAC.Babushkin, kappa, symt, 0., Complex(0.)) )  
                            elseif string(mp)[1] == 'M'                                
                                push!(channels, RadiativeAuger.Channel(mp, JAC.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                            end 
                        end
                    end
                end
            end
            push!(sharings, RadiativeAuger.Sharing(pEnergy, eEnergy, EmProperty(0., 0.), true, channels) )
        end
        return( sharings )  
    end


    """
    `JAC.RadiativeAuger.displayLines(lines::Array{RadiativeAuger.Line,1})`  ... to display a list of lines, sharings and channels that have been 
         selected due to the prior settings. A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{RadiativeAuger.Line,1})
        println(" ")
        println("  Selected radiative-Auger lines:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(170))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(34, "Energy  " * JAC.TableStrings.inUnits("energy"); na=4);              
        sb = sb * JAC.TableStrings.center(34, "i -- f        omega     e_Auger  "; na=4)
        sa = sa * JAC.TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * JAC.TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(170)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.5e", JAC.convert("energy: from atomic", energy)) * "  "
            #
            for  sharing  in  line.sharings
                sb =      @sprintf("%.4e", JAC.convert("energy: from atomic", sharing.photonEnergy))   * "  "
                sb = sb * @sprintf("%.4e", JAC.convert("energy: from atomic", sharing.electronEnergy)) * "    "
                kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
                for  i in 1:length(sharing.channels)
                    push!( kappaMultipoleSymmetryList, (sharing.channels[i].kappa, sharing.channels[i].multipole, sharing.channels[i].gauge, 
                                                        sharing.channels[i].symmetry) )
                end
                wa = JAC.TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
                sc = sa * sb * wa[1];    println( sc )  
                for  i = 2:length(wa)
                    sc = JAC.TableStrings.hBlank( length(sa*sb) ) * wa[i];    println( sc )
                end
            end
        end
        println("  ", JAC.TableStrings.hLine(170))
        #
        return( nothing )
    end


    """
    `JAC.RadiativeAuger.displayResults(lines::Array{RadiativeAuger.Line,1})`  ... to list all results, energies, rates, etc. of the selected 
         lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(lines::Array{RadiativeAuger.Line,1})
        println(" ")
        println("  Radiative-Auger rates:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(148))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(36, "Energy  " * JAC.TableStrings.inUnits("energy"); na=4);              
        sb = sb * JAC.TableStrings.center(36, "i -- f        omega     e_Auger  "; na=4)
        sa = sa * JAC.TableStrings.center(30, "Cou -- differ. rate -- Bab"; na=3)      
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("rate") * "          " * JAC.TableStrings.inUnits("rate"); na=3)
        sa = sa * JAC.TableStrings.center(30, "Cou -- total rate -- Bab"; na=3)      
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("rate") * "          " * JAC.TableStrings.inUnits("rate"); na=3)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(148)) 
        #  
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.5e", JAC.convert("energy: from atomic", energy)) * "    "
            #
            first = true
            for  sharing  in  line.sharings
                sb =      @sprintf("%.4e", JAC.convert("energy: from atomic", sharing.photonEnergy))   * "  "
                sb = sb * @sprintf("%.4e", JAC.convert("energy: from atomic", sharing.electronEnergy)) * "    "
                sb = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", sharing.differentialCs.Coulomb))     * "    "
                sb = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", sharing.differentialCs.Babushkin))   * "    "
                if  first                  first = false
                sb = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.totalRate.Coulomb))     * "    "
                sb = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.totalRate.Babushkin))   * "    "
                end 
                println(sa*sb)
            end
        end
        println("  ", JAC.TableStrings.hLine(148))
        #
        return( nothing )
    end

end # module
