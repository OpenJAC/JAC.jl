
"""
`module  JAC.PhotoIonization`  
    ... a submodel of JAC that contains all methods for computing photoionization properties between some initial 
        and final-state multiplets.
"""
module PhotoIonization

    using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..Radial, ..Nuclear, ..ManyElectron, ..PhotoEmission, 
                  ..PlasmaShift, ..TableStrings

    """
    `struct  PhotoIonization.Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing photoionization lines.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}             ... List of photon energies [in user-selected units].  
        + electronEnergies        ::Array{Float64,1}             ... List of electron energies; usually only one of these lists are utilized. 
        + calcAnisotropy          ::Bool                         ... True, if the beta anisotropy parameters are to be calculated and false otherwise (o/w). 
        + calcPartialCs           ::Bool                         ... True, if partial cross sections are to be calculated and false otherwise.  
        + calcTensors             ::Bool                         ... True, if statistical tensors of the excited atom are to be calculated and false o/w. 
        + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
        + stokes                  ::ExpStokes                    ... Stokes parameters of the incident radiation.
    """
    struct Settings  <:  AbstractProcessSettings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        photonEnergies            ::Array{Float64,1} 
        electronEnergies          ::Array{Float64,1} 
        calcAnisotropy            ::Bool 
        calcPartialCs             ::Bool 
        calcTensors               ::Bool 
        printBefore               ::Bool
        lineSelection             ::LineSelection
        stokes                    ::ExpStokes
    end 


    """
    `PhotoIonization.Settings()`  ... constructor for the default values of photoionization line computations
    """
    function Settings()
        Settings(Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], Float64[], false, false, false, false, 
                 LineSelection(), Basics.ExpStokes())
    end


    """
    `PhotoIonization.Settings(set::PhotoIonization.Settings;`
    
            multipoles=..,          gauges=..,                  photonEnergies=..,          electronEnergies=..,     
            calcAnisotropy=..,      calcPartialCs0..,           calcTensors=..,             printBefore=..,             
            lineSelection=..,       stokes=..)
                        
        ... constructor for modifying the given PhotoIonization.Settings by 'overwriting' the previously selected parameters.
    """
    function Settings(set::PhotoIonization.Settings;    
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
        photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                electronEnergies::Union{Nothing,Array{Float64,1}}=nothing,       
        calcAnisotropy::Union{Nothing,Bool}=nothing,                            calcPartialCs::Union{Nothing,Bool}=nothing,      
        calcTensors::Union{Nothing,Bool}=nothing,                               printBefore::Union{Nothing,Bool}=nothing,  
        lineSelection::Union{Nothing,LineSelection}=nothing,                    stokes::Union{Nothing,ExpStokes}=nothing)  
        
        if  multipoles      == nothing   multipolesx      = set.multipoles        else  multipolesx      = multipoles       end 
        if  gauges          == nothing   gaugesx          = set.gauges            else  gaugesx          = gauges           end 
        if  photonEnergies  == nothing   photonEnergiesx  = set.photonEnergies    else  photonEnergiesx  = photonEnergies   end 
        if  electronEnergies== nothing   electronEnergiesx= set.electronEnergies  else  electronEnergiesx= electronEnergies end 
        if  calcAnisotropy  == nothing   calcAnisotropyx  = set.calcAnisotropy    else  calcAnisotropyx  = calcAnisotropy   end 
        if  calcPartialCs   == nothing   calcPartialCsx   = set.calcPartialCs     else  calcPartialCsx   = calcPartialCs    end 
        if  calcTensors     == nothing   calcTensorsx     = set.calcTensors       else  calcTensorsx     = calcTensors      end 
        if  printBefore     == nothing   printBeforex     = set.printBefore       else  printBeforex     = printBefore      end 
        if  lineSelection   == nothing   lineSelectionx   = set.lineSelection     else  lineSelectionx   = lineSelection    end 
        if  stokes          == nothing   stokesx          = set.stokes            else  stokesx          = stokes           end 

        Settings( multipolesx, gaugesx, photonEnergiesx, electronEnergiesx, calcAnisotropyx, calcPartialCsx, calcTensorsx, 
                  printBeforex, lineSelectionx, stokesx)
    end


    # `Base.show(io::IO, settings::PhotoIonization.Settings)`  ... prepares a proper printout of the variable settings::PhotoIonization.Settings.
    function Base.show(io::IO, settings::PhotoIonization.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "electronEnergies:         $(settings.electronEnergies)  ")
        println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
        println(io, "calcPartialCs:            $(settings.calcPartialCs)  ")
        println(io, "calcTensors:              $(settings.calcTensors)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "lineSelection:            $(settings.lineSelection)  ")
        println(io, "stokes:                   $(settings.stokes)  ")
    end


    """
    `struct  PhotoIonization.Channel`  
        ... defines a type for a photoionization channel to help characterize a single multipole and scattering (continuum) state 
            of many electron-states with a single free electron.

        + multipole      ::EmMultipole          ... Multipole of the photon absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Photoionization amplitude associated with the given channel.
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
    `struct  Line`  ... defines a type for a photoionization line that may include the definition of channels.

        + initialLevel   ::Level                  ... initial-(state) level
        + finalLevel     ::Level                  ... final-(state) level
        + electronEnergy ::Float64                ... Energy of the (outgoing free) electron.
        + photonEnergy   ::Float64                ... Energy of the absorbed photon.
        + crossSection   ::EmProperty             ... Cross section for this photoionization.
         + channels       ::Array{PhotoIonization.Channel,1}  ... List of PhotoIonization.Channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        electronEnergy   ::Float64
        photonEnergy     ::Float64
        crossSection     ::EmProperty
         channels         ::Array{PhotoIonization.Channel,1}
    end


    """
    `PhotoIonization.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
        ... constructor for an photoionization line between a specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, crossSection::EmProperty)
        Line(initialLevel, finalLevel, totalRate, 0., 0., crossSection, false, PhotoChannel[] )
    end


    # `Base.show(io::IO, line::PhotoIonization.Line)`  ... prepares a proper printout of the variable line::PhotoIonization.Line.
    function Base.show(io::IO, line::PhotoIonization.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "electronEnergy:    $(line.electronEnergy)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
    end


    """
    `PhotoIonization.amplitude(kind::String, channel::PhotoIonization.Channel, omega::Float64, continuumLevel::Level, 
                                   initialLevel::Level, grid::Radial.Grid)`  
        ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || O^(photoionization) || alpha_i J_i>  
            due to the electron-photon interaction for the given final and initial level, the partial wave of the outgoing 
            electron as well as the given multipole and gauge. A value::ComplexF64 is returned.
    """
    function amplitude(kind::String, channel::PhotoIonization.Channel, omega::Float64, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid)
        if      kind in [ "photoionization"]
        #-----------------------------------
            amplitude = PhotoEmission.amplitude("absorption", channel.multipole, channel.gauge, omega, continuumLevel, initialLevel, grid, 
                                                display=false, printout=false)
            ## amplitude = im^Basics.subshell_l(Subshell(101, channel.kappa)) * exp( -im*channel.phase ) * amplitude
            amplitude = exp( -im*channel.phase ) * amplitude
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `PhotoIonization.computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                     settings::PhotoIonization.Settings; printout::Bool=false)`  
        ... to compute all amplitudes and properties of the given line; a line::PhotoIonization.Line is returned for which the amplitudes and 
            properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                          settings::PhotoIonization.Settings; printout::Bool=false)
        newChannels = PhotoIonization.Channel[];;   contSettings = Continuum.Settings(false, nrContinuum);    csC = 0.;    csB = 0.
        for channel in line.channels
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
            newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            #
            #==  Additional printout for Kabachnik/Pfeiffer, June 2023
            printSummary, iostream = Defaults.getDefaults("summary flag/stream")
            if  printSummary 
                en = Defaults.convertUnits("energy: from atomic to eV", cOrbital.energy)
                subsh = string(cOrbital.subshell) * "     "
                println(iostream, "   channel = $channel ")    
                println(iostream, "   shell = $(subsh[4:10])  energy = $en eV   phase = $phase  gauge = $(channel.gauge)")     
            end ==#
            #
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
            amplitude  = PhotoIonization.amplitude("photoionization", newChannel, line.photonEnergy, newcLevel, newiLevel, grid)
            push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry, 
                                                        newChannel.phase, amplitude) )
            if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
            elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
            elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
            end
        end
        Ji2 = Basics.twice(line.initialLevel.J)
        ##x csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
        ##x csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") / line.photonEnergy / (Ji2 + 1)
        ##x csFactor     = 4 * pi^2 / Defaults.getDefaults("alpha") / line.photonEnergy / (Ji2 + 1)
        csFactor     = 8 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
        csFactor     = csFactor / 2.   # Not fully clear, arises likely from the Rydberg normalization
        ##  Correct for energy normalization 
        ##  if  line.electronEnergy < 2.0   csFactor = csFactor * (line.electronEnergy/2.0)^1.5     end
        crossSection = EmProperty(csFactor * csC, csFactor * csB)
        newLine = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, 
                                        crossSection, newChannels)

        return( newLine )
    end



    """
    `PhotoIonization.computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                                           settings::PlasmaShift.PhotoSettings)`  
        ... to compute all amplitudes and properties of the given line but for the given plasma model; 
            a line::PhotoIonization.Line is returned for which the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, settings::PlasmaShift.PhotoSettings)
        newChannels = PhotoIonization.Channel[];;   contSettings = Continuum.Settings(false, grid.NoPoints-50);    csC = 0.;    csB = 0.
        for channel in line.channels
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel)
            newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel)
            @warn "Adapt a proper continuum orbital for the plasma potential"
            cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
            @warn "Adapt a proper Auger amplitude for the plasma e-e interaction"
            amplitude = 1.0
            # amplitude  = PhotoIonization.amplitude("photoionization", channel, line.photonEnergy, newcLevel, newiLevel, grid)
            push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry, 
                                                        newChannel.phase, amplitude) )
            if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
            elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
            elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
            end
        end
        Ji2 = Basics.twice(line.initialLevel.J)
        csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
        crossSection = EmProperty(csFactor * csC, csFactor * csB)
        println("plasma-photo cs = $crossSection")
        newline = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, crossSection, true, newChannels)
        return( newline )
    end


    """
    `PhotoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                      settings::PhotoIonization.Settings; output::Bool=true)`  
        ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{PhotoIonization.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::PhotoIonization.Settings; output::Bool=true)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        
        println("")
        printstyled("PhotoIonization.computeLines(): The computation of photo-ionization and properties starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    PhotoIonization.displayLines(lines)    end
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = PhotoIonization.Line[]
        for  line in lines
            println("\n>> Calculate photoionization amplitudes and properties for line: $(line.initialLevel.index) - $(line.finalLevel.index) " *
                    "for the photon energy $(Defaults.convertUnits("energy: from atomic", line.photonEnergy)) " * Defaults.GBL_ENERGY_UNIT)
            newLine = PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings) 
            push!( newLines, newLine)
            #== Additional printout for each line, if required ... due to collaboration with Kabachnik/Pfeiffer, June 2023
            if  true  && printSummary   PhotoIonization.displayResultsDetailed(stdout,   newLine, settings)     
                                        PhotoIonization.displayResultsDetailed(iostream, newLine, settings)  end  ==#
        end
        # Print all results to screen
        PhotoIonization.displayResults(stdout, newLines, settings)
        PhotoIonization.displayTimeDelay(stdout, newLines, settings)
        if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)     
                           PhotoIonization.displayTimeDelay(iostream, newLines, settings)    end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end



    """
    `PhotoIonization.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                         settings::PhotoIonization.Settings; output=true, printout::Bool=true)`  
        ... to compute the photoionization transition amplitudes and all properties as requested by the given settings. The computations
            and printout is adapted for large cascade computations by including only lines with at least one channel and by sending
            all printout to a summary file only. A list of lines::Array{PhotoIonization.Lines} is returned.
    """
    function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                  settings::PhotoIonization.Settings; output=true, printout::Bool=true)
        
        lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        # if  settings.printBefore    PhotoIonization.displayLines(lines)    end  
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  en in settings.photonEnergies     maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
                          for  en in settings.electronEnergies   maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = PhotoIonization.Line[]
        for  (i,line)  in  enumerate(lines)
            if  rem(i,10) == 0    println("> Photo line $i:  ... calculated ")    end
            newLine = PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings, printout=printout) 
            ##x if  rem(i,10) == 0    println("> Photo line $i:  ... not calculated ")    end
            ##x newLine = PhotoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, Basics.EmProperty(1.0, 0.), 
            ##x                                false, PhotoIonization.Channel[] )
            push!( newLines, newLine)
        end
        # Print all results to a summary file, if requested
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)   end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end


    """
    `PhotoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                            settings::PlasmaShift.PhotoSettings; output::Bool=true)`  
        ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{PhotoIonization.Lines} is returned.
    """
    function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                 settings::PlasmaShift.PhotoSettings; output::Bool=true)
        println("")
        printstyled("PhotoIonization.computeLinesPlasma(): The computation of photo-ionization cross sections starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        photoSettings = PhotoIonization.Settings(settings.multipoles, settings.gauges, settings.photonEnergies, settings.electronEnergies, 
                            false, false, false, settings.printBefore, settings.selectLines, settings.selectedLines, Basics.ExpStokes() )
        
        lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, photoSettings)
        # Display all selected lines before the computations start
        if  settings.printBefore    PhotoIonization.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoIonization.Line[]
        for  line in lines
            newLine = PhotoIonization.computeAmplitudesPropertiesPlasma(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        PhotoIonization.displayResults(stdout, lines, photoSettings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   PhotoIonization.displayResults(iostream, lines, photoSettings)     end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end



    """
    `PhotoIonization.computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)`  
        ... to compute the partial photoionization cross section for initially unpolarized atoms by unpolarized plane-wave photons.
            A value::Float64 is returned.
    """
    function  computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)
        # Define an internal Racah expression for the summation over t, lambda
        function Racahexpr(kappa::Int64, Ji::AngularJ64, Jf::AngularJ64, Mf::AngularM64, J::AngularJ64, Jp::AngularJ64, 
                           L::Int64, Lp::Int64, p::Int64, pp::Int64)
            # Determine the allowed values of t
            t1 = Basics.oplus( AngularJ64(Lp), Jf);    t2 = Basics.oplus( AngularJ64(L), Jf);    tList = intersect(t1, t2)
            wb = 0.
            for  t  in tList
                for  lambda = -1:2:1
                    j = AngularMomentum.kappa_j(kappa);   Mf_lambda = Basics.add(AngularM64(lambda), Mf)
                    wb = wb + (1.0im * lambda)^p * (-1.0im * lambda)^pp * 
                         AngularMomentum.ClebschGordan( AngularJ64(Lp), AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                         AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                         AngularMomentum.Wigner_9j(j, Jp, Jf, J, Ji, AngularJ64(L), Jf, AngularJ64(Lp), t)
                end
            end
            
            return( wb )
        end
        
         wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J;
        kappaList = PhotoIonization.getLineKappas(line)
        for  kappa in kappaList
            for  cha  in line.channels
                if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)   continue    end
                J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
                for  chp  in line.channels  
                    if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)    continue    end
                    Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                    wa = wa + 1.0im^(L - Lp) * (-1)^(L + Lp) * AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) *  
                         Racahexpr(kappa, Ji, Jf, Mf, J, Jp, L, Lp, p, pp) * cha.amplitude * conj(chp.amplitude)
                end
            end
        end
        csFactor = 8 * pi^3 * Defaults.getDefaults("alpha") / (2*line.photonEnergy * (Basics.twice(Ji) + 1))
        wa       = csFactor * wa

        return( wa )
    end



    """
    `PhotoIonization.computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line, 
                                                             settings::PhotoIonization.Settings)`  
        ... to compute the statistical tensor of the photoion in its final level after the photoionization of initially unpolarized atoms 
            by plane-wave photons with given Stokes parameters (density matrix). A value::ComplexF64 is returned.
    """
    function  computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line, 
                                                  settings::PhotoIonization.Settings)
        wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J   
        kappaList = PhotoIonization.getLineKappas(line);    P1 = settings.stokes.P1;   P2 = settings.stokes.P2;   P3 = settings.stokes.P3
        for  kappa in kappaList
            j = AngularMomentum.kappa_j(kappa)
            for  cha  in line.channels
                if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)   continue    end
                J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
                for  chp  in line.channels  
                    if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)    continue    end
                    Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                    #
                    for  lambda = -1:2:1
                        for  lambdap = -1:2:1
                            if  lambda == lambdap   wb = (1.0 + 0.0im +lambda*P3)    else    wb = P1 - lambda * P2 * im    end
                            wa = wa + wb * 1.0im^(L - Lp + p - pp) * lambda^p * lambdap^pp *
                                 sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) ) *  
                                 AngularMomentum.phaseFactor([J, +1, Jp, +1, Jf, +1, Ji, +1, j, +1, AngularJ64(1)]) *
                                 AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), AngularJ64(Lp),  AngularM64(-lambda), 
                                                                    AngularJ64(k),  AngularM64(q)) *
                                 AngularMomentum.Wigner_6j(Jf, j, Jp, J, AngularJ64(k), Jf) * 
                                 AngularMomentum.Wigner_6j(Jp, Ji, AngularJ64(Lp), AngularJ64(L), AngularJ64(k), J) * 
                                 cha.amplitude * conj(chp.amplitude)
                        end
                    end
                end
            end
        end
        
        wa = pi / (Basics.twice(Ji) + 1) * wa
        return( wa )
    end

    

    """
    `PhotoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)`  
        ... to determine a list of photoionization Channel for a transitions from the initial to final level and by taking into account 
            the particular settings of for this computation; an Array{PhotoIonization.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)
        channels = PhotoIonization.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        if  Basics.UseCoulomb  in  settings.gauges   gaugeM = Basics.UseCoulomb    else   gaugeM = Basics.UseBabushkin    end
        for  mp in settings.multipoles
            ### for  gauge in settings.gauges
                symList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                ##x println("mp = $mp   symi = $symi   symList = $symList")
                for  symt in symList
                    kappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
                    for  kappa in kappaList
            for  gauge in settings.gauges
                        # Include further restrictions if appropriate
                        if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                            push!(channels, PhotoIonization.Channel(mp, Basics.Coulomb,   kappa, symt, 0., Complex(0.)) )
                        elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                            push!(channels, PhotoIonization.Channel(mp, Basics.Babushkin, kappa, symt, 0., Complex(0.)) )  
                        elseif string(mp)[1] == 'M'  &&   gauge == gaugeM                               
                            push!(channels, PhotoIonization.Channel(mp, Basics.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                        end 
                    end
                end
            end
        end
        return( channels )  
    end


    """
    `PhotoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)`  
        ... to determine a list of PhotoIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
            and  by taking into account the particular selections and settings for this computation; an Array{PhotoIonization.Line,1} 
            is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)
        lines = PhotoIonization.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    # Add lines for all photon energies
                    for  omega in settings.photonEnergies
                        # Photon energies are still in 'pre-defined' units; convert to Hartree
                        omega_au = Defaults.convertUnits("energy: to atomic", omega)
                        energy   = omega_au - (fLevel.energy - iLevel.energy)
                        if  energy < 0.    continue   end  
                        channels = PhotoIonization.determineChannels(fLevel, iLevel, settings) 
                        push!( lines, PhotoIonization.Line(iLevel, fLevel, energy, omega_au, EmProperty(0., 0.), channels) )
                    end
                    # Add lines for all electron energies
                    for  en in settings.electronEnergies
                        # Electron energies are still in 'pre-defined' units; convert to Hartree
                        energy_au = Defaults.convertUnits("energy: to atomic", en)
                        omega     = energy_au + (fLevel.energy - iLevel.energy)
                        if  energy_au < 0.    continue   end  
                        channels = PhotoIonization.determineChannels(fLevel, iLevel, settings) 
                        push!( lines, PhotoIonization.Line(iLevel, fLevel, energy_au, omega, EmProperty(0., 0.), channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `PhotoIonization.displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})`  
        ... to display the calculated data, ordered by the initial levels and the photon energies involved.
            Neat tables of all initial levels and photon energies as well as all associated cross sections are printed
            but nothing is returned otherwise.
    """
    function  displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})
        # Extract and display all initial levels by their total energy and symmetry
        energies = Float64[]
        for  line  in  lines    push!(energies, line.initialLevel.energy)  end     
        energies = unique(energies);    energies = sort(energies)
        println(stream, "\n  Initial levels, available in the given photoionization line data:")
        println(stream, "\n  ", TableStrings.hLine(42))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(12, "Level energy"   ; na=3);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(42))
        for  energy in energies
            isNew = true
            for  line in lines
                if  line.initialLevel.energy == energy  &&  isNew
                    sa  = "  ";    sym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity )
                    sa = sa * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
                    sa = sa * TableStrings.center(10, string(sym); na=4)
                    sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * "    "
                    println(stream, sa);   isNew = false
                end
            end
        end
        #
        # Extract and display photon energies for which line data are available
        omegas = Float64[]
        for  line  in  lines    push!(omegas, line.photonEnergy)   end;    
        omegas = unique(omegas);    omegas = sort(omegas)
        sa  = "    ";   
        for   omega in omegas  sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "   "      end
        #
        println(stream, "\n  Photon energies $(TableStrings.inUnits("energy")) for which given photoionization line data are given:")
        println(stream, "\n" * sa)
        #
        # Extract and display total cross sections, ordered by the initial levels and photon energies, for which line data are available
        println(stream, "\n  Initial level, photon energies and total cross sections for given photoionization line data: \n")
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(12, "Omega"   ; na=3);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(26, "Cou--Total CS--Bab"   ; na=5);               
        sb = sb * TableStrings.center(26,TableStrings.inUnits("cross section"); na=5)
        sa = sa * TableStrings.center(18, "Electron energies"   ; na=3);               
        sb = sb * TableStrings.center(18,TableStrings.inUnits("energy"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(107))
        #
        if  length(lines) == 0    println("  >>> No photoionization data in lines here ... return ");    return(nothing)    end 
        #
        wLine = lines[1]
        for  energy in energies
            for  omega in omegas
                cs = Basics.EmProperty(0.);   electronEnergies = Float64[]
                for  line in lines
                    if  line.initialLevel.energy == energy  &&   line.photonEnergy == omega     cs = cs + line.crossSection
                        push!(electronEnergies, line.electronEnergy); wLine = line   end
                end
                sa  = "  ";    sym = LevelSymmetry( wLine.initialLevel.J, wLine.initialLevel.parity )
                sa = sa * TableStrings.center(10, TableStrings.level(wLine.initialLevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Coulomb))   * "  " *
                          @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Babushkin)) * "     "
                for en in electronEnergies  sa = sa * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", en)) * "  " end
                println(stream, sa)
            end
        end
        
        
        return( nothing )
    end


    """
    `PhotoIonization.displayLines(lines::Array{PhotoIonization.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoIonization.Line,1})
        nx = 175
        println(" ")
        println("  Selected photoionization lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(10, "Energy_fi"; na=3);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "omega"; na=3);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        nchannels = 0
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            energy = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "   "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaMultipoleSymmetryList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge, 
                                                    line.channels[i].symmetry) )
                nchannels = nchannels + 1
            end
            ##x println("PhotoIonization-diplayLines-ad: kappaMultipoleSymmetryList = ", kappaMultipoleSymmetryList)
            wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
            sb = sa * wa[1];    println( sb )  
            for  i = 2:length(wa)
                sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx), "\n")
        println("  A total of $nchannels channels need to be calculated. \n")
        #
        return( nothing )
    end


    """
    `PhotoIonization.displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
        nx = 130
        println(stream, " ")
        println(stream, "  Total photoionization cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " * 
                                              TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        # 
        wx = 1.0 ## wx = 1.0  (Schippers, August'23; wx = 2.0; wx = pi/2)
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Babushkin))   * "                 "
            ##x sa = sa * @sprintf("%.6e", line.crossSection.Coulomb)     * "    "
            ##x sa = sa * @sprintf("%.6e", line.crossSection.Babushkin)   * "    "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        #==
        # Print, if useful, the sum of all cross sections (Schippers, August'23)
        tcs = Basics.EmProperty(0.);    for  line in lines   tcs = tcs + line.crossSection   end 
        sa = repeat(" ", 102)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Coulomb))     * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Babushkin))   * "                 "
        println(stream, sa)  ==#
        #
        if  settings.calcPartialCs  
            nx = 144 
            println(stream, " ")
            println(stream, "  Partial cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
            println(stream, " ")
            println(stream, "  ", TableStrings.hLine(nx))
                sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
            sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
            sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
            sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(12, "omega"     ; na=4)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
            sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
            sa = sa * TableStrings.center( 7, "M_f"; na=1);                                     sb = sb * TableStrings.hBlank(11)
            sa = sa * TableStrings.center(30, "Cou -- Partial cross section -- Bab"; na=3)      
            sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " * 
                                              TableStrings.inUnits("cross section"); na=3)
            println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
            #   
            for  line in lines
                sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                             fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
                en = line.finalLevel.energy - line.initialLevel.energy
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
                sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
                println(stream, sa)
                MfList = Basics.projections(line.finalLevel.J)
                for  Mf in MfList
                    sb  = TableStrings.hBlank(97)
                    wac = PhotoIonization.computePartialCrossSectionUnpolarized(Basics.Coulomb, Mf, line)
                    wab = PhotoIonization.computePartialCrossSectionUnpolarized(Basics.Babushkin, Mf, line)
                    sb  = sb * TableStrings.flushright( 8, string(Mf))                             * "       "
                    sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wac.re))     * "    "
                    sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wab.re))     * "    "
                    println(stream, sb)
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        #
        if  settings.calcTensors  
            nx = 144 
            println(stream, " ")
            println(stream, "  Reduced statistical tensors of the photoion in its final level after the photoionization ")
            println(stream, "  of initially unpolarized atoms by plane-wave photons with given Stokes parameters (density matrix):")
            println(stream, "\n     + tensors are printed for k = 0, 1, 2 and if non-zero only.")
            println(stream,   "     + Stokes parameters are:  P1 = $(settings.stokes.P1),  P2 = $(settings.stokes.P2),  P3 = $(settings.stokes.P3) ")
            println(stream, " ")
            println(stream, "  ", TableStrings.hLine(nx))
                sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
            sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
            sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
            sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(12, "omega"     ; na=4)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
            sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
            sa = sa * TableStrings.center(10, "k    q"; na=4);                                  sb = sb * TableStrings.hBlank(11)
            sa = sa * TableStrings.center(30, "Cou --  rho_kq (J_f)  -- Bab"; na=3)      
            println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
            #   
            for  line in lines
                sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                             fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
                en = line.finalLevel.energy - line.initialLevel.energy
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
                sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
                println(stream, sa)
                for  k = 0:2
                    for q = -k:k
                        sb   = TableStrings.hBlank(102)
                        rhoc = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, Basics.Coulomb,   line, settings)
                        rhob = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, Basics.Babushkin, line, settings)
                        sb   = sb * string(k) * " " * TableStrings.flushright( 4, string(q))             * "       "
                        sb   = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", rhoc.re))     * "    "
                        sb   = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", rhob.re))     * "    "
                        println(stream, sb)
                    end
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        return( nothing )
    end


    """
    `PhotoIonization.displayResultsDetailed(stream::IO, line::PhotoIonization.Line, settings::PhotoIonization.Settings)`  
        ... to list the detailed results, energies, etc. for the given line. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResultsDetailed(stream::IO, line::PhotoIonization.Line, settings::PhotoIonization.Settings)
        symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
        symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)
        aeffC = ComplexF64(0.);     aeffB = ComplexF64(0.);    
        beffC = Float64(0.);        beffB = Float64(0.);       deffC = Float64(0.);     deffB = Float64(0.)
        for  channel in line.channels   
            if      channel.gauge == Basics.Coulomb     aeffC = aeffC + channel.amplitude
                                                        beffC = beffC + abs(channel.amplitude)^2
                                                        deffC = deffC + abs(channel.amplitude)^2 * atan(channel.amplitude.im, channel.amplitude.re)
                                                        ## deffC = deffC + abs(channel.amplitude)^2 * angle(channel.amplitude)
            elseif  channel.gauge == Basics.Babushkin   aeffB = aeffB + channel.amplitude  
                                                        beffB = beffB + abs(channel.amplitude)^2
                                                        deffB = deffB + abs(channel.amplitude)^2 * atan(channel.amplitude.im, channel.amplitude.re)
                                                        ## deffB = deffB + abs(channel.amplitude)^2 * angle(channel.amplitude)
            else
            end
        end
        deffC = deffC / beffC;    deffB = deffB / beffB  
        #
        sa   = "\n  Results for PI line from the transition  $(line.initialLevel.index) -  $(line.finalLevel.index):  " *
               "  $symi  - $symf "
        
        println(stream, sa, "\n  ", TableStrings.hLine(length(sa)-3))
        println(stream, "\n  Photon energy               = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy)) *
                          "  " * TableStrings.inUnits("energy") )
        println(stream,   "  Electron energy             = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) *
                          "  " * TableStrings.inUnits("energy") )
        println(stream,   "  Total cross section         = " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb)) *
                          "  " * TableStrings.inUnits("cross section") * " (Coulomb gauge)" )
        println(stream,   "                                " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin)) *
                          "  " * TableStrings.inUnits("cross section") * " (Babushkin gauge)" )
        println(stream,   "  A^eff (Re, Im, abs, phi)    = " * @sprintf("%.4e", aeffC.re)   * "   " * @sprintf("%.4e", aeffC.im) * "   " * 
                                                               @sprintf("%.4e", abs(aeffC)) * "   " * @sprintf("%.4e", angle(aeffC)) * "  (Coulomb gauge)" )
        println(stream,   "                              = " * @sprintf("%.4e", aeffB.re)   * "   " * @sprintf("%.4e", aeffB.im) * "   " * 
                                                               @sprintf("%.4e", abs(aeffB)) * "   " * @sprintf("%.4e", angle(aeffB)) * "  (Babushkin gauge)" )
        println(stream,   "  A^r_eff, del_eff (non-coh)  = " * @sprintf("%.4e", sqrt(beffC))* "   " * @sprintf("%.4e", deffC)        * "  (Coulomb gauge)" )
        println(stream,   "                              = " * @sprintf("%.4e", sqrt(beffB))* "   " * @sprintf("%.4e", deffB)        * "  (Babushkin gauge)" )
        println(stream,   "  A^r_eff, del_eff (coherent) = " * @sprintf("%.4e", abs(aeffC)) * "   " * @sprintf("%.4e", atan(aeffC.im, aeffC.re)) * "  (Coulomb gauge)" )
        println(stream,   "                              = " * @sprintf("%.4e", abs(aeffB)) * "   " * @sprintf("%.4e", atan(aeffB.im, aeffB.re)) * "  (Babushkin gauge)" )
        println(stream, "\n  Kappa    Total J^P  Mp     Gauge               Amplitude          Real-Amplitude  Cross section (b)   Phase    atan()" *
                        "\n  " * TableStrings.hLine(112) )
                    
        for  channel in line.channels 
            Ji2 = Basics.twice(line.initialLevel.J)
            csFactor = 4 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
            # csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
            cs = csFactor * abs(channel.amplitude)^2
            sg = string(channel.gauge) * "      "
            sk = "    " * string(channel.kappa)
            sb = "    " * @sprintf("%.4e", channel.amplitude.re)
            sc = "    " * @sprintf("%.4e", channel.amplitude.im)
            sp = "    " * @sprintf("%.4e", channel.phase)
            sx = "    " * @sprintf("%.4e", atan(channel.amplitude.im, channel.amplitude.re) )
            sa = " " * sk[end-3:end] * "          " * string(channel.symmetry) * "    " * string(channel.multipole) * "     " * sg[1:11] * 
                 sb[end-12:end] * sc[end-12:end] * "    " * @sprintf("%.4e", abs(channel.amplitude)) *
                 "       " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs)) * "   " * sp[end-12:end] * 
                 "       " * sx[end-12:end]
            println(stream, sa) 
        end
        
        return( nothing )
    end


    """
    `PhotoIonization.displayTimeDelay(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`  
        ... to display the time delay for "two" lines which just refer to two neighbored energies but to the same initial
            and final level. The procedure terminates if more than two PhotoIonization.Lines occur or if they are not
            identical with regard the number, coupling and gauges of the channels. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayTimeDelay(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
        if        length(lines) != 2                                                        error("stop a")
        elseif    length(lines[1].channels)  !=  length(lines[1].channels)                  error("stop b")
        elseif    !(0 <= lines[2].electronEnergy - lines[1].electronEnergy <= 0.1)          error("stop c")
        elseif    LevelSymmetry(lines[1].initialLevel.J, lines[1].initialLevel.parity) != 
                  LevelSymmetry(lines[2].initialLevel.J, lines[2].initialLevel.parity)      error("stop d")
        elseif    LevelSymmetry(lines[1].finalLevel.J, lines[1].finalLevel.parity) != 
                  LevelSymmetry(lines[2].finalLevel.J, lines[2].finalLevel.parity)          error("stop e")
        end 
        #
        deltaTau = lines[2].electronEnergy - lines[1].electronEnergy;   
        meanTauNC = ComplexF64(0.);   meanTauNB = ComplexF64(0.);    meanTauNxC = ComplexF64(0.);   meanTauNxB = ComplexF64(0.); 
        meanTauDC = ComplexF64(0.);   meanTauDB = ComplexF64(0.);    meanTauDxC = ComplexF64(0.);   meanTauDxB = ComplexF64(0.); 
        #==
        # Calculate non-coherent amplitudes
        for  ch = 1:length(lines[1].channels)
            phase1 = lines[1].channels[ch].phase # + lines[1].channels[ch].kappa
            phase2 = lines[2].channels[ch].phase # + lines[2].channels[ch].kappa
            w1 = abs(lines[1].channels[ch].amplitude) * exp(-im*phase1)
            w2 = abs(lines[1].channels[ch].amplitude) * exp(-im*phase2)
            if      lines[1].channels[ch].gauge != lines[2].channels[ch].gauge   ||
                    lines[1].channels[ch].kappa != lines[2].channels[ch].kappa              error("stop f") 
            elseif  lines[1].channels[ch].gauge == Basics.Coulomb     
                    meanTauNC  = meanTauNC  + conj(lines[1].channels[ch].amplitude) / deltaTau *
                                              abs(lines[1].channels[ch].amplitude) *
                                              exp(im*( angle(lines[2].channels[ch].amplitude) - angle(lines[1].channels[ch].amplitude) ) )
                    meanTauDC  = meanTauDC  + conj(lines[1].channels[ch].amplitude) * lines[1].channels[ch].amplitude 
                    meanTauNxC = meanTauNxC + conj(w1) /  deltaTau * (w2 - w1)
                    meanTauDxC = meanTauDxC + conj(w1) * w1
            elseif  lines[1].channels[ch].gauge == Basics.Babushkin    
                    meanTauNB  = meanTauNB  + conj(lines[1].channels[ch].amplitude) / deltaTau *
                                              abs(lines[1].channels[ch].amplitude) *
                                              exp(im*( angle(lines[2].channels[ch].amplitude) - angle(lines[1].channels[ch].amplitude) ) )
                                                  ## (lines[2].channels[ch].amplitude - lines[1].channels[ch].amplitude )
                    meanTauDB  = meanTauDB  + conj(lines[1].channels[ch].amplitude) * lines[1].channels[ch].amplitude 
                    meanTauNxB = meanTauNxB + conj(w1) /  deltaTau * (w2 - w1)
                    meanTauDxB = meanTauDxB + conj(w1) * w1
            else
            end
        end
        @show "\naa", meanTauNC, meanTauDC, meanTauNB, meanTauDB
        @show "\nbb", meanTauNxC, meanTauDxC, meanTauNxB, meanTauDxB
        meanTauC  = -im * meanTauNC / meanTauDC;        meanTauB  = -im * meanTauNB / meanTauDB
        meanTauxC = -im * meanTauNxC / meanTauDxC;      meanTauxB = -im * meanTauNxB / meanTauDxB
        #
        ## For second energy
        ## meanTauxC = meanTauxC;   meanTauxB = meanTauxB   ==#
        #  
        # Calculate coherent amplitudes
        amplitude1C = ComplexF64(3.0e-6);   amplitude2C = ComplexF64(0.)
        amplitude1B = ComplexF64(0.);       amplitude2B = ComplexF64(0.)
        for  ch = 5:6 ## length(lines[1].channels)
            if      lines[1].channels[ch].gauge != lines[2].channels[ch].gauge   ||
                    lines[1].channels[ch].kappa != lines[2].channels[ch].kappa              error("stop f") 
            elseif  lines[1].channels[ch].gauge == Basics.Coulomb
                ## if  abs( angle(lines[1].channels[ch].amplitude) )  < 1.45  &&  abs( angle(lines[1].channels[ch].kappa) )  != 3
                ## if abs(lines[1].channels[ch].amplitude - lines[2].channels[ch].amplitude)  >  1.0e-4
                ##     amplitude1C = amplitude1C + lines[1].channels[ch].amplitude
                ##     amplitude2C = amplitude2C - lines[2].channels[ch].amplitude
                ## else
                    amplitude1C = amplitude1C + lines[1].channels[ch].amplitude
                    amplitude2C = amplitude2C + lines[2].channels[ch].amplitude
                ## end
            elseif  lines[1].channels[ch].gauge == Basics.Babushkin     
                ## if  abs( angle(lines[1].channels[ch].amplitude) )  < 1.45  &&  abs( angle(lines[1].channels[ch].kappa) )  != 3
                ## if abs(lines[1].channels[ch].amplitude - lines[2].channels[ch].amplitude)  >  1.0e-4
                ##     amplitude1B = amplitude1B + lines[1].channels[ch].amplitude
                ##     amplitude2B = amplitude2B - lines[2].channels[ch].amplitude
                ## else
                    amplitude1B = amplitude1B + lines[1].channels[ch].amplitude
                    amplitude2B = amplitude2B + lines[2].channels[ch].amplitude
                ## end
            else
            end
        end
        #
        @show amplitude2C, amplitude2B
        delta_deffC = angle(amplitude2C) - angle(amplitude1C);      meanTauyC = abs(delta_deffC) / deltaTau / 300.
        delta_deffB = angle(amplitude2B) - angle(amplitude1B);      meanTauyB = abs(delta_deffB) / deltaTau / 300.
        #
        symi = LevelSymmetry(lines[1].initialLevel.J, lines[1].initialLevel.parity)
        symf = LevelSymmetry(lines[1].finalLevel.J, lines[1].finalLevel.parity)
        sa   = "\n  Mean time delay of the transition  $(lines[1].initialLevel.index) -  $(lines[1].finalLevel.index):  " *
               "  $symi  - $symf "
        
        println(stream, sa, "\n  ", TableStrings.hLine(length(sa)-3))
        println(stream, "\n  Electron energies   = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", lines[1].electronEnergy)) * "   " *
                                                       @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", lines[2].electronEnergy)) *
                          "  " * TableStrings.inUnits("energy") )
        ##x println(stream,   "  < tau(E) >          = " * @sprintf("%.4e", meanTauC.re)   * "   " * @sprintf("%.4e", meanTauC.im) * "   " * "  (Coulomb gauge)" )
        ##x println(stream,   "                      = " * @sprintf("%.4e", meanTauB.re)   * "   " * @sprintf("%.4e", meanTauB.im) * "   " * "  (Babushkin gauge)" )
        ##  println(stream,   "  < tau(E) ... >      = " * @sprintf("%.4e", meanTauxC.re)   * "   " * @sprintf("%.4e", meanTauxC.im) * "   " * "  (Coulomb gauge)" )
        ##  println(stream,   "  abs(A_s)*e^-i*phi_s = " * @sprintf("%.4e", meanTauxB.re)   * "   " * @sprintf("%.4e", meanTauxB.im) * "   " * "  (Babushkin gauge)" )
        println(stream,   "  < tau(E) ... >      = " * @sprintf("%.4e", meanTauyC)   * "   " * "  (Coulomb gauge)" )
        println(stream,   "  d delta_eff/dE      = " * @sprintf("%.4e", meanTauyB)   * "   " * "  (Babushkin gauge)" )
                          
        
        return( nothing )
    end



    """
    `PhotoIonization.extractPhotonEnergies(lines::Array{PhotoIonization.Line,1})`  
        ... to extract all photon energies for which photoionization data and cross sections are provided by lines;
            an list of energies::Array{Float64,1} is returned.
    """
    function  extractPhotonEnergies(lines::Array{PhotoIonization.Line,1})
        pEnergies = Float64[]
        for  line  in  lines    push!(pEnergies, line.photonEnergy)  end     
        pEnergies = unique(pEnergies)
        pEnergies = sort(pEnergies)
        
        return( pEnergies )
    end


    """
    `PhotoIonization.extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)`  
        ... to extract from lines all those that refer to the given omega;
            a reduced list rLines::Array{PhotoIonization.Line,1} is returned.
    """
    function  extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)
        rLines = PhotoIonization.Line[]
        for  line  in  lines    
            if  line.photonEnergy == omega   @show  "aa", omega;    push!(rLines, line)  end  
        end
        
        return( rLines )
    end


    """
    `PhotoIonization.extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`  
        ... to extract from lines the total PI cross section that refer to the given omega and initial level;
            a cross section cs::EmProperty is returned.
    """
    function  extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
        cs = Basics.EmProperty(0.)
        for  line  in  lines    
            if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  &&  
                line.photonEnergy       == omega 
                cs = cs + line.crossSection
            end  
        end
       
        return( cs )
    end


    """
    `PhotoIonization.extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, shell::Shell, initialLevel)`  
        ... to extract from lines the total PI cross section that refer to the given omega and initial level and to 
            to the ionization of an electron from shell; a cross section cs::EmProperty is returned.
    """
    function  extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, shell::Shell, initialLevel)
        cs = Basics.EmProperty(0.)
        for  line  in  lines    
            if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  &&  
                line.photonEnergy       == omega 
                # Now determined of whether the photoionization refers to the given shell
                confi     = Basics.extractLeadingConfiguration(line.initialLevel)
                conff     = Basics.extractLeadingConfiguration(line.finalLevel)
                shellOccs = Basics.extractShellOccupationDifference(confi::Configuration, conff::Configuration)
                if  length(shellOccs) > 1  ||  shellOccs[1][2] < 0      error("stop a")   end
                if  shellOccs[1][1] == shell   cs = cs + line.crossSection      end  
            end  
        end
        
        return( cs )
    end


    """
    `PhotoIonization.interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`  
        ... to interpolate (or extrapolate) from lines the total PI cross section for any given omega and 
            initial level. The procedure applies a linear interpolation/extrapolation by just using the 
            cross sections from the two nearest (given) omega points; a cross section cs::EmProperty is returned.
    """
    function  interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
        # First determine for which omegas cross sections are available and which associated cross sections
        # are to be applied in the interpolation
        omegas     = PhotoIonization.extractPhotonEnergies(lines)
        oms, diffs = Basics.determineNearestPoints(omega, 2, omegas)
        if  oms[1] < oms[2]     om1 = oms[1];   om2 = oms[2];   diff = - diffs[1]  
        else                    om1 = oms[2];   om2 = oms[1];   diff = - diffs[2] 
        end
        cs1        = PhotoIonization.extractCrossSection(lines, om1, initialLevel)
        cs2        = PhotoIonization.extractCrossSection(lines, om2, initialLevel)
        if  omega < om1  &&  omega < om2   ||  omega > om1  &&  omega > om2
            @warn("No extrapolation of cross sections; cs = 0.");       
            return( Basics.EmProperty(0.) )
        end
        # Interpolate/extrapolate linearly with the cross section data for the two omegas
        wm         = 1 / (om2 - om1) * (cs2 - cs1)
        cs         = cs1  +  diff * wm
        ## @show om1, om2, diff, cs1.Coulomb, cs2.Coulomb, cs.Coulomb
        
        return( cs )
    end


    """
    `PhotoIonization.getLineKappas(line::PhotoIonization.Line)`  
        ... returns a list of kappa-values (partial waves) which contribute to the given line, to which one or several channels are 
            assigned. An kappaList::Array{Int64,1} is returned.
    """
    function getLineKappas(line::PhotoIonization.Line)
        kappaList = Int64[]
        for  ch  in line.channels
            kappaList = union(kappaList, ch.kappa)
        end
        
        return( kappaList )
    end

end # module
