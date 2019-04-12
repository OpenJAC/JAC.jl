
"""
`module  JAC.PhotoIonization`  
    ... a submodel of JAC that contains all methods for computing photoionization properties between some initial and final-state 
        multiplets; it is using JAC, JAC.ManyElectron, JAC.Radial.
"""
module PhotoIonization

    using Printf, JAC, JAC.ManyElectron, JAC.Radial
    
    global JAC_counter = 0


    """
    `struct  PhotoIonization.Settings`  ... defines a type for the details and parameters of computing photoionization lines.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}             ... List of photon energies.  
        + calcAnisotropy          ::Bool                         ... True, if the beta anisotropy parameters are to be calculated and false otherwise. 
        + calcPartialCs           ::Bool                         ... True, if partial cross sections are to be calculated and false otherwise.  
        + calcTensors             ::Bool                         ... True, if the statistical tensors of the excited atom are to be calculated and 
                                                                     false otherwise. 
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
        + stokes                  ::ExpStokes                    ... Stokes parameters of the incident radiation.
    """
    struct Settings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        photonEnergies            ::Array{Float64,1} 
        calcAnisotropy            ::Bool 
        calcPartialCs             ::Bool 
        calcTensors               ::Bool 
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
        stokes                    ::ExpStokes
    end 


    """
    `JAC.PhotoIonization.Settings()`  ... constructor for the default values of photoionization line computations
    """
    function Settings()
        Settings(EmMultipole[E1], UseGauge[JAC.UseCoulomb], Float64[], false, false, false, false, false, Tuple{Int64,Int64}[])
    end


    """
    `Base.show(io::IO, settings::PhotoIonization.Settings)`  ... prepares a proper printout of the variable settings::PhotoIonization.Settings.
    """
    function Base.show(io::IO, settings::PhotoIonization.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
        println(io, "calcPartialCs:            $(settings.calcPartialCs)  ")
        println(io, "calcTensors:              $(settings.calcTensors)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
        println(io, "stokes:                   $(settings.stokes)  ")
    end


    """
    `struct  PhotoIonization.Channel`  ... defines a type for a photoionization channel to help characterize a single multipole and scattering 
                                           (continuum) state of many electron-states with a single free electron.

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
        + hasChannels    ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      free-electron energy, kappa, multipole, etc., or not.
        + channels       ::Array{PhotoIonization.Channel,1}  ... List of PhotoIonization.Channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        electronEnergy   ::Float64
        photonEnergy     ::Float64
        crossSection     ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{PhotoIonization.Channel,1}
    end


    """
    `JAC.PhotoIonization.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
        ... constructor for an photoionization line between a specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, crossSection::EmProperty)
        Line(initialLevel, finalLevel, totalRate, 0., 0., crossSection, false, PhotoChannel[] )
    end


    """
    `Base.show(io::IO, line::PhotoIonization.Line)`  ... prepares a proper printout of the variable line::PhotoIonization.Line.
    """
    function Base.show(io::IO, line::PhotoIonization.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "electronEnergy:    $(line.electronEnergy)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
    end


    """
    `JAC.PhotoIonization.amplitude(kind::String, channel::PhotoIonization.Channel, energy::Float64, continuumLevel::Level, 
                                   initialLevel::Level, grid::Radial.Grid)`  
        ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || O^(photoionization) || alpha_i J_i>  
            due to the electron-photon interaction for the given final and initial level, the partial wave of the outgoing electron as well as 
            the given multipole and gauge. A value::ComplexF64 is returned.
    """
    function amplitude(kind::String, channel::PhotoIonization.Channel, energy::Float64, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid)
        if      kind in [ "photoionization"]
        #-----------------------------------
            amplitude = JAC.Radiative.amplitude("absorption", channel.multipole, channel.gauge, energy, continuumLevel, initialLevel, grid)
            amplitude = im^JAC.subshell_l(Subshell(101, channel.kappa)) * exp( -im*channel.phase ) * amplitude
        else    error("stop b")
        end
        
        ##x println("photo amplitude = $amplitude")
        return( amplitude )
    end



    """
    `JAC.PhotoIonization.computeAmplitudesProperties(line::PhotoIonization.Line, nm::JAC.Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                     settings::PhotoIonization.Settings)`  
        ... to compute all amplitudes and properties of the given line; a line::PhotoIonization.Line is returned for which the amplitudes and 
            properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoIonization.Line, nm::JAC.Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                          settings::PhotoIonization.Settings)
        newChannels = PhotoIonization.Channel[];;   contSettings = JAC.Continuum.Settings(false, nrContinuum);    csC = 0.;    csB = 0.
        for channel in line.channels
            newiLevel = JAC.generateLevelWithSymmetryReducedBasis(line.initialLevel)
            newiLevel = JAC.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            newfLevel = JAC.generateLevelWithSymmetryReducedBasis(line.finalLevel)
            cOrbital, phase  = JAC.Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel  = JAC.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
            amplitude  = JAC.PhotoIonization.amplitude("photoionization", channel, line.photonEnergy, newcLevel, newiLevel, grid)
            push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry, 
                                                        newChannel.phase, amplitude) )
            if       channel.gauge == JAC.Coulomb     csC = csC + abs(amplitude)^2
            elseif   channel.gauge == JAC.Babushkin   csB = csB + abs(amplitude)^2
            elseif   channel.gauge == JAC.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
            end
        end
        Ji2 = JAC.AngularMomentum.twoJ(line.initialLevel.J)
        csFactor     = 4 * pi^2 * JAC.give("alpha") * line.photonEnergy / (2*(Ji2 + 1))
        crossSection = EmProperty(csFactor * csC, csFactor * csB)
        ##x println("photo cs = $crossSection")
        newLine = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, 
                                        crossSection, true, newChannels)
        ##x println("photo newLine = $newLine")
        return( newLine )
    end



    """
    `JAC.PhotoIonization.computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                                           settings::PlasmaShift.PhotoSettings)`  
        ... to compute all amplitudes and properties of the given line but for the given plasma model; 
            a line::PhotoIonization.Line is returned for which the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::PlasmaShift.PhotoSettings)
        newChannels = PhotoIonization.Channel[];;   contSettings = JAC.Continuum.Settings(false, grid.nr-50);    csC = 0.;    csB = 0.
        for channel in line.channels
            newiLevel = JAC.generateLevelWithSymmetryReducedBasis(line.initialLevel)
            newiLevel = JAC.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            newfLevel = JAC.generateLevelWithSymmetryReducedBasis(line.finalLevel)
            @warn "Adapt a proper continuum orbital for the plasma potential"
            cOrbital, phase  = JAC.Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel  = JAC.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
            @warn "Adapt a proper Auger amplitude for the plasma e-e interaction"
            amplitude = 1.0
            # amplitude  = JAC.PhotoIonization.amplitude("photoionization", channel, line.photonEnergy, newcLevel, newiLevel, grid)
            push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry, 
                                                        newChannel.phase, amplitude) )
            if       channel.gauge == JAC.Coulomb     csC = csC + abs(amplitude)^2
            elseif   channel.gauge == JAC.Babushkin   csB = csB + abs(amplitude)^2
            elseif   channel.gauge == JAC.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
            end
        end
        Ji2 = JAC.AngularMomentum.twoJ(line.initialLevel.J)
        csFactor     = 4 * pi^2 * JAC.give("alpha") * line.photonEnergy / (2*(Ji2 + 1))
        crossSection = EmProperty(csFactor * csC, csFactor * csB)
        println("plasma-photo cs = $crossSection")
        newline = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, crossSection, true, newChannels)
        return( newline )
    end


    """
    `JAC.PhotoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                      settings::PhotoIonization.Settings; output::Bool=true)`  
        ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{PhotoIonization.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                           settings::PhotoIonization.Settings; output::Bool=true)
        println("")
        printstyled("JAC.PhotoIonization.computeLines(): The computation of photo-ionization and properties starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        lines = JAC.PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.PhotoIonization.displayLines(lines)    end
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
        nrContinuum = JAC.Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = PhotoIonization.Line[]
        for  line in lines
            newLine = JAC.PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.PhotoIonization.displayResults(stdout, newLines, settings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.PhotoIonization.displayResults(iostream, newLines, settings)     end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `JAC.PhotoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                            settings::PlasmaShift.PhotoSettings; output::Bool=true)`  
        ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{PhotoIonization.Lines} is returned.
    """
    function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                 settings::PlasmaShift.PhotoSettings; output::Bool=true)
        println("")
        printstyled("JAC.PhotoIonization.computeLinesPlasma(): The computation of photo-ionization cross sections starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        photoSettings = JAC.PhotoIonization.Settings(settings.multipoles, settings.gauges, settings.photonEnergies, false, false, false,
                                                     settings.printBeforeComputation, settings.selectLines, settings.selectedLines)
        
        lines = JAC.PhotoIonization.determineLines(finalMultiplet, initialMultiplet, photoSettings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.PhotoIonization.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoIonization.Line[]
        for  line in lines
            newLine = JAC.PhotoIonization.computeAmplitudesPropertiesPlasma(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.PhotoIonization.displayResults(stdout, lines, photoSettings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.PhotoIonization.displayResults(iostream, lines, photoSettings)     end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end



    """
    `JAC.PhotoIonization.computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)`  
        ... to compute the partial photoionization cross section for initially unpolarized atoms by unpolarized plane-wave photons.
            A value::Float64 is returned.
    """
    function  computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)
        # Define an internal Racah expression for the summation over t, lambda
        function Racahexpr(kappa::Int64, Ji::AngularJ64, Jf::AngularJ64, Mf::AngularM64, J::AngularJ64, Jp::AngularJ64, 
                           L::Int64, Lp::Int64, p::Int64, pp::Int64)
            # Determine the allowed values of t
            t1 = JAC.oplus( AngularJ64(Lp), Jf);    t2 = JAC.oplus( AngularJ64(L), Jf);    tList = intersect(t1, t2)
            wb = 0.
            for  t  in tList
                for  lambda = -1:2:1
                    j = JAC.AngularMomentum.kappa_j(kappa);   Mf_lambda = JAC.add(AngularM64(lambda), Mf)
                    wb = wb + (1.0im * lambda)^p * (-1.0im * lambda)^pp * 
                         JAC.AngularMomentum.ClebschGordan( AngularJ64(Lp), AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                         JAC.AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                         JAC.AngularMomentum.Wigner_9j(j, Jp, Jf, J, Ji, AngularJ64(L), Jf, AngularJ64(Lp), t)
                end
            end
            
            return( wb )
        end
        
        if  !line.hasChannels   error("No channels are defined for the given PhotoIonization.line.")         end
        wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J;
        kappaList = JAC.PhotoIonization.getLineKappas(line)
        for  kappa in kappaList
            for  cha  in line.channels
                if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != JAC.Magnetic)   continue    end
                J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
                for  chp  in line.channels  
                    if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != JAC.Magnetic)    continue    end
                    Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                    wa = wa + 1.0im^(L - Lp) * (-1)^(L + Lp) * JAC.AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) *  
                         Racahexpr(kappa, Ji, Jf, Mf, J, Jp, L, Lp, p, pp) * cha.amplitude * conj(chp.amplitude)
                end
            end
        end
        csFactor = 8 * pi^3 * JAC.give("alpha") / (2*line.photonEnergy * (JAC.AngularMomentum.twoJ(Ji) + 1))
        wa       = csFactor * wa

        return( wa )
    end



    """
    `JAC.PhotoIonization.computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line, 
                                                             settings::PhotoIonization.Settings)`  
        ... to compute the statistical tensor of the photoion in its final level after the photoionization of initially unpolarized atoms 
            by plane-wave photons with given Stokes parameters (density matrix). A value::ComplexF64 is returned.
    """
    function  computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line, 
                                                  settings::PhotoIonization.Settings)
        if  !line.hasChannels   error("No channels are defined for the given PhotoIonization.line.")         end
        wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J   
        kappaList = JAC.PhotoIonization.getLineKappas(line);    P1 = settings.stokes.P1;   P2 = settings.stokes.P2;   P3 = settings.stokes.P3
        for  kappa in kappaList
            j = JAC.AngularMomentum.kappa_j(kappa)
            for  cha  in line.channels
                if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != JAC.Magnetic)   continue    end
                J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
                for  chp  in line.channels  
                    if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != JAC.Magnetic)    continue    end
                    Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                    #
                    for  lambda = -1:2:1
                        for  lambdap = -1:2:1
                            if  lambda == lambdap   wb = (1.0 + 0.0im +lambda*P3)    else    wb = P1 - lambda * P2 * im    end
                            wa = wa + wb * 1.0im^(L - Lp + p - pp) * lambda^p * lambdap^pp *
                                 sqrt( JAC.AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) ) *  
                                 JAC.AngularMomentum.phaseFactor([J, +1, Jp, +1, Jf, +1, Ji, +1, j, +1, AngularJ64(1)]) *
                                 JAC.AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), AngularJ64(Lp),  AngularM64(-lambda), 
                                                                    AngularJ64(k),  AngularM64(q)) *
                                 JAC.AngularMomentum.Wigner_6j(Jf, j, Jp, J, AngularJ64(k), Jf) * 
                                 JAC.AngularMomentum.Wigner_6j(Jp, Ji, AngularJ64(Lp), AngularJ64(L), AngularJ64(k), J) * 
                                 cha.amplitude * conj(chp.amplitude)
                        end
                    end
                end
            end
        end
        
        wa = pi / (JAC.AngularMomentum.twoJ(Ji) + 1) * wa
        return( wa )
    end

    

    """
    `JAC.PhotoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)`  
        ... to determine a list of photoionization Channel for a transitions from the initial to final level and by taking into account 
            the particular settings of for this computation; an Array{PhotoIonization.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)
        channels = PhotoIonization.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        if  JAC.UseCoulomb  in  settings.gauges   gaugeM = JAC.UseCoulomb    else   gaugeM = JAC.UseBabushkin    end
        for  mp in settings.multipoles
            for  gauge in settings.gauges
                symList = JAC.AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                ##x println("mp = $mp   symi = $symi   symList = $symList")
                for  symt in symList
                    kappaList = JAC.AngularMomentum.allowedKappaSymmetries(symt, symf)
                    for  kappa in kappaList
                        # Include further restrictions if appropriate
                        if     string(mp)[1] == 'E'  &&   gauge == JAC.UseCoulomb      
                            push!(channels, PhotoIonization.Channel(mp, JAC.Coulomb,   kappa, symt, 0., Complex(0.)) )
                        elseif string(mp)[1] == 'E'  &&   gauge == JAC.UseBabushkin    
                            push!(channels, PhotoIonization.Channel(mp, JAC.Babushkin, kappa, symt, 0., Complex(0.)) )  
                        elseif string(mp)[1] == 'M'  &&   gauge == gaugeM                               
                            push!(channels, PhotoIonization.Channel(mp, JAC.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                        end 
                    end
                end
            end
        end
        return( channels )  
    end


    """
    `JAC.PhotoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)`  
        ... to determine a list of PhotoIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
            and  by taking into account the particular selections and settings for this computation; an Array{PhotoIonization.Line,1} 
            is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)
        if    settings.selectLines    selectLines   = true;                         
                       selectedLines = JAC.determineSelectedLines(settings.selectedLines, initialMultiplet, finalMultiplet)
        else                          selectLines   = false
        end
    
        lines = PhotoIonization.Line[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !((i,f) in selectedLines )    continue   end
                for  omega in settings.photonEnergies
                    # Photon energies are still in 'pre-defined' units; convert to Hartree
                    omega_au = JAC.convert("energy: to atomic", omega)
                    energy   = omega_au - (finalMultiplet.levels[f].energy - initialMultiplet.levels[i].energy)
                    if  energy < 0    continue   end  

                    channels = JAC.PhotoIonization.determineChannels(finalMultiplet.levels[f], initialMultiplet.levels[i], settings) 
                    push!( lines, PhotoIonization.Line(initialMultiplet.levels[i], finalMultiplet.levels[f], energy, omega_au, 
                                                       EmProperty(0., 0.), true, channels) )
                end
            end
        end
        return( lines )
    end


    """
    `JAC.PhotoIonization.displayLines(lines::Array{PhotoIonization.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoIonization.Line,1})
        println(" ")
        println("  Selected photoionization lines:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(175))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * JAC.TableStrings.hBlank(18)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(10, "Energy_fi"; na=3);              
        sb = sb * JAC.TableStrings.center(10, JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.center(10, "omega"; na=3);              
        sb = sb * JAC.TableStrings.center(10, JAC.TableStrings.inUnits("energy"); na=2)
        sa = sa * JAC.TableStrings.center(12, "Energy e_p"; na=3);              
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * JAC.TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(175)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=3)
            energy = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", energy))              * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", line.photonEnergy))   * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", line.electronEnergy)) * "    "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaMultipoleSymmetryList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge, 
                                                    line.channels[i].symmetry) )
            end
            ##x println("PhotoIonization-diplayLines-ad: kappaMultipoleSymmetryList = ", kappaMultipoleSymmetryList)
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
    `JAC.PhotoIonization.displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
        println(stream, " ")
        println(stream, "  Total photoionization cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(130))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * JAC.TableStrings.hBlank(18)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * JAC.TableStrings.center(12,JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * JAC.TableStrings.hBlank(13)
        sa = sa * JAC.TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)      
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("cross section") * "          " * 
                                              JAC.TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(130)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
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
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.crossSection.Babushkin))   * "                 "
            sa = sa * @sprintf("%.6e", line.crossSection.Coulomb)     * "    "
            sa = sa * @sprintf("%.6e", line.crossSection.Babushkin)   * "    "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(130))
        #
        #
        if  settings.calcPartialCs   
            println(stream, " ")
            println(stream, "  Partial cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(144))
                sa = "  ";   sb = "  "
            sa = sa * JAC.TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * JAC.TableStrings.hBlank(18)
            sa = sa * JAC.TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * JAC.TableStrings.hBlank(22)
            sa = sa * JAC.TableStrings.center(12, "f--Energy--i"; na=4)               
            sb = sb * JAC.TableStrings.center(12,JAC.TableStrings.inUnits("energy"); na=4)
            sa = sa * JAC.TableStrings.center(12, "omega"     ; na=4)             
            sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=4)
            sa = sa * JAC.TableStrings.center(12, "Energy e_p"; na=3)             
            sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=3)
            sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * JAC.TableStrings.hBlank(13)
            sa = sa * JAC.TableStrings.center( 7, "M_f"; na=1);                                     sb = sb * JAC.TableStrings.hBlank(11)
            sa = sa * JAC.TableStrings.center(30, "Cou -- Partial cross section -- Bab"; na=3)      
            sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("cross section") * "          " * 
                                                  JAC.TableStrings.inUnits("cross section"); na=3)
            println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(144)) 
            #   
            for  line in lines
                sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                             fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=3)
                en = line.finalLevel.energy - line.initialLevel.energy
                sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", en))                  * "    "
                sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", line.photonEnergy))   * "    "
                sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", line.electronEnergy)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
                sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=2)
                println(stream, sa)
                MfList = JAC.projections(line.finalLevel.J)
                for  Mf in MfList
                    sb  = JAC.TableStrings.hBlank(97)
                    wac = JAC.PhotoIonization.computePartialCrossSectionUnpolarized(JAC.Coulomb, Mf, line)
                    wab = JAC.PhotoIonization.computePartialCrossSectionUnpolarized(JAC.Babushkin, Mf, line)
                    sb  = sb * JAC.TableStrings.flushright( 8, string(Mf))                             * "       "
                    sb  = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", wac.re))     * "    "
                    sb  = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", wab.re))     * "    "
                    println(stream, sb)
                end
            end
            println(stream, "  ", JAC.TableStrings.hLine(144))
        end
        #
        #
        if  settings.calcTensors   
            println(stream, " ")
            println(stream, "  Reduced statistical tensors of the photoion in its final level after the photoionization ")
            println(stream, "  of initially unpolarized atoms by plane-wave photons with given Stokes parameters (density matrix):")
            println(stream, "\n     + tensors are printed for k = 0, 1, 2 and if non-zero only.")
            println(stream,   "     + Stokes parameters are:  P1 = $(settings.stokes.P1),  P2 = $(settings.stokes.P2),  P3 = $(settings.stokes.P3) ")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(144))
                sa = "  ";   sb = "  "
            sa = sa * JAC.TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * JAC.TableStrings.hBlank(18)
            sa = sa * JAC.TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * JAC.TableStrings.hBlank(22)
            sa = sa * JAC.TableStrings.center(12, "f--Energy--i"; na=4)               
            sb = sb * JAC.TableStrings.center(12,JAC.TableStrings.inUnits("energy"); na=4)
            sa = sa * JAC.TableStrings.center(12, "omega"     ; na=4)             
            sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=4)
            sa = sa * JAC.TableStrings.center(12, "Energy e_p"; na=3)             
            sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("energy"); na=3)
            sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * JAC.TableStrings.hBlank(13)
            sa = sa * JAC.TableStrings.center(10, "k    q"; na=4);                                  sb = sb * JAC.TableStrings.hBlank(11)
            sa = sa * JAC.TableStrings.center(30, "Cou --  rho_kq (J_f)  -- Bab"; na=3)      
            println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(144)) 
            #   
            for  line in lines
                sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                             fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=3)
                en = line.finalLevel.energy - line.initialLevel.energy
                sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", en))                  * "    "
                sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", line.photonEnergy))   * "    "
                sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", line.electronEnergy)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
                sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=2)
                println(stream, sa)
                for  k = 0:2
                    for q = -k:k
                        sb   = JAC.TableStrings.hBlank(102)
                        rhoc = JAC.PhotoIonization.computeStatisticalTensorUnpolarized(k, q, JAC.Coulomb,   line, settings)
                        rhob = JAC.PhotoIonization.computeStatisticalTensorUnpolarized(k, q, JAC.Babushkin, line, settings)
                        sb   = sb * string(k) * " " * JAC.TableStrings.flushright( 4, string(q))             * "       "
                        sb   = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", rhoc.re))     * "    "
                        sb   = sb * @sprintf("%.6e", JAC.convert("cross section: from atomic", rhob.re))     * "    "
                        println(stream, sb)
                    end
                end
            end
            println(stream, "  ", JAC.TableStrings.hLine(144))
        end
        #
        return( nothing )
    end


    """
    `JAC.PhotoIonization.getLineKappas(line::PhotoIonization.Line)`  
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
