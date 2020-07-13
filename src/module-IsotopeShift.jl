
"""
`module  JAC.IsotopeShift`  
    ... a submodel of JAC that contains all methods for computing isotope-dependent properties for some level(s).
"""
module IsotopeShift

    using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings
    global JAC_counter = 0


    """
    `struct  IsotopeShift.Settings`  
        ... defines a type for the details and parameters of computing isotope-shift M and F parameters.

        + calcNMS                  ::Bool             ... True if mass-shift parameters M_nmn need to be calculated, and false otherwise.
        + calcSMS                  ::Bool             ... True if mass-shift parameters M_sms need to be calculated, and false otherwise.
        + calcF                    ::Bool             ... True if the field-shift parameter need to be calculated, and false otherwise.
        + printBefore   ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Level,1}   ... List of selected levels.
        + methodF                  ::String           ... Method to calculate the field-shift parameter F.
    """
    struct Settings 
        calcNMS                    ::Bool
        calcSMS                    ::Bool
        calcF                      ::Bool
        printBefore     ::Bool
        selectLevels               ::Bool
        selectedLevels             ::Array{Level,1}
        methodF                    ::String
    end 


    """
    `IsotopeShift.Settings()`  
        ... constructor for an `empty` instance of IsotopeSettings for the computation of isotope M and F parameters.
    """
    function Settings()
        Settings(false, false, false, false, false, Level[], "")
    end


    # `Base.show(io::IO, settings::IsotopeShift.Settings)`  ... prepares a proper printout of the variable settings::IsotopeShift.Settings.
    function Base.show(io::IO, settings::IsotopeShift.Settings) 
        println(io, "calcNMS:                  $(settings.calcNMS)  ")
        println(io, "calcSMS:                  $(settings.calcSMS)  ")
        println(io, "calcF:                    $(settings.calcF)  ")
        println(io, "printBefore:   $(settings.printBefore)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
        println(io, "methodF:                  $(settings.methodF)  ")
    end


    
    """
    `struct  IsotopeShift.Outcome`  
        ... defines a type to keep the outcome of a isotope-shift computation, such as the K and F parameters as well other 
            results.

        + level                     ::Level              ... Atomic level to which the outcome refers to.
        + Knms                      ::Float64            ... K_nms parameter
        + Ksms                      ::Float64            ... K_sms parameter
        + F                         ::Float64            ... F parameter
        + amplitudeKnms             ::Complex{Float64}   ... K_nms amplitude
        + amplitudeKsmsA            ::Complex{Float64}   ... K_sms,A amplitude
        + amplitudeKsmsB            ::Complex{Float64}   ... K_sms,B amplitude
        + amplitudeKsmsC            ::Complex{Float64}   ... K_sms,C amplitude
    """
    struct Outcome 
        level                       ::Level 
        Knms                        ::Float64
        Ksms                        ::Float64
        F                           ::Float64
        amplitudeKnms               ::Complex{Float64}
        amplitudeKsmsA              ::Complex{Float64}
        amplitudeKsmsB              ::Complex{Float64}
        amplitudeKsmsC              ::Complex{Float64}
    end 


    """
    `IsotopeShift.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
    """
    function Outcome()
        Outcome(Level(), 0., 0., 0.,   0., 0., 0., 0.)
    end


    # `Base.show(io::IO, outcome::IsotopeShift.Outcome)`  ... prepares a proper printout of the variable outcome::IsotopeShift.Outcome.
    function Base.show(io::IO, outcome::IsotopeShift.Outcome) 
        println(io, "level:                   $(outcome.level)  ")
        println(io, "Knms:                    $(outcome.Knms)  ")
        println(io, "Ksms:                    $(outcome.Ksms)  ")
        println(io, "F:                       $(outcome.F)  ")
        println(io, "amplitudeKnms:           $(outcome.amplitudeKnms)  ")
        println(io, "amplitudeKsmsA:          $(outcome.amplitudeKsmsA)  ")
        println(io, "amplitudeKsmsB:          $(outcome.amplitudeKsmsB)  ")
        println(io, "amplitudeKsmsC:          $(outcome.amplitudeKsmsC)  ")
    end


    """
    `IsotopeShift.amplitude(kind::String, rLevel::Level, sLevel::Level, nm::Nuclear.Model, grid::Radial.Grid)` 
        ... to compute either the  H^(NMS),  H^(SMS,A),  H^(SMS,B) or  H^(SMS,C) normal and specific mass-shift amplitudes
            <alpha_r J_r || H^(A)) || alpha_s J_s>  for a given pair of levels. A value::ComplexF64 is returned.
    """
    function  amplitude(kind::String, rLevel::Level, sLevel::Level, nm::Nuclear.Model, grid::Radial.Grid)
        nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
        printstyled("Compute mass-shift $(kind[1:9]) matrix of dimension $nr x $ns ... ", color=:light_green)
        #
        if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
        #
        for  r = 1:nr
            for  s = 1:ns
                me = 0.
                if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                    rLevel.parity != sLevel.parity    continue    
                end 
                #
                if      kind == "H^(NMS)   amplitude"
                #----------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[1]
                        ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = InteractionStrength.hamiltonian_nms(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], nm, grid)
                        me = me + coeff.T  * sqrt( ja + 1) * tamp  
                    end 
                #
                elseif  kind == "H^(SMS,A) amplitude"
                #------------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[2]
                        if  coeff.nu != 1  continue   end
                        ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = InteractionStrength.X1_smsA(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                            sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                        me = me + coeff.V * sqrt( ja + 1) * tamp  
                    end
                #
                elseif  kind == "H^(SMS,B) amplitude"
                #------------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[2]
                        if  coeff.nu != 1  continue   end
                        ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = InteractionStrength.X1_smsB(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                            sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                        me = me + coeff.V * sqrt( ja + 1) * tamp  
                    end
                #
                elseif  kind == "H^(SMS,C) amplitude"
                #------------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[2]
                        if  coeff.nu != 1  continue   end
                        ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = InteractionStrength.X1_smsC(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                            sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                        me = me + coeff.V * sqrt( ja + 1) * tamp  
                    end
                #
                else    error("stop a")
                end
                #
                matrix[r,s] = me
            end
        end
        printstyled("done.\n", color=:light_green)
        amplitude = transpose(rLevel.mc) * matrix * sLevel.mc 
        #
        return( amplitude )
    end



    """
    `IsotopeShift.computeAmplitudesProperties(outcome::IsotopeShift.Outcome, nm::Nuclear.Model, grid::Radial.Grid, 
                                                  settings::IsotopeShift.Settings) 
        ... to compute all amplitudes and properties for a given level; an outcome::IsotopeShift.Outcome is returned for 
            which the amplitudes and properties are now evaluated explicitly.
    """
    function  computeAmplitudesProperties(outcome::IsotopeShift.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings)
        global JAC_counter
        Knms  = 0.0;    Ksms  = 0.0;    F  = 0.0;    
        amplitudeKnms  = 0.0;    amplitudeKsmsA  = 0.0;   amplitudeKsmsB  = 0.0;  amplitudeKsmsC  = 0.0;
        if  settings.calcNMS
            amplitudeKnms  = IsotopeShift.amplitude("H^(NMS)   amplitude", outcome.level, outcome.level, nm, grid)
        end
        
        if  settings.calcSMS
            amplitudeKsmsA = IsotopeShift.amplitude("H^(SMS,A) amplitude", outcome.level, outcome.level, nm, grid)
            amplitudeKsmsB = IsotopeShift.amplitude("H^(SMS,B) amplitude", outcome.level, outcome.level, nm, grid)
            amplitudeKsmsC = IsotopeShift.amplitude("H^(SMS,C) amplitude", outcome.level, outcome.level, nm, grid)
            Knms = amplitudeKnms / nm.mass
            Ksms = (amplitudeKsmsA + amplitudeKsmsB + amplitudeKsmsC) / nm.mass
        end
        
        if  settings.calcF
            ##x F = Radial.amplitude("Delta V^nuc amplitude/(R-R')", nm, outcome.level, outcome.level, grid)
            F = 1.0
        end
        newOutcome = IsotopeShift.Outcome( outcome.level, Knms, Ksms, F, 
                                           amplitudeKnms, amplitudeKsmsA, amplitudeKsmsB, amplitudeKsmsC )
        return( newOutcome )
    end


    """
    `IsotopeShift.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings; output=true)` 
        ... to compute (as selected) the isotope-shift M and F parameters for the levels of the given multiplet and as 
            specified by the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings; output=true)
        println("")
        printstyled("IsotopeShift.computeOutcomes(): The computation of the isotope-shift parameters starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = IsotopeShift.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBefore    IsotopeShift.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = IsotopeShift.Outcome[]
        for  outcome in outcomes
            newOutcome = IsotopeShift.computeAmplitudesProperties(outcome, nm, grid, settings) 
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        IsotopeShift.displayResults(stdout, newOutcomes)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    IsotopeShift.displayResults(iostream, newOutcomes)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `IsotopeShift.determineOutcomes(multiplet::Multiplet, settings::IsotopeShift.Settings)`  
        ... to determine a list of Outcomes's for the computation of the isotope-shift M and F parameters for the given multiplet. 
            It takes into account the particular selections and settings. An Array{IsotopeShift.Outcome,1} is returned. Apart from 
            the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::IsotopeShift.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        outcomes = IsotopeShift.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, IsotopeShift.Outcome(multiplet.levels[i], 0., 0., 0.,   0., 0., 0., 0.) )
        end
        return( outcomes )
    end


    """
    `IsotopeShift.displayOutcomes(outcomes::Array{IsotopeShift.Outcome,1})`  
        ... to display a list of levels that have been selected for the computations. A small neat table of all 
            selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{IsotopeShift.Outcome,1})
        nx = 43
        println(" ")
        println("  Selected IsotopeShift levels:")
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
            sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
            println( sa )
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `IsotopeShift.displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1})`  
        ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{IsotopeShift.Outcome,1})
        nx = 102
        println(stream, " ")
        println(stream, "  IsotopeShift parameters and amplitudes:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  ";   sc = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5);  sc = sc * TableStrings.hBlank(45)
        sa = sa * TableStrings.center(11, "K_nms"; na=4)              
        sb = sb * TableStrings.center(11, "[Hz]" ; na=4)
        sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
        sa = sa * TableStrings.center(11, "K_sms"; na=4)             
        sb = sb * TableStrings.center(11, "[Hz]" ; na=4)
        sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
        sa = sa * TableStrings.center(11, "K_ms "; na=4)              
        sb = sb * TableStrings.center(11, "[Hz]" ; na=4)
        sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
        sa = sa * TableStrings.center(11, " F ";   na=4)              
        sb = sb * TableStrings.center(11, "[..]" ; na=4)
        println(stream, sa);    println(stream, sb);    println(stream, sc);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))               * "    "
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.Knms))                 * "    "
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.Ksms))                 * "    "
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.Knms + outcome.Ksms))  * "    "
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.F))                    * "    "
            println(stream, sa )
            sb = TableStrings.hBlank(47)
            sb = sb * @sprintf("%.5e", outcome.Knms * 1822.888 * 2)                                            * "    "
            sb = sb * @sprintf("%.5e", outcome.Ksms * 1822.888 * 2)                                            * "    "
            sb = sb * @sprintf("%.5e", (outcome.Knms + outcome.Ksms) * 1822.888 * 2)                           * "    "
            sb = sb * @sprintf("%.5e", outcome.F)                                                              * "    "
            println(stream, sb )
        end
        println(stream, "  ", TableStrings.hLine(nx), "\n\n")
        #
        # Now printout the individual amplitudes in an extra table
        nx = 158
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(26, "H^(nms) amplitude"    ; na=3);             sb = sb * TableStrings.hBlank(35)
        sa = sa * TableStrings.center(26, "H^(sms,A) amplitude"  ; na=3);             sb = sb * TableStrings.hBlank(35)
        sa = sa * TableStrings.center(26, "H^(sms,B) amplitude"  ; na=3);             sb = sb * TableStrings.hBlank(35)
        sa = sa * TableStrings.center(26, "H^(sms,C) amplitude"  ; na=3);             sb = sb * TableStrings.hBlank(35)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))          * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKnms.re,  " ", outcome.amplitudeKnms.im)  * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsA.re, " ", outcome.amplitudeKsmsA.im) * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsB.re, " ", outcome.amplitudeKsmsB.im) * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsC.re, " ", outcome.amplitudeKsmsC.im) * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))

        return( nothing )
    end

end # module
