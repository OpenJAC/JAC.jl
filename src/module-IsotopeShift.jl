
"""
`module  JAC.IsotopeShift`  ... a submodel of JAC that contains all methods for computing isotope-dependent properties for some level or 
                                between some initial and final-state multiplets; it is using JAC, JAC.ManyElectron, JAC.Radial.
"""
module IsotopeShift

    using Printf, JAC, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  IsotopeShift.Settings`  ... defines a type for the details and parameters of computing isotope-shift M and F parameters.

        + calcNMS                  ::Bool             ... True if mass-shift parameters M_nmn need to be calculated, and false otherwise.
        + calcSMS                  ::Bool             ... True if mass-shift parameters M_sms need to be calculated, and false otherwise.
        + calcF                    ::Bool             ... True if the field-shift parameter need to be calculated, and false otherwise.
        + printBeforeComputation   ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Level,1}   ... List of selected levels.
        + methodF                  ::String           ... Method to calculate the field-shift parameter F.
    """
    struct Settings 
        calcNMS                    ::Bool
        calcSMS                    ::Bool
        calcF                      ::Bool
        printBeforeComputation     ::Bool
        selectLevels               ::Bool
        selectedLevels             ::Array{Level,1}
        methodF                    ::String
    end 


    """
    `JAC.IsotopeShift.Settings()`  ... constructor for an `empty` instance of IsotopeSettings for the computation of isotope M and F parameters.
    """
    function Settings()
        Settings(false, false, false, false, false, Level[], "")
    end


    """
    `Base.show(io::IO, settings::IsotopeShift.Settings)`  ... prepares a proper printout of the variable settings::IsotopeShift.Settings.
   """
    function Base.show(io::IO, settings::IsotopeShift.Settings) 
        println(io, "calcNMS:                  $(settings.calcNMS)  ")
        println(io, "calcSMS:                  $(settings.calcSMS)  ")
        println(io, "calcF:                    $(settings.calcF)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
        println(io, "methodF:                  $(settings.methodF)  ")
    end


    
    """
    `struct  IsotopeShift.Outcome`  ... defines a type to keep the outcome of a isotope-shift computation, such as the K and F parameters
                                        as well other results.

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
    `JAC.IsotopeShift.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
    """
    function Outcome()
        Outcome(Level(), 0., 0., 0.,   0., 0., 0., 0.)
    end


    """
    `Base.show(io::IO, outcome::IsotopeShift.Outcome)`  ... prepares a proper printout of the variable outcome::IsotopeShift.Outcome.
    """
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
    `JAC.IsotopeShift.amplitude(kind::String, rLevel::Level, sLevel::Level, nm::JAC.Nuclear.Model, grid::Radial.Grid)` 
         ... to compute either the  H^(NMS),  H^(SMS,A),  H^(SMS,B) or  H^(SMS,C) normal and specific mass-shift amplitudes
             <alpha_r J_r || H^(A)) || alpha_s J_s>  for a given pair of levels. A value::ComplexF64 is returned.
    """
    function  amplitude(kind::String, rLevel::Level, sLevel::Level, nm::JAC.Nuclear.Model, grid::Radial.Grid)
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
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.hamiltonian_nms(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], nm, grid)
                        me = me + coeff.T  * sqrt( ja + 1) * tamp  
                    end 
                #
                elseif  kind == "H^(SMS,A) amplitude"
                #------------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[2]
                        if  coeff.nu != 1  continue   end
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.X1_smsA(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                                sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                        me = me + coeff.V * sqrt( ja + 1) * tamp  
                    end
                #
                elseif  kind == "H^(SMS,B) amplitude"
                #------------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[2]
                        if  coeff.nu != 1  continue   end
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.X1_smsB(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                                sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                        me = me + coeff.V * sqrt( ja + 1) * tamp  
                    end
                #
                elseif  kind == "H^(SMS,C) amplitude"
                #------------------------------------
                    wa = compute("angular coefficients: e-e, Ratip2013", rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa[2]
                        if  coeff.nu != 1  continue   end
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.X1_smsC(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
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
    `JAC.IsotopeShift.computeAmplitudesProperties(outcome::IsotopeShift.Outcome, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                                  settings::IsotopeShift.Settings) 
         ... to compute all amplitudes and properties for a given level; an outcome::IsotopeShift.Outcome is returned for which the amplitudes 
         and properties are now evaluated explicitly.
    """
    function  computeAmplitudesProperties(outcome::IsotopeShift.Outcome, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings)
        global JAC_counter
        Knms  = 0.0;    Ksms  = 0.0;    F  = 0.0;    
        amplitudeKnms  = 0.0;    amplitudeKsmsA  = 0.0;   amplitudeKsmsB  = 0.0;  amplitudeKsmsC  = 0.0;
        if  settings.calcNMS
            amplitudeKnms  = JAC.IsotopeShift.amplitude("H^(NMS)   amplitude", outcome.level, outcome.level, nm, grid)
        end
        
        if  settings.calcSMS
            amplitudeKsmsA = JAC.IsotopeShift.amplitude("H^(SMS,A) amplitude", outcome.level, outcome.level, nm, grid)
            amplitudeKsmsB = JAC.IsotopeShift.amplitude("H^(SMS,B) amplitude", outcome.level, outcome.level, nm, grid)
            amplitudeKsmsC = JAC.IsotopeShift.amplitude("H^(SMS,C) amplitude", outcome.level, outcome.level, nm, grid)
            Knms = amplitudeKnms / nm.mass
            Ksms = (amplitudeKsmsA + amplitudeKsmsB + amplitudeKsmsC) / nm.mass
        end
        
        if  settings.calcF
            ##x F = JAC.Radial.amplitude("Delta V^nuc amplitude/(R-R')", nm, outcome.level, outcome.level, grid)
            F = 1.0
        end
        newOutcome = IsotopeShift.Outcome( outcome.level, Knms, Ksms, F, 
                                           amplitudeKnms, amplitudeKsmsA, amplitudeKsmsB, amplitudeKsmsC )
        return( newOutcome )
    end


    #==
    """
    `JAC.IsotopeShift.computeMatrixMnms_old(basis::Basis, grid::Radial.Grid, settings::IsotopeShift.Settings)`  ... to compute the matrix 
         M_nms = (<csf_r|| H_nms ||csf_s>) of the normal mass shift for the given CSF basis; a (length(basis.csfs) x length(basis.csfs)}-dimensional
         matrix::Array{Float64,2) is returned. See Naze et al., CPC 184 (2013) 2187, section 4.1 for further details on the decomposition of 
         H_nms.  
    """
    function  computeMatrixMnms_old(basis::Basis, grid::Radial.Grid, settings::IsotopeShift.Settings) 
        n = length(basis.csfs)
  
        print("Compute M_nms matrix of dimension $n x $n for the given basis ...")
        matrix = zeros(Float64, n, n)
        for  r = 1:n
            for  s = 1:n
                wa = compute("angular coefficients: e-e, Ratip2013", basis.csfs[r], basis.csfs[s])
                me = 0.
                for  coeff in wa[1]
                    jj = JAC.subshell_2j(basis.orbitals[coeff.a].subshell)
                    me = me + coeff.T * sqrt( jj + 1) * 
                                        JAC.InteractionStrength.hamiltonian_nms(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                end
                matrix[r,s] = me
            end
        end 
        println("   ... done.")
        return( matrix )
    end



    """
    `JAC.IsotopeShift.computeMatrixMsms_old(k::Int64, basis::Basis, grid::Radial.Grid, settings::IsotopeShift.Settings)` ... to compute the matrix 
         M_sms (k=1,2,3) = (<csf_r|| H_sms,k ||csf_s>) of the specific mass shift for the given CSF basis; a (length(basis.csfs) x 
         length(basis.csfs)}-dimensional matrix::Array{Float64,2) is returned. See Naze et al., CPC 184 (2013) 2187, section 4.2 for further 
         details on the decomposition of H_sms.  
    """
    function  computeMatrixMsms_old(k::Int64, basis::Basis, grid::Radial.Grid, settings::IsotopeShift.Settings) 
        n = length(basis.csfs)
  
        print("Compute M_sms (k=$k) matrix of dimension $n x $n for the given basis ...")
        matrix = zeros(Float64, n, n)
        for  r = 1:n
            for  s = 1:n
                wa = compute("angular coefficients: e-e, Ratip2013", basis.csfs[r], basis.csfs[s])
                me = 0.
                for  coeff in wa[2]
                    jj = JAC.subshell_2j(basis.orbitals[coeff.a].subshell)
                    if       k == 1    wa = JAC.InteractionStrength.X1_sms1(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b], 
                                                                                      basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                    elseif   k == 2    wa = JAC.InteractionStrength.X1_sms2(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b], 
                                                                                      basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid) 
                    elseif   k == 3    wa = JAC.InteractionStrength.X1_sms3(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b], 
                                                                                      basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                    else   error("stop a")
                    end
                    me = me + coeff.V * sqrt( jj + 1) * wa
                end
                matrix[r,s] = me
            end
        end 
        println("   ... done.")
        return( matrix )
    end
    ==#


    """
    `JAC.IsotopeShift.computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings; output=true)` 
         ... to compute (as selected) the isotope-shift M and F parameters for the levels of the given multiplet and as specified by the given settings. 
         The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings; output=true)
        println("")
        printstyled("JAC.IsotopeShift.computeOutcomes(): The computation of the isotope-shift parameters starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = JAC.IsotopeShift.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    JAC.IsotopeShift.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = IsotopeShift.Outcome[]
        for  outcome in outcomes
            newOutcome = JAC.IsotopeShift.computeAmplitudesProperties(outcome, nm, grid, settings) 
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        JAC.IsotopeShift.displayResults(stdout, newOutcomes)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    JAC.IsotopeShift.displayResults(iostream, newOutcomes)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `JAC.IsotopeShift.determineOutcomes(multiplet::Multiplet, settings::IsotopeShift.Settings)`  ... to determine a list of Outcomes's for 
         the computation of the isotope-shift M and F parameters for the given multiplet. It takes into account the particular selections and 
         settings. An Array{IsotopeShift.Outcome,1} is returned. Apart from the level specification, all physical properties are set to zero 
         during the initialization process.
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
    `JAC.IsotopeShift.displayOutcomes(outcomes::Array{IsotopeShift.Outcome,1})`  ... to display a list of levels that have been selected 
         for the computations. A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{IsotopeShift.Outcome,1})
        println(" ")
        println("  Selected IsotopeShift levels:")
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
    `JAC.IsotopeShift.displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1})`  ... to display the energies, M_ms and F-parameters, etc. for the 
         selected levels. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{IsotopeShift.Outcome,1})
        println(stream, " ")
        println(stream, "  IsotopeShift parameters and amplitudes:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(102))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=5)              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=5)
        sa = sa * JAC.TableStrings.center(11, "K_nms"; na=4)              
        sb = sb * JAC.TableStrings.center(11, "[Hz]" ; na=4)
        sa = sa * JAC.TableStrings.center(11, "K_sms"; na=4)             
        sb = sb * JAC.TableStrings.center(11, "[Hz]" ; na=4)
        sa = sa * JAC.TableStrings.center(11, "K_ms "; na=4)              
        sb = sb * JAC.TableStrings.center(11, "[Hz]" ; na=4)
        sa = sa * JAC.TableStrings.center(11, " F ";   na=4)              
        sb = sb * JAC.TableStrings.center(11, "[..]" ; na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(102)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.level.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", outcome.level.energy))               * "    "
            sa = sa * @sprintf("%.5e", JAC.convert("energy: from atomic to Hz", outcome.Knms))                 * "    "
            sa = sa * @sprintf("%.5e", JAC.convert("energy: from atomic to Hz", outcome.Ksms))                 * "    "
            sa = sa * @sprintf("%.5e", JAC.convert("energy: from atomic to Hz", outcome.Knms + outcome.Ksms))  * "    "
            sa = sa * @sprintf("%.5e", JAC.convert("energy: from atomic to Hz", outcome.F))                    * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(102), "\n\n")
        #
        # Now printout the individual amplitudes in an extra table
        println(stream, "  ", JAC.TableStrings.hLine(158))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(26, "H^(nms) amplitude"    ; na=3);             sb = sb * JAC.TableStrings.hBlank(35)
        sa = sa * JAC.TableStrings.center(26, "H^(sms,A) amplitude"  ; na=3);             sb = sb * JAC.TableStrings.hBlank(35)
        sa = sa * JAC.TableStrings.center(26, "H^(sms,B) amplitude"  ; na=3);             sb = sb * JAC.TableStrings.hBlank(35)
        sa = sa * JAC.TableStrings.center(26, "H^(sms,C) amplitude"  ; na=3);             sb = sb * JAC.TableStrings.hBlank(35)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(158)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.level.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", outcome.level.energy))          * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKnms.re,  " ", outcome.amplitudeKnms.im)  * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsA.re, " ", outcome.amplitudeKsmsA.im) * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsB.re, " ", outcome.amplitudeKsmsB.im) * "    "
            sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsC.re, " ", outcome.amplitudeKsmsC.im) * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(158))

        return( nothing )
    end

end # module
