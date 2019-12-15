
"""
`module  JAC.BasicsQZ`  
    ... a submodel of JAC that contains methods for setting-up and performing atomic structure computations,
        i.e. self-consistent-fields, configuration-interaction, etc.
"""
module BasicsQZ

    using Printf, Interact, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..Einstein, ..Hfs, 
          ..ManyElectron, ..Nuclear, ..PhotoEmission, ..Radial, ..TableStrings
    
    export dummyQZ

    """
    `Basics.read()`  ... reads in data from different files and sources.

    + `("CSF list: Grasp92", "cslFilename")`  
        ... reads in the CSF list as given by the (existing) .csl file cslFilename; an basis::Basis is returned with 
            basis.isDefined = true and with a proper sequence of orbitals, but with basis.orbitals = Orbital[].

    + `("orbital list: Grasp92", "orbFilename")`  
        ... read in the orbitals from a (formatted) .rwf file orbFilename; a list of orbitals::Array{OrbitalsR,1} 
            is returned with all subfields specified.
    """
    function Basics.read(sa::String, filename::String)
        if      sa == "CSF list: Grasp92"                          wa = readCslFileGrasp92(filename)
        elseif  sa == "orbital list: Grasp92"                      wa = readOrbitalFileGrasp92(filename)
        elseif  sa == "energies & mixing coefficients: Grasp92"    wa = readMixingFileGrasp92(filename)
        else    error("Unsupported keystring = $sa ")
        end

        return( wa )
    end


    """
    Basics.readCslFileGrasp92(filename::String)`  
        ... reads in the CSF list from a Grasp92 .csl file; a basis::Basis is returned.
    """
    function Basics.readCslFileGrasp92(filename::String)
        coreSubshells = Subshell[];    peelSubshells =  Subshell[]

        f  = open(filename)
        sa = readline(f);   sa[1:14] != "Core subshells"   &&   error("Not a Grasp92 .cls file.")
        sa = readline(f);   if  length(sa) > 3    go = true    else  go = false    end
        i  = -1
        while go
            i = i + 1;    sh = Basics.subshellGrasp( strip(sa[5i+1:5i+5]) );    push!(coreSubshells, sh)
            if  5i + 6 >  length(sa)    break    end
        end
        sa = readline(f);   println("sa = $sa")
        sa = readline(f);   if  length(sa) > 3    go = true    else  go = false    end
        i  = -1
        while go
            i = i + 1;    sh = Basics.subshellGrasp( strip(sa[5i+1:5i+5]) );    push!(peelSubshells, sh)
            if  5i + 6 >  length(sa)    break    end
        end
        # Prepare a list of all subshells to define the common basis
        subshells = deepcopy(coreSubshells);    append!(subshells, peelSubshells)

        # Now, read the CSF in turn until the EOF
        readline(f)
        NoCSF = 0;    csfs = CsfR[]
        Defaults.setDefaults("relativistic subshell list", subshells)

        while true
            sa = readline(f);      if length(sa) == 0   break    end
            sb = readline(f);      sc = readline(f)
            NoCSF = NoCSF + 1
            push!(csfs, ManyElectron.CsfRGrasp92(subshells, coreSubshells, sa, sb, sc) )
        end

        println("  ... $NoCSF CSF read in from Grasp92 file  $filename)") 
        NoElectrons = sum( csfs[1].occupation )

        Basis(true, NoElectrons, subshells, csfs, coreSubshells, Dict{Subshell,Radial.Orbital}() )
    end


    """
    `Basics.readOrbitalFileGrasp92(filename::String, grid.Radial.Grid)`  
        ... reads in the orbitals list from a Grasp92 .rwf or .out file; a dictionary 
            orbitals::Dict{Subshell,Radial.Orbital} is returned.
    """
    function Basics.readOrbitalFileGrasp92(filename::String, grid::Radial.Grid)
        # using FortranFiles
        f = FortranFile(filename)
        
        orbitals = Dict{Subshell,Radial.Orbital}();    first = true

        if (String(Base.read(f, FString{6})) != "G92RWF")
        close(f)
        error("File \"", filename, "\" is not a proper Grasp92 RWF file")
        end

        Defaults.setDefaults("standard grid", grid);   

        while true
        try
            n, kappa, energy, mtp = Base.read(f, Int32, Int32, Float64, Int32)
            pz, P, Q = Base.read(f, Float64, (Float64, mtp), (Float64, mtp))
            ra = Base.read(f, (Float64, mtp))
                # Now place the data into the right fields of Orbital
                subshell = Subshell( Int64(n), Int64(kappa) )   
                if  energy < 0   isBound = true    else    isBound = false   end 
                useStandardGrid = true
                println("Basics.readOrbitalFileGrasp92-aa: warning ... Radial functions still need to be interpolated upon the given grid.")
                Px = P
                Qx = Q 

                orbitals = Base.merge( orbitals, Dict( subshell => Orbital(subshell, isBound, useStandardGrid, energy, Px, Qx, Radial.Grid() ) ))
        catch ex
            if     ex isa EOFError  break
            else   throw(ex)
            end
        end
        end
        close(f)

        return( orbitals )
    end


    """
    `Basics.readMixFileRelci(filename::String)`  
        ... reads in the mixing coefficients from a RELCI .mix file; a multiplet::Multiplet is returned.
    """
    function Basics.readMixFileRelci(filename::String, basis::Basis)
        levels = Level[];    name =  "from Relci"

        f  = open(filename)
        sa = readline(f);   sa[1:14] != "G92MIX (format"   &&   error("Not a formatted relci.mix file.")
        sa = readline(f);   NoElectrons = Base.parse(Int64, strip(sa[1:6]) );    NoCsf = Base.parse(Int64, strip(sa[7:12]) ) 
                            NoSubshells = Base.parse(Int64, strip(sa[13:18]) ) 
        sa = readline(f);   NoLevels    = Base.parse(Int64, strip(sa[1:6]) )
        #
        sa = readline(f)
        sa = readline(f);   list2J = AngularJ64[AngularJ64(0)  for  i = 1:NoLevels];    listp = Basics.Parity[Basics.plus  for  i = 1:NoLevels]
        for  i = 1:NoLevels
            jj = Base.parse(Int64, strip(sa[(i-1)*6+1:(i-1)*6+5]) );   list2J[i] = AngularJ64(jj//2)
            pp = string( sa[(i-1)*6+6] );                              listp[i]  = Basics.Parity( pp )
        end
        #
        sa   = readline(f);   listEnergy = Float64[0.  for  i = 1:NoLevels]
        enav = Base.parse(Float64, strip(sa[1:26]) )
        for  i = 1:NoLevels
            en = Base.parse(Float64, strip(sa[26i+1:26i+26]) );   listEnergy[i] = en
        end
        #
        eigenvectors = zeros(NoCsf,NoLevels)
        for  j = 1:NoCsf
            sa   = readline(f)
            for  i = 1:NoLevels
            cm = Base.parse(Float64, strip(sa[16*(i-1)+1:16(i-1)+16]) );   eigenvectors[j,i] = cm
            end
        end
        #
        for  j = 1:NoLevels
            mc = Float64[0.  for  i = 1:NoCsf]
            for  i = 1:NoCsf    mc[i] = eigenvectors[i,j]    end
            Jx = list2J[j];   Mx = AngularM64( Jx.num//Jx.den )
            push!( levels, Level(Jx, Mx, listp[j], j, enav+listEnergy[j], 0., true, basis, mc) )
        end

        multiplet = Multiplet(name, levels)
        return( multiplet )
    end


    """
    `Basics.readMixingFileGrasp92(filename::String, basis::Basis)`  
        ... reads in the mixing coefficients from a Grasp92 .mix file by using the given basis; A multiplet::Multiplet 
            is returned.
    """
    function Basics.readMixingFileGrasp92(filename::String, basis::Basis)
        name = "from Grasp"
        f = FortranFile(filename)

        if (String(Base.read(f, FString{6})) != "G92MIX")
            close(f)
        error("File \"", filename, "\" does not seem to be a Grasp92 MIX file")
        end

        no_electrons, no_csf, no_subshells = read(f, Int32, Int32, Int32)
        no_asf = Base.read(f, Int32)
        level_numbers = Base.read(f, (Int32, no_asf))
        ijip = Base.read(f, (Int32, 2, no_asf))

        mjs = ijip[1,:] # 2j+1
        mparities = ijip[2,:] # 1=='+', -1=='-'

        mixfile_to_j(i::Int32) = J((i-1)//2)
        mixfile_to_p(i::Int32) = i==1 ?  Basics.plus : Basics.minus

        js = mixfile_to_j.(mjs)
        parities = mixfile_to_p.(mparities)

        avg_energy, eval = read(f, Float64, (Float64, no_asf))
        energies = avg_energy + eval
        eigenvectors = read(f, (Float64, no_csf, no_asf))

        ##x @show no_electrons
        ##x @show no_csf
        ##x @show no_subshells
        ##x @show no_asf
        ##x println("angular momenta = ", js)
        ##x @show parities
        ##x @show energies
        ##x println("first eigenvector = ", eigenvectors[:,1])

        close(f)
        multiplet = Multiplet(name, Level[] )
        return( multiplet )
    end
    


    """
    `Basics.recast()`  
        ... recasts some data from one number/representation into another one; cf. the supported keystrings and return values.

    + `("rate: radiative, to decay width", value::Float64)`  
        ... to recast a given radiative rate (Einstein A in atomic units) into a decay withs, taking the selected energy unit 
            into account. A Float64 is returned.

    + `("rate: radiative, to Einstein A", value::Float64)`  
        ... to recast a given spontaneous radiative rate (= Einstein A-coefficient), taking the selected unit into account. 
            A Float64 is returned.

    + `("rate: radiative, to Einstein B", value::Float64)`  
        ... to recast a given radiative rate (Einstein A in atomic units) into an Einstein B-coefficient, taking the selected 
            unit into account. A Float64 is returned.

    + `("rate: radiative, to g_f", value::Float64)`  
        ... to recast a given radiative rate (Einstein A in atomic units) into an oscillator strength g_f; 
            a Float64 is returned.
    """
    function Basics.recast(sa::String, line::Union{Einstein.Line,PhotoEmission.Line}, wa::Float64)
        ##x global  CONVERT_ENERGY_AU_TO_EV, CONVERT_ENERGY_AU_TO_KAYSERS, CONVERT_ENERGY_AU_TO_PER_SEC, CONVERT_TIME_AU_TO_SEC,
        ##x         CONVERT_CROSS_SECTION_AU_TO_BARN,  CONVERT_RATE_AU_TO_PER_SEC 
    
        if       sa == "rate: radiative, to decay width"
            width = Defaults.convertUnits("energy: from atomic", wa)
            return( width )

        elseif   sa == "rate: radiative, to Einstein A"
            einsteinA = Defaults.convertUnits("rate: from atomic", wa)
            return( einsteinA )

        elseif   sa == "rate: radiative, to Einstein B"
            einsteinB = (AngularMomentum.twoJ(line.initialLevel.J) + 1) / (AngularMomentum.twoJ(line.finalLevel.J) + 1) * pi^2 *
                        Defaults.getDefaults("speed of light: c")^3 / line.omega^3  * wa
            einsteinB = Defaults.convertUnits("Einstein B: from atomic", einsteinB)
            return( einsteinB )

        elseif   sa == "rate: radiative, to g_f"
            gf = (AngularMomentum.twoJ(line.initialLevel.J) + 1) / (AngularMomentum.twoJ(line.finalLevel.J) + 1) / 2. * 
                Defaults.getDefaults("speed of light: c")^3 / line.omega^2 * wa     
            return( gf )

        else     error("Unsupported keystring = $sa")
        end
    end
    
    

    """
    `Basics.sortByEnergy(multiplet::Multiplet)`  
        ... to sort all levels in the multiplet into a sequence of increasing energy; a (new) multiplet::Multiplet is returned.
    """
    function Basics.sortByEnergy(multiplet::Multiplet)
        sortedLevels = Base.sort( multiplet.levels , lt=Base.isless)
        newLevels = Level[];   index = 0
        for lev in sortedLevels
            index = index + 1
            push!(newLevels, Level(lev.J, lev.M, lev.parity, index, lev.energy, lev.relativeOcc, lev.hasStateRep, lev.basis, lev.mc) )
        end
        
        newMultiplet = Multiplet(multiplet.name, newLevels)
        
        return( newMultiplet )  
    end

        
    """
    `Basics.sortByEnergy(multiplet::Hfs.IJF_Multiplet)`  
        ... to sort all hyperfine levels in the multiplet into a sequence of increasing energy; a (new) multiplet::Hfs.IJF_Multiplet 
            is returned.
    """
    function Basics.sortByEnergy(multiplet::Hfs.IJF_Multiplet)
        sortedLevels = Base.sort( multiplet.levelFs , lt=Base.isless)
        newLevels = Hfs.IJF_Level[];   index = 0
        for lev in sortedLevels
            index = index + 1
            push!(newLevels, Hfs.IJF_Level(lev.I, lev.F, lev.M, lev.parity, lev.energy, lev.hasStateRep, lev.basis, lev.mc) )
        end
        
        newMultiplet = Hfs.IJF_Multiplet(multiplet.name, newLevels)
        
        return( newMultiplet )  
    end


    """
    `Basics.tabulate(stream::IO, sa::String, multiplet::Multiplet, levelNos::Array{Int64,1})`  
        ... tabulates the energies from the multiplet with level numbers in levelNos due to different criteria.

    + `(stream::IO, "multiplet: energies", multiplet::Multiplet)`  
        ... to tabulate the energies of all levels of the given multiplet into a neat format; nothing is returned.
    + `(stream::IO, "multiplet: energy relative to immediately lower level", multiplet::Multiplet, levelNos::Array{Int64,1})`  
        ... to tabulate the energy splitting between neighboured levels of all levels of the given multiplet into a neat
            format; nothing is returned.
    + `(stream::IO, "multiplet: energy of each level relative to lowest level", multiplet::Multiplet, levelNos::Array{Int64,1}`  
        ... to tabulate the energy splitting of all levels with regard to the lowest level of the given multiplet into 
            a neat format; nothing is returned.
    """
    function Basics.tabulate(stream::IO, sa::String, multiplet::Multiplet, levelNos::Array{Int64,1})
        if        sa == "multiplet: energies"
            println(stream, "\n  Eigenenergies:")
            sb = "  Level  J Parity          Hartrees       " * "             eV                   " *  TableStrings.inUnits("energy")     
            println(stream, "\n", sb, "\n")
            for  i = 1:length(multiplet.levels)
                if  i in levelNos
                    lev = multiplet.levels[i]
                    en  = lev.energy;    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                    sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                    @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
                end
            end

        elseif    sa == "multiplet: energy relative to immediately lower level"
            println(stream, "\n  Energy of each level relative to immediately lower level:")
            sb = "  Level  J Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")  
            println(stream, "\n", sb, "\n")
            for  i = 2:length(multiplet.levels)
                if  i in levelNos
                    lev = multiplet.levels[i]
                    en    = lev.energy - multiplet.levels[i-1].energy;    
                    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                    sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                    @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
                end
            end

        elseif    sa == "multiplet: energy of each level relative to lowest level"
            println(stream, "\n  Energy of each level relative to lowest level:")
            sb = "  Level  J Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")      
            println(stream, "\n", sb, "\n")
            for  i = 2:length(multiplet.levels)
                if  i in levelNos
                    lev = multiplet.levels[i]
                    en    = lev.energy - multiplet.levels[1].energy;    
                    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                    sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                    @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
                end
            end
        else
            error("Unsupported keystring.")
        end

        return( nothing )  
    end



    """
    `Basics.tabulate(sa::String, multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)`  
        ... tabulates the energies from the multiplet due to different criteria.

    + `("multiplet: energies", multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)`  
        ... to tabulate the energies of all hyperfine levels of the given multiplet into a neat format; nothing is returned.
    + `("multiplet: energy of each level relative to lowest level", multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)`  
        ... to tabulate the energy splitting of all levels with regard to the lowest level of the given multiplet into 
            a neat format; nothing is returned.
    """
    function Basics.tabulate(sa::String, multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)
        if        sa == "multiplet: energies"
            println(stream, "\n  Eigenenergies for nuclear spin I = $(multiplet.levelFs[1].I):")
            sb = "  Level  F Parity          Hartrees       " * "             eV                   " *  TableStrings.inUnits("energy")     
            println(stream, "\n", sb, "\n")
            for  i = 1:length(multiplet.levelFs)
                lev = multiplet.levelFs[i]
                en  = lev.energy;    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
            end

        elseif    sa == "multiplet: energy of each level relative to lowest level"
            println(stream, "\n  Energy of each level relative to lowest level for nuclear spin I = $(multiplet.levelFs[1].I):")
            sb = "  Level  F Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")      
            println(stream, "\n", sb, "\n")
            for  i = 2:length(multiplet.levelFs)
                lev = multiplet.levelFs[i]
                en    = lev.energy - multiplet.levelFs[1].energy;    
                en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity))  * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
            end
        else
            error("Unsupported keystring.")
        end

        return( nothing )  
    end

    
    """
    `Basics.tabulateKappaSymmetryEnergiesDirac(kappa::Int64, evalues::Array{Float64,1}, ns::Int64, nuclearModel::Nuclear.Model)`  
        ... tabulates the eigenenergies for a given symmetry block kappa together with the corresponding Dirac energies for a 
            point-like nucleus. The index ns tells the number of 'negative-continnum' energies in the given evalues. nothing is 
            returned.
    """
    function Basics.tabulateKappaSymmetryEnergiesDirac(kappa::Int64, evalues::Array{Float64,1}, ns::Int64, nuclearModel::Nuclear.Model)
        Z = nuclearModel.Z
        # Determine the allowed principal quantum numbers n
        l = Basics.subshell_l(Subshell(101,kappa))
        println("  ", TableStrings.hLine(77))
        sa = "  "
        sa = sa * TableStrings.center( 7, "Index";            na=2)
        sa = sa * TableStrings.center(10, "Subshell";         na=3)
        sa = sa * TableStrings.center(17, "Energies [a.u.]";  na=2)
        sa = sa * TableStrings.center(17, "Dirac-E  [a.u.]";  na=2)
        sa = sa * TableStrings.center(17, "Delta-E / |E|";    na=2)
        println(sa)
        println("  ", TableStrings.hLine(77))
        for  i = ns+1:ns+8
            sa = " " * TableStrings.center( 6, TableStrings.level(i-ns); na=2)
            sa = sa *  TableStrings.flushright(10, string(Subshell(i-ns+l, kappa)); na=6)
            en = Basics.computeDiracEnergy(Subshell(i-ns+l, kappa), Z)
            if  evalues[i] >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", evalues[i])                       * "    "
            if  en         >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", en)                               * "    "
            if  evalues[i]-en >= 0.   sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", (evalues[i]-en)/abs(evalues[i])) * "    "
            println(sa)
        end
        println("      :       :    ")
        for  i = length(evalues)-5:length(evalues)
            sa = " " * TableStrings.center( 6, TableStrings.level(i-ns); na=2)
            sa = sa *  TableStrings.flushright(10, string(Subshell(i-ns+l, kappa)); na=6)
            en = Basics.computeDiracEnergy(Subshell(i-ns+l, kappa), Z)
            if  evalues[i] >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", evalues[i])                       * "    "
            if  en         >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", en)                               * "    "
            if  evalues[i]-en >= 0.   sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", (evalues[i]-en)/abs(evalues[i])) * "    "
            println(sa)
        end
        println("  ", TableStrings.hLine(77))

        return( nothing )
    end

    """
    `Basics.tools()`  ... select different tools from a menu for which no further input is required.
    """
    function Basics.tools(x::Int64)
        println("Tools.select(): A menu if nothing comes in")
        #
        t1 = "Simple tools JAC without input: "
        b1 = dropdown(["Start task", "Grid calculator", "Another task"])
        b2 = slider(1:100, label = "To what extent?", value = 33)
        update = button("Update")
        ui = vbox( hbox( pad(1em, t1) ),
                hbox( pad(1em, b1), pad(1em, b2), pad(1em, update) )
                )
        Interact.display(ui)
        output = Interact.@map (&update;  (observe(b1)[], observe(b2)[]) ) 
        return( output )
    end


    """
    `Basics.tools(dict::Dict)`  
        ... select different tools from a menu if proper results::dict are given in terms of a dictionary as obtained
            usually from an Atomic.Computation.
    """
    function Basics.tools(dict::Dict)
        println("Tools.select(): A menu if a dict comes in; $dict ")
    end


    """
    `Basics.yesno(question::String, sa::String)`  
        ... Returns true if the answer 'yes' or 'y' is found, and false otherwise; sa = {"Y", "N"} determines how the zero-String 
            "" is interpreted. The given question is repeated until a proper answer is obtained.
    """
    function yesno(question::String, sa::String)
        while  true
            print(question * "  ");    reply = strip( readline(STDIN) )
            if      sa == "Y"
                if      reply in ["", "Y", "Yes", "y", "yes"]   return( true )
                elseif  reply in [    "N", "No",  "n", "no" ]   return( false )
                        println("... answer not recognized ... redo:")
                end
            elseif  sa == "N"
                if      reply in [    "Y", "Yes", "y", "yes"]   return( true )
                elseif  reply in ["", "N", "No",  "n", "no" ]   return( false )
                        println("... answer not recognized ... redo:")
                end
            end
        end
    end


end # module
