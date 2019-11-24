
"""
`module  JAC.BasicsAG`  
    ... a submodel of JAC that contains methods that support tasks related to spectroscopic computation.
"""
module BasicsAG

    using Printf,  LinearAlgebra, ..AngularMomentum, ..Atomic, ..Basics, ..Continuum, ..Defaults, ..Einstein, ..LSjj, ..ManyElectron, 
                  ..PhotoEmission, ..Radial
    
    export recastAG

    """
    `Basics.add(pota::Radial.Potential, potb::Radial.Potential)`  
        ... to add two radial potentials together that are defined on the same grid. A potential::RadialPotential is returned 
            that inherits its radial size from the potential that is defined in a larger range of r-values.
    """
    function  Basics.add(pota::Radial.Potential, potb::Radial.Potential)
        if  pota.grid.NoPoints != potb.grid.NoPoints  ||   pota.grid.rnt != potb.grid.rnt  ||  pota.grid.h != potb.grid.h  ||  
            pota.grid.hp != potb.grid.hp                   error("stop a")   end
        name = pota.name * "+" * potb.name;   nmax = max(length(pota.Zr), length(potb.Zr));   nmin = min(length(pota.Zr), length(potb.Zr))
        Zr   = zeros(nmax);   
        nx = length(pota.Zr);    Zr[1:nx] = Zr[1:nx] + pota.Zr[1:nx] 
        nx = length(potb.Zr);    Zr[1:nx] = Zr[1:nx] + potb.Zr[1:nx]
        Zr[nmin+1:nmax] .= Zr[nmin] 

        potential = Radial.Potential(name, Zr, pota.grid)
        return( potential )
    end



    """
    `Basics.addZerosToCsfR(nz::Int64, csf::CsfR)`  
        ... 'adds' a number of zeros to the occupation, seniority, etc. of a given CsfR.
    """
    function  Basics.addZerosToCsfR(nz::Int64, csf::CsfR)
        !csf.useStandardSubshells   &&   error("Zeros can only be added to CSF with standard subshell order.")
        
        println("nz = $nz")
        occupation = csf.occupation;   seniority = csf.seniority;    subshellJ = csf.subshellJ;    subshellX = csf.subshellX 
        
        for  i = 1:nz
            push!(occupation, 0);    push!(seniority, 0);    push!(subshellJ, AngularJ64(0));    push!(subshellX, subshellJ[end])
        end
        
        newCsf = CsfR( csf.useStandardSubshells, csf.J, csf.parity, occupation, seniority, subshellJ, subshellX, Subshell[])
        return(newCsf)
    end

    
    
    """
    `Basics.analyze("level decomposition: % of NR configurations", level::Level)`  
        ... to analyze and list the non-relativistic configurations with a weight larger than 5%. A list of NR configurations with their 
            corresponding weights are printed, but nothing is returned.
    """
    function Basics.analyze(sa::String, level::Level)
        !(sa == "level decomposition: % of NR configurations")   &&   error("Unsupported keystring = $sa")
        !(level.hasBasis)                                        &&   error("Levels without a basis cannot be analyzed.")
        
        confList   = generate("configuration list: NR, from basis", level.basis)
        percentage = zeros( length(confList) )
        
        # Now generate for each CSF in this basis a Configuration, compare with confList and add contribution to the corresponding percentage
        for  k = 1:length(level.basis.csfs)
            csf = level.basis.csfs[k]
            shellList = Dict{Shell,Int64}[]
            for  i = 1:length(csf.occupation)
                n = csf.subshell[i].n;    l = Basics.subshell_l( csf.subshell[i] );    sh = Shell(n, l)
                if     haskey(shellList, sh)    shellList(sh) = shellList(sh) + csf.occupation[i]
                else   shellList = merge( shellList, Dict(sh => csf.occupation[i]) )
                end
            end
            confNew = Configuration( shellList, NoElectrons )
            # Compare and add the weight to the right configuration
            for  i = 1:length(confList)
                if   confList[i] == confNew    percentage[i] = percentage[i] + level.mc[k]^2;    break    end 
            end
        end
    
        # Check that all weights sum up properly
        wa = sum( percentage );   abs(1. - wa) > 1.0e-3  &&    error("Total percentage $wa must add to 1.")
        # Sort the percentage and print both, configurations and percentage
        wb = sortperm( percentage )
        for i in wb
            println("Configuration: $(confList[i])  ... with  $(percentage[i]) %")
        end

        return( nothing )
    end



    """
    `Basics.analyze("level decomposition: % of jj-coupled CSF", level::Level, N::Int64)`  
        ... to anaylze and list (up to) N relativistic CSF, together with their weight |c_n|^2 in the expansion of the given level. 
            A list of CSF and their corresponding weights are printed, but nothing is returned otherwise.  **Not yet implemented !**
    """
    function Basics.analyze(sa::String, level::Level, N::Int64)
        !(sa == "level decomposition: % of jj-coupled CSF")   &&   error("Unsupported keystring = $sa")
        error("Not yet implemented !")

        return( nothing )
    end


        
    """
    `Basics.analyzeConvergence(blocka::Basics.Eigen, blockb::Basics.Eigen)`  
        ... to analyze the convergence for two blocks of eigenvalues und eigenvectors; the value > 0 indicates the relative maximum 
            difference abs( a-b/a+b ) between two eigenvalues. The value is set arbitrarely set to 1.0e+2 if the number of eigenvalues
            differ in the two blocks.
    """
    function Basics.analyzeConvergence(blocka::Basics.Eigen, blockb::Basics.Eigen)
        va = blocka.values;   vb = blockb.values
        if       length(va) !=  length(vb)    wa = 1.0e+2   
        else     wa = 0.;       for  i = 1:length(va)   wa = max(wa, abs( (va[i]-vb[i]) / (va[i]+vb[i]) ))   end
        end
        return( wa )
    end


        
    """
    `Basics.analyzeConvergence(a::Radial.Orbital, b::Radial.Orbital)`  
        ... to analyze the convergence for two radial orbitals; the value > 0 indicates the relative maximum difference abs( a-b/a+b ) 
            between the energies or (large or small) components of the orbitals. The function assumes that the orbitals are represented 
            on the same grid and that the components are zero if only one of the orbitals is defined beyond some maximum grid point.
    """
    function Basics.analyzeConvergence(a::Radial.Orbital, b::Radial.Orbital)
        if       a.subshell !=  b.subshell                                             error("stop a")   
        elseif   a.useStandardGrid !=  b.useStandardGrid   ||  !(a.useStandardGrid)    error("stop b")   
        else     wa = 0.;   wa = max(wa, abs( (a.energy-b.energy) / (a.energy+b.energy) ));   nx = min(length(a.P), length(b.P))
                #=== for  i = 1:nx 
                    wp = abs(a.P[i])+abs(b.P[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.P[i]-b.P[i]) / wp ))   end 
                    ##x if  abs( (a.P[i]-b.P[i]) / wp ) > 1.01    Defaults.warn(AddWarning, "analyzeConv: P, i = $i Pa=$(a.P[i]) Pb=$(b.P[i])")   end
                    wp = abs(a.Q[i])+abs(b.Q[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.Q[i]-b.Q[i]) / wp ))   end   
                    ##x if  abs( (a.Q[i]-b.Q[i]) / wp ) > 1.01    Defaults.warn(AddWarning, "analyzeConv: Q, i = $i Qa=$(a.Q[i]) Qb=$(b.Q[i]) ")   end
                end ===#
        end
        return( wa )
    end

        

    """
    `Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)`  
        ... to determine the NoEnergySharings sharings by using the Gauss-Legendre integration points in the interval
            [0., energy]; an Array{Tuple{Float64,Float64},1} is returned. This methods assumes that the package 
            GaussQuadrature is 'used'.
    """
    function Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)
        xx, ww = GaussQuadrature.legendre(NoEnergySharings)
        sharings = Tuple{Float64,Float64}[]
        for  x in xx
            xrel = (x+1.) / 2
            push!(sharings, (xrel*energy, energy - xrel*energy) )
        end
        return( sharings )
    end



    """
    `Basics.determineHoleShells(conf::Configuration)`  
        ... to determine the hole-shells of a given non-relativistic configuration; a shellList::Array{Shell,1} is returned.
    """
    function Basics.determineHoleShells(conf::Configuration)

        shellList = Shell[]
        for  (k,v) in conf.shells
        if   conf.shells[k] < 2*(2k.l+1)   push!(shellList, k)    end 
        end

        return( shellList )
    end



    """
    `Basics.determineParity(conf::Configuration)`  ... to determine the parity of a given non-relativistic configuration.
    """
    function Basics.determineParity(conf::Configuration)

        par = 1
        for  (k,v) in conf.shells
        if   iseven(k.l)    p = 1   else   p = -1    end    
        par = par * (p^v)
        end

        if       par == 1    return( Basics.plus )  
        elseif   par == -1   return( Basics.minus )
        else     error("stop b")
        end  
    end


    """
    `Basics.determineParity(conf::ConfigurationR)`  ... to determine the parity of a given relativistic configuration.
    """
    function Basics.determineParity(conf::ConfigurationR)

        par = 1
        for  (k,v) in conf.subshells
        l = Basics.subshell_l(k)
        if   iseven(l)    p = 1   else   p = -1    end    
        par = par * (p^v)
        end

        if       par == 1    return( Basics.plus )  
        elseif   par == -1   return( Basics.minus )
        else     error("stop b")
        end  
    end


    """
    `Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)`  
        ... to determine the specified lines as tupels of initial- and final levels. A level index is 0 in lineList always refers to
            all levels from the corresponding multiplet. A (unique) newList::Array{Tuple{Int64,Int64},1} is returned.
    """
    function Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)
        newList = Tuple{Int64,Int64}[]
        for tupl in lineList
            if      tupl[1] >  length(initialMultiplet.levels)   error("Initial level index $(tupl[1]) must be <= No. initial levels.")
            elseif  tupl[2] >  length(finalMultiplet.levels)     error("Final level index $(tupl[2]) must be <= No. final levels.")
            elseif  tupl[1] == 0  &&   tupl[2] == 0    error("Tupel (0,0) is not allowed.")
            elseif  tupl[1] == 0       for i = 1:length(initialMultiplet.levels)   push!(newList, (i, tupl[2]))   end
            elseif  tupl[2] == 0       for i = 1:length(finalMultiplet.levels)     push!(newList, (tupl[1], i))   end
            else                                                                   push!(newList, (tupl[1], tupl[2]))
            end
        end
        newList = unique(newList)
        return( newList)
    end     


    """
    `Basics.determineSelectedPathways(pathwayList::Array{Tuple{Int64,Int64,Int64},1}, initialMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                      finalMultiplet::Multiplet)`  
        ... to determine the specified pathways as tupels of initial-, intermediate- and final levels. A level index 0 in pathwayList always 
            refers to all levels from the corresponding multiplet. A (unique) newList::Array{Tuple{Int64,Int64,Int64},1} is returned.
    """
    function Basics.determineSelectedPathways(pathwayList::Array{Tuple{Int64,Int64,Int64},1}, initialMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                    finalMultiplet::Multiplet)
        newList = Tuple{Int64,Int64,Int64}[]
        for tupl in pathwayList
            if      tupl[1] >  length(initialMultiplet.levels)       error("Initial level index $(tupl[1]) must be <= No. initial levels.")
            elseif  tupl[2] >  length(intermediateMultiplet.levels)  error("Intermediate level index $(tupl[2]) must be <= No. intermediate levels.")
            elseif  tupl[3] >  length(finalMultiplet.levels)         error("Final level index $(tupl[3]) must be <= No. final levels.")
            elseif  (tupl[1] + tupl[2] == 0  || tupl[1] + tupl[3] == 0  ||  tupl[2] + tupl[3] == 0)   &&  error("Only tupels with one 0-index are allowed.")
            elseif  tupl[1] == 0       for i = 1:length(initialMultiplet.levels)          push!(newList, (i, tupl[2], tupl[3]))   end
            elseif  tupl[2] == 0       for i = 1:length(intermediateMultiplet.levels)     push!(newList, (tupl[1], i, tupl[3]))   end
            elseif  tupl[3] == 0       for i = 1:length(finalMultiplet.levels)            push!(newList, (tupl[1], tupl[2], i))   end
            else                                                                          push!(newList, (tupl[1], tupl[2], tupl[3]))
            end
        end
        newList = unique(newList)
        return( newList)
    end     

    
    
    """
    `Basics.diagonalize("matrix: Julia, eigfact", matrix::Array{Float64,2})`  
        ... to apply the standard eigfact() method from Julia for a quadratic, full matrix; nothing about the symmetry of 
            the matrix is assumed here; an eigen::Basics.Eigen is returned.
    """
    function Basics.diagonalize(sa::String, matrix::Array{Float64,2})
        if       sa == "matrix: Julia, eigfact" 
            # Use the standard eigfact() method from Julia for a quadratic, full matrix   
            wa = eigen( matrix )
            ##x println("diagonalize-ab: values  = ", wa[:values] )
            ##x println("diagonalize-ac: vectors = ", wa[:vectors] )
            vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
            for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d])    end
            wc = Basics.Eigen( wa.values, vectors )
            return( wc )
        elseif   sa == "matrix: Julia method B"     error("Not yet implemented")
        else     error("Unsupported keystring = $sa")
        end
    end


    """
    `Basics.diagonalize("generalized eigenvalues: Julia, eigfact", matrixA::Array{Float64,2}, matrixB::Array{Float64,2})`  
        ... to apply the standard eigfact() method from Julia for a generalized eigenvalue problem with two quadratic   
            (full) matrices; nothing about the symmetry of the matrix is assumed here; an eigen::JAC.Eigen is returned.
    """
    function Basics.diagonalize(sa::String, matrixA::Array{Float64,2}, matrixB::Array{Float64,2})
        if       sa == "generalized eigenvalues: Julia, eigfact" 
            # Use the standard eigfact() method from Julia for two quadratic, full matrices   
            wa = eigen( matrixA, matrixB )
            ##x println("diagonalize-ab: values  = ", wa[:values] )
            ##x println("diagonalize-ac: vectors = ", wa[:vectors] )
            ##x vectors = Vector{Float64}[];    wb = wa[:vectors];    d = size(wb)[1]
            vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
            # for  i = 0:d-1    push!(vectors, real(wb[i*d+1:i*d+d]) )    end
            for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d] )    end
            ##x wc = Basics.Eigen( wa[:values], vectors )
            # wc = Basics.Eigen( real(wa.values), vectors )
            wc = Basics.Eigen( wa.values, vectors )
            return( wc )
        else     error("Unsupported keystring = $sa")
        end
    end

    """
    `Basics.display("constants")`  or  `("physical constants")`  
        ... to display (all) currently defined physical constants; nothing is returned if not indicated otherwise. 
            Cf. Defaults.setDefaults().

    `Basics.display("settings")`     ... to display (all) currently defined settings of the JAC module.
    """
    function Basics.display(sa::String)

        if        sa in ["constants", "physical constants"]
            println("Physical constants are defines as follows:  \n", 
                    "------------------------------------------  \n")
            sb = "  + Fine-structure constant:    " * string( Defaults.getDefaults("alpha") );                println(sb)        
            sb = "  + Electron mass [kg]:         " * string( Defaults.getDefaults("electron mass: kg") );    println(sb)        
            sb = "  + Electron mass [amu]:        " * string( Defaults.getDefaults("electron mass: amu") );   println(sb)        
            println()

        elseif   sa == "settings"
            println("Current settings of the JAC module:  \n", 
                    "-----------------------------------  \n")
            sb = "  + Framework:                              " * string( Defaults.getDefaults("framework") );                println(sb)        
            sb = "  + Energy unit:                            " * string( Defaults.getDefaults("unit: energy") );             println(sb)        
            sb = "  + Rate and transition probability unit:   " * string( Defaults.getDefaults("unit: rate") );               println(sb)        
            sb = "  + Cross section unit:                     " * string( Defaults.getDefaults("unit: cross section") );      println(sb) 
            sb = "  + Time unit:                              " * string( Defaults.getDefaults("unit: time") );               println(sb) 
            println()       
            
            if      Defaults.getDefaults("standard grid") != false    println("  + A standard grid has been defined; cf. Defaults.getDefaults()" )
            elseif  Defaults.getDefaults("standard grid") == false    println("  + No standard grid has yet been defined; cf. Defaults.setDefaults()" )
            end
            println()       
    
        else    error("Unsupported keystring:: $sa")
        end

        return( nothing )
    end


    """
    `Basics.excludeDoubles(confList::Array{Configuration,1})`  
        ... to exlude from the (non-relativistic) confList all 'doubles', ie. configurations that occured before in the list; 
            a new confList::Array{Configuration,1} is returned.
    """
    function Basics.excludeDoubles(confList::Array{Configuration,1})
        confListNew = [ deepcopy(confList[1]) ]

        for  ic  in 2:length(confList)
            add = true
            for  id in 1:ic-1    if   confList[ic] == confList[id]    add = false;   break   end    end
            if  add   push!( confListNew, confList[ic])   end
        end    
        
        return( confListNew )
    end
        
        
    """
    `Basics.excludeDoubles(csfList::Array{CsfR,1})`  
        ... to exlude from the (relativistic) csfList all 'doubles', ie. CSF that occured already before in the list; 
            a new csfList::Array{CsfR,1} is returned.
    """
    function Basics.excludeDoubles(confList::Array{CsfR,1})
        error("Not yet implemented !")
        
        csfListNew = [ deepcopy(csfList[1]) ]
        
        for  ic  in 2:length(csfList)
            add = true
            for  id in 1:ic-1
                if   csfList[ic] == csfList[id]    add = false;   break   end
            end
            if  add   push!( csfListNew, csfList[ic])	end
        end    
        return( csfListNew )
    end
    
        
    """
    `Basics.extractNoOpenShells(conf::Configuration)`  
        ... determine the number of open (nonrelativistic) shells in conf; a singleton of type AbstractOpenShell is returned.
    """
    function Basics.extractNoOpenShells(conf::Configuration)
        ##x println("conf = $conf")
        ns = 0
        for  (shell, occ)  in conf.shells
            if      occ == 0  ||  occ == 2*(2*shell.l + 1)
            else    ns = ns + 1
            end
        end
        ##x println("conf = $conf  ns = $ns")
        
        if      ns == 0   return ( LSjj.ZeroOpenShell() )
        elseif  ns == 1   return ( LSjj.OneOpenShell() )
        elseif  ns == 2   return ( LSjj.TwoOpenShells() )
        elseif  ns == 3   return ( LSjj.ThreeOpenShells() )
        else    error("stop a")
        end
        
    end


    """
    `Basics.extractNonrelativisticShellList(subshells::Array{Subshell,1})  
        ... extract a nonrelativistic shell list from the given subshell list; the shells appear in the same order
            as the corresponding subshells. A shellList::Array{Subshell,1} is returned.
    """
    function Basics.extractNonrelativisticShellList(subshells::Array{Subshell,1}) 
        shellList = Shell[]

        for  subsh  in subshells
            l = Basics.subshell_l(subsh);   shell = Shell( subsh.n, l)
            if  shell in shellList      else  push!(shellList, shell)   end
        end 
        ##x println("Basics.extractNonrelativisticShellList()::  shellList = $shellList")
        
        return( shellList )
    end


    """
    `Basics.extractNonrelativisticConfigurations(basis::Basis)  
        ... extract all nonrelativistic configurations that contribute to the given set of CSF in basis.csfs. 
            A confList::Array{Configuration,1} is returned.
    """
    function Basics.extractNonrelativisticConfigurations(basis::Basis)
        confList  = Configuration[]
        subshells = basis.subshells
        
        for  csf  in basis.csfs
            newShells = Dict{Shell,Int64}();    NoElectrons = 0
            for  i = 1:length(subshells)
                n = subshells[i].n;    l = Basics.subshell_l(subshells[i]);    occ = csf.occupation[i]
                shell = Shell(n,l)
                if    haskey(newShells, shell)    newShells[shell] = newShells[shell] + occ;   NoElectrons = NoElectrons + occ
                else  
                    if   occ > 0  newShells = Base.merge( newShells, Dict( shell => occ));     NoElectrons = NoElectrons + occ   end
                end
            end
            if  basis.NoElectrons != NoElectrons    error("stop a")
            else
                conf = Configuration(newShells, NoElectrons)
                if  conf in confList    continue;    else    push!( confList,  conf)    end
            end
        end 
        ##x println("Basics.extractNonrelativisticConfigurations()::  confList = $confList")
        
        return( confList )
    end



    """
    `Basics.extractNonrelativisticConfigurationFromCsfR(csf::CsfR,  basis::Basis)`  
        ... extract the nonrelativistic configuration from the given CSF that is defined in basis.csfs. 
            A conf::Configuration is returned.
    """
    function  Basics.extractNonrelativisticConfigurationFromCsfR(csf::CsfR,  basis::Basis)
        subshells = basis.subshells
        
        newShells = Dict{Shell,Int64}();    NoElectrons = 0
        for  i = 1:length(subshells)
            n = subshells[i].n;    l = Basics.subshell_l(subshells[i]);    occ = csf.occupation[i]
            shell = Shell(n,l)
            if    haskey(newShells, shell)    newShells[shell] = newShells[shell] + occ;   NoElectrons = NoElectrons + occ
            else  
                if   occ > 0  newShells = Base.merge( newShells, Dict( shell => occ));     NoElectrons = NoElectrons + occ   end
            end
        end
        
        if      basis.NoElectrons != NoElectrons    error("stop a")
        else    conf = Configuration(newShells, NoElectrons)
        end 
        ##x println("Basics.extractNonrelativisticConfigurationFromCsfR()::  conf = $conf")
        
        return( conf )
    end


    """
    `Basics.extractOpenShellQNfromCsfR(csfR::CsfR, basis::Basis)`  
        ... extracts the quantum numbers of the nonrelativistic open shells from the given csfr. Here, we assume that csfR is 
            defined in the given basis, and it is checked that the subshells of each nonrelativistic shell occur subsequently
            in basis and in the right sequence j = l - 1/2, j = l + 1/2, respectively. 
            An rQN::Dict{Shell, Tuple{Tuple{Int64,Int64,Int64,Int64,Int64}, Tuple{Int64,Int64,Int64,Int64,Int64}}}
                    = Dict(shell => (rQNm, rQNp) ) is returned with rQNm = (idxm, Nm, Qm, Jm, Jmx) and 
            rQNp = (idxp, Np, Qp, Jp, Jpx). The two tuples rQNm, rQNp contains all the subshell quantum numbers and their
            coupling in csfR.
    """
    function Basics.extractOpenShellQNfromCsfR(csfR::CsfR, basis::Basis)
        wa = Dict{Shell,Tuple{NTuple{5,Int64},NTuple{5,Int64}}}()
        for  s = 1:length(basis.subshells)
            sx = s + 1;     subsh = basis.subshells[s];     n = subsh.n;        occ   = csfR.occupation[s]
            l  = Basics.subshell_l(subsh);    j2 = Basics.subshell_2j(subsh);   shell = Shell(n,l)
            # Get l,j from next subshell ... if still defined
            if  sx <= length(basis.subshells)       subshx = basis.subshells[sx]
                nx = subshx.n;    lx = Basics.subshell_l(subshx);   j2x = Basics.subshell_2j(subshx);   occx   = csfR.occupation[sx]
            end
            # Now collect the quantum numbers into tuples
            if  l == 0   &&   0 < occ < 2
                ## wa = Base.merge( wa, Dict( shell => ( (s, csfR.occupation[s], csfR.seniority[s], 
                ##                                        Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
                ##                                        (-9, 0, -9, 0, Basics.twice(csfR.subshellX[s]) ) ) ) )
                wa = Base.merge( wa, Dict( shell => ( (-9, 0, -9, 0, 0 ),
                                                    (s, csfR.occupation[s], csfR.seniority[s], 
                                                    Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ) ) ) )
            elseif  2l + 1 == j2    continue
            elseif  2l - 1 == j2    &&    (n != nx  ||  l != lx  ||  j2x != j2 + 2)
                    error("stop a")
            elseif  0 < occ + occx < 2* (2l + 1)
                wa = Base.merge( wa, Dict(shell => ( (s, csfR.occupation[s], csfR.seniority[s], 
                                                    Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
                                                    (sx, csfR.occupation[sx], csfR.seniority[sx], 
                                                    Basics.twice(csfR.subshellJ[sx]), Basics.twice(csfR.subshellX[sx]) ) ) ) )
            end
        end
        
        ##x println("Basics.extractOpenShellQNfromCsfR()::  csfR = $csfR  \n  wa = $wa")
        
        return( wa )
    end


    """
    `Basics.extractRelativisticSubshellList(rep::Atomic.Representation)`  
        ... extract all (relativistic) subshells that (will) contribute to the given RAS expansion, and for which (start) 
            orbitals will be needed in course of the computation. A subshellList::Array{Subshell,1} is returned.
    """
    function Basics.extractRelativisticSubshellList(rep::Atomic.Representation)
        # First extract all shells, then unique this list and, finally, convert in corresponding subshells
        shellList = Shell[]
        
        for  conf in rep.refConfigs
            for  (k,v) in conf.shells   push!( shellList, k)    end
        end
        
        for step in rep.repType.steps
            append!( shellList, step.seFrom);    append!( shellList, step.seTo)    
            append!( shellList, step.deFrom);    append!( shellList, step.deTo)   
            append!( shellList, step.teFrom);    append!( shellList, step.teTo)    
            append!( shellList, step.qeFrom);    append!( shellList, step.qeTo)   
        end
        
        shellList = unique(shellList)
        ##x println("ee shellList = $shellList ")
        subshellList = Subshell[]
        
        for  sh in shellList    
            if  sh.l == 0    push!( subshellList, Subshell(sh.n, -1))
            else             push!( subshellList, Subshell(sh.n, sh.l));    push!( subshellList, Subshell(sh.n, -sh.l -1))
            end
        end
        
        ##x println("ee subshellList = $subshellList ")
        return( subshellList )
    end


    """
    `Basics.extractShellOccupationFromCsfR(csfR::CsfR, basis::Basis)`  
        ... extracts the nonrelativistic shell list as well as the (nonrelativistic) occupation of the given CsfR. 
            Here, we assume that csfR is defined in the given basis and that all subshells occur in standard sequence
            for a jjJ-LSJ transformation. A tuple tp::Tuple{Array{Shell,1},Array{Int64,1}}  of two arrays of equal 
            length is returned that contains the nonrelativistic shells and their occupation in the specification of 
            csfR. 
    """
    function Basics.extractShellOccupationFromCsfR(csfR::CsfR, basis::Basis)
        shellList = Basics.extractNonrelativisticShellList(basis.subshells) 
        occList   = Int64[]
        
        for  s = 1:length(basis.subshells)
            l = Basics.subshell_l(basis.subshells[s]);   j2 = Basics.subshell_2j(basis.subshells[s])
            if      l == 0          push!( occList, csfR.occupation[s])
            elseif  2*l + 1 == j2   continue    # contributions from j = l + 1/2 are collected together with j = l - 1/2 below
            else    push!( occList, csfR.occupation[s] + csfR.occupation[s+1])
            end
        end 
        ##x println("Basics.extractShellOccupationFromCsfR()::  shellList = $shellList,  occList = $occList")
        
        return( (shellList, occList) )
    end

end # module
