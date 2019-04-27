
"""
`Basics.read()`  ... reads in data from different files and sources.

  + `("CSF list: Grasp92", "cslFilename")`  ... reads in the CSF list as given by the (existing) .csl file cslFilename; an basis::Basis
                 is returned with basis.isDefined = true and with a proper sequence of orbitals, but with basis.orbitals = Orbital[].

  + `("orbital list: Grasp92", "orbFilename")`  ... read in the orbitals from a (formatted) .rwf file orbFilename; a list of 
                 orbitals::Array{OrbitalsR,1} is returned with all subfields specified.
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
Basics.readCslFileGrasp92(filename::String)`  ... reads in the CSF list from a Grasp92 .csl file; a basis::Basis is returned.
"""
function Basics.readCslFileGrasp92(filename::String)
    coreSubshells = Subshell[];    peelSubshells =  Subshell[]

    f  = open(filename)
    sa = readline(f);   sa[1:14] != "Core subshells"   &&   error("Not a Grasp92 .cls file.")
    sa = readline(f);   if  length(sa) > 3    go = true    else  go = false    end
    i  = -1
    while go
        i = i + 1;    sh = JAC.subshellGrasp( strip(sa[5i+1:5i+5]) );    push!(coreSubshells, sh)
        if  5i + 6 >  length(sa)    break    end
    end
    sa = readline(f);   println("sa = $sa")
    sa = readline(f);   if  length(sa) > 3    go = true    else  go = false    end
    i  = -1
    while go
        i = i + 1;    sh = JAC.subshellGrasp( strip(sa[5i+1:5i+5]) );    push!(peelSubshells, sh)
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
        push!(csfs, JAC.ManyElectron.CsfRGrasp92(subshells, coreSubshells, sa, sb, sc) )
    end

    println("  ... $NoCSF CSF read in from Grasp92 file  $filename)") 
    NoElectrons = sum( csfs[1].occupation )

    Basis(true, NoElectrons, subshells, csfs, coreSubshells, Dict{JAC.Subshell,JAC.Radial.Orbital}() )
end


"""
`Basics.readOrbitalFileGrasp92(filename::String, grid.JAC.Radial.Grid)`  ... reads in the orbitals list from a Grasp92 .rwf or .out file; 
                            a dictionary orbitals::Dict{JAC.Subshell,JAC.Radial.Orbital} is returned.
"""
function Basics.readOrbitalFileGrasp92(filename::String, grid::JAC.Radial.Grid)
    # using FortranFiles
    f = FortranFile(filename)
    
    orbitals = Dict{JAC.Subshell,JAC.Radial.Orbital}();    first = true

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
            println("JAC.Basics.readOrbitalFileGrasp92-aa: warning ... Radial functions still need to be interpolated upon the given grid.")
            Px = P
            Qx = Q 

            orbitals = Base.merge( orbitals, Dict( subshell => Orbital(subshell, isBound, useStandardGrid, energy, Px, Qx, JAC.Radial.Grid() ) ))
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
`Basics.readMixFileRelci(filename::String)`  ... reads in the mixing coefficients from a RELCI .mix file; a multiplet::Multiplet is returned.
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
    sa = readline(f);   list2J = AngularJ64[AngularJ64(0)  for  i = 1:NoLevels];    listp = JAC.Parity[JAC.plus  for  i = 1:NoLevels]
    for  i = 1:NoLevels
        jj = Base.parse(Int64, strip(sa[(i-1)*6+1:(i-1)*6+5]) );   list2J[i] = AngularJ64(jj//2)
        pp = string( sa[(i-1)*6+6] );                              listp[i]  = JAC.Parity( pp )
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

    multiplet = JAC.Multiplet(name, levels)
    return( multiplet )
end


"""
`Basics.readMixingFileGrasp92(filename::String, basis::Basis)`  ... reads in the mixing coefficients from a Grasp92 .mix file by using
                           the given basis; A multiplet::Multiplet is returned.
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
    mixfile_to_p(i::Int32) = i==1 ?  JAC.plus : JAC.minus

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
    multiplet = JAC.Multiplet(name, JAC.Level[] )
    return( multiplet )
end
 

