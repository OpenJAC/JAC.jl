
"""
`module  JAC.Basics`  
    ... a submodel of JAC that contains many basic types/struct that are specific to the JAC module; this module also defines
        a number of basic functions/methods that are later extended. We here provide proper docstrings for all abstract and
        conrete types (struct) in order to allow the user to easy understand and access the individual fields of each type/struct 
        definition.
"""
module Basics

    using Printf

    export  AbstractScField, AbstractPotential, AbstractLevelProperty, AbstractProcess, add, analyze, AngularJ, AngularJ64, AngularM, 
            AngularM64, Auger,
            CartesianVector, compute, CorePolarization, diagonalize,
            estimate, EmMultipole, E1, M1, E2, M2, E3, M3, E4, M4, EmProperty, ExpStokes, EmStokes, 
            generate,
            interpolate, integrate, 
            LevelSymmetry, LevelKey, LevelSelection, LineSelection, 
            modify, 
            oplus,
            Photo, Parity, perform, provide, PathwaySelection,
            Radiative,
            SolidAngle, Shell, Subshell, subshell_2j,
            tabulate, tools, TensorComp, 
            UseCoulomb, UseBabushkin,
            WeightedCartesian


    abstract type  AngularJ  end


    """
    `struct  AngularJ64 <: AngularJ`  ... defines a type for angular momenta j.
    
        + num  ::Int64              ... numerator
        + den  ::Int64              ... denominator, must be 1 or 2
    """
    struct  AngularJ64 <: AngularJ
        num      ::Int64
        den      ::Int64
    end


    """
    `Basics.AngularJ64(m::Integer)`  ... constructor for a given integer (numerator).
    """
    function AngularJ64(j::Integer)
        j < 0   &&   error("j must be positive.")
        AngularJ64(j, 1)
    end


    """
    `Basics.AngularJ64(rational::Rational{Int64})`  ... constructor for a given  Rational{Int64}.
    """
    function AngularJ64(rational::Rational{Int64})
        !(rational.den in [1,2])   &&   error("Denominator must be 1 or 2.")
        rational.num < 0           &&   error("j must be positive.")
        AngularJ64(rational.num, rational.den)
    end


    # `Base.show(io::IO, j::AngularJ64)`  ... prepares a proper printout of the variable j::AngularJ64.
    function Base.show(io::IO, j::AngularJ64) 
        if      j.den == 1    print(io, j.num)
        elseif  j.den == 2    print(io, j.num, "/2")
        else    error("stop a")
        end
    end
    

    """
    `Basics.oplus(ja::AngularJ64, jb::AngularJ64)`  
        ... adds the angular momenta ja `oplus` jb and returns a list::Array{AngularJ64,1} of j-valuescin the interval 
            |ja - jb| <= j <= ja + jb.
    """
    function  oplus(ja::AngularJ64, jb::AngularJ64)
        if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
        if  jb.den == 1   jb2 = 2jb.num   else   jb2 = jb.num   end
        jList = AngularJ64[]
        for  j = abs(ja2 - jb2):2:ja2+jb2    push!(jList, AngularJ64(j//2) )    end
        return( jList )
    end
    
    oplus(ja::AngularJ64, jb::Int64) = oplus(ja, AngularJ64(jb))
    oplus(ja::Int64, jb::AngularJ64) = oplus(AngularJ64(ja), jb)
    oplus(ja::Int64, jb::Int64)      = oplus(AngularJ64(ja), AngularJ64(jb))
    
    Base.Float64(ja::AngularJ64) = ja.num / ja.den


    """
    `Basics.projections(ja::AngularJ64)`  
        ... returns all allowed projections of a given angular momenta ja as a list::Array{AngularM64,1} of m-values, i.e. -ja, -ja+1, ..., j.
    """
    function  projections(ja::AngularJ64)
        if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
        mList = AngularM64[]
        for  m = -ja2:2:ja2    push!(mList, AngularM64(m//2) )    end
        return( mList )
    end


    abstract type  AngularM  end

    """
    `struct  Basics.AngularM64 <: AngularM`  ... defines a type for angular momentum projections m.

        + num  ::Int64              ... numerator
        + den  ::Int64              ... denominator, must be 1 or 2
    """
    struct  AngularM64 <: AngularM
        num      ::Int64
        den      ::Int64
    end


    """
    `Basics.AngularM64(m::Integer)`  ... constructor for a given integer (numerator).
    """
    function AngularM64(m::Integer)
        AngularM64(m, 1)
    end


    """
    `Basics.AngularM64(m::Integer, j::AngularJ64)`  
        ... constructor for a given integer (numerator) that is consistent with a given j-value.
    """
    function AngularM64(m::Integer, j::AngularJ64)
        !(j.den == 1)      &&  error("m must be integer for j = $(j).")
        j.num < abs(m)     &&  error("abs(m) must be <= j = $(j).")
        AngularM64(m, 1)
    end


    """
    `Basics.AngularM64(rational::Rational{Int64})`  ... constructor for a given  Rational{Int64}.
    """
    function AngularM64(rational::Rational{Int64})
        !(rational.den in [1,2])  &&  error("Denominator must be 1 or 2.")
        AngularM64(rational.num, rational.den)
    end


    """
    `Basics.AngularM64(rational::Rational{Int64}, j::AngularJ64)`  
        ... constructor for a given  Rational{Int64} that is consistent with a given m-value.
    """
    function AngularM64(rational::Rational{Int64}, j::AngularJ64)
        !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
        !(j.den == rational.den)      &&   error("j,m must be both integer or half-integer.")
        j.num < abs(rational.num)     &&   error("abs(m) must be <= j = $(j).")
        AngularM64(rational.num, rational.den)
    end


    """
    `Basics.AngularM64(j::AngularJ64)`  ... constructor for a given j::AngularJ64 to define m = j.
    """
    function AngularM64(j::AngularJ64)
        AngularM64(j.num, j.den)
    end


    # `Base.show(io::IO, m::AngularM64)`  ... prepares a proper printout of the variable m::AngularM64.
    function Base.show(io::IO, m::AngularM64) 
        if      m.den == 1    print(io, m.num)
        elseif  m.den == 2    print(io, m.num, "/2")
        else    error("stop a")
        end
    end


    """
    `Basics.add(ma::AngularM64, mb::AngularM64)`  
        ... adds the projections of the angular momenta ma + mb and returns a mc::AngularM64.
    """
    function  add(ma::AngularM64, mb::AngularM64)
        if  ma.den == 1   ma2 = 2ma.num   else   ma2 = ma.num   end
        if  mb.den == 1   mb2 = 2mb.num   else   mb2 = mb.num   end
        return( AngularM64( (ma2+mb2)//2 ) )
    end



    # Conversion between HalfInt and AngularJ64, AngularM64

    twice(x::Union{AngularJ64,AngularM64}) = ifelse(isone(x.den), twice(x.num), x.num)
    twice(x) = x + x



    """
    `struct  Basics.Eigen`  
        ... defines a simple struct to communicate eigenvalues and eigenvectors if different diagonalization procedures are used.

        + values   ::Array{Float64,1}            ... List of eigenvalues.
        + vectors  ::Array{Vector{Float64},1}    ... List of eigenvectors.
    """
    struct  Eigen
        values     ::Array{Float64,1}
        vectors    ::Array{Vector{Float64},1}
    end


    """
    `Basics.isEqual(a::Eigen, b::Eigen, accuracy::Float64)`  
        ... determines whether a == b with the requested accuracy; a value::Bool is returned. The test is made element-wise on 
            abs(a.number - b.number) > accuracy.
    """
    function isEqual(a::Eigen, b::Eigen, accuracy::Float64)
        # Check the dimensions
        if  accuracy <= 0.                          error("Requires accuracy > 0.")                                    end
        if  size(a.values) != size(b.values)   ||   size(a.vectors) != size(b.vectors)                return( false )  end
        # Check the eigenvalues and vectors
        for  i = 1:length(a.values)    
            if  abs(a.values[i] - b.values[i]) > accuracy                                             return( false )  end    
            for  j = 1:length(a.vectors[i])   if  abs(a.vectors[i][j] - b.vectors[i][j]) > accuracy   return( false )  end    end     
        end
        return( true )
    end 



    """
    `struct  Basics.SolidAngle`  ... defines a type for a solid angle Omega = (theta, phi).

        + theta  ::Float64       ... polar angle theta with regard to the z-axis; 0 <= theta <= pi.
        + phi    ::Float64       ... azimuthal angle phi (with regard to the x-z plane in mathematical positive direction); 0 <= phi <= 2 pi.
    """
    struct  SolidAngle
        theta    ::Float64
        phi      ::Float64
    end


    # `Base.show(io::IO, Omega::SolidAngle)`  ... prepares a proper printout of the variable Omega::SolidAngle.
    function Base.show(io::IO, Omega::SolidAngle)
        sa = "Omega == (" * string(Omega.theta) * ", " * string(Omega.phi) * ")"
        print(io, sa)
    end


    """
    `@enum   Parity`  ... defines a enumeration for the allowed values of parity.

        + plus, minus               ... with obvious meaning
    """
    @enum   Parity    plus   minus
    export  Parity,   plus,  minus


    """
    `Basics.Parity(sa::String)`  ... constructor for a given String.
    """
    function Parity(sa::String)
        if       sa == "+"    wa = plus
        elseif   sa == "-"    wa = minus
        else     error("stop a")
        end

        return( wa )
    end
    
    
    """
    `Basics.invertParity(p::Parity)`  ... inverts the given parity plus <--> minus.
    """
    function invertParity(p::Parity)
        if       p == plus    return( minus )
        elseif   p == minus   return( plus )
        else     error("stop a")
        end
    end


    # `Base.string(p::Parity)`  ... provides a proper printout of the variable p::Parity.
    function Base.string(p::Parity) 
        if      p == plus     return( "+" )
        elseif  p == minus    return( "-" )  
        else    error("stop a")
        end
    end


    """
    `struct  Basics.LevelSymmetry`  ... defines a struct for defining the overall J^P symmetry of a level.

        + J          ::AngularJ64  ... total angular momentum of a level
        + parity     ::Parity      ... total parity of the level
    """
    struct  LevelSymmetry
        J            ::AngularJ64
        parity       ::Parity  
    end


    """
    `Basics.LevelSymmetry(rational::Rational{Int64}, parity::Parity)`  ... constructor for a given (Rational,Parity).
    """
    function  LevelSymmetry(rational::Rational{Int64}, parity::Parity)
        !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
        LevelSymmetry( AngularJ64(rational.num, rational.den), parity )    
    end


    """
    `Basics.LevelSymmetry(i::Int64, parity::Parity)`  ... constructor for a given (Int64,Parity).
    """
    function  LevelSymmetry(i::Int64, parity::Parity)
        LevelSymmetry( AngularJ64(i), parity )    
    end


    """
    `Basics.LevelSymmetry(rational::Rational{Int64}, sa::String)`  ... constructor for a given (Rational,String).
    """
    function  LevelSymmetry(rational::Rational{Int64},sa::String)
        !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
        LevelSymmetry( AngularJ64(rational.num, rational.den), Parity(sa) )    
    end


    """
    `Basics.LevelSymmetry(i::Int64, sa::String)`  ... constructor for a given (Int64,String).
    """
    function  LevelSymmetry(i::Int64,sa::String)
        LevelSymmetry( AngularJ64(i), Parity(sa) )    
    end


    # `Base.show(io::IO, sym::LevelSymmetry)`  ... prepares a proper printout of the variable sym::LevelSymmetry.
    function Base.show(io::IO, sym::LevelSymmetry) 
        print(io, string(sym) )
    end


    # `Base.string(sym::LevelSymmetry)`  ... provides a proper printout of the variable sym::LevelSymmetry.
    function Base.string(sym::LevelSymmetry) 
        if      sym.parity == plus     return( "$(sym.J) +" )
        elseif  sym.parity == minus    return( "$(sym.J) -" )  
        else    error("stop a")
        end
    end


    """
    `struct  Basics.LevelKey`  ... defines a struct for identifying a level by its symmetry, energy, etc.

        + sym          ::LevelSymmetry    ... Level symmetry.
        + index        ::Int64            ... Index of this level in its original multiplet, or 0.
        + energy       ::Float64          ... Total energy.
        + relativeOcc  ::Float64          ... Relative occupation of this level, if part of some multiplet/cascade.
    """
    struct  LevelKey
        sym            ::LevelSymmetry
        index          ::Int64    
        energy         ::Float64 
        relativeOcc    ::Float64 
    end


    """
    `Basics.LevelKey()`  ... constructor for an empty instance of a LevelKey.
    """
    function  LevelKey()
        LevelKey( LevelSymmetry(0,Basics.plus), 0, 0., 0.)    
    end


    # `Base.show(io::IO, key::LevelKey)`  ... prepares a proper printout of the variable key::LevelKey.
    function Base.show(io::IO, key::LevelKey) 
        print(io, string(key) )
    end


    # `Base.string(key::LevelKey)`  ... provides a proper printout of the variable key::LevelKey.
    function Base.string(key::LevelKey) 
        return( "LevelKey[$(key.sym) ($(key.index)), en=$(key.energy), occ=$(key.relativeOcc)]")
    end


    """
    `struct  Basics.LineKey`  
        ... defines a struct for identifying a line by the keys of the initial and final level, its transition energy,
            decay strength, etc.

        + iLevelKey    ::LevelKey         ... Key of the initial level.
        + fLevelKey    ::LevelKey         ... Key of the final level.
        + energy       ::Float64          ... Total energy.
        + strength     ::Float64          ... Strength of the given line.
    """
    struct  LineKey
        iLevelKey      ::LevelKey  
        fLevelKey      ::LevelKey
        energy         ::Float64   
        strength       ::Float64 
    end


    """
    `Basics.LineKey()`  ... constructor for an empty instance of a LineKey.
    """
    function  LineKey()
        LineKey( LevelKey(), LevelKey(), 0., 0.)    
    end


    # `Base.show(io::IO, key::LineKey)`  ... prepares a proper printout of the variable key::LineKey.
    function Base.show(io::IO, key::LineKey) 
        print(io, string(key) )
    end


    # `Base.string(key::LineKey)`  ... provides a proper printout of the variable key::LevelKey.
    function Base.string(key::LineKey) 
        return( "LineKey[i=$(key.iLevelKey) --> f=$(key.fLevelKey), en=$(key.energy), strength=$(key.strength)]")
    end


    
    """
    `struct  Basics.LevelSelection`  
        ... defines a struct to specify a list of levels by means of their (level) indices or level symmetries.

        + active       ::Bool                     ... true, if some selection has been made.
        + indices      ::Array{Int64,1}           ... List of selected indices.
        + symmetries   ::Array{LevelSymmetry,1}   ... List of selected symmetries
    """
    struct  LevelSelection
        active         ::Bool  
        indices        ::Array{Int64,1}
        symmetries     ::Array{LevelSymmetry,1}
    end


    """
    `Basics.LevelSelection()`  ... constructor for an inactive LevelSelection.
    """
    function  LevelSelection()
        LevelSelection( false, Int64[], LevelSymmetry[])    
    end


    """
    `Basics.LevelSelection(active::Bool; indices::Array{Int64,1}=Int64[], symmetries::Array{LevelSymmetry,1}=LevelSymmetry[])`  
        ... constructor for specifying the details of a LevelSelection.
    """
    function  LevelSelection(active::Bool; indices::Array{Int64,1}=Int64[], symmetries::Array{LevelSymmetry,1}=LevelSymmetry[])
        if  active   LevelSelection( true, indices, symmetries)  
        else         LevelSelection()
        end
    end

    function Base.show(io::IO, selection::LevelSelection) 
        print(io, string(selection) )
    end

    function Base.string(selection::LevelSelection) 
        if  selection.active   sa = "LevelSelection:  indices = $(selection.indices);    symmetries = $(selection.symmetries)."
        else                   sa = "Inactive LevelSelection."
        end
        return( sa )
    end


    """
    `struct  Basics.LineSelection`  
        ... defines a struct to specify a list of level pair by means of their (level) indices or level symmetries.

        + active        ::Bool                                          ... true, if some selection has been made.
        + indexPairs    ::Array{Tuple{Int64,Int64},1}                   ... List of selected index pairs.
        + symmetryPairs ::Array{Tuple{LevelSymmetry,LevelSymmetry},1}   ... List of selected symmetry pairs.
    """
    struct  LineSelection
        active          ::Bool  
        indexPairs      ::Array{Tuple{Int64,Int64},1}
        symmetryPairs   ::Array{Tuple{LevelSymmetry,LevelSymmetry},1}
    end


    """
    `Basics.LineSelection()`  ... constructor for an inactive LineSelection.
    """
    function  LineSelection()
        LineSelection( false, Tuple{Int64,Int64}[], Tuple{LevelSymmetry,LevelSymmetry}[])    
    end


    """
    `Basics.LineSelection(active::Bool; indexPairs::Array{Tuple{Int64,Int64},1}=Tuple{Int64,Int64}[],
                                        symmetryPairs::Array{Tuple{LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry}[])`  
        ... constructor for specifying the details of a LineSelection.
    """
    function  LineSelection(active::Bool; indexPairs::Array{Tuple{Int64,Int64},1}=Tuple{Int64,Int64}[],
                                          symmetryPairs::Array{Tuple{LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry}[])
        if  active   LineSelection( true, indexPairs, symmetryPairs)  
        else         LineSelection()
        end
    end

    function Base.show(io::IO, selection::LineSelection) 
        print(io, string(selection) )
    end

    function Base.string(selection::LineSelection) 
        if  selection.active   sa = "LineSelection:  indexPairs = $(selection.indexPairs);    symmetryPairs = $(selection.symmetryPairs)."
        else                   sa = "Inactive LineSelection."
        end
        return( sa )
    end


    """
    `struct  Basics.PathwaySelection`  
        ... defines a struct to specify a list of level triple (pathways) by means of their (level) indices or level symmetries.

        + active          ::Bool                                          ... true, if some selection has been made.
        + indexTriples    ::Array{Tuple{Int64,Int64,Int64},1}             ... List of selected index triples.
        + symmetryTriples ::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}  ... List of selected symmetry triples.
    """
    struct  PathwaySelection
        active            ::Bool  
        indexTriples      ::Array{Tuple{Int64,Int64,Int64},1}  
        symmetryTriples   ::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}
    end


    """
    `Basics.PathwaySelection()`  ... constructor for an inactive PathwaySelection.
    """
    function  PathwaySelection()
        PathwaySelection( false, Tuple{Int64,Int64,Int64}[], Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry}[])    
    end


    """
    `Basics.PathwaySelection(active::Bool; indexTriples::Array{Tuple{Int64,Int64,Int64},1}=Tuple{Int64,Int64,Int64}[],
                   symmetryTriples::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry}[])`  
        ... constructor for specifying the details of a PathwaySelection.
    """
    function  PathwaySelection(active::Bool; indexTriples::Array{Tuple{Int64,Int64,Int64},1}=Tuple{Int64,Int64,Int64}[],
                symmetryTriples::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry}[])
        if  active   PathwaySelection( true, indexTriples, symmetryTriples)  
        else         PathwaySelection()
        end
    end

    function Base.show(io::IO, selection::PathwaySelection) 
        print(io, string(selection) )
    end

    function Base.string(selection::PathwaySelection) 
        if  selection.active   sa = "PathwaySelection:  indexTriples = $(selection.indexTriples);    symmetryTriples = $(selection.symmetryTriples)."
        else                   sa = "Inactive PathwaySelection."
        end
        return( sa )
    end

    
    
    """
    `struct  Basics.Shell`  ... defines a enumeration for the allowed values of a non-relativistic shell.  

        + n     ::Int64  ... principal quantum number 
        + l     ::Int64  ... orbital angular quantum number
    """
    struct  Shell
        n       ::Int64 
        l       ::Int64 
    end


    """
    `Basics.Shell(sa::Union{String,SubString{String}})`  ... constructor for a given String, such as 1s, 2s, 2p, ... .
    """
    function Shell(sa::Union{String,SubString{String}}) 
        wa = strip( Base.convert(String, sa) )
        n  = parse(Int64, wa[1:end-1]);    l = shellNotation( string(wa[end]) )
        Shell(n, l)   
    end


    # `Base.show(io::IO, sh::Shell)`  ... prepares a proper printout of the variable sh::Shell.
    function Base.show(io::IO, sh::Shell) 
        print(io, string(sh) )
    end


    # `Base.string(sh::Shell)`  ... provides a proper printout of the variable sh::Shell.
    function Base.string(sh::Shell) 
        return( string(sh.n) * shellNotation(sh.l) )
    end


    """
    `Basics.shellNotation(l::Int64)`  ... returns the corresponding spectroscopy letter for a given orbital angular quantum number l.
    """
    function shellNotation(l::Int64)
        !(0 <= l <= 25)      &&   error("Orbital QN 0 <= l <= 25; l = $l")
        wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
        return( wa[l+1] )   
    end


    """
    `Basics.shellNotation(sa::String)`  ... returns for a given spectroscopy letter the orbital angular quantum number l.
    """
    function shellNotation(sa::String)
        wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
        wb = findall(in([sa]), wa)
        !(length(wb)==1)      &&   error("Unrecognized orbital letter = $sa")
        return( wb[1]-1)   
    end


    """
    `Basics.shellSplitIntoSubshells(sh::Shell)'  ... to convert a non-relativistic shell into a list of Subshell[..]
    """
    function shellSplitIntoSubshells(sh::Shell)
        subshells = Subshell[]
        if     sh.l == 0   push!(subshells, Subshell(sh.n, -1) )
        else   j2 = 2sh.l - 1;    kappa = Int64(  (j2+1)/2 );   push!(subshells, Subshell(sh.n, kappa) )
            j2 = 2sh.l + 1;    kappa = Int64( -(j2+1)/2 );   push!(subshells, Subshell(sh.n, kappa) )
        end
        return( subshells )
    end


    """
    `Basics.shellSplitOccupation(sh::Shell, occ::Int64)'  
        ... to split the occupation of a shell into proper subshell occupations; an  Array{Dict{Basics.Subshell,Int64},1} is returned.
    """
    function shellSplitOccupation(sh::Shell, occ::Int64)

        subshells = shellSplitIntoSubshells(sh)
        wa        = Dict{Subshell,Int64}[]

        for  k1 in 0:occ
            for  k2 in 0:occ
                if  k1 + k2 == occ
                    if       length( subshells ) == 1  
                        !(0 <= occ <= Basics.subshell_2j( subshells[1] ) + 1)   &&    error("Wrong pair of subshell = $(subshells[1]) and occupation = $occ") 
                        wb = Dict( subshells[1] => occ )
                        push!(wa, wb)
                        return( wa )

                        elseif  length( subshells ) == 2
                            j2_1 = Basics.subshell_2j(subshells[1])
                        j2_2 = Basics.subshell_2j(subshells[2])
                        if  k1 <= j2_1 + 1   &&    k2 <= j2_2 + 1    
                        wb = Dict( subshells[1] => k1, subshells[2] => k2 )
                        # Delete contributions with zero occupation
                        ## if  k1 == 0    delete!(wb, subshells[1])   end
                        ## if  k2 == 0    delete!(wb, subshells[2])   end
                        push!(wa, wb)
                        end
                    else    error("stop a") 
                    end
                end
            end
        end

        return( wa )
    end


    """
    `struct  Basics.Subshell`  ... defines a type for the allowed values of a relativistic subhell.  

        + n        ::Int64  ... principal quantum number 
        + kappa    ::Int64  ... relativistic angular quantum number
    """
    struct  Subshell
        n          ::Int64 
        kappa      ::Int64 
    end


    """
    `Basics.Subshell(sa::String)`  ... constructor for a given String, such as 1s_1/2, 2s_1/2, 2p_3/2, ... .
    """
    function Subshell(sa::String) 
        wa = strip( sa );   
        wb = findnext("_", wa, 1);    if  wb == nothing   error("No underscore in sa = $sa")           else   wb = wb[1]   end
        wc = findnext("/", wa, wb);   if  wc == nothing   error("No / in wa[] = $(wa[wb[1]+1:end])")   else   wc = wc[1]   end
        wd = string( wa[wb+1:wc-1] )
        n  = parse(Int64, wa[1:wb-2]);                  l  = shellNotation( string(wa[ wb-1 ]) )
        j2 = parse(Int64, wd );                         !(string(wa[end]) == "2")    &&  error("Unrecognized subshell string sa = $sa") 
        if  l - j2/2 > 0     kappa = Int64( (j2+1)/2 )    else    kappa = Int64( -(j2+1)/2 )   end
        Subshell(n, kappa)   
    end


    """
    `Basics.Subshell(n::Int64, symmetry::LevelSymmetry)`  ... constructor for a given principal quantum number n and (level) symmetry.
    """
    function Subshell(n::Int64, symmetry::LevelSymmetry) 
        if  symmetry.parity == Basics.plus
            if      symmetry.J == AngularJ64(1//2)   kappa = -1
            elseif  symmetry.J == AngularJ64(3//2)   kappa =  2
            elseif  symmetry.J == AngularJ64(5//2)   kappa = -3
            elseif  symmetry.J == AngularJ64(7//2)   kappa =  4
            elseif  symmetry.J == AngularJ64(9//2)   kappa = -5
            elseif  symmetry.J == AngularJ64(11//2)  kappa =  6
            elseif  symmetry.J == AngularJ64(13//2)  kappa = -7
            else    error("stop a")
            end
        else
            if      symmetry.J == AngularJ64(1//2)   kappa =  1
            elseif  symmetry.J == AngularJ64(3//2)   kappa = -2
            elseif  symmetry.J == AngularJ64(5//2)   kappa =  3
            elseif  symmetry.J == AngularJ64(7//2)   kappa = -4
            elseif  symmetry.J == AngularJ64(9//2)   kappa =  5
            elseif  symmetry.J == AngularJ64(11//2)  kappa = -6
            elseif  symmetry.J == AngularJ64(13//2)  kappa =  7
            else    error("stop b")
            end
        end
        
        return( Subshell(n,kappa) )
    end


    # `Base.show(io::IO, sh::Subshell)`  ... prepares a proper printout of the variable sh::Subshell.
    function Base.show(io::IO, sh::Subshell) 
        print(io, string(sh) )
    end


    # `Base.string(sh::Subshell)`  ... provides a proper printout of the variable sh::Subshell.
    function Base.string(sh::Subshell) 
        sa = string(sh.n) * shellNotation( Basics.subshell_l(sh) ) * "_" * string( Basics.subshell_j(sh) )
        return( sa )
    end


    """
    `Basics.subshellGrasp(sa::Stringl)`  ... returns a sh::Subshell if sa denotes a (relativistic) Subshell in Grasp-notation.
    """
    function subshellGrasp(sa::String)
            if   sa == "1s"   return( Subshell("1s_1/2")  )   
        elseif   sa == "2s"   return( Subshell("2s_1/2")  )   elseif   sa == "2p-"   return( Subshell("2p_1/2") )
        elseif   sa == "2p"   return( Subshell("2p_3/2")  ) 
        elseif   sa == "3s"   return( Subshell("3s_1/2")  )   elseif   sa == "3p-"   return( Subshell("3p_1/2") )
        elseif   sa == "3p"   return( Subshell("3p_3/2")  )   elseif   sa == "3d-"   return( Subshell("3d_3/2") )
        elseif   sa == "3d"   return( Subshell("3d_5/2")  ) 
        elseif   sa == "4s"   return( Subshell("4s_1/2")  )   elseif   sa == "4p-"   return( Subshell("4p_1/2") )
        elseif   sa == "4p"   return( Subshell("4p_3/2")  )   elseif   sa == "4d-"   return( Subshell("4d_3/2") )
        elseif   sa == "4d"   return( Subshell("4d_5/2")  )   elseif   sa == "4f-"   return( Subshell("4f_5/2") )
        elseif   sa == "4f"   return( Subshell("4f_7/2")  ) 
        elseif   sa == "5s"   return( Subshell("5s_1/2")  )   elseif   sa == "5p-"   return( Subshell("5p_1/2") )
        elseif   sa == "5p"   return( Subshell("5p_3/2")  )   elseif   sa == "5d-"   return( Subshell("5d_3/2") )
        elseif   sa == "5d"   return( Subshell("5d_5/2")  )   elseif   sa == "5f-"   return( Subshell("5f_5/2") )
        elseif   sa == "5f"   return( Subshell("5f_7/2")  )   elseif   sa == "5g-"   return( Subshell("5g_7/2") )
        elseif   sa == "5g"   return( Subshell("5g_9/2")  ) 
        elseif   sa == "6s"   return( Subshell("6s_1/2")  )   elseif   sa == "6p-"   return( Subshell("6p_1/2") )
        elseif   sa == "6p"   return( Subshell("6p_3/2")  )   elseif   sa == "6d-"   return( Subshell("6d_3/2") )
        elseif   sa == "6d"   return( Subshell("6d_5/2")  )   elseif   sa == "6f-"   return( Subshell("6f_5/2") )
        elseif   sa == "6f"   return( Subshell("6f_7/2")  )   elseif   sa == "6g-"   return( Subshell("6g_7/2") )
        elseif   sa == "6g"   return( Subshell("6g_9/2")  ) 
        elseif   sa == "7s"   return( Subshell("7s_1/2")  )   elseif   sa == "7p-"   return( Subshell("7p_1/2") )
        elseif   sa == "7p"   return( Subshell("7p_3/2")  )   elseif   sa == "7d-"   return( Subshell("7d_3/2") )
        else     error("Grasp notation ... not yet implemented; sa = $sa")
        end
    end


    """
    `Basics.subshell_l(sh::Subshell)`  ... returns the orbital angular quantum number l for a given sh::Subshell; an Int64 is returned
    """
    function subshell_l(sh::Subshell)
        if   sh.kappa < 0   l = abs(sh.kappa) -1   else   l = sh.kappa   end
        return( l )   
    end


    """
    `Basics.subshell_j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
    """
    function subshell_j(sh::Subshell)
        return( AngularJ64(abs(sh.kappa) - 1//2) )   
    end


    """
    `Basics.subshell_2j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
    """
    function subshell_2j(sh::Subshell)
        return( 2* abs(sh.kappa) - 1 )   
    end


    """
    `Basics.subshellsFromClosedShellConfiguration("[Ne]")`  
        ... to provide a list of (relativistic) subshells for the given closed-shell configuration.
    """
    function subshellsFromClosedShellConfiguration(sa::String)
        if       sa == "[He]"    wa = [ Subshell("1s_1/2")]    
        elseif   sa == "[Ne]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]   
        elseif   sa == "[Mg]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3s_1/2")]   
        elseif   sa == "[Ar]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                            Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2")]   
        elseif   sa == "[Kr]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                            Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"),   
                                                            Subshell("3d_3/2"), Subshell("3d_5/2"), 
                                        Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2")]   
        elseif   sa == "[Xe]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                            Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"),   
                                                            Subshell("3d_3/2"), Subshell("3d_5/2"), 
                                        Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2"), 
                                                            Subshell("4d_3/2"), Subshell("4d_5/2"), 
                                        Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2")]   
        else
            error("Unsupported keystring = $sa ")
        end
        
        return( wa )
    end


    """
    `struct   Basics.SubshellStateR`  
        ... defines a struct for the relativistic antisymmetric subshell states in the seniority scheme; this struct can be accessed 
            only internally and, therefore, the only the standard constructor is supported.

        + subshell     ::Subshell        ... to refer to a particular subshell. 
        + occ          ::Int64           ... occupation of the subshell
        + nu           ::Int64           ... seniority of the subshell state, often called nu
        + Jsub2        ::Int64           ... 2*Jsub
    """
    struct  SubshellStateR
        subshell       ::Subshell
        occ            ::Int64
        nu             ::Int64
        Jsub2          ::Int64
    end     


    # `Base.show(io::IO, state::SubshellStateR)`  ... prepares a proper printout of the variable state::SubshellStateR.
    function Base.show(io::IO, state::SubshellStateR) 
        print(io, string(state) )
    end


    # `Base.string(state::SubshellStateR)`  ... provides a proper printout of the variable state::SubshellStateR.
    function Base.string(state::SubshellStateR)
        sa = "SubshellState: [" * string(state.subshell) * "^$(state.occ), nu=$(state.nu), 2*J_sub=$(state.Jsub2)]"
    end


    """
    `Basics.subshellStateString(subshell::String, occ::Int64, seniorityNr::Int64, Jsub::AngularJ64, X::AngularJ64)`  
        ... to provide a string of a given subshell state in the form '[2p_1/2^occ]_(seniorityNr, J_sub), X=Xo' ... .
    """
    function subshellStateString(subshell::String, occ::Int64, seniorityNr::Int64, Jsub::AngularJ64, X::AngularJ64)
        sa = "[" * subshell * "^$occ]_($seniorityNr, " * string(Jsub) * ") X=" * string(X)
        return( sa )
    end



    """
    `struct  Basics.EmMultipole`  ... defines a struct for the em multipoles, E1, M1, E2, ...

        + L           ::Int64    ... multipolarity
        + electric    ::Bool     ... True if the multipole is an electric one, and false otherwise.
    """  
    struct  EmMultipole
        L             ::Int64
        electric      ::Bool   
    end


    const  E1 = EmMultipole(1, true);     const  M1 = EmMultipole(1, false)
    const  E2 = EmMultipole(2, true);     const  M2 = EmMultipole(2, false)
    const  E3 = EmMultipole(3, true);     const  M3 = EmMultipole(3, false)
    const  E4 = EmMultipole(4, true);     const  M4 = EmMultipole(4, false)


    """
    `Basics.EmMultipole(sa::String)`  ... constructor for a given String.
    """
    function EmMultipole(sa::String)
        if       sa == "E1"    L = 1;   electric = true
        elseif   sa == "M1"    L = 1;   electric = false
        elseif   sa == "E2"    L = 2;   electric = true
        elseif   sa == "M2"    L = 2;   electric = false
        elseif   sa == "E3"    L = 3;   electric = true
        elseif   sa == "M3"    L = 3;   electric = false
        elseif   sa == "E4"    L = 4;   electric = true
        elseif   sa == "M4"    L = 4;   electric = false
        else     error("stop a")
        end

        EmMultipole(L, electric)
    end


    # `Base.show(io::IO, mp::EmMultipole)`  ... prepares a proper printout of the variable mp::EmMultipole.
    function Base.show(io::IO, mp::EmMultipole) 
        if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
        print(io, string(mp) )
    end


    # `Base.string(mp::EmMultipole)`  ... provides a proper printout of the variable mp::EmMultipole.
    function Base.string(mp::EmMultipole)
        mp.L < 1    &&    error("Improper multipolarity L = $(mp.L)") 
        if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
        return( sa )
    end


    """
    `Basics.multipole_p(mp::EmMultipole)`  ... returns the kind integer p=0 (magnetic) or p=1 (electric) of the given  mp::EmMultipole.
    """
    function multipole_p(mp::EmMultipole)
        if  mp.electric   return(1)    else    return(0)   end
    end


    """
    `@enum   Basics.EmGauge`  
        ... defines a enumeration for the (allowed) gauges to deal with the radiation field; this differs from the UseGauge which can 
            be selected for explicit computations.
                        
            + NoGauge      ... No gauge define (yet).
            + Coulomb      ... Coulomb gauge (= velocity gauge in the non-relativistic limit; Grant, 1974).
            + Babushkin    ... Babushkin gauge (= length gauge in the non-relativistic limit; Grant, 1974). 
            + Magnetic     ... Magnetic gauge (transverse velocity gauge).
            + Velocity     ... Transverse velocity gauge (Johnson, 2007).
            + Length       ... Length gauge (Johnson, 2007).
    """
    @enum   EmGauge    NoGauge    Coulomb   Babushkin   Magnetic   Velocity   Length
    export  EmGauge,   NoGauge,   Coulomb,  Babushkin,  Magnetic,  Velocity,  Length


    """
    `Basics.EmGauge(sa::String)`  ... constructor for a given String.
    """
    function EmGauge(sa::String)
        if       sa in ["none", "NoGauge"]         wa = NoGauge
        elseif   sa == "Coulomb"                   wa = Coulomb
        elseif   sa == "Babushkin"                 wa = Babushkin
        elseif   sa in ["Magnetic", "magnetic"]    wa = Magnetic
        else     error("stop a")
        end

        return( wa )
    end


    # `Base.string(gauge::EmGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
    function Base.string(gauge::EmGauge) 
        if      gauge == NoGauge     return( "no gauge" )
        elseif  gauge == Coulomb     return( "Coulomb" )   
        elseif  gauge == Babushkin   return( "Babushkin" )   
        elseif  gauge == Magnetic    return( "Magnetic" )   
        else    error("stop a")
        end
    end


    """
    `@enum   Basics.UseGauge`  ... defines a enumeration for the (allowed) gauges to be selected in explicit computations
    """
    @enum   UseGauge    UseCoulomb   UseBabushkin
    export  UseGauge,   UseCoulomb,  UseBabushkin


    """
    `Basics.UseGauge(sa::String)`  ... constructor for a given String.
    """
    function UseGauge(sa::String)
        if       sa == "Coulomb"                   wa = UseCoulomb
        elseif   sa == "Babushkin"                 wa = UseBabushkin
        else     error("stop a")
        end

        UseGauge(wa)
    end


    # `Base.string(gauge::UseGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
    function Base.string(gauge::UseGauge) 
        if      gauge == UseCoulomb     return( "Coulomb" )   
        elseif  gauge == UseBabushkin   return( "Babushkin" )   
        else    error("stop a")
        end
    end



    """
    `struct  EmProperty`  ... defines a type to maintain two gauge forms of a computed result that depends on the radiation field.

        + Coulomb         ::Float64    ... Value for the Coulomb gauge of the radiation field.
        + Babushkin       ::Float64    ... Value for the Coulomb gauge of the radiation field.
    """
    struct  EmProperty 
        Coulomb           ::Float64
        Babushkin         ::Float64
    end 


    """
    `Basics.EmProperty(wa::Float64)`  ... constructor for an `constant` instance of EmProperty.
    """
    function EmProperty(wa::Float64)
        EmProperty(wa, wa)
    end

    
    Base.:+(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb + b.Coulomb, a.Babushkin + b.Babushkin)
    Base.:+(a::EmProperty, b) = EmProperty(a.Coulomb + b, a.Babushkin + b)
    Base.:+(a, b::EmProperty) = b + a
    Base.:*(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb * b.Coulomb, a.Babushkin * b.Babushkin)
    Base.:*(a, b::EmProperty) = EmProperty(a * b.Coulomb, a * b.Babushkin)
    Base.:/(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb / b.Coulomb, a.Babushkin / b.Babushkin)
    Base.:/(a, b::EmProperty) = EmProperty(a / b.Coulomb, a / b.Babushkin)


    # `Base.show(io::IO, property::EmProperty)`  ... prepares a proper printout of the variable property::EmProperty.
    function Base.show(io::IO, property::EmProperty) 
        sa = Base.string(property);                print(io, sa)
    end


    # `Base.string(property::EmProperty)`  ... provides a String notation for the variable property::EmProperty.
    function Base.string(property::EmProperty)
        sa = "$(property.Coulomb) [Coulomb],  $(property.Babushkin) [Babushkin]"
        return( sa )
    end


    """
    `struct  Basics.ExpStokes`  
        ... defines a type to maintain the (experimentally) given Stokes parameter for the polarization of incoming 
            photons or electrons.

        + P1              ::Float64    ... Stokes P1 parameter.
        + P2              ::Float64    ... Stokes P2 parameter.
        + P3              ::Float64    ... Stokes P3 parameter.
    """
    struct  ExpStokes 
        P1                ::Float64
        P2                ::Float64
        P3                ::Float64
    end 


    """
    `Basics.ExpStokes()`  ... constructor for an (empty or unpolarized) Stokes vector.
    """
    function ExpStokes()
        ExpStokes(0., 0., 0.)
    end


    # `Base.show(io::IO, ps::ExpStokes)`  ... prepares a proper printout of the variable ps::ExpStokes.
    function Base.show(io::IO, ps::ExpStokes) 
        sa = Base.string(ps);                print(io, sa)
    end


    # `Base.string(ps::ExpStokes)`  ... provides a String notation for the variable ps::ExpStokes.
    function Base.string(ps::ExpStokes)
        sa = "exp. Stokes parameters P1 = $(ps.P1),  P2 = $(ps.P2),  P3 = $(ps.P3)"
        return( sa )
    end


    """
    `struct  Basics.EmStokes`  
        ... defines a type to maintain the (computed) given Stokes parameter for the polarization of emitted photons or electrons.

        + P1              ::EmProperty    ... Stokes P1 parameter in Coulomb and Babushkin gauge
        + P2              ::EmProperty    ... Stokes P2 parameter.
        + P3              ::EmProperty    ... Stokes P3 parameter.
    """
    struct  EmStokes 
        P1                ::EmProperty
        P2                ::EmProperty
        P3                ::EmProperty
    end 



    # `Base.show(io::IO, ps::EmStokes)`  ... prepares a proper printout of the variable ps::EmStokes.
    function Base.show(io::IO, ps::EmStokes) 
        sa = Base.string(ps);                print(io, sa)
    end


    # `Base.string(ps::EmStokes)`  ... provides a String notation for the variable ps::EmStokes.
    function Base.string(ps::EmStokes)
        sa = "Stokes parameter P1 = $(ps.P1),  P2 = $(ps.P2),  P3 = $(ps.P3)  "
        return( sa )
    end

    
    
    """
    `abstract type Basics.AbstractPolarization` 
        ... defines an abstract type to comprise various polarizations of light and electron beams.

        + LinearPolarization        ... to specify a linearly-polarized pulse/beam.
        + LeftCircular              ... to specify a left-circularly polarized pulse/beam.
        + RightCircular             ... to specify a right-circularly polarized pulse/beam.
        + LeftElliptical            ... to specify an elliptically polarized pulse/beam.
        + RightElliptical           ... to specify an elliptically polarized pulse/beam.
        + NonePolarization          ... to specify an upolarized pulse/beam.
        + DensityMatrixPolarization ... to specify the polarization of a pulse/beam by its (2x2) 
                                        density matrix (not yet).
    """
    abstract type  AbstractPolarization  end
    
    struct         LinearPolarization      <:  AbstractPolarization   end
    struct         LeftCircular            <:  AbstractPolarization   end
    struct         RightCircular           <:  AbstractPolarization   end
    
    
    """
    `struct     Basics.LeftElliptical          <:  Basics.AbstractPolarization`   
    
            + ellipticity      ::Float64     ... Ellipticity of the beam in the range 0...1.
    """
    struct         LeftElliptical          <:  AbstractPolarization
            ellipticity      ::Float64
    end
    
    
    """
    `struct     Basics.RightElliptical          <:  Basics.AbstractPolarization`   
    
            + ellipticity      ::Float64     ... Ellipticity of the beam in the range 0...1.
    """
    struct         RightElliptical          <:  AbstractPolarization
            ellipticity      ::Float64
    end
    
    
    struct         NonePolarization        <:  AbstractPolarization   end

    function Base.string(pol::LinearPolarization)   return( "linearly-polarized" )            end
    function Base.string(pol::LeftCircular)         return( "left-circularly polarized" )     end
    function Base.string(pol::RightCircular)        return( "right-circularly polarized" )    end
    function Base.string(pol::LeftElliptical)       return( "left-elliptically polarized with ellipticity $(pol.ellipticity)" )   end
    function Base.string(pol::RightElliptical)      return( "right-elliptically polarized with ellipticity $(pol.ellipticity)" )  end
    function Base.string(pol::NonePolarization)     return( "unpolarized" )                   end

    
    """
    `struct  Basics.TensorComp`  ... defines a type for a component of the statistical tensor as associated with an atomic or ionic level.

        + k               ::Int64          ... rank of the tensor component
        + q               ::Int64          ... q-component
        + comp            ::ComplexF64     ... value of the tensor component
    """
    struct  TensorComp 
        k                 ::Int64
        q                 ::Int64 
        comp              ::ComplexF64
    end 



    # `Base.show(io::IO, st::TensorComp)`  ... prepares a proper printout of the variable st::TensorComp.
    function Base.show(io::IO, st::TensorComp) 
        sa = Base.string(st);                print(io, sa)
    end


    # `Base.string(st::TensorComp)`  ... provides a String notation for the variable st::TensorComp.
    function Base.string(st::TensorComp)
        sa = "Statistical tensor component rho_($(st.k),$(st.q)) = $(st.comp)  "
        return( sa )
    end


    """
    `Basics.tensorComp(k::Int64, q::Int64, tcomps::Array{TensorComp,1};  withZeros::Bool = false)`  
        ... returns the value of the statistical tensor component rho_kq if it is contained in tcomps.
    """
    function tensorComp(k::Int64, q::Int64, tcomps::Array{TensorComp,1};  withZeros::Bool = false)
        for  tc in  tcomps
            if   k == tc.k    &&    q == tc.q    return( tc.comp )    end
        end
        
        if   withZeros    return( ComplexF64(0.) )    else    error("Statistical tensor component not found for k = $k and q = $q ")    end
    end


    """
    `struct  Basics.CartesianVector`  
        ... defines a type to represent a point in Cartesian{Type} coordinates.

        + x              ::Type       ... x-coordinate.
        + y              ::Type       ... y-coordinate.
        + z              ::Type       ... z-coordinate.
    """
    struct  CartesianVector{Type} 
        x                ::Type
        y                ::Type
        z                ::Type
    end 



    # `Base.show(io::IO, vector::CartesianVector)`  ... prepares a proper printout of the variable vector::CartesianVector.
    function Base.show(io::IO, vector::CartesianVector) 
        sa = Base.string(vector);                print(io, sa)
    end


    # `Base.string(vector::CartesianVector)`  ... provides a String notation for the variable vector::CartesianVector.
    function Base.string(vector::CartesianVector)
        sa = "Cartesian vector V.x = $(vector.x),  V.y = $(vector.y),  V.z = $(vector.z)  "
        return( sa )
    end


    """
    `struct  Basics.WeightedCartesian{Type}`  
        ... defines a type to represent a point in Cartesian{Type} coordinates with weight.

        + x              ::Type       ... x-coordinate.
        + y              ::Type       ... y-coordinate.
        + z              ::Type       ... z-coordinate.
        + w              ::Float64    ... (real) weight.
    """
    struct  WeightedCartesian{Type}
        x                ::Type
        y                ::Type
        z                ::Type
        w                ::Float64
    end 



    # `Base.show(io::IO, point::WeightedCartesian)`  ... prepares a proper printout of the variable point::WeightedCartesian.
    function Base.show(io::IO, vector::WeightedCartesian) 
        sa = Base.string(vector);                print(io, sa)
    end


    # `Base.string(point::WeightedCartesian)`  ... provides a String notation for the variable point::WeightedCartesian.
    function Base.string(vector::WeightedCartesian)
        sa = "Cartesian vector V.x = $(vector.x),  V.y = $(vector.y),  V.z = $(vector.z) with weight w = $(vector.w) "
        return( sa )
    end


    """
    `struct  Basics.CorePolarization`  
        ... defines a type to describe the influence of core-polarization upon the potential or interaction operator.

        + doApply          ::Bool                ... True if core-polarization is to be applied, and false otherwise.
        + coreAlpha        ::Float64             ... Polarizibility of the ionic core.
        + coreRadius       ::Float64             ... Radius of the ionic core.
        + coreShells       ::Array{Shell,1}      ... Shells of the ionic core
    """
    struct  CorePolarization
        doApply            ::Bool     
        coreAlpha          ::Float64 
        coreRadius         ::Float64
        coreShells         ::Array{Shell,1} 
    end 


    """
    `Basics.CorePolarization()`  ... constructor for the default values of core-polarization contributions.
    """
    function CorePolarization()
        CorePolarization(false, 0., 0., Shell[])
    end



    # `Base.show(io::IO, cp::CorePolarization)`  ... prepares a proper printout of the variable cp::CorePolarization.
    function Base.show(io::IO, cp::CorePolarization) 
        sa = Base.string(cp);                print(io, sa)
    end


    # `Base.string(cp::CorePolarization)`  ... provides a String notation for the variable cp::CorePolarization.
    function Base.string(cp::CorePolarization)
        if  cp.doApply  sa = "Apply "   else    sa = "Do not apply "     end
        sa = sa * "core-polarization with alpha_c = $(cp.coreAlpha) a.u., r_c = $(cp.coreRadius) a.u. " *
                  "core shells = $(cp.coreShells)."
        return( sa )
    end

    
    """
    `abstract type Basics.AbstractContinuumNormalization` 
        ... defines an abstract and a number of singleton types for dealing with the normalization of continuum orbitals.

        + PureSineNorm       ... normalize with regard to an (asymtotic) pure sine funtion, sin(kr).
        + CoulombSineNorm    ... normalize with regard to an (asymtotic) Coulombic-sine funtion, sin(kr + ...).
        + OngRussekNorm      ... normalize by following Ong & Russek (1973).
    """
    abstract type  ContinuumNormalization                          end
    struct     PureSineNorm         <:  ContinuumNormalization     end
    struct     CoulombSineNorm      <:  ContinuumNormalization     end
    struct     OngRussekNorm        <:  ContinuumNormalization     end

    export  AbstractContinuumNormalization,   PureSineNorm,   CoulombSineNorm,   OngRussekNorm

    
    """
    `abstract type Basics.AbstractContinuumSolutions` 
        ... defines an abstract and a number of singleton types for solving the continuum orbitals in a given potential.

        + ContBessel              ... generate a pure Bessel function for the large component together with kinetic balance.
        + ContSine                ... generate a pure Sine function for the large component together with kinetic balance.
        + NonrelativisticCoulomb  ... generate a non-relativistic Coulomb function for the large component together with kinetic balance.
        + AsymptoticCoulomb       ... generate a pure (asymptotic) Coulombic function for both components.
        + BsplineGalerkin         ... generate a continuum orbital with the Galerkin method.dealing with warnings that are made during a run or REPL session.
    """
    abstract type  AbstractContinuumSolutions                            end
    struct     ContBessel             <:  AbstractContinuumSolutions     end
    struct     ContSine               <:  AbstractContinuumSolutions     end
    struct     NonrelativisticCoulomb <:  AbstractContinuumSolutions     end
    struct     AsymptoticCoulomb      <:  AbstractContinuumSolutions     end
    struct     BsplineGalerkin        <:  AbstractContinuumSolutions     end

    export  AbstractContinuumSolutions, ContBessel, ContSine, NonrelativisticCoulomb, AsymptoticCoulomb, BsplineGalerkin 

    
    """
    `abstract type Basics.AbstractWarning` 
        ... defines an abstract and a number of singleton types for dealing with warnings that are made during a run or REPL session.
            Cf. Defaults.warn().

        + AddWarning        ... add a Warning to a warningList.
        + PrintWarnings     ... print all warnings into a jac-warn.report file.
        + ResetWarnings     ... reset (empty) the warningList, usually at the beginning of a new run.to distinguish between different warnings
    """
    abstract type  AbstractWarning                          end
    struct     AddWarning           <:  AbstractWarning     end
    struct     PrintWarnings        <:  AbstractWarning     end
    struct     ResetWarnings        <:  AbstractWarning     end

    export  AbstractWarning, AddWarning, PrintWarnings, ResetWarnings 
    
    
    """
    `abstract type Basics.AbstractExcitationScheme` 
        ... defines an abstract and a number of singleton types to distinguish between different schemes for
            generating configuration lists as they frequently occur in Green function and cascade computations.

      + struct NoExcitationScheme        
        ... dummy scheme for (unsupported) initialization of this abstract tpye.

      + struct DeExciteSingleElectron        
        ... generates all excitations and de-excitations of a single electron from a given list of bound
            electron configurations. The number of electrons of the generated configurations is the same as 
            for the given bound configurations.

      + struct DeExciteTwoElectrons       
        ... generates all excitations and de-excitations of one or two electrons from a given list of bound
            electron configurations. The number of electrons of the generated configurations is the same as 
            for the given bound configurations.
            
      + struct AddSingleElectron             
        ... generates configurations by just adding a single electrons to a given list of bound
            electron configurations. The number of electrons of the generated configurations is N+1.
            
      + struct ExciteByCapture             
        ... generates all excitations and de-excitations of one or more electron from a given list of bound
            electron configurations, together with an capture of an additional electron. 
            The number of electrons of the generated configurations is N+1.
    """
    abstract type  AbstractExcitationScheme                               end
    struct         NoExcitationScheme      <:  AbstractExcitationScheme   end
    struct         DeExciteSingleElectron  <:  AbstractExcitationScheme   end
    struct         DeExciteTwoElectrons    <:  AbstractExcitationScheme   end
    struct         AddSingleElectron       <:  AbstractExcitationScheme   end
    struct         ExciteByCapture         <:  AbstractExcitationScheme   end

    
    """
    `abstract type Basics.AbstractScField` 
        ... defines an abstract and a number of singleton types to distinguish between different self-consistent fields

        + struct ALField          ... to represent an average-level field.        
        + struct EOLField         ... to represent an (extended) optimized-level field.        
        + struct DFSField         ... to represent an mean Dirac-Fock-Slater field.        
        + struct DFSwCPField      ... to represent an mean Dirac-Fock-Slater with core-polarization field.        
        + struct HSField          ... to represent an mean Hartree-Slater field.        
        + struct NuclearField     ... to represent a pure nuclear (potential) field.        
    """
    abstract type  AbstractScField                          end
    struct     ALField              <:  AbstractScField     end
    struct     EOLField             <:  AbstractScField     end
    struct     DFSField             <:  AbstractScField     end
    struct     HSField              <:  AbstractScField     end
    struct     NuclearField         <:  AbstractScField     end

    """
    `struct  Basics.DFSwCPField          <:  AbstractScField`  
        ... defines a type to describe a mean Dirac-Fock-Slater field with core-polarization.

        + corePolarization  ::CorePolarization   ... Parametrization of the core-polarization potential/contribution.
    """
    struct     DFSwCPField          <:  AbstractScField   
        corePolarization    ::CorePolarization
    end

    export  AbstractScField, ALField, EOLField, DFSField, DFSwCPField, HSField, NuclearField

    
    
    """
    `abstract type Basics.AbstractEeInteraction` 
        ... defines an abstract and a number of singleton types for specifying the electron-electron interaction.

        + struct CoulombInteraction   ... to represent the Coulomb part of the electron-electron interaction.        
        + struct BreitInteraction     ... to represent the Breit part of the electron-electron interaction.        
        + struct CoulombBreit         ... to represent the Coulomb+Breit part of the electron-electron interaction.        
    """
    abstract type  AbstractEeInteraction                          end
    struct     CoulombInteraction   <:  AbstractEeInteraction     end
    struct     BreitInteraction     <:  AbstractEeInteraction     end
    struct     CoulombBreit         <:  AbstractEeInteraction     end
    
    export  AbstractEeInteraction, CoulombInteraction, BreitInteraction, CoulombBreit

    
    #==
    """
    `abstract type Basics.AbstractCImethod` 
        ... defines an abstract and a number of singleton types to determine CI matrix solutions.

        + struct FullCIjulia          ... to diagonalize the full CI matrix with the internal Julia function.        
    """
    abstract type  AbstractCImethod                          end
    struct     FullCIeigen          <:  AbstractCImethod     end
    
    export  AbstractCImethod, FullCIeigen  ==#

    """
    `abstract type Basics.AbstractPotential` 
        ... defines an abstract and a number of singleton types to distinguish between different (electronic) atomic potentials.

        + struct DFSpotential     ... to represent a Dirac-Fock-Slater potential.        
        + struct CoreHartree      ... to represent a core-Hartree potential.        
        + struct KohnSham         ... to represent a Kohn-Sham potential.        
        + struct HartreeSlater    ... to represent a Hartree-Slater potential.        
    """
    abstract type  AbstractPotential                        end
    struct    DFSpotential          <:  AbstractPotential   end
    struct    CoreHartree           <:  AbstractPotential   end
    struct    KohnSham              <:  AbstractPotential   end
    struct    HartreeSlater         <:  AbstractPotential   end
    
    export  AbstractPotential, DFSpotential, CoreHartree, KohnSham, HartreeSlater

    
    """
    `abstract type Basics.AbstractPropertySettings` 
        ... defines an abstract type to distinguish between different settings of atomic level properties.
    """
    abstract type  AbstractPropertySettings                 end

    
    """
    `abstract type Basics.AbstractProcessSettings` 
        ... defines an abstract type to distinguish between different settings of atomic processes.
    """
    abstract type  AbstractProcessSettings                  end
    
    export  AbstractPropertySettings, AbstractProcessSettings

    
    #======================================     Remove abstract idenitfier and only use settings to specify properties and processes
    """
    `abstract type Basics.AbstractLevelProperty` 
        ... defines an abstract and a number of singleton types to distinguish between different (electronic) atomic potentials.

        + struct NoProperty       ... No level property defined.        
        + struct HFS              ... Hyperfine A and B parameters.     
        + struct Isotope          ... Isotope shift M and F parameters.      
        + struct LandeJ           ... Lande g_J factors.
        + struct LandeF           ... Lande g_F factors.   
        + struct Polarizibility   ... static and dynamic polarizibilities. 
        + struct Plasma           ... CI computations including interactions from various plasma models. 
        + struct Zeeman           ... Zeeman splitting of fine-structure levels.
        + struct EinsteinX        
            ... Einstein A, B coefficients and oscillator strength; although not a 'level property', in the Einstein module these computations 
                are treated within a single basis and without all relaxation effects, etc. The `Einstein' property therefore helps to obtain a 
                quick overview about transition probabilities of a transition arry or if many of these Einstein coefficients need to be 
                calculated for a cascade.   
    """
    abstract type  AbstractLevelProperty                         end
    struct    NoProperty            <:  AbstractLevelProperty    end
    struct    AlphaX                <:  AbstractLevelProperty    end
    struct    FormF                 <:  AbstractLevelProperty    end
    struct    HFS                   <:  AbstractLevelProperty    end
    struct    Isotope               <:  AbstractLevelProperty    end
    struct    LandeJ                <:  AbstractLevelProperty    end
    struct    LandeF                <:  AbstractLevelProperty    end
    struct    Plasma                <:  AbstractLevelProperty    end
    struct    Polarizibility        <:  AbstractLevelProperty    end
    struct    Yields                <:  AbstractLevelProperty    end
    struct    Zeeman                <:  AbstractLevelProperty    end
    struct    EinsteinX             <:  AbstractLevelProperty    end

    
    export  AbstractProperty,  NoProperty,  AlphaX,  EinsteinX,  FormF,  HFS,  Isotope,  LandeJ,  LandeF, Plasma,  
            Polarizibility,  Yields,  Zeeman

    function Base.string(prop::NoProperty)      return( "no property" )                                     end
    function Base.string(prop::AlphaX)          return( "alpha variation" )                                 end
    function Base.string(prop::EinsteinX)       return( "Einstein A, B coefficients" )                      end  
    function Base.string(prop::FormF)           return( "form & scattering factors" )                       end
    function Base.string(prop::HFS)             return( "HFS A, B parameters" )                             end
    function Base.string(prop::Isotope)         return( "isotope shift M and F parameters" )                end
    function Base.string(prop::LandeJ)          return( "Lande g_J factors" )                               end
    function Base.string(prop::LandeF)          return( "Lande g_F factors" )                               end
    function Base.string(prop::Plasma)          return( "CI computations with plasma-type interactions" )   end
    function Base.string(prop::Polarizibility)  return( "multipole polarizibility" )                        end
    function Base.string(prop::Yields)          return( "fluorescence & Auger yields" )                     end
    function Base.string(prop::Zeeman)          return( "Zeeman splitting in an external magnetic field" )  end 
    ================================================#    

    
    """
    `abstract type Basics.AbstractProcess` 
        ... defines an abstract and a number of singleton types to distinguish different atomic processes.

        + struct Auger            ... Auger transitions, i.e. single autoionization or the emission of a single free electron into the continuum.
        + struct AugerInPlasma    ... Auger transitions but calculated for a specified plasma model.
        + struct Compton          ... Rayleigh-Compton scattering cross sections.
        + struct Coulex           ... Coulomb-excitation of target or projeticle electrons by fast, heavy ions.
        + struct Coulion          ... Coulomb-ionization of target or projeticle electrons by fast, heavy ions.
        + struct Dierec           ... di-electronic recombination, i.e. the dielectronic capture of a free electron and the subsequent emission of a photon.
        + struct DoubleAuger      ... Double Auger rates.
        + struct ImpactExcAuto    ... di-electronic recombination, i.e. the dielectronic capture of a free electron and the subsequent emission of a photon.
        + struct MultiPhotonDE    ... multi-photon excitation and decay rates, including 2-photon, etc. processes.
        + struct MultiPI          ... multi-photon (single-electron) ionization.
        + struct MultiPDI         ... multi-photon (single-electron) double ionization.
        + struct Photo            ... Photoionization processes, i.e. the emission of a single free electron into the continuum due to an external light field.
        + struct PhotoDouble      ... Photo-double ionization rates.
        + struct PhotoExc         ... Photoexcitation rates.
        + struct PhotoExcFluor    ... photoexcitation fluorescence rates and cross sections.
        + struct PhotoExcAuto     ... photoexcitation autoionization cross sections and collision strengths.
        + struct PhotoInPlasma    ... Photoionization processes but calculated for a specified plasma model.
        + struct PhotoIonFluor    ... photoionization fluorescence rates and cross sections.
        + struct PhotoIonAuto     ... photoionization autoionization cross sections and collision strengths.
        + struct Radiative        ... Radiative (multipole) transitions between bound-state levels of the same charge state.
        + struct Rec              ... radiative electron capture, i.e. the capture of a free electron with the simultaneous emission of a photon.
        + struct Eimex            ... electron-impact excitation cross sections and collision strengths.
        + struct RAuger           ... Radiative Auger rates.
    """
    abstract type  AbstractProcess                          end
    struct    NoProcess             <:  AbstractProcess     end
    struct    Auger                 <:  AbstractProcess     end
    struct    AugerInPlasma         <:  AbstractProcess     end
    struct    Compton               <:  AbstractProcess     end
    struct    Coulex                <:  AbstractProcess     end
    struct    Coulion               <:  AbstractProcess     end
    struct    Dierec                <:  AbstractProcess     end
    struct    DoubleAuger           <:  AbstractProcess     end
    struct    ElecCapture           <:  AbstractProcess     end
    struct    ImpactExcAuto         <:  AbstractProcess     end
    struct    InternalConv          <:  AbstractProcess     end
    struct    MultiPhotonDE         <:  AbstractProcess     end
    struct    MultiPI               <:  AbstractProcess     end
    struct    MultiPDI              <:  AbstractProcess     end
    struct    Photo                 <:  AbstractProcess     end
    struct    PhotoDouble           <:  AbstractProcess     end
    struct    PhotoExc              <:  AbstractProcess     end
    struct    PhotoExcFluor         <:  AbstractProcess     end
    struct    PhotoExcAuto          <:  AbstractProcess     end
    struct    PhotoInPlasma         <:  AbstractProcess     end
    struct    PhotoIonFluor         <:  AbstractProcess     end
    struct    PhotoIonAuto          <:  AbstractProcess     end
    struct    Radiative             <:  AbstractProcess     end
    struct    Rec                   <:  AbstractProcess     end
    struct    Eimex                 <:  AbstractProcess     end
    struct    RAuger                <:  AbstractProcess     end
    struct    PairA1P               <:  AbstractProcess     end

    export  NoProcess, Auger, AugerInPlasma, Compton, Coulex, Coulion, Dierec, DoubleAuger, Eimex, ElecCapture, ImpactExcAuto, InternalConv, 
            MultiPhotonDE, MultiPI, MultiPDI, Photo, PhotoDouble, PhotoExc, PhotoExcAuto, PhotoExcFluor, PhotoInPlasma, PhotoIonAuto, 
            PhotoIonFluor, Radiative, RAuger, Rec, PairA1P, Coulion

    function Base.string(propc::NoProcess)          return( "no process" )                         end
    function Base.string(propc::Auger)              return( "Auger" )                              end
    function Base.string(propc::AugerInPlasma)      return( "Auger in plasma" )                    end
    function Base.string(propc::Compton)            return( "Rayleigh-Compton" )                   end
    function Base.string(propc::Dierec)             return( "Dielectronic recombination" )         end
    function Base.string(propc::DoubleAuger)        return( "Double Auger" )                       end
    function Base.string(propc::ElecCapture)        return( "Electron capture" )                   end
    function Base.string(propc::ImpactExcAuto)      return( "ImpactExcAuto" )                      end
    function Base.string(propc::InternalConv)       return( "InternalConv" )                       end
    function Base.string(propc::MultiPhotonDE)      return( "multi-photon excitation & decay" )    end
    function Base.string(propc::Photo)              return( "Photo" )                              end
    function Base.string(propc::PhotoDouble)        return( "single-photon double ionization" )    end
    function Base.string(propc::PhotoExcFluor)      return( "Photo-Excitation-Fluoresence" )       end
    function Base.string(propc::PhotoExcAuto)       return( "Photo-Excitation-Autoionization" )    end
    function Base.string(propc::PhotoInPlasma)      return( "Photo in plasma" )                    end
    function Base.string(propc::PhotoIonFluor)      return( "Photo-Ionization-Fluoresence" )       end
    function Base.string(propc::PhotoIonAuto)       return( "Photo-Ionization-Autoionization" )    end
    function Base.string(propc::Radiative)          return( "Radiative" )                          end
    function Base.string(propc::RAuger)             return( "Radiative Auger" )                    end
    function Base.string(propc::Rec)                return( "Rec" )                                end
    
    """
    `abstract type Basics.AbstractSpectrumKind` 
        ... defines an abstract and a number of singleton types for specifying different kinds of spectra that need to be displayed

        + struct BarIntensities                         
            ... to display just intensity-bars at given x-positions.  
        + struct DiscreteLines                         
            ... to display lines (for guiding the eyes) but which are only defined at discrete x-points.
        + struct DiscretePoints                         
            ... to display discrete points that are defined at discrete x-points.
        + struct LorentzianIntensitiesConstantWidths    
            ... to display the superposition of Lorentzians with given position and intensity but constant widths.        
        + struct LorentzianIntensities  
            ... to display the superposition of Lorentzians with given position, intensity and individual widths.     
    """
    abstract type  AbstractSpectrumKind                                          end
    struct     BarIntensities                       <:  AbstractSpectrumKind     end
    struct     DiscreteLines                        <:  AbstractSpectrumKind     end
    struct     DiscretePoints                       <:  AbstractSpectrumKind     end
    struct     LorentzianIntensitiesConstantWidths  <:  AbstractSpectrumKind     end
    struct     LorentzianIntensities                <:  AbstractSpectrumKind     end
    
    export  AbstractSpectrumKind, BarIntensities, DiscreteLines, DiscretePoints, LorentzianIntensitiesConstantWidths, LorentzianIntensities

            
            
    # Functions/methods that are later added to the module Basics
    function add                                                    end
    function addZerosToCsfR                                         end
    function analyze                                                end
    function analyzeConvergence                                     end
    function compute                                                end
    function computeDensity                                         end
    function computeDiracEnergy                                     end
    function computeMeanSubshellOccupation                          end
    function computeMultipletForGreenApproach                       end
    function computePotentialCoreHartree                            end
    function computePotentialHartree                                end
    function computePotentialHartreeSlater                          end
    function computePotentialKohnSham                               end
    function computePotentialDFS                                    end
    function computePotentialDFSwCP                                 end
    function computePotentialExtendedHartree                        end
    function computeScfCoefficients                                 end
    function determineEnergySharings                                end
    function determineHoleShells                                    end
    function determineMeanEnergy                                    end
    function determineParity                                        end
    function determinePolarizationLambda                            end
    function determinePolarizationVector                            end
    function determineSelectedLines                                 end
    function determineSelectedPathways                              end
    function diagonalize                                            end
    function diracDelta                                             end
    function display                                                end
    function displayLevels                                          end
    function displayOrbitalOverlap                                  end
    function excludeDoubles                                         end
    function extractLeadingConfiguration                            end
    function extractNoOpenShells                                    end
    function extractNonrelativisticShellList                        end
    function extractNonrelativisticConfigurations                   end
    function extractNonrelativisticConfigurationFromCsfR            end
    function extractOpenShellQNfromCsfNR                            end
    function extractOpenShellQNfromCsfR                             end
    function extractRelativisticSubshellList                        end
    function extractShellOccupationFromCsfR                         end
    function generate()                                             end
    function generateBasis                                          end
    function generateConfigurations                                 end
    function generateConfigurationsForExcitationScheme              end
    function generateConfigurationsWithElectronCapture              end
    function generateLevelWithExtraElectron                         end
    function generateLevelWithExtraSubshell                         end
    function generateLevelWithSymmetryReducedBasis                  end
    function generateOrbitalsForBasis                               end
    function generateOrbitalsForPotential                           end
    function generateOrbitalSuperposition                           end
    function generateShellList                                      end
    function generateSubshellList                                   end
    function generateSpectrumLorentzian                             end
    function generateSpectrumGaussian                               end
    function integrate                                              end
    function integrateOnGridNewtonCotes                             end
    function integrateOnGridSimpsonRule                             end
    function integrateOnGridTrapezRule                              end
    function interpolateOnGridGrasp92                               end
    function interpolateOnGridTrapezRule                            end
    function isSimilar                                              end
    function isSymmetric                                            end
    function isStandardSubshellList                                 end
    function isViolated                                             end
    function isZero                                                 end
    function lastPoint                                              end
    function merge                                                  end
    function modifyLevelEnergies                                    end
    function perform                                                end
    function performSCF                                             end
    function performCI                                              end
    function plot                                                   end
    function read                                                   end
    function readCslFileGrasp92                                     end
    function readOrbitalFileGrasp92                                 end
    function readMixFileRelci                                       end
    function readMixingFileGrasp92                                  end
    function recast                                                 end
    function selectLevel                                            end
    function selectLevelPair                                        end
    function selectLevelTriple                                      end
    function selectSymmetry                                         end
    function shiftTotalEnergies                                     end
    function sortByEnergy                                           end
    function tabulate                                               end
    function tabulateKappaSymmetryEnergiesDirac                     end
    function tools                                                  end
    function yesno                                                  end
    
end # module

