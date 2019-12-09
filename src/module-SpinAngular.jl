
"""
`module  JAC.SpinAngular`  
    ... a submodel of JAC that contains all methods for computing the spin-angular coefficients for one- and two-particle operators 
        of given rank and for symmetry-adapted CSF. At the present, these spin-angular can be obtained for all scalar (and symmetric) 
        one and two-particle interactions and they are sometimes refered to as 'pure angular coefficients' in order to distinguish 
        them from those used in GRASP (and its derivatives) where part of the `physical interaction' were included orginal 
        into the angular integrals. The spin-angular coefficients usally appear in the form:
        
            <CSF_l | operator | CSF_r> = sum_{t} coeff_t (L,abcd) * X^L (abcd)
            
        where CSF is a standard (jj-coupled) configuration state function, t a summation index, and where X^L(abcd) denotes the 
        effective interaction strength. 
"""
module SpinAngular

    using Printf, JAC, ..Basics,  ..Defaults, ..ManyElectron, ..Radial
    
    export  Xronecker

    
    """
    `abstract type SpinAngular.AbstractAngularType` 
        ... defines an abstract type and a number of data types to work with one- and two-particle operators of given rank, see also:
        
        + struct OneParticleOperator    ... to represent a single-particle operator with well-defined spherical tensor rank.
        + struct TwoParticleOperator    ... to represent a two-particle operator with well-defined spherical tensor rank.
    """
    abstract type  AbstractAngularType                           end

    
    """
    `struct  SpinAngular.OneParticleOperator  <:  AbstractAngularType`  
        ... a struct for defining the spherial tensor structure of a single-particle operator.

        + rank            ::Int64             ... Rank of the operator.
        + parity          ::Basics.Parity     ... Parity of the operator (if needed ??)
        + sameOrbitalSet  ::Bool              ... True if orbitals for both CSF are taken from the same orbital set (if needed ??) 
    """
    struct  OneParticleOperator  <:  AbstractAngularType 
        rank              ::Int64         
        parity            ::Basics.Parity 
        sameOrbitalSet    ::Bool
    end

    
    """
    `SpinAngular.OneParticleOperator()`  ... constructor for setting the default values.
    """
    function OneParticleOperator()
    	OneParticleOperator( 0, Basics.plus, false)
    end
    
    
    # `Base.show(io::IO, op::OneParticleOperator)`  ... prepares a proper printout of the op::OneParticleOperator.
    function Base.show(io::IO, op::OneParticleOperator)
          sa = "One-particle operator O^($(op.rank) $(op.parity))"
          if op.sameOrbitalSet  sa = sa * ", defined for CSF from the same orbital set."
          else                  sa = sa * ", defined for CSF from different orbital sets."
          end
    	  println(io, sa)
    end

    
    """
    `struct  SpinAngular.TwoParticleOperator`  
        ... a struct for defining the spherial tensor structure of a two-particle operator.

        + rank            ::Int64             ... Rank of the operator (if needed ??).
        + parity          ::Basics.Parity     ... Parity of the operator (if needed ??)
        + sameOrbitalSet  ::Bool              ... True if orbitals for both CSF are taken from the same orbital set (if needed ??) 
    """
    struct  TwoParticleOperator
        rank              ::Int64         
        parity            ::Basics.Parity 
        sameOrbitalSet    ::Bool
    end

    
    """
    `SpinAngular.TwoParticleOperator()`  ... constructor for setting the default values.
    """
    function TwoParticleOperator()
    	TwoParticleOperator( 0, Basics.plus, false)
    end

    
    # `Base.show(io::IO, op::TwoParticleOperator)`  ... prepares a proper printout of the op::TwoParticleOperator.
    function Base.show(io::IO, op::TwoParticleOperator)
          sa = "Two-particle operator O^($(op.rank) $(op.parity))"
          if op.sameOrbitalSet  sa = sa * ", defined for CSF from the same orbital set."
          else                  sa = sa * ", defined for CSF from different orbital sets."
          end
    	  println(io, sa)
    end

    
    """
    `struct  SpinAngular.Coefficient1p`  
        ... a struct for defining a single spin-angular coefficient for a reduced one-particle matrix element <a || o^(L) || b>.

        + L        ::Int64             ... Rank (L or nu) of the single-particle interaction strength.
        + a        ::Subshell          ... Left-hand subshell (orbital).
        + b        ::Subshell          ... Right-hand subshell (orbital).
        + v        ::Float64           ... (Value of) spin-angular coefficient.
    """
    struct  Coefficient1p
        L          ::Int64    
        a          ::Subshell 
        b          ::Subshell
        v          ::Float64
    end

    
    # `Base.show(io::IO, coeff::Coefficient1p)`  ... prepares a proper printout of the coeff::Coefficient1p.
    function Base.show(io::IO, coeff::Coefficient1p)
          sa = "V^($(coeff.L)) [$(coeff.a), $(coeff.b)] = $(coeff.v)"
    	  println(io, sa)
    end

    
    """
    `struct  SpinAngular.Coefficient2p`  
        ... a struct for defining a single spin-angular coefficient for a reduced two-particle matrix element <ab || o^(L) || cd>,
            such as the Coulomb interaction strength.

        + L        ::Int64             ... Rank (L or nu) of the single-particle interaction strength.
        + a        ::Subshell          ... Left-hand subshell (orbital).
        + b        ::Subshell          ... Left-hand subshell (orbital).
        + c        ::Subshell          ... Right-hand subshell (orbital).
        + d        ::Subshell          ... Right-hand subshell (orbital).
        + v        ::Float64           ... (Value of) spin-angular coefficient.
    """
    struct  Coefficient2p
        L          ::Int64    
        a          ::Subshell 
        b          ::Subshell
        c          ::Subshell 
        d          ::Subshell  
        v          ::Float64
    end

    
    # `Base.show(io::IO, coeff::Coefficient2p)`  ... prepares a proper printout of the coeff::Coefficient2p.
    function Base.show(io::IO, coeff::Coefficient2p)
          sa = "V^($(coeff.L)) [$(coeff.a), $(coeff.b)| $(coeff.c), $(coeff.d)] = $(coeff.v)"
    	  println(io, sa)
    end


    """
    `SpinAngular.computeCoefficients(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR)`  
        ... computes the spin-angular coefficients for the reduced one-particle matrix element
            <leftCsf || op || rightCsf >; a list of one-particle coefficients coeffs::Array{Coefficient1p,1} is returned.
    """
    function  computeCoefficients(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR)
        a = Subshell("1s_1/2");     b = Subshell("2p_1/2")
        coeff  = Coefficient1p( 1, a, b, 1.0 )
        coeffs = [coeff]
        return( coeffs )
    end


    """
    `SpinAngular.computeCoefficients(op::SpinAngular.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR)`  
        ... computes the spin-angular coefficients for the reduced otwo-particle matrix element
            <leftCsf || op || rightCsf >; a Tuple of two lists with one- and two-particle coefficients 
            tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}}  is returned.
    """
    function  computeCoefficients(op::SpinAngular.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR)
        a = Subshell("1s_1/2");     b = Subshell("2p_1/2")
        coeff1p  = Coefficient1p( 0, a, b, 2.0 );           coeffs1p = [coeff1p]
        coeff2p  = Coefficient2p( 0, a, b, a, b, 2.0 );     coeffs2p = [coeff2p]
        return( (coeffs1p, coeffs2p) )
    end

    
    

    """
    `SpinAngular.computeRecouplingMatrixFor2Shells(..., leftCsf::CsfR, rightCsf::CsfR)`  
        ... computes the recoupling matrix for two open shells ... 
            all arguments should be clearly explained, please
    """
    function  computeCoefficients(leftCsf::CsfR, rightCsf::CsfR)
        return( 1.0 )
    end


end # module


#========================================================
Comments & thoughts on the implementation of this module
--------------------------------------------------------

++  I looked briefly through your ANO/RATIP implementation to remember the basic ingredients of such an implementation 
    and to suggest some first functions which you likely need in 'this' or a rather similar form.
++  Overall, however, I belief that the overall code should be a factor 3-5 shorter and more easy to check.
++  Above you find the basic data structures which will allow me/JAC to communicate with this SpinAngular module.

++  In the first instance, I'm interested only in the spin-angular coefficients that belong to a single reduced
    many-electron matrix element ... I do currently not plan to benefit from the simultaneous computation of a whole
    raw or matrix of (matrix) elements. 
   
++  I will be happy (and grateful) if we can agree about proper names of functions and if we make use of the 
    standard Julia features; names like  anco_case_6 or similar are no longer appropriate and will make you own
    life difficult. 
    Let's introduce a AbstractAngularCase if really necessary.

++  Please, use ? AngularMomentum.Wigner_6j(a, b, c, d, e, f) and similar for Wigner_3j and _9j symbols to get the Wigner symbols
    if needed.
    
    

++  Please, use CamelCase notation an start all variables with a lowercase letter; arrays are usually given the same name
    with an additional s:  coeff::Coefficient1p  vs.  coeffs::Array{Coefficient1p,1}
++  All angular momenta should be of type j::AngularJ64 (for j-values) and of type m::AngularM64 (for m-values, 
    and if they appear at all in your implementation), please.
++  For a given subhsell sh::Subshell, get n and kappa simply by: sh.n, sh.kappa
++  For a given subhsell sh::Subshell, get l and j simply by:     Basics.subshell_l(sh), Basics.subshell_j(sh)

++  Look for ?CsfR ... to get the underlying definition of a CSF in JAC; it allows a direct access to J, parity,
    the occupation, seniority and all angular momenta from the coupling. 
    
++  if csf.useStandardSubshells  subshells = Defaults.GBL_STANDARD_SUBSHELL_LIST
    else                         subshells = csf.subshells    
    end
    In the latter case, these subshells could be different for the leftCsf and rightCsf.

========================================================#
