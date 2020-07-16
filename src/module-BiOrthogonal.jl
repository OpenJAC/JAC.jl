
"""
`module  JAC.BiOrthogonal`  
    ... a submodel of JAC that contains all methods for a bi-orthogonal transformation of two non-orthogonal orbitals set.
"""
module BiOrthogonal

    using Printf, JAC, ..Basics,  ..Defaults, ..ManyElectron, ..Radial



    """
    `BiOrthogonal.computeTransformation(leftMp::Multiplet, rightMp::Multiplet)`  
        ... computes the bi-orthogonal transformation of the two given multiplets by rotating the radial functions 
            by counter-rotating the corresponding CI coefficients. A tuple of two multiplets
            tpl(newLeftMp::::Multiplet, newRightMp::Multiplet) is returned.
    """
    function  computeTransformation(leftMp::Multiplet, rightMp::Multiplet)
        return( 1.0 )
    end
    
    

    """
    `BiOrthogonal.computeTransformationMatrices()`  
        ... computes the transformation matrices for rotating the radial functions and for counter-rotating CI coefficients.
    """
    function  computeTransformationMatrices()
        return( 1.0 )
    end
    
    

    """
    `BiOrthogonal.generateBiorthogonalShellMatrices()`  
        ... generates the biorthogonal forms of initial and final shells.
    """
    function  generateBiorthogonalShellMatrices()
        return( 1.0 )
    end
    
    

    """
    `BiOrthogonal.generateCounterRotatingCiMatrices()`  
        ... generates the counter-rotating matrices for the CI coefficients of the initial and final states.
    """
    function  generateCounterRotatingCiMatrices()
        return( 1.0 )
    end

end # module


#========================================================
Comments & thoughts on the implementation of this module
--------------------------------------------------------

++  I looked briefly through your rbiotransform implementation in GRASP2018 to get a very first impression which basic ingredients 
    are needed for such an implementation. A large fraction of code refers to 'internal data handling of quantum numbers and 
    functions' which might become obsolete or can be coded much more compact in Julia.
++  Moreover, for the manipulation of matrices (multiplication, inversion, ...), etc., we shall use the internal Julia features 
    as much as possible. This will likely reduced the size of the code considerably since some of the procedures in rbiotransform
    can be written as a single line.

++  Overall, therefore, I belief that the code could be made more compact by a factor 3-7, and also more easy to check.

++  A first goal should be to write-down the basic formulas for and in the notation of the UserGuide of JAC.
    This will enable us to bring the notation and internal use much closer to each other; I noticed that large parts in the
    rbiotransform implementation in GRASP2018 are highly technical with very little use of the underlying physics language
    which is required to formulate the transformation itself.
   
++  I will be happy (and grateful) if we can agree about proper names of functions and if we make use of the 
    standard Julia features.

++  Please, use ? AngularMomentum.Wigner_6j(a, b, c, d, e, f) and similar for Wigner_3j and _9j symbols to get the Wigner symbols
    if needed.

    
++  Let's use CamelCase notation an start all variables with a lowercase letter; arrays are usually given the same name
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

