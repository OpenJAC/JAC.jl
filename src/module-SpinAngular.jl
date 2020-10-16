
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

    using  Printf, ..Basics,  ..Defaults, ..ManyElectron, ..Radial
    
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
          sa = "\n V^($(coeff.L)) [$(coeff.a), $(coeff.b)] = $(coeff.v)"
    	  print(io, sa)
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
          sa = "\n V^($(coeff.L)) [$(coeff.a), $(coeff.b)| $(coeff.c), $(coeff.d)] = $(coeff.v)"
    	  print(io, sa)
    end

    
    """
    `@enum   SchemeEta`  
        ... defines a enumeration for the allowed values of the Eta scheme to distinguish the coupling and selection of the
            tensor operators.

        + Eta_a         ... a^(qj)_mq
        + Eta_W         ... W^(k_12) = [a x a]
        + Eta_aW        ... W^(k_12, k_2) = [a1 x [a2 x a3]]
        + Eta_Wa        ... W^(k_12, k_2) = [[a1 x a2] x a3]
        + Eta_WW        ... W^(kk0) = [[a1 x a2]^k x [a3 x a4]^k]^0
    """
    @enum   SchemeEta    Eta_none   Eta_a   Eta_W   Eta_aW   Eta_Wa   Eta_WW
    
    
    """
    `struct  SpinAngular.SchemeGamma??`  
        ... a singleton struct to distinguish between different coupling schemes of the matrix elements.
    """
    struct  SchemeGamma01   end
    struct  SchemeGamma02   end
    struct  SchemeGamma03   end

    
    
    """
    `struct  SpinAngular.SubshellTerm`  
        ... a struct for defining a subshell term/state  |j (nu) alpha Q J> == |j (nu) Q J Nr> for a subshell with well-defined j.

        + j          ::AngularJ64      ... subshell j
        + Q          ::AngularJ64      ... quasi-spin
        + seniority  ::Int64           ... seniority
        + J          ::AngularJ64      ... total J of subshell term
        + Nr         ::Int64           ... Additional quantum number Nr = 0,1,2.
    """
    struct  SubshellTerm
        j            ::AngularJ64
        Q            ::AngularJ64
        seniority    ::Int64
        J            ::AngularJ64
        Nr           ::Int64
    end

    
    # `Base.show(io::IO, term::SubshellTerm)`  ... prepares a proper printout of term::SubshellTerm.
    function Base.show(io::IO, term::SubshellTerm)
        if  term.Nr == 0   sa = "|$(term.j) ($(term.seniority)) $(term.Q) $(term.J)>"
        else               sa = "|$(term.j) ($(term.seniority)) $(term.Q) $(term.J); Nr=$(term.Nr)>"
        end
        println(io, sa)
    end

    
    """
    `struct  SpinAngular.Lambda`  
        ... a struct for defining a subshell term/state  |j (nu) alpha Q J> == |j (nu) Q J Nr> for a subshell with well-defined j.

        + Ji         ::AngularJ64      ... J_i
        + Jj         ::AngularJ64      ... J_j
        + Jip        ::AngularJ64      ... J_i'
        + Jjp        ::AngularJ64      ... J_j'
    """
    struct  Lambda
        Ji           ::AngularJ64
        Jj           ::AngularJ64
        Jip          ::AngularJ64
        Jjp          ::AngularJ64
    end

    
    # `Base.show(io::IO, lambda::Lambda)`  ... prepares a proper printout of lambda::Lambda.
    function Base.show(io::IO, lambda::Lambda)
        sa = "Lambda($(lambda.Ji), $(lambda.Jj); lambda.Jip), $(lambda.Jjp))"
        println(io, sa)
    end
    
    include("module-SpinAngular-inc-reducedcoeffs.jl")
    
    
    
    #######################################################################################################################
    #######################################################################################################################


    """
    `SpinAngular.computeCoefficients(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
        ... computes the spin-angular coefficients for the reduced one-particle matrix element
            <leftCsf || op^(k) || rightCsf >  if both CSF refer to the same list of subshells; a list of one-particle 
            coefficients coeffs::Array{Coefficient1p,1} is returned.
    """
    function  computeCoefficients(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
        if      op.rank == 0   coeffs = computeCoefficientsScalar(op,    leftCsf, rightCsf, subshells)
        elseif  op.rank >  0   coeffs = computeCoefficientsNonScalar(op, leftCsf, rightCsf, subshells)
        else    error("stop a")
        end
        
        return (coeffs)
    end


    """
    `SpinAngular.computeCoefficients(op::SpinAngular.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
        ... computes the spin-angular coefficients for the reduced two-particle matrix element
            <leftCsf || op^(k) || rightCsf >   if both CSF refer to the same list of subshells; a Tuple of two lists 
            with one- and two-particle coefficients 
            tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}}  is returned. 
    """
    function  computeCoefficients(op::SpinAngular.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
        if      op.rank == 0   coeffs = computeCoefficientsScalar(op, leftCsf, rightCsf, subshells)
        else    error("stop a")
        end
        
        return (coeffs)
    end


    """
    `SpinAngular.computeCoefficientsNonScalar(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
        ... computes the spin-angular coefficients for the reduced (scalar) one-particle matrix element
            <leftCsf || op^(0) || rightCsf >  if both CSF refer to the same list of subshells; 
            a list of one-particle coefficients coeffs::Array{Coefficient1p,1} is returned.
    """
    function  computeCoefficientsNonScalar(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
        coeffs1p = Coefficient1p[]
        # Cycle through all subshells of the bra- and ket-CSF
        for  (i,shi)  in  enumerate(subshells)
            for  (j,shj)  in  enumerate(subshells)
            
                push!(coeffs1p, Coefficient1p(op.rank, shi, shj, 0.))
            end
        end
        
        return( coeffs1p )
    end


    """
    `SpinAngular.computeCoefficientsScalar(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
        ... computes the spin-angular coefficients for the reduced (scalar) one-particle matrix element
            <leftCsf || op^(0) || rightCsf >  if both CSF refer to the same list of subshells; 
            a list of one-particle coefficients coeffs::Array{Coefficient1p,1} is returned.
    """
    function  computeCoefficientsScalar(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
        coeffs1p = Coefficient1p[]
        # Here, I'm following your notations in Eq. (6) in Gaigalas et al. (CPC, 2001)
        # Cycle through all subshells of the bra- and ket-CSF
        for  (i,shi)  in  enumerate(subshells)
            for  (j,shj)  in  enumerate(subshells)
                if  shi.kappa != shj.kappa      continue    end
                Ni  = leftCsf.occupation[i];   Nj  = rightCsf.occupation[j]
                ji  = Basics.subshell_j(shi);  jj  = Basics.subshell_j(shj)
                Ji  = leftCsf.subshellJ[i];    Jj  = rightCsf.subshellJ[j]
                si  = leftCsf.seniority[i];    sj  = rightCsf.seniority[j]
                Qi  = AngularJ64(1);           Qj  = AngularJ64(1)                          # !! This is yet not correct
                Nri = 0;                       Nrj = 0                                      # !! This is yet not correct
                iTerm     = SpinAngular.SubshellTerm(ji, Qi, si, Ji, Nri)
                jTerm     = SpinAngular.SubshellTerm(jj, Qj, sj, Jj, Nrj)
                braLambda = SpinAngular.Lambda(ji, ji, ji, ji)                              # !! This is yet not correct
                ketLambda = SpinAngular.Lambda(ji, ji, ji, ji)                              # !! This is yet not correct
                wa        = (-1)^SpinAngular.phaseDelta(op, Ni, Nj) * sqrt(Basics.twice(ji)+1) *
                            SpinAngular.recouplingMatrix(op, ji, jj, braLambda, ketLambda)
                #
                if  shi.n == shj.n    wb = SpinAngular.irreducibleTensor(iTerm, Ni, AngularM64(1//2), AngularM64(1//2), iTerm, Ni)     
                else                  wb = 0.  
                end
                if  shi.n != shj.n    wb = wb + SpinAngular.irreducibleTensor(iTerm, Ni, AngularM64(1//2),  jTerm, Nj) *
                                                SpinAngular.irreducibleTensor(iTerm, Ni, AngularM64(-1//2), jTerm, Nj)      end 
                push!(coeffs1p, Coefficient1p(op.rank, shi, shj, wa*wb))
            end
        end
        
        return( coeffs1p )
    end


    """
    `SpinAngular.computeCoefficientsScalar(op::SpinAngular.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
        ... computes the spin-angular coefficients for the reduced (scalar) two-particle matrix element
            <leftCsf || op^(0) || rightCsf >; a Tuple of two lists with one- and two-particle coefficients 
            tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}}  is returned.
    """
    function  computeCoefficientsScalar(op::SpinAngular.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
        coeffs1p   = Coefficient1p[];       coeffs2p    = Coefficient2p[]
        # Here, I'm following your notations in Eq. (13) in Gaigalas et al. (CPC, 2001)
        # Cycle through all subshells of the bra- and ket-CSF
        for  (i,shi)  in  enumerate(subshells)
            for  (j,shj)  in  enumerate(subshells)
                for  (ip,ship)  in  enumerate(subshells)
                    for  (jp,shjp)  in  enumerate(subshells)
                        schemeGamma = SpinAngular.determineSchemeGamma(shi, shj, ship, shjp)
                        schemeEta   = SpinAngular.determineSchemeEta(shi, shj, ship, shjp)
                        
                        Ni  = leftCsf.occupation[i];   Nj  = leftCsf.occupation[j];   Nip = rightCsf.occupation[ip];  Njp  = rightCsf.occupation[jp]
                        ji  = Basics.subshell_j(shi);  jj  = Basics.subshell_j(shj);  jip = Basics.subshell_j(ship);  jjp  = Basics.subshell_j(shjp)
                        Ji  = leftCsf.subshellJ[i];    Jj  = leftCsf.subshellJ[j];    Jip = rightCsf.subshellJ[ip];   Jjp  = rightCsf.subshellJ[jp] 
                        si  = leftCsf.seniority[i];    sj  = leftCsf.seniority[j];    sip = leftCsf.seniority[ip];    sjp  = leftCsf.seniority[jp]
                        Qi  = AngularJ64(1);           Qj  = AngularJ64(1);           Qip = AngularJ64(1);            Qjp  = AngularJ64(1)                          # !! This is yet not 
                        Nri = 0;                       Nrj = 0;                       Nrip= 0;                        Nrjp= 0                                        
                        # !! This is yet not correct
                        iTerm     = SpinAngular.SubshellTerm(ji, Qi, si, Ji, Nri)
                        jTerm     = SpinAngular.SubshellTerm(jj, Qj, sj, Jj, Nrj)
                        ipTerm    = SpinAngular.SubshellTerm(jip, Qip, sip, Jip, Nrip)
                        jpTerm    = SpinAngular.SubshellTerm(jjp, Qjp, sjp, Jjp, Nrjp)
                        #
                        wa = 0.
                        for  k12 = 0:10                                                                 # !! This is yet not correct
                            braLambda = SpinAngular.Lambda(ji, ji, ji, ji)                              # !! This is yet not correct
                            ketLambda = SpinAngular.Lambda(ji, ji, ji, ji)                              # !! This is yet not correct
                            k         = 1                                                               # !! This is yet not correct
                            wa        = wa + (-1)^SpinAngular.phaseDelta(op, Ni, Nj, Nip, Njp) * 
                                        SpinAngular.coefficientThetap(schemeEta, k, shi, shj, ship, shjp) *
                                        SpinAngular.irreducibleTensorT(schemeEta, ji, jj, jip, jjp, braLambda, ketLambda, schemeGamma) *
                                        SpinAngular.recouplingMatrix(op, ji, jj, jip, jjp, braLambda, ketLambda, schemeGamma)
                        end
                        #
                        if  length(coeffs2p) < 7    push!(coeffs2p, Coefficient2p(op.rank, shi, shj, ship, shjp, wa))   end
                        #  This condition is of course arbitrary; from where do we get the one-particle coefficients 
                        #  from kinetic energy, like for a Hamiltonian ? !!
                    end
                end
            end
        end
        
        return( coeffs2p )
    end


    """
    `SpinAngular.coefficientThetap(eta::SchemeEta, k::Int64, shi::Subshell, shj::Subshell, ship::Subshell, shjp::Subshell)`  
        ... computes the coefficient Theta' (ni li ji, nj lj jj; ni' li' ji', nj lj jj', Eta) apart from the effective 
            interaction strength that is set X^k (shi, shj, ship, shjp) = 1 in the procedure; a value::Float64 is returned.
    """
    function  coefficientThetap(eta::SchemeEta, k::Int64, shi::Subshell, shj::Subshell, ship::Subshell, shjp::Subshell)
        # Here likely call: AngularMomentum.
        if      eta == Eta_none
        elseif  eta == Eta_a
        elseif  eta == Eta_W
        else    error("stop a")
        end
        
        return( 1.0 )
    end


    """
    `SpinAngular.determineSchemeEta(shi::Subshell, shj::Subshell, ship::Subshell, shjp::Subshell)`  
        ... determines the coupling scheme Eta for a set of subshells (shi, shj, ship, shjp); 
            a scheme::SchemeEta is returned.
            !! Here, I do not really understand how you read off Eta = Eta_none, ... from the shells and or their occupation;
               this may require to provide further input to this procedure. !!
    """
    function  determineSchemeEta(shi::Subshell, shj::Subshell, ship::Subshell, shjp::Subshell)
        # Determine the correct scheme here
        return( Eta_none )
    end


    """
    `SpinAngular.determineSchemeGamma(shi::Subshell, shj::Subshell, ship::Subshell, shjp::Subshell)`  
        ... determines the coupling scheme Gamma for a set of subshells (shi, shj, ship, shjp); 
            a scheme::SchemeGamma?? is returned.
            !! Here, I do not really understand how you read off Gamma = 1, ..., 42  from the shells and or their occupation;
               this may require to provide further input to this procedure. !!
    """
    function  determineSchemeGamma(shi::Subshell, shj::Subshell, ship::Subshell, shjp::Subshell)
        # Determine the correct scheme here
        return( SchemeGamma01() )
    end


    """
    `SpinAngular.irreducibleTensor(term::SubshellTerm, N::Int64, mq::AngularM64, termp::SubshellTerm, Np::Int64)`  
        ... computes the submatrix elements (j^N alpha Q J || a^(qj)_{m_q} || j'^N' alpha' Q' J')
            as defined by Eq. (7) in Gaigalas et al. (CPC, 2001); a value::Float64 is returned.
            !! We should check of whether we better include the occupation N into SubshellTerm; this would make the communication
               simpler. !!
    """
    function  irreducibleTensor(term::SubshellTerm, N::Int64, mq::AngularM64, termp::SubshellTerm, Np::Int64)
        # Here likely call: AngularMomentum.ClebschGordan, AngularMomentum.Wigner_3j
        
        return( 1.0 )
    end


    """
    `SpinAngular.irreducibleTensor(term::SubshellTerm, N::Int64, mLeft::AngularM64, mRight::AngularM64, termp::SubshellTerm, Np::Int64)`  
        ... computes the submatrix elements (j^N alpha Q J || [a^(qj)_{m_q1} x a^(qj)_{m_q2}]^kj || j'^N' alpha' Q' J')
            as defined by Eq. (8) in Gaigalas et al. (CPC, 2001); a value::Float64 is returned.
            !! We should check of whether we better include the occupation N into SubshellTerm; this would make the communication
               simpler. !!
    """
    function  irreducibleTensor(term::SubshellTerm, N::Int64, mLeft::AngularM64, mRight::AngularM64, termp::SubshellTerm, Np::Int64)
        # Here likely call: AngularMomentum.ClebschGordan, AngularMomentum.Wigner_3j
        
        return( 1.0 )
    end
    

    """
    `SpinAngular.irreducibleTensorT(eta::SchemeEta, ja::AngularJ64, jb::AngularJ64, jap::AngularJ64, jbp::AngularJ64,  
                                    braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda, scheme::SpinAngular.SchemeGamma01)`  
        ... computes matrix elements T(Eta, ja, jb, ja', jb', Lambda^bra, Lambda^ket, Gamma=1) recouplings matrix for a two-particle
            operator; a value::Float64 is returned.
    """
    function  irreducibleTensorT(eta::SchemeEta, ja::AngularJ64, jb::AngularJ64, jap::AngularJ64, jbp::AngularJ64,  
                                 braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda, scheme::SpinAngular.SchemeGamma01)
        # Here likely call functions:  AngularMomentum.
        if      eta == Eta_none
        elseif  eta == Eta_a
        elseif  eta == Eta_W
        else    error("stop a")
        end
        
        return( 1.0 )
    end


    """
    `SpinAngular.phaseDelta(op::SpinAngular.OneParticleOperator, Ni::Int64, Nj::Int64)`  
        ... computes the phase Delta for a one-particle operator from the occupation number of the shells involved;
            a phase::Int64 is returned
    """
    function  phaseDelta(op::SpinAngular.OneParticleOperator, Ni::Int64, Nj::Int64)
        phase = 0                           # !! This is yet not correct
        return( phase )
    end


    """
    `SpinAngular.phaseDelta(op::SpinAngular.TwoParticleOperator, Ni::Int64, Nj::Int64, Nip::Int64, Njp::Int64)`  
        ... computes the phase Delta for a two-particle operator from the occupation number of the shells involved;
            a phase::Int64 is returned
    """
    function  phaseDelta(op::SpinAngular.TwoParticleOperator, Ni::Int64, Nj::Int64, Nip::Int64, Njp::Int64)
        phase = 0                           # !! This is yet not correct
        return( phase )
    end

    
    """
    `SpinAngular.recouplingMatrix(op::SpinAngular.OneParticleOperator, ja::AngularJ64, jb::AngularJ64, 
                                  braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda)`  
        ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) recouplings matrix for a one-particle
            operator; a value::Float64 is returned.
    """
    function  recouplingMatrix(op::SpinAngular.OneParticleOperator, ja::AngularJ64, jb::AngularJ64, 
                               braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda)
        return( 1.0 )
    end


    """
    `SpinAngular.recouplingMatrix(op::SpinAngular.TwoParticleOperator, ja::AngularJ64, jb::AngularJ64, jap::AngularJ64, jbp::AngularJ64, 
                                  braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda, scheme::SpinAngular.SchemeGamma01)`  
        ... computes the recoupling matrix R(ja, jb, ja', jb', Lambda^bra, Lambda^ket, Gamma=1) recouplings matrix for a two-particle
            operator; a value::Float64 is returned.
    """
    function  recouplingMatrix(op::SpinAngular.TwoParticleOperator, ja::AngularJ64, jb::AngularJ64, jap::AngularJ64, jbp::AngularJ64, 
                               braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda, scheme::SpinAngular.SchemeGamma01)
        return( 1.0 )
    end


    """
    `SpinAngular.recouplingMatrix(op::SpinAngular.TwoParticleOperator, ja::AngularJ64, jb::AngularJ64, jap::AngularJ64, jbp::AngularJ64, 
                                  braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda, scheme::SpinAngular.SchemeGamma02)`  
        ... computes the recoupling matrix R(ja, jb, ja', jb', Lambda^bra, Lambda^ket, Gamma=2) recouplings matrix for a two-particle
            operator; a value::Float64 is returned.
    """
    function  recouplingMatrix(op::SpinAngular.TwoParticleOperator, ja::AngularJ64, jb::AngularJ64, jap::AngularJ64, jbp::AngularJ64, 
                               braLambda::SpinAngular.Lambda, ketLambda::SpinAngular.Lambda, scheme::SpinAngular.SchemeGamma02)
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
