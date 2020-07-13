
"""
`module JAC.InteractionStrength`  
    ... a submodel of JAC that contains all methods for evaluating the interaction strength (reduced matrix elements) 
        for various atomic interactions.
"""
module InteractionStrength

    using  GSL, JAC, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals


    """
    `struct  InteractionStrength.XLCoefficient`  ... defines a type for coefficients of the two-electron (Breit) interaction

        + mu        ::Int64      ... Kind of integral, presently fixed to mu = 5 to keep similarity with RATIP during development.
        + nu        ::Int64      ... Rank of the integral.
        + a         ::Orbital    ... Orbitals a, b, c, d.
        + b         ::Orbital
        + c         ::Orbital
        + d         ::Orbital
        + coeff     ::Float64    ... corresponding coefficient.
    """
    struct XLCoefficient 
        mu          ::Int64
        nu          ::Int64
        a           ::Orbital
        b           ::Orbital
        c           ::Orbital
        d           ::Orbital
        coeff       ::Float64 
    end 


    """
    `InteractionStrength.dipole(a::Orbital, b::Orbital, grid::Radial.Grid)`  
        ... computes the  <a|| d ||b>  reduced matrix element of the dipole operator for orbital functions a, b. 
            A value::Float64 is returned. 
    """
    function dipole(a::Orbital, b::Orbital, grid::Radial.Grid)
        wa = AngularMomentum.CL_reduced_me(a.subshell, 1, b.subshell) * RadialIntegrals.rkDiagonal(1, a, b, grid)
        return( wa )
    end


    """
    `InteractionStrength.hamiltonian_nms(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
        ... computes the  <a|| h_nms ||b>  reduced matrix element of the normal-mass-shift Hamiltonian for orbital 
            functions a, b. A value::Float64 is returned.  For details, see Naze et al., CPC 184 (2013) 2187, Eq. (37).
    """
    function hamiltonian_nms(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
        if  a.subshell.kappa != b.subshell.kappa   return( 0. )   end
        wa = RadialIntegrals.isotope_nms(a, b, nm.Z, grid) / (2 * nm.mass)
        println("**  <$(a.subshell) || h^nms || $(b.subshell)>  = $wa" )
        return( wa )
    end


    """
    `InteractionStrength.hfs_t1(a::Orbital, b::Orbital, grid::Radial.Grid)`  
        ... computes the <a|| t^(1) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment 
            of the nucleus for orbital functions a, b. A value::Float64 is returned.  
    """
    function hfs_t1(a::Orbital, b::Orbital, grid::Radial.Grid)
        # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
        minusa = Subshell(1, -a.subshell.kappa)
        wb =   - (a.subshell.kappa + b.subshell.kappa) * AngularMomentum.CL_reduced_me_rb(minusa, 1, b.subshell)
        wc =   RadialIntegrals.rkNonDiagonal(-2, a, b, grid)
        wa =   wb * wc
        #
        ##x println("**  <$(a.subshell) || t1 || $(b.subshell)>  = $wa   = $wb * $wc" )
        return( wa )
    end


    """
    `InteractionStrength.hfs_t2(a::Orbital, b::Orbital, grid::Radial.Grid)`  
        ... computes the <a|| t^(2) ||b> reduced matrix element for the HFS coupling to the electric-quadrupole moment of 
            the nucleus for orbital functions a, b. A value::Float64 is returned.  
     """
    function hfs_t2(a::Orbital, b::Orbital, grid::Radial.Grid)
        # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
        wb = - AngularMomentum.CL_reduced_me_rb(a.subshell, 2, b.subshell)
        ## wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
        wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
        wa =   wb * wc
        #
        ##x println("**  <$(a.subshell) || t2 || $(b.subshell)>  = $wa   = $wb * $wc" )
        return( wa )
    end


    """
    `InteractionStrength.MbaAbsorptionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
        ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp, absorption) || a> 
            for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
            within the given gauge. A value::Float64 is returned.  
    """
    function MbaAbsorptionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
        wa = MbaEmissionCheng(mp, gauge, omega, b, a, grid)
        wa = conj(wa)
        return( wa )
    end


    """
    `InteractionStrength.MbaEmissionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
        ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp,emission) || a> 
            for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
            within the given gauge. A value::Float64 is returned.  
    """
    function MbaEmissionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
        kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
        #
        if       gauge == Basics.Magnetic
            ChengI = AngularMomentum.ChengI(-kapa, kapb, AngularJ64(mp.L))
            wa     = -1.0im / sqrt(mp.L*(mp.L+1)) * ChengI * (kapb + kapa) * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
        #
        elseif   gauge == Basics.Babushkin
            ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L))
            wr     = (mp.L+1) * RadialIntegrals.GrantJL(mp.L, q, a, b, grid::Radial.Grid)
            wr     = wr +  (kapb-kapa-mp.L-1) * RadialIntegrals.GrantIL0(mp.L+1, q, a, b, grid::Radial.Grid)
            wr     = wr +  (kapb-kapa+mp.L+1) * RadialIntegrals.GrantIL0(mp.L+1, q, b, a, grid::Radial.Grid)
            wa     = 1.0im / sqrt(mp.L*(mp.L+1)) * ChengI * wr
        #
        elseif   gauge == Basics.Coulomb
            ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L))
            wr     = (mp.L+1) * (kapb-kapa+mp.L) / (2mp.L+1) * RadialIntegrals.GrantIL0(mp.L-1, q, a, b, grid::Radial.Grid)
            wr     = wr +  ((mp.L+1) * (kapb-kapa+mp.L) - (2mp.L+1) * (kapb-kapa)) / (2mp.L+1) * 
                                                               RadialIntegrals.GrantIL0(mp.L+1, q, a, b, grid::Radial.Grid)
            wr     = wr + (mp.L+1) * (kapb-kapa-mp.L) / (2mp.L+1) * RadialIntegrals.GrantIL0(mp.L-1, q, b, a, grid::Radial.Grid)
            wr     = wr +  ((mp.L+1) * (kapb-kapa-mp.L) - (2mp.L+1) * (kapb-kapa)) / (2mp.L+1) * 
                                                               RadialIntegrals.GrantIL0(mp.L+1, q, b, a, grid::Radial.Grid) 
            wa     = 1.0im / sqrt(mp.L*(mp.L+1)) * ChengI * wr
        else     error("stop a")
        end
        wa = conj(wa)
 
        return( wa )
    end


    """
    `InteractionStrength.MbaEmissionAndrey(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`  
        ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp, emission) || a>  
            for the interaction with the Mp multipole component of the radiation field and the transition frequency omega,
            and within the given gauge. A value::Float64 is returned. At present, only the magnetic matrix elements are 
            implemented. 
    """
    function MbaEmissionAndrey(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
        kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
        ja   = Basics.subshell_j(a.subshell);   jb   = Basics.subshell_j(b.subshell);   
        #
        if       gauge == Basics.Magnetic
            wa = -1.0im * sqrt( (2mp.L+1) * (Basics.subshell_2j(a.subshell)+1) ) / sqrt(4pi) / sqrt(mp.L * (mp.L+1)) * (-1)^mp.L * (kapb + kapa)
            wa = wa * AngularMomentum.ClebschGordan(ja, AngularM64(1//2), AngularJ64(mp.L), AngularM64(0), jb, AngularM64(1//2)) 
            wa = wa * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
            wa = conj(wa)
        #
        elseif   gauge == Basics.Babushkin
            wa = 0.
        #
        elseif   gauge == Basics.Coulomb
            wa = 0.
        else     error("stop a")
        end
         
        return( wa )
    end


    """
    `InteractionStrength.multipoleTransition(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
        ... to compute the (single-electron reduced matrix element) multipole-transition interaction strength 
            <b || T^(Mp, absorption) || a> due to Johnson (2007) for the interaction with the Mp multipole component of the
            radiation field and the transition frequency omega, and within the given gauge. A value::Float64 is returned.  
    """
    function multipoleTransition(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
        function besselPrime_jl(L::Int64, x::Float64)    return( GSL.sf_bessel_jl(L-1, x) - (L+1)/x * GSL.sf_bessel_jl(L, x) )       end

        kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
        mtp  = min(size(a.P, 1), size(b.P, 1))
        !(grid.meshType == Radial.MeshGL())  &&  error("Only for Radial.MeshGL() implemented so far.")
        #
        if       gauge == Basics.Magnetic
            ChengI = AngularMomentum.ChengI(-kapa, kapb, AngularJ64(mp.L));   if  abs(ChengI) < 1.0e-10  return( 0. )   end
            wa = Complex(0.)
            for  i = 2:mtp
                wa = wa + (kapa+kapb) / (mp.L+1) * GSL.sf_bessel_jl(mp.L, q * grid.r[i]) * (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * grid.wr[i]  
            end
            wa = ChengI * wa
            #
        elseif   gauge == Basics.Velocity
            ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L));    if  abs(ChengI) < 1.0e-10  return( 0. )   end
            wa = Complex(0.)
            for  i = 2:mtp
                wa = wa - (kapa-kapb) / (mp.L+1) * 
                          ( besselPrime_jl(mp.L, q * grid.r[i]) + GSL.sf_bessel_jl(mp.L, q * grid.r[i]) / (q * grid.r[i]) ) * 
                          (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * grid.wr[i]
                wa = wa + mp.L * GSL.sf_bessel_jl(mp.L, q * grid.r[i]) / (q * grid.r[i]) * (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) * grid.wr[i]  
            end
            wa = ChengI * wa
            #
        elseif   gauge == Basics.Length
            ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L));    if  abs(ChengI) < 1.0e-10  return( 0. )   end
            wa = Complex(0.)
            for  i = 2:mtp
                wa = wa + GSL.sf_bessel_jl(mp.L, q * grid.r[i]) * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * grid.wr[i] 
                        + GSL.sf_bessel_jl(mp.L+1, q * grid.r[i]) * 
                          ( (kapa-kapb) / (mp.L+1) * (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) +
                            (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) ) * grid.wr[i]
            end
            wa = ChengI * wa
            #
        else     error("stop a")
        end
 
        return( wa )
    end


    """
    `InteractionStrength.schiffMoment(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
        ... computes the  <a|| h^(Schiff-moment) ||b>  reduced matrix element of the Schiff-moment Hamiltonian for orbital 
            functions a, b and for the nuclear density as given by the nuclear model. A value::Float64 is returned.  
    """
    function schiffMoment(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
        printstyled("\nWarning -- InteractionStrength.schiffMoment():: Not yet implemented.", color=:cyan)
        wb = 1.0 + 2.0im
        return( wb )
    end


    """
    `InteractionStrength.weakCharge(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
        ... computes the  <a|| h^(weak-charge) ||b>  reduced matrix element of the weak-charge Hamiltonian for orbital functions 
            a, b and for the nuclear density as given by the nuclear model. A value::Float64 is returned.  
    """
    function weakCharge(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
        printstyled("\nWarning -- InteractionStrength.weakCharge():: Not yet implemented.", color=:cyan)
        wb = 1.0 + 2.0im
        return( wb )
    end


    """
    `InteractionStrength.XL_Breit_WO(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
        ... computes the the effective Breit interaction strengths X^L_Breit (abcd) for given rank L and orbital functions 
            a, b, c and d  at the given grid but without optimization. A value::Float64 is returned. At present, only the zero-frequency 
            Breit interaction is taken into account.
    """
    function XL_Breit_WO(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
        ja2 = Basics.subshell_2j(a.subshell)
        jb2 = Basics.subshell_2j(b.subshell)
        jc2 = Basics.subshell_2j(c.subshell)
        jd2 = Basics.subshell_2j(d.subshell)
        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||  L == 0  
            return( 0. )
        end
        #
        xcList   = XL_Breit0_coefficients(L,a,b,c,d)
        ##x println("\nBreit strength: L=$L, nx = $(length(xcList)) ")
        XL_Breit = XL_Breit0_densities(xcList, grid)
        #
        return( XL_Breit )
    end
    
    
    """
    `InteractionStrength.XL_Breit_reset_storage(keep::Bool)`  
        ... resets the global storage of XL_Breit interaction strength; nothing is returned.
    """
    function XL_Breit_reset_storage(keep::Bool)
        if  keep
            println("  reset GBL_Storage_XL_Breit storage ...")
            global GBL_Storage_XL_Breit = Dict{String, Float64}()
        else
        end
        return( nothing )      
    end


    """
    `InteractionStrength.XL_Breit(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)`  
        ... computes the the effective Breit interaction strengths X^L_Breit (abcd) for given rank L and orbital functions 
            a, b, c and d  at the given grid. For keep=true, the procedure looks up the (global) directory GBL_Storage_XL_Coulomb
            and returns the corresponding value without re-calculation of the interaction strength; it also 'stores' the calculated
            value if not yet included. For keep=false, the interaction strength is always computed on-fly. A value::Float64 is returned. 
            At present, only the zero-frequency Breit interaction is taken into account.
    """
    function XL_Breit(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)
        global GBL_Storage_XL_Breit
        ja2 = Basics.subshell_2j(a.subshell)
        jb2 = Basics.subshell_2j(b.subshell)
        jc2 = Basics.subshell_2j(c.subshell)
        jd2 = Basics.subshell_2j(d.subshell)
        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||  L == 0  
            return( 0. )
        end
        
        # Now distiguish due to the optional argument keep
        if  keep
            sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
            if haskey(GBL_Storage_XL_Breit, sa )
                XL_Breit = GBL_Storage_XL_Breit[sa]
            else
                xcList   = XL_Breit0_coefficients(L,a,b,c,d)
                XL_Breit = XL_Breit0_densities(xcList, grid)
                global GBL_Storage_XL_Breit = Base.merge(GBL_Storage_XL_Breit, Dict( sa => XL_Breit))
            end
        else
            xcList   = XL_Breit0_coefficients(L,a,b,c,d)
            XL_Breit = XL_Breit0_densities(xcList, grid)
        end
        #

        return( XL_Breit )
    end


    """
    `InteractionStrength.XL_BreitDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
        ... computes the the effective Breit interaction strengths X^L_Breit (abcd) for given rank L and orbital functions 
            a, b, c and d  at the given grid. A value::Float64 is returned. At present, only the zero-frequency Breit 
            interaction is taken into account.
    """
    function XL_BreitDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
        error("stop a")
    end


    """
    `InteractionStrength.XL_Breit0_coefficients(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital)`  
        ... evaluates the combinations and pre-coefficients for the zero-frequency Breit interaction  
            X^L_Breit (omega=0.; abcd) for given rank L and orbital functions a, b, c and d. A list of coefficients 
            xcList::Array{XLCoefficient,1} is returned.
    """
    function XL_Breit0_coefficients(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital)
        xcList = XLCoefficient[]
        
        la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
        lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
        lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
        ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

        xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
        if   rem(L,2) == 1    xc = - xc                end 
        if   abs(xc)  <  1.0e-10    return( xcList )   end

        # Consider the individual contributions from sum_nu and sum_mu. First, take T^(nu,L)_mu = R^(nu,L)_mu
        nu = L - 1

        if  rem(la+lc+nu,2) == 1   &&   rem(lb+ld+nu,2) == 1   &&   L != 0
            wa = (L+1) / ( L*(L+L-1)*(L+L+1) ) 
            # mu = 1
            xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa+L) * (d.subshell.kappa-b.subshell.kappa+L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, a, b, c, d, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, b, a, d, c, xcc) )   end
            # mu = 2
            xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa-L) * (d.subshell.kappa-b.subshell.kappa-L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, c, d, a, b, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, d, c, b, a, xcc) )   end
            # mu = 3
            xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa+L) * (d.subshell.kappa-b.subshell.kappa-L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, a, d, c, b, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, d, a, b, c, xcc) )   end
            # mu = 4
            xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa-L) * (d.subshell.kappa-b.subshell.kappa+L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, c, b, a, d, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, b, c, d, a, xcc) )   end
        end
        #
        nu = L

        if  rem(la+lc+nu,2) == 1   &&   rem(lb+ld+nu,2) == 1   &&   L != 0
            wa = - (a.subshell.kappa + c.subshell.kappa) * (b.subshell.kappa + d.subshell.kappa) / (L*(L+1))
            # mu = 1
            xcc = xc * wa
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, a, b, c, d, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, b, a, d, c, xcc) )   end
            # mu = 2
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, c, d, a, b, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, d, c, b, a, xcc) )   end
            # mu = 3
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, a, d, c, b, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, d, a, b, c, xcc) )   end
            # mu = 4
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, c, b, a, d, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, b, c, d, a, xcc) )   end
        end
        #
        nu = L + 1

        if  rem(la+lc+nu,2) == 1   &&   rem(lb+ld+nu,2) == 1   &&   L != 0
            wa =  L / ( (L+1)*(L+L+1)*(L+L+3) )
            # mu = 1
            xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa - L - 1) * (d.subshell.kappa - b.subshell.kappa - L - 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, a, b, c, d, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, b, a, d, c, xcc) )   end
            # mu = 2
            xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa + L + 1) * (d.subshell.kappa - b.subshell.kappa + L + 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, c, d, a, b, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, d, c, b, a, xcc) )   end
            # mu = 3
            xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa - L - 1) * (d.subshell.kappa - b.subshell.kappa + L + 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, a, d, c, b, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, d, a, b, c, xcc) )   end
            # mu = 4
            xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa + L + 1) * (d.subshell.kappa - b.subshell.kappa - L - 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, c, b, a, d, xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, nu, b, c, d, a, xcc) )   end
        end

        # Add contributions of the S^k_mu integrals
        if  rem(la+lc+L-1,2) == 1   &&   rem(lb+ld+L+1,2) == 1
            # mu = 1
            wb =  1 / ( (L+L+1)*(L+L+1) )
            xcc = xc * wb * (c.subshell.kappa - a.subshell.kappa + L) * (d.subshell.kappa - b.subshell.kappa - L - 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, b, a, d, c,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, b, a, d, c,  -(L+L+1) / 2 * xcc) )   end
            # mu = 2
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L) * (c.subshell.kappa - a.subshell.kappa - L - 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, a, b, c, d,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, a, b, c, d,  -(L+L+1) / 2 * xcc) )   end
            # mu = 3
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L + 1) * (c.subshell.kappa - a.subshell.kappa - L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, d, c, b, a,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, d, c, b, a,  -(L+L+1) / 2 * xcc) )   end
            # mu = 4
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L) * (c.subshell.kappa - a.subshell.kappa + L + 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, c, d, a, b,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, c, d, a, b,  -(L+L+1) / 2 * xcc) )   end
            # mu = 5
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L + 1) * (c.subshell.kappa - a.subshell.kappa + L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, d, a, b, c,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, d, a, b, c,  -(L+L+1) / 2 * xcc) )   end
            # mu = 6
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L) * (c.subshell.kappa - a.subshell.kappa - L - 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, a, d, c, b,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, a, d, c, b,  -(L+L+1) / 2 * xcc) )   end
            # mu = 7
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L - 1) * (c.subshell.kappa - a.subshell.kappa - L)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, b, c, d, a,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, b, c, d, a,  -(L+L+1) / 2 * xcc) )   end
            # mu = 8
            xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L) * (c.subshell.kappa - a.subshell.kappa + L + 1)
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L+1, c, b, a, d,   (L+L+1) / 2 * xcc) )   end
            if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient(5, L-1, c, b, a, d,  -(L+L+1) / 2 * xcc) )   end
        end

        return( xcList )
    end


    """
    `InteractionStrength.XL_Breit0_densities(xcList::Array{XLCoefficient,1}, grid::Radial.Grid)`  
        ... computes the the effective Breit interaction strengths X^L,0_Breit (abcd) for given rank L and a list of 
            orbital functions a, b, c, d and angular coefficients at the given grid. A value::Float64 is returned. 
            At present, only the zero-frequency Breit interaction is taken into account.
    """
    function XL_Breit0_densities(xcList::Array{XLCoefficient,1}, grid::Radial.Grid)
        
        if  grid.meshType == Radial.MeshGL()
            wa = 0.
            for  xc  in  xcList  ## [end:end]
                # Use the minimal extent of any involved orbitals; this need to be improved
                mtp_ac = min(size(xc.a.P, 1), size(xc.c.P, 1));   mtp_bd = min(size(xc.b.P, 1), size(xc.d.P, 1))
                for  r = 2:mtp_ac
                    for  s = 2:mtp_bd
                        wc = xc.coeff * grid.wr[r] * grid.wr[s] 
                        if      s > r   continue
                        elseif  s == r
                            wa = wa + wc * (xc.a.P[r] * xc.c.Q[r]) * (grid.r[s]^xc.nu) / (grid.r[r]^(xc.nu+1)) * (xc.b.P[s] * xc.d.Q[s]) / 2.0
                        else 
                            wa = wa + wc * (xc.a.P[r] * xc.c.Q[r]) * (grid.r[s]^xc.nu) / (grid.r[r]^(xc.nu+1)) * (xc.b.P[s] * xc.d.Q[s])
                        end
                    end
                end
            end
            return( wa )
        else
            error("stop a")
        end
    end


    """
    `InteractionStrength.matrixL_Coulomb(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, primitives::Bsplines.Primitives)`  
        ... computes the partly-contracted (effective) Coulomb interaction matrices M^L_Coulomb (abcd) for given rank L and orbital functions 
            a, b, c and d at the given grid. The matrix M^L is defined for the primitives and contracted over the two orbitals
            b, d (for a=c) or  b, c (for a=d).  An error message is issued if a != c && a != d. A matrix::Array{Float64,2} is returned.
    """
    function matrixL_Coulomb(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, primitives::Bsplines.Primitives)
        grid = primitives.grid;   nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        wm = zeros( nsL+nsS, nsL+nsS )
        # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
        la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
        lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
        lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
        ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
            rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
            return( wm )
        end
        xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
        if   rem(L,2) == 1    xc = - xc    end 

        if  a.subshell == c.subshell
            # Direct interaction; contract the full interaction array over the orbitals b and d
            for  i = 1:nsL
                for  k = 1:nsL 
                    Ba = primitives.bsplinesL[i].bs;    Bc = primitives.bsplinesL[k].bs
                    wm[i,k] = RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.P, Bc, d.P, grid) + 
                              RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.Q, Bc, d.Q, grid)
                end
            end
            for  i = 1:nsS
                for  k = 1:nsS
                    Ba = primitives.bsplinesS[i].bs;    Bc = primitives.bsplinesS[k].bs
                    wm[nsL+i,nsL+k] = RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.P, Bc, d.P, grid) + 
                                      RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.Q, Bc, d.Q, grid)
                end
            end
        elseif true 
            println("Skip exchange integrals")
            return( wm )
        elseif  a.subshell == d.subshell
            # Exchange interaction; contract the full interaction array over the orbitals b and c
            for  i = 1:nsL
                for  k = 1:nsL 
                    Ba = primitives.bsplinesL[i].bs;    Bd = primitives.bsplinesL[k].bs
                    wm[i,k] = RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.P, c.P, Bd, grid)
                end
                for  k = 1:nsS 
                    Ba = primitives.bsplinesL[i].bs;    Bd = primitives.bsplinesS[k].bs
                    wm[i,nsL+k] = RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.P, c.Q, Bd, grid)
                end
            end
            for  i = 1:nsS
                for  k = 1:nsL 
                    Ba = primitives.bsplinesS[i].bs;    Bd = primitives.bsplinesL[k].bs
                    wm[nsL+i,k] = RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.Q, c.P, Bd, grid)
                end
                for  k = 1:nsS 
                    Ba = primitives.bsplinesS[i].bs;    Bd = primitives.bsplinesS[k].bs
                    wm[nsL+i,nsL+k] = RadialIntegrals.SlaterRkComponent_2dim(L, Ba, b.Q, c.Q, Bd, grid)
                end
            end
        else    error("stop d")
        end

        return( wm )
    end


    """
    `InteractionStrength.XS_Coulomb(large::Bool, largep::Bool, coeffs::Array{JAC.AngularCoefficientsRatip2013.AngularVcoeff,1}, 
                                    orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid; exchange::Bool=false)`  
        ... computes the effective XS_Coulomb interaction function for given components (large, largep) and a set of (Coulomb-like) angular coefficients.
            The procedures assumes that all coefficients V^K (abcd) fulfill  a==c || a == d, and a error is issued if this is not the case.
            For the remaining two orbitals b, d (or c), the correct components are added due to boolean selectors.
            The XS_Coulomb interaction function is defined by
            
                S(r) = sum_k  Int ds  U(r,s)  V^k(abcd) b.P(large) cd.P(largep)
                
            where b.P(large) refers to b.P of orbital b if large=true and b.Q if large= false, and analogue for cd.P. A Sfunc::Float64[] is returned.
    """
    function XS_Coulomb(large::Bool, largep::Bool, coeffs::Array{JAC.AngularCoefficientsRatip2013.AngularVcoeff,1}, 
                        orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid; exchange::Bool=false)
        function ul(k::Int64, r::Float64, s::Float64) :: Float64
            if     r <= s    return( r^k/s^(k+1) )
            elseif r > s     return( s^k/r^(k+1) )
            end
        end
        #
        Sfunc = zeros( grid.nr );   na = 1
        for  r = 2:grid.nr
            # Determine integrant over s
            ws = zeros( grid.nr )
            for  coeff in  coeffs
                a = orbitals[coeff.a];   b = orbitals[coeff.b];   c = orbitals[coeff.c];   d = orbitals[coeff.d];   L = coeff.nu
                # Determine the weight xc of this coefficient
                la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
                lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
                lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
                ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

                if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
                    rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1      error("stop a")
                end
                xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
                if   rem(L,2) == 1    xc = - xc    end 
                # Divide total coefficient by occupation w_a of the considered shell (DF equations vs. energy functional)
                xc = xc / (ja2 + 1)
                #
                ##x if na < 6   na = na + 1;    @show L, coeff.a, coeff.b, coeff.c, coeff.d, coeff.V * xc  end
                #
                # Now add the contributions for this coefficient
                if      coeff.a == coeff.c  &&  large  &&  largep
                    mtp_bd = min(size(b.P, 1), size(d.P, 1));    for  s = 2:mtp_bd   ws[s] = ws[s] + coeff.V * xc * b.P[s] * d.P[s]* ul(L, grid.r[r], grid.r[s])    end 
                elseif  coeff.a == coeff.c  &&  !(large)  &&  !(largep)
                    mtp_bd = min(size(b.Q, 1), size(d.Q, 1));    for  s = 2:mtp_bd   ws[s] = ws[s] + coeff.V * xc * b.Q[s] * d.Q[s]* ul(L, grid.r[r], grid.r[s])    end 
                elseif  coeff.a == coeff.d  &&  exchange  &&  large  &&  largep         
                    mtp_bc = min(size(b.P, 1), size(c.P, 1));    for  s = 2:mtp_bc   ws[s] = ws[s] + coeff.V * xc * b.P[s] * c.P[s]* ul(L, grid.r[r], grid.r[s])    end 
                elseif  coeff.a == coeff.d  &&  exchange  &&  !(large)  &&  largep 
                    mtp_bc = min(size(b.Q, 1), size(c.P, 1));    for  s = 2:mtp_bc   ws[s] = ws[s] + coeff.V * xc * b.Q[s] * c.P[s]* ul(L, grid.r[r], grid.r[s])    end 
                elseif  coeff.a == coeff.d  &&  exchange  &&  large     &&  !(largep)
                     mtp_bc = min(size(b.P, 1), size(c.Q, 1));    for  s = 2:mtp_bc   ws[s] = ws[s] + coeff.V * xc * b.P[s] * c.Q[s]* ul(L, grid.r[r], grid.r[s])    end 
                elseif  coeff.a == coeff.d  &&  exchange  &&  !(large)  &&  !(largep)
                    mtp_bc = min(size(b.Q, 1), size(c.Q, 1));    for  s = 2:mtp_bc   ws[s] = ws[s] + coeff.V * xc * b.Q[s] * c.Q[s]* ul(L, grid.r[r], grid.r[s])    end 
                elseif  coeff.a == coeff.d  &&  !(exchange)
                elseif  coeff.a == coeff.c  &&  (large     &&  !(largep)   ||   !(large)  &&  largep)
                else    @show L, coeff.a, coeff.b, coeff.c, coeff.d, coeff.V * xc
                end
            end
            #
            for  s = 2:grid.nr   Sfunc[r] = Sfunc[r]  +  ws[s] * grid.wr[s]  end
        end
        
        return( Sfunc )
    end


    """
    `InteractionStrength.XL_Coulomb_WO(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
        ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
            a, b, c and d at the given grid but without optimization. A value::Float64 is returned.
    """
    function XL_Coulomb_WO(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
        # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
        la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
        lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
        lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
        ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
            rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
            return( 0. )
        end
        xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
        if   rem(L,2) == 1    xc = - xc    end 
        
        XL_Coulomb = xc * RadialIntegrals.SlaterRk_2dim_WO(L, a, b, c, d, grid)
        ##x XL_Coulomb = xc * RadialIntegrals.SlaterRk_new(L, a, b, c, d, grid)
        return( XL_Coulomb )
    end
    
    
    """
    `InteractionStrength.XL_Coulomb_reset_storage(keep::Bool; printout::Bool=false)`  
        ... resets the global storage of XL_Coulomb interaction strength; nothing is returned.
    """
    function XL_Coulomb_reset_storage(keep::Bool; printout::Bool=false)
        if  keep
            if printout     println(">> Reset GBL_Storage_XL_Coulomb storage.")     end
            global GBL_Storage_XL_Coulomb = Dict{String, Float64}()
        else
        end
        return( nothing )      
    end


    """
    `InteractionStrength.XL_Coulomb(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)`  
        ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
            a, b, c and d at the given grid. For keep=true, the procedure looks up the (global) directory GBL_Storage_XL_Coulomb
            and returns the corresponding value without re-calculation of the interaction strength; it also 'stores' the calculated
            value if not yet included. For keep=false, the interaction strength is always computed on-fly. A value::Float64 is 
            returned.
    """
    function XL_Coulomb(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)
        global GBL_Storage_XL_Coulomb
        # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
        la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
        lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
        lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
        ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
            rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
            return( 0. )
        end
        
        # Now distiguish due to the optional argument keep
        if  keep
            sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
            if haskey(GBL_Storage_XL_Coulomb, sa )
                XL_Coulomb = GBL_Storage_XL_Coulomb[sa]
            else
                XL_Coulomb = InteractionStrength.XL_Coulomb(L::Int64, a, b, c, d, grid)
                ## global GBL_Storage_XL_Coulomb = Base.merge(GBL_Storage_XL_Coulomb, Dict( sa => XL_Coulomb))
                global GBL_Storage_XL_Coulomb[sa] = XL_Coulomb
            end
        else
            xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
            if   rem(L,2) == 1    xc = - xc    end 
            XL_Coulomb = xc * RadialIntegrals.SlaterRk_2dim(L, a, b, c, d, grid)
        end
        
        return( XL_Coulomb )
    end


    """
    `InteractionStrength.XL_CoulombDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
        ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
            a, b, c and d at the given grid. A value::Float64 is returned.
    """
    function XL_CoulombDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
        # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
        la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
        lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
        lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
        ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
            rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
            return( 0. )
        end
        xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
        if   rem(L,2) == 1    xc = - xc    end 
        
        XL_Coulomb = xc * RadialIntegrals.SlaterRk_2dim_Damped(tau::Float64, L, a, b, c, d, grid)
        ##x XL_Coulomb = xc * RadialIntegrals.SlaterRk_new(L, a, b, c, d, grid)
        return( XL_Coulomb )
    end


    """
    `InteractionStrength.XL_Coulomb_DH(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid, lambda::Float64)`  
        ... computes the the effective Coulomb-Debye-Hückel interaction strengths X^L_Coulomb_DH (abcd) for given rank L and 
            orbital functions a, b, c and d at the given grid and for the given screening parameter lambda. A value::Float64 is 
            returned.
    """
    function XL_Coulomb_DH(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid, lambda::Float64)
        # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
        la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
        lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
        lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
        ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
            rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
            return( 0. )
        end
        xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
        if   rem(L,2) == 1    xc = - xc    end 

        XL_Coulomb_DH = xc * RadialIntegrals.SlaterRk_DebyeHueckel_2dim(L, a, b, c, d, grid, lambda)
        return( XL_Coulomb_DH )
    end


    """
    `InteractionStrength.XL_plasma_ionSphere(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, lambda::Float64)`  
        ... computes the effective interaction strengths X^L_ion-sphere (abcd) for given rank L and orbital functions 
            a, b, c and d and for the plasma parameter lambda. A value::Float64 is returned.  **Not yet implemented !**
    """
    function XL_plasma_ionSphere(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, lambda::Float64)
        error("Not yet implemented")
    end


    """
    `InteractionStrength.X1_smsA(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
        ... computes the the effective interaction strengths X^1_sms,A (abcd) for fixed rank 1 and orbital functions 
            a, b, c and d at the given grid. A value::Float64 is returned.
    """
    function X1_smsA(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
        wa = AngularMomentum.CL_reduced_me(a.subshell, 1, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, 1, d.subshell) *
             RadialIntegrals.Vinti(a, c, grid) * RadialIntegrals.Vinti(b, d, grid) / (2 * nm.mass)
        ## println("**  <$(a.subshell) || Vinti || $(c.subshell)>  = $(RadialIntegrals.Vinti(a, c, grid)) " )
        ## println("**  <$(b.subshell) || Vinti || $(d.subshell)>  = $(RadialIntegrals.Vinti(b, d, grid)) " )
        return( wa )
    end


    """
    `InteractionStrength.X1_smsB(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
        ... computes the the effective interaction strengths X^1_sms,B (abcd) for fixed rank 1 and orbital functions 
            a, b, c and d at the given grid. A value::Float64 is returned.
    """
    function X1_smsB(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
        wa = - AngularMomentum.CL_reduced_me(b.subshell, 1, d.subshell) * RadialIntegrals.Vinti(b, d, grid) *
               RadialIntegrals.isotope_smsB(a, c, nm.Z, grid) / (2 * nm.mass) 
        return( wa )
    end


    """
    `InteractionStrength.X1_smsC(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)` 
        ... computes the the effective interaction strengths X^1_sms,C (abcd) for fixed rank 1 and orbital functions 
            a, b, c and d at the given grid. A value::Float64 is returned.
    """
    function X1_smsC(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
        wa = - AngularMomentum.CL_reduced_me(b.subshell, 1, d.subshell) * AngularMomentum.CL_reduced_me(a.subshell, 1, c.subshell) * 
               RadialIntegrals.Vinti(b, d, grid) * RadialIntegrals.isotope_smsC(a, c, nm.Z, grid) / (2 * nm.mass) 
        return( wa )
    end


    """
    `InteractionStrength.zeeman_Delta_n1(a::Orbital, b::Orbital, grid::Radial.Grid)`  
        ... computes the <a|| Delta n^(1) ||b> reduced matrix element for the Zeeman-Schwinger contribution to the coupling 
            to an external magnetic field for orbital functions a, b. A value::Float64 is returned.
    """
    function zeeman_Delta_n1(a::Orbital, b::Orbital, grid::Radial.Grid)
        # Use Andersson, Jönson (2008), CPC ... test for the proper definition of the C^L tensors.
        minusa = Subshell(1, -a.subshell.kappa)
        wa = (Defaults.getDefaults("electron g-factor") - 2) / 2. * (a.subshell.kappa + b.subshell.kappa + 1) * 
             AngularMomentum.CL_reduced_me(minusa, 1, b.subshell) * RadialIntegrals.rkDiagonal(0, a, b, grid)
        return( wa )
    end


    """
    `InteractionStrength.zeeman_n1(a::Orbital, b::Orbital, grid::Radial.Grid)`  
        ... computes the <a|| n^(1) ||b> reduced matrix element for the Zeeman coupling to an external magnetic field for 
            orbital functions a, b. A value::Float64 is returned. 
    """
    function zeeman_n1(a::Orbital, b::Orbital, grid::Radial.Grid)
        # Use Andersson, Jönson (2008), CPC ... test for the proper definition of the C^L tensors.
        minusa = Subshell(1, -a.subshell.kappa)
        wa = - AngularMomentum.CL_reduced_me(minusa, 1, b.subshell) / (2 * Defaults.getDefaults("alpha")) *
               RadialIntegrals.rkNonDiagonal(1, a, b, grid)
        return( wa )
    end

end # module

