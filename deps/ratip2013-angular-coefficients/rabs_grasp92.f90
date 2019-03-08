module rabs_grasp92
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module contains a number of `original' GRASP92 procedures which are needed at some intermediate level of the 
! development of RABS. They are (nearly) the same like in the GRASP92 package but have been adapted to the present standard 
! of Fortran 90. In particular, no common arrays or other rather obsolete features occur in the present context. 
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_dirac_orbital
   use rabs_functions_math
   use rabs_nucleus
   implicit none
   !
   private :: Bessel_grasp92
                 ! Calculates one of the functions phi_k(x) or psi_k(x) of the transverse Breit interactions which are 
                 ! closely related to the spherical Bessel functions.
   public  :: debye_integral_explicit
                 ! Calculates the relativistic Debye integral (as a generalization of the Slater integral) for four given
		 ! radial orbital functions. 
   private :: dpbdt_grasp92
                 ! Computes H times the derivative of the large and small components of a given radial wave function.
   public  :: hydrogen_point_grasp92
                 ! Return a hydrogenic function for a point nucleus with given nuclear charge.
   public  :: hydrogenic_bound_orbital
                 ! Computes the  Dirac-Coulomb  bound-state orbital radial wavefunction following Akhiezer and Berestetskii 
                 ! but modified to ensure positive slope at the origin.
   public  :: I_ab_grasp92
                 ! Returns the one-electron integral I (ab) for two given orbital functions a and b.
   public  :: interpolate_rwf_grasp92
                 ! Interpolates and renormalizes the large and small component of a Dirac orbital on a given grid.
   public  :: is_bound_orbital
                 ! Returns .true. if a bound orbital is given as argument and .false. otherwise.
   public  :: load_rwf_file_grasp92
                 ! Reads in a set of GRASP92 orbital functions into an appropriate structure.
   public  :: ncharg_grasp92
                 ! Evaluates the nuclear charge density, and stores it in the array  zdist_grasp92(:)           
   public  :: quad_grasp92
                 ! Performs the quadrature for a given function along the radial grid rp_grasp92 by using a five-point 
                 ! Closed Newton-Cotes formula. 
   public  :: radgrd_grasp92
                 ! Sets up the radial grid  r_grasp92(:)  and the associated arrays  rp_grasp92(:)  and  rpor_grasp92(:).
   public  :: readwrite_grasp92_grid
                 ! Reads or writes the GRASP-92 grid data structure from or to a file.
   public  :: readwrite_grasp92_orbital
                 ! Reads or writes a type(grasp92_orbital) data structure from or to a file.
   private :: Rk_bar_acbd_integral_grasp92
                 ! Evaluates the transverse interaction integral R-bar^k[ac,bd].
   public  :: rk_integral_grasp92
                 ! Calculates the r^k matrix element for two radial orbital function from the same set of (orthonormal) 
                 ! orbitals. 
   public  :: rk_integral_grasp92_ab
                 ! Calculates the r^k matrix element for radial orbital functions which may belong to different orbital sets.
   public  :: rk_integral_grasp92_cd
                 ! Calculates the r^k matrix element for radial orbital functions which are given explicitly.
   private :: schmidt_orthogonalize_grasp92
                 ! Orthogonalizes the radial wavefunctions in a given wave function array.
   public  :: selfenergy_ratio_grasp92
                 ! Returns a self-energy estimate for the a given orbital function.
   public  :: setqic_grasp92
                 ! This  subroutine sets up the coefficients for SUBROUTINEs DPBDT, QUAD, RINTI, START, YZK, ZKF. 
   private :: Sk_acbd_integral_grasp92
                 ! Evaluates the transverse interaction integral S^k[ac,bd].
   public  :: slater_integral_grasp92
                 ! Calculates the relativistic Slater integral for four given radial orbital functions. 
   public  :: slater_integral_explicit
                 ! Calculates the relativistic Slater integral for four given radial orbital functions. 
   public  :: V_ab_grasp92
                 ! Returns the one-electron integral V (ab) for two given orbital functions a and b as appropriate for 
                 ! specific mass-shift calculations.
   public  :: vpintf_grasp92
                 ! Computes nuclear vacuum polarization integrals.            
   public  :: W_integral_grasp92
                 ! Calculates one of the prediefined integrals for four given radial orbital functions. 
   public  :: yz_k_grasp92
                 ! Calculates the Hartree Y- and Z-functions. 
   public  :: zf_k_grasp92
                 ! Calculates the Hartree Z-functionals.
   !
   real(kind=dp), dimension(:), allocatable, public ::    r_grasp92, rp_grasp92, rpor_grasp92, zdist_grasp92, &
                                                          vacpol_2_grasp92, vacpol_4_grasp92
   !
   real(kind=dp), dimension(:,:), allocatable, public ::  rtk_grasp92, rprtk_grasp92, rpbrtk_grasp92, rtki_grasp92
   !
   real(kind=dp), public :: accy_grasp92 = 1.0e-8_dp
   !
   real(kind=dp), public :: c1_grasp92, c2_grasp92, c3_grasp92, c4_grasp92
   real(kind=dp), dimension(2:5,2:4), public :: cnc5c_grasp92
   real(kind=dp), dimension(1:6,2:6), public :: cnc6c_grasp92
   real(kind=dp), dimension(13,13),   public :: a13_grasp92(13,13)
   real(kind=dp), dimension(6),       public :: c_grasp92(6)
   !
   ! Define an appropriate data structure for keeping the radial orbital wave functions 
   type, public :: grasp92_orbital
      integer :: number_of_rwf
      type(orbital_function), dimension(:), pointer :: rwf
   end type grasp92_orbital
   !
   type, public :: orbital_function
      type(nkappa)  :: orbital
      integer       :: mtp
      real(kind=dp) :: energy, gamma, pz, phase
      real(kind=dp), dimension(:), pointer :: P, Q   
   end type orbital_function
   !
   ! Define several counter
   integer, public :: grasp92_Rk_slater                 = 0, &
                      grasp92_Rkbar_integral            = 0, &
                      grasp92_Sk_integral               = 0
   !
   ! Define some logical flags for debugging
   logical, public :: debug_I_ab_grasp92         = .true.,  &
                      debug_radgrd_grasp92       = .false., &
                      debug_load_rwf_grasp92     = .false., &
                      debug_schmidt_grasp92      = .false., &
                      debug_interpolate_grasp92  = .false., &
                      debug_vpintf_grasp92       = .true.,  &
                      debug_V_ab_grasp92         = .true.,  &
                      debug_Rk_bar_acbd_grasp92  = .true.,  &
                      debug_Sk_acbd_grasp92      = .true.
   !
contains
   !
   subroutine Bessel_grasp92(k,omega,bessel_j,bessel_n)  
   !--------------------------------------------------------------------------------------------------------------------
   ! Evaluates the functions 
   !
   !                 (2k+1)!!                                           
   !   bessel_j  =   --------  j  (w * r) - 1  =  phi  (w * r) - 1      
   !                        k   k                    k                
   !                 (w * r)                                             
   !                                                                
   ! and                                                                
   !                      k+1                                          
   !                  (w * r)                                             
   !   bessel_n  =  - --------  n  (w * r) - 1  =  psi  (w * r) - 1       
   !                  (2k-1)!!   k                    k             
   !                                                                  
   ! where j(x) and n(x) are spherical Bessel functions, and omega is the frequency factor. The writeup (B J  McKenzie, 
   ! I P Grant, and P H Norrington, Computer Phys Commun 21 (1980) 233-246) is incorrect in its description of the output 
   ! of this routine.                                                      
   !                                                                  
   ! The routine uses equations given in M Abramowitz and I A STegun to evaluate the functions. Devices are used to reduce 
   ! the number of actual evaluations of these functions. 
   !
   ! This routine is similar to Bessel of GRASP92 (written by F A Parpia) but has been modified and adapted to the 
   ! Fortran 90/95 standard.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                                :: k
      real(kind=dp), intent(in)                          :: omega
      real(kind=dp), dimension(1:n_grasp92), intent(out) :: bessel_j, bessel_n
      !
      integer       :: i, iswap, irem, j, jchan, nn
      real(kind=dp) :: b, cn, dfnm, dfn, obwa, s1, s2, skeep, sn, scn, ssn, wa, xbess1, xbess2
      real(kind=dp), parameter :: eps = 1.0e-6_dp
      !
      bessel_j(:) = zero;   bessel_n(:) = zero
      nn = k
      !
      ! If omega < eps10  neglect the calculation
      if (omega < eps10) return
      !
      ! Use a four-term power series for low w*r
      do  5   j = 2,n_grasp92
         wa     = - half *(r_grasp92(j) * omega)**2
         xbess1 = one;  s1 = zero 
         xbess2 = one;  s2 = zero
         do  i = 1,4
            xbess1 = xbess1 * wa  / (i*(2*(nn+i)+one))
            xbess2 = xbess2 * wa /  (i*(2*(i-nn)-one))
            s1 = s1 + xbess1;   s2 = s2 + xbess2
            if (abs(xbess1) < abs(s1)*eps  .and.  abs(xbess2) < abs(s2)*eps) then
               bessel_j(j) = s1
               bessel_n(j) = s2
               goto 5
            end if
         end do
         jchan = j
         goto 6
    5 continue
      !
      ! If here then calculated whole array using four-term power series. 
      return
      !
    6 continue
      !
      ! Use sin/cos expansion when power series requires more than four terms terms to converge
      if (nn == 0) then
         dfnm = one;   dfn = one
      else
         dfnm = one
         do  i = 3,2*nn-1,2
            dfnm = dfnm * i
         end do
         dfn = dfnm * (2*nn+one)
      end if
      dfnm = one / dfnm
      irem = mod(nn,4)
      !
      select case(irem)
      case(1)  ! nn = 1, 5, 9, ...
         ssn   = -one;   scn =  one;   iswap = 1
      case(2)  ! nn = 2, 6, 10, ...
         ssn   = -one;   scn = -one;   iswap = 0
      case(3)  ! nn = 3, 7, 11, ...
         ssn   =  one;   scn = -one;   iswap = 1
      case default ! nn = 0, 4, 8,...
         ssn   =  one;   scn =  one;   iswap = 0
      end select
      !
      do  j = jchan,n_grasp92
         wa = omega * r_grasp92(j)
         if (iswap == 0) then
            sn = ssn * sin(wa);   cn = scn * cos(wa)
         else
            sn = ssn * cos(wa);   cn = scn * sin(wa)
         endif
         obwa = one / wa
         b    = obwa
         s1   = b * sn;   s2 = b*cn
         do  i = 1,nn
            skeep = sn
            sn    = cn
            cn    = -skeep
            b     =  b * obwa * ((nn+i)*(nn-i+one)) / (two*i)
            s1    = s1 + b*sn
            s2    = s2 + b*cn
         end do
         s1 = s1 * dfn/ (wa**nn) - one
         s2 = s2 * (wa**(nn+1))* dfnm - one
         bessel_j(j) = s1
         bessel_n(j) = s2
      end do
      !
   end subroutine Bessel_grasp92
   !
   !
   function debye_integral_explicit(nu,lambda,rwf_a,rwf_b,rwf_c,rwf_d)      result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is the Debye-Slater integral
   !
   !                            k
   !                           R_Debye (abcd)
   !
   ! as defined in PLASMA. In contrast to the slater_integral_grasp92(), here an explicit 2-dim integration is carried out. 
   ! This is much more time consuming but allows for cases where the radial functions do not allow the solution of an ODE.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: nu
      real(kind=dp), intent(in)          :: lambda
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: value
      !
      integer       :: i, j, mac, mbd
      real(kind=dp) :: rr, rac, rbd
      real(kind=dp), dimension(1:n_grasp92+10) :: rho_ac, rho_bd, ta, tb, tc
      !
      rho_ac(:) = zero;  mac = min(rwf_a%mtp, rwf_c%mtp) - 2
      rho_ac(2:mac) = rwf_a%P(2:mac)*rwf_c%P(2:mac) + rwf_a%Q(2:mac)*rwf_c%Q(2:mac)
      !
      rho_bd(:) = zero;  mbd = min(rwf_b%mtp, rwf_d%mtp) - 2
      rho_bd(2:mbd) = rwf_b%P(2:mbd)*rwf_d%P(2:mbd) + rwf_b%Q(2:mbd)*rwf_d%Q(2:mbd)
      !
      ta(:) = zero
      do  i = 2,min(mac,mbd)
         rr  = r_grasp92(i)
         rac = rho_ac(i);   rbd = rho_bd(i)
         !
         do  j = 2,i
            tb(j)  = rac * U_L_rs(nu,lambda,r_grasp92(j),rr) * rho_bd(j)*rp_grasp92(j)
            tc(j)  = rbd * U_L_rs(nu,lambda,r_grasp92(j),rr) * rho_ac(j)*rp_grasp92(j)
         end do
         !
         tb(i)    = rac * U_L_rs(nu,lambda,r_grasp92(i),rr) * rho_bd(i) * rp_grasp92(i) * half
         tb(i+1:) = zero
         !
         tc(i)    = rbd * U_L_rs(nu,lambda,r_grasp92(i),rr) * rho_ac(i) * rp_grasp92(i) * half
         tc(i+1:) = zero
         !
         ta(i) = quad_grasp92(tb,i) + quad_grasp92(tc,i)
      end do
      !
      ta(1:mac) = rp_grasp92(1:mac) * ta(1:mac)  
      value = quad_grasp92(ta,mac)
      !
      !
      contains
         !
         function U_L_rs(L,lambda,s,r)                           result(U)
         !--------------------------------------------------------------------------------------------------------------
         ! This internal function calculates the U_L (r,s) functions for s <= r.
         !--------------------------------------------------------------------------------------------------------------
         !
         integer, intent(in)       :: L
         real(kind=dp), intent(in) :: lambda, r, s
         real(kind=dp)             :: U, sum, suma
         !
         integer                   :: p, q
         !
         !! U = (s**L) / (r**(L+1))
         !
         !!x lambda = 0.01
         !
         sum = zero
         do  p = 0,2
            do  q = 0,L
               sum = sum + (two**(L-q)) * (lambda**(L+p+p-q)) *  factorial_dp(L+q) * factorial_dp(L+p) /          &
                           ( factorial_dp(L+L+p+p+1) * factorial_dp(L-q) * factorial_dp(p) * factorial_dp(q)) *   & 
                           (s**(L+p+p)) * exp(-lambda*r) / (r**(q+1))
            end do
            if (p == 2) suma = sum
         end do
         !
         U =  (L+L+one) * sum
         !
         end function U_L_rs
         !
   end function debye_integral_explicit
   !
   !
   subroutine dpbdt_grasp92(rwf,P_prime,Q_prime)
   !--------------------------------------------------------------------------------------------------------------------
   ! Computes H times the derivative, with respect to the internal grid, of the large and small components of a the 
   ! radial wave function rwf. These are returned in the arrays P_prime(:) and Q_prime which must have a dimension of at 
   ! least rwf%mtp points. A  thirteen-point  Lagrange formaula is used for the calculation of derivatives.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(in)       :: rwf
      real(kind=dp), dimension(:), intent(out) :: P_prime, Q_prime
      !
      integer       :: i, irow, k, loc
      real(kind=dp) :: a1, a2, a3, a4, a5, a6, aik, hdpbdt, hdqbdt
      !
      a1 = a13_grasp92(7,1);   a2 = a13_grasp92(7,2);   a3 = a13_grasp92(7,3)   
      a4 = a13_grasp92(7,4);   a5 = a13_grasp92(7,5);   a6 = a13_grasp92(7,6)
      !
      ! Compute derivative in three separate regions; first, points 1 to 6
      do  i = 1,6
         hdpbdt = zero;   hdqbdt = zero
         do  k = 1,13
            aik = a13_grasp92(i,k)
            hdpbdt = hdpbdt + a13_grasp92(i,k) * rwf%P(k)
            hdqbdt = hdqbdt + a13_grasp92(i,k) * rwf%Q(k)
         end do
	 P_prime(i) = hdpbdt;   Q_prime(i) = hdqbdt
      end do
      !
      ! Next, points 7 to mtp-6; special treatment for this region because of the symmetry of the differentiation formula
      do  i = 7,rwf%mtp-6
         P_prime(i) =  a1 * (rwf%P(i-6) - rwf%P(i+6))  +  a2 * (rwf%P(i-5) - rwf%P(i+5)) &
                     + a3 * (rwf%P(i-4) - rwf%P(i+4))  +  a4 * (rwf%P(i-3) - rwf%P(i+3)) &
                     + a5 * (rwf%P(i-2) - rwf%P(i+2))  +  a6 * (rwf%P(i-1) - rwf%P(i+1)) 
         Q_prime(i) =  a1 * (rwf%Q(i-6) - rwf%Q(i+6))  +  a2 * (rwf%Q(i-5) - rwf%Q(i+5)) &
                     + a3 * (rwf%Q(i-4) - rwf%Q(i+4))  +  a4 * (rwf%Q(i-3) - rwf%Q(i+3)) &
                     + a5 * (rwf%Q(i-2) - rwf%Q(i+2))  +  a6 * (rwf%Q(i-1) - rwf%Q(i+1)) 
      end do
      !
      ! Last, points mtp-5 to mtp
      do  i = rwf%mtp-5,rwf%mtp
         irow = i-rwf%mtp+13
         hdpbdt = zero;   hdqbdt = zero
         do  k = 1,13
            loc    = rwf%mtp - 13 + k
            hdpbdt = hdpbdt + a13_grasp92(irow,k) * rwf%P(loc)
            hdqbdt = hdqbdt + a13_grasp92(irow,k) * rwf%Q(loc)
         end do
	 P_prime(i) = hdpbdt;   Q_prime(i) = hdqbdt
      end do
      !
   end subroutine dpbdt_grasp92
   !
   !
   subroutine hydrogen_point_grasp92(rwf_hydrogen,Z)   
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates a hydrogen orbital function for a point nucleus with charge Z. This subroutine computes the Dirac-Coulomb 
   ! bound-state orbital radial Berestetskii modified to ensure positive slope at the origin wavefunction. Equations (13.5) 
   ! and (13.5') of  Akhiezer and Berestetskii modified to ensure positive slope at the origin for RG are used.
   ! This routine has been taken from GRASP92, DCBSRW which was written by F A Parpia; it has beed adapted to the Fortran 
   ! 90/95 standard.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(inout) :: rwf_hydrogen
      real(kind=dp), intent(in)             :: Z
      !
      integer       :: i,iordr1, iordr2, k, mtp, nr, nrfac
      real(kind=dp) :: a, alfa, an1, an2, argr, argi, b, bn, bign, bignmk, eps, fac, facn, fden, ff, fg, fk, f1, f2,   &
                       gamma, ovlfac, rgamm1, rgamm2, rg0, cutoff, rho, rhon, twogp1, zalfa
      real(kind=dp), dimension(1:n_grasp92) :: ta, tb
      !
      ! Ensure that the principal quantum number, kappa, and the
      ! nuclear charge is physical
      if (rabs_use_stop) then
      if (rwf_hydrogen%orbital%n < 0) then
         print *, " "
         print *, "Principal quantum number is ",rwf_hydrogen%orbital%n
         stop  "hydrogen_point_grasp92(): program stop A."
      else if (rwf_hydrogen%orbital%kappa == 0   .or.                                         &
               angular_momentum_l(rwf_hydrogen%orbital%kappa) >  rwf_hydrogen%orbital%n-1) then
         print *, " "
         print *, "Principal quantum number and kappa are ", rwf_hydrogen%orbital%n, rwf_hydrogen%orbital%kappa
         stop  "hydrogen_point_grasp92(): program stop B."
      else if (Z < zero   .or.   Z > c) then
         print *, " "
         print *, "Nuclear charge and c are ",Z,c
         stop  "hydrogen_point_grasp92(): program stop C."
      end if
      end if
      !
      ! Atomic units
      alfa = one / c
      !
      ! Now determine all the parameters
      nr     = rwf_hydrogen%orbital%n-abs(rwf_hydrogen%orbital%kappa)
      fk     = abs(rwf_hydrogen%orbital%kappa)
      zalfa  = Z*alfa
      gamma  = sqrt(rwf_hydrogen%orbital%kappa*rwf_hydrogen%orbital%kappa - zalfa*zalfa)
      twogp1 = gamma + gamma + one
      bign   = sqrt(rwf_hydrogen%orbital%n*rwf_hydrogen%orbital%n &
                    -two*(rwf_hydrogen%orbital%n-fk)*(fk-gamma))
      eps    = one / sqrt(one + (zalfa/(gamma+(rwf_hydrogen%orbital%n-fk)))**2)
      !
      ! EPS is the total energy divided by C*C; this must be converted to the units and reference energy of GRASP
      rwf_hydrogen%energy = (one-eps)*c*c
      !
      ! Now the normalization constants
      nrfac  = 1
      do  i = 1,nr
         nrfac = nrfac * i
      end do
      argr   = twogp1 + nr
      argi   = zero
      rgamm1 = gamma_complex(cmplx(argr,argi,kind=dp))
      argr   = twogp1
      rgamm2 = gamma_complex(cmplx(argr,argi,kind=dp))
      !
      fac = - sqrt(rgamm1) / (rgamm2*sqrt(dble(nrfac))) * sqrt (Z/(two*bign*bign*(bign-rwf_hydrogen%orbital%kappa)))
      !
      ! Ensure that the slope of the large-component function is positive
      ! at the origin
      if (rwf_hydrogen%orbital%kappa > 0) then
         fac = -fac
      end if
      !
      fg = fac * sqrt(one + eps)
      ff = fac * sqrt(one - eps)
      !
      ! Now set up the coefficients of the confluent hypergeometric
      ! functions  F (-NR+1,2*GAMMA+1;RHO)  and  F (-NR,2*GAMMA+1;RHO)
      ! in the workspace arrays  TA  and  TB , respectively
      if (nr == 0) then
         iordr1 = 0;        iordr2 = 0
      else
         iordr1 = nr - 1;   iordr2 = nr
      end if
      !
      fac  = one;   facn = one
      a    = -nr
      an1  = a + one;   an2  = a
      b    = twogp1;    bn   = b
      !
      k = 0
    2 k = k+1
      fden = one / (facn*bn)
      if (k <= iordr1) then
         ta(k) = an1 * fden
      end if
      !
      if (k <= iordr2) then
         tb(k) = an2 * fden
         a     = a + one
         an1   = an1 * (a+one)
         an2   = an2 * a
         b     = b + one
         bn    = bn * b
         fac   = fac + one
         facn  = facn * fac
         goto 2
      end if
      !
      ! Now tabulate the function over the entire grid
      rwf_hydrogen%P(1)  = zero
      rwf_hydrogen%Q(1)  = zero
      fac    = (Z+Z)/bign
      bignmk = bign - rwf_hydrogen%orbital%kappa
      do  i   = 2,n_grasp92
         rho  = fac * r_grasp92(i)
         rhon = rho
         k = 0
         f1 = one;   f2 = one
    3    k = k+1
         if (k <= iordr1) then
            f1  = f1 + ta(k) * rhon
         end if
         if (k  <= iordr2) then
            f2   = f2 + tb(k) * rhon
            rhon = rhon * rho
            goto 3
         end if
         f1 = nr * f1
         f2 = bignmk * f2
         ovlfac = exp(-half*rho) * (rho**gamma)
         rwf_hydrogen%P(i)  = fg * ovlfac * (f1-f2)
         rwf_hydrogen%Q(i)  = ff * ovlfac * (f1+f2)
      end do
      !
      ! Determine the effective maximum tabulation point based on the
      ! cutoff; define the cutoff conservatively
      cutoff = accy_grasp92/ten
      !
      mtp = n_grasp92+1
    5 mtp = mtp-1
      if (abs(rwf_hydrogen%P(mtp)) < cutoff) then
         rwf_hydrogen%P(mtp) = zero
         rwf_hydrogen%Q(mtp) = zero
         goto 5
      end if
      rwf_hydrogen%mtp = mtp
      !
      if (mtp == n_grasp92) then
         print *, "Warning in hydrogen_point_grasp92(): Radial grid of insufficient extent; P(",n_grasp92,") = ",              &
                  rwf_hydrogen%P(n_grasp92)
         print *, "Cutoff (",cutoff,") exceeded."
      end if
      !
      ! Compute the coefficient of r**gamma at the origin
      rg0 = fg * (fac**gamma) * (nr-bignmk)
      !
   end subroutine hydrogen_point_grasp92
   !
   !
   subroutine  hydrogenic_bound_orbital(n,kappa, Znuc,energy, rg0,rg,rf, mtp)              
   !--------------------------------------------------------------------------------------------------------------------
   ! This subroutine computes the  Dirac-Coulomb  bound-state orbital radial wavefunction. Equations (13.5) and 
   ! (13.5') of  Akhiezer and Berestetskii modified to ensure positive slope at the origin for RG are used. 
   !
   ! The arguments are as follows: 
   !
   !    N    : (Input)   The (usual) principal quantum number 
   !    KAPPA: (Input)   The relativistic angular quantum number 
   !    Z    : (Input)   The nuclear charge
   !    E    : (Output)  The Dirac-Coulomb Eigenenergy  
   !    RG0  : (Output)  Coefficient of the leading term in the series expansion of the large component near the origin
   !    RG   : (Output)  r times the large component wavefunction of Akhiezer and Berestetskii 
   !    RF   : (Output)  r times the small component wavefunction of Akhiezer and Berestetskii 
   !    MTP  : (Output)  Maximum tabulation point  
   !
   ! This procedures was originally written by Farid A Parpia but has been adapted to the present standard.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)        :: n, kappa
      integer, intent(out)       :: mtp
      real(kind=dp), intent(in)  :: Znuc
      real(kind=dp), intent(out) :: energy, rg0
      real(kind=dp), dimension(:), intent(out) :: rg, rf
      !
      integer          :: i, k, nr, iordr1, iordr2, nrfac
      real(kind=dp)    :: alpha, bign, eps, fk, fn, fkappa, fnr, gamma, twogp1, zalfa, rgamm1, rgamm2, fac, facn, fg, ff, &
                          a, an1, an2, b, bn, fden, bignmk, rho, rhon, f1, f2, ovlfac, cutoff, accy_grasp92
      complex(kind=dp)                    :: za, zb
      real(kind=dp), dimension(n_grasp92) :: ta, tb
      !
      ! Ensure that the principal and angular quantum numbers are physical
      if (n <= 0) then
         print *, "Principal quantum number is n = ",n
         stop "hydrogenic_bound_orbital(): program stop A."
      else if (kappa == 0  .or.  kappa == n  .or. abs(kappa) > n) then
         print *, "Improper kappa quantum number; n, kappa = ",n,kappa
         stop "hydrogenic_bound_orbital(): program stop B."
      else if (Znuc <= zero  .or.  Znuc > c) then
         print *, "Nuclear charge is to small or exceeds c ;Z, c = ",Znuc, c
         stop "hydrogenic_bound_orbital(): program stop C."
      end if
      !
      ! Use atomic units and determine all parameters
      alpha = one / c
      !
      fn  = n;   fkappa = kappa;   k = abs(kappa);   fk = k;   nr = n - k
      fnr = nr;  zalfa = Znuc * alpha
      gamma  = sqrt(fk*fk - zalfa*zalfa)
      twogp1 = gamma + gamma + one
      bign   = sqrt (fn*fn - two*fnr*(fk-gamma))
      eps    = one / sqrt(one + (zalfa/(gamma+fnr))**2 )
      !
      ! eps  is the total energy divided by C*C; this must be converted to the units and reference energy of GRASP
      energy = (one - eps) *c*c
      !
      ! Now the normalization constants
      nrfac = 1
      do  i = 1,nr
         nrfac = nrfac * i
      end do
      !
      za      = cmplx(twogp1+fnr, zero)
      zb      = gamma_complex(za)
      rgamm1  = zb
      !
      za      = cmplx(twogp1, zero)
      zb      = gamma_complex(za)
      rgamm2  = zb
      ! 
      fac = - sqrt (rgamm1)/ (rgamm2 *sqrt( dble(nrfac) )) * sqrt (Znuc/(two*bign*bign*(bign-fkappa)))
      !
      ! Ensure that the slope of the large-component function is positive at the origin
      if (kappa > 0) fac = -fac
      fg = fac * sqrt(one+eps)
      ff = fac * sqrt(one-eps)
      !
      ! Now set up the coefficients of the confluent hypergeometric functions  F (-NR+1,2*GAMMA+1;RHO)  and  
      ! F (-NR,2*GAMMA+1;RHO) in the workspace arrays  TA  and  TB , respectively
      if (nr == 0) then
         iordr1 = 0;    iordr2 = 0
      else
         iordr1 = nr-1; iordr2 = nr
      endif
      !
      fac  = one
      facn = one
      a    = -fnr
      an1  = a + one
      an2  = a
      b    = twogp1
      bn   = b
      !
      k = 0
    2 k = k+1
      fden = one/(facn*bn)
      if (k <= iordr1) then
         ta(k) = an1 * fden
      endif
      !
      if (k <= iordr2) then
         tb(k) = an2 * fden
         a    = a + one
         an1  = an1 * (a+one)
         an2  = an2 * a
         b    = b + one
         bn   = bn * b
         fac  = fac + one
         facn = facn * fac
         goto 2
      endif
      !
      ! Now tabulate the function over the entire grid
      rg(1) = zero
      rf(1) = zero
      fac   = (Znuc + Znuc)/bign
      bignmk = bign-fkappa
      do  i = 2,n_grasp92
         rho  = fac * r_grasp92(i)
         rhon = rho
         k    = 0
         f1   = one
         f2   = one
    3    k = k+1
         if (k <= iordr1) then
            f1 = f1 + ta(k) * rhon
         endif
         if (k <= iordr2) then
            f2   = f2 + tb(k) * rhon
            rhon = rhon * rho
            goto 3
         endif
         f1 = fnr * f1
         f2 = bignmk * f2
         ovlfac = exp (-half*rho)*(rho**gamma)
         rg(i) = fg * ovlfac * (f1-f2)
         rf(i) = ff * ovlfac * (f1+f2)
      end do
      !
      ! Determine the effective maximum tabulation point based on the cutoff; define the cutoff conservatively
      accy_grasp92 = h_grasp92**6
      cutoff       = accy_grasp92 / ten
      !
      mtp = n_grasp92 + 1
    5 mtp = mtp-1
      if (abs (rg(mtp)) .lt. cutoff) then
         rg(mtp) = zero
         rf(mtp) = zero
         goto 5
      endif
      !
      if (mtp == n_grasp92) then
         print *, "*** hydrogenic_bound_orbital(): Radial grid of insufficient extent; n, kappa = ",n, kappa
         print *, "    n_grasp92, rg(n_grasp92), cutoff = ",n_grasp92, rg(n_grasp92), cutoff
      end if
      !
      ! Compute the coefficient of R**GAMMA at the origin
      rg0 = fg * (fac**gamma) * (fnr-bignmk)
      
      !
   end subroutine hydrogenic_bound_orbital
   !
   !
   function I_ab_grasp92(rwf_a,rwf_b,mode)                  result(I_ab)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the one-particle integral I (ab) for two given radial orbital functions. The analytical expression for 
   ! this quantity is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C Pyper,  
   ! Computer  Phys Commun 21 (1980) 211 .
   ! 
   ! The optional argument may be mode = "kinetic" to return only the kinetic part of the integral.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(in)     :: rwf_a, rwf_b
      real(kind=dp)                          :: I_ab
      character(len=*), intent(in), optional :: mode
      !
      integer       :: i, mtp
      real(kind=dp), dimension(n_grasp92)    :: P_prime_b, Q_prime_b
      real(kind=dp), dimension(n_grasp92+10) :: ta
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Code added to support the Ce@C_82 calculations of Galya Kashenock
      logical           :: first = .true.
      integer           :: ii
      real(kind=dp)  :: ri
      real(kind=dp), dimension(12000), save :: ta_galya
      !
      !!x	if (first) then
      !!x	   first = .false.
      !!x	   write(47,*) "Radial grid to define the additional potential"
      !!x	   write(47,*) "----------------------------------------------"
      !!x	   write(47,*) " "
      !!x	   do  i = 1,n_grasp92
      !!x	      write(47,*) i, r_grasp92(i)
      !!x	   end do
      !!x	end if
      !!x	!
      !!x	if (first) then
      !!x	   first = .false.
      !!x	   do  i = 1,n_grasp92
      !!x	      read(48,*) ii, ri, ta_galya(i)
      !!x	   end do
      !!x	  print *, "********"
      !!x	  ii = 3200
      !!x	  print *, "ii, ta_galya(ii), ta_galya(ii+400), ta_galya(ii+800) ", &
      !!x		 ii, ta_galya(ii), ta_galya(ii+400), ta_galya(ii+800)
      !!x	end if
      !!x	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Stop if orbitals rwf_a and rwf_b have different kappa values
      if (rabs_use_stop                      .and.  &
          rwf_a%orbital%kappa /= rwf_b%orbital%kappa) then
         print *, "I_ab_grasp92(): Attempt to calculate I("//              &
	          orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)//","// &
		  orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//")."
         stop "I_ab_grasp92(): program stop A."
      end if
      !
      ! Kinetic energy contribution; first, piece involving derivatives
      mtp = min( rwf_a%mtp, rwf_b%mtp )
      call dpbdt_grasp92(rwf_b,P_prime_b,Q_prime_b)
      ta(1)     = zero
      ta(2:mtp) = rwf_a%Q(2:mtp) * P_prime_b(2:mtp) - &
                  rwf_a%P(2:mtp) * Q_prime_b(2:mtp)
      I_ab = c* quad_grasp92(ta,mtp) / h_grasp92
      !
      ! Pieces not involving derivatives
      ta(1)     = zero
      ta(2:mtp) = rp_grasp92(2:mtp) * rwf_a%Q(2:mtp) * rwf_b%Q(2:mtp)
      I_ab = I_ab - two * c*c * quad_grasp92(ta,mtp)
      !
      ta(1)     = zero
      ta(2:mtp) = rpor_grasp92(2:mtp) * (rwf_a%P(2:mtp) * rwf_b%Q(2:mtp)  &
                                       + rwf_a%Q(2:mtp) * rwf_b%P(2:mtp))
      I_ab = I_ab + c * rwf_b%orbital%kappa * quad_grasp92(ta,mtp)
      !
      if (present(mode)) then
         if (mode(1:7) == "kinetic") then
         if (debug_I_ab_grasp92) then
          !  write(99,*)                                                   &
	  !  " K_ab("//orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)// &
	  !  ","//orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//") = ",I_ab
         end if
         return
         end if
      end if
      !
      ! Contribution from nuclear potential
      ta(1) = zero
      do  i = 2,mtp
         ta(i) = rp_grasp92(i) * nuclear_potential(r_grasp92(i)) * (rwf_a%P(i) * rwf_b%P(i) + rwf_a%Q(i) * rwf_b%Q(i))
      end do
      !
      !!x	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!x	! Code added to support the Ce@C_82 calculations of Galya Kashenock
      !!x	do  i = 2,mtp
      !!x	   ta(i) = rp_grasp92(i) *					  &
      !!x		 (nuclear_potential(r_grasp92(i))+three*ta_galya(i)) *    &
      !!x		 (rwf_a%P(i) * rwf_b%P(i) + rwf_a%Q(i) * rwf_b%Q(i))
      !!x	end do
      !!x	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      I_ab = I_ab + quad_grasp92(ta,mtp)
      !
      if (debug_I_ab_grasp92) then
         write(99,*) " I_ab("//orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)// &
	             ","//orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//") = ",I_ab
      end if
      !
   end function I_ab_grasp92
   !
   !
   subroutine interpolate_rwf_grasp92(Pa,Qa,ma,ra,P_new,Q_new,mtp_new)
   !--------------------------------------------------------------------------------------------------------------------
   ! This  subprogram interpolates the arrays  Pa(1:ma), Qa(1:ma), tabulated on grid ra(1:ma) into the arrays P_new(), 
   ! Q_new(). (Aitken's algorithm is used. See F B Hildebrand, Introduction  to Numerical  Analysis, 2nd ed., McGraw-Hill, 
   ! New York, NY, 1974.) The orbital is renormalized.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: ma
      integer, intent(out) :: mtp_new
      real(kind=dp), dimension(ma), intent(in) :: Pa, Qa, ra
      real(kind=dp), dimension(:), pointer     :: P_new, Q_new
      character(len=256) :: buffer      
      !
      integer       :: i, ilirok, ildiag, ilothr, irow, k, kount, llo, lhi, locnxt, mfj, nrstlo, nrsthi
      real(kind=dp) :: diff, difft, dpbp, dqbq, dxkmn1, dxirow, factor, pestl, pestt, qestl, qestt, rama, rn, xbar
      logical       :: set
      logical, dimension(n_grasp92+10000) :: used
      !
      ! mxord is the maximum order of the interpolation
      integer, parameter :: mxord = 13
      real(kind=dp), dimension(mxord)               :: x, dx
      real(kind=dp), dimension((mxord*(mxord+1))/2) :: polyp, polyq
      !
      ! Initialization
      rama = ra(ma)
      rn   = r_grasp92(n_grasp92)
      !
      !!x print *, "***interpolate_rwf_grasp92*** rama,rn = ",rama,rn
      ! Checks
      if (rama > rn) then
         write (*,1) rn, rama
         stop
       1 format(/"interpolate_rwf_grasp92(): Grid of insufficient extent:", &
                /" present grid has r(n) = ",1p,1d19.12," Bohr radii",      &
                /"          require r(n) = ",   1d19.12," Bohr radii" )
      end if
      !
      ! Determine mtp_new
      i = n_grasp92
    2 i = i - 1
      if (r_grasp92(i) <= rama) then
         mfj = i
      else
         goto 2
      end if
      mtp_new = mfj
      !
      ! Allocate the arrays P_new and Q_new
      allocate( P_new(1:mtp_new), Q_new(1:mtp_new) )
      !
      P_new(1) = zero   ! this is always true in GRAP92
      Q_new(1) = zero
      !
      ! Overall initialization for interpolation
      nrstlo = 0
      kount  = 0
      !
      ! Perform interpolation
      do  i = 2,mfj
         !
         ! Initialization for interpolation
         xbar  = r_grasp92(i);   irow  = 0
         pestl = zero;           qestl = zero
         !
         ! Determine the nearest two grid points bounding the present grid point
    3    k = nrstlo + 1
         if (ra(k) < xbar) then
            nrstlo = k
            goto 3
         else
            nrsthi = k
         end if
         !
         ! Clear relevant piece of use-indicator array
         llo = max(nrstlo-mxord,  1)
         lhi = min(nrsthi+mxord, ma)
         do  k = llo,lhi
            used(k) = .false.
         end do
         !
         ! Determine next-nearest grid point
    4    irow = irow + 1
         llo  = max (nrstlo-irow+1, 1)
         lhi  = min (nrsthi+irow-1,ma)
         set  = .false.
         do  k = llo,lhi
            if (.not. used(k)) then
               if (.not. set) then
                  diff   = ra(k) - xbar
                  locnxt = k
                  set    = .true.
               else
                  difft  = ra(k) - xbar
                  if (abs(difft) < abs(diff)) then
                     diff   = difft
                     locnxt = k
                  end if
               end if
            end if
         end do
         used(locnxt) = .true.
         x(irow)      = ra(locnxt)
         dx(irow)     = diff
         !
         ! Fill table for this row
         do  k = 1,irow
            ilirok = iloc(irow,k)
            if (k == 1) then
               polyp(ilirok) = pa(locnxt)
               polyq(ilirok) = qa(locnxt)
            else
               ildiag = iloc(k-1,k-1)
               ilothr = iloc(irow,k-1)
               dxkmn1 = dx(k-1)
               dxirow = dx(irow)
               factor = one / (x(irow) - x(k-1))
               polyp(ilirok) =  ( polyp(ildiag) * dxirow  - polyp(ilothr) * dxkmn1 ) * factor
               polyq(ilirok) =  ( polyq(ildiag) * dxirow  - polyq(ilothr) * dxkmn1 ) * factor
            end if
         end do
         !
         ! Check for convergence
         ildiag = iloc (irow,irow)
         pestt  = polyp(ildiag)
         qestt  = polyq(ildiag)
         if (pestt == zero   .or.   qestt == zero) then
            if (irow < mxord) then
               goto 4
            else
               P_new(i) = pestt
               Q_new(i) = qestt
            end if
         else
            dpbp = abs((pestt-pestl)/pestt)
            dqbq = abs((qestt-qestl)/qestt)
            if (dqbq < accy_grasp92  .and.  dpbp < accy_grasp92) then
               P_new(i) = pestt
               Q_new(i) = qestt
            else
               pestl    = pestt
               qestl    = qestt
               if (irow .lt. mxord) then
                  goto 4
               else
                  P_new(i) = pestt
                  Q_new(i) = qestt
                  kount    = kount+1
               end if
            end if
         end if
      end do
      !
      if (kount > 0) then
         write (99,5) accy_grasp92,kount,mfj
       5 format (/"interpolate_rwf_grasp92(): Interpolation procedure not converged to",    &
                 1p,1d19.12," for ",1i3," of ",1i3," tabulation points")
      end if
      !
      !
      contains
         !
         function iloc (ind1,ind2)                           result(loc)
         !--------------------------------------------------------------------------------------------------------------
         ! This internal function dispenses with the need for a two-dimensional array for the interpolation. It replaces 
         ! astatement function in the original code.
         !--------------------------------------------------------------------------------------------------------------
         !
         integer, intent(in) :: ind1, ind2
         integer             :: loc
         !
         loc = (ind1*(ind1-1)) / 2 + ind2
         !
         end function iloc
         !
   end subroutine interpolate_rwf_grasp92
   !
   !
   function is_bound_orbital(rwf)                      result(yes)
   !--------------------------------------------------------------------------------------------------------------------
   ! This logical function returns .true. if rwf represents a bound orbital and .false. otherwise.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(in) :: rwf
      logical :: yes
      !
      integer       :: i
      real(kind=dp) :: P, Q
      !
      P = zero;   Q = zero
      do  i = rwf%mtp,rwf%mtp-10,-1
         P = P + abs(rwf%P(i));   Q = Q + abs(rwf%Q(i));   
      end do
      !
      if (P < 1.0e-5   .and.   Q < 1.0e-5) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function is_bound_orbital
   !
   !
   subroutine load_rwf_file_grasp92(wave,is_formatted,ierr, notall_needed)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in a set of GRASP92 orbital functions from stream 21 into a derived data structure of type(grasp92_orbital).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(out)                 :: ierr
      logical, intent(in)                  :: is_formatted     
      logical, intent(in), optional        :: notall_needed   
      type(grasp92_orbital), intent(inout) :: wave
      !
      integer       :: i, ios, j, nwin, n, kappa, mtp
      real(kind=dp) :: energy, norm, pz
      real(kind=dp), dimension(:), pointer :: P, Q, ra 
      !
      ! Read orbital information from read orbitals file; write summary to  .dbg  file if option set
      if (debug_load_rwf_grasp92) then
         write (99,2)
       2 format(/,"From subroutine load_rwf_file_grasp92():", &
                /," Orbital",8x,"Eigenvalue",19x,"Norm")
      end if
      nwin = 0
    1 if (is_formatted) then
         read (21,*,iostat = ios) n,kappa,energy,mtp
      else
         read (21,iostat = ios)   n,kappa,energy,mtp
      end if
      !
      !! print *, "n,kappa,energy,mtp = ",n,kappa,energy,mtp
      !
      if (ios == 0) then
         allocate( P(1:mtp), Q(1:mtp), ra(1:mtp) )
         if (is_formatted) then
            read (21,*) pz
            do  i = 1,mtp
               read (21,*) ra(i),P(i),Q(i)
            end do
         else
            read (21)   pz, (P(i),i = 1,mtp),(Q(i),i = 1,mtp)
            read (21)   (ra(i),i = 1,mtp)
         end if
         !
         do  j = 1,wave%number_of_rwf
            !!x print *, "*** dd, j = ",j
            if (wave%rwf(j)%mtp == 0  .and.  wave%rwf(j)%orbital%n == n  .and.  wave%rwf(j)%orbital%kappa == kappa) then
               wave%rwf(j)%pz     = pz 
               wave%rwf(j)%energy = energy
               call interpolate_rwf_grasp92(P,Q,mtp,ra, wave%rwf(j)%P,wave%rwf(j)%Q,wave%rwf(j)%mtp)
               !!x print *, "*** dd"
               !
               ! Normalization
               norm = rk_integral_grasp92(wave,0,j,j)
               norm = sqrt(norm)
               if(.not. norm .eq. 0) then
                  do  i = 1,wave%rwf(j)%mtp
                      wave%rwf(j)%P(i) = wave%rwf(j)%P(i) / norm 
                      wave%rwf(j)%Q(i) = wave%rwf(j)%Q(i) / norm 
                  end do
               end if
               if (debug_load_rwf_grasp92) then
                  write (99, "(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)") wave%rwf(j)%orbital%n,   &
                             orbital_symmetry(wave%rwf(j)%orbital%kappa), wave%rwf(j)%energy,norm
               end if
               !!x print *, "*** ee"
               nwin = nwin+1
            endif
         end do
         deallocate( P, Q, ra )
         !!x print *, "*** ff"
         goto 1
      end if
      if (debug_load_rwf_grasp92) write (99,*) " Orbitals renormalised;"
      !
      ! Return if not all orbitals are needed
      if (present(notall_needed)  .and.  notall_needed) return
      !
      !
      ! Stop with an error message if all orbitals are not known
      if (nwin < wave%number_of_rwf) then
         print *, "load_rwf_file_grasp92(): All required orbitals not found."
         ierr = 1
         return
      else
         ierr = 0
      end if
      !
      ! Schmidt orthogonalise the orbitals
!       call schmidt_orthogonalize_grasp92(wave)
      if (debug_load_rwf_grasp92) then
         write (99,*) " orbitals orthogonalised and renormalised;"
      end if
      print *, " ... load complete;"
      !
   end subroutine load_rwf_file_grasp92
   !
   !
   subroutine ncharg_grasp92(mtp)
   !--------------------------------------------------------------------------------------------------------------------
   ! Evaluates the nuclear charge density, and stores it in the array zdist_grasp92(:).
   !
   ! Calls: nuclear_potential_es
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(out) :: mtp
      !
      integer              :: i
      logical              :: form1, form2
      real(kind=dp)        :: abc, abc2, abc3, cba, en, extrm, pi2, &
                              s2mcba, s3mcba, zdisti, znorm
      !
      ! Allocate the array zdist_grasp92(:)
      allocate( zdist_grasp92(1:n_grasp92) )
      zdist_grasp92(:) = zero
      !
      ! Fermi charge distribution
      if (nuclear_model == "fermi") then
         cba   = fermi_c_parameter / fermi_a_parameter
         abc   = fermi_a_parameter / fermi_c_parameter
         abc2  = abc*abc
         abc3  = abc2*abc
         pi2   = pi*pi
         call nuclear_potential_es(-cba,s2mcba,s3mcba)
         en    = one + pi2*abc2 - 6.0_dp*abc3*s3mcba
         znorm = three * nuclear_charge / (four*pi*en*fermi_c_parameter**3)
         form1 = .true.
         form2 = .false.
         do  i = 1,n_grasp92
            if (form1) then
               extrm = exp ((r_grasp92(i)-fermi_c_parameter)/fermi_a_parameter)
               zdist_grasp92(i) = znorm /(one + extrm)
               if (one/extrm <= eps10*eps10*eps10) then
                  form1 = .false.
                  form2 = .true.
               end if
            else if (form2) then
               zdisti = znorm * exp( -(r_grasp92(i)-fermi_c_parameter) /fermi_a_parameter )
               if (abs(zdisti) > zero) then
                  zdist_grasp92(i) = zdisti
               else
                  mtp = i
                  return
               end if
            end if
         end do
      end if
      !
   end subroutine ncharg_grasp92
   !
   !
   function quad_grasp92(ta,mtp)                          result(result)
   !--------------------------------------------------------------------------------------------------------------------
   ! The argument result is an approximation  to the integral of F(R) from  zero to infinity, where the values of  
   ! RP(I)*F(R(I))  are tabulated in the array  ta(1:mtp). The integral in the interval zero to r_grasp92(j) is computed 
   ! by use of an analytical fit
   !
   !                              SIGMA
   !                   F(R) = A R 
   !
   ! A five-point Closed Newton-Cotes formula (cf. F B Hildebrand, Introduction to Numerical Analysis, second edition, 
   ! McGraw-Hill, New York, 1974, p. 93)  is  used  to  compute the integral in the interval  r(j:mtp). The  contribution  
   ! from  the  tail  of  the function beyond the last  tabular  point (mtp) is assumed to be negligible. The method 
   ! uses  mtp+3  tabulation points. Array ta should therefore be dimensioned to at least  N+4 .
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                        :: mtp
      real(kind=dp), dimension(:), intent(inout) :: ta
      real(kind=dp)                              :: result
      !
      integer       :: i, ip1, loc, mtpm1
      real(kind=dp) :: fri, frip1, quott, ratio, ri, rip1, sigma, tai, taip1
      !
      ! Find first values that will permit computation of exponent
      mtpm1 = mtp - 1
      do  i = 2,mtpm1
         tai = ta(i)
         if (abs(tai) > zero) then
            ip1   = i + 1
            taip1 = ta(ip1)
            quott = taip1 / tai
            if (quott > zero) then
               !
               ! Exponent from fit
               frip1 = taip1 / rp_grasp92(ip1)
               fri   = tai   / rp_grasp92(i  )
               ratio = frip1 / fri
               rip1  = r_grasp92(ip1)
               ri    = r_grasp92(i  )
               sigma = log(ratio) / log(rip1/ri)
               !
               ! Analytical integration and error estimate for interval r(1:i)
               fri    = ri * fri
               result = fri / (sigma+one)
               !
               ! Set the tail to zero
               do  loc = 1,3
                  ta(mtp+loc) = zero
               end do
               !
               ! Newton-Cotes quadature for the remainder
               result = result +c1_grasp92 * tai
               do  loc = ip1,mtp,4
                  result = result + c2_grasp92 *(ta(loc) + ta(loc+2)) + c3_grasp92 * ta(loc+1) + c4_grasp92 * ta(loc+3)
               end do
               if (mod(mtp-i,4) == 0) result = result - c1_grasp92 * ta(mtp)
               !
               ! Test of result's accuracy; `decomment' to activate
               !! estder = 10.0d 00*ri*fri
               !! ratio = abs (estder/result)
               !! if (ratio .gt. accy) then 
               !!    print (*,300) ratio
               !!    300 format (/"quad_grasp92(): Estimated accuracy is ",
               !!               es10.3,/, " decrease rnt or improve input ",
               !!                 "data conditioning to ameliorate."/)
               !!
               return
            end if
         end if
      end do
      !
      ! No value which will permit computation of exponent
      result = zero
      !
   end function quad_grasp92
   !
   !
   subroutine radgrd_grasp92()
   !--------------------------------------------------------------------------------------------------------------------
   ! This routine sets up the radial grid  r_grasp92  and the associated arrays  rp_grasp92  and  rpor_grasp92. Different 
   ! grids are generated depending on whether or not  hp_grasp92 is zero. It has been taken from GRASP92 (Parpia, 
   ! Froese Fischer & Grant, 1996 CPC  and adapted to Fortran 90.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer       :: i, ii, ii1, ii2, kk, np10, nb2, nrows
      real(kind=dp) :: a, delr, epslon, eph, ett, ettm1, fofr, fpri, rest, rests, rlast, rpow, t
      !
      np10 = n_grasp92 + 10
      allocate( r_grasp92(1:np10), rp_grasp92(1:np10), rpor_grasp92(1:np10) )
      !
      ! rpor_grasp92(1) is never used in the program: it is arbitrarily
      ! set to zero
      r_grasp92(1)    = zero
      rpor_grasp92(1) = zero
      !
      ! Now set up the grids
      if (hp_grasp92 == zero) then
         !
         ! Exponential grid if hp_grasp92 is zero
         rp_grasp92(1) = rnt_grasp92
         eph = exp (h_grasp92)
         ett = one
         !
         ! Set up the arrays r_grasp92, rp_grasp92, rpor_grasp92
         do  i = 2,np10
            ett   = eph * ett
            ettm1 = ett - one
            r_grasp92(i)    = rnt_grasp92 * ettm1
            rp_grasp92(i)   = rnt_grasp92 * ett
            rpor_grasp92(i) = ett/ettm1
         end do
      else
         !
         ! Asymptotically-linear exponential grid otherwise:
         epslon = 1.0e3 * epsilon(one)
         a      = h_grasp92 / hp_grasp92
         rlast  = zero
         rest   = zero
         rp_grasp92(1) = rnt_grasp92 / (a * rnt_grasp92 + one)
         !
         ! Set up the arrays r, rp, rpor
         do  i = 2,np10
            t = h_grasp92 * (i-1)
            !
            ! Solve the implicit equation for r using the Newton-Raphson method
    2       rests = rest + rnt_grasp92
            fofr  = log(rests/rnt_grasp92) + a * rest - t
            fpri  = rests /(a * rests + one)
            delr  = -fofr * fpri
            rest  =  rlast + delr
            !
            if (abs(delr/rest) < epslon) then
               r_grasp92(i) = rest
               rests = rest + rnt_grasp92
               fpri  = rests /(a * rests + one)
               rp_grasp92(i)   = fpri
               rpor_grasp92(i) = fpri/rest
            else
               rlast = rest
               goto 2
            end if
         end do
      endif
      !
      ! Initialize some further arrays for calculating Hartree Y- and Z-functions; these arrays store r' * r ; 
      ! they increase the storage required by GRASP92, but reduce the CPU effort
      allocate( rpbrtk_grasp92(2:n_grasp92,1:kmp1_grasp92),  rprtk_grasp92(2:n_grasp92,1:kmp1_grasp92),   &
                rtk_grasp92(2:n_grasp92,1:kmp1_grasp92),     rtki_grasp92(2:n_grasp92,1:kmp1_grasp92) )
      !
      do  kk = 1,kmp1_grasp92
         do  ii = 2,n_grasp92
            rpow                  = r_grasp92(ii)**kk
            rtk_grasp92(ii,kk)    = rpow
            rprtk_grasp92(ii,kk)  = rp_grasp92(ii) * rpow
            rtki_grasp92(ii,kk)   = one / rpow
            rpbrtk_grasp92(ii,kk) = rpor_grasp92(ii) / rpow
         end do
      end do
      !
      ! Debug printout
      if (debug_radgrd_grasp92) then
         write (99,300)
         nb2 = n_grasp92/2
         if (2*nb2 .eq. n_grasp92) then
            nrows = nb2
         else
            nrows = nb2 + 1
         end if
         do  ii = 1,nrows
            ii1 = ii;   ii2 = ii1+nrows
            if (ii2 <= n_grasp92) then
               write(99,301) r_grasp92(ii1),rp_grasp92(ii1),rpor_grasp92(ii1),        &
                             r_grasp92(ii2),rp_grasp92(ii2),rpor_grasp92(ii2)
            else if (ii1 <= n_grasp92) then
               write(99,301) r_grasp92(ii1),rp_grasp92(ii1),rpor_grasp92(ii1)
            end if
         end do
         !
         300 format(/"From subroutine radgrd_grasp92():"                              &
                    /2(" -------- r -------- -------- r' ------- ------- r'/r ------"))
         301 format(1p,6(1x,1d19.12))
      end if
      !
   end subroutine radgrd_grasp92
   !
   !
   subroutine readwrite_grasp92_grid(stream,read_from)             
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads or writes the GRASP-92 grid structure from or to a file on stream.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                  :: stream
      logical, intent(in)                  :: read_from
      !
      if (read_from) then
         read(stream,*)  n_grasp92,rnt_grasp92,h_grasp92,hp_grasp92, accy_grasp92
         ! Set up the coefficients for the numerical procedures
         call setqic_grasp92()
         !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp92()
         !
      else
         write(stream,*) n_grasp92,rnt_grasp92,h_grasp92,hp_grasp92, accy_grasp92
      end if
      !
   end subroutine readwrite_grasp92_grid
   !
   !
   subroutine readwrite_grasp92_orbital(stream,read_from,wave)             
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads or writes a type(grasp92_orbital) data structure from or to a file on stream.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                  :: stream
      logical, intent(in)                  :: read_from
      type(grasp92_orbital), intent(inout) :: wave
      !
      integer       :: i, j
      !
      if (read_from) then
         read(stream,*) wave%number_of_rwf
	 allocate( wave%rwf(1:wave%number_of_rwf) )
      else
         write(stream,*) wave%number_of_rwf
      end if
      !
      do  i = 1,wave%number_of_rwf
         if (read_from) then
	    read(stream,*) wave%rwf(i)%orbital%n,  wave%rwf(i)%orbital%kappa,  wave%rwf(i)%mtp,   wave%rwf(i)%energy, &
                           wave%rwf(i)%gamma,      wave%rwf(i)%pz
	    allocate( wave%rwf(i)%P(1:wave%rwf(i)%mtp), wave%rwf(i)%Q(1:wave%rwf(i)%mtp) )
            do  j = 1,wave%rwf(i)%mtp
               read(stream,*) wave%rwf(i)%P(j), wave%rwf(i)%Q(j)
            end do
	 else
	    write(stream,*) wave%rwf(i)%orbital%n,  wave%rwf(i)%orbital%kappa, wave%rwf(i)%mtp,   wave%rwf(i)%energy, &
                            wave%rwf(i)%gamma,      wave%rwf(i)%pz
            do  j = 1,wave%rwf(i)%mtp
               write(stream,*) wave%rwf(i)%P(j), wave%rwf(i)%Q(j)
            end do
	 end if
      end do
      !
   end subroutine readwrite_grasp92_orbital
   !
   !
   function Rk_bar_acbd_integral_grasp92(nu,rwf_a,rwf_c,rwf_b,rwf_d, Bessel_j,Bessel_n)   result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Evaluates the transverse interaction integral R bar (k; a c | b d ; w). 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: nu
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp), dimension(1:n_grasp92), intent(in) :: Bessel_j, Bessel_n
      !
      real(kind=dp)                      :: value
      !
      integer :: mtp, mtp_ac, mtp_bd
      real(kind=dp),dimension(1:n_grasp92+10) :: ta, tb
      !
      mtp_ac = min( rwf_a%mtp, rwf_c%mtp )
      mtp_bd = min( rwf_b%mtp, rwf_d%mtp )
      !
      ta(1:mtp_ac) = rwf_a%P(1:mtp_ac) * rwf_c%Q(1:mtp_ac) * (one + Bessel_j(1:mtp_ac)) 
      call zf_k_grasp92(nu,ta,tb,mtp_ac)
      mtp   = min( mtp_ac, mtp_bd )
      ta(1) = zero
      ta(2:mtp) =  rwf_b%P(2:mtp)*rwf_d%Q(2:mtp) * (one + Bessel_n(2:mtp)) * tb(2:mtp) * rpor_grasp92(2:mtp)
      !
      value = quad_grasp92(ta,mtp)
      !      
      grasp92_Rkbar_integral = grasp92_Rkbar_integral + 1
      !
      if (debug_Rk_bar_acbd_grasp92) then
         write(99,*)  " Rk-bar(k=",nu,";"// orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)//","// &
	                                    orbital_name(rwf_c%orbital%n,rwf_c%orbital%kappa)//"|"// &
	                                    orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//","// &
	                                    orbital_name(rwf_d%orbital%n,rwf_d%orbital%kappa)//") = ",value
      end if
      !
   end function Rk_bar_acbd_integral_grasp92
   !
   !
   function rk_integral_grasp92(wave,kk,i,j)               result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is an approximation to:
   !
   !           k 
   !      I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)
   !                   i     j       i     j
   !
   ! where I ( G(r) ; range )  denotes  the  integral  of G(r) over range. The radial functions are taken from the wave 
   ! functions array.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)               :: kk, i, j
      type(grasp92_orbital), intent(in) :: wave
      real(kind=dp)                     :: value
      !
      integer       :: l, mtp
      real(kind=dp), dimension(:), allocatable :: ta
      !
      ! Determine the maximum tabulation point for the integrand
      mtp = min(wave%rwf(i)%mtp, wave%rwf(j)%mtp)
      allocate( ta(1:mtp+10) )
      !
      ! Tabulate the integrand as required for subroutine quad-grasp92; 
      ! the value at the first tabulation point is arbitrary
      ta(1) = zero
      do  l = 2,mtp
         ta(l) = (r_grasp92(l)**kk) * (wave%rwf(i)%P(l) * wave%rwf(j)%P(l) + wave%rwf(i)%Q(l) * wave%rwf(j)%Q(l)) &
                                    * rp_grasp92(l)
      end do
      !
      ! Perform the quadrature
      value = quad_grasp92(ta,mtp)
      deallocate( ta )
      !
   end function rk_integral_grasp92
   !
   !
   function rk_integral_grasp92_ab(wave_a,wave_b,kk,i,j)   result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is an approximation to:
   !
   !           k 
   !      I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)
   !                   i     j       i     j
   !
   ! where I ( G(r) ; range )  denotes  the  integral  of G(r) over range. The radial functions are taken from the 
   ! corresponding wave functions array, i.e. P_i and Q_i from wave_a and P_j and Q_j from wave_b.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)               :: kk, i, j
      type(grasp92_orbital), intent(in) :: wave_a, wave_b
      real(kind=dp)                     :: value
      !
      integer       :: l, mtp
      real(kind=dp), dimension(:), allocatable :: ta
      !
      ! Determine the maximum tabulation point for the integrand
      mtp = min(wave_a%rwf(i)%mtp, wave_b%rwf(j)%mtp)
      allocate( ta(1:mtp+10) )
      !
      ! Tabulate the integrand as required for subroutine quad-grasp92; 
      ! the value at the first tabulation point is arbitrary
      ta(1) = zero
      do  l = 2,mtp
         ta(l) = (r_grasp92(l)**kk) * (wave_a%rwf(i)%P(l) * wave_b%rwf(j)%P(l) + wave_a%rwf(i)%Q(l) * wave_b%rwf(j)%Q(l)) &
                                    * rp_grasp92(l)
      end do
      !
      ! Perform the quadrature
      value = quad_grasp92(ta,mtp)
      deallocate( ta )
      !
   end function rk_integral_grasp92_ab
   !
   !
   function rk_integral_grasp92_cd(rwf_a,rwf_b,kk)         result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is an approximation to:
   !
   !           k 
   !      I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)
   !                   a     b       a     b
   !
   ! where I ( G(r) ; range )  denotes  the  integral  of G(r) over range. The radial functions are taken from the 
   ! corresponding wave functions array, i.e. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: kk
      type(orbital_function), intent(in) :: rwf_a, rwf_b
      real(kind=dp)                      :: value
      !
      integer       :: l, mtp
      real(kind=dp), dimension(:), allocatable :: ta
      !
      ! Determine the maximum tabulation point for the integrand
      mtp = min(rwf_a%mtp, rwf_b%mtp)
      allocate( ta(1:mtp+10) )
      !
      ! Tabulate the integrand as required for subroutine quad-grasp92; 
      ! the value at the first tabulation point is arbitrary
      ta(1) = zero
      do  l = 2,mtp
         ta(l) = (r_grasp92(l)**kk) * (rwf_a%P(l) * rwf_b%P(l) + rwf_a%Q(l) * rwf_b%Q(l)) * rp_grasp92(l)
      end do
      !
      ! Perform the quadrature
      value = quad_grasp92(ta,mtp)
      deallocate( ta )
      !
   end function rk_integral_grasp92_cd
   !
   !
   subroutine schmidt_orthogonalize_grasp92(wave)
   !--------------------------------------------------------------------------------------------------------------------
   ! This routine Schmidt orthogonalises radial wavefunctions in the given wave function array.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(grasp92_orbital), intent(inout) :: wave
      !
      integer       :: i, counter, k, kappa, l, mtp, mtp_min, novl
      real(kind=dp) :: eps, factor, norm, overlap
      integer, dimension(:), allocatable       :: j
      real(kind=dp), dimension(:), allocatable :: ovlap, P, Q
      !
      ! Set tabulated values of the radial wavefunction to zero if they are 
      ! less than eps
      eps = 0.1_dp * accy_grasp92
      !
      ! Determine the number of interesting overlaps
      novl = 0
      do  k = 1,wave%number_of_rwf-1
         kappa = wave%rwf(k)%orbital%kappa
         do  l = k+1,wave%number_of_rwf
            if (kappa == wave%rwf(l)%orbital%kappa) novl = novl + 1
         end do
      end do
      !
      if (novl == 0) return
      allocate( j(1:novl), ovlap(1:novl) )
      !
    1 do  l = 2,wave%number_of_rwf 
         kappa   = wave%rwf(l)%orbital%kappa
         counter = 0
    2    do  k = 1,l-1
            if (kappa == wave%rwf(k)%orbital%kappa) then
               !
               ! Compute overlap
               overlap = rk_integral_grasp92(wave,0,l,k)
               !
               ! Schmidt orthogonalise
               counter        = counter + 1
               j(counter)     = k
               ovlap(counter) = overlap
               wave%rwf(l)%pz = wave%rwf(l)%pz - overlap * wave%rwf(k)%pz
               mtp = max(wave%rwf(l)%mtp, wave%rwf(k)%mtp)
               mtp_min = min(wave%rwf(l)%mtp, wave%rwf(k)%mtp)
               allocate( P(1:mtp), Q(1:mtp) )
               do  i = 1,mtp_min
                  P(i) = wave%rwf(l)%P(i) - overlap*wave%rwf(k)%P(i)
                  Q(i) = wave%rwf(l)%Q(i) - overlap*wave%rwf(k)%Q(i)
               end do
               if (wave%rwf(l)%mtp < wave%rwf(k)%mtp) then
                  P(mtp_min+1:mtp) = - overlap*wave%rwf(k)%P(mtp_min+1:mtp)
                  Q(mtp_min+1:mtp) = - overlap*wave%rwf(k)%Q(mtp_min+1:mtp)
               else if (wave%rwf(l)%mtp > wave%rwf(k)%mtp) then
                  P(mtp_min+1:mtp) = wave%rwf(l)%P(mtp_min+1:mtp)
                  Q(mtp_min+1:mtp) = wave%rwf(l)%Q(mtp_min+1:mtp)
               end if
               !
               ! Reallocate storage and normalise
               wave%rwf(l)%mtp = mtp
               deallocate( wave%rwf(l)%P, wave%rwf(l)%Q )
               allocate(   wave%rwf(l)%P(1:mtp), wave%rwf(l)%Q(1:mtp) )
               wave%rwf(l)%P(:) = P(:);   wave%rwf(l)%Q(:) = Q(:)
               deallocate( P, Q )
               !
               norm    = rk_integral_grasp92(wave,0,l,l)
               !! print *, "schmidt_orthogonalize_grasp92(): norm = ",norm
               factor  = one / sqrt(norm)
               wave%rwf(l)%pz = factor * wave%rwf(l)%pz
               do  i = 2,mtp
                  wave%rwf(l)%P(i) = factor * wave%rwf(l)%P(i)
                  wave%rwf(l)%Q(i) = factor * wave%rwf(l)%Q(i)
               end do
               !
               ! Find new wave%rwf(l)%mtp
               mtp = mtp + 1
    3          mtp = mtp - 1
               if (abs(wave%rwf(l)%P(mtp)) < eps) then
                  wave%rwf(l)%P(mtp) = zero
                  wave%rwf(l)%Q(mtp) = zero
                  goto 3
               else
                  wave%rwf(l)%mtp = mtp
               end if
            end if
         end do
         !
         ! Print overlap information
         if (debug_schmidt_grasp92   .and. counter > 0) then
            write (99,4) (ovlap(i),wave%rwf(l)%orbital%n, orbital_symmetry(wave%rwf(l)%orbital%kappa),         &
                         wave%rwf(j(i))%orbital%n, orbital_symmetry(wave%rwf(j(i))%orbital%kappa),i = 1,counter)
          4 format(1p,5(2x,1d10.3," = <",1i2,1a2,"|",1i2,1a2,">"))
         end if
      end do
      !
   end subroutine schmidt_orthogonalize_grasp92
   !
   !
   function selfenergy_ratio_grasp92(rwf)                  result(ratio)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the QED estimate for the orbital rwf.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(in)    :: rwf
      real(kind=dp)                         :: ratio, R0_qed
      !
      integer                               :: i, mtp
      real(kind=dp)                         :: value, value_hyd
      real(kind=dp), dimension(1:n_grasp92) :: ta
      type(orbital_function)                :: rwf_hydrogen
      !
      R0_qed = 2.0e-5_dp * nuclear_charge
      !
      ! Determine mtp
      do  i = 1,n_grasp92
         if (r_grasp92(i) > R0_qed) then
            mtp = i
            goto 1
         end if
      end do
      stop "selfenergy_shell_grasp92(): program stop A." 
      !
    1 ta(1) = zero
      ta(2:mtp) = (rwf%P(2:mtp)**2 + rwf%Q(2:mtp)**2) * rp_grasp92(2:mtp)
      value = quad_grasp92(ta,mtp)
      !
      rwf_hydrogen%orbital = nkappa(rwf%orbital%n,rwf%orbital%kappa)
      allocate( rwf_hydrogen%P(1:n_grasp92), rwf_hydrogen%Q(1:n_grasp92) )
      call hydrogen_point_grasp92(rwf_hydrogen,nuclear_charge)
      !
      ta(1) = zero
      ta(2:mtp) = (rwf_hydrogen%P(2:mtp)**2 + rwf_hydrogen%Q(2:mtp)**2) * rp_grasp92(2:mtp) 
      value_hyd = quad_grasp92(ta,mtp)
      !
      ratio = value / value_hyd
      !
      deallocate( rwf_hydrogen%P, rwf_hydrogen%Q )
      !
   end function selfenergy_ratio_grasp92
   !
   !
   subroutine setqic_grasp92
   !--------------------------------------------------------------------------------------------------------------------
   ! This  subroutine sets up the coefficients for SUBROUTINEs DPBDT, QUAD, RINTI, START, YZK, ZKF. No subroutines are 
   ! called. It has been taken from GRASP92 (Parpia, Froese Fischer & Grant, 1996 CPC  and adapted to Fortran 90.
   !--------------------------------------------------------------------------------------------------------------------
      !
      logical :: first = .true.
      real(kind=dp), save                     :: b13den, denom, c5den, c6den
      real(kind=dp), dimension(13,13), save   :: b13
      real(kind=dp), dimension(6), save       :: cg
      real(kind=dp), dimension(1:5,2:5), save :: c5num
      real(kind=dp), dimension(1:6,2:6), save :: c6num
      !
      integer :: i, j
      real(kind=dp) :: factor
      !
      if (first) then
      !
      ! Thirteen-point lagrange interpolation coefficients for first
      ! derivative 
      b13( 1,1:13) = (/             -1486442880.0d00,      5748019200.0d00, -15807052800.0d00,  &
                 35126784000.0d00, -59276448000.0d00,     75873853440.0d00, -73766246400.0d00,  &
                 54195609600.0d00, -29638224000.0d00,     11708928000.0d00,  -3161410560.0d00,  &
                   522547200.0d00,    -39916800.0d00  /)
      b13( 2,1:13) = (/               -39916800.0d00,      -967524480.0d00,   2634508800.0d00,  &
                 -4390848000.0d00,   6586272000.0d00,     -7903526400.0d00,   7376624640.0d00,  &
                 -5269017600.0d00,   2822688000.0d00,     -1097712000.0d00,    292723200.0d00,  &
                   -47900160.0d00,      3628800.0d00  /)
      b13( 3,1:13) = (/                 3628800.0d00,       -87091200.0d00,   -684478080.0d00,  &
                  1596672000.0d00,  -1796256000.0d00,      1916006400.0d00,  -1676505600.0d00,  &
                  1149603840.0d00,   -598752000.0d00,       228096000.0d00,    -59875200.0d00,  &
                     9676800.0d00,      -725760.0d00  /)
      b13( 4,1:13) = (/                 -725760.0d00,        13063680.0d00,   -143700480.0d00,  &
                  -476910720.0d00,   1077753600.0d00,      -862202880.0d00,    670602240.0d00,  &
                  -431101440.0d00,    215550720.0d00,       -79833600.0d00,     20528640.0d00,  &
                    -3265920.0d00,       241920.0d00  /)
      b13( 5,1:13) = (/                  241920.0d00,        -3870720.0d00,     31933440.0d00,  &
                  -212889600.0d00,   -303937920.0d00,       766402560.0d00,   -447068160.0d00,  &
                   255467520.0d00,   -119750400.0d00,        42577920.0d00,    -10644480.0d00,  &
                     1658880.0d00,      -120960.0d00  /)
      b13( 6,1:13) = (/                 -120960.0d00,         1814400.0d00,    -13305600.0d00,  &
                    66528000.0d00,   -299376000.0d00,      -148262400.0d00,    558835200.0d00,  &
                  -239500800.0d00,     99792000.0d00,       -33264000.0d00,      7983360.0d00,  &
                    -1209600.0d00,        86400.0d00  /)
      b13( 7,1:13) = (/                   86400.0d00,        -1244160.0d00,      8553600.0d00,  &
                   -38016000.0d00,    128304000.0d00,      -410572800.0d00,            0.0d00,  &
                   410572800.0d00,   -128304000.0d00,        38016000.0d00,     -8553600.0d00,  &
                     1244160.0d00,       -86400.0d00  /)
      b13( 8,1:13) = (/                  -86400.0d00,         1209600.0d00,     -7983360.0d00,  &
                    33264000.0d00,    -99792000.0d00,       239500800.0d00,   -558835200.0d00,  &
                   148262400.0d00,    299376000.0d00,       -66528000.0d00,     13305600.0d00,  &
                    -1814400.0d00,       120960.0d00  /)
      b13( 9,1:13) = (/                  120960.0d00,        -1658880.0d00,     10644480.0d00,  &
                   -42577920.0d00,    119750400.0d00,      -255467520.0d00,    447068160.0d00,  &
                  -766402560.0d00,    303937920.0d00,       212889600.0d00,    -31933440.0d00,  &
                     3870720.0d00,      -241920.0d00  /)
      b13(10,1:13) = (/                 -241920.0d00,         3265920.0d00,    -20528640.0d00,  &
                    79833600.0d00,   -215550720.0d00,       431101440.0d00,   -670602240.0d00,  &
                   862202880.0d00,  -1077753600.0d00,       476910720.0d00,    143700480.0d00,  &
                   -13063680.0d00,       725760.0d00  /)
      b13(11,1:13) = (/                  725760.0d00,        -9676800.0d00,     59875200.0d00,  &
                  -228096000.0d00,    598752000.0d00,     -1149603840.0d00,   1676505600.0d00,  &
                 -1916006400.0d00,   1796256000.0d00,     -1596672000.0d00,    684478080.0d00,  &
                    87091200.0d00,     -3628800.0d00  /)
      b13(12,1:13) = (/                -3628800.0d00,        47900160.0d00,   -292723200.0d00,  &
                  1097712000.0d00,  -2822688000.0d00,      5269017600.0d00,  -7376624640.0d00,  &
                  7903526400.0d00,  -6586272000.0d00,      4390848000.0d00,  -2634508800.0d00,  &
                   967524480.0d00,     39916800.0d00  /)
      b13(13,1:13) = (/                39916800.0d00,      -522547200.0d00,   3161410560.0d00,  &
                -11708928000.0d00,  29638224000.0d00,    -54195609600.0d00,  73766246400.0d00,  &
                -75873853440.0d00,  59276448000.0d00,    -35126784000.0d00,  15807052800.0d00,  &
                 -5748019200.0d00,   1486442880.0d00  /)
      !
      b13den = 479001600.0_dp
      !
      ! Coefficients for Sienkiewicz-Baylis formula
      cg(1:6) = (/ 1771.0d00,   9235.0d00,   5890.0d00,   4610.0d00,     35.0d00,     59.0d00 /)
      !
      denom = 5760.0d00
      !
      ! Five-point Newton-Cotes coefficients for closed integration. Expressed as rational numbers
      c5num(1:5,2) = (/ 251.0d00, 646.0d00,   -264.0d00, 106.0d00, -19.0d00 /)
      c5num(1:5,3) = (/ 232.0d00, 992.0d00,    192.0d00,  32.0d00,  -8.0d00 /)
      c5num(1:5,4) = (/ 243.0d00, 918.0d00,    648.0d00, 378.0d00, -27.0d00 /)
      c5num(1:5,5) = (/ 224.0d00,1024.0d00,    384.0d00,1024.0d00, 224.0d00 /)
      !
      c5den = 720.0d00
      !
      ! Six-point newton-cotes coefficients for closed integration.
      ! Expressed as rational numbers
      c6num(1:6,2) = (/ 475.0d00,1427.0d00,    -798.0d00, 482.0d00,   -173.0d00,  27.0d00  /)
      c6num(1:6,3) = (/ 448.0d00,2064.0d00,     224.0d00, 224.0d00,    -96.0d00,  16.0d00  /)
      c6num(1:6,4) = (/ 459.0d00,1971.0d00,    1026.0d00,1026.0d00,   -189.0d00,  27.0d00  /)
      c6num(1:6,5) = (/ 448.0d00,2048.0d00,     768.0d00,2048.0d00,    448.0d00,   0.0d00  /)
      c6num(1:6,6) = (/ 475.0d00,1875.0d00,    1250.0d00,1250.0d00,   1875.0d00, 475.0d00  /)
      !
      c6den = 1440.0d00
      !
      end if
      !
      ! Lagrange interpolation coefficients, do this initialization once per run only
      if (first) then
      !
      ! Thirteen-point coefficients for dpbdt
         factor = one/b13den
         do  j = 1,13
            do  i = 1,13
               a13_grasp92(i,j) = b13(i,j)*factor
            end do
         end do
         !
         first = .false.
      end if
      !
      !
      ! Sienkiewicz-Baylis coefficients for sbstep
      c_grasp92(1) = cg(1)/denom
      factor       = h_grasp92/denom
      do  i = 2,6
         c_grasp92(i) = cg(i)*factor
      end do
      !
      ! Newton-Cotes coefficients for yzk and quad
      factor = h_grasp92/c5den
      do  j = 2,4
         do  i = 2,5
            cnc5c_grasp92(i,j) = factor*c5num(i,j)
         end do
      end do
      !
      c1_grasp92 = factor*c5num(1,5)
      c2_grasp92 = factor*c5num(2,5)
      c3_grasp92 = factor*c5num(3,5)
      c4_grasp92 = c1_grasp92 + c1_grasp92
      !
      ! Newton-cotes coefficients for start
      factor = h_grasp92/c6den
      do  j = 2,6
         do  i = 1,6
            cnc6c_grasp92(i,j) = factor*c6num(i,j)
         end do
      end do
      !
   end subroutine setqic_grasp92
   !
   !
   function Sk_acbd_integral_grasp92(nu,rwf_a,rwf_c,rwf_b,rwf_d,omega, Bessel_j_m1,Bessel_n_m1,Bessel_j_p1,Bessel_n_p1)  &
                                                                                                              result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Evaluates transverse interaction integrals
   !                       (k) 
   !                      S   (a,c;b,d;w). 
   !
   ! See I P Grant and B J McKenzie, J Phys B: At Mol Phys, 13 (1980) 2671-2681.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: nu
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp), intent(in)          :: omega
      real(kind=dp), dimension(1:n_grasp92), intent(in) :: Bessel_j_m1,Bessel_n_m1, Bessel_j_p1,Bessel_n_p1
      !
      real(kind=dp)                                     :: value, wk
      !
      integer                                           :: mtp, mtp_ac, mtp_bd
      real(kind=dp),dimension(1:n_grasp92+10)           :: ta, tb, tc
      !
      value = zero
      ta(:) = zero;    tb(:) = zero;   tc(:) = zero
      !
      mtp_ac = min( rwf_a%mtp, rwf_c%mtp )
      mtp_bd = min( rwf_b%mtp, rwf_d%mtp )
      wk = nu + nu + one
      !          (k-1)
      ! Compute Z     (rho  ; s)
      !                   ac
      ta(1:mtp_ac) = rwf_a%P(1:mtp_ac) * rwf_c%Q(1:mtp_ac)
      call zf_k_grasp92(nu-1,ta,tc,mtp_ac)
      !
      if (abs(omega) < eps10) then
         call zf_k_grasp92(nu+1,ta,tb,mtp_ac)
         mtp = min( mtp_ac, mtp_bd )
         !
         ta(1) = zero
         ta(2:mtp) =  rwf_b%P(2:mtp) * rwf_d%Q(2:mtp) * rpor_grasp92(2:mtp) * (tb(2:mtp) - tc(2:mtp))
         !
         value = half * wk * quad_grasp92(ta,mtp)
      else
         ta(1:mtp_ac) = -ta(1:mtp_ac) *  Bessel_j_m1(1:mtp_ac)
         call zf_k_grasp92(nu-1,ta,tb,mtp_ac)
         mtp = min( mtp_ac, mtp_bd )
         !
         ta(1) = zero
         ta(2:mtp) =  ((one + Bessel_n_p1(2:mtp)) * tb(2:mtp) - tc(2:mtp) * Bessel_n_p1(2:mtp))   &
                    * rwf_b%P(2:mtp) * rwf_d%Q(2:mtp) / (r_grasp92(2:mtp)**2) * rpor_grasp92(2:mtp)
         !
         value = ((wk/omega)**2) * quad_grasp92(ta,mtp)
         return
         !
         ta(1:mtp_bd) =  rwf_b%P(1:mtp_bd) * rwf_d%Q(1:mtp_bd) * (one + Bessel_j_p1(1:mtp_bd))
         call zf_k_grasp92(nu+1,ta,tb,mtp_bd)
         mtp = min( mtp_bd, mtp_ac )
         !
         ta(1) = zero
         ta(2:mtp) =  rwf_a%P(2:mtp) * rwf_c%Q(2:mtp) &
                      *  (one + Bessel_n_m1(2:mtp)) * tb(2:mtp) - tb(2:mtp) * r_grasp92(2:mtp) * rp_grasp92(2:mtp)
         value = value - quad_grasp92(ta,mtp)*omega*omega / ((nu+nu+three)*(nu+nu-one))
      end if
      !      
      grasp92_Sk_integral = grasp92_Sk_integral + 1
      !
      if (debug_Sk_acbd_grasp92) then
         write(99,*)  " Sk(k=",nu,";"// orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)//","// &
	                                orbital_name(rwf_c%orbital%n,rwf_c%orbital%kappa)//"|"// &
	                                orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//","// &
	                                orbital_name(rwf_d%orbital%n,rwf_d%orbital%kappa)//") = ",value
      end if
      !
   end function Sk_acbd_integral_grasp92
   !
   !
   function slater_integral_grasp92(nu,rwf_a,rwf_b,rwf_c,rwf_d,is_bound)         result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is the Slater integral
   !
   !                            k
   !                           R (abcd)
   !
   ! as defined in GRASP92. The radial orbital functions rwf_a, ..., however, may belong to different orbitals sets so 
   ! that 'relaxation effects' could easily be included.
   !
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: nu
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: value, xx
      logical, optional                  :: is_bound
      !
      integer                            :: mtp
      real(kind=dp), dimension(1:n_grasp92+10) :: ta, tb
      !
      if (is_bound  .or.   &
         (is_bound_orbital(rwf_b)  .and.   is_bound_orbital(rwf_d)) ) then
         call yz_k_grasp92(nu,rwf_b,rwf_d,tb,mtp)
	 mtp = min(mtp, rwf_a%mtp, rwf_c%mtp)     
         ta(1) = zero
         ta(2:mtp) = (rwf_a%P(2:mtp)*rwf_c%P(2:mtp) + rwf_a%Q(2:mtp)*rwf_c%Q(2:mtp))  *  rpor_grasp92(2:mtp) * tb(2:mtp)
         !!x xx = slater_integral_explicit(nu,rwf_a,rwf_b,rwf_c,rwf_d)
      else if (is_bound_orbital(rwf_a)  .and.   is_bound_orbital(rwf_c)) then
         call yz_k_grasp92(nu,rwf_a,rwf_c,tb,mtp)
	 mtp = min(mtp, rwf_b%mtp, rwf_d%mtp)     
         ta(1) = zero
         ta(2:mtp) = (rwf_b%P(2:mtp)*rwf_d%P(2:mtp) + rwf_b%Q(2:mtp)*rwf_d%Q(2:mtp))  *  rpor_grasp92(2:mtp) * tb(2:mtp)
      else if (rabs_use_stop) then 
         !! stop "slater_integral_grasp92(): program stop A."
         !! print *, "slater_integral_grasp92(): program stop A."
         !! value = slater_integral_explicit(nu,rwf_a,rwf_b,rwf_c,rwf_d)
         !! return
         call yz_k_grasp92(nu,rwf_a,rwf_c,tb,mtp)
	 mtp = min(mtp, rwf_b%mtp, rwf_d%mtp)     
         ta(1) = zero
         ta(2:mtp) = (rwf_b%P(2:mtp)*rwf_d%P(2:mtp) + rwf_b%Q(2:mtp)*rwf_d%Q(2:mtp))  *  rpor_grasp92(2:mtp) * tb(2:mtp)
      end if
      !
      ! Perform the quadrature
      value = quad_grasp92(ta,mtp)
      !
      grasp92_Rk_slater = grasp92_Rk_slater + 1
      !
      ! Print the Slater integral out, if requested
      if (.true.) then
         write(62,1) nu, orbital_name(rwf_c%orbital%n,rwf_c%orbital%kappa), &
                         orbital_name(rwf_d%orbital%n,rwf_d%orbital%kappa), & 
	                 orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa), & 
	                 orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa), &
			 value*2.3_dp, rwf_c%energy*27.21_dp, rwf_d%energy*27.21_dp, &
                                       rwf_a%energy*27.21_dp, rwf_b%energy*27.21_dp
	 !
         write(63,2) nu, rwf_c%orbital%n,rwf_c%orbital%kappa, rwf_d%orbital%n,rwf_d%orbital%kappa, & 
	                 rwf_a%orbital%n,rwf_a%orbital%kappa, rwf_b%orbital%n,rwf_b%orbital%kappa, &
			 value*2.3_dp,  rwf_b%energy*27.21_dp
	 !
       1 format("*** R_Slater(", i1, ";", a4,1x, a4,1x, a4,1x, a4,1x, ") = ",es14.6, 10x, 4(es14.6,1x) )             
       2 format(i3, i5,i3,  i5,i3,  i5,i3,  i5,i3, 3x, es14.6, 10x, es14.6 ) 
      end if
      !
   end function slater_integral_grasp92
   !
   !
   function slater_integral_explicit(nu,rwf_a,rwf_b,rwf_c,rwf_d)   result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is the Slater integral
   !
   !                            k
   !                           R (abcd)
   !
   ! as defined in GRASP92. In contrast to the slater_integral_grasp92(), here an explicit 2-dim integration is carried 
   ! out. This is much more time consuming but allows for cases where the radial functions do not allow the solution of 
   ! an ODE.
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: nu
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: value
      !
      integer       :: i, j, mac, mbd
      real(kind=dp) :: rr, rac, rbd
      real(kind=dp), dimension(1:n_grasp92+10) :: rho_ac, rho_bd, ta, tb, tc
      !
      rho_ac(:) = zero;  mac = min(rwf_a%mtp, rwf_c%mtp) - 2  !!! - 121
      rho_ac(2:mac) = rwf_a%P(2:mac)*rwf_c%P(2:mac) + rwf_a%Q(2:mac)*rwf_c%Q(2:mac)
      !
      rho_bd(:) = zero;  mbd = min(rwf_b%mtp, rwf_d%mtp) - 2  !!! - 121
      rho_bd(2:mbd) = rwf_b%P(2:mbd)*rwf_d%P(2:mbd) + rwf_b%Q(2:mbd)*rwf_d%Q(2:mbd)
      !
      ta(:) = zero
      do  i = 2,min(mac,mbd)
         rr  = r_grasp92(i)
         rac = rho_ac(i)/(rr**(nu+1));   rbd = rho_bd(i)/(rr**(nu+1))
         !
         tb(2:i)  = rac * (r_grasp92(2:i)**nu) * rho_bd(2:i) * rp_grasp92(2:i)
         tb(i)    = rac * (r_grasp92(i)**nu) * rho_bd(i) * rp_grasp92(i) * half
         tb(i+1:) = zero
         tc(2:i)  = rbd * (r_grasp92(2:i)**nu) * rho_ac(2:i) * rp_grasp92(2:i)
         tc(i)    = rbd * (r_grasp92(i)**nu) * rho_ac(i) * rp_grasp92(i) * half
         tc(i+1:) = zero
         !
         ta(i) = quad_grasp92(tb,i) + quad_grasp92(tc,i)
      end do
      !
      ta(1:mac) = rp_grasp92(1:mac) * ta(1:mac)  
      value = quad_grasp92(ta,mac)
      !
   end function slater_integral_explicit
   !
   !
   function V_ab_grasp92(rwf_a,rwf_b)                     result(V_ab)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the one-particle integral V (ab) for two given radial orbital functions. The analytical expression for 
   ! this quantity is given as eq (3.23) in  F A  Parpia, M. Tong and C F Fischer.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(in)     :: rwf_a, rwf_b
      real(kind=dp)                          :: V_ab
      !
      integer       :: mtp
      real(kind=dp) :: factor1, factor2
      real(kind=dp), dimension(n_grasp92)    :: P_prime_b, Q_prime_b
      real(kind=dp), dimension(n_grasp92+10) :: ta
      !
      ! Kinetic energy contribution; first, piece involving derivatives
      mtp = min( rwf_a%mtp, rwf_b%mtp )
      call dpbdt_grasp92(rwf_b,P_prime_b,Q_prime_b)
      ta(1)     = zero
      ta(2:mtp) = rwf_a%P(2:mtp) * P_prime_b(2:mtp) + rwf_a%Q(2:mtp) * Q_prime_b(2:mtp)
      V_ab = quad_grasp92(ta,mtp) / h_grasp92
      !
      ! Pieces not involving derivatives
      factor1 = half * ( rwf_a%orbital%kappa *( rwf_a%orbital%kappa+one) - rwf_b%orbital%kappa *( rwf_b%orbital%kappa+one))  
      factor2 = half * (-rwf_a%orbital%kappa *(-rwf_a%orbital%kappa+one) + rwf_b%orbital%kappa *(-rwf_b%orbital%kappa+one))
      ta(1)     = zero
      ta(2:mtp) = rpor_grasp92(2:mtp) * (factor1 * rwf_a%P(2:mtp) * rwf_b%P(2:mtp)  +    &
                                         factor2 * rwf_a%Q(2:mtp) * rwf_b%Q(2:mtp)) 
      V_ab = V_ab - quad_grasp92(ta,mtp)
      !
      if (debug_V_ab_grasp92) then
         write(99,*) " V_ab("//orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)//       &
	             ","//orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//") = ",V_ab
      end if
      !
   end function V_ab_grasp92
   !
   !
   function vpintf_grasp92(rwf_a,rwf_b)                    result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Computes nuclear vacuum polarization integrals.
   !
   ! Calls: quad_grasp92().
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(orbital_function), intent(in) :: rwf_a, rwf_b
      real(kind=dp)                      :: value
      !
      integer :: mtp
      real(kind=dp), dimension(1:n_grasp92+10) :: ta
      !
      mtp   = min(rwf_a%mtp, rwf_b%mtp)
      ta(1) = zero
      ta(2:mtp) = (rwf_a%P(2:mtp)*rwf_b%P(2:mtp)  +  rwf_a%Q(2:mtp)*rwf_b%Q(2:mtp)) * zdist_grasp92(2:mtp)
      !
      value = quad_grasp92(ta,mtp)
      !
      if (debug_vpintf_grasp92) then
        ! write(99,*)                                                       &
	!  " VP_ab("//orbital_name(rwf_a%orbital%n,rwf_a%orbital%kappa)//   &
	!  ","//orbital_name(rwf_b%orbital%n,rwf_b%orbital%kappa)//") = ",value
      end if
      !
   end function vpintf_grasp92
   !
   !
   function W_integral_grasp92(mu,nu,rwf_a,rwf_c,rwf_b,rwf_d)                   result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! The value of this function is one of the W_mu^nu(abcd) integrals as defined in the basis set program. The radial 
   ! orbital functions rwf_a, ... are given however as GRASP92 orbitals. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: mu, nu
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: value
      !
      integer       :: mtp, mtp_ac, mtp_bd
      real(kind=dp) :: omega
      real(kind=dp), dimension(1:n_grasp92+10) :: ta, tb
      real(kind=dp), dimension(1:n_grasp92   ) :: Bessel_j,Bessel_n, Bessel_j_m1,Bessel_n_m1, Bessel_j_p1,Bessel_n_p1
      !
      mtp_ac = min( rwf_a%mtp, rwf_c%mtp )   
      mtp_bd = min( rwf_b%mtp, rwf_d%mtp ) 
      !  
      if (rabs_use_stop   .and.   mu == 1) then
         stop "W_integral_grasp92: program stop A."
      else if (mu == 5) then
         ta(1:mtp_bd) = rwf_b%P(1:mtp_bd) * rwf_d%Q(1:mtp_bd)
         call zf_k_grasp92(nu,ta,tb,mtp_bd)
	 mtp = min(mtp_bd, mtp_ac)
         !!x print *, "mtp, tb(2:5) = ",mtp, tb(2:5)     
         ta(1) = zero
         ta(2:mtp) =  rwf_a%P(2:mtp) * rwf_c%Q(2:mtp) * tb(2:mtp) *  rpor_grasp92(2:mtp)
         !
         value = quad_grasp92(ta,mtp)
      else if (mu == 6) then
         !
         ! Calculate all parts independently
         omega = abs((rwf_a%energy - rwf_c%energy)/c)
         call Bessel_grasp92(nu,omega,bessel_j,bessel_n) 
         value = Rk_bar_acbd_integral_grasp92(nu,rwf_a,rwf_c,rwf_b,rwf_d,Bessel_j,Bessel_n)
         value = value + Rk_bar_acbd_integral_grasp92(nu,rwf_b,rwf_d,rwf_a,rwf_c,Bessel_j,Bessel_n)
         !
         omega = abs((rwf_b%energy - rwf_d%energy)/c)
         call Bessel_grasp92(nu,omega,bessel_j,bessel_n)
         value = value + Rk_bar_acbd_integral_grasp92(nu,rwf_a,rwf_c,rwf_b,rwf_d,Bessel_j,Bessel_n)
         value = value + Rk_bar_acbd_integral_grasp92(nu,rwf_b,rwf_d,rwf_a,rwf_c,Bessel_j,Bessel_n)
         value = half * value
      else if (mu == 7) then
         !
         ! Calculate both parts independently
         omega = abs((rwf_a%energy - rwf_c%energy)/c)
         call Bessel_grasp92(nu-1,omega,bessel_j_m1,bessel_n_m1)
         call Bessel_grasp92(nu+1,omega,bessel_j_p1,bessel_n_p1)
         value = half * Sk_acbd_integral_grasp92(nu,rwf_a,rwf_c,rwf_b,rwf_d, omega,Bessel_j_m1,Bessel_n_m1,  &
                                                                                   Bessel_j_p1,Bessel_n_p1)
         !
         omega = abs((rwf_b%energy - rwf_d%energy)/c)
         call Bessel_grasp92(nu-1,omega,bessel_j_m1,bessel_n_m1)
         call Bessel_grasp92(nu+1,omega,bessel_j_p1,bessel_n_p1)
         value = value + half * Sk_acbd_integral_grasp92(nu,rwf_a,rwf_c,rwf_b,rwf_d,&
                                omega,Bessel_j_m1,Bessel_n_m1,Bessel_j_p1,Bessel_n_p1)
      else if (rabs_use_stop) then
         stop "W_integral_grasp92: program stop B."
      end if
      !
   end function W_integral_grasp92
   !
   !
   subroutine yz_k_grasp92(k,rwf_a,rwf_b,yk,mtp)
   !--------------------------------------------------------------------------------------------------------------------
   ! This routine evaluates Hartree Y- and Z-functions:
   !
   !            (K)            (K)           (K)
   !           Y   (I,J;r) =  Z   (I,J;r) + W   (I,J;r)
   !
   ! where
   !
   !  (K)
   ! Z   (I,J;r) =  I ( (s/r)   (P (s)*P (s) + Q (s)*Q (s)) ; 0 - r )
   !                              I     J       I     J 
   !
   ! and 
   !
   !  (K)                    K+1              
   ! W   (I,J;r) =  I ( (r/s)   (P (s)*P (s) + Q (s)*Q (s)) ; r - 
   !                              I     J       I     J    INFINITY )
   !
   ! where  I ( G(r,s) ; RANGE )  denotes the integral of G(r,s) over a range in  s.  The Y-function is tabulated in the 
   ! array  TB, the Z-function in array TA, but this function is not returned.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                       :: k
      type(orbital_function), intent(in)        :: rwf_a, rwf_b
      real(kind=dp), dimension(:), intent(out)  :: yk
      integer, intent(out)                      :: mtp
      !
      integer                                   :: ii, kk
      real(kind=dp)                             :: dif, sum, zklim
      real(kind=dp), dimension(1:n_grasp92+10)  :: rho, temp, wk, zk
      !
      if (k+1 > kmp1_grasp92) then
         print *, "yz_k_grasp92(): Revise the value of parameter kmp1_grasp92 to at least ",k+1,";"
	 print *, " terminating execution ..."
         stop
      end if
      !
      ! Determine maximum tabulation point as location beyond which rho(:)   (see comment below) would be zero
      mtp = min(rwf_a%mtp,rwf_b%mtp)
      !
      ! Compute p (s)*p (s)+q (s)*q (s) and store in rho
      !          i     j     i     j
      !
      rho(2:mtp) = rwf_a%P(2:mtp)*rwf_b%P(2:mtp) + rwf_a%Q(2:mtp)*rwf_b%Q(2:mtp)
      ! 
      !                            k
      ! Fill array temp(:) with r' * r  * RHO; set additional four points to 0
      temp(1) = zero
      if (k == 0) then
         temp(2:mtp) = rp_grasp92(2:mtp) * rho(2:mtp)
      else
         temp(2:mtp) = rprtk_grasp92(2:mtp,k) * rho(2:mtp)
      end if
      temp(mtp+1:mtp+4) = zero
      !
      !                              k
      ! Compute the first values of r * ZK using semi-open Newton-Cotes formula
      zk(1) = zero
      do  ii = 2,4
         sum = zero
         do  kk = 2,5
            sum = sum + cnc5c_grasp92(kk,ii) * temp(kk)
         end do
         zk(ii) = sum
      end do
      !                       k
      ! Compute remainder of r * zk up to mtp+3; use closed Newton-Cotes formula
      do  ii = 5,mtp+3
         zk(ii) = zk(ii-4) + c1_grasp92*(temp(ii-4)+temp(ii  )) + c2_grasp92*(temp(ii-3)+temp(ii-1)) &
                           + c3_grasp92* temp(ii-2)
      end do
      !                                     k   (k)
      ! Determine the asymptotic value of  r * Z   ; apply a correction to Z^(0) in the manner of  C Froese Fischer, 
      ! The Hartree-Fock Method for Atoms, John Wiley & Sons, New York, 1977, p 235.
      if (k == 0) then
         if (rwf_a%mtp == rwf_b%mtp  .and.  rwf_a%energy == rwf_b%energy) then
            zklim = one
         else
            zklim = zero
         end if
         do  kk = mtp+3,mtp,-1
            dif = zk(kk) - zklim
            do  ii = kk,2,-4
               zk(ii) = zk(ii) - dif
            end do
         end do
      else
         zklim = zk(mtp+3)
      end if
      !
      ! Tabulate zk(:) for entire internal grid
      if (k == 0) then
         zk(mtp+4:n_grasp92) = zklim
      else
         zk(2:mtp+3)         = zk(2:mtp+3) * rtki_grasp92(2:mtp+3,k)
         zk(mtp+4:n_grasp92) = zklim * rtki_grasp92(mtp+4:n_grasp92,k)
      end if
      !
      !                   k+1
      ! Start array wk / r      
      do  ii = n_grasp92+4,mtp+1,-1
         wk(ii) = zero
      end do
      !                           k+1
      ! Fill array temp with r'/ r    * rho; set temp(1) = 0 to avoid 0/0
      temp(1) = zero
      if (k == 0) then
         temp(2:mtp) = rpor_grasp92(2:mtp) * rho(2:mtp)
      else
         temp(2:mtp) = rpbrtk_grasp92(2:mtp,k) * rho(2:mtp)
      end if
      !                            k+1
      ! Compute remainder of WK / r   : march in to the origin
      do  ii = mtp,2,-1
         wk(ii) = wk(ii+4) + c1_grasp92*(temp(ii  )+temp(ii+4)) + c2_grasp92*(temp(ii+1)+temp(ii+3))   &
                           + c3_grasp92*(temp(ii+2))
      end do
      wk(1) = zero
      !
      ! Compute wk
      if (k == 0) then
         wk(2:mtp) = wk(2:mtp) * r_grasp92(2:mtp)
      else
         wk(2:mtp) = wk(2:mtp) * rtk_grasp92(2:mtp,k+1)
      end if
      !
      ! Assemble solution
      yk(1) = zero
      yk(2:n_grasp92) = zk(2:n_grasp92) + wk(2:n_grasp92)
      !
   end subroutine yz_k_grasp92
   !
   !
   subroutine zf_k_grasp92(k,ta,zk,mtp)
   !--------------------------------------------------------------------------------------------------------------------
   ! This routine evaluates Hartree Z-functionals:
   !                                                                   
   !           (k)                     k                               
   !          Z   [f(r);r] =  I ( (s/r)   f(s) ; 0 - r )               
   !                                                                   
   ! where  I ( g(r,s) ; range )  denotes the integral of g(r,s) over range  in  s .    The Z-functional is tabulated 
   ! at return in  zk(:). The f-function is assumed to be tabulated in array  ta(:).  
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                      :: k
      real(kind=dp), dimension(:), intent(in)  :: ta
      real(kind=dp), dimension(:), intent(out) :: zk
      integer, intent(in)                      :: mtp
      !
      integer                                  :: ii, kk
      real(kind=dp)                            :: sum, zklim
      real(kind=dp), dimension(1:n_grasp92+10) :: temp
      !
      !                       k
      ! For  k > 0  compute  r   and store in  rttk
!       if (k > 0) then
!          rttk(2:n_grasp92) = r_grasp92(2:n_grasp92)**k
!       end if
!       !
!       ! Compute  rp(s)*f(s)  and store it in  rhop
!       rhop(2:mtp) = rp_grasp92(2:mtp) * ta(2:mtp)
      !
      ! Fill array temp with r**k * rhop
      if (k+1 > kmp1_grasp92) then
         print *, "zf_k_grasp92(): Revise the value of parameter "//&
                    " kmp1_grasp92 to at least ",k+1,";"
         print *, " terminating execution ..."
         stop
      end if

      temp(1) = zero
      if (k == 0) then
         temp(2:mtp) = ta(2:mtp) * rp_grasp92(2:mtp)
      else
         temp(2:mtp) = rprtk_grasp92(2:mtp,k) * ta(2:mtp)
      end if
      temp(mtp+1:mtp+4) = zero
      !
      !                                    k
      !  Compute the first few values of  r  * zk(:)  using semi-open Newton-Cotes formulae
      zk(1) = zero
      do  ii = 2,4
         sum = zero
         do  kk = 2,5
            sum = sum + cnc5c_grasp92(kk,ii) * temp(kk)
         end do
         zk(ii) = sum
      end do
      !                       k
      ! Compute remainder of r  * zk: march out to mtp+3
      do  ii = 5,mtp+3
         zk(ii) = zk(ii-4) + c1_grasp92*(temp(ii-4)+temp(ii  )) + c2_grasp92*(temp(ii-3)+temp(ii-1)) &
                           + c3_grasp92* temp(ii-2)
      end do
      !                                     k   (k)
      ! Determine the asymptotic value of  r * Z
      ! Compute ZK
      zklim = zk(mtp+3)
      if (k == 0) then
         zk(mtp+4:n_grasp92) = zklim
      else !if(k > 0) then
         zk(2:mtp+3) = zk(2:mtp+3) / rtk_grasp92(2:mtp+3,k) ! rttk(2:mtp+3)
         zk(mtp+4:n_grasp92) = zklim / rtk_grasp92(mtp+4:n_grasp92,k) ! rttk(mtp+4:n_grasp92)
      end if
      !
   end subroutine zf_k_grasp92
   !
end module rabs_grasp92
