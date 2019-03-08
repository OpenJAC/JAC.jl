module rabs_functions_math
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module contains a number of mathematical and related functions as they appear frequently in the computation of 
! atomic properties and structures.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_naglib
   implicit none
   private
   !
   private :: arctan_0_2pi            
                 ! Returns the tan^(-1) (x/y) in the range 0 ... 2*pi.
   public  :: beta_func             
                 ! Returns B(a,b) for integer and half-arguments.
   public  :: calculate_determinant             
                 ! Calculates the determinant of a given quadratic matrix.
   public  :: complete_contraction4 
                 ! Calculates the complete contraction of a 4-dimensional with four vectors.
   public  :: coul_cern_kabachnik
                 ! Calculates the non-relativistic coulomb function for a free electron with given energy and orbital 
                 ! momentum.
   public  :: distance
                 ! Calculates the distance between two (3-dim) points, given either in "cartesian" or "spherical" 
                 ! coordinates.       
   public  :: double_factorial        
                 ! Returns the value n!! for given integer n.
   public  :: factorial        
                 ! Returns the value n! for given integer n.
   public  :: factorial_dp        
                 ! Returns the approximate value n! for given integer n as real(kind=dp).
   public  :: gamma_complex            
                 ! Returns Gamma(arg) for complex arguments.
   public  :: gamma_func            
                 ! Returns Gamma(n) for integer arguments.
   public  :: Gammax_over_ax        
                 ! Calculates Gamma(x)/(a^x)
   public  :: get_maximum_minimum
                 ! Returns the value and location of a local maximum or minimum by interpolation
   public  :: hypergeometric_2F1    
                 ! Calculates Gauss' hypergeometric function 2F_1 = F(a,b,c;x) for real x.
   public  :: hypergeometric_2F1_complex_b    
                 ! Calculates Gauss' hypergeometric function 2F_1 = F(a,b,c;x) for complex a, b, c and x.
   public  :: hypergeometric_2F1_complex_g    
                 ! Calculates Gauss' hypergeometric function 2F_1 = F(a,b,c;x) for complex a, b, c and x.
   public  :: hypergeometric_Phi    
                 ! Calculates the degenerate hypergeometric function 1F_1 = Phi(a,b;x) for real x.
   public  :: incomplete_beta_func  
                 ! Returns the incomplete Beta(a,b,x) function.
   public  :: integral_one_dim      
                 ! Calculates a simple estimate of an 1-dim integral in a given interval.
   public  :: interpolation_aitken
                 ! Interpolates a value by using Aitken's algorithm.   
   public  :: interpolation_lagrange
                 ! Interpolates a value by using Lagrange's algorithm.   
   public  :: is_triangle
                 ! Returns .true. if three integers fulfill the triangular condition and .false. otherwise.       
   public  :: legendre_associated
                 ! Returns the (complex) value of the associated Legendre function P_n^m (z) for complex argument.      
   public  :: legendre_polynomial
                 ! Returns the (complex) value of the Legendre polynomial P_n (z) for complex argument. 
   public  :: LU_decomposition_of_matrix       
                 ! Carries out a LU decomposition of a given matrix.
   public  :: pochhammer             
                 ! Calculates the value of a pochhammer symbol
   public  :: set_beta_gamma_arrays 
                 ! Initializes array for beta_func(), gamma_func(), and incomplete_beta_func().
   public  :: set_gauss_zeros       
                 ! Returns the zeros and weights for a Gauss-Legendre integration.
   private :: set_gauss_zeros_aux  
                 ! Auxiliarity routine for set_gauss_zeros().
   public  :: spherical_Bessel_jL   
                 ! Returns the value of the spherical Bessel function j_L(x) for real values of x.
   public  :: spherical_Legendre_Pmn      
                 ! Returns the value of the spherical (associated) Legendre function P^m_n(x) for real x.
   public  :: transform_cart_to_cylinder  
                 ! Transforms cartesian coordinates (x,y,z) into (rho,phi,z).
   public  :: transform_cart_to_spherical 
                 ! Transforms cartesian coordinates (x,y,z) into (r,theta,phi).
   public  :: transform_cylinder_to_cart  
                 ! Transforms cylindrical coordinates (rho,phi,z) into (x,y,z).
   public  :: transform_spherical_to_cart 
                 ! Transforms spherical coordinates (r,theta,phi) into (x,y,z).
   public  :: Whittaker_M           
                 ! Returns the value of the Whittaker function M_a,b(x) for real values of x.
   public  :: Whittaker_W           
                 ! Returns the value of the Whittaker function W_a,b(x) for real values of x.
   !
   !
   integer, parameter, private :: ndbeta = 32, ndmx = 19
   real(kind=dp), pointer, dimension(:),     private :: gamma_store
   real(kind=dp), pointer, dimension(:,:),   private :: beta_store
   real(kind=dp), pointer, dimension(:,:,:), private :: beta_frac_coeff
   !
   integer, dimension(0:19), private, parameter :: dfak = &
      (/       1,            1,            2,               3,       &
               8,           15,           48,             105,       &
             384,          945,         3840,           10395,       &
           46080,       135135,       645120,         2027025,       &
        10321920,     34459425,    185794560,       654729075       /)
   !  3715891200,  13749310575,  81749606400,    316234143225       /)
   !
   integer, private    :: ifail
   !
   !
contains
   !
   !
   function arctan_0_2pi(arg1,arg2)                        result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns 
   !                   -1                                         
   !       value = tan   (arg1/arg2),          0 <= value < 2*pi.
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: arg1, arg2
      real(kind=dp)             :: value
      !
      logical, save :: first = .true., intrin = .true.
      !
      ! Determine whether the FORTRAN intrinsic function ATAN2 always
      ! returns a positive value
      if (first) then
         value = atan2(-one,-one)
         if (value > zero) then
            intrin = .true.
         else
            intrin = .false.
         end if
         first = .false.
      end if
      !
      ! Use the intrinsic function if it passes the above test; otherwise add 2*pi to the negative values returned by 
      ! the intrinsic function
      if (intrin) then
         value = atan2(arg1,arg2)
      else
         if (arg1 >= zero) then
            value = atan2(arg1,arg2)
         else
            value = pi + pi + atan2(arg1,arg2)
         end if
      end if
      !
   end function arctan_0_2pi
   !
   !
   function beta_func(a2,b2)                          result(beta)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Beta(a,b) = Beta(a2/2,b2/2) function for integer and half-integer arguments using the array 
   ! beta_store(:,:) which needs to be initialized before.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: a2, b2
      real(kind=dp)       :: beta
      !
      if (rabs_use_stop           .and.  &
         (a2 > ndbeta  .or.  b2 > ndbeta)) then
         stop "beta_func(): program stop A."
      end if
      !
      beta = beta_store(a2,b2)
      !
   end function beta_func
   !
   !
   subroutine calculate_determinant(matrix,determinant,n)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the determinant of the n x n matrix; ndim is the dimension in the declaration of matrix.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)        :: n
      real(kind=dp), intent(out) :: determinant
      real(kind=dp), dimension(:,:), intent(inout) :: matrix
      !
      integer       :: j, p
      real(kind=dp) :: d, wa
      !! integer       :: ier
      !! real(kind=dp) :: d1, d2, wa
      !! real(kind=dp), dimension(n) :: worka, workb     
      !
      ! Return zero if one 'raw' is identically zero
      determinant = zero
      do  j = 1,n
        wa = zero
	do  p = 1,n
	   wa = max(wa, abs(matrix(j,p)))
	end do
	if (wa < eps20) return
      end do
      !
      call LU_decomposition_of_matrix(matrix,n,d)
      determinant = d
      do  j = 1,n
         determinant = determinant * matrix(j,j)
      end do
      !
      ! LU-decomposition
      !
      ! The following three lines show how to use the subroutine LUDATN from
      ! the IMSL Library (IMSL Inc, 1984, 7500 Bellaire Boulevard, Houston 
      ! TX 77036, USA). This invokation has to be replaced if some other
      ! routine of numerical library is used.
      ! The only relevant OUTPUT parameter from this routine, however, is
      ! determinant = the value of the determinant det(DWD).
      !
      !! d1 = one;   ier = 0
      !! call ludatn(matrix,ndim,n,matrix,ndim,0,d1,d2,worka,workb,wa,ier)
      !! determinant = d1  *  two ** d2
      !! if (abs(determinant) > eps10) print *, "determinant = ",determinant
      !
      ! NAG library function f03aaf
      !
      !! call f03aaf(matrix,ndim,5,determinant,work,ifail)
      !! if (ifail /= 0) then
      !!    print *, "calculate_determinant(), f03aaf: ifail = ",ifail
      !! end if
      !! print *, "determinant, ifail = ",determinant, ifail
      !
      ! NAG library function f03aff
      !
      !! eps = eps10
      !! call f03aff(n,eps,matrix,ndim,d1,id2,work,ifail)
      !! if (ifail /= 0) then
      !!    print *, "calculate_determinant(): ifail = ",ifail
      !! end if
      !! deter = d1  *  two ** id2
      !
   end subroutine calculate_determinant
   !
   !
   function complete_contraction4(wa,va,vc,vb,vd,na,nc,nb,nd)                       result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of a complete contraction of the 4-dimensional array wa with the vectors va, vc, vb, vd in this 
   ! sequence of indices.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                         :: na, nc, nb, nd
      real(kind=dp), dimension(1:na), intent(in)  :: va
      real(kind=dp), dimension(1:nb), intent(in)  :: vb
      real(kind=dp), dimension(1:nc), intent(in)  :: vc
      real(kind=dp), dimension(1:nd), intent(in)  :: vd
      real(kind=dp), dimension(1:na,1:nc,1:nb,1:nd), intent(in)  :: wa
      real(kind=dp)                                              :: value
      !
      integer :: a, b, c, d
      !
      value = zero
      do  d = 1,nd
         do  b = 1,nb
            do  c = 1,nc
               do  a = 1,na
                  !x print *, "complete_c: a,c,b,d,value = ",a,c,b,d,value
                  value = value + wa(a,c,b,d)*va(a)*vc(c)*vb(b)*vd(d)
               end do
            end do
         end do
      end do
      !
   end function complete_contraction4
   !
   !
   subroutine coul_cern_kabachnik(rho,eta,minl,maxl,fc,fcp,gc,gcp,accur,step)
   !--------------------------------------------------------------------------------------------------------------------
   ! Coulomb wavefunctions calculated at r = rho by the continued-fraction method of steed.   minl,maxl are actual 
   ! l-values. See Barnett et al Comp. Phys. Commun. v.8 p.377 (1974).
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp)               :: K,K1,K2,K3,K4,M1,M2,M3,M4
      real(kind=dp), dimension(:) :: FC,FCP,GC,GCP
      !
      integer       :: minl, maxl, ktr, lmax, lmin1, ktrp, l, lp, i2
      real(kind=dp) :: rho, eta, accur, step, pace, acc, r, xll1, eta2, turn, tf, f, tfp, fp, etar, rho2, pl, pmx, dk,   &
                       del, d, h, p, q, ar, ai, br, bi, wi, dr, di, dq, t, g, gp, w, r3, h2, dp_x, etah, h2ll, s, rh2, rh
      !
      !
      PACE = STEP
      ACC = ACCUR
      IF(PACE.LT.100.0) PACE = 100.0
      IF(ACC.LT.1.0E-15.OR.ACC.GT.1.0E-6) ACC = 1.0E-6
      R = RHO
      KTR = 1
      LMAX = MAXL
      LMIN1 = MINL + 1
      XLL1 = FLOAT(MINL*LMIN1)
      ETA2 = ETA*ETA
      TURN = ETA + SQRT(ETA2 + XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.0E-6) KTR = -1
      KTRP = KTR
      GO TO 2
1     R = TURN
      TF = F
      TFP = FP
      LMAX = MINL
      KTRP = 1
2     ETAR = ETA*R
      RHO2 = R*R
      PL = FLOAT(LMAX + 1)
      PMX = PL + 0.5
      !
      ! *** CONTINUED FRACTION FOR FP(MAXL)/F(MAXL) XL IS F  XLPRIME IS FP **
      !
      FP = ETA/PL + PL/R
      DK = ETAR*2.0
      DEL = 0.0
      D = 0.0
      F = 1.0
      K = (PL*PL - PL + ETAR)*(2.0*PL - 1.0)
      IF(PL*PL+PL+ETAR.NE.0.0) GO TO 3
      R = R + 1.0E-6
      GO TO 2
3     H = (PL*PL + ETA2)*(1.0 - PL*PL)*RHO2
      K = K + DK + PL*PL*6.0
      D = 1.0/(D*H + K)
      DEL = DEL*(D*K - 1.0)
      IF(PL.LT.PMX) DEL = -R*(PL*PL + ETA2)*(PL +1.0)*D/PL
      PL = PL + 1.0
      FP = FP + DEL
      IF(D.LT.0.0) F=-F
      IF(PL.GT.20000.) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP = F*FP
      IF(LMAX.EQ.MINL) GO TO 5
      FC(LMAX+1) = F
      FCP(LMAX+1) = FP
      !
      ! *** DOWNWARD RECURSION TO MINL FOR F AND FP, ARRAYS GC,GCP ARE STORAGE
      !
      L = LMAX
          DO 4 LP = LMIN1,LMAX
          PL = FLOAT(L)
          GC(L+1) = ETA/PL + PL/R
          GCP(L+1) = SQRT(ETA2 + PL*PL)/PL
          FC(L) = (GC(L+1)*FC(L+1) + FCP(L+1))/GCP(L+1)
          FCP(L) = GC(L+1)*FC(L) - GCP(L+1)*FC(L+1)
4         L = L - 1
      F = FC(LMIN1)
      FP = FCP(LMIN1)
5     IF(KTRP.EQ.-1) GO TO 1
      !
      ! *** REPEAT FOR R = TURN IF RHO LT TURN
      ! *** NOW OBTAIN P + 1.Q FOR MINL FROM CONTINUED FRACTION (32)
      ! *** REAL ARITHMETIC TO FACILITATE CONVERSION TO IBM USING REAL*8
      !
      P = 0.0
      Q = R - ETA
      PL = 0.0
      AR = -(ETA2 + XLL1)
      AI = ETA
      BR = 2.0*Q
      BI = 2.0
      WI = 2.0*ETA
      DR = BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP_x = -(AR*DI + AI*DR)
      DQ = (AR*DR - AI*DI)
6     P = P + DP_x
      Q = Q + DQ
      PL = PL +2.0
      AR = AR + PL
      AI = AI + WI
      BI = BI + 2.0
      D = AR*DR - AI*DI + BR
      DI = AI*DR + AR*DI + BI
      T = 1.0/(D*D + DI*DI)
      DR = T*D
      DI = -T*DI
      H = BR*DR - BI*DI - 1.0
      K = BI*DR + BR*DI
      T = DP_x*H - DQ*K
      DQ = DP_x*K + DQ*H
      DP_x = T
      IF(PL.GT.46000.) GO TO 11
      IF(ABS(DP_x)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6
      P = P/R
      Q = Q/R
      !
      ! *** SOLVE FOR FP,G,GP AND NORMALISE F  AT L=MINL
      !
      G = (FP - P*F)/Q
      GP = P*G - Q*F
      W = 1.0/SQRT(FP*G - F*GP)
      G = W*G
      GP = W*GP
      IF(KTR.EQ.1) GO TO 8
      F = TF
      FP = TFP
      LMAX = MAXL
      !
      ! *** RUNGE-KUTTA INTEGRATION OF G(MINL) AND GP(MINL) INWARDS FROM TURN
      ! *** SEE FOX AND MAYERS 1968 PG 202
      !
      IF(RHO.LT.0.2*TURN) PACE = 999.0
      R3 = 1.0/3.0D0
      H = (RHO - TURN)/(PACE + 1.0)
      H2 = 0.5*H
      !!!   I2 = IFIX(PACE + 0.001)
      I2 = INT(PACE + 0.001)
      ETAH = ETA*H
      H2LL = H2*XLL1
      S = (ETAH + H2LL/R)/R - H2
7     RH2 = R + H2
      T = (ETAH + H2LL/RH2)/RH2 - H2
      K1 = H2*GP
      M1 = S*G
      K2 = H2*(GP + M1)
      M2 = T*(G + K1)
      K3 = H*(GP + M2)
      M3 = T*(G + K2)
      M3 = M3 + M3
      K4 = H2*(GP + M3)
      RH = R + H
      S = (ETAH + H2LL/RH)/RH - H2
      M4 = S*(G + K3)
      G = G + (K1 + K2 + K2 + K3 + K4)*R3
      GP = GP + (M1 + M2 + M2 + M3 + M4)*R3
      R = RH
      I2 = I2 - 1
      IF(ABS(GP).GT.1.0E37) GO TO 11
      IF(I2.GE.0) GO TO 7 
      W = 1.0/(FP*G - F*GP)
      !
      ! *** UPWARD RECURSION FROM GC(MINL) AND GCP(MINL),STORED VALUES ARE R,S
      ! *** RENORMALISE FC, FCP FOR EACH L-VALUE
      !
8     GC (LMIN1) = G
      GCP(LMIN1) = GP
      IF (LMAX.EQ.MINL) GO TO 10
            DO 9 L=LMIN1,LMAX 
            T = GC(L+1)
            GC(L+1) = (GC(L)*GC(L+1) - GCP(L))/GCP(L+1)
            GCP(L+1) = GC(L)*GCP(L+1) - GC(L+1)*T
            FC(L+1) = W*FC(L+1)
9           FCP(L+1) = W*FCP(L+1)
      FC(LMIN1) = FC(LMIN1)*W
      FCP(LMIN1) = FCP(LMIN1)*W
      RETURN
10    FC(LMIN1) = W*F
      FCP(LMIN1) = W*FP
      RETURN
11    W = 0.0
      G = 0.0
      GP = 0.0
      GO TO 8
      !
   end subroutine coul_cern_kabachnik
   !
   !
   function distance(coordinates,x1,y1,z1,x2,y2,z2)                                  result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the "distance" between two 3-dimensional points; these points can be given in either "cartesian" or 
   ! "spherical" coordinates.
   !
   ! Calls: transform_spherical_to_cart().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=9), intent(in)  :: coordinates
      real(kind=dp), intent(in)     :: x1, y1, z1, x2, y2, z2
      real(kind=dp)                 :: value, x1w, y1w, z1w, x2w, y2w, z2w
      !
      if (coordinates(1:9) == "cartesian") then
         value = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) )
      else if (coordinates(1:9) == "spherical") then
         call transform_spherical_to_cart(x1,y1,z1,x1w,y1w,z1w)
         call transform_spherical_to_cart(x2,y2,z2,x2w,y2w,z2w)
         value = sqrt( (x2w-x1w)*(x2w-x1w) + (y2w-y1w)*(y2w-y1w) + (z2w-z1w)*(z2w-z1w) )
      else if (rabs_use_stop) then
         stop "distance(): program stop A."
      end if
      !
   end function distance
   !
   !
   function double_factorial(n)                                            result(dfactorial)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value n!! for given positive integer n.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: n
      integer             :: dfactorial
      !
      if (n < 0  .or.  n > 23) then
         print *, "n = ",n
         stop "double_factorial(): program stop A."
      else
         dfactorial = dfak(n)
      end if
      !
   end function double_factorial
   !
   !
   function factorial(n)                                                    result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value n! for given positive integer n.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: n
      integer             :: value, i
      !
      if (rabs_use_stop   .and.   n < 0) then
         print *, "n = ",n
         stop "factorial(): program stop A."
      else if (n== 0) then
         value = 1 
      else
         value = 1
	 do  i = n,2,-1
	    value = value * i
	 end do 
      end if
      !
   end function factorial
   !
   !
   function factorial_dp(n)                                   result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value n! for given positive integer n.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: n
      integer             :: i
      real(kind=dp)       :: value
      !
      if (rabs_use_stop   .and.   n < 0) then
         print *, "n = ",n
         stop "factorial_dp(): program stop A."
      else if (n== 0) then
         value = one 
      else
         value = one
	 do  i = n,2,-1
	    value = value * i
	 end do 
      end if
      !
   end function factorial_dp
   !
   !
   function gamma_complex(arg)                            result(cgamma)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns in cgamma the complex Gamma function of the complex argument arg. Only real(cgamma) is nonzero if aimag(arg) 
   ! is zero. The arctan function required must return angles (in radians) in the range  [0,2*\pi).
   ! This routine was originally written by F A Parpia and has been adapted to the Fortran 90/95 standard.
   !--------------------------------------------------------------------------------------------------------------------
      !
      complex(kind=dp), intent(in)  :: arg
      complex(kind=dp)              :: cgamma
      !
      logical, save       :: first = .true.
      real(kind=dp), save :: hlntpi, twoi
      logical             :: negarg
      !
      ! These are the Bernoulli numbers B02, B04, ..., B14, expressed as rational numbers. From Abramowitz and Stegun, p. 810. 
      real(kind=dp), dimension(7), save ::                                                         &
         fn = (/   1.0_dp,   -1.0_dp,    1.0_dp,   -1.0_dp,    5.0_dp,  -691.0_dp,    7.0_dp  /),  &
         fd = (/   6.0_dp,   30.0_dp,   42.0_dp,   30.0_dp,   66.0_dp,  2730.0_dp,    6.0_dp  /)
      !
      integer       :: i
      real(kind=dp) :: argum, argur, argur2, argui, argui2, clngr, clngi, diff, fac, facneg, obasq, obasqi, obasqr,  &
                       ovlfac, ovlfr, ovlfi, termr, termi, zfacr, zfaci
      !
      ! On the first entry to this routine, set up the constants required for the reflection formula (cf. Abramowitz and 
      ! Stegun 6.1.17) and Stirling's approximation (cf. Abramowitz and Stegun 6.1.40).
      if (first) then
         hlntpi = half * log(pi+pi)
         do  i = 1,7
            fn(i) = fn(i) / fd(i)
            twoi  = dble(i+i)
            fn(i) = fn(i) / (twoi*(twoi-one))
         end do
         first = .false.
      end if
      !
      ! Cases where the argument is real
      if (aimag(arg) == zero) then
         !
         ! Cases where the argument is real and negative
         if (real(arg) <= zero) then
         !
         ! Stop with an error message if the argument is too near a pole
         ! diff <= precis+precis, originally
            diff = abs( nint(real(arg)) - real(arg) )
            if (diff <= eps10*eps10) then
               print *, "Argument (",arg,") is too close to a pole."
               stop     "gamma_complex(): program stop A."
            else
               !
               ! Otherwise use the reflection formula (Abramowitz and Stegun 6.1.17) to ensure that the argument is 
               ! suitable for Stirling's formula
               argum = pi / (-real(arg) * sin(pi*real(arg)))
               if (argum < zero) then
                  argum = -argum
                  clngi = pi
               else
                  clngi = zero
               end if
               facneg = log (argum)
               argur  = -real(arg)
               negarg = .true.
            end if
            ! 
            ! Cases where the argument is real and positive
         else
            clngi  = zero
            argur  = real(arg)
            negarg = .false.
         end if
         !
         ! Use abramowitz and stegun formula 6.1.15 to ensure that the argument in stirling's formula is greater than 10
         ovlfac = one
       2 if (argur < ten) then
            ovlfac = ovlfac * argur
            argur  = argur + one
            goto 2
         end if
         !
         ! Use stirling's formula to compute log (gamma (argum))
         clngr = (argur-half) * log(argur) - argur + hlntpi
         fac   = argur
         obasq = one / (argur*argur)
         do  i = 1,7
            fac   = fac * obasq
            clngr = clngr + fn(i) * fac
         end do
         !
         ! Include the contributions from the recurrence and reflection formulae
         clngr = clngr - log(ovlfac)
         if (negarg) clngr = facneg - clngr
      else
         !
         ! Cases where the argument is complex
         argur  = real(arg);   argui = aimag(arg)
         argui2 = argui * argui
         !
         ! Use the recurrence formula (Abramowitz and Stegun 6.1.15) to ensure that the magnitude of the argument in 
         ! Stirling's formula is greater than 10
         ovlfr = one;   ovlfi = zero
    4    argum = sqrt(argur*argur + argui2)
         if (argum < 10) then
            termr = ovlfr*argur - ovlfi*argui
            termi = ovlfr*argui + ovlfi*argur
            ovlfr = termr
            ovlfi = termi
            argur = argur + one
            goto 4
         end if
         !
         ! Use stirling's formula to compute log (gamma (argum))
         argur2 = argur * argur
         termr  = half * log(argur2 + argui2)
         termi  = arctan_0_2pi(argui,argur)
         clngr  = (argur-half) * termr - argui*termi - argur+hlntpi
         clngi  = (argur-half) * termi + argui*termr - argui
         fac    = (argur2+argui2)**(-2)
         obasqr = (argur2-argui2) * fac
         obasqi = -two * argur * argui * fac
         zfacr  = argur;   zfaci  = argui
         do  i = 1,7
            termr = zfacr*obasqr - zfaci*obasqi
            termi = zfacr*obasqi + zfaci*obasqr
            fac   = fn(i)
            clngr = clngr + termr*fac
            clngi = clngi + termi*fac
            zfacr = termr;   zfaci = termi
         end do
         !
         ! Add in the relevant pieces from the recurrence formula
         clngr = clngr - half * log(ovlfr*ovlfr + ovlfi*ovlfi)
         clngi = clngi - arctan_0_2pi (ovlfi,ovlfr)
      end if
      !
      ! Now exponentiate the complex log gamma function to get the complex gamma function
      if (clngr <= maxexponent(one)   .and.  clngr >= minexponent(one)) then
         fac = exp (clngr)
      else
         print *, "Argument to exponential function",clngr,"out of range."
         stop     "gamma_complex(): program stop A."
      end if
      !
      cgamma = cmplx(fac *cos(clngi), fac *sin(clngi))
      !
   end function gamma_complex
   !
   !
   function gamma_func(a2)                                 result(gamma)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Gamma(a) = Gamma(a2/2) function for integer and half-integer arguments using the array 
   ! gamma_store(:,:) which needs to be initialized before.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: a2
      real(kind=dp)       :: gamma
      !
      if (rabs_use_stop   .and.  a2 > ndbeta) then
         stop "gamma_func(): program stop A."
      end if
      gamma = gamma_store(a2)
      !
   end function gamma_func
   !
   !
   function Gammax_over_ax(x,a)                           result(Gammax)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the value Gamma(x) / a^x by calling the function nag_nag_s14aaf() from the NAG library in order to compute 
   ! Gamma function for non-integer arguments.
   !
   ! Calls: nag_s14aaf() [from NAG library].
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: a, x
      real(kind=dp)             :: Gammax
      !
      if (rabs_use_naglib) then
         Gammax = nag_s14aaf(x,ifail) / ( a ** x )
      else 
         stop "Gammax_over_ax(): program stop A."
      end if
      !
   end function Gammax_over_ax
   !
   !
   subroutine get_maximum_minimum(keyword,x,y,n,xloc,yval)                      
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the location xloc and the corresponding value yval of either a local maximum or minimum by interpolation. 
   ! The 'function' must be contained in y(n) = f(x(n)) and must be either convex or concav according to the keyword. 
   ! It is also assumed that the maximum/minimum is not located 'too close' of one of the boundaries x(1) and x(n).
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=7), intent(in)              :: keyword
      integer, intent(in)                       :: n
      real(kind=dp), dimension(1:n), intent(in) :: x, y
      real(kind=dp), intent(out)                :: xloc, yval
      !
      real(kind=dp) :: x1, x2, x3, dx, y1, y2, y3
      !
      dx = (x(n) - x(1)) / (ten*two)
      x1 = x(1) - dx/two
      !
      do
         x1 = x1 + dx;   x2 = x1 + dx;   x3 = x2 + dx
         call interpolation_aitken(x,y,n,x1,y1)
         call interpolation_aitken(x,y,n,x2,y2)
         call interpolation_aitken(x,y,n,x3,y3)
         !
         select case(keyword)
         case("maximum")
            if (y1 < y2   .and.   y2 > y3) then
               if (abs(x1-x2) < 1.0e-6) then
                  xloc = x2;   yval = y2
                  return
               end if
               dx = -dx / ten;   x1 = x3 - dx
            else if (x1 < x(1)  .or.  x3 < x(1)  .or.  x1 > x(n)  .or.  x3 > x(n)) then
               stop "get_maximum_minimum(): program stop A."
            end if
         case("minimum")
            if (y1 > y2   .and.   y2 < y3) then
               if (abs(x1-x2) < 1.0e-6) then
                  xloc = x2;   yval = y2
                  return
               end if
               dx = -dx / ten;   x1 = x3 - dx
            else if (x1 < x(1)  .or.  x3 < x(1)  .or.  x1 > x(n)  .or.  x3 > x(n)) then
               stop "get_maximum_minimum(): program stop B."
            end if
         case default
            stop "get_maximum_minimum(): program stop C."
         end select
      end do
      !
   end subroutine get_maximum_minimum
   !
   !
   function hypergeometric_2F1(a,b,c,x)                    result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Gauss' hypergeometric function 2F_1 = F(a,b,c;x) for real parameter a, b, and c, and for 
   ! real argument x. The hypergeometric function is calculated by its power expansion as given by I. S. Gradshteyn and 
   ! I. M. Ryzhik, Tables of Integrals, Series, and Products (Academic Press, New York, London a.o. 1980).
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: a, b, c, x
      real(kind=dp)             :: value
      !
      integer                   :: i
      real(kind=dp)             :: wa, wb, wc, wd
      !  
      if (rabs_use_stop   .and.  c == zero) then
         stop "hypergeometric_2F1(): program stop A."
      end if
      !
      wa = a*b/c;   wb = x;   wc = one;  
      !
      value = one + wa*wb
      !
      do  i = 2,400
         wa = wa * (a+i-1) * (b+i-1) / (c+i-1)
         wb = wb * x
         wc = wc * i
         wd = wa*wb/wc
         value = value + wd
	 !
         if (abs(wd) < epsilon(value)) then
            return
         end if
      end do
      !
      print *, "i, value, wd = ",i, value, wd
      stop     "hypergeometric_2F1(): program stop B."
      !
   end function hypergeometric_2F1
   !
   !
   function hypergeometric_2F1_complex_b(a,b,c,z)                                       result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Gauss' hypergeometric function 2F_1 = F(a,b,c;x) for real parameter a, and c, and for 
   ! complex b and argument z. The hypergeometric function is calculated by its power expansion as given by 
   ! I. S. Gradshteyn and I. M. Ryzhik, Tables of Integrals, Series, and Products (Academic Press, New York, 
   ! London a.o. 1980).
   !--------------------------------------------------------------------------------------------------------------------
      !
      complex(kind=dp), intent(in) :: a, b, c, z
      complex(kind=dp)             :: value
      !
      integer                      :: i
      complex(kind=dp)             :: wa, wb, wc, wd
      !  
      if (rabs_use_stop   .and.  c == zero) then
         stop "hypergeometric_2F1_complex_b(): program stop A."
      end if
      !
      wa = a*b/c;   wb = z;   wc = one;  
      !
      value = one + wa*wb
      !
      do  i = 2,400
         wa = wa * (a+i-1) * (b+i-1) / (c+i-1)
         wb = wb * z
         wc = wc * i
         wd = wa*wb/wc
         value = value + wd
	 !
         if (abs(wd) < epsilon(abs(value))) then
            return
         end if
      end do
      !
      print *, "i, value, wd = ",i, value, wd
      stop     "hypergeometric_2F1_complex_b(): program stop B."
      !
   end function hypergeometric_2F1_complex_b
   !
   !
   function hypergeometric_2F1_complex_g(a,b,c,z)                            result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Gauss' hypergeometric function 2F_1 = F(a,b,c;x) for real parameter a, and c, and for 
   ! complex b and argument z. The hypergeometric function is calculated by its power expansion as given by 
   ! I. S. Gradshteyn and I. M. Ryzhik, Tables of Integrals, Series, and Products (Academic Press, New York, 
   ! London a.o. 1980).
   !--------------------------------------------------------------------------------------------------------------------
      !
      complex(kind=dp), intent(in) :: a, b, c, z
      complex(kind=dp)             :: value
      !
      integer                      :: i
      complex(kind=dp)             :: wa, wb, wc, wd, F1, F2
      !  
      if (rabs_use_stop   .and.  c == zero) then
         stop "hypergeometric_2F1_complex_g(): program stop A."
      else if (abs(z) <= one) then
         print *, "Standard way"
         value = hypergeometric_2F1_complex_b(a,b,c,z)
      !! else if (abs(one-one/z) < one  	     .and.  &
      !!	  abs(atan2(aimag(z),real(z))) < Pi  .and.  &
      !!	  abs(atan2(aimag(-z),real(one-z))) < Pi) then
      !!    print *, "3rd extended way"
      !!    value = gamma_complex(c)*gamma_complex(c-a-b)		       &
      !!	  /(gamma_complex(c-a)*gamma_complex(c-b)) *(z**a)	       &
      !!	  * hypergeometric_2F1_complex_b(a,a-c+one,a+b-c+one,one-one/z)&
      !!	  + (one-z)**(c-a-b) * z**(a-c) 		               &
      !!	  * gamma_complex(c)*gamma_complex(a+b-c)		       &
      !!	  /(gamma_complex(a)*gamma_complex(b))  		       &
      !!	  * hypergeometric_2F1_complex_b(c-a,one-a,c-a-b+one,one-one/z)
      else if (abs(one-z) < one  .and.  &
               abs(atan2(aimag(-z),real(one-z))) < Pi) then
         F1 = gamma_complex(c)*gamma_complex(c-a-b) /(gamma_complex(c-a)*gamma_complex(c-b))                &
               * hypergeometric_2F1_complex_b(a,b,a+b-c+one,one-z)
         F2 = (one-z)**(c-a-b) * gamma_complex(c)*gamma_complex(a+b-c)                                      &
               /(gamma_complex(a)*gamma_complex(b)) * hypergeometric_2F1_complex_b(c-a,c-b,c-a-b+one,one-z)
         print *, "2nd extended way; F1, F2 = ",F1,F2
         !
         value = gamma_complex(c)*gamma_complex(c-a-b) /(gamma_complex(c-a)*gamma_complex(c-b))             &
                 * hypergeometric_2F1_complex_b(a,b,a+b-c+one,one-z)                                        &
               + (one-z)**(c-a-b) * gamma_complex(c)*gamma_complex(a+b-c)                                   &
                 /(gamma_complex(a)*gamma_complex(b)) * hypergeometric_2F1_complex_b(c-a,c-b,c-a-b+one,one-z)
      else if (abs(z) > one  .and.  &
               abs(atan2(aimag(-z),real(-z))) < Pi) then
         print *, "1st extended way"
         value = gamma_complex(c)*gamma_complex(b-a) /(gamma_complex(b)*gamma_complex(c-a)) * (-z)**(-a)    &
                 * hypergeometric_2F1_complex_b(a,one-c+a,one-b+a,one/z)                                    &
               + gamma_complex(c)*gamma_complex(a-b) /(gamma_complex(c-b)*gamma_complex(a)) * (-z)**(-b)    &
                 * hypergeometric_2F1_complex_b(b,one-c+b,one-a+b,one/z)
      else
         stop "hypergeometric_2F1_complex_g(): program stop B."
      end if
      !
      !
   end function hypergeometric_2F1_complex_g
   !
   !
   function hypergeometric_Phi(a,b,x)                                                result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the degenerate hypergeometric function 1F_1 = Phi(a,b;x) for real parameters a and b, and for 
   ! real argument x. The degenerate hypergeometric function is calculated by applying its power expansion as given by 
   ! I. S. Gradshteyn and I. M. Ryzhik, Tables of Integrals, Series, and Products (Academic Press, New York, 
   ! London a.o. 1980).
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: a, b, x
      real(kind=dp)             :: value
      !
      integer                   :: i
      real(kind=dp)             :: wa, wb, wc, wd
      !  
      if (rabs_use_stop   .and.   b == zero) then
         stop "hypergeometric_Phi(): program stop A."
      end if
      !
      wa = a/b;   wb = x;   wc = one;  
      !
      value = one + wa*wb
      !
      do  i = 2,400
         wa = wa * (a+i-1) / (b+i-1)
         wb = wb * x / i
         wd = wa * wb
         value = value + wd
         if (abs(wd) < epsilon(value)) then
            return
         end if
      end do
      !
      !! print *, "i, a, b, x, value, wd = ", i, a, b, x, value, wd
      print *, "hypergeometric_Phi(): program stop B."
      stop
      !
   end function hypergeometric_Phi
   !
   !
   function incomplete_beta_func(a2,b2,x)                                result(incbeta)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the incomplete Beta(a,b,x) = Beta(a2/2,b2/2,x) function for integer and half-integer arguments 
   ! a, b by using the array beta_frac_coeff(:,:) which needs to be initialized before.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: a2, b2
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: incbeta
      !
      integer                   :: m
      real(kind=dp)             :: a, b, bt, wn, wx
      !
      if (rabs_use_stop           .and.  &
         (a2 > ndbeta  .or.  b2 > ndbeta)) then
         stop "incomplete_beta_func(): program stop A."
      end if
      !
      a  = a2 / two
      b  = b2 / two
      !
      if (rabs_use_stop   .and.  (x < zero  .or.  x > one)) then
         stop "incomplete_beta_func(): program stop B."
      else if (x == zero  .or.  x == one ) then
         bt = zero
      else
         bt = (x**a) * ( (one - x)**b )
      end if
      !
      if (x < (a + one)/(a + b + two) ) then
         wx = x
         wn = beta_frac_coeff(ndmx,a2,b2) * wx + one
         do  m = (ndmx-1),1,-1
            wn = beta_frac_coeff(m,a2,b2) * wx / wn  +  one
         end do
         incbeta = bt / (a * wn)
      else
         wx = one - x
         wn = beta_frac_coeff(ndmx,b2,a2) * wx + one
         do  m = (ndmx-1),1,-1
            wn = beta_frac_coeff(m,b2,a2) * wx / wn  +  one
         end do
         incbeta = beta_store(a2,b2) - bt / (b * wn)
      endif
      !
   end function incomplete_beta_func
   !
   !
   function integral_one_dim(points,func,delta_x)         result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a simple estimate of the integral value  int func(x) dx  with an equidistant step size delta_x. points is 
   ! the number of grid points and func(i) contains the corresponding values of the function.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer,       intent(in) :: points
      real(kind=dp), intent(in) :: delta_x
      real(kind=dp), dimension(:), intent(in) :: func
      real(kind=dp)             :: value
      !
      integer                   :: i
      !  
      value = zero
      do  i = 2,points
        value = value + func(i) * delta_x 
      end do
      !
   end function integral_one_dim
   !
   !
   subroutine interpolation_aitken(xarr,yarr,narr,xval,yval,accy)
   !--------------------------------------------------------------------------------------------------------------------
   ! This routine returns  yval  as the functions value of  xval  by interpolating a pair of arrays xarr(1:narr), 
   ! yarr(1:narr), that tabulate a  function.  Aitken's algorithm is used. See, for  instance, F B Hildebrand, Introduction 
   ! to Numerical Analysis, 2nd ed., McGraw-Hill, New York, NY, 1974. accy is the  desired accuracy of the estimate: 
   ! a warning message is issued if this is not achieved.  A warning is also issued when the routine is extrapolating. 
   ! This procedures is adapted to Fortran 90/95 from the routine interp() of GRASP92 which originally was written by 
   ! F A Parpia.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer,       intent(in)  :: narr
      real(kind=dp), intent(in)  :: xval
      real(kind=dp), intent(out) :: yval
      real(kind=dp), dimension(narr), intent(in) :: xarr, yarr
      real(kind=dp), intent(in), optional        :: accy
      !
      ! mxord is the maximum order of the interpolation
      integer, parameter                            :: mxord = 11
      logical, dimension(2*mxord+2)                 :: used
      real(kind=dp), dimension(mxord)               :: dx, x, est
      real(kind=dp), dimension((mxord*(mxord+1))/2) :: poly
      !
      logical       :: set
      integer       :: ibest, ilirok, ildiag, ilothr, irow, k, lhi, llo, llr, locnxt, nrsthi, nrstlo
      real(kind=dp) :: debe, debeb, diff, difft
      !
      ! Determine the nearest two XARR entries bounding XVAL
      if (xval < xarr(1)) then
         nrstlo = 1;   nrsthi = 1
         print *, "interpolation_aitken(): Extrapolating, not interpolating."
      elseif (xval > xarr(narr)) then
         nrstlo = narr;   nrsthi = narr
         print *, "interpolation_aitken(): Extrapolating, not interpolating."
      else
         k = 0
       1 k = k+1
         if (xarr(k) < xval) then
            nrstlo = k
            goto 1
         else
            nrsthi = k
         end if
      end if
      !
      ! Clear relevant piece of use-indicator array
      llo = max(nrstlo-mxord,   1)
      lhi = min(nrsthi+mxord,narr)
      llr = llo - 1
      do  k = llo,lhi
         used(k-llr) = .false.
      end do
      !
      ! Determine next-nearest XARR entry
      do  irow = 1,mxord
         llo = max(nrstlo-irow+1,   1)
         lhi = min(nrsthi+irow-1,narr)
         set = .false.
         do  k = llo,lhi
            if (.not.used(k-llr)) then
               if (.not.set) then
                  diff = xarr(k) - xval
                  locnxt = k
                  set = .true.
               else
                  difft = xarr(k) - xval
                  if (abs(difft) < abs(diff)) then
                     diff = difft
                     locnxt = k
                  end if
               end if
            end if
         end do
         used(locnxt-llr) = .true.
         x(irow)  = xarr(locnxt)
         dx(irow) = diff
         !
         ! Fill table for this row
         do  k = 1,irow
            ilirok = iloc(irow,k)
            if (k == 1) then
               poly(ilirok) = yarr(locnxt)
            else
               ildiag       = iloc(k-1,k-1)
               ilothr       = iloc(irow,k-1)
               poly(ilirok) = (poly(ildiag)*dx(irow) - poly(ilothr)*dx(k-1)) / (x(irow)-x(k-1))
            endif
         end do
         !
         ! Pick off the diagonal element
         ildiag    = iloc(irow,irow)
         est(irow) = poly(ildiag)
      end do
      !
      ! Now the estimate vector is filled in, so obtain the best estimate
      debeb = abs((est(2) - est(1)) / est(2))
      ibest = 2
      do  irow = 3,mxord
         debe = abs((est(irow) - est(irow-1)) / est(irow))
         if (debe < debeb) then
            debeb = debe
            ibest = irow
         end if
      end do
      yval = est(ibest)
      !
      if (present(accy)) then
         if (debeb > accy) then
            write(*,2) debeb, accy
          2 format( "interpolation_aitken(): Accuracy of interpolation (", e10.3,") is below input criterion (",e10.3,").")
         end if
      end if
      !
      contains
         !
         function iloc (ind1,ind2)                           result(loc)
         !--------------------------------------------------------------------------------------------------------------
         ! This internal function dispenses with the need for a two-dimensional array for the interpolation. It replaces a
         ! statement function in the original code.
         !--------------------------------------------------------------------------------------------------------------
         !
         integer, intent(in) :: ind1, ind2
         integer             :: loc
         !
         loc = (ind1*(ind1-1)) / 2 + ind2
         !
         end function iloc
         !
   end subroutine interpolation_aitken
   !
   !
   subroutine interpolation_lagrange(arg,val,x,y,n)
   !--------------------------------------------------------------------------------------------------------------------
   ! Interpolates the value val(x) from the values arg(i),val(i), i=1,n by using a Lagrange-interpolation formulae. 
   !--------------------------------------------------------------------------------------------------------------------
   !
   integer, intent(in)                     :: n
   real(kind=dp), intent(in)               :: x
   real(kind=dp), intent(out)              :: y
   real(kind=dp), dimension(:), intent(in) :: arg, val
   !
   integer       :: j, l
   real(kind=dp) :: pl
   !
   y = zero
   do  l = 1,n
      pl = one
      do  j = 1,n
         if (l-j /= 0) then
            pl = (x - arg(j)) * pl / (arg(l) - arg(j))
         end if
      end do
      y = y + pl * val(l)
   end do
   !
   end subroutine interpolation_lagrange
   !
   !
   function is_triangle(i1,i2,i3)                            result(yes)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns .true. if the lengths i1, i2, and i3 may form a triangle  and .false. otherwise.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: i1, i2, i3
      logical             :: yes
      !
      if (i1 <= i2+i3   .and.   i2 <= i3+i1   .and.   i3 <= i1+i2) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function is_triangle
   !
   !
   function legendre_associated(n,m,z)                     result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the (complex) value of the Legendre function P_n^m (z) for complex argument z.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)              :: n, m
      complex(kind=dp), intent(in)     :: z
      complex(kind=dp)                 :: value
      !
      integer                          :: i
      complex(kind=dp), dimension(0:n) :: Pm
      !
      if (rabs_use_stop   .and. &
         (n < 0   .or.   m /= 1)) then
         stop "legendre_associated(): program stop A."
      else if (n == 0) then
         value = zero
      else if (n == 1) then
         value = sqrt(one - z*z)
      else
         Pm(0) = zero;   Pm(1) = sqrt(one - z*z)
	 do  i = 1,n-1
	    Pm(i+1) = ((two*i + one) *z* Pm(i) - (i+m)*Pm(i-1)) / (i - m + one)
	 end do
	 value = Pm(n)
      end if
      !
   end function legendre_associated
   !
   !
   function legendre_polynomial(n,z)                                         result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the (complex) value of the Legendre polynomial P_n (z) for complex argument z.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)              :: n
      complex(kind=dp), intent(in)     :: z
      complex(kind=dp)                 :: value
      !
      integer                          :: i
      complex(kind=dp), dimension(0:n) :: P
      !
      if (rabs_use_stop   .and.   n < 0) then
         stop "legendre_polynomial(): program stop A."
      else if (n == 0) then
         value = one
      else if (n == 1) then
         value = z
      else
         P(0) = one;   P(1) = z
	 do  i = 1,n-1
	    P(i+1) = ((two*i + one) * z * P(i) - i * P(i-1)) / (i + one)
	 end do
	 value = P(n)
      end if
      !
   end function legendre_polynomial
   !
   !
   subroutine LU_decomposition_of_matrix(matrix,n,d)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the LU decomposition of the n x n matrix; ndim is the dimension in the declaration of matrix. This routine 
   ! has been adopted from an algorithmus as discussed in the 'Numerical Recipies'. The real, intent(out) parameter 
   ! d = 1.0 or d = -1.0 denotes the number of permutations and may be used to calculate the determinant as d times the 
   ! product of all diagonal matrix elements (at output time).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)        :: n
      real(kind=dp), intent(out) :: d
      real(kind=dp), dimension(:,:), intent(inout) :: matrix
      !
      integer                              :: i, imax, j, k
      real(kind=dp)                        :: a_max, dummy, sum
      integer, dimension(1:size(matrix,1)) :: index
      real(kind=dp), dimension(n)          :: vv
      !
      d = one
      do  i = 1,n
         a_max = zero
         do  j = 1,n;   a_max = max(a_max, abs(matrix(i,j)));   end do
         if (a_max == zero) then
            print *, "LU_decomposition_of_matrix(): singular matrix."
            ! stop
         else
            vv(i) = one / a_max
         end if
      end do
      !
      do  j = 1,n
         do  i = 1,j-1
            sum = matrix(i,j)
            do  k = 1,i-1;   sum = sum - matrix(i,k)*matrix(k,j);   end do
            matrix(i,j) = sum
         end do
         a_max = zero
         do  i = j,n
            sum = matrix(i,j)
            do  k = 1,j-1;   sum = sum - matrix(i,k)*matrix(k,j);   end do
            matrix(i,j) = sum
            dummy = vv(i) * abs(sum)
            if (dummy >= a_max) then
               imax = i;   a_max = dummy
            end if
         end do
         if (j /= imax) then
            do  k = 1,n
               dummy = matrix(imax,k);   matrix(imax,k) = matrix(j,k)
               matrix(j,k) = dummy
            end do
            d = -d;   vv(imax) = vv(j)
         end if
         index(j) = imax
         if (matrix(j,j) == zero)  matrix(j,j) = eps20
         if (j /= n) then
            dummy = one / matrix(j,j)
            do  i = j+1,n;   matrix(i,j) = matrix(i,j) * dummy;   end do
         end if
      end do
      !
   end subroutine LU_decomposition_of_matrix
   !      
   !
   function pochhammer(x,n)                                   result(ph)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of Pochhammer's symbol (z)_n for real arguments of x.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: n
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: ph
      !
      if (n == 0) then
         ph = one
      else if (rabs_use_naglib) then
         ph = nag_s14aaf(x+n,ifail) / nag_s14aaf(x,ifail)
      else 
         stop "pochhammer(): program stop A."
      end if
      !
   end function pochhammer
   !
   !
   subroutine set_beta_gamma_arrays()
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram allocates and initializes some arrays for the calculation of the Beta function, the incomplete Beta 
   ! function, as well as the Gamma function.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer          :: a, b, m
      real(kind=dp)    :: wa, wb, wm, wa2m
      !! real(kind=dp) :: nag_s14aaf
      !
      ! Allocate the arrays
      allocate ( gamma_store(2*ndbeta) )
      allocate ( beta_store(ndbeta,ndbeta) )
      allocate ( beta_frac_coeff(ndmx,ndbeta,ndbeta) )
      !
      ! The array gamma_store(n) contains the value Gamma(n/2)
      do  a = 1, 2*ndbeta
         wa = a / two
         !! if (rabs_use_naglib) then
            gamma_store(a) = nag_s14aaf(wa,ifail)
         !! else if (rabs_use_stop) then
         !!    stop "set_beta_gamma_arrays(): program stop A."
         !! end if
      end do
      !
      ! The array beta_store(a,b) contains the value Beta(a/2,b/2)
      do  b = 1, ndbeta
         do  a = 1, ndbeta
            beta_store(a,b) = gamma_store(a) * gamma_store(b) / gamma_store(a+b)
         end do
      enddo
      !
      ! The array beta_frac_coeff(n,a,b) contains fractional coefficients for the calculation of the incomplete 
      ! Beta(x,a,b) function with integer and half-integer arguments a,b
      !
      do b = 1,ndbeta
         do  a = 1,ndbeta
            wa = a / two
            wb = b / two
            beta_frac_coeff(1,a,b) = -(wa + wb) / (wa + one)
            do  m = 1,(ndmx/2)
               wm   = m
               wa2m = wa + wm + wm
               beta_frac_coeff(m+m,a,b)   = wm * (wb - wm) /  ((wa2m - one) * wa2m)
               beta_frac_coeff(m+m+1,a,b) = - (wa + wm) * (wa + wb + wm) / (wa2m * (wa2m + one) )
            end do
         end do
      end do
      !
   end subroutine set_beta_gamma_arrays
   !
   !
   subroutine set_gauss_zeros(ax,bx,nx,zw,w)
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram set n-point gauss zeros and weights for the interval (ax,bx) into the arrays z(1:nx) and w(1:nx) 
   ! respectively.
   !
   ! Calls: set_gauss_zeros_aux
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer                                     :: nx
      real(kind=dp), intent(in)                   :: ax, bx
      real(kind=dp), dimension(1:nx), intent(out) :: zw, w
      !
      integer                         :: k, m, n, j, jp, jmid, jtab
      integer, dimension(1:96)        :: ltab 
      real(kind=dp)                   :: alpha, beta, delta, wtemp
      real(kind=dp), dimension(1:273) :: b, y 
      !
      call set_gauss_zeros_aux(b,y,ltab)
      !
      ! Test n
      n = nx
      alpha = half *(bx+ax)
      beta  = half *(bx-ax)
      if (n < 1) then
         print *, "set_gauss_zeros(): n has a non-permissible value of ",n 
         stop
      else if (n /= 1) then
         if ( n <= 16  .or.  n == 20  .or.  n == 24  .or.  n == 32  .or.              &
              n == 40  .or.  n == 48  .or.  n == 64  .or.  n == 80  .or. n == 96 ) then  
            goto 2
         else
            print *, "set_gauss_zeros(): n has a non-permissible value of ",n 
            stop
         end if
      else
         zw(1) = alpha
         w(1)  = bx - ax
         goto 3
      end if
      !
      ! Set  k  equal to initial subscript and store results
    2 k = ltab(n)
      m = n / 2
      !
      do  j = 1,m
         jtab  = k-1+j
         wtemp = beta * b(jtab)
         delta = beta * y(jtab)
         zw(j) = alpha - delta
         w(j)  = wtemp
         jp    = n+1-j
         zw(jp)= alpha + delta
         w(jp) = wtemp
      end do
      !
      if (n-m-m == 0) then
      else
         zw(m+1) = alpha
         jmid    = k+m
         w(m+1)  = beta * b(jmid)
      end if
      ! 
    3 continue
      !
   end subroutine set_gauss_zeros
   !
   !
   subroutine set_gauss_zeros_aux(b,y,ltab)
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram is an auxilarity routine for set_gauss_zeros to set n-point Gauss zeros and weights.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, dimension(1:96), intent(out)        :: ltab 
      real(kind=dp), dimension(1:273), intent(out) :: b(273), y(273) 
      !
      integer                               :: i
      integer, dimension(1:96), save        :: ktab 
      real(kind=dp), dimension(1:273), save :: a(273), x(273) 
      !
      ! Table of initial subscripts for n= 2 (1) 16 (4) 96
      data ktab(2)/1/;      data ktab(3)/2/;      data ktab(4)/4/
      data ktab(5)/6/;      data ktab(6)/9/;      data ktab(7)/12/
      data ktab(8)/16/;     data ktab(9)/20/;     data ktab(10)/25/
      data ktab(11)/30/;    data ktab(12)/36/;    data ktab(13)/42/
      data ktab(14)/49/;    data ktab(15)/56/;    data ktab(16)/64/
      data ktab(20)/72/;    data ktab(24)/82/;    data ktab(28)/82/
      data ktab(32)/94/;    data ktab(36)/94/;    data ktab(40)/110/
      data ktab(44)/110/;   data ktab(48)/130/;   data ktab(52)/130/
      data ktab(56)/130/;   data ktab(60)/130/;   data ktab(64)/154/
      data ktab(68)/154/;   data ktab(72)/154/;   data ktab(76)/154/
      data ktab(80)/186/;   data ktab(84)/186/;   data ktab(88)/186/
      data ktab(92)/186/;   data ktab(96)/226/
      !
      ! Table of abscissae (x) and weights (a) for interval (-1,+1).
      !
      ! n = 2
      data x(1)/0.577350269189626_dp  /, a(1)/1.000000000000000_dp  /
      ! n = 3
      data x(2)/0.774596669241483_dp  /, a(2)/0.555555555555556_dp  /
      data x(3)/0.000000000000000_dp  /, a(3)/0.888888888888889_dp  /
      ! n = 4
      data x(4)/0.861136311594053_dp  /, a(4)/0.347854845137454_dp  /
      data x(5)/0.339981043584856_dp  /, a(5)/0.652145154862546_dp  /
      ! n = 5
      data x(6)/0.906179845938664_dp  /, a(6)/0.236926885056189_dp  /
      data x(7)/0.538469310105683_dp  /, a(7)/0.478628670499366_dp  /
      data x(8)/0.000000000000000_dp  /, a(8)/0.568888888888889_dp  /
      ! n = 6
      data x(9)/0.932469514203152_dp  /, a(9)/0.171324492379170_dp  /
      data x(10)/0.661209386466265_dp /, a(10)/0.360761573048139_dp /
      data x(11)/0.238619186083197_dp /, a(11)/0.467913934572691_dp /
      ! n = 7
      data x(12)/0.949107912342759_dp /, a(12)/0.129484966168870_dp /
      data x(13)/0.741531185599394_dp /, a(13)/0.279705391489277_dp /
      data x(14)/0.405845151377397_dp /, a(14)/0.381830050505119_dp /
      data x(15)/0.000000000000000_dp /, a(15)/0.417959183673469_dp /
      ! n = 8
      data x(16)/0.960289856497536_dp /, a(16)/0.101228536290376_dp /
      data x(17)/0.796666477413627_dp /, a(17)/0.222381034453374_dp /
      data x(18)/0.525532409916329_dp /, a(18)/0.313706645877887_dp /
      data x(19)/0.183434642495650_dp /, a(19)/0.362683783378362_dp /
      ! n = 9
      data x(20)/0.968160239507626_dp /, a(20)/0.081274388361574_dp /
      data x(21)/0.836031107326636_dp /, a(21)/0.180648160694857_dp /
      data x(22)/0.613371432700590_dp /, a(22)/0.260610696402935_dp /
      data x(23)/0.324253423403809_dp /, a(23)/0.312347077040003_dp /
      data x(24)/0.000000000000000_dp /, a(24)/0.330239355001260_dp /
      ! n = 10
      data x(25)/0.973906528517172_dp /, a(25)/0.066671344308688_dp /
      data x(26)/0.865063366688985_dp /, a(26)/0.149451349150581_dp /
      data x(27)/0.679409568299024_dp /, a(27)/0.219086362515982_dp /
      data x(28)/0.433395394129247_dp /, a(28)/0.269266719309996_dp /
      data x(29)/0.148874338981631_dp /, a(29)/0.295524224714753_dp /
      ! n = 11
      data x(30)/0.978228658146057_dp /, a(30)/0.055668567116174_dp /
      data x(31)/0.887062599768095_dp /, a(31)/0.125580369464905_dp /
      data x(32)/0.730152005574049_dp /, a(32)/0.186290210927734_dp /
      data x(33)/0.519096129206812_dp /, a(33)/0.233193764591990_dp /
      data x(34)/0.269543155952345_dp /, a(34)/0.262804544510247_dp /
      data x(35)/0.000000000000000_dp /, a(35)/0.272925086777901_dp /
      ! n = 12
      data x(36)/0.981560634246719_dp /, a(36)/0.047175336386512_dp /
      data x(37)/0.904117256370475_dp /, a(37)/0.106939325995318_dp /
      data x(38)/0.769902674194305_dp /, a(38)/0.160078328543346_dp /
      data x(39)/0.587317954286617_dp /, a(39)/0.203167426723066_dp /
      data x(40)/0.367831498998180_dp /, a(40)/0.233492536538355_dp /
      data x(41)/0.125233408511469_dp /, a(41)/0.249147045813403_dp /
      ! n = 13
      data x(42)/0.984183054718588_dp /, a(42)/0.040484004765316_dp /
      data x(43)/0.917598399222978_dp /, a(43)/0.092121499837728_dp /
      data x(44)/0.801578090733310_dp /, a(44)/0.138873510219787_dp /
      data x(45)/0.642349339440340_dp /, a(45)/0.178145980761946_dp /
      data x(46)/0.448492751036447_dp /, a(46)/0.207816047536889_dp /
      data x(47)/0.230458315955135_dp /, a(47)/0.226283180262897_dp /
      data x(48)/0.000000000000000_dp /, a(48)/0.232551553230874_dp /
      ! n = 14
      data x(49)/0.986283808696812_dp /, a(49)/0.035119460331752_dp /
      data x(50)/0.928434883663574_dp /, a(50)/0.080158087159760_dp /
      data x(51)/0.827201315069765_dp /, a(51)/0.121518570687903_dp /
      data x(52)/0.687292904811685_dp /, a(52)/0.157203167158194_dp /
      data x(53)/0.515248636358154_dp /, a(53)/0.185538397477938_dp /
      data x(54)/0.319112368927890_dp /, a(54)/0.205198463721296_dp /
      data x(55)/0.108054948707344_dp /, a(55)/0.215263853463158_dp /
      ! n = 15
      data x(56)/0.987992518020485_dp /, a(56)/0.030753241996117_dp /
      data x(57)/0.937273392400706_dp /, a(57)/0.070366047488108_dp /
      data x(58)/0.848206583410427_dp /, a(58)/0.107159220467172_dp /
      data x(59)/0.724417731360170_dp /, a(59)/0.139570677926154_dp /
      data x(60)/0.570972172608539_dp /, a(60)/0.166269205816994_dp /
      data x(61)/0.394151347077563_dp /, a(61)/0.186161000015562_dp /
      data x(62)/0.201194093997435_dp /, a(62)/0.198431485327111_dp /
      data x(63)/0.000000000000000_dp /, a(63)/0.202578241925561_dp /
      ! n = 16
      data x(64)/0.989400934991650_dp /, a(64)/0.027152459411754_dp /
      data x(65)/0.944575023073233_dp /, a(65)/0.062253523938648_dp /
      data x(66)/0.865631202387832_dp /, a(66)/0.095158511682493_dp /
      data x(67)/0.755404408355003_dp /, a(67)/0.124628971255534_dp /
      data x(68)/0.617876244402644_dp /, a(68)/0.149595988816577_dp /
      data x(69)/0.458016777657227_dp /, a(69)/0.169156519395003_dp /
      data x(70)/0.281603550779259_dp /, a(70)/0.182603415044924_dp /
      data x(71)/0.095012509837637_dp /, a(71)/0.189450610455069_dp /
      ! n = 20
      data x(72)/0.993128599185094_dp /, a(72)/0.017614007139152_dp /
      data x(73)/0.963971927277913_dp /, a(73)/0.040601429800386_dp /
      data x(74)/0.912234428251325_dp /, a(74)/0.062672048334109_dp /
      data x(75)/0.839116971822218_dp /, a(75)/0.083276741576704_dp /
      data x(76)/0.746331906460150_dp /, a(76)/0.101930119817240_dp /
      data x(77)/0.636053680726515_dp /, a(77)/0.118194531961518_dp /
      data x(78)/0.510867001950827_dp /, a(78)/0.131688638449176_dp /
      data x(79)/0.373706088715419_dp /, a(79)/0.142096109318382_dp /
      data x(80)/0.227785851141645_dp /, a(80)/0.149172986472603_dp /
      data x(81)/0.076526521133497_dp /, a(81)/0.152753387130725_dp /
      ! n = 24
      data x(82)/0.995187219997021_dp /, a(82)/0.012341229799987_dp /
      data x(83)/0.974728555971309_dp /, a(83)/0.028531388628933_dp /
      data x(84)/0.938274552002732_dp /, a(84)/0.044277438817419_dp /
      data x(85)/0.886415527004401_dp /, a(85)/0.059298584915436_dp /
      data x(86)/0.820001985973902_dp /, a(86)/0.073346481411080_dp /
      data x(87)/0.740124191578554_dp /, a(87)/0.086190161531953_dp /
      data x(88)/0.648093651936975_dp /, a(88)/0.097618652104113_dp /
      data x(89)/0.545421471388839_dp /, a(89)/0.107444270115965_dp /
      data x(90)/0.433793507626045_dp /, a(90)/0.115505668053725_dp /
      data x(91)/0.315042679696163_dp /, a(91)/0.121670472927803_dp /
      data x(92)/0.191118867473616_dp /, a(92)/0.125837456346828_dp /
      data x(93)/0.064056892862605_dp /, a(93)/0.127938195346752_dp /
      ! n = 32
      data x(94)/0.997263861849481_dp /, a(94)/0.007018610009470_dp /
      data x(95)/0.985611511545268_dp /, a(95)/0.016274394730905_dp /
      data x(96)/0.964762255587506_dp /, a(96)/0.025392065309262_dp /
      data x(97)/0.934906075937739_dp /, a(97)/0.034273862913021_dp /
      data x(98)/0.896321155766052_dp /, a(98)/0.042835898022226_dp /
      data x(99)/0.849367613732569_dp /, a(99)/0.050998059262376_dp /
      data x(100)/0.794483795967942_dp/, a(100)/0.058684093478535_dp/
      data x(101)/0.732182118740289_dp/, a(101)/0.065822222776361_dp/
      data x(102)/0.663044266930215_dp/, a(102)/0.072345794108848_dp/
      data x(103)/0.587715757240762_dp/, a(103)/0.078193895787070_dp/
      data x(104)/0.506899908932229_dp/, a(104)/0.083311924226946_dp/
      data x(105)/0.421351276130635_dp/, a(105)/0.087652093004403_dp/
      data x(106)/0.331868602282127_dp/, a(106)/0.091173878695763_dp/
      data x(107)/0.239287362252137_dp/, a(107)/0.093844399080804_dp/
      data x(108)/0.144471961582796_dp/, a(108)/0.095638720079274_dp/
      data x(109)/0.048307665687738_dp/, a(109)/0.096540088514727_dp/
      ! n = 40
      data x(110)/0.998237709710559_dp/, a(110)/0.004521277098533_dp/
      data x(111)/0.990726238699457_dp/, a(111)/0.010498284531152_dp/
      data x(112)/0.977259949983774_dp/, a(112)/0.016421058381907_dp/
      data x(113)/0.957916819213791_dp/, a(113)/0.022245849194166_dp/
      data x(114)/0.932812808278676_dp/, a(114)/0.027937006980023_dp/
      data x(115)/0.902098806968874_dp/, a(115)/0.033460195282547_dp/
      data x(116)/0.865959503212259_dp/, a(116)/0.038782167974472_dp/
      data x(117)/0.824612230833311_dp/, a(117)/0.043870908185673_dp/
      data x(118)/0.778305651426519_dp/, a(118)/0.048695807635072_dp/
      data x(119)/0.727318255189927_dp/, a(119)/0.053227846983936_dp/
      data x(120)/0.671956684614179_dp/, a(120)/0.057439769099391_dp/
      data x(121)/0.612553889667980_dp/, a(121)/0.061306242492928_dp/
      data x(122)/0.549467125095128_dp/, a(122)/0.064804013456601_dp/
      data x(123)/0.483075801686178_dp/, a(123)/0.067912045815233_dp/
      data x(124)/0.413779204371605_dp/, a(124)/0.070611647391286_dp/
      data x(125)/0.341994090825758_dp/, a(125)/0.072886582395804_dp/
      data x(126)/0.268152185007253_dp/, a(126)/0.074723169057968_dp/
      data x(127)/0.192697580701371_dp/, a(127)/0.076110361900626_dp/
      data x(128)/0.116084070675255_dp/, a(128)/0.077039818164247_dp/
      data x(129)/0.038772417506050_dp/, a(129)/0.077505947978424_dp/
      ! n = 48
      data x(130)/0.998771007252426_dp/, a(130)/0.003153346052305_dp/
      data x(131)/0.993530172266350_dp/, a(131)/0.007327553901276_dp/
      data x(132)/0.984124583722826_dp/, a(132)/0.011477234579234_dp/
      data x(133)/0.970591592546247_dp/, a(133)/0.015579315722943_dp/
      data x(134)/0.952987703160430_dp/, a(134)/0.019616160457355_dp/
      data x(135)/0.931386690706554_dp/, a(135)/0.023570760839324_dp/
      data x(136)/0.905879136715569_dp/, a(136)/0.027426509708356_dp/
      data x(137)/0.876572020274247_dp/, a(137)/0.031167227832798_dp/
      data x(138)/0.843588261624393_dp/, a(138)/0.034777222564770_dp/
      data x(139)/0.807066204029442_dp/, a(139)/0.038241351065830_dp/
      data x(140)/0.767159032515740_dp/, a(140)/0.041545082943464_dp/
      data x(141)/0.724034130923814_dp/, a(141)/0.044674560856694_dp/
      data x(142)/0.677872379632663_dp/, a(142)/0.047616658492490_dp/
      data x(143)/0.628867396776513_dp/, a(143)/0.050359035553854_dp/
      data x(144)/0.577224726083972_dp/, a(144)/0.052890189485193_dp/
      data x(145)/0.523160974722233_dp/, a(145)/0.055199503699984_dp/
      data x(146)/0.466902904750958_dp/, a(146)/0.057277292100403_dp/
      data x(147)/0.408686481990716_dp/, a(147)/0.059114839698395_dp/
      data x(148)/0.348755886292160_dp/, a(148)/0.060704439165893_dp/
      data x(149)/0.287362487355455_dp/, a(149)/0.062039423159892_dp/
      data x(150)/0.224763790394689_dp/, a(150)/0.063114192286254_dp/
      data x(151)/0.161222356068891_dp/, a(151)/0.063924238584648_dp/
      data x(152)/0.097004699209462_dp/, a(152)/0.064466164435950_dp/
      data x(153)/0.032380170962869_dp/, a(153)/0.064737696812683_dp/
      ! n = 64
      data x(154)/0.999305041735772_dp/, a(154)/0.001783280721696_dp/
      data x(155)/0.996340116771955_dp/, a(155)/0.004147033260562_dp/
      data x(156)/0.991013371476744_dp/, a(156)/0.006504457968978_dp/
      data x(157)/0.983336253884625_dp/, a(157)/0.008846759826363_dp/
      data x(158)/0.973326827789910_dp/, a(158)/0.011168139460131_dp/
      data x(159)/0.961008799652053_dp/, a(159)/0.013463047896718_dp/
      data x(160)/0.946411374858402_dp/, a(160)/0.015726030476024_dp/
      data x(161)/0.929569172131939_dp/, a(161)/0.017951715775697_dp/
      data x(162)/0.910522137078502_dp/, a(162)/0.020134823153530_dp/
      data x(163)/0.889315445995114_dp/, a(163)/0.022270173808383_dp/
      data x(164)/0.865999398154092_dp/, a(164)/0.024352702568710_dp/
      data x(165)/0.840629296252580_dp/, a(165)/0.026377469715054_dp/
      data x(166)/0.813265315122797_dp/, a(166)/0.028339672614259_dp/
      data x(167)/0.783972358943341_dp/, a(167)/0.030234657072402_dp/
      data x(168)/0.752819907260531_dp/, a(168)/0.032057928354851_dp/
      data x(169)/0.719881850171610_dp/, a(169)/0.033805161837141_dp/
      data x(170)/0.685236313054233_dp/, a(170)/0.035472213256882_dp/
      data x(171)/0.648965471254657_dp/, a(171)/0.037055128540240_dp/
      data x(172)/0.611155355172393_dp/, a(172)/0.038550153178615_dp/
      data x(173)/0.571895646202634_dp/, a(173)/0.039953741132720_dp/
      data x(174)/0.531279464019894_dp/, a(174)/0.041262563242623_dp/
      data x(175)/0.489403145707052_dp/, a(175)/0.042473515123653_dp/
      data x(176)/0.446366017253464_dp/, a(176)/0.043583724529323_dp/
      data x(177)/0.402270157963991_dp/, a(177)/0.044590558163756_dp/
      data x(178)/0.357220158337668_dp/, a(178)/0.045491627927418_dp/
      data x(179)/0.311322871990210_dp/, a(179)/0.046284796581314_dp/
      data x(180)/0.264687162208767_dp/, a(180)/0.046968182816210_dp/
      data x(181)/0.217423643740007_dp/, a(181)/0.047540165714830_dp/
      data x(182)/0.169644420423992_dp/, a(182)/0.047999388596458_dp/
      data x(183)/0.121462819296120_dp/, a(183)/0.048344762234802_dp/
      data x(184)/0.072993121787799_dp/, a(184)/0.048575467441503_dp/
      data x(185)/0.024350292663424_dp/, a(185)/0.048690957009139_dp/
      ! n = 80
      data x(186)/0.999553822651630_dp/, a(186)/0.001144950003186_dp/
      data x(187)/0.997649864398237_dp/, a(187)/0.002663533589512_dp/
      data x(188)/0.994227540965688_dp/, a(188)/0.004180313124694_dp/
      data x(189)/0.989291302499755_dp/, a(189)/0.005690922451403_dp/
      data x(190)/0.982848572738629_dp/, a(190)/0.007192904768117_dp/
      data x(191)/0.974909140585727_dp/, a(191)/0.008683945269260_dp/
      data x(192)/0.965485089043799_dp/, a(192)/0.010161766041103_dp/
      data x(193)/0.954590766343634_dp/, a(193)/0.011624114120797_dp/
      data x(194)/0.942242761309872_dp/, a(194)/0.013068761592401_dp/
      data x(195)/0.928459877172445_dp/, a(195)/0.014493508040509_dp/
      data x(196)/0.913263102571757_dp/, a(196)/0.015896183583725_dp/
      data x(197)/0.896675579438770_dp/, a(197)/0.017274652056269_dp/
      data x(198)/0.878722567678213_dp/, a(198)/0.018626814208299_dp/
      data x(199)/0.859431406663111_dp/, a(199)/0.019950610878141_dp/
      data x(200)/0.838831473580255_dp/, a(200)/0.021244026115782_dp/
      data x(201)/0.816954138681463_dp/, a(201)/0.022505090246332_dp/
      data x(202)/0.793832717504605_dp/, a(202)/0.023731882865930_dp/
      data x(203)/0.769502420135041_dp/, a(203)/0.024922535764115_dp/
      data x(204)/0.744000297583597_dp/, a(204)/0.026075235767565_dp/
      data x(205)/0.717365185362099_dp/, a(205)/0.027188227500486_dp/
      data x(206)/0.689637644342027_dp/, a(206)/0.028259816057276_dp/
      data x(207)/0.660859898986119_dp/, a(207)/0.029288369583267_dp/
      data x(208)/0.631075773046871_dp/, a(208)/0.030272321759557_dp/
      data x(209)/0.600330622829751_dp/, a(209)/0.031210174188114_dp/
      data x(210)/0.568671268122709_dp/, a(210)/0.032100498673487_dp/
      data x(211)/0.536145920897131_dp/, a(211)/0.032941939397645_dp/
      data x(212)/0.502804111888784_dp/, a(212)/0.033733214984611_dp/
      data x(213)/0.468696615170544_dp/, a(213)/0.034473120451753_dp/
      data x(214)/0.433875370831756_dp/, a(214)/0.035160529044747_dp/
      data x(215)/0.398393405881969_dp/, a(215)/0.035794393953416_dp/
      data x(216)/0.362304753499487_dp/, a(216)/0.036373749905835_dp/
      data x(217)/0.325664370747701_dp/, a(217)/0.036897714638276_dp/
      data x(218)/0.288528054884511_dp/, a(218)/0.037365490238730_dp/
      data x(219)/0.250952358392272_dp/, a(219)/0.037776364362001_dp/
      data x(220)/0.212994502857666_dp/, a(220)/0.038129711314477_dp/
      data x(221)/0.174712291832646_dp/, a(221)/0.038424993006959_dp/
      data x(222)/0.136164022809143_dp/, a(222)/0.038661759774076_dp/
      data x(223)/0.097408398441584_dp/, a(223)/0.038839651059051_dp/
      data x(224)/0.058504437152420_dp/, a(224)/0.038958395962769_dp/
      data x(225)/0.019511383256793_dp/, a(225)/0.039017813656306_dp/
      ! n = 96
      data x(226)/0.999689503883230_dp/, a(226)/0.000796792065552_dp/
      data x(227)/0.998364375863181_dp/, a(227)/0.001853960788946_dp/
      data x(228)/0.995981842987209_dp/, a(228)/0.002910731817934_dp/
      data x(229)/0.992543900323762_dp/, a(229)/0.003964554338444_dp/
      data x(230)/0.988054126329623_dp/, a(230)/0.005014202742927_dp/
      data x(231)/0.982517263563014_dp/, a(231)/0.006058545504235_dp/
      data x(232)/0.975939174585136_dp/, a(232)/0.007096470791153_dp/
      data x(233)/0.968326828463264_dp/, a(233)/0.008126876925698_dp/
      data x(234)/0.959688291448742_dp/, a(234)/0.009148671230783_dp/
      data x(235)/0.950032717784437_dp/, a(235)/0.010160770535008_dp/
      data x(236)/0.939370339752755_dp/, a(236)/0.011162102099838_dp/
      data x(237)/0.927712456722308_dp/, a(237)/0.012151604671088_dp/
      data x(238)/0.915071423120898_dp/, a(238)/0.013128229566961_dp/
      data x(239)/0.901460635315852_dp/, a(239)/0.014090941772314_dp/
      data x(240)/0.886894517402420_dp/, a(240)/0.015038721026994_dp/
      data x(241)/0.871388505909296_dp/, a(241)/0.015970562902562_dp/
      data x(242)/0.854959033434601_dp/, a(242)/0.016885479864245_dp/
      data x(243)/0.837623511228187_dp/, a(243)/0.017782502316045_dp/
      data x(244)/0.819400310737931_dp/, a(244)/0.018660679627411_dp/
      data x(245)/0.800308744139140_dp/, a(245)/0.019519081140145_dp/
      data x(246)/0.780369043867433_dp/, a(246)/0.020356797154333_dp/
      data x(247)/0.759602341176647_dp/, a(247)/0.021172939892191_dp/
      data x(248)/0.738030643744400_dp/, a(248)/0.021966644438744_dp/
      data x(249)/0.715676812348967_dp/, a(249)/0.022737069658329_dp/
      data x(250)/0.692564536642171_dp/, a(250)/0.023483399085926_dp/
      data x(251)/0.668718310043916_dp/, a(251)/0.024204841792364_dp/
      data x(252)/0.644163403784967_dp/, a(252)/0.024900633222483_dp/
      data x(253)/0.618925840125468_dp/, a(253)/0.025570036005349_dp/
      data x(254)/0.593032364777572_dp/, a(254)/0.026212340735672_dp/
      data x(255)/0.566510418561397_dp/, a(255)/0.026826866725591_dp/
      data x(256)/0.539388108324357_dp/, a(256)/0.027412962726029_dp/
      data x(257)/0.511694177154667_dp/, a(257)/0.027970007616848_dp/
      data x(258)/0.483457973920596_dp/, a(258)/0.028497411065085_dp/
      data x(259)/0.454709422167743_dp/, a(259)/0.028994614150555_dp/
      data x(260)/0.425478988407300_dp/, a(260)/0.029461089958167_dp/
      data x(261)/0.395797649828908_dp/, a(261)/0.029896344136328_dp/
      data x(262)/0.365696861472313_dp/, a(262)/0.030299915420827_dp/
      data x(263)/0.335208522892625_dp/, a(263)/0.030671376123669_dp/
      data x(264)/0.304364944354496_dp/, a(264)/0.031010332586313_dp/
      data x(265)/0.273198812591049_dp/, a(265)/0.031316425596861_dp/
      data x(266)/0.241743156163840_dp/, a(266)/0.031589330770727_dp/
      data x(267)/0.210031310460567_dp/, a(267)/0.031828758894411_dp/
      data x(268)/0.178096882367618_dp/, a(268)/0.032034456231992_dp/
      data x(269)/0.145973714654896_dp/, a(269)/0.032206204794030_dp/
      data x(270)/0.113695850110665_dp/, a(270)/0.032343822568575_dp/
      data x(271)/0.081297495464425_dp/, a(271)/0.032447163714064_dp/
      data x(272)/0.048812985136049_dp/, a(272)/0.032516118713868_dp/
      data x(273)/0.016276744849602_dp/, a(273)/0.032550614492363_dp/
      !
      do i=1,273
         b(i) = a(i)
         y(i) = x(i)
      end do
      do  i=1,96
         ltab(i) = ktab(i)
      end do
      !
   end subroutine set_gauss_zeros_aux
   !
   !
   function spherical_Bessel_jL(L,x)                       result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the spherical Bessel function j_L(x) = sqrt(pi/2x) J_(L+1/2) (x) for real values of x by 
   ! calling a proper function from the NAG library.
   !
   ! Call(s): nag_s17def() [from NAG library].
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: L
      real(kind=qp), intent(in) :: x
      !
      integer                   :: nz
      real(kind=dp)             :: nu, value, test
      complex(kind=qp)          :: zz
      complex(kind=qp), dimension(1) :: cy
      !
      value = one
      !
      if (rabs_use_naglib  .and.   x < 30000.0_dp  .and.   .false.) then
         nu = L + half 
         zz = cmplx(x,zero,kind=qp)
	 print *, "nu, zz = ",nu, zz
         call nag_s17def(nu,zz,1,"s",cy(1),nz,ifail)
         !!x print *, "Bessel, L, x, cy = ", L, x, cy
         !
         if (ifail /= 0  .or.  nz /= 0) then
           print *, "spherical_Bessel_jL(): program stop A."
         end if
         !
         value = sqrt( pi/(two*x) ) * cy(1)
         !
         ! Return here in all standard calculations
         return
      end if
      !
      ! Calculate directly
      select case(L)
      case(0)
         if (x < 1.0e-3) then
            test = one - x*x/6.0_dp + x**4/120.0_dp - x**6/5040.0_dp
         else
            test = sin(x) / x
         end if
      case(1)
         if (x < 1.0e-3) then
            test = x/three - x*x*x/30.0_dp + x**5/840.0_dp - x**7/45360.0_dp
         else
            test = sin(x) / (x*x) - cos(x)/x
         end if
      case(2)
         if (x < 1.0e-2) then
            test = x*x/15.0_dp - x**4/210.0_dp + x**6/7560.0_dp - x**8/498960.0_dp
         else
            test = three*sin(x)/(x**3) - three*cos(x)/(x*x) - sin(x)/x
         end if
      case(3)
         if (x < 9.0e-2) then
            test = x**3/105.0_dp - x**5/1890.0_dp + x**7/83160.0_dp - x**9/6486480.0_dp
         else
            test = cos(x)/x - 6.0_dp*sin(x)/(x*x) - 15.0_dp*cos(x)/(x**3) + 15.0_dp*sin(x)/(x**4)
         end if
      case(4)
         if (x < 3.0e-1) then
            test = x**4/945.0_dp - x**6/20790.0_dp + x**8/1081080.0_dp - x**10/97297200.0_dp
         else
            test = sin(x)/x + 10.0_dp*cos(x)/(x*x) - 45.0_dp*sin(x)/(x**3) - 105.0_dp*cos(x)/(x**4) + 105.0_dp*sin(x)/(x**5) 
         end if
      case(5)
         if (x < 8.0e-1) then
            test = x**5/10395.0_dp - x**7/270270.0_dp + x**9/16216200.0_dp - x**11/1654052400.0_dp
         else
            test = -cos(x)/x + 15.0_dp*sin(x)/(x*x) + 105.0_dp*cos(x)/(x**3) -            &
                   420.0_dp*sin(x)/(x**4) - 945.0_dp*cos(x)/(x**5) + 945.0_dp*sin(x)/(x**6)
         end if
      case default
         stop "spherical_Bessel_jL():"
      end select
      !
      value = test
      !
      !! ! Compare results
      !! if (abs((value - test)/value) > 1.0e-6_dp) then
      !!    print *, "spherical_Bessel_jL(): L, x, nag, direct = ", L, x, value, test
      !! end if
      !
   end function spherical_Bessel_jL
   !
   !
   function spherical_Bessel_jL_old(L,x)                    result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the spherical Bessel function j_L(x) = sqrt(pi/2x) J_(L+1/2) (x) for real values of x by calling
   ! a proper function from the NAG library.
   !
   ! Call(s): nag_s17def() [from NAG library].
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: L
      real(kind=dp), intent(in) :: x
      !
      integer          :: nz
      real(kind=dp)    :: value, test
      complex*16       :: zz, cy
      character*1      :: scale
      DOUBLE PRECISION :: nu
      !
      if (rabs_use_naglib  .and.   x < 30000.0_dp) then
         scale = "s"
         nu = L + half
         zz = cmplx(x,zero,kind=dp)
	 !
	 !! print *, "a:: nu,zz,1,scale,cy,nz,ifail = ", nu,zz,1,scale,cy,nz,ifail
         call nag_s17def(nu,zz,1,scale,cy,nz,ifail)
         !
         if (ifail /= 0  .or.  nz /= 0) then
           print *, "spherical_Bessel_jL(): program stop A."
         end if
         !
         value = sqrt( pi/(two*x) ) * cy
         !
         ! Return here in all standard calculations
         return
      end if
      !
      ! Calculate directly
      select case(L)
      case(0)
         if (x < 1.0e-3) then
            test = one - x*x/6.0_dp + x**4/120.0_dp - x**6/5040.0_dp
         else
            test = sin(x) / x
         end if
      case(1)
         if (x < 1.0e-3) then
            test = x/three - x*x*x/30.0_dp + x**5/840.0_dp - x**7/45360.0_dp
         else
            test = sin(x) / (x*x) - cos(x)/x
         end if
      case(2)
         if (x < 1.0e-2) then
            test = x*x/15.0_dp - x**4/210.0_dp + x**6/7560.0_dp - x**8/498960.0_dp
         else
            test = three*sin(x)/(x**3) - three*cos(x)/(x*x) - sin(x)/x
         end if
      case(3)
         if (x < 9.0e-2) then
            test = x**3/105.0_dp - x**5/1890.0_dp + x**7/83160.0_dp - &
                   x**9/6486480.0_dp
         else
            test = cos(x)/x - 6.0_dp*sin(x)/(x*x) - 15.0_dp*cos(x)/(x**3) + 15.0_dp*sin(x)/(x**4)
         end if
      case(4)
         if (x < 3.0e-1) then
            test = x**4/945.0_dp - x**6/20790.0_dp + x**8/1081080.0_dp - x**10/97297200.0_dp
         else
            test = sin(x)/x + 10.0_dp*cos(x)/(x*x) - 45.0_dp*sin(x)/(x**3) - 105.0_dp*cos(x)/(x**4) + 105.0_dp*sin(x)/(x**5) 
         end if
      case(5)
         if (x < 8.0e-1) then
            test = x**5/10395.0_dp - x**7/270270.0_dp + x**9/16216200.0_dp - x**11/1654052400.0_dp
         else
            test = -cos(x)/x + 15.0_dp*sin(x)/(x*x) + 105.0_dp*cos(x)/(x**3) -            &
                   420.0_dp*sin(x)/(x**4) - 945.0_dp*cos(x)/(x**5) + 945.0_dp*sin(x)/(x**6)
         end if
      case default
         stop "spherical_Bessel_jL():"
      end select
      !
      !! print *, "Bessel, L, x, test = ", L, x, test
      !
      !! ! Compare results
      !! if (abs((value - test)/value) > 1.0e-6_dp) then
      !!    print *, "spherical_Bessel_jL(): L, x, nag, direct = ", L, x, value, test
      !! end if
      !
   end function spherical_Bessel_jL_old
   !
   !
   function spherical_Legendre_Pmn(m,n,x)                  result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the spherical (associated) Legendre function P^m_n(x) as given in terms of the Gauss' 
   ! hypergeometric function F(a,b,c;x) for real parameters a, b, c and real values of x.
   !
   ! Call(s): hypergeometric_2F1, nag_s14aaf (from NAG-Libarary).
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: m, n, x
      real(kind=dp)             :: value
      !
      if (rabs_use_stop   .and.  x == 0) then
         stop "spherical_Legendre_Pmn(): program stop A."
      end if
      !
      if (rabs_use_naglib) then
         value = ((1+x)/(1-x))**(m*half) / nag_s14aaf(1-m,ifail)  * hypergeometric_2F1(-n,n+1,1-m,half*(one-x))
      else if (rabs_use_stop) then
         stop "spherical_Legendre_Pmn(): program stop B."
      end if
      !
   end function spherical_Legendre_Pmn
   !
   !
   subroutine transform_cart_to_cylinder(x,y,z,rho,phi,zc)
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram transforms the cartesian coordinates (x,y,z) into cylindrical coordinates (rho,phi,z). It terminates 
   ! with a proper message if r = zero.
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in)  :: x, y, z
      real(kind=dp), intent(out) :: rho, phi, zc
      !
      real(kind=dp)              :: xn
      intrinsic sqrt, atan
      !
      rho = sqrt( x*x + y*y )
      if (rabs_use_stop   .and.   rho <= zero) then
         stop "transform_cart_to_cylinder(): program stop A."
      end if
      !
      if (x == zero) then
         xn = 1.0e-99_dp
         print *, "transform_cart_to_cylinder(): x, xn = ", x, xn
      else
         xn = x
      end if
      phi   = atan(y/xn)
      !
      zc = z
      !
   end subroutine transform_cart_to_cylinder
   !
   !
   subroutine transform_cart_to_spherical(x,y,z,r,theta,phi)
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram transforms the cartesian coordinates (x,y,z) into spherical coordinates (r,theta,phi). It terminates 
   ! with a proper message if r = zero.
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in)  :: x, y, z
      real(kind=dp), intent(out) :: r,theta,phi
      !
      real(kind=dp)              :: xn
      intrinsic sqrt, atan, acos
      !
      r = sqrt( x*x + y*y + z*z )
      if (rabs_use_stop   .and.   r <= zero) then
         stop "transform_cart_to_spherical(): program stop A."
      end if
      !
      if (x == zero) then
         xn = 1.0e-99_dp
         print *, "xn = ",xn
      else
         xn = x
      end if
      phi   = atan(y/xn)
      !
      theta = acos( z/r )
      !
   end subroutine transform_cart_to_spherical
   !
   !
   subroutine transform_cylinder_to_cart(rho,phi,zc,x,y,z)
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram transforms the cylindrical coordinates (rho,phi,z) into cartesian coordinates (x,y,z). 
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in)  :: rho, phi, zc
      real(kind=dp), intent(out) :: x, y, z
      !
      intrinsic sin, cos
      !
      x = rho * cos(phi)
      y = rho * sin(phi)
      z = zc
      !
   end subroutine transform_cylinder_to_cart
   !
   !
   subroutine transform_spherical_to_cart(r,theta,phi,x,y,z)
   !--------------------------------------------------------------------------------------------------------------------
   ! This subprogram transforms the spherical coordinates (r,theta,phi) into cartesian coordinates (x,y,z). 
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in)  :: r,theta,phi
      real(kind=dp), intent(out) :: x, y, z
      !
      intrinsic sin, cos
      !
      x = r * sin(theta) * cos(phi)
      y = r * sin(theta) * sin(phi)
      z = r * cos(theta)
      !
   end subroutine transform_spherical_to_cart
   !
   !
   function Whittaker_M(a,b,x)                             result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Whittaker function M_a,b(x) for real parameters a and b and for real argument of x. 
   !
   ! Calls: hypergeometric_Phi().
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in)  :: a, b, x
      real(kind=dp)              :: value
      !
      value = x**(b+half) * exp(-x*half) * hypergeometric_Phi(b-a+half,b+b+1,x)
      !
   end function Whittaker_M
   !
   !
   function Whittaker_W(a,b,x)                             result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Whittaker function W_a,b(x) for real parameters a and b and for real argument of x. 
   !
   ! Calls: nag_s14aaf() [from NAG], Whittaker_M().
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in)  :: a, b, x
      real(kind=dp)              :: value
      !! real(kind=dp)           :: value, nag_s14aaf
      !
      print *, "WARNING: function Whittaker_W() has not been tested."
      print *, "WARNING: function Whittaker_W() has not been tested."
      !
      if (rabs_use_naglib) then
         value = nag_s14aaf(-b-b,ifail) / nag_s14aaf(half-b-a,ifail) * Whittaker_M(a,b,x)  +     &
                 nag_s14aaf( b+b,ifail) / nag_s14aaf(half+b-a,ifail) * Whittaker_M(a,-b,x)
      else if (rabs_use_stop) then
         stop "Whittaker_W(): program stop A."
      end if
      !
   end function Whittaker_W
   !
end module rabs_functions_math
