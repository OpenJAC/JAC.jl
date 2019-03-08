module rabs_angular
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module provides a set of procedures which support the analytic integration over angular coordinates or which are 
! related to this integration. This includes the Wigner n-j symbols but also a set of other coefficients for closed- and 
! open-shell SCF iterations and for the electron-electron interaction.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_dirac_orbital
   use rabs_functions_math
   use rabs_io_dialog
   use rabs_naglib
   implicit none
   private
   !
   public  :: angular_DF_Cabk      
                 ! Calculates the C(abk) angular coefficient of the direct Coulomb term in the closed-shell SCF iteration 
                 ! scheme.
   public  :: angular_DF_Dabk      
                 ! Calculates the D(abk) angular coefficient of the exchange Coulomb term in the closed-shell SCF 
                 ! iteration scheme.
   public  :: angular_E_beta_nu      
                 ! Calculates the E_beta^nu(abL) angular coefficient as defined by Grant (1988) for the Gaunt interaction 
                 ! and for other interactions.
   public  :: angular_nu      
                 ! Calculates the nu_nu,L angular coefficient of the Gaunt interaction as defined by Grant and Pyper (1976).
   public ::  bipolar_spherical_Ylm 
                 ! Returns the value of the bipolar spherical harmonic for angles theta1, phi1, theta2, phi2.
   public  :: cfp_coefficient             
                 ! Returns a coefficient of fractional parentage by looking up an appropriate table.
   private :: cfp_three_half             
                 ! Returns a coefficient of fractional parentage for equivalent electrons with j=3/2 by looking up an 
                 ! appropriate table.
   private :: cfp_five_half             
                 ! Returns a coefficient of fractional parentage for equivalent electrons with j=5/2 by looking up an 
                 ! appropriate table.
   private :: cfp_seven_half             
                 ! Returns a coefficient of fractional parentage for equivalent electrons with j=7/2 by looking up an 
                 ! appropriate table.
   private :: cfp_dummy             
                 ! Returns a coefficient of fractional parentage for 1 or 2 equivalent electrons with j > 7/2.
   public  :: CL_reduced_me             
                 ! Calculates the reduced matrix element of the C^L tensor  <kapa || C^(L) || kapb>.
   public  :: CL_reduced_me_mod             
                 ! Calculates the reduced matrix element of the C^L tensor  [kapa || C^(L) || kapb] due to 
                 ! Erikas&Gediminas (2009).
   public  :: CL_reduced_me_hfs             
                 ! Calculates the reduced matrix element of the C^L tensor  <kapa || C^(L) || kapb> due to the 
                 ! definition in HFS.
   public  :: Clebsch_Gordan            
                 ! Calculates the Clebsch-Gordan coefficient  <ja,ma,jb,mb;Jab,Mab>  by using the phase convention due
                 ! to Condon and Shortley.
   public  :: sigma_TtL_reduced_me
                 ! Calculates the reduced matrix element of the tensor  <kapa || sigma dot T_{tL} || kapb>.
   public  :: sigma_1mp_reduced_me
                 ! Calculates the (modified) reduced matrix element  [-kapa||sigma^1||kapb] due to Erikas&Gediminas (2009).   
   public  :: sigma_1pm_reduced_me
                 ! Calculates the (modified) reduced matrix element  [kapa||sigma^1||-kapb] due to Erikas&Gediminas (2009).   
   public  :: triangle             
                 ! Calculates the tringular factor Delta(ja,jb,jc).
   public  :: wigner_3j_symbol     
                 ! Calculates the value of a Wigner 3-j symbol for given quantum numbers.
   private :: wigner_6j_Delta      
                 ! Calculates the value of the Delta(a,b,c) symbol for given quantum numbers.
   public  :: wigner_6j_symbol     
                 ! Calculates the value of a Wigner 6-j symbol for given quantum numbers.
   public  :: wigner_6j_triangle
                 ! Calculate the triangular factors for 6j-symbol.
   public  :: wigner_9j_symbol     
                 ! Calculates the value of a Wigner 9-j symbol for given quantum numbers.
   public  :: wigner_9j_triangle
                 ! Calculate the triangular factors for 9j-symbol.
   public  :: wigner_eckardt_geometry     
                 ! Returns the "geometrical factor" from a complete matrix element owing to the Wigner-Eckardt theorem.
   public  :: wigner_rotation_D_matrix    
                 ! Returns the value of the Wigner rotation matrix D^j_pq(alpha, beta, gamma).
   public  :: wigner_rotation_d_small   
                 ! Returns the value of the Wigner rotation matrix d^j_pq(beta).
                 !
   integer, private :: ifail
   !
   ! Define some global variables to store the logarithmic Gamma function
   integer, parameter, private :: number_of_ln_gamma = 500
   logical, private            :: ln_gamma_not_initialized = .true.
   real(kind=dp), private      :: ln_gamma(number_of_ln_gamma)
   !
contains
   !
   function angular_DF_Cabk(kappa,kappa_p,k,q_p,equivalent)         result(Cabk)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the C(abk) angular coefficient of the direct Coulomb term in the DF energy functional within the 
   ! configuration-average approach. See I. Lindgren and A. Rosen (1974) for more details. Use  equivalent = .true. 
   ! if kappa and kappa_p arise from the same subshell and .false. otherwise. q_p is the occupation of the primed
   ! subshell (kappa_p). The parameter kappa is not needed within this procedure but is kept for similiarity with other 
   ! procedures.
   !
   ! Calls: angular_momentum_j().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: kappa, kappa_p, k, q_p
      logical, intent(in)  :: equivalent
      real(kind=dp)        :: Cabk
      !
      Cabk = kappa      ! kappa is not really needed for the calculation but kept for similiarity with other procedures
      if (k /= 0) then
         Cabk = zero
         return
      else
         Cabk = angular_momentum_j(kappa_p) + one
      end if
      !
      if (equivalent) then
         Cabk = Cabk * (q_p - one) / angular_momentum_j(kappa_p) 
      else
         Cabk = Cabk * q_p / (angular_momentum_j(kappa_p) + one)
      end if 
      !
   end function angular_DF_Cabk
   !
   !
   function angular_DF_Dabk(kappa,kappa_p,k,q_p,equivalent) result(Dabk)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the D(abk) angular coefficient of the exchange Coulomb term in the DF energy functional within the 
   ! configuration-average approach. See I. Lindgren and A. Rosen (1974) for more details. Use  equivalent = .true. if 
   ! kappa and kappa_p arise from the same subshell and .false. otherwise. q_p is the occupation of the primed
   ! subshell (kappa_p)
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: kappa, kappa_p, k, q_p
      logical, intent(in)  :: equivalent
      real(kind=dp)        :: Dabk
      !
      integer              :: l, l_p, j, j_p
      !
      l   = angular_momentum_l(kappa)
      l_p = angular_momentum_l(kappa_p)
      j   = angular_momentum_j(kappa)
      j_p = angular_momentum_j(kappa_p)
      !
      if (mod(l+l_p+k,2) == 1) then
         Dabk = zero
         return
      else
         Dabk = -(j_p + one) * wigner_3j_symbol(j,2*k,j_p,-1,0,1)**2
      end if
      !
      if (equivalent) then
         Dabk = Dabk * (q_p - one) / angular_momentum_j(kappa_p) 
      else
         Dabk = Dabk * q_p / (angular_momentum_j(kappa_p) + one)
      end if 
      !
   end function angular_DF_Dabk
   !
   !
   function angular_E_beta_nu(beta,nu,kappa,kappa_p,L)     result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the E_beta^nu (kappa,kappa',L) angular coefficient of the electron--electron interaction as defined by 
   ! I. P. Grant, in "Relativistic Atomic Structure Calculations" (1989), Eq.(66a-c).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: beta, kappa, kappa_p, L, nu
      real(kind=dp)        :: value
      !
      if (rabs_use_stop   .and.  beta /= -1   .and.   beta /= 1) then
         stop "angular_E_beta_nu(): program stop A."
      end if
      !
      if (nu == L-1   .and.   L > 0) then
         value = (nu + one - beta*(kappa_p - kappa)) / sqrt( (nu+one)*(nu+nu+one) )
      else if (nu == L   .and.   L >= 0) then
         value = - beta*(kappa + kappa_p) / sqrt( (nu+one)*nu )
      else if (nu == L+1   .and.   L >= 0) then
         value = (-nu - beta*(kappa_p - kappa)) / sqrt( nu*(nu+nu+one) )
      else if (rabs_use_stop) then
         stop "angular_E_beta_nu(): program stop B."
      end if
      !
   end function angular_E_beta_nu
   !
   !
   function angular_nu(nu,L)                               result(value)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the nu_nu,L angular coefficient of the Gaunt interaction as defined by I. P. Grant and N. C. Pyper, 
   ! J. Phy. B9 (1976) 761, Eq. (8).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: nu, L
      real(kind=dp)        :: value
      !
      if (nu == L-1   .and.   L > 0) then
         value = - (L+one) / (L+l+one)
      else if (nu == L   .and.   L >= 0) then
         value = one
      else if (nu == L+1   .and.   L >= 0) then
         value = - L / (L+L+one)
      else if (rabs_use_stop) then
         stop "angular_nu(): program stop A."
      end if
      !
   end function angular_nu
   !
   !
   function bipolar_spherical_Ylm(l1,l2,L,M,theta1,phi1,theta2,phi2)       result(bipolar)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the bipolar spherical harmonic   {Y_l1 (theta1,phi1) otimes Y_l2 (theta2,phi2)}_LM
   ! for given angles theta1, phi1, theta2, and phi2.
   !
   ! Calls: spherical_Ylm()
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: l1, l2, L, M
      real(kind=dp), intent(in) :: theta1, phi1, theta2, phi2
      real(kind=dp)             :: bipolar
      !
      integer                   :: m1, m2
      real(kind=dp)             :: wa
      !
      bipolar = zero
      !
      do  m1 = -l1,l1
         do  m2 = -l2,l2
	    wa = Clebsch_Gordan(l1+l1,m1+m1,l2+l2,m2+m2,L+L,M+M)
	    if (abs(wa) > eps10) then
	       bipolar = bipolar + wa * spherical_Ylm(l1,m1,theta1,phi1) * spherical_Ylm(l2,m2,theta2,phi2)
	    end if
	 end do
      end do
      !
   end function bipolar_spherical_Ylm
   !
   !
   subroutine cfp_coefficient(lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,coefp)
   !--------------------------------------------------------------------------------------------------------------------
   ! Selects a cfp coefficient from an appropriate table of fractional parentage  coefficients in jj-coupling. It has 
   ! been taken from the GRASP package and slightly adopted to the format of the RABS program.
   !
   ! Input variables:
   !
   !   lock     : + or - (2*j + 1).
   !   nel      : number of equivalent electrons in shell.
   !   ijd/ijp  : total j of daughter/parent state.
   !   ivd/ivp  : seniority of daughter/parent states.
   !   iwd/iwp  : other quantum number (if needed).
   !
   ! Output variable: 
   !
   !   coefp    : numerical result
   !
   ! This control routine does not check the input variables for consistency, except the trivial case of j = 1/2. 
   ! All other checks are performed at a lower level. The package will return correct results for j = 3/2, 5/2, 7/2. 
   ! Higher values of j return a value 1.0 if NEL = 1 or 2; otherwise 0 with an error signal. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: lock, nel, ijd, ivd, iwd, ijp, ivp, iwp
      real(kind=dp), intent(out):: coefp
      !
      integer :: k
      !
      k = abs(lock) / 2
      !
      select case(k)
      case(1)
         write (*,1)
         stop "cfp_coefficient(): program stop A."
      case(2)
         call cfp_three_half(nel,ijd,ijp,coefp)
      case(3)
         call cfp_five_half (nel,ijd,ivd,ijp,ivp,coefp)
      case(4)
         call cfp_seven_half(nel,ijd,ivd,ijp,ivp,coefp)
      case(5:)
         call cfp_dummy(lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,coefp)
      end select
      !
      1 format ("cfp_coefficient(): Unnecessary attempt to form a cfp for an electron with j = 1/2.")
      !
   end subroutine cfp_coefficient
   !
   !
   subroutine cfp_three_half(nel,ijd,ijp,coefp)
   !--------------------------------------------------------------------------------------------------------------------
   ! Table look-up for fractional parentage coefficients of  equivalent electrons with j = 3/2. See listing of CFP for 
   ! argument list.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: nel, ijd, ijp
      real(kind=dp), intent(out):: coefp
      !
      real(kind=dp), parameter :: c1 = 6.0_dp, c2 = 5.0_dp
      !
      if (nel <= 0   .or. nel > 4) goto 7
      !
      select case(nel)
      case(1)
         if (ijd /= 3   .or.   ijp /= 0) goto 7
         coefp = one
         return
      case(2)
         if (ijp /= 3) goto 7
         if (ijd == 0   .or.   ijd == 4) then
            coefp = one
            return
         end if
         goto 7
      case(3)
         if (ijd /= 3) goto 7
         if (ijp == 0) then
            coefp = sqrt(one/c1)
            return
         else if(ijp /= 4) then
            goto 7
         else
            coefp = -sqrt (c2/c1)
            return
         end if
      case(4)
         if (ijd /= 0   .or.   ijp /= 3) goto 7
         coefp = one
         return
      case(:0,5:)
         stop "cfp_three_half(): program stop A."
      end select
      !
    7 write (*,1) nel,ijd,ijp
      stop
      !
    1 format ("cfp3: error in trying to compute a cfp for a state with ",1i2," electrons with j = 3/2;",  &
            / " ijd = ",1i2,", ijp = ",1i2,".")
      !
   end subroutine cfp_three_half
   !
   !
   subroutine cfp_five_half(nel,ijd,ivd,ijp,ivp,coefp)
   !--------------------------------------------------------------------------------------------------------------------
   ! Table look-up for fractional parentage coefficients of  equivalent electrons with j = 3/2. See listing of CFP for 
   ! argument list.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: nel, ijd, ivd, ijp, ivp
      real(kind=dp), intent(out):: coefp
      !
      integer       :: ij1, iv1, ij2, iv2, is, k, kd, kp, n
      real(kind=dp) :: denom, dnel, fact
      !
      ! Set 1.0 tables of data
      !
      integer, dimension(3)    :: norm = (/  18,7,14 /)
      integer, dimension(3,3)  :: ij   = reshape(source=(/  5,0,5,0,4,3,0,8,9 /),   shape=(/3,3/)), &
                                  iv   = reshape(source=(/  1,0,1,8,2,3,8,2,3 /),   shape=(/3,3/)), &
                                  num3 = reshape(source=(/ -4,0,0,5,-5,3,9,2,-11/), shape=(/3,3/))
      !
      real(kind=dp), parameter :: d0 = 0.0_dp, d1 = 1.0_dp, d7 = 7.0_dp
      !
      ! 2.0 Locate entry in cfp table.
      if (nel <= 0) goto 12
      if (nel >= 4) goto 1
      !
      n = nel
      ij1 = ijd;   iv1 = ivd;   ij2 = ijp;   iv2 = ivp;   goto 2
      !
    1 if (nel > 6) goto 12
      n = 7 - nel
      ij1 = ijp;   iv1 = ivp;   ij2 = ijd;   iv2 = ivd
      !
      ! 2.1 Find 'daughter' index.
    2 k = 0
    3 k = k+1
      if (k > 3) then
         print *, "Find 'daughter' index.n = ",n
         goto 12
      end if
      if (ij(n,k) /= ij1) goto 3
      if (iv(n,k) /= iv1) goto 3
      kd = k
      !
      ! 2.2 Find 'parent' index.
      if (n /= 1) goto 4
      if (iv2 /= 0) goto 12
      if (ij2 == 0) goto 6
      goto 12
    4 k = 0
    5 k = k+1
      if (k > 3) then
         print *, "Find 'parent' index."
         goto 12
      end if
      if (ij(n-1,k) /= ij2) goto 5
      if (iv(n-1,k) /= iv2) goto 5
      kp = k
      !
      ! 3.0 compute coefficients.
      ! 3.1 table look-up
      !!x goto (6,6,7), n
      select case (n)
      case(1,2)
        goto 6
      case(3)
        goto 7
      end select
      !
    6 coefp = d1;   goto 10
    7 continue
      coefp = dble (num3(kd,kp))
      denom = dble (norm(kd))
      if (coefp < zero) then  ! Arithmetic if (coefp) 9,11,8 is obsolet
         goto 9
      else if (coefp == zero) then
         goto 11
      else
         goto 8
      end if
    8 continue
      coefp = sqrt (coefp/denom)
      goto 10
    9 continue
      coefp = -sqrt(-coefp/denom)
      !
      ! 3.2 insert additional factors for hole states
   10 if (nel <= 3) goto 11
      dnel  = nel
      fact  = ((d7-dnel)/dnel) * (d1+dble(ijp)) / (d1+dble(ijd))
      coefp = coefp * sqrt(fact)
      is    = abs((ijd-ijp-ivd+ivp)/2)
      if (mod(is,2) == 0) goto 11
      coefp = -coefp
   11 return
      !
      ! 4.0 fault mode section.
   12 coefp = d0
      write (*,15) nel, ijd, ivd, ijp, ivp
      !x stop
      !
   15 format ("cfp_five_half(): Error in trying to compute a cfp for a state with ",1i2," electrons with j = 5/2;",  &
            / " ijd = ",1i2,", ivd = ",1i2,", ijp = ",1i2,", ivp = ",1i2,".")
      !
   end subroutine cfp_five_half
   !
   !
   subroutine cfp_seven_half(nel,ijd,ivd,ijp,ivp,coefp)
   !--------------------------------------------------------------------------------------------------------------------
   ! Table look-up for fractional parentage coefficients of  equivalent electrons with j = 7/2. See listing of CFP for 
   ! argument list.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: nel, ijd, ivd, ijp, ivp
      real(kind=dp), intent(out):: coefp
      !
      integer                   :: ij1, iv1, ij2, iv2, is, k, kd, kp, n
      real(kind=dp)             :: denom, dnel, fact
      !
      ! Set 1.0 tables of data
      !
      integer, dimension(4,8)  ::   ij   = reshape(source=(/  7,0,7,0,0,4,3,4,0,8,5,8,0,12,9,12,0,0,   &
                                                             11,4,0,0,15,8,0,0,0,10,0,0,0,16  /),   shape=(/4,8/)), &
                                    iv   = reshape(source=(/  1,0,1,0,10,2,3,2,10,2,3,2,10,2,3,2,10,   &
                                                             10,3,4,10,10,3,4,10,10,10,4,10,10,10,4 /), shape=(/4,8/))
      integer, dimension(6,4)  ::   num3 = reshape(source=(/  9,0,0,0,0,0,-5,3,121,143,-55,0,-9,-11,   &
                                                             12,-900,39,5,-13,0,-65,343,104,-17  /),  shape=(/6,4/))
      integer, dimension(6)    ::  norm3 = (/  36,14,198,1386,198,22  /) 
      integer, dimension(8,6)  ::  num4  = reshape(source=(/  1,280,308,1144,0,0,0,0,0,54,-121,0,-968, &
                                                              169,462,0,0,-231,-14,195,-77,2366,-343,0,0,-65,250,     &
                                                             -245,-1755,90,-945,140,0,-210,91,624,280,2275,650,      &
                                                              234,0,0,140,-1224,0,-560,680,627 /),  shape=(/8,6/))
      integer, dimension(8)    ::  norm4 = (/1,840,924,3432,3080,5460,3080,1001/) 
      !
      real(kind=dp), parameter :: d0 = 0.0_dp, d1 = 1.0_dp, d9 = 9.0_dp
      !
      ! 2.0 locate entry in cfp table.
      if (nel <= 0) goto 14
      if (nel >= 5) goto 1
      !
      n = nel
      ij1 = ijd;   iv1 = ivd;   ij2 = ijp;   iv2 = ivp;   goto 2
    1 if (nel > 8) goto 14
      n = 9-nel
      ij1 = ijp;   iv1 = ivp;   ij2 = ijd;   iv2 = ivd
      !
      ! 2.1 find 'daughter' index.
    2 k = 0
    3 k = k+1
      if (k > 8) goto 14
      if (ij(n,k) /= ij1) goto 3
      if (iv(n,k) /= iv1) goto 3
      kd = k
      !
      ! 2.2 find 'parent' index.
      if (n /= 1) goto 4
      if (iv2 /= 0) goto 14
      if (ij2 == 0) goto 6
      goto 14
    4 k = 0
    5 k = k+1
      if (k > 8) goto 14
      if (ij(n-1,k) /= ij2) goto 5
      if (iv(n-1,k) /= iv2) goto 5
      kp = k
      !
      ! 3.0 compute coefficients.
      ! 3.1 table look-up
      !!x goto (6,6,7,11), n
      select case (n)
      case(1,2)
        goto 6
      case(3)
        goto 7
      case(4)
         goto 11
      end select
      !
    6 coefp = d1;   goto 12
    7 continue
      coefp = dble (num3(kd,kp));   denom = dble (norm3(kd))
    8 if (coefp < zero) then  ! Arithmetic if (coefp) 10,13,9 is obsolet
         goto 10
      else if (coefp == zero) then
         goto 13
      else
         goto 9
      end if
    9 continue
      coefp = sqrt(coefp/denom);     goto 12
   10 continue
      coefp = -sqrt(-coefp/denom);   goto 12
   11 continue
      coefp = dble (num4(kd,kp));    denom = dble (norm4(kd));   goto 8
      !
      ! 3.2 insert additional factors for hole states
   12 if (nel <= 4) goto 13
      dnel  = dble (nel)
      fact  = ((d9-dnel)/dnel) * (d1+dble(ijp)) / (d1+dble(ijd))
      coefp = coefp * sqrt(fact)
      is = abs((ijd-ijp-ivd+ivp)/2-3)
      if (mod(is,2) == 0) goto 13
      coefp = -coefp
   13 return
      !
      ! 4.0 fault mode section.
   14 write (*,15) nel, ijd, ivd, ijp, ivp
      stop
      !
   15 format ("cfp7: error in trying to compute a cfp for a state with ',1i2,' electrons with j = 7/2;",   &
            / " ijd = ",1i2,", ivd = ",1i2,  ", ijp = ",1i2,", ivp = ",1i2,".")
      !
   end subroutine cfp_seven_half
   !
   !
   subroutine cfp_dummy(lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,coefp)
   !--------------------------------------------------------------------------------------------------------------------
   ! This is a dummy subroutine. It returns correct values for 1 or 2 particle or single hole states, and signals an 
   ! error otherwise.   
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: lock, nel, ijd, ivd, iwd, ijp, ivp, iwp
      real(kind=dp), intent(out):: coefp
      !
      integer :: locj
      !
      if (nel == 1   .or.   nel == 2   .or.   nel == abs(lock)) then
         coefp = one
      else
         if (ijd == 99    .and.   ijp == 99   .and. &
             ivd == iwd   .and. ivp == iwp)   then
            ! Just to avoid warnings because of unused arguments 
            print *, " " 
         end if
         !
         locj = abs(lock) - 1
         write (*,1) locj
         stop
      end if
      !
      1 format ("cfp_dummy(): Inadmissable attempt to obtain a cfp for a state of a shell with j = ",1i2,"/2;",  &
              / " new subprogram required.")
      !
   end subroutine cfp_dummy
   !
   !
   function CL_reduced_me(kapa,L,kapb)                     result(redme)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the reduced matrix element of the C^L spherical tensor  <kapa || C^(L) || kapb>.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kapa, kapb, L
      real(kind=dp)       :: redme
      !
      integer :: ja, jb, la, lb
      !
      redme = zero
      !
      ja = angular_momentum_j(kapa);   jb = angular_momentum_j(kapb)
      la = angular_momentum_l(kapa);   lb = angular_momentum_l(kapb)
      !
      if (rabs_use_stop   .and.   mod(ja+1,2) /= 0) then
         stop "CL_reduced_me(): program stop A."
      end if
      !
      redme = ((-1)**((ja+1)/2)) * sqrt( (ja+one)*(jb+one) ) * wigner_3j_symbol(ja,L+L,jb,1,0,-1)
      !
   end function CL_reduced_me
   !
   !
   function CL_reduced_me_mod(kapa,L,kapb)                  result(redme)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the modified form [kapa || C^(L) || kapb] of the reduced matrix element of the C^L spherical tensor
   ! <kapa || C^(L) || kapb> due to Erikas&Gediminas (2009).
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kapa, kapb, L
      real(kind=dp)       :: redme
      !
      integer             :: ja, jb, la, lb
      !
      redme = zero
      !
      ja = angular_momentum_j(kapa);   jb = angular_momentum_j(kapb)
      la = angular_momentum_l(kapa);   lb = angular_momentum_l(kapb)
      !
      if (rabs_use_stop   .and.   mod(ja+1,2) /= 0) then
         stop "CL_reduced_me_mod(): program stop A."
      end if
      !
      ! Implement the parity selection
      if (mod(la+L+lb,2) /= 0) then
         print *, "CL_reduced_me_mod(): parity rule requests pi(...) = 0."
         return
      end if
      !
      redme = ((-1)**((ja+1)/2)) * sqrt( (jb+one) ) * wigner_3j_symbol(ja,L+L,jb,1,0,-1)
      !
   end function CL_reduced_me_mod
   !
   !
   function CL_reduced_me_hfs(kapa,L,kapb)                     result(redme)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the reduced matrix element of the C^L spherical tensor  <kapa || C^(L) || kapb>   due to the definition 
   ! in the HFS component.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kapa, kapb, L
      real(kind=dp)       :: redme
      !
      integer :: ja, jb, la, lb
      !
      redme = zero
      !
      ja = angular_momentum_j(kapa);   jb = angular_momentum_j(kapb)
      la = angular_momentum_l(kapa);   lb = angular_momentum_l(kapb)
      !
      if (rabs_use_stop   .and.   mod(jb+L+L+1,2) /= 0) then
         stop "CL_reduced_me_hfs(): program stop A."
      end if
      !
      if (mod(la+L+lb,2) == 0) then
         redme = ((-1)**((jb-L-L-1)/2)) * sqrt( (jb+one) ) * wigner_3j_symbol(ja,jb,L+L,1,-1,0)
      else
         redme = zero
      end if
      !
   end function CL_reduced_me_hfs
   !
   !
   function Clebsch_Gordan(ja,ma,jb,mb,Jab,Mab)               result(CG)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the Clebsch-Gordan coefficient   <ja',ma',jb',mb';Jab',Mab'>   by using the phase convention due to
   ! Condon and Shortley; the dummy arguments of this procedure are the double of the actual angular momentum quantum 
   ! numbers, i.e. ja = ja' + ja', ma = ma' + ma', ...
   !
   ! Calls: wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: ja, jb, Jab, ma, mb, Mab
      real(kind=dp)       :: CG
      !
      CG = zero
      !
      if (rabs_use_stop   .and.   mod(ja-jb+Mab,2) /= 0) then
         stop "Clebsch_Gordan(): program stop A."
      end if
      !
      CG = ((-1)**((ja-jb+Mab)/2)) * sqrt( (Jab+one) ) * wigner_3j_symbol(ja,jb,Jab,ma,mb,-Mab)
      !
   end function Clebsch_Gordan
   !
   !
   function sigma_TtL_reduced_me(kapa,L,t,kapb)              result(redme)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the reduced matrix element of the C^L spherical tensor  <kapa || sigma dot T_{tL} || kapb>.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kapa, kapb, L, t
      real(kind=dp)       :: redme
      !
      redme = zero
      !
      if (L == t + 1) then
         redme = sqrt( (t+one)/(two*two*pi) ) * (one + (kapa+kapb)/(t+one)) * CL_reduced_me(-kapa,t,kapb) 
      else if (L == t) then
         redme = sqrt( (t+t+one)/(two*two*pi*t*(t+one)) ) * (kapb-kapa) * CL_reduced_me(kapa,t,kapb) 
      else if (L == t - 1) then
         redme = sqrt( t/(two*two*pi) ) * (-one + (kapa+kapb)/t ) * CL_reduced_me(-kapa,t,kapb) 
      else
         stop "sigma_TtL_reduced_me(): program stop A."
      end if
      !
   end function sigma_TtL_reduced_me
   !
   !
   function sigma_1mp_reduced_me(kapa,kapb)                 result(redme)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the (modified) reduced matrix element of the sigma^1 spherical tensor [-kapa||sigma^1||kapb] due to 
   ! Erikas&Gediminas (2009). Note that kapa must be given without the minus sign.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kapa, kapb
      real(kind=dp)       :: redme
      !
      integer :: ja, jb, la, lb
      !
      redme = zero
      !
      ja = angular_momentum_j(kapa);   jb = angular_momentum_j(kapb)
      la = angular_momentum_l(kapa);   lb = angular_momentum_l(kapb)
      !
      redme = (  Clebsch_Gordan(lb+lb,0,1,1,ja,1)  * Clebsch_Gordan(lb+lb,0,1,1,jb,1)                                  & 
	       - Clebsch_Gordan(lb+lb,2,1,-1,ja,1) * Clebsch_Gordan(lb+lb,2,1,-1,jb,1) ) / Clebsch_Gordan(jb,1,2,0,ja,1)
      !
   end function sigma_1mp_reduced_me
   !
   !
   function sigma_1pm_reduced_me(kapa,kapb)                 result(redme)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the (modified) reduced matrix element of the sigma^1 spherical tensor [kapa||sigma^1||-kapb] due to 
   ! Erikas&Gediminas (2009). Note that kapb must be given without the minus sign.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kapa, kapb
      real(kind=dp)       :: redme
      !
      integer :: ja, jb, la, lb
      !
      redme = zero
      !
      ja = angular_momentum_j(kapa);   jb = angular_momentum_j(kapb)
      la = angular_momentum_l(kapa);   lb = angular_momentum_l(kapb)
      !
      redme = (  Clebsch_Gordan(la+la,0,1,1,ja,1)  * Clebsch_Gordan(la+la,0,1,1,jb,1)                                  & 
	       - Clebsch_Gordan(la+la,2,1,-1,ja,1) * Clebsch_Gordan(la+la,2,1,-1,jb,1) )  / Clebsch_Gordan(jb,1,2,0,ja,1)
      !
   end function sigma_1pm_reduced_me
   !
   !
   function triangle(i2a,i2b,i2c)                          result(Delta)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the tringular factor Delta(ja,jb,jc). The arguments in this integer function are i2a = 2*ja+1, ... 
   ! The result is 0 if the triangular condition failes and 1 otherwise. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: i2a, i2b, i2c
      integer              :: Delta, i
      !
      i = i2b - i2c
      if (i2a >= abs(i) + 1  .and.   i2a <= i2b + i2c - 1) then
         Delta = 1
      else
         Delta = 0
      end if
      !
   end function triangle
   !
   !
   function wigner_3j_symbol(ja,jb,jc,ma,mb,mc)              result(w3j)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the Wigner 3-j symbol by its algebraic formulae as displayed in many texts on the theory of angular 
   ! momentum (see R. D. Cowan, The Theory of Atomic Structure and Spectra;  University of California Press, 1981, p. 142).
   ! The integer arguments ja, ... of this function must be the double of the corresponding quantum numbers in the 3-j 
   ! symbol, i.e. jk = jk' + jk', mk = mk' + mk' in  
   !
   !                       (  ja'   jb'   jc'  ) 
   !                       (                   ) 
   !                       (  ma'   mb'   mc'  )
   !
   ! Calls:  nag_s14aaf() [from NAG], triangle().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)      :: ja, jb, jc, ma, mb, mc
      real(kind=dp)            :: w3j
      !
      integer                  :: i, k, kmin, kmax
      integer, dimension(1:15) :: ik 
      real(kind=qp)            :: delta, qsum, sumk
      !
      ! Test the triangular condition and that for magnetic quantum numbers
      !
      if (ma+mb+mc /= 0) then
         w3j = zero
	 return
      else if (triangle(ja+1,jb+1,jc+1) == 0) then
         w3j = zero
	 return
      else if (abs(ma) > ja  .or.  abs(mb) > jb  .or.  abs(mc) > jc) then
         w3j = zero
	 return
      else if (rabs_use_stop                 .and. &
              (mod(ma+ja+ja,2) /= mod(ja,2)  .or.  mod(mb+jb+jb,2) /= mod(jb,2)  .or.  mod(mc+jc+jc,2) /= mod(jc,2)) ) then
         stop "wigner_3j_symbol(): program stop A."
      endif
      !
      !
      if (ja > 75 .or. jb > 75 .or. jc > 75)  then
         call set_warning("modified wigner_3j_symbol(), November 2010; arbitrary value")
         print *, "modified wigner_3j_symbol(), November 2010"
         w3j = one  ! arbitrary value
	 return
      end if
      !
      !
      ik(1)  =  ja + jb - jc;   ik(2)  =  ja - jb + jc
      ik(3)  = -ja + jb + jc;   ik(4)  =  ja + jb + jc + two
      !
      ik(5)  =  ja - ma;   ik(6)  =  ja + ma;   ik(7)  =  jb - mb
      ik(8)  =  jb + mb;   ik(9)  =  jc - mc;   ik(10) =  jc + mc
      !
      ik(11) =  jb - jc - ma;   ik(12) =  ja - jc + mb
      ik(13) =  jc - jb + ma;   ik(14) =  ja - jb - mc
      !
      do  i = 1,14
	 if (rabs_use_stop   .and.   mod(ik(i),2) == 1) then
            stop "wigner_3j_symbol(): program stop B."
	 end if
	 ik(i) = ik(i) / 2
      end do
      !
      ! Calculate the 3-j delta factor
      if (rabs_use_naglib)  then
         delta =   1.0_dp * nag_s14aaf(ik(1)+one,ifail) * nag_s14aaf(ik(2)+one,ifail) &
                          * nag_s14aaf(ik(3)+one,ifail) * nag_s14aaf(ik(5)+one,ifail) &
                          * nag_s14aaf(ik(6)+one,ifail) * nag_s14aaf(ik(7)+one,ifail) &
                          * nag_s14aaf(ik(8)+one,ifail) * nag_s14aaf(ik(9)+one,ifail) &
                          * nag_s14aaf(ik(10)+one,ifail) / nag_s14aaf(ik(4)+one,ifail) 
      else
         delta =   1.0_dp * factorial_dp(ik(1))  * factorial_dp(ik(2)) * factorial_dp(ik(3))  * factorial_dp(ik(5)) &
                          * factorial_dp(ik(6))  * factorial_dp(ik(7)) * factorial_dp(ik(8))  * factorial_dp(ik(9)) &
                          * factorial_dp(ik(10)) / factorial_dp(ik(4)) 
      end if
      !
      ! Find out the intervall of summation  k  and sum up
      kmin = max(0,ik(11),ik(12))
      kmax = min(ik(1),ik(5),ik(8))
      !
      qsum = 0.0_dp
      do  k = kmin,kmax
         if (rabs_use_naglib)  then
	    sumk = 1.0_dp * &
	           nag_s14aaf(k+one,ifail)        * nag_s14aaf(ik(1)-k+one,ifail) &
	         * nag_s14aaf(ik(5)-k+one,ifail)  * nag_s14aaf(ik(8)-k+one,ifail) &
	         * nag_s14aaf(ik(13)+k+one,ifail) * nag_s14aaf(k-ik(12)+one,ifail) 
         else
	    sumk = 1.0_dp * factorial_dp(k)        * factorial_dp(ik(1)-k)  * factorial_dp(ik(5)-k)  &
                          * factorial_dp(ik(8)-k)  * factorial_dp(ik(13)+k) * factorial_dp(k-ik(12)) 
         end if
	 sumk = 1.0_dp / sumk 
	 if (mod(k,2) == 0) then
	    qsum = qsum + sumk
	 else
	    qsum = qsum - sumk
	 end if
      end do
      !
      if (mod(ik(14),2) /= 0 ) then
	 w3j = - sqrt(delta) * qsum
      else
	 w3j = sqrt(delta) * qsum
      end if
      !
      !! print *, "ja,jb,jc,ma,mb,mc,w3j = ",ja,"/2 ",jb,"/2 ",jc,"/2 ", &
      !!                                     ma,"/2 ",mb,"/2 ",mc,"/2 ",w3j
      !
   end function wigner_3j_symbol
   !
   !
   function wigner_6j_Delta(a,b,c)                         result(Delta)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the Delta factor of three angular momenta a', b', c' as given in R. D. Cowan, "The Theory of Atomic 
   ! Structure and Spectra", University of California Press, Berkeley, Los Angeles, London, 1981, Eq. (5.24). It is 
   ! assumed that a, b, c represent the double of angular momenta (i.e. a = a' + a', ...) and that they fulfill the
   ! triangular condition triangle(a+1,b+1,c+1) = 1. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: a, b, c
      real(kind=dp)       :: Delta
      !
      if (rabs_use_stop) then
      if (mod(a+b-c,2)  /= 0  .or.  mod(a-b+c,2)   /= 0    .or.  mod(-a+b+c,2) /= 0  .or.  mod(a+b+c+2,2) /= 0  ) then
         print *, "a, b, c = ", a, b, c
         stop "wigner_6j_Delta(): program stop A."
      else if ((a+b-c)/2 < 0  .or.  (a-b+c)/2 < 0  .or.  (-a+b+c)/2 < 0  .or.  (a+b+c+2)/2 < 0) then
         print *, "a, b, c = ", a, b, c
         stop "wigner_6j_Delta(): program stop B."
      end if
      end if
      !
      if (rabs_use_naglib)  then
         Delta = nag_s14aaf((a+b-c)/2+one,ifail)  * nag_s14aaf((a-b+c)/2+one,ifail) * &
                 nag_s14aaf((-a+b+c)/2+one,ifail) / nag_s14aaf((a+b+c+2)/2+one,ifail)
      else
         Delta = factorial((a+b-c)/2)  * factorial((a-b+c)/2) * factorial((-a+b+c)/2) / factorial((a+b+c+2)/2)
      end if
      Delta = sqrt(Delta)
      !
   end function wigner_6j_Delta
   !
   !
   function wigner_6j_symbol(j1,j2,j3,l1,l2,l3,no_triangles) result(w6j)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the numerical value of a 6j-symbol (given in w6j) using the algebraic formulae in R.D. Cowan, 
   ! "The Theory of Atomic Structure and Spectra", University of California Press, Berkeley, Los Angeles, London, 
   ! 1981, Eq. (5.23). A straightforward numerical evaluation is performed. The integer arguments ja, ... of this 
   ! function must be the double of the corresponding quantum numbers in the 6-j symbol, i.e. jk = jk' + jk' > in  
   !
   !                       (  j1'   j2'   j3'  ) 
   !                      {(                   )} 
   !                       (  l1'   l2'   l3'  )
   !
   ! Calls:  triangle().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)           :: j1, j2, j3, l1, l2, l3
      logical, intent(in), optional :: no_triangles
      real(kind=dp)                 :: w6j, x
      !  
      integer                       :: n1, n2, n3, n4, n5, n6, n7, k, kmin, kmax, icount, ki
      !
      if (ln_gamma_not_initialized) then
         ln_gamma(1) = one;   ln_gamma(2) = one;   x = two
         do k = 3, 30
            ln_gamma(k) = ln_gamma(k-1) * x;   x = x + one
         end do
         !
         do k = 1, 30
            ln_gamma(k) = log(ln_gamma(k))
         end do
         x = three * ten
         do k = 31, number_of_ln_gamma
            ln_gamma(k) = ln_gamma(k-1) + log(x);   x = x + one
         end do
         ln_gamma_not_initialized = .false.
      end if
      !
      w6j = zero
      !
      if (present(no_triangles)) then
         if (no_triangles) goto 1
      end if
      if (triangle(j1+1,j2+1,j3+1) == 0  .or.  triangle(j1+1,l2+1,l3+1) == 0  .or. triangle(l1+1,j2+1,l3+1) == 0  .or. &
          triangle(l1+1,l2+1,j3+1) == 0     ) then
         return
      else if (rabs_use_stop			   .and. &
           (j1 < 0  .or.  j2 < 0  .or.  j3 < 0  .or.  l1 < 0  .or.  l2 < 0  .or.  l3 < 0)) then
         print *, "j1, j2, j3, l1 ,l2, l3 = ",j1, j2, j3, l1 ,l2, l3
         stop  "wigner_6j_symbol(): program stop A."
      end if
    1 continue  
      !  
      ! Evaluate the algebraic formulae
      !
      n1 = (j1 + j2 + j3)/2;         n2 = (l2 + l1 + j3)/2;       n3 = (j1 + l2 + l3)/2;       n4 = (j2 + l1 + l3)/2
      n5 = (j1 + j2 + l2 + l1)/2;    n6 = (j1 + l1 + j3 + l3)/2;  n7 = (j2 + l2 + j3 + l3)/2
      kmin = max(n1,n2,n3,n4) + 1;   kmax = min(n5,n6,n7) + 1
      !
      w6j = one;   icount = 0
      if (kmin /= kmax) then
         do k = kmin+1, kmax
            ki = kmax - icount
	    w6j = one -(w6j * ki*(n5-ki+two)*(n6-ki+two)*(n7-ki+two))/         &
                  ((ki-one-n1)*(ki-one-n2)*(ki-one-n3)*(ki-one-n4))
            icount = icount + 1
         end do
      end if
      w6j = w6j*exp(                                                           &
        ( ln_gamma(kmin+1)   -ln_gamma(kmin-n1)-ln_gamma(kmin-n2)   -          &
	  ln_gamma(kmin-n3)  -ln_gamma(kmin-n4)-ln_gamma(n5+2-kmin) -          &
	  ln_gamma(n6+2-kmin)-ln_gamma(n7+2-kmin)) +                           &
	((ln_gamma(n1+1-j1)+ln_gamma(n1+1-j2)+ln_gamma(n1+1-j3)-ln_gamma(n1+2)+&
	  ln_gamma(n2+1-l2)+ln_gamma(n2+1-l1)+ln_gamma(n2+1-j3)-ln_gamma(n2+2)+&
	  ln_gamma(n3+1-j1)+ln_gamma(n3+1-l2)+ln_gamma(n3+1-l3)-ln_gamma(n3+2)+&
          ln_gamma(n4+1-j2)+ln_gamma(n4+1-l1)+ln_gamma(n4+1-l3)-ln_gamma(n4+2))&
	  /two))
      if(mod(n5+kmin,2) == 0) w6j=-w6j;   if(mod(j1+j2+l2+l1,4) /= 0) w6j=-w6j
      !
   end function wigner_6j_symbol
   !
   !
   function wigner_6j_triangle(j1,j2,j3,l1,l2,l3)             result(Delta)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculate the triangular factors for 6j-symbol. The arguments in this integer function are j1 = 2*ja, ... 
   ! The result is 0 if the triangular condition failes and 1 itherwise.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)      :: j1, j2, j3, l1, l2, l3
      integer                  :: Delta
      !
      Delta = 0
      if (triangle(j1+1,j2+1,j3+1) == 0) return   
      if (triangle(j1+1,l2+1,l3+1) == 0) return
      if (triangle(l1+1,j2+1,l3+1) == 0) return    
      if (triangle(l1+1,l2+1,j3+1) == 0) return
      Delta = 1
      !
   end function wigner_6j_triangle 
   !
   !
   function wigner_9j_symbol(j11,j12,j13,j21,j22,j23,j31,j32,j33,no_triangels)        result(w9j)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the numerical value of a Wigner 9-j symbol by using the expansion and summation over a triple product
   ! of 6j-symbols as listed, for instance, in R.D. Cowan, "The Theory of Atomic Structure and Spectra", University of 
   ! California Press, Berkeley, Los Angeles, London, 1981, Eq. (5.23). A straightforward numerical evaluation is 
   ! performed. The integer arguments ja, ... of this function must be the double of the corresponding quantum numbers 
   ! in the 9-j symbol, i.e. jk = jk' + jk > 0 in  
   !
   !                       (  j11'   j12'   j13'  ) 
   !                       (  j21'   j22'   j23'  ) 
   !                       (  j31'   j32'   j33'  ) 
   !
   ! Calls:  triangle(), wigner_6j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)            :: j11, j12, j13, j21, j22, j23, j31, j32, j33
      logical, intent(in), optional  :: no_triangels
      real(kind=dp)                  :: w9j
      !
      integer :: j, jmin, jmax
      !
      w9j = zero
      !
      if (present(no_triangels)  .and. no_triangels) goto 1
      if (triangle(j11+1,j12+1,j13+1) == 0  .or. triangle(j21+1,j22+1,j23+1) == 0 .or. triangle(j31+1,j32+1,j33+1) == 0 .or.&
          triangle(j11+1,j21+1,j31+1) == 0  .or. triangle(j12+1,j22+1,j32+1) == 0 .or. triangle(j13+1,j23+1,j33+1) == 0) then
         return
      end if
    1 continue
      !  
      ! Evaluate the algebraic expansion formulae
      !
      jmin = max(abs(j11-j33),abs(j12-j23),abs(j32-j21))
      jmax = min(j11+j33,j12+j23,j32+j21)
      !
      if (jmin <= jmax) then
         do  j = jmin,jmax,2
            w9j = w9j + (-1)**j * (j+one) * wigner_6j_symbol(j11,j21,j31,j32,j33,j) *  &
                                            wigner_6j_symbol(j12,j22,j32,j21,j,j23) *  &
                                            wigner_6j_symbol(j13,j23,j33,j,j11,j12)
         end do
      end if
      !
   end function wigner_9j_symbol
   !
   !
   function wigner_9j_triangle(j11,j12,j13,j21,j22,j23,j31,j32,j33)      result(Delta)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the triangular factor for a Wigner 9-j symbol. The integer arguments ja, ... of this function must be 
   ! the double of the corresponding quantum numbers in the 9-j symbol, i.e. jk = jk' + jk > 0 in  
   !
   !                       (  j11'   j12'   j13'  ) 
   !                       (  j21'   j22'   j23'  ) 
   !                       (  j31'   j32'   j33'  ) 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
      integer             :: Delta
      !
      Delta = 0
      if (triangle(j11+1,j12+1,j13+1) == 0) return
      if (triangle(j21+1,j22+1,j23+1) == 0) return
      if (triangle(j31+1,j32+1,j33+1) == 0) return
      if (triangle(j11+1,j21+1,j31+1) == 0) return
      if (triangle(j12+1,j22+1,j32+1) == 0) return
      if (triangle(j13+1,j23+1,j33+1) == 0) return
      Delta = 1
      !
   end function wigner_9j_triangle
   !
   !
   function wigner_eckardt_geometry(totalJ_f,MM_f,J_t,MM_t,totalJ_i,MM_i)       result(geom)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the "geometrical factor" of the Wigner-Eckardt theorem by splitting a (complete) matrix element in a 
   ! "reduced" matrix element and a geometrical cofactor. Here, the definition of Grant, J. Phys. B7 (1974), 1458 is used.
   !
   !         <complete ME>  =  geom  *  <reduced ME> 
   !
   ! The given quantum numbers totalJf, MM_f .... are the doubles of the original ones which characterize the complete 
   ! matrix element.
   !
   ! Calls: is_triangle(), wigner_3j_symbol().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: totalJ_f, MM_f, J_t, MM_t, totalJ_i, MM_i
      !
      integer              :: phase
      real(kind=dp)        :: geom, gfactor
      ! 
      geom = zero
      if (.not.is_triangle(totalJ_f,J_t,totalJ_i)   .or.  MM_t + MM_i - MM_f /= 0) then
         return
      end if
      !
      gfactor = sqrt(totalJ_f + one)
      phase   = totalJ_f - MM_f + 2 + 32
      if (mod(phase,4) == 0) then
         gfactor = - gfactor
      else if (rabs_use_stop       .and. (mod(phase,4) == 1   .or.   mod(phase,4) == 3)) then
         stop "Wigner_Eckardt_geometry(): program stop A."
      end if
      geom = gfactor * wigner_3j_symbol(totalJ_f,J_t,totalJ_i,-MM_f,MM_t,MM_i)
      !
   end function wigner_eckardt_geometry
   !
   !
   function wigner_rotation_D_matrix(j,p,q,alpha,beta,gamma)                      result(wa)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Wigner rotation matrix D^j_pq(alpha, beta, gamma).
   !
   ! Calls: wigner_rotation_d_small().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: j, p, q
      real(kind=dp), intent(in) :: alpha, beta, gamma
      complex(kind=dp)          :: wa
      !
      if (abs(p) > j  .or.  abs(q) > j) then
         stop "wigner_rotation_D_matrix(): program stop A."
      else if (alpha < zero  .or.  alpha > two*pi) then
         stop "wigner_rotation_D_matrix(): program stop B."
      else if (beta  < zero  .or.  beta  > pi)     then
         stop "wigner_rotation_D_matrix(): program stop C."
      else if (gamma < zero  .or.  gamma > two*pi) then
         stop "wigner_rotation_D_matrix(): program stop D."
      end if
      !
      wa = exp(-i_cmplx * p * alpha) * wigner_rotation_d_small(j,p,q,beta) * exp(-i_cmplx * q * gamma)      
      ! 
   end function wigner_rotation_D_matrix
   !
   !
   function wigner_rotation_d_small(j,p,q,beta)                      result(wa)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the Wigner rotation matrix d^j_pq(beta).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: j, p, q
      real(kind=dp), intent(in) :: beta
      real(kind=dp)             :: wa
      !
      integer                   :: k, kmin1, kmax1, m, n
      !
      if (abs(p) > j  .or.  abs(q) > j) then
         stop "wigner_rotation_d_small(): program stop A."
      else if (beta  < zero  .or.  beta  > pi)     then
         stop "wigner_rotation_d_small(): program stop B."
      else if (phase_am(j-p) == one  .or.  phase_am(j-q) == one) then
         stop "wigner_rotation_d_small(): program stop C."
      end if
      !
      wa = zero 
      ! 
      kmin1 = max(0,-p-q)/2;   kmax1 =  min(j-p,j-q)/2
      !
      do  k = kmin1, kmax1
         m  = (p+q) / 2 + k+k;   n = (j+j-p-q) / 2 - k-k
         wa = wa + phase_am(k+k) * (cos(beta/two)**m) * (sin(beta/two)**n) /                                        &
                   (factorial_dp(k) * factorial_dp((j-p)/2 -k) * factorial_dp((j-q)/2 -k) * factorial_dp((p+q)/2 +k))
      end do;
      ! 
   end function wigner_rotation_d_small
   !
end module rabs_angular
