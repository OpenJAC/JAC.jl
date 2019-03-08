module rabs_dirac_orbital
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module defines all variables and arrays which are related to the symmetry blocks and to the quantum numbers of 
! the relativistic Dirac orbitals. It includes also several internal procedures to transform equivalent information 
! about these orbitals in a suitable form.
!
! This module stores the (one-electron) spinor coefficients and matrix elements for all symmetry blocks. It also keeps 
! the selected type of atomic basis set, a string which often determines the call of individual subprocedures for the 
! various basis sets.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_functions_math
   use rabs_naglib
   use rabs_nucleus
   implicit none
   !
   public :: angular_momentum_l 
                ! Returns the orbital angular momentum l for given kappa.
   public :: angular_momentum_j 
                ! Returns the angular momentum j for given kappa. 
   public :: angular_momentum_kappa 
                ! Returns the kappa angular momentum for given l and 2j.
   public :: angular_momentum_string
                ! Returns a string for printing an angular momentum as a character with given length (standard len = 5).
   public :: dirac_alpha_matrix  
                ! Returns one of the Dirac alpha matrices (alpha_x, ...).
   public :: dirac_orbital_energy_point  
                ! Returns the relativistic one-electron Dirac orbital energy for a point nucleus. 
   public :: dirac_orbital_rk_moment  
                ! Returns the <r^k> moment for a one-electron Dirac orbital and a point nucleus. 
   public :: get_kappa_channels
                ! Returns a list of kappa values which can 'couple' two total J^P symmetries together. 
   public :: get_kappa_from_name
                ! Returns the kappa value for a given orbital name string of (len=2).
   public :: get_subshell_occupation 
                ! Returns the subshell occupation for a given pqn and kappa.
   public :: legendre_Plm      
                ! Returns the value of the associate Legendre polynomial P_lm(theta). 
   public :: orbital_energy       
                ! Returns the one-electron orbital energy for a given orbital.
   public :: orbital_name       
                ! Returns the name of an orbital as a character(len=4) string.
   public :: orbital_symmetry   
                ! Returns the symmetry as a (len=2) string.
   public :: orbital_X_vector     
                ! Returns the X vector for a given orbital.
   public :: spherical_Clm      
                ! Returns the value of the normalized spherical harmonic C_lm(theta,phi).
   public :: spherical_Ylm     
                ! Returns the value of the spherical harmonic Y_lm(theta,phi).
   public :: spherical_Ylm_norm 
                ! Returns the value of the normalization integral over Y_lm and Y_l'm'.
   public :: spherical_Omega_km 
                ! Returns the values of the two-component spherical Dirac spinor Omega_kappa,m.
   public :: Xk_name       
                ! Returns the effective interaction strength X^k(a,b,c,d) as a character string using proper orbital 
                ! names.
   !
   integer, public                  :: kappa_min_block, kappa_max_block
   integer, public                  :: kappa_min_block_scf, kappa_max_block_scf
   integer, public                  :: number_of_orbitals
   integer, public, parameter       :: nw = 49
   integer, dimension(1:nw), public :: kappa_qn, principal_qn
   integer, dimension(1:nw), public :: orbital_occ
   !
   !
   ! Set a maximal |kappa| value to dimension several arrays
   integer, parameter, public :: kappamax = 10
   !
   integer,       dimension(-kappamax:kappamax), public :: bs_size
   real(kind=dp), dimension(-kappamax:kappamax), public :: alpha_N, beta_N
   !
   type(matrix_dp), dimension(:), pointer, public :: spinor_coeff
   type(matrix_dp), dimension(:), pointer, public :: spinor_overlap_ll
   type(matrix_dp), dimension(:), pointer, public :: spinor_overlap_ss
   type(matrix_dp), dimension(:), pointer, public :: spinor_nuclear_me_ll
   type(matrix_dp), dimension(:), pointer, public :: spinor_nuclear_me_ss
   type(matrix_dp), dimension(:), pointer, public :: spinor_pi_me_ls
   type(matrix_dp), dimension(:), pointer, public :: spinor_pi_me_sl
   type(vector_dp), dimension(:), pointer, public :: one_electron_energy
   !
   type(matrix_dp), private  :: fock_matrix, s_matrix, x_matrix 
   !
   integer,       dimension(:), pointer, private  :: iteration 
   real(kind=dp), dimension(:), pointer, private  :: alfr, alfi, beta 
   !
   character(len=12), public :: rabs_type
   !
   !
   !
   character(len=2), dimension(-10:70), private, parameter :: pqn_string = &
      (/ "xx", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1",   &
         " 0",                                                         &
         " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", "10",   &
         "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",   &
         "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",   &
         "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",   &
         "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",   &
         "51", "52", "53", "54", "55", "56", "57", "58", "59", "60",   &
         "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"    /)
   character(len=6), dimension(-120:120), private, parameter :: two_jm_string =                              &
      (/ "  -60 ", "-119/2", "  -59 ", "-117/2", "  -58 ", "-115/2", "  -57 ", "-113/2", "  -56 ", "-111/2", &
         "  -55 ", "-109/2", "  -54 ", "-107/2", "  -53 ", "-105/2", "  -52 ", "-103/2", "  -51 ", "-101/2", &
         "  -50 ", " -99/2", "  -49 ", " -97/2", "  -48 ", " -95/2", "  -47 ", " -93/2", "  -46 ", " -91/2", &
         "  -45 ", " -89/2", "  -44 ", " -87/2", "  -43 ", " -85/2", "  -42 ", " -83/2", "  -41 ", " -81/2", &
         "  -40 ", " -79/2", "  -39 ", " -77/2", "  -38 ", " -75/2", "  -37 ", " -73/2", "  -36 ", " -71/2", &
         "  -35 ", " -69/2", "  -34 ", " -67/2", "  -33 ", " -65/2", "  -32 ", " -63/2", "  -31 ", " -61/2", &
         "  -30 ", " -59/2", "  -29 ", " -57/2", "  -28 ", " -55/2", "  -27 ", " -53/2", "  -26 ", " -51/2", &
         "  -25 ", " -49/2", "  -24 ", " -47/2", "  -23 ", " -45/2", "  -22 ", " -43/2", "  -21 ", " -41/2", &
         "  -20 ", " -39/2", "  -19 ", " -37/2", "  -18 ", " -35/2", "  -17 ", " -33/2", "  -16 ", " -31/2", &
         "  -15 ", " -29/2", "  -14 ", " -27/2", "  -13 ", " -25/2", "  -12 ", " -23/2", "  -11 ", " -21/2", &
         "  -10 ", " -19/2", "   -9 ", " -17/2", "   -8 ", " -15/2", "   -7 ", " -13/2", "   -6 ", " -11/2", &
         "   -5 ", "  -9/2", "   -4 ", "  -7/2", "   -3 ", "  -5/2", "   -2 ", "  -3/2", "   -1 ", "  -1/2", &
         "    0 ",                                                                                           &
         "   1/2", "    1 ", "   3/2", "    2 ", "   5/2", "    3 ", "   7/2", "    4 ", "   9/2", "    5 ", &
         "  11/2", "    6 ", "  13/2", "    7 ", "  15/2", "    8 ", "  17/2", "    9 ", "  19/2", "   10 ", &
         "  21/2", "   11 ", "  23/2", "   12 ", "  25/2", "   13 ", "  27/2", "   14 ", "  29/2", "   15 ", &
         "  31/2", "   16 ", "  33/2", "   17 ", "  35/2", "   18 ", "  37/2", "   19 ", "  39/2", "   20 ", &
         "  41/2", "   21 ", "  43/2", "   22 ", "  45/2", "   23 ", "  47/2", "   24 ", "  49/2", "   25 ", &
         "  51/2", "   26 ", "  53/2", "   27 ", "  55/2", "   28 ", "  57/2", "   29 ", "  59/2", "   30 ", &
         "  61/2", "   31 ", "  63/2", "   32 ", "  65/2", "   33 ", "  67/2", "   34 ", "  69/2", "   35 ", &
         "  71/2", "   36 ", "  73/2", "   37 ", "  75/2", "   38 ", "  77/2", "   39 ", "  79/2", "   40 ", &
         "  81/2", "   41 ", "  83/2", "   42 ", "  85/2", "   43 ", "  87/2", "   44 ", "  89/2", "   45 ", &
         "  91/2", "   46 ", "  93/2", "   47 ", "  95/2", "   48 ", "  97/2", "   49 ", "  99/2", "   50 ", &
         " 101/2", "   51 ", " 103/2", "   52 ", " 105/2", "   53 ", " 107/2", "   54 ", " 109/2", "   55 ", &
         " 111/2", "   56 ", " 113/2", "   57 ", " 115/2", "   58 ", " 117/2", "   59 ", " 119/2", "   60 " /)    
   !!! character(len=2), dimension(-80:80), private, parameter :: kappa_string =  &
   !!!    (/ "80", "79", "78", "77", "76", "75", "74", "73", "72", "71",   &
   !!!       "70", "69", "68", "67", "66", "65", "64", "63", "62", "61",   &
   !!!       "60", "59", "58", "57", "56", "55", "54", "53", "52", "51",   &
   !!!       "50", "49", "48", "47", "46", "45", "44", "43", "42", "41",   &
   !!!       "40", "39", "38", "37", "36", "35", "34", "33", "32", "31",   &
   !!!       "30", "29", "28", "27", "26", "25", "24", "23", "22", "21",   &
   !!!       "y ", "x ", "w ", "v ", "u ", "t ", "r ", "q ", "o ", "n ",   &
   !!!       "m ", "l ", "k ", "i ", "h ", "g ", "f ", "d ", "p ", "s ",   &
   !!!       "xx",                                                         &
   !!!       "p-", "d-", "f-", "g-", "h-", "i-", "k-", "l-", "m-", "n-",   &
!!! 	 "o-", "q-", "r-", "t-", "u-", "v-", "w-", "x-", "y-", "z-",   &
!!! 	 "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",   &
!!! 	 "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",   &
!!! 	 "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",   &
!!! 	 "51", "52", "53", "54", "55", "56", "57", "58", "59", "60",   &
!!! 	 "61", "62", "63", "64", "65", "66", "67", "68", "69", "70",   &
!!! 	 "71", "72", "73", "74", "75", "76", "77", "78", "79", "80"    /)
   character(len=2), dimension(-80:80), private, parameter :: kappa_string =  &
      (/ "80", "79", "78", "77", "76", "75", "74", "73", "72", "71",   &
         "70", "69", "68", "67", "66", "65", "64", "63", "62", "61",   &
         "60", "59", "58", "57", "56", "55", "54", "53", "52", "51",   &
         "50", "49", "48", "47", "46", "45", "44", "43", "42", "41",   &
         "40", "39", "38", "37", "36", "35", "34", "33", "32", "31",   &
         "30", "29", "28", "27", "26", "25", "24", "23", "22", "21",   &
         "x ", "w ", "v ", "u ", "t ", "r ", "q ", "o ", "n ", "m ",   &
         "l ", "k ", "j ", "i ", "h ", "g ", "f ", "d ", "p ", "s ",   &
         "xx",                                                         &
         "p-", "d-", "f-", "g-", "h-", "i-", "j-", "k-", "l-", "m-",   &
	 "n-", "o-", "q-", "r-", "t-", "u-", "v-", "w-", "x-", "y-",   &
	 "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",   &
	 "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",   &
	 "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",   &
	 "51", "52", "53", "54", "55", "56", "57", "58", "59", "60",   &
	 "61", "62", "63", "64", "65", "66", "67", "68", "69", "70",   &
	 "71", "72", "73", "74", "75", "76", "77", "78", "79", "80"    /)
   !
contains
   !
   function angular_momentum_l(kappa)                          result(l)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the orbital angular momentum l(kappa).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kappa
      integer             :: l
      !
      if (kappa < 0) then
         l = -kappa -1
      elseif (kappa > 0) then
         l = kappa
      else if (rabs_use_stop) then 
         stop "angular_momentum_l(): program stop A."
      end if 
      !
   end function angular_momentum_l
   !
   !
   function angular_momentum_j(kappa)                          result(j)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the angular momentum 2*j(kappa).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kappa
      integer             :: j
      !
      if (kappa < 0  .or.  kappa > 0) then
         j = 2 * abs(kappa) - 1
      else if (rabs_use_stop) then 
         stop "angular_momentum_j(): program stop A."
      end if 
      !
   end function angular_momentum_j
   !
   !
   function angular_momentum_kappa(l,j2)                   result(kappa)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the kappa quantum number kappa(l,2*j).
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: l, j2
      integer              :: kappa
      !
      if (rabs_use_stop     .and.                      &
         (abs(l+l-j2) /= 1  .or.  l < 0   .or.  j2 < 1)) then
         stop "angular_momentum_kappa(): program stop A."
      else if (l+l < j2) then
         kappa = -l - 1
      else
         kappa = l
      end if
      !
   end function angular_momentum_kappa
   !
   !
   function angular_momentum_string(two_jm,length)               result(string)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a string for the angular momentum two_jm with given (optional) length; the standard is len = 6.
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)           :: two_jm
      integer, optional, intent(in) :: length
      character(len=6)              :: string
      !
      character(len=6) :: stringa
      !
      stringa = two_jm_string(two_jm)
      if (present(length)) then
         if (length > 0  .and.  length < 7) then
            stringa = adjustl(stringa)
            string(1:length)    = stringa(1:length)
            string(length+1:6)  = "     "
         else
           stop "angular_momentum_string(): program stop A."
         end if
      else 
         string = stringa
      end if
      !
   end function angular_momentum_string
   !
   !
   subroutine dirac_alpha_matrix(coordinate,alpha_i) 
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns one of the (4x4) Dirac alpha matrices as specified by the parameter 
   ! coordinate = { "0", "1", "2", "3", "x", "y", "z" }. For coordinate = "0" the beta-matrix is returned.
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      character(len=1), intent(in)                  :: coordinate 
      complex(kind=dp), dimension(4,4), intent(out) :: alpha_i 
      !
      alpha_i(1:4,1:4) = zero
      !
      if (coordinate == "0") then
         alpha_i(1,4) = one;   alpha_i(2,3) = one
         alpha_i(3,2) = one;   alpha_i(4,1) = one
      else if (coordinate == "1"  .or.  coordinate == "x") then
         alpha_i(1,4) = one;   alpha_i(2,3) = one
         alpha_i(3,2) = one;   alpha_i(4,1) = one
      else if (coordinate == "2"  .or.  coordinate == "y") then
         alpha_i(1,4) = cmplx(zero,-one,kind=dp)
         alpha_i(2,3) = cmplx(zero, one,kind=dp)
         alpha_i(3,2) = cmplx(zero,-one,kind=dp)
         alpha_i(4,1) = cmplx(zero, one,kind=dp)
      else if (coordinate == "3"  .or.  coordinate == "z") then
         alpha_i(1,3) = one;   alpha_i(2,4) = -one
         alpha_i(3,1) = one;   alpha_i(4,2) = -one
      else if (rabs_use_stop) then   
         stop "dirac_alpha_matrix(): program stop A."
      end if
      !
   end subroutine dirac_alpha_matrix
   !
   !
   function dirac_orbital_energy_point(pqn,kappa,Z)       result(energy)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the relativistic one-electron Dirac orbital energy for a point nucleus. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)       :: pqn, kappa
      real(kind=dp), intent(in) :: Z
      !
      integer                   :: nr
      real(kind=dp)             :: energy, gamma, N
      !
      nr     = pqn - abs(kappa)
      gamma  = sqrt( kappa*kappa - Z*Z / (c*c) )
      N      = sqrt( pqn*pqn - 2*nr*(abs(kappa) - gamma) )
      energy = c*c *( sqrt( one - Z*Z / (N*N*c*c) ) - one )
      !
   end function dirac_orbital_energy_point
   !
   !
   function dirac_orbital_rk_moment(pqn,kappa,Z,k)     result(rk_moment)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the <r^k> moment for a hydrogenic one-electron Dirac orbital in the case of a point nucleus. The formulas 
   ! below are taken from I.P. Grant, "Quantum electrodynamics of atoms", table 2.3, p.48 (unpublished).
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)       :: pqn, kappa, k
      real(kind=dp), intent(in) :: Z
      !
      integer                   :: nr
      real(kind=dp)             :: gamma, N, rk_moment
      !
      nr     = pqn - abs(kappa)
      gamma  = sqrt( kappa*kappa - Z*Z / (c*c) )
      N      = sqrt( pqn*pqn - 2*nr*(abs(kappa) - gamma) )
      !
      select case(k)
      case(2)
         rk_moment = 2 *( N*N*(5*N*N - 2*kappa*kappa)*(1 - Z*Z/(N*N*c*c)) + N*N*(1-gamma*gamma)                       &
                        - 3*kappa*N*N*sqrt(1 - Z*Z/(N*N*c*c)) ) / (4*Z*Z)
      case(1)
         rk_moment = ((3*N*N - kappa*kappa)*sqrt(1 - Z*Z/(N*N*c*c)) - kappa) /  (2*Z)     
      case(0)
         rk_moment = 1
      case(-1)
         rk_moment = ( pqn*gamma + (abs(kappa) - gamma)*abs(kappa) ) / ( 2*gamma*N*N*N ) * (2*Z)
      case(-2)
         rk_moment = ( kappa*kappa*sqrt(1 - Z*Z/(N*N*c*c)) ) / ( 2*gamma*gamma*N*N*N*(2*gamma + sign(kappa,kappa)) ) * &
                     (4*Z*Z)
      case(-3)
         rk_moment = ( -3*N*N*kappa*sqrt(1 - Z*Z/(N*N*c*c)) + N*N +2*gamma*gamma*kappa*kappa ) /                       &
                     ( 4*N*N*N*N*N*(gamma-1)*gamma*(gamma+1) * (2*gamma-1)*(2*gamma+1) ) * (8*Z*Z*Z)
      case default
         stop "dirac_orbital_rk_moment(): program stop A."
      end select
      !
   end function dirac_orbital_rk_moment
   !
   !
   subroutine get_kappa_channels(ja,parity_a,jb,parity_b,kappa_list,number)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a list of kappa values which can 'couple' two total J^P symmetries together, i.e. from 
   ! (ja,parity_a) --> (jb,parity_b).
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)               :: ja, jb
      character(len=1), intent(in)      :: parity_a, parity_b
      integer, intent(out)              :: number
      integer, dimension(:),intent(out) :: kappa_list
      !
      integer :: ipa, j, kappa, l
      ! 
      if (parity_a == "+") then
         ipa =  1
      else if (parity_a == "-") then
         ipa = -1
      else if (rabs_use_stop) then
         stop "get_kappa_channels(): program stop A."
      end if
      !
      number = 0
      do  kappa = 1,30
	 l = angular_momentum_l(-kappa)
	 j = angular_momentum_j(-kappa)
	 if ( ((ipa * ((-1)**l) ==  1 .and. parity_b == "+")  .or.  (ipa * ((-1)**l) == -1 .and. parity_b == "-")) .and. & 
              abs(ja - j) <= jb   .and.  ja + j >= jb ) then
	    number             = number + 1
	    kappa_list(number) = -kappa
	 end if
	 !
	 l = angular_momentum_l(kappa)
	 j = angular_momentum_j(kappa)
	 if ( ((ipa * ((-1)**l) ==  1 .and. parity_b == "+")  .or.  (ipa * ((-1)**l) == -1 .and. parity_b == "-")) .and. & 
              abs(ja - j) <= jb   .and.  ja + j >= jb ) then
	    number             = number + 1
	    kappa_list(number) = kappa
	 end if
      end do
      !
   end subroutine get_kappa_channels
   !
   !
   function get_kappa_from_name(string,fail)               result(kappa)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the kappa value of a given orbital name string of (len=2).
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      character(len=2), intent(in)   :: string
      logical, optional, intent(out) :: fail 
      integer                        :: i, kappa
      ! 
      do  i = -80,80
         if (string == kappa_string(i)) then
            kappa = i
            if (present(fail)) fail = .false.
            return
         end if
      end do
      !
      print *, "string = '"//string//"'"
      if (present(fail)) then
         fail = .true.
      else if (rabs_use_stop) then
         stop "get_kappa_from_name(): program stop A."
      end if
      !
   end function get_kappa_from_name
   !
   !
   function get_subshell_occupation(pqn,kappa)               result(occ)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the occupation of a given subshell pqn, kappa in the present run of the program.
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: pqn, kappa
      integer             :: i, occ
      ! 
      ! If subshell is not found
      occ = 0
      !
      do  i = 1,nw
         if (pqn == principal_qn(i)   .and.   kappa == kappa_qn(i)) then
            occ = orbital_occ(i)
            exit
         end if
      end do
      !
      if (rabs_use_stop   .and. (occ < 0   .or.   occ > angular_momentum_j(kappa) + 1)) then
         stop "get_subshell_occupation(): program stop A."
      end if 
      !
   end function get_subshell_occupation
   !
   !
   function legendre_Plm(l,m,theta)  result(Plm)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the associate Legendre polynomial P_lm (theta) for given quantum numbers l,m and the angle 
   ! theta. This function uses a simple tabulation of the associate Legendre polynomials for some low values of l and 
   ! stops otherwise. It does not deal with any general representation of these polynomials. This procedure has not yet 
   ! been implemented properly.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: l, m
      real(kind=dp), intent(in) :: theta
      real(kind=dp)             :: Plm
      !
      ! |m| <= l must apply
      if ( abs(m) .gt. l ) then
         stop "legendre_Plm(): Quantum numbers must fullfill |m| <= l."
      end if
      !
      select case(l)
      case(0)
         select case(m)
         case(0)
            Plm = sin(theta)
         case default
            stop "legendre_Plm(): program stop A."
         end select
      case(1)
         select case(m)
         case(-1)
         case(0)
         case(1)
         case default
            stop "legendre_Plm(): program stop B."
         end select
      case(2)
         select case(m)
         case(-2)
         case(-1)
         case(0)
         case(1)
         case(2)
         case default
            stop "legendre_Plm(): program stop C."
         end select
      case(3)
         select case(m)
         case(-3)
         case(-2)
         case(-1)
         case(0)
         case(1)
         case(2)
         case(3)
         case default
            stop "legendre_Plm(): program stop D."
         end select
      case default
         stop "legendre_Plm(): program stop E."
      end select
      !
   end function legendre_Plm
   !
   !
   function orbital_energy(pqz,kappa)                    result(epsilon)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the one-electron orbital energy E(pqz,kappa) where pqz is the principal quantum number. 
   !
   ! Calls: angular_momentum_l().
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)        :: pqz, kappa
      real(kind=dp)              :: epsilon
      !
      integer :: n, l
      !
      n     = bs_size(kappa)
      l     = angular_momentum_l(kappa)
      if (n > 0) then
         epsilon = one_electron_energy(kappa)%vector(n + pqz - l)
      else if (n < 0) then 
         epsilon = one_electron_energy(kappa)%vector(n + pqz + l + 1)
      else if (rabs_use_stop) then
         stop "orbital_energy(): program stop A."
      end if
      !
   end function orbital_energy
   !
   !
   function orbital_name(pqn,kappa)               result(orbital_string)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the name of an orbital as a character(len=4) string like  1s ,  2p-, etc from the principal quantum number 
   ! pqn and the kappa value; the first two character denote the principal qn and the last two the corresponding subshell.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: pqn, kappa
      character(len=4)    :: orbital_string
      !
      orbital_string(1:2) = pqn_string(pqn)
      orbital_string(3:4) = kappa_string(kappa)
      !
   end function orbital_name
   !
   !
   function orbital_symmetry(kappa)                  result(orbital_sym)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the orbital symmetry as a character(len=2) string like  s ,  p-, etc from the kappa value. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kappa
      character(len=2)    :: orbital_sym
      !
      if (abs(kappa) <= 80) then
         orbital_sym(1:2) = kappa_string(kappa)
      else
         orbital_sym(1:2) = "##"
      end if
      !
   end function orbital_symmetry
   !
   !
   subroutine orbital_X_vector(pqz,kappa,large,n,x_vector)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the vector of X coefficients for the orbital component P(pqz,kappa) where pqz is the principal quantum 
   ! number. The coefficients for the large component is chosen for large = .true. and those for the small component 
   ! otherwise. n is the number of coefficients and the dimension of the (output) x_vector.
   !
   ! Calls: angular_momentum_l().
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: pqz, kappa,n
      logical, intent(in) :: large
      real(kind=dp), dimension(:), intent(out) :: x_vector 
      !
      integer :: index, l
      !
      l     = angular_momentum_l(kappa)
      if (pqz > 0) then
         index = n + pqz - l
      else if (pqz < 0) then
         index = n + pqz + l + 1
      else if (rabs_use_stop) then
         stop "orbital_X_vector(): program stop A."
      end if
      if (large) then
         x_vector(1:n) = spinor_coeff(kappa)%matrix(1:n,index)
      else
         x_vector(1:n) = spinor_coeff(kappa)%matrix(n+1:n+n,index)
      endif
      !
   end subroutine orbital_X_vector
   !
   !
   function spherical_Clm(l,m,theta,phi)  result(Clm)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the normalized spherical harmonic C_lm(theta,phi) for given spherical angles theta and phi.
   !
   ! Calls: nag_s14aaf() [from NAG].
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)       :: l, m
      real(kind=dp), intent(in) :: theta, phi
      complex(kind=dp)          :: Clm, imphi
      !
      integer                   :: ifail, mp
      !! real(kind=dp)          :: nag_s14aaf
      !
      if (rabs_use_stop   .and.   abs(m) > l) then
         stop "spherical_Clm(): program stop A."
      else
         mp = abs(m)
      end if 
      !
      imphi = cmplx(zero,one,dp)
      imphi = imphi * m * phi
      Clm   = exp(imphi)
      if (rabs_use_naglib) then
         Clm   = Clm * sqrt( nag_s14aaf(real(l-m+1,dp),ifail) / nag_s14aaf(real(l+m+1,dp),ifail) )
      else
         Clm   = Clm * sqrt( one * factorial(l-m) / factorial(l+m) )
      end if
      Clm   = (-1)**m * Clm * legendre_Plm(l,m,theta)
      !
      if (m == -mp) then
         Clm = (-1)**mp * conjg(Clm)
      end if
      !
   end function spherical_Clm
   !
   !
   function spherical_Ylm(l,m,theta,phi)                     result(Ylm)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the spherical harmonic Y_lm(theta,phi) as usually defined for given spherical angles theta and 
   ! phi. The definition of the spherical harmonics in this procedure is taken from the monograph of Varshalovich 
   ! et al. (1988).
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)       :: l, m
      real(kind=dp), intent(in) :: theta, phi
      complex(kind=dp)          :: Ylm, iphi
      !
      if (rabs_use_stop   .and.  l < 0) then
         stop "spherical_Ylm(): program stop A."
      else if (abs(m) > l) then
         Ylm = zero
         return
      end if 
      !
      ! Ylm = sqrt( (2*l+1) / (two*two*pi) ) * spherical_Clm(l,m,theta,phi)
      !
      iphi = cmplx(zero,phi,dp)
      select case(l)
      case(0)
         select case(m)
         case(0)
            Ylm = one / (two*sqrt(pi))
         case default
            stop "spherical_Ylm(): program stop B."
         end select
      case(1)
         select case(m)
         case(1)
            Ylm = - sqrt(three/(two*pi)) / two * sin(theta) * exp(iphi)
         case(0)
            Ylm = sqrt(three/pi) / two * cos(theta)
         case(-1)
            Ylm =   sqrt(three/(two*pi)) / two * sin(theta) * exp(-iphi)
         case default
            stop "spherical_Ylm(): program stop C."
         end select
      case(2)
         select case(m)
         case(2)
            Ylm = sqrt(3*5/(two*pi)) / (two*two) * sin(theta) * sin(theta) * exp(iphi*two)
         case(1)
            Ylm = - sqrt(3*5/(two*pi)) / two * cos(theta) * sin(theta) * exp(iphi)
         case(0)
            Ylm = sqrt(5/pi) / (two*two) * (three * cos(theta) * cos(theta) - one)
         case(-1)
            Ylm = sqrt(3*5/(two*pi)) / two * cos(theta) * sin(theta) * exp(-iphi)
         case(-2)
            Ylm = sqrt(3*5/(two*pi)) / (two*two) * sin(theta) * sin(theta) * exp(-iphi*two)
         case default
            stop "spherical_Ylm(): program stop D."
         end select
      case(3)
         select case(m)
         case(3)
            Ylm = - sqrt(5*7/pi) / (two*two*two) * sin(theta) * sin(theta) * sin(theta) * exp(iphi*three)
         case(2)
            Ylm = sqrt(3*5*7/(two*pi)) / (two*two) * cos(theta) * sin(theta) * sin(theta) * exp(iphi*two)
         case(1)
            Ylm = - sqrt(3*7/pi) / (two*two*two) * (5*cos(theta) * cos(theta) - one) * sin(theta) * exp(iphi)
         case(0)
            Ylm = sqrt(7/pi) / (two*two) * (5*cos(theta) * cos(theta) - three) * cos(theta)
         case(-1)
            Ylm = sqrt(3*7/pi) / (two*two*two) * (5*cos(theta) * cos(theta) - one) * sin(theta) * exp(-iphi)
         case(-2)
            Ylm = sqrt(3*5*7/(two*pi)) / (two*two) * cos(theta) * sin(theta) * sin(theta) * exp(-iphi*two)
         case(-3)
            Ylm = sqrt(5*7/pi) / (two*two*two) * sin(theta) * sin(theta) * sin(theta) * exp(-iphi*three)
         case default
            stop "spherical_Ylm(): program stop E."
         end select
      case(4)
         select case(m)
         case(4)  
            Ylm = three * sqrt(5*7/(two*pi)) / (two*two*two*two) * (sin(theta)**4) * exp(iphi*two*two)
         case(3)
            Ylm = -three * sqrt(5*7/pi) / (two*two*two) * (sin(theta)**3) * cos(theta) * exp(iphi*three)
         case(2)
            Ylm = three * sqrt(5/(two*pi)) / (two*two*two) * sin(theta) * sin(theta) * (7*cos(theta)*cos(theta) - one)   &
                  * exp(iphi*two)
         case(1)
            Ylm = -three * sqrt(5/pi) / (two*two*two) * sin(theta) * (7*(cos(theta)**3) - three*cos(theta)) * exp(iphi)
         case(0)
            Ylm = three * sqrt(1/pi) / (two*two*two*two) * (35*(cos(theta)**4) - 30*cos(theta)*cos(theta) + three)   
         case(-1)
            Ylm = three * sqrt(5/pi) / (two*two*two) * sin(theta) * (7*(cos(theta)**3) - three*cos(theta))* exp(-iphi)
         case(-2)
            Ylm = three * sqrt(5/(two*pi)) / (two*two*two) * sin(theta) * sin(theta) * (7*cos(theta)*cos(theta) - one)   &
                  * exp(-iphi*two)
         case(-3)
            Ylm = three * sqrt(5*7/pi) / (two*two*two) * (sin(theta)**3) * cos(theta) * exp(-iphi*three)
         case(-4)
            Ylm = three * sqrt(5*7/(two*pi)) / (two*two*two*two) * (sin(theta)**4) * exp(-iphi*two*two)
         case default
            stop "spherical_Ylm(): program stop F."
         end select
      case(5)
         select case(m)
         case(-5)
            Ylm = three/32.0_dp * sqrt(77.0_dp/pi)* (sin(theta)**5) * exp(-5.0_dp*iphi)
         case(-4)
            Ylm = three/32.0_dp * sqrt(770.0_dp/pi)* (sin(theta)**4) * cos(theta) * exp(-two*two*iphi)
         case(-3)
            Ylm = one/32.0_dp * sqrt(385.0_dp/pi)* (sin(theta)**3) * (-one+9.0_dp*(cos(theta)**2)) * exp(-three*iphi)
         case(-2)
            Ylm = one/16.0_dp * sqrt(2310.0_dp/pi) * (sin(theta)**2) * cos(theta) * (-one+three*(cos(theta)**2)) * &
                  exp(-two*iphi)
         case(-1)
            Ylm = one/32.0_dp * sqrt(330.0_dp/pi) * sin(theta)* (one-14.0_dp*(cos(theta)**2)+21.0_dp*(cos(theta)**4)) * &
                  exp(-iphi)
         case(0)
            Ylm = one/16.0_dp * sqrt(11.0_dp/pi) * cos(theta) * (15.0_dp-70.0_dp*(cos(theta)**2)+63.0_dp*(cos(theta)**4))
         case(1)
            Ylm = - one/32.0_dp * sqrt(330.0_dp/pi) * sin(theta) * (one-14.0_dp*(cos(theta)**2)+21*(cos(theta)**4)) *   &
                  exp(iphi)
         case(2)
            Ylm = - one/16.0_dp * sqrt(2310.0_dp/pi) * (sin(theta)**2) * cos(theta) * (-one+three*(cos(theta)**2)) *    &
                  exp(two*iphi)
         case(3)
            Ylm = - one/32.0_dp * sqrt(385.0_dp/pi) * (sin(theta)**3) * (-one+9.0_dp*(cos(theta)**2)) * exp(three*iphi)
         case(4)
            Ylm = three/32.0_dp * sqrt(770.0_dp/pi) * (sin(theta)**4) * cos(theta) * exp(two*two*iphi)
         case(5)
            Ylm = -three/32.0_dp * sqrt(77.0_dp/pi) * (sin(theta)**5) * exp(5.0_dp*iphi)
         case default
            stop "spherical_Ylm(): program stop G."
         end select
      case(6)
         select case(m)
         case(-6)
            Ylm = one/64.0_dp * sqrt(3003.0_dp/pi) * (sin(theta)**6) * exp(-6.0_dp*iphi)
         case(-5)
            Ylm = three/32.0_dp * sqrt(1001.0_dp/pi) * (sin(theta)**5) * cos(theta) * exp(-5.0_dp*iphi)
         case(-4)
            Ylm = three/64.0_dp * sqrt(182.0_dp/pi) * (sin(theta)**4) * (-one+11.0_dp *(cos(theta)**2)) * &
                  exp(-two*two*iphi)
         case(-3)
            Ylm = one/32.0_dp * sqrt(1365.0_dp/pi) * (sin(theta)**3) * cos(theta) * (-three+11.0_dp*(cos(theta)**2)) *  &
                  exp(-three*iphi)
         case(-2)
            Ylm = one/64.0_dp * sqrt(1365.0_dp/pi) * (sin(theta)**2) *                &
                  (one-18.0_dp*(cos(theta)**2)+33.0_dp*(cos(theta)**4)) * exp(-two*iphi)
         case(-1)
            Ylm = one/32.0_dp * sqrt(546.0_dp/pi) * sin(theta)*cos(theta) *           &
                  (5.0_dp-30.0_dp*(cos(theta)**2)+33.0_dp*(cos(theta)**4))* exp(iphi)
         case(0)
            Ylm = one/32.0_dp * sqrt(13.0_dp/pi) * (-5.0_dp+  105.0_dp*(cos(theta)**2)-315.0_dp*(cos(theta)**4)+   &
                  231.0_dp*(cos(theta)**6))
         case(1)
            Ylm = -one/32.0_dp * sqrt(546.0_dp/pi) * sin(theta)*cos(theta)*           &
                  (5.0_dp-30.0_dp*(cos(theta)**2)+33.0_dp*(cos(theta)**4))* exp(iphi)
         case(2)
            Ylm = one/64.0_dp * sqrt(1365.0_dp/pi) * (sin(theta)**2) *                &
                  (one-18.0_dp*(cos(theta)**2)+33.0_dp*(cos(theta)**4)) * exp(two*iphi)
         case(3)
            Ylm = -one/32.0_dp * sqrt(1365.0_dp/pi) * (sin(theta)**3) *       &
                  cos(theta) * (-three+11.0_dp*(cos(theta)**2)) * exp(three*iphi)
         case(4)
            Ylm = three/64.0_dp * sqrt(182.0_dp/pi) * (sin(theta)**4) * (-one+11.0_dp*(cos(theta)**2)) * exp(two*two*iphi)
         case(5)
            Ylm = -three/32.0_dp * sqrt(1001.0_dp/pi) * (sin(theta)**5) * cos(theta) * exp(5.0_dp*iphi)
         case(6)
            Ylm = one/64.0_dp * sqrt(3003.0_dp/pi) * (sin(theta)**6) * exp(6.0_dp*iphi)
         case default
            stop "spherical_Ylm(): program stop H."
         end select
      case(7)
         select case(m)
         case(-7)
            Ylm = three/128.0_dp * sqrt(1430.0_dp/pi) * (sin(theta)**7) * exp(-7.0_dp*iphi)
         case(-6)
            Ylm = three/64.0_dp * sqrt(5005.0_dp/pi) * (sin(theta)**6) * cos(theta) * exp(-6.0_dp*iphi)
         case(-5)
            Ylm = three/128.0_dp * sqrt(770.0_dp/pi) * (sin(theta)**5) * (-one+13.0_dp*(cos(theta)**2)) * exp(-5.0_dp*iphi)
         case(-4)
            Ylm = three/64.0_dp * sqrt(770.0_dp/pi) * (sin(theta)**4) *       &
                  cos(theta) * (-three+13.0_dp*(cos(theta)**2)) * exp(-4.0_dp*iphi)
         case(-3)
            Ylm = three/128.0_dp * sqrt(70.0_dp/pi) * (sin(theta)**3) *       &
                  (three-66.0_dp*(cos(theta)**2)+143.0_dp*(cos(theta)**4)) * exp(-three*iphi)
         case(-2)
            Ylm = three/64.0_dp * sqrt(35.0_dp/pi) * (sin(theta)**2) * cos(theta) * (15.0_dp-110.0_dp*(cos(theta)**2)+  &
                  143.0_dp*(cos(theta)**4)) * exp(-two*iphi)
         case(-1)
            Ylm = one/128.0_dp * sqrt(210.0_dp/pi) * sin(theta)* (-5.0_dp+    &
                  135.0_dp*(cos(theta)**2)-495.0_dp*(cos(theta)**4)+ 429.0_dp*(cos(theta)**6)) * exp(-iphi)
         case(0)
            Ylm = one/32.0_dp * sqrt(15.0_dp/pi) * cos(theta) * (-35.0_dp+    &
                  315.0_dp*(cos(theta)**2)-693.0_dp*(cos(theta)**4)+ 429.0_dp*(cos(theta)**6))
         case(1)
            Ylm = -one/128.0_dp * sqrt(210.0_dp/pi) * sin(theta) * (-5.0_dp+  &
                  135.0_dp*(cos(theta)**2)-495.0_dp*(cos(theta)**4)+ 429.0_dp*(cos(theta)**6)) * exp(iphi)
         case(2)
            Ylm = three/64.0_dp * sqrt(35.0_dp/pi) * (sin(theta)**2) *        &
                  cos(theta) * (15.0_dp-110.0_dp*(cos(theta)**2)+ 143.0_dp*(cos(theta)**4)) * exp(two*iphi)
         case(3)
            Ylm = -three/128.0_dp * sqrt(70.0_dp/pi) * (sin(theta)**3) *      &
                  (three-66.0_dp*(cos(theta)**2)+143.0_dp*(cos(theta)**4)) * exp(three*iphi)
         case(4)
            Ylm = three/64.0_dp* sqrt(770.0_dp/pi) * (sin(theta)**4) * cos(theta) * (-three+13.0_dp*(cos(theta)**2)) *  &
                  exp(two*two*iphi)
         case(5)
            Ylm = -three/128.0_dp * sqrt(770.0_dp/pi) * (sin(theta)**5) * (-one+13.0_dp*(cos(theta)**2)) * exp(5.0_dp*iphi)
         case(6)
            Ylm = three/64.0_dp * sqrt(5005.0_dp/pi) * (sin(theta)**6) * cos(theta) * exp(6.0_dp*iphi)
         case(7)
            Ylm = -three/128.0_dp * sqrt(1430.0_dp/pi) * (sin(theta)**7) * exp(7.0_dp*iphi)
         case default
            stop "spherical_Ylm(): program stop I."
         end select
      case(8)
         select case(m)
         case(-8)
            Ylm = three/512.0_dp * sqrt(24310.0_dp/pi) * (sin(theta)**8) * exp(-8.0_dp*iphi)
         case(-7)
            Ylm = three/128.0_dp * sqrt(24310.0_dp/pi) * (sin(theta)**7) * cos(theta) * exp(-7.0_dp*iphi)
         case(-6)
            Ylm = one/128.0_dp * sqrt(7293.0_dp/pi) * (sin(theta)**6) * (-one+15.0_dp*(cos(theta)**2)) * exp(-6.0_dp*iphi)
         case(-5)
            Ylm = three/128.0_dp * sqrt(34034.0_dp/pi) * (sin(theta)**5) *    &
                  cos(theta) * (-one+5.0_dp*(cos(theta)**2))*exp(-5.0_dp*iphi)
         case(-4)
            Ylm = three/256.0_dp * sqrt(2618.0_dp/pi) * (sin(theta)**4) *     &
                  (one-26.0_dp*(cos(theta)**2)+65.0_dp*(cos(theta)**4)) * exp(-4.0_dp*iphi)
         case(-3)
            Ylm = one/128.0_dp * sqrt(39270.0_dp/pi) * (sin(theta)**3) *      &
                  cos(theta) * (three-26.0_dp*(cos(theta)**2)+ 39.0_dp*(cos(theta)**4)) * exp(-three*iphi)
         case(-2)
            Ylm = three/128.0_dp * sqrt(595.0_dp/pi) * (sin(theta)**2) *      &
                  (-one+33.0_dp*(cos(theta)**2)-143.0_dp*(cos(theta)**4)+ 143.0_dp*(cos(theta)**6)) * exp(-two*iphi)
         case(-1)
            Ylm = three/128.0_dp * sqrt(34.0_dp/pi) * sin(theta) *cos(theta)* &
                (-35.0_dp+385.0_dp*(cos(theta)**2)-1001.0_dp*(cos(theta)**4)+ 715.0_dp*(cos(theta)**6)) * exp(-iphi)
         case(0)
            Ylm = one/256.0_dp * sqrt(17.0_dp/pi) * (35.0_dp-                 &
                  1260.0_dp*(cos(theta)**2)+6930.0_dp*(cos(theta)**4)- 12012.0_dp*(cos(theta)**6)+6435.0_dp*(cos(theta)**8))
         case(1)
            Ylm = -three/128.0_dp * sqrt(34.0_dp/pi) * sin(theta)*cos(theta)* &
                (-35.0_dp+385.0_dp*(cos(theta)**2)-1001.0_dp*(cos(theta)**4)+ 715.0_dp*(cos(theta)**6)) * exp(iphi)
         case(2)
            Ylm = three/128.0_dp * sqrt(595.0_dp/pi) * (sin(theta)**2) *      &
                  (-one+33.0_dp*(cos(theta)**2)-143.0_dp*(cos(theta)**4)+ 143.0_dp*(cos(theta)**6)) *exp(two*iphi)
         case(3)
            Ylm = -one/128.0_dp * sqrt(39270.0_dp/pi) * (sin(theta)**3) *     &
                  cos(theta) * (three-26.0_dp*(cos(theta)**2)+ 39.0_dp*(cos(theta)**4)) * exp(three*iphi)
         case(4)
            Ylm = three/256.0_dp * sqrt(2618.0_dp/pi) * (sin(theta)**4) *     &
                  (one-26.0_dp*(cos(theta)**2)+65.0_dp*(cos(theta)**4)) * exp(4.0_dp*iphi) 
         case(5)
            Ylm = -three/128.0_dp * sqrt(34034.0_dp/pi) * (sin(theta)**5) *   &
                  cos(theta) * (-one+5.0_dp*(cos(theta)**2))*exp(5.0_dp*iphi)
         case(6)
            Ylm = one/128.0_dp * sqrt(7293.0_dp/pi) * (sin(theta)**6) *  (-one+15.0_dp*(cos(theta)**2)) * exp(6.0_dp*iphi)
         case(7)
            Ylm = -three/128.0_dp * sqrt(24310.0_dp/pi) * (sin(theta)**7) * cos(theta) * exp(7.0_dp*iphi)
         case(8)
            Ylm = three/512.0_dp * sqrt(24310.0_dp/pi) * (sin(theta)**8) * exp(8.0_dp*iphi)
         case default
            stop "spherical_Ylm(): program stop J."
         end select
      case(9)
         select case(m)
         case(-9)
            Ylm = one/512.0_dp * sqrt(230945.0_dp/pi) * (sin(theta)**9) * exp(-9.0_dp*iphi)
         case(-8)
            Ylm = three/512.0_dp * sqrt(461890.0_dp/pi) * (sin(theta)**8) * cos(theta) * exp(-8.0_dp*iphi)
         case(-7)
            Ylm = three/512.0_dp * sqrt(13585.0_dp/pi) * (sin(theta)**7) * (-one+17.0_dp*(cos(theta)**2)) * &
                  exp(-7.0_dp*iphi)
         case(-6)
            Ylm = one/128.0_dp * sqrt(40755.0_dp/pi) * (sin(theta)**6) *      &
                  cos(theta)*(-three+17.0_dp*(cos(theta)**2))*exp(-6.0_dp*iphi)
         case(-5)
            Ylm = three/256.0_dp * sqrt(2717.0_dp/pi) * (sin(theta)**5) *     &
                  (one-30.0_dp*(cos(theta)**2)+85.0_dp*(cos(theta)**4)) * exp(-5.0_dp*iphi)
         case(-4)
            Ylm = three/256.0_dp * sqrt(190190.0_dp/pi) * (sin(theta)**4) *   &
                  cos(theta) * (one-10.0_dp*(cos(theta)**2)+ 17.0_dp*(cos(theta)**4)) *exp(-4.0_dp*iphi)
         case(-3)
            Ylm = one/256.0_dp * sqrt(21945.0_dp/pi) * (sin(theta)**3) *      &
                  (-one+39.0_dp*(cos(theta)**2)-195.0_dp*(cos(theta)**4)+ 221.0_dp*(cos(theta)**6)) * exp(-three*iphi)
         case(-2)
            Ylm = three/128.0_dp * sqrt(1045.0_dp/pi) * (sin(theta)**2) *     &
                  cos(theta) * (-7.0_dp+91.0_dp*(cos(theta)**2)- 273.0_dp*(cos(theta)**4)+221.0_dp*(cos(theta)**6)) *        &
                  exp(-two*iphi)
         case(-1)
            Ylm = three/512.0_dp * sqrt(190.0_dp/pi) * sin(theta) *           &
                  (7.0_dp-308.0_dp*(cos(theta)**2)+2002.0_dp*(cos(theta)**4)- &
                4004.0_dp*(cos(theta)**6)+2431.0_dp*(cos(theta)**8))*exp(-iphi)
         case(0)
            Ylm = one/256.0_dp * sqrt(19.0_dp/pi) * cos(theta) *              &
                  (315.0_dp-4620.0_dp*(cos(theta)**2)+ 18018.0_dp*(cos(theta)**4)-25740.0_dp*(cos(theta)**6)+      &
                  12155.0_dp*(cos(theta)**8))
         case(1)
            Ylm = -three/512.0_dp * sqrt(190.0_dp/pi) * sin(theta) *          &
                  (7.0_dp-308.0_dp*(cos(theta)**2)+2002.0_dp*(cos(theta)**4)- &
                   4004.0_dp*(cos(theta)**6)+2431.0_dp*(cos(theta)**8))*exp(iphi)
         case(2)
            Ylm = three/128.0_dp * sqrt(1045.0_dp/pi) * (sin(theta)**2) *     &
                  cos(theta) * (-7.0_dp+91.0_dp*(cos(theta)**2)- 273.0_dp*(cos(theta)**4)+221.0_dp*(cos(theta)**6)) *  &
                  exp(2*iphi)
         case(3)
            Ylm = -one/256.0_dp * sqrt(21945.0_dp) * (sin(theta)**3) *        &
                  (-one+39.0_dp*(cos(theta)**2)-195.0_dp*(cos(theta)**4)+ 221.0_dp*(cos(theta)**6)) * exp(three*iphi)
         case(4)
            Ylm = three/256.0_dp * sqrt(190190.0_dp/pi) * (sin(theta)**4) *   &
                  cos(theta) * (one-10.0_dp*(cos(theta)**2)+ 17.0_dp*(cos(theta)**4)) * exp(4.0_dp*iphi)
         case(5)
            Ylm = -three/256.0_dp * sqrt(2717.0_dp/pi) * (sin(theta)**5) *    &
                  (one-30.0_dp*(cos(theta)**2)+85*(cos(theta)**4)) * exp(5.0_dp*iphi)
         case(6)
            Ylm = one/128.0_dp * sqrt(40755.0_dp/pi) * (sin(theta)**6) *      &
                  cos(theta) * (-three+17.0_dp*(cos(theta)**2)) * exp(6.0_dp*iphi) 
         case(7)
            Ylm = -three/512.0_dp * sqrt(13585.0_dp/pi) * (sin(theta)**7) *   &
                  (-one+17.0_dp*(cos(theta)**2)) * exp(7.0_dp*iphi)
         case(8)
            Ylm = three/512.0_dp * sqrt(461890.0_dp/pi) * (sin(theta)**8) * cos(theta) * exp(8.0_dp*iphi)
         case(9)
            Ylm = -one/512.0_dp * sqrt(230945.0_dp/pi) * (sin(theta)**9) * exp(9.0_dp*iphi)
         case default
            stop "spherical_Ylm(): program stop K."
         end select
      case(10)
         select case(m)
         case(-10)
            Ylm = one/1024.0_dp * sqrt(969969.0_dp/pi) * (sin(theta)**10) * exp(-10.0_dp*iphi)
         case(-9)
            Ylm = one/512.0_dp * sqrt(0.4849845e7_dp/pi) * (sin(theta)**9) * cos(theta) * exp(-9.0_dp*iphi)
         case(-8)
            Ylm = one/1024.0_dp * sqrt(510510.e0_dp/pi) * (sin(theta)**8) *   &
                  (-one+19.0_dp*(cos(theta)**2)) * exp(-8.0_dp*iphi)
         case(-7)
            Ylm = three/512.0_dp * sqrt(85085.e0_dp/pi) * (sin(theta)**7) *   &
                  cos(theta) * (-three+19.0_dp*(cos(theta)**2)) * exp(-7.0_dp*iphi)
         case(-6)
            Ylm = three/1024.0_dp * sqrt(5005.e0_dp/pi) * (sin(theta)**6) *   &
                  (three-102.0_dp*(cos(theta)**2)+323.0_dp*(cos(theta)**4)) * exp(-6.0_dp*iphi)
         case(-5)
            Ylm = three/256.0_dp * sqrt(1001.e0_dp/pi) * (sin(theta)**5) * cos(theta)  *(15.0_dp-170.0_dp*(cos(theta)**2)+ &
                  323.0_dp*(cos(theta)**4)) * exp(-5.0_dp*iphi)
         case(-4)
            Ylm = three/512.0_dp * sqrt(10010.e0_dp/pi) * (sin(theta)**4) *   &
                  (-one+45.0_dp*(cos(theta)**2)-255.0_dp*(cos(theta)**4)+ 323.0_dp*(cos(theta)**6))*exp(-4.0_dp*iphi)
         case(-3)
            Ylm = three/256.0_dp * sqrt(5005.e0_dp/pi) * (sin(theta)**3) * cos(theta) * &
                  (-7.0_dp+105.0_dp*(cos(theta)**2)- 357.0_dp*(cos(theta)**4)+323.0_dp*(cos(theta)**6)) * exp(-three*iphi)
         case(-2)
            Ylm = three/1024.0_dp * sqrt(770.0_dp/pi) * (sin(theta)**2) *     &
                  (7.0_dp-364.0_dp*(cos(theta)**2)+2730.0_dp*(cos(theta)**4)- &
                  6188.0_dp*(cos(theta)**6)+4199.0_dp*(cos(theta)**8)) * exp(-two*iphi)
         case(-1)
            Ylm = one/512.0_dp * sqrt(2310.0_dp/pi) * sin(theta)*cos(theta)*  &
               (63.0_dp-1092.0_dp*(cos(theta)**2)+4914.0_dp*(cos(theta)**4)-  &
               7956.0_dp*(cos(theta)**6)+4199.0_dp*(cos(theta)**8)) * exp(-iphi)
         case(0)
            Ylm = one/512.0_dp * sqrt(21.0_dp/pi) * (-63.0_dp+  3465.0_dp*(cos(theta)**2)-  30030.0_dp*(cos(theta)**4) &
                  +90090.0_dp*(cos(theta)**6) -  109395.0_dp*(cos(theta)**8)+46189.0_dp*(cos(theta)**10))
         case(1)
            Ylm = -one/512.0_dp * sqrt(2310.0_dp/pi) * sin(theta)*cos(theta)* &
                (63.0_dp-1092.0_dp*(cos(theta)**2)+4914.0_dp*(cos(theta)**4)- &
                7956.0_dp*(cos(theta)**6)+4199.0_dp*(cos(theta)**8)) * exp(iphi)
         case(2)
            Ylm = three/1024.0_dp * sqrt(770.0_dp/pi) * (sin(theta)**2) *     &
                  (7.0_dp-364.0_dp*(cos(theta)**2)+2730.0_dp*(cos(theta)**4)- &
                  6188.0_dp*(cos(theta)**6)+4199.0_dp*(cos(theta)**8)) * exp(two*iphi)
         case(3)
            Ylm = -three/256.0_dp * sqrt(5005.0_dp/pi) * (sin(theta)**3) *  cos(theta) *                         &
                  (-7.0_dp+105.0_dp*(cos(theta)**2)- 357.0_dp*(cos(theta)**4)+323.0_dp*(cos(theta)**6)) * exp(three*iphi)
         case(4)
            Ylm = three/512.0_dp * sqrt(10010.0_dp/pi) * (sin(theta)**4) * (-one+45.0_dp*(cos(theta)**2)-255.0_dp* &
                  (cos(theta)**4)+  323.0_dp*(cos(theta)**6)) * exp(4.0_dp*iphi)
         case(5)
            Ylm = -three/256.0_dp * sqrt(1001.0_dp/pi) * (sin(theta)**5) * cos(theta) * (15.0_dp-170.0_dp*(cos(theta)**2) &
                  +323.0_dp*(cos(theta)**4)) * exp(5.0_dp*iphi)
         case(6)
            Ylm = three/1024.0_dp * sqrt(5005.0_dp/pi) * (sin(theta)**6) *    &
                  (three-102.0_dp*(cos(theta)**2)+323.0_dp*(cos(theta)**4)) * exp(6.0_dp*iphi)
         case(7)
            Ylm = -three/512.0_dp * sqrt(85085.0_dp/pi) * (sin(theta)**7) *   &
                  cos(theta) * (-three+19.0_dp*(cos(theta)**2))*exp(7.0_dp*iphi)
         case(8)
            Ylm = one/1024.0_dp * sqrt(510510.0_dp/pi) * (sin(theta)**8) * (-one+19.0_dp*(cos(theta)**2))*exp(8.0_dp*iphi)
         case(9)
            Ylm = -one/512.0_dp * sqrt(0.4849845e7_dp/pi) * (sin(theta)**9) * cos(theta) * exp(9.0_dp*iphi)
         case(10)
            Ylm = one/1024.0_dp * sqrt(969969.0_dp/pi) * (sin(theta)**10) * exp(10.0_dp*iphi)
         case default
            stop "spherical_Ylm(): program stop L."
         end select
      case default
         stop "spherical_Ylm(): program stop M."
      end select
   end function spherical_Ylm
   !
   !
   function spherical_Ylm_norm(la,ma,lb,mb,ngl)             result(norm)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the value of the normalization and orthogonalization relation of the spherical harmonics, i.e.
   !
   !    integral{ d_theta d_phi sin(theta) Y_la,ma(theta,phi)^* Y_lb,mb(theta,phi) }
   !
   ! using a two--dimensional Gauss-Legendre integration scheme over (theta,phi) and ngl zeros in each dimension.
   !
   ! Calls: set_gauss_zeros().
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: la, ma, lb, mb, ngl
      complex(kind=dp)    :: norm
      !
      integer             :: i, j
      real(kind=dp)       :: phi, theta
      real(kind=dp), dimension(:), pointer :: xtheta, xphi, wtheta, wphi
      !
      if (rabs_use_stop  .and.              &
         (abs(ma) > la   .or.  abs(mb) > lb)) then
         stop "spherical_Ylm_norm(): program stop A."
      end if 
      !
      allocate( xtheta(1:ngl), wtheta(1:ngl) )
      allocate( xphi(1:ngl),   wphi(1:ngl) )
      !
      norm = cmplx(zero,zero,dp)
      !
      call set_gauss_zeros(-pi, pi,ngl,xphi,wphi)
      call set_gauss_zeros(zero,pi,ngl,xtheta,wtheta)
      !
      ! Sum over all contributions
      do  i = 1,ngl
         do  j = 1,ngl
            phi   = xphi(i)
            theta = xtheta(j)
            norm  = norm + wphi(i) * wtheta(j) * sin(theta) * conjg(spherical_Ylm(la,ma,theta,phi)) *       &
                           spherical_Ylm(lb,mb,theta,phi)
         end do
      end do
      !
      deallocate( xtheta, wtheta, xphi, wphi )
      !
   end function spherical_Ylm_norm
   !
   !
   subroutine spherical_Omega_km(kappa,mm,theta,phi,omega1,omega2)  
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the values of the two-component spherical Dirac spinor Omega_kappa,m(theta,phi) for given spherical angles 
   ! theta and phi. The cases kappa < 0 and kappa > 0 are treated separately.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(), spherical_Ylm().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)           :: kappa, mm
      real(kind=dp), intent(in)     :: theta, phi
      complex(kind=dp), intent(out) :: omega1, omega2
      !
      integer                       :: j, l
      !
      j = angular_momentum_j(kappa)
      l = angular_momentum_l(kappa)
      !
      ! Treat kappa < 0 (a = +1)
      if (kappa < 0) then
         omega1 =   sqrt( real(j+mm,dp) / real(j+j,dp) ) * spherical_Ylm((j-1)/2,(mm-1)/2,theta,phi)
         omega2 =   sqrt( real(j-mm,dp) / real(j+j,dp) ) * spherical_Ylm((j-1)/2,(mm+1)/2,theta,phi)
      !
      ! Treat kappa > 0 (a = -1)
      elseif (kappa > 0) then
         omega1 = - sqrt( real(j+2-mm,dp) / real(j+j+4,dp) ) * spherical_Ylm((j+1)/2,(mm-1)/2,theta,phi)
         omega2 =   sqrt( real(j+2+mm,dp) / real(j+j+4,dp) ) * spherical_Ylm((j+1)/2,(mm+1)/2,theta,phi)
      end if
      !
   end subroutine spherical_Omega_km
   !
   !
   function Xk_name(k,an,akappa,bn,bkappa,cn,ckappa,dn,dkappa) result(Xk_string)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the effective interaction strength X^k(a,b,c,d) as a character(len=*) string using proper orbital names 
   ! such as 1s ,  2p-, etc from the principal quantum numbers an, bn, ... and the kappa value akappa, bkappa, ...; 
   ! the first two character denote 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: k, an, akappa, bn, bkappa, cn, ckappa, dn, dkappa
      character(len=50)   :: Xk_string
      !
      Xk_string = "X^"//trim(adjustl(pqn_string(k)))//" ("// orbital_name(an,akappa)//","//orbital_name(bn,bkappa)//";"// &
                                                             orbital_name(cn,ckappa)//","//orbital_name(dn,dkappa)//")"
      !
   end function Xk_name
   !
end module rabs_dirac_orbital
