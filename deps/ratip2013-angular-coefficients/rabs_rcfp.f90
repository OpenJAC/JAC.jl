module rabs_rcfp
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module provides a set of procedures which support the analytic integration over the angular coordinates for 
! matrix elements with symmetry-adapted subshell states. This includes the definition and set-up of the reduced 
! coefficients of fractional parentage as well as of reduced matrix elements for the tensor operators W^(k_q k_j) and
! T^{(k). Many of these coefficients are given by tables as defined in the header of this module. Currently, coefficients
! are provided for subshells with j <= 9/2 in jj-coupling.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_functions_string
   implicit none
   !
   public  :: rcfp_get_coefficient
                 ! Returns the value of the reduced coefficient of fractional parentage (rcfp) for subshells with 
                 ! j = 1/2, 3/2, 5/2, 7/2.
   public  :: rcfp_get_reduced_W
                 ! Returns the value of the reduced matrix elements W^(k_q k_j) for subshells with 
                 ! j = 1/2, 3/2, 5/2, 7/2, 9/2. 
   private :: rcfp_get_reduced_W_from_table
                 ! Returns the value of the reduced matrix elements W^(k_q k_j) for subshells with 
                 ! j = 1/2, 3/2, 5/2, 7/2 by looking up appropriate tables.
   private :: rcfp_get_reduced_W_all
                 ! Auxiliarity routine for rcfp_get_reduced_W. Returns the values of the reduced matrix elements 
                 ! of W^(00), W^(10), and W^(01) for j= 1/2, 3/2, 5/2, 7/2 from tables.
   private :: rcfp_get_reduced_W_three_half 
                 ! Auxiliary  routine for rcfp_get_reduced_W. Returns the the value of the reduced matrix elements 
                 ! of W^(12) and W^(03) for j= 3/2 from tables.
   private :: rcfp_get_reduced_W_five_half
                 ! Auxiliary routine for rcfp_get_reduced_W. Returns the value of the reduced matrix elements of 
                 ! W^(12), W^(03), W^(14), and W^(05) for j= 5/2 from tables.   
   private :: rcfp_get_reduced_W_seven_half
                 ! Auxiliary routine for rcfp_get_reduced_W. Returns the value of the reduced matrix elements of 
                 ! W^(12), W^(03), W^(14), W^(05), W^(16), and W^(07) for j= 7/2 from tables.
   public  :: rcfp_calculate_Wk_me
                 ! Calculate the value of the matrix element W^(k_j) subshells with j = 1/2, 3/2, 5/2, 7/2, 9/2. 
   public  :: rcfp_calculate_Wk_times_Wk_0_me
                 ! Calculate the value of the matrix element [ W^(k_j) * W^(k_j) ]^(0) for subshells with 
                 ! j= 1/2, 3/2, 5/2, 7/2, 9/2.
   public  :: rcfp_calculate_a_times_Wk_me
                 ! Calculates the value of the matrix element [ a^(q j)_m_q * W^(k_j) ]^(0) for subshells with 
		 ! j = 1/2, 3/2, 5/2, 7/2, 9/2.
   public  :: rcfp_calculate_Wk_times_a_me
                 ! Calculates the value of the matrix element [ W^(k_j) * a^(q j)_m_q ]^(0) for subshells with  
		 ! j = 1/2, 3/2, 5/2, 7/2, 9/2.
   public  :: rcfp_calculate_a_me
                 ! Calculates the value of the matrix element (j^N QJ :: a^(q j)_m_q :: j^N' QJ) for subshells with
                 ! j= 1/2, 3/2, 5/2, 7/2, 9/2. 
   public  :: rcfp_Q_space_delta
                 ! Calculate trivial delta factors in the Q space.
   public  :: rcfp_input
                 ! Collects and proceeds all input for the calculations of coefficients of fractional parentage 
                 ! (cfp and rcfp) and further reduced matrix elements for symmetry-adapted subshell states. 
   private :: rcfp_analyse
                 ! Analyzes a given input string for different angular momenta and other quantum numbers.
   public  :: rcfp_get_term_number
                 ! Returns an (internal) index for a  given subshell state.
   public  :: rcfp_Clebsch_Gordan_qusispin
                 ! Calculates specific the Clebsch-Gordan coefficient which need in quasispin formalism.
   !
   ! Define several derived data structures to keep information about subshell states and terms together.
   type, public  :: subshell_state
      integer :: state         ! State number of the subshell.
      integer :: n             ! Principal quantum number n.
      integer :: nq            ! Number of electrons in the subshell.
      integer :: subshellMQ    ! Subshell quasispin projection 2*M_Q.
   end type subshell_state
   !
   type, public :: subshell_term
      integer :: j             ! Angular momentum 2*j.
      integer :: Q             ! Subshell total quasispin 2*Q.
      integer :: nu            ! Seniority number.
      integer :: subshellJ     ! Subshell total angular momentum 2*J.
      integer :: Nr            ! State identifier Nr.
   end type subshell_term
   !
   type, public :: reduced_coeff
      integer :: phase         ! Weight factor.
      integer :: nom     
      integer :: denom   
   end type reduced_coeff
   !  
   ! Define all possible antisymetric subshell terms in jj-coupling
   type(subshell_term), dimension(1:63), parameter, public ::      &
      terms_jj = (/                                                &
      !
      ! j = 1/2; these terms have indices terms_jj(1:2)  
      subshell_term(1, 0, 1, 1, 0), subshell_term(1, 1, 0, 0, 0),  &
      !
      ! j = 3/2; these terms have indices terms_jj(3:5)
      subshell_term(3, 1, 1, 3, 0), subshell_term(3, 2, 0, 0, 0),  &
      subshell_term(3, 0, 2, 4, 0),                                &                
      !
      ! j = 5/2; these terms have indices terms_jj(6:11)
      subshell_term(5, 2, 1, 5, 0), subshell_term(5, 0, 3, 3, 0),  & 
      subshell_term(5, 0, 3, 9, 0), subshell_term(5, 3, 0, 0, 0),  &
      subshell_term(5, 1, 2, 4, 0), subshell_term(5, 1, 2, 8, 0),  &
      !
      ! j = 7/2; these terms have indices terms_jj(12:25)
      subshell_term(7, 1, 3, 3, 0), subshell_term(7, 1, 3, 5, 0),  &
      subshell_term(7, 3, 1, 7, 0), subshell_term(7, 1, 3, 9, 0),  &
      subshell_term(7, 1, 3,11, 0), subshell_term(7, 1, 3,15, 0),  &
      subshell_term(7, 4, 0, 0, 0), subshell_term(7, 2, 2, 4, 0),  &
      subshell_term(7, 0, 4, 4, 0), subshell_term(7, 2, 2, 8, 0),  &
      subshell_term(7, 0, 4, 8, 0), subshell_term(7, 0, 4,10, 0),  &
      subshell_term(7, 2, 2,12, 0), subshell_term(7, 0, 4,16, 0),  &
      !
      ! j = 9/2; these terms have indices terms_jj(26:63)
      subshell_term(9, 0, 5, 1, 0), subshell_term(9, 2, 3, 3, 0),  &
      subshell_term(9, 2, 3, 5, 0), subshell_term(9, 0, 5, 5, 0),  &
      subshell_term(9, 2, 3, 7, 0), subshell_term(9, 0, 5, 7, 0),  &
      subshell_term(9, 4, 1, 9, 0), subshell_term(9, 2, 3, 9, 0),  &
      subshell_term(9, 0, 5, 9, 0), subshell_term(9, 2, 3,11, 0),  &
      subshell_term(9, 0, 5,11, 0), subshell_term(9, 2, 3,13, 0),  &
      subshell_term(9, 0, 5,13, 0), subshell_term(9, 2, 3,15, 0),  &
      subshell_term(9, 0, 5,15, 0), subshell_term(9, 2, 3,17, 0),  &
      subshell_term(9, 0, 5,17, 0), subshell_term(9, 0, 5,19, 0),  &
      subshell_term(9, 2, 3,21, 0), subshell_term(9, 0, 5,25, 0),  &
      subshell_term(9, 5, 0, 0, 0), subshell_term(9, 1, 4, 0, 0),  &
      subshell_term(9, 3, 2, 4, 0), subshell_term(9, 1, 4, 4, 0),  &
      subshell_term(9, 1, 4, 6, 0), subshell_term(9, 3, 2, 8, 0),  &
      subshell_term(9, 1, 4, 8, 1), subshell_term(9, 1, 4, 8, 2),  &
      subshell_term(9, 1, 4,10, 0), subshell_term(9, 3, 2,12, 0),  &
      subshell_term(9, 1, 4,12, 1), subshell_term(9, 1, 4,12, 2),  &
      subshell_term(9, 1, 4,14, 0), subshell_term(9, 3, 2,16, 0),  &
      subshell_term(9, 1, 4,16, 0), subshell_term(9, 1, 4,18, 0),  &
      subshell_term(9, 1, 4,20, 0), subshell_term(9, 1, 4,24, 0) /)
   !	  	  
   ! Set-up the reduced coefficients of fractional parentage (rcfp)
   ! j = 1/2			 
   type(reduced_coeff), dimension(1,1), parameter, public ::              &
      rcfp_one_half = reshape(source = (/                                 &
      reduced_coeff(-1,      4,     1)                  /), shape=(/1,1/))
   ! j = 3/2
   type(reduced_coeff), dimension(1,2), parameter, public ::              &
      rcfp_three_half = reshape(source = (/                               &
      reduced_coeff(-1,     12,     1), reduced_coeff(-1,     20,     1)  &
	                                                /), shape=(/1,2/))
   ! j = 5/2	
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      rcfp_five_half = reshape(source = (/                                &
      reduced_coeff(-1,     24,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,     30,     1), &
      reduced_coeff( 1,    120,     7), reduced_coeff(-1,     90,     7), &
      reduced_coeff(-1,     54,     1), reduced_coeff(-1,     48,     7), &
      reduced_coeff( 1,    330,     7)                  /), shape=(/3,3/))
   ! j = 7/2	
   type(reduced_coeff), dimension(6,8), parameter, public ::              &
      rcfp_seven_half = reshape(source = (/                               &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     40,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     54,     7), reduced_coeff(-1,     33,     1), &
      reduced_coeff(-1,     40,     1), reduced_coeff(-1,     65,     7), &
      reduced_coeff( 1,     30,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 1,     88,     7), reduced_coeff(-1,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   1755,    77), &
      reduced_coeff(-1,     40,    11), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 1,    198,     7), reduced_coeff(-1,     36,    11), &
      reduced_coeff(-1,     72,     1), reduced_coeff( 1,   4500,    77), &
      reduced_coeff(-1,    234,    11), reduced_coeff(-1,    360,    11), &
      reduced_coeff(-1,     78,    35), reduced_coeff( 1,    156,     5), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    108,    91), &
      reduced_coeff(-1,     30,     1), reduced_coeff( 1,     96,    13), &
      reduced_coeff( 1,     66,     5), reduced_coeff( 1,     49,     5), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,     27,     1), &
      reduced_coeff( 1,    130,     7), reduced_coeff( 1,    136,     7), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    195,    11), &
      reduced_coeff(-1,    104,     1), reduced_coeff(-1,    245,    11), &
      reduced_coeff(-1,    624,    11), reduced_coeff( 1,   1224,    11), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   2720,   143), &
      reduced_coeff(-1,   2448,    77), reduced_coeff(-1,   7752,    91)  &
   	                                                /), shape=(/6,8/))
   ! j = 9/2		          
   type(reduced_coeff), dimension(20,3), parameter, public ::             &
      rcfp_nine_half_ket_46_48 = reshape(source = (/                      &
      ! rcfp for term_jj(46)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     60,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(47)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,     12,     1), &
      reduced_coeff(-1,      8,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(48)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     20,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,   1664,    33), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     50,     1), reduced_coeff(-1,    130,    33), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,    272,    11), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,    560,    11), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)  &
   	                                               /), shape=(/20,3/))
   !      
   type(reduced_coeff), dimension(20,3), parameter, public ::             &
      rcfp_nine_half_ket_49_51 = reshape(source = (/                      &
      ! rcfp for term_jj(49)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,    130,     7), reduced_coeff(-1,     48,     7), &
      reduced_coeff( 1,     44,    21), reduced_coeff(-1,    340,    77), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    140,    33), &
      reduced_coeff(-1,     70,    11), reduced_coeff( 1,   3808,   143), &
      reduced_coeff(-1,   2280,   143), reduced_coeff(-1,    110,    13), &
      reduced_coeff(-1,    918,   143), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(50)
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,     63,     5), &
      reduced_coeff(-1,    312,    55), reduced_coeff(-1,      4,    11), &
      reduced_coeff( 1,      9,    11), reduced_coeff( 1,    255,    11), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   5040,   143), &
      reduced_coeff( 1,    840,   143), reduced_coeff(-1,   1428,   143), &
      reduced_coeff(-1,     95,   143), reduced_coeff( 1,    504,   715), &
      reduced_coeff(-1,   2856,   143), reduced_coeff( 1,  13566,   715), &
      reduced_coeff(-1,    850,   143), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(51)
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,    384,    11), &
      reduced_coeff( 1,    156,    11), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,   1920,   143), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     90,     1), reduced_coeff( 1,   7350,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   8160,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   1512,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   3648,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   9000,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)  &
   	                                               /), shape=(/20,3/))
   !      
   type(reduced_coeff), dimension(20,3), parameter, public ::             &
      rcfp_nine_half_ket_52_54 = reshape(source = (/                      &
      ! rcfp for term_jj(52)
      reduced_coeff(-1,   2184,   253), reduced_coeff(-1,     63,    23), &
      reduced_coeff(-1,  59904, 19481), reduced_coeff(-1, 302460, 19481), &
      reduced_coeff( 1,   5265,   161), reduced_coeff( 1, 848691,253253), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1, 145152, 36179), &
      reduced_coeff( 1, 217728, 36179), reduced_coeff( 1,1049580, 36179), &
      reduced_coeff( 1, 287337, 36179), reduced_coeff( 1,   5184,   299), &
      reduced_coeff( 1, 261120, 36179), reduced_coeff( 1, 691866, 36179), &
      reduced_coeff(-1,    750, 36179), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,  76608,  3289), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(53)
      reduced_coeff( 1,   1224,  1265), reduced_coeff( 1,   2652,   115), &
      reduced_coeff( 1, 188598, 13915), reduced_coeff(-1,  31824,  2783), &
      reduced_coeff( 1,    204,    23), reduced_coeff(-1,  12996, 13915), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,  25500,  2783), &
      reduced_coeff(-1,  38250,  2783), reduced_coeff(-1,    768,  2783), &
      reduced_coeff( 1,  81396, 13915), reduced_coeff( 1,   3213,   115), &
      reduced_coeff(-1,  60543,  2783), reduced_coeff(-1,3066144,236555), &
      reduced_coeff( 1, 727776, 47311), reduced_coeff( 1,    207,    17), &
      reduced_coeff(-1,  41553, 21505), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(54)
      reduced_coeff( 1,     52,     5), reduced_coeff(-1,     36,     5), &
      reduced_coeff( 1,     84,     5), reduced_coeff(-1,     56,    13), &
      reduced_coeff( 1,    126,    13), reduced_coeff(-1,   3570,   325), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    360,    13), &
      reduced_coeff( 1,     60,    13), reduced_coeff( 1,    102,    13), &
      reduced_coeff(-1,   1064,    65), reduced_coeff(-1,   2142,    65), &
      reduced_coeff( 1,     42,    13), reduced_coeff( 1,    228,    65), &
      reduced_coeff( 1,    252,    13), reduced_coeff( 1,    342,    13), &
      reduced_coeff(-1,     66,    65), reduced_coeff( 1,    230,    13), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)  &
   	                                               /), shape=(/20,3/))
   !      
   type(reduced_coeff), dimension(20,3), parameter, public ::             &
      rcfp_nine_half_ket_55_57 = reshape(source = (/                      &
      ! rcfp for term_jj(55)
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    144,    11), &
      reduced_coeff( 1,    416,    11), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,     32,   165), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,    130,     1), reduced_coeff(-1,   1922,    33), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    896,    55), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    476,    11), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   1344,    11), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   2052,    55), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,    308,     5), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(56)
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    132,     5), &
      reduced_coeff(-1,  52728,  6655), reduced_coeff( 1,    196,  1331), &
      reduced_coeff(-1,     50,    11), reduced_coeff(-1,  24990,  1331), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,  21160,  1331), &
      reduced_coeff( 1,  31740,  1331), reduced_coeff( 1,  26250,  1331), &
      reduced_coeff( 1,  12920,  1331), reduced_coeff(-1,    357,    55), &
      reduced_coeff( 1,   9583,  1331), reduced_coeff(-1, 344988,  6655), &
      reduced_coeff( 1,   5700,  1331), reduced_coeff( 1,    171,    11), &
      reduced_coeff(-1,     15,   121), reduced_coeff(-1,   4830,   121), &
      reduced_coeff(-1,     84,    11), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(57)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1, 209950,  9317), reduced_coeff( 1,  77520,  9317), &
      reduced_coeff( 1,  12920,   231), reduced_coeff(-1,  25688,  9317), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   4522,  3993), &
      reduced_coeff( 1,   2261,  1331), reduced_coeff(-1,  48640, 22627), &
      reduced_coeff(-1, 285144,  9317), reduced_coeff( 1,    931,    44), &
      reduced_coeff( 1, 273885, 90508), reduced_coeff(-1, 112908,429913), &
      reduced_coeff(-1,2138580,158389), reduced_coeff( 1, 137781,  3740), &
      reduced_coeff( 1,6654375,156332), reduced_coeff(-1,  59616, 39083), &
      reduced_coeff( 1, 284089, 17765), reduced_coeff( 0,      1,     1)  &
   	                                               /), shape=(/20,3/)) 
   !     
   type(reduced_coeff), dimension(20,3), parameter, public ::             &
      rcfp_nine_half_ket_58_60 = reshape(source = (/                      &
      ! rcfp for term_jj(58)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 1,   1530,    77), reduced_coeff( 1,  13056,  1001), &
      reduced_coeff( 1,  29376,  1001), reduced_coeff( 1,    720,  1001), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   1890,   143), &
      reduced_coeff(-1,    315,   143), reduced_coeff( 1, 101124,  2431), &
      reduced_coeff(-1,   4560,  1001), reduced_coeff(-1,  13965,   572), &
      reduced_coeff(-1,  35131,  9724), reduced_coeff( 1,  13500,  2431), &
      reduced_coeff(-1, 685900, 17017), reduced_coeff(-1,   1197,   884), &
      reduced_coeff(-1,  28875,   884), reduced_coeff(-1,   5060,   221), &
      reduced_coeff( 1,    759,    17), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(59)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 1,  22848,   715), reduced_coeff( 0,      1,     1), &
      reduced_coeff(-1,    170,     1), reduced_coeff( 1,    918,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,  32832,   715), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   9044,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,    576,    13), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   7524,    65), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 1,   1012,     5), reduced_coeff( 0,      1,     1), &
      ! rcfp for term_jj(60)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,   2128,   143), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   1938,   143), &
      reduced_coeff( 1,   2907,   143), reduced_coeff(-1,   4860,   143), &
      reduced_coeff(-1,  17136,  2717), reduced_coeff(-1,   1309,    52), &
      reduced_coeff(-1,   8505,   572), reduced_coeff( 1,  15876,   247), &
      reduced_coeff( 1,    420,    13), reduced_coeff( 1,   1287,    20), &
      reduced_coeff(-1,   6075,   988), reduced_coeff(-1,    132,    19), &
      reduced_coeff(-1,    253,    95), reduced_coeff(-1,    650,    19)  &
   	                                               /), shape=(/20,3/)) 
   !     
   type(reduced_coeff), dimension(20,3), parameter, public ::             &
      rcfp_nine_half_ket_61_63 = reshape(source = (/                      &
      ! rcfp for term_jj(61)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 1,    570,    13), &
      reduced_coeff( 1,     95,    13), reduced_coeff( 1,   4104,   221), &
      reduced_coeff( 1,    504,    65), reduced_coeff(-1,   1463,   260), &
      reduced_coeff(-1,  39501,   884), reduced_coeff(-1,  60516,  1105), &
      reduced_coeff(-1,   1596,   221), reduced_coeff( 1,   3933,    68), &
      reduced_coeff( 1,    621,  3740), reduced_coeff(-1,    840,    17), &
      reduced_coeff(-1,    805,    17), reduced_coeff( 1,    390,    11), &
      ! rcfp for term_jj(62)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,   5796,   221), &
      reduced_coeff( 1,  17664,  1235), reduced_coeff(-1,   5313,    65), &
      reduced_coeff( 1,   1771,   221), reduced_coeff(-1,  16632,  1615), &
      reduced_coeff(-1,     88,    17), reduced_coeff( 1,    693,    17), &
      reduced_coeff(-1,  94269,  1615), reduced_coeff( 1, 192500,  7429), &
      reduced_coeff( 1,  30030,   323), reduced_coeff( 1,  24570,   437), &
      ! rcfp for term_jj(63)
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), &
      reduced_coeff( 0,      1,     1), reduced_coeff(-1,  15000,   323), &
      reduced_coeff( 1,    280,    17), reduced_coeff(-1,   1170,    17), &
      reduced_coeff(-1,  48750,  3553), reduced_coeff(-1,  15600,   437), &
      reduced_coeff( 1,   3510,    19), reduced_coeff(-1,  33930,   253)  &
   	                                               /), shape=(/20,3/))
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(00) :: j QJ)  for j = 1/2, 3/2, 5/2, 7/2.
   ! j = 1/2						   
   type(reduced_coeff), dimension(1:2), parameter, public ::              &
      W_00_one_half =    (/                                               &
      reduced_coeff(-1,     2,     1), reduced_coeff(-1,       2,     1)/)
   ! j = 3/2						   
   type(reduced_coeff), dimension(1:3), parameter, public ::              &
      W_00_three_half =  (/                                               &
      reduced_coeff(-1,    16,     1), reduced_coeff(-1,       6,     1), &     
      reduced_coeff(-1,    10,     1)                                   /)
   ! j = 5/2						   
   type(reduced_coeff), dimension(1:6), parameter, public ::              &
      W_00_five_half =   (/                                               &
      reduced_coeff(-1,    54,     1), reduced_coeff(-1,      12,     1), &     
      reduced_coeff(-1,    30,     1), reduced_coeff(-1,      12,     1), &     
      reduced_coeff(-1,    30,     1), reduced_coeff(-1,      54,     1)/)
   ! j = 7/2						   
   type(reduced_coeff), dimension(1:14), parameter, public ::             &
      W_00_seven_half =  (/                                               &
      reduced_coeff(-1,    32,     1), reduced_coeff(-1,      48,     1), &     
      reduced_coeff(-1,   128,     1), reduced_coeff(-1,      80,     1), &
      reduced_coeff(-1,    96,     1), reduced_coeff(-1,     128,     1), &	  
      reduced_coeff(-1,    20,     1), reduced_coeff(-1,      60,     1), &	  
      reduced_coeff(-1,    20,     1), reduced_coeff(-1,     108,     1), &
      reduced_coeff(-1,    36,     1), reduced_coeff(-1,      44,     1), &	       
      reduced_coeff(-1,   156,     1), reduced_coeff(-1,      68,     1)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(10) :: j QJ)  for j = 1/2, 3/2, 5/2, 7/2.
   ! j = 1/2						    
   type(reduced_coeff), dimension(1:2), parameter, public ::              &
      W_10_one_half =    (/                                               &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,       6,     1)/)
   ! j = 3/2						   
   type(reduced_coeff), dimension(1:3), parameter, public ::              &
      W_10_three_half =  (/                                               &
      reduced_coeff(-1,    12,     1), reduced_coeff(-1,      12,     1), &     
      reduced_coeff( 0,     1,     1)                                   /)
   ! j = 5/2						   
   type(reduced_coeff), dimension(1:6), parameter, public ::              &
      W_10_five_half =   (/                                               &
      reduced_coeff(-1,    48,     1), reduced_coeff( 0,       1,     1), &     
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      20,     1), &     
      reduced_coeff(-1,    10,     1), reduced_coeff(-1,      18,     1)/)
   ! j = 7/2						   
   type(reduced_coeff), dimension(1:14), parameter, public ::             &
      W_10_seven_half =  (/                                               &
      reduced_coeff(-1,     6,     1), reduced_coeff(-1,       9,     1), &     
      reduced_coeff(-1,   120,     1), reduced_coeff(-1,      15,     1), &
      reduced_coeff(-1,    18,     1), reduced_coeff(-1,      24,     1), &	  
      reduced_coeff(-1,    30,     1), reduced_coeff(-1,      30,     1), &	  
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      54,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &	       
      reduced_coeff(-1,    78,     1), reduced_coeff( 0,       1,     1)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(01) :: j QJ)  for j = 1/2, 3/2, 5/2, 7/2.
   ! j = 1/2						   
   type(reduced_coeff), dimension(1:2), parameter, public ::              &
      W_01_one_half =    (/ &
      reduced_coeff(-1,     6,     1), reduced_coeff( 0,       1,     1)/)
   ! j = 3/2						   
   type(reduced_coeff), dimension(1:3), parameter, public ::              &
      W_01_three_half =  (/                                               &
      reduced_coeff(-1,    12,     1), reduced_coeff( 0,       1,     1), &     
      reduced_coeff(-1,    12,     1)                                   /)
   ! j = 5/2						   
   type(reduced_coeff), dimension(1:6), parameter, public ::              &
      W_01_five_half =   (/                                               &
      reduced_coeff(-1,   126,     7), reduced_coeff(-1,      12,     7), &     
      reduced_coeff(-1,   198,     7), reduced_coeff( 0,       1,     1), &     
      reduced_coeff(-1,    48,     7), reduced_coeff(-1,     288,     7)/)
   ! j = 7/2						   
   type(reduced_coeff), dimension(1:14), parameter, public ::             &
      W_01_seven_half =  (/                                               &
      reduced_coeff(-1,    10,     7), reduced_coeff(-1,       5,     1), &
      reduced_coeff(-1,   168,     7), reduced_coeff(-1,     165,     7), &
      reduced_coeff(-1,   286,     7), reduced_coeff(-1,     680,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      30,     7), &
      reduced_coeff(-1,    10,     7), reduced_coeff(-1,     180,     7), &
      reduced_coeff(-1,    60,     7), reduced_coeff(-1,     110,     7), &
      reduced_coeff(-1,   546,     7), reduced_coeff(-1,     408,     7)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(12) :: j Q'J')  for j = 3/2, 5/2, 7/2.
   ! j = 3/2  for odd seniority			     
   type(reduced_coeff), dimension(1,1), parameter, public ::              &
      W_12_three_half_odd = reshape(source =   (/                         &
      reduced_coeff( 1,    60,     1)                   /), shape=(/1,1/))
   ! j = 3/2  for even seniority			     
   type(reduced_coeff), dimension(2,2), parameter, public ::              &
      W_12_three_half_even = reshape(source =  (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,      30,     1), &
      reduced_coeff(-1,    30,     1), reduced_coeff( 0,       1,     1)  &
                                                        /), shape=(/2,2/))
   ! j = 5/2  for odd seniority			     
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_12_five_half_odd = reshape(source =    (/                         &
      reduced_coeff( 1,   420,     7), reduced_coeff( 1,     360,     7), &
      reduced_coeff( 1,   270,     7), reduced_coeff( 1,     360,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,   270,     7), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1)                   /), shape=(/3,3/))
   ! j = 5/2  for even seniority			     
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_12_five_half_even = reshape(source =   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,    1960,    49), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    1960,    49), &
      reduced_coeff(-1,  1000,    49), reduced_coeff( 1,    2430,    49), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,    2430,    49), &
      reduced_coeff( 1,  1980,    49)                   /), shape=(/3,3/))
   ! j = 7/2  for odd seniority			     
   type(reduced_coeff), dimension(21), parameter, public ::               &
      W_12_seven_half_odd =                    (/                         &
      reduced_coeff(-1,   252,    10), reduced_coeff( 1,    1056,    70), &
      reduced_coeff( 1,   144,     7), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   507,    10), reduced_coeff(-1,      88,     1), &
      reduced_coeff(-1,   325,    14), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     200,     3), &
      reduced_coeff(-1,   520,    21), reduced_coeff( 1,      80,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    6125,   462), &
      reduced_coeff(-1,   560,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   390,   539), reduced_coeff(-1,    3840,    49), &
      reduced_coeff( 1,  2040,    49)                                   /)
   ! j = 7/2  for even seniority			     
   type(reduced_coeff), dimension(36), parameter, public ::               &
      W_12_seven_half_even =                   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      50,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,  1280,    49), reduced_coeff( 1,     990,    49), &
      reduced_coeff( 1,  2640,    49), reduced_coeff( 1,    1950,    49), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,   480,    49), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,     360,   539), &
      reduced_coeff( 1,  1872,    49), reduced_coeff( 1,      42,     1), &
      reduced_coeff( 1,   390,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,    48,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     234,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,    1040,    11), &
      reduced_coeff( 1,   340,     7), reduced_coeff( 0,       1,     1)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(03) :: j Q'J')  for j = 3/2, 5/2, 7/2.
   ! j = 3/2  for odd seniority
   type(reduced_coeff), dimension(1,1), parameter, public ::              &
      W_03_three_half_odd = reshape(source =   (/                         &
      reduced_coeff(-1,    28,     1)                   /), shape=(/1,1/))
   ! j = 3/2  for even seniority			     
   type(reduced_coeff), dimension(2,2), parameter, public ::              &
      W_03_three_half_even = reshape(source =  (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,      28,     1)  &
                                                        /), shape=(/2,2/))
   ! j = 5/2  for odd seniority			         
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_03_five_half_odd = reshape(source =    (/                         &
      reduced_coeff(-1,   882,    21), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   384,    21), reduced_coeff(-1,     400,    21), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     400,    21), &    
      reduced_coeff( 1,   286,    21)                   /), shape=(/3,3/)) 
   ! j = 5/2  for even seniority			     
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_03_five_half_even = reshape(source =   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   162,     7), reduced_coeff(-1,     300,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,     300,     7), &
      reduced_coeff( 1,    22,     7)                   /), shape=(/3,3/)) 
   ! j = 7/2  for odd seniority			        
   type(reduced_coeff), dimension(21), parameter, public ::               &
      W_03_seven_half_odd =                    (/                         & 
      reduced_coeff(-1,  1188,    70), reduced_coeff(-1,     196,    10), & 
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,     234,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   189,   110), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,  1911,   242), reduced_coeff( 1,    1470,   121), & 
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      56,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,  394805, 22022), &
      reduced_coeff( 1,  5250,   121), reduced_coeff( 1,   53760,  1573), &
      reduced_coeff(-1,    78,   847), reduced_coeff( 1,   17408,   847), &
      reduced_coeff( 1, 12920,  1001)                                   /) 
   ! j = 7/2  for even seniority			       
   type(reduced_coeff), dimension(36), parameter, public ::               &
      W_03_seven_half_even =                   (/                         & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   110,     7), reduced_coeff( 0,       1,     1), & 
      reduced_coeff(-1,   240,     7), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    1920,    77), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,      52,     7), &
      reduced_coeff(-1,   224,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,   32490,   847), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,  7644,   121), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,    12,    70), reduced_coeff( 1,     224,    10), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,    52,    70), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,  2040,    77), reduced_coeff(-1,     364,   121), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,    1292,    77)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(14) :: j Q'J')  for j = 5/2, 7/2.
   ! j = 5/2  for odd seniority			         
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_14_five_half_odd = reshape(source =    (/                         &
      reduced_coeff( 1,  3780,    35), reduced_coeff(-1,     720,    35), &
      reduced_coeff(-1,  4950,    35), reduced_coeff(-1,     720,    35), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  4950,    35), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1)                   /), shape=(/3,3/))
   ! j = 5/2  for even seniority			     
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_14_five_half_even = reshape(source =   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  3528,    49), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  2430,    49), reduced_coeff( 1,    1980,    49), &
      reduced_coeff(-1,  3528,    49), reduced_coeff( 1,    1980,    49), &
      reduced_coeff(-1,  7722,    49)                   /), shape=(/3,3/))
   ! j = 7/2  for odd seniority			         
   type(reduced_coeff), dimension(21), parameter, public ::               &
      W_14_seven_half_odd =                    (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,      54,     7), &
      reduced_coeff(-1,   528,     7), reduced_coeff( 1,     546,    11), &
      reduced_coeff(-1,   378,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,  4335,   110), reduced_coeff(-1,      96,    11), &
      reduced_coeff( 1, 11271,  1694), reduced_coeff( 1,    3822,   121), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     120,     1), &
      reduced_coeff( 1, 12000,    77), reduced_coeff(-1,     624,    11), &
      reduced_coeff(-1,   960,    11), reduced_coeff(-1,   30345,  3146), &
      reduced_coeff(-1,   210,   121), reduced_coeff(-1,  228480,  1573), &
      reduced_coeff( 1,580476,  5929), reduced_coeff( 1,  146880,  5929), &
      reduced_coeff(-1,627912,  7007)                                   /)
      ! j = 7/2  for even seniority			       
   type(reduced_coeff), dimension(36), parameter, public ::               &
      W_14_seven_half_even =                   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      90,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  2640,    49), reduced_coeff( 1,     480,    49), &
      reduced_coeff(-1,   360,   539), reduced_coeff( 1,   20592,   539), &
      reduced_coeff(-1,    42,     1), reduced_coeff( 1,     390,    11), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,468180,  5929), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,   21840,   847), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,  359424,  5929), &
      reduced_coeff(-1,  6750,   539), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,917280,  5929), reduced_coeff(-1,   10710,   121), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,    36,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,     858,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    5304,   121), &
      reduced_coeff(-1, 69768,   847), reduced_coeff( 0,       1,     1)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(05) :: j Q'J')  for j= 5/2, 7/2.
   ! j = 5/2  for odd seniority			         
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_05_five_half_odd = reshape(source =    (/                         &
      reduced_coeff(-1,  1386,    21), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,     440,    21), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     440,    21), & 
      reduced_coeff(-1,  1430,    21)                   /), shape=(/3,3/))
      ! j = 5/2  for even seniority			     
   type(reduced_coeff), dimension(3,3), parameter, public ::              &
      W_05_five_half_even = reshape(source =   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     330,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     330,     7), &
      reduced_coeff( 1,   572,     7)                   /), shape=(/3,3/))
      ! j = 7/2  for odd seniority			         
   type(reduced_coeff), dimension(21), parameter, public ::               &
      W_05_seven_half_odd =                    (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     144,     7), &
      reduced_coeff( 1,    14,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   975,    22), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1, 50421,  2002), reduced_coeff( 1,     336,    11), & 
      reduced_coeff(-1,  7000,   143), reduced_coeff(-1,      88,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    6845, 26026), &
      reduced_coeff(-1,  3360,   143), reduced_coeff( 1,   14280,  1859), &
      reduced_coeff(-1,  1836,    77), reduced_coeff(-1,  103360,  1001), &
      reduced_coeff(-1, 28424,143143)                                   /) 
      ! j = 7/2  for even seniority			       
   type(reduced_coeff), dimension(36), parameter, public ::               &
      W_05_seven_half_even =                   (/                         & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &  	      
      reduced_coeff( 1,   390,     7), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,      70,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,      32,     7), &
      reduced_coeff( 1,    14,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,     576,    77), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,   210,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1, 63888,  1183), reduced_coeff(-1,     154,    13), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    8568,   169), &
      reduced_coeff(-1,   176,     7), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  1938,    91), reduced_coeff( 1,    1088,    11), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,   28424,  1183)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(16) :: j Q'J')  for j= 7/2.    
   ! j = 7/2  for odd seniority			       
   type(reduced_coeff), dimension(21), parameter, public ::               &
      W_16_seven_half_odd =                    (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     576,    11), &
      reduced_coeff( 1,   390,    77), reduced_coeff(-1,     144,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,     520,    11), &
      reduced_coeff(-1,    49,   121), reduced_coeff(-1,   12480,   121), &  	      
      reduced_coeff(-1,   408,    11), reduced_coeff( 1,     520,     3), & 
      reduced_coeff(-1,  1960,    33), reduced_coeff(-1,    1664,    11), &
      reduced_coeff( 1,  3264,    11), reduced_coeff( 1,  552250,  4719), &
      reduced_coeff(-1, 43520,   847), reduced_coeff(-1,   38760, 11011), &
      reduced_coeff(-1,  2652,   121), reduced_coeff( 1,   15504,   121), &
      reduced_coeff( 1, 38760,   143)                                   /) 
   ! j = 7/2  for even seniority			       
   type(reduced_coeff), dimension(36), parameter, public ::               &
      W_16_seven_half_even =                   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1,   130,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &  	      
      reduced_coeff( 1,   390,    11), reduced_coeff(-1,      48,     1), & 
      reduced_coeff(-1,   234,     7), reduced_coeff( 1,    1040,    11), & 
      reduced_coeff( 1,   340,     7), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  3120,   121), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,    1170,   121), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,   18720,   121), &
      reduced_coeff(-1,    36,    11), reduced_coeff( 1,     858,     7), &
      reduced_coeff(-1,  5304,   121), reduced_coeff(-1,   69768,   847), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  1020,    11), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff(-1,   33592,   121), &
      reduced_coeff( 1, 31654,   121), reduced_coeff( 0,       1,     1)/)
   !
   ! Define the value of reduced matrix elements 
   ! ( j QJ :: W(07) :: j Q'J')  for j= 7/2.    
   ! j = 7/2 for odd seniority			         
   type(reduced_coeff), dimension(21), parameter, public ::               &
      W_07_seven_half_odd =                    (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,   162,     7), reduced_coeff( 1,     272,     7), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1, 11025,  1573), reduced_coeff( 1,    4624,   121), & 
      reduced_coeff(-1,  1632,   143), reduced_coeff(-1,     120,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1, 1224510, 20449), &
      reduced_coeff( 1,306000, 11011), reduced_coeff( 1,12558240,143143), &
      reduced_coeff(-1,  6460,   121), reduced_coeff( 1,   77520,  1573), &
      reduced_coeff(-1,297160,  1859)                                   /) 
   ! j = 7/2  for even seniority		       
   type(reduced_coeff), dimension(36), parameter, public ::               &
      W_07_seven_half_even =                   (/                         &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,      60,     1), & 
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  1600,    77), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  2040,    77), reduced_coeff( 1,    4410,   121), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1, 18360,   121), reduced_coeff( 0,       1,     1), &
      reduced_coeff(-1, 18816,   845), reduced_coeff(-1,   11016,   455), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,   11628,  1183), &
      reduced_coeff( 1,    34,     5), reduced_coeff( 0,       1,     1), &
      reduced_coeff( 1,  7752,   143), reduced_coeff( 1,    9690,   121), &
      reduced_coeff( 0,     1,     1), reduced_coeff( 1,  222870,  1859)/)
   !
   ! Define the minimal and maximal limits of the subshell terms for 
   ! odd number operators in second quantization 
   integer, dimension(63), parameter, private :: rcfp_min_odd = (/        &
       2, 1, 4, 3, 3, 9, 9, 9, 6, 6, 6,18,18,18,18,18,18,12,12,12,        &
      12,12,12,12,12,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,        &
      46,46,46,46,46,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,        &
      26,26,26 /) 
   ! 
   integer, dimension(63), parameter, private :: rcfp_max_odd = (/        &
       2, 1, 5, 3, 3,11,11,11, 8, 8, 8,25,25,25,25,25,25,17,17,17,        &
      17,17,17,17,17,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,        &
      63,63,63,63,63,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,        &
      45,45,45 /)
   !
   ! Define the minimal and maximal limits of the subshell terms for 
   ! even number operators in second quantization 
   integer, dimension(63), parameter, private :: rcfp_min_even =(/        &
       1, 2, 3, 4, 4, 6, 6, 6, 9, 9, 9,12,12,12,12,12,12,18,18,18,        &
      18,18,18,18,18,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,        &
      26,26,26,26,26,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,        &
      46,46,46 /)
   !
   integer, dimension(63), parameter, private :: rcfp_max_even =(/        &
       1, 2, 3, 5, 5, 8, 8, 8,11,11,11,17,17,17,17,17,17,25,25,25,        &
      25,25,25,25,25,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,        &
      45,45,45,45,45,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,        &
      63,63,63 /)
   !
contains
   !
   function rcfp_get_coefficient(no_bra,no_ket)                   result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the reduced coefficient of fractional parentage 
   ! ( j QJ :: a^(qj) :: j Q'J') for subshells with j = 1/2, 3/2, 5/2, 7/2.
   !
   ! Calls: 
   !-----------------------------------------------------------------------
      !
      integer, intent(in) :: no_bra, no_ket
      real(kind=dp)       :: coeff
      !   
      integer             :: no_a, no_b, phase, nom, denom
      !
      coeff = zero
      !   
      if (rcfp_min_odd(no_bra)  /=  rcfp_min_even(no_ket)               .or. &
          triangle(terms_jj(no_ket)%Q +1,2,terms_jj(no_bra)%Q +1) == 0  .or. &
          triangle(terms_jj(no_bra)%subshellJ +1,terms_jj(no_bra)%j +1,      &
                   terms_jj(no_ket)%subshellJ +1) == 0) return
      !
      if (no_bra <= no_ket) then
         no_a = no_bra;   no_b = no_ket
      else
         no_a = no_ket;   no_b = no_bra
      end if
      !
      select case (terms_jj(no_bra)%j)
      case (1) 
         ! j=1/2 
         phase = rcfp_one_half(no_a, no_b-1)%phase
	 nom   = rcfp_one_half(no_a, no_b-1)%nom
	 denom = rcfp_one_half(no_a, no_b-1)%denom
      case (3)
         ! j=3/2 
         phase = rcfp_three_half(no_a-2, no_b-3)%phase
	 nom   = rcfp_three_half(no_a-2, no_b-3)%nom
	 denom = rcfp_three_half(no_a-2, no_b-3)%denom
      case (5)
         ! j=5/2 
         phase = rcfp_five_half(no_a-5, no_b-8)%phase
	 nom   = rcfp_five_half(no_a-5, no_b-8)%nom
	 denom = rcfp_five_half(no_a-5, no_b-8)%denom
      case (7)
         ! j=7/2 
         phase = rcfp_seven_half(no_a-11, no_b-17)%phase
	 nom   = rcfp_seven_half(no_a-11, no_b-17)%nom
	 denom = rcfp_seven_half(no_a-11, no_b-17)%denom
      case (9)
         ! j=9/2 
	 if (no_b > 45 .and. no_b < 49) then
            phase = rcfp_nine_half_ket_46_48(no_a-25,no_b-45)%phase
	    nom   = rcfp_nine_half_ket_46_48(no_a-25,no_b-45)%nom
	    denom = rcfp_nine_half_ket_46_48(no_a-25,no_b-45)%denom
 	 else if (no_b > 48 .and. no_b < 52) then
            phase = rcfp_nine_half_ket_49_51(no_a-25,no_b-48)%phase
	    nom   = rcfp_nine_half_ket_49_51(no_a-25,no_b-48)%nom
	    denom = rcfp_nine_half_ket_49_51(no_a-25,no_b-48)%denom
 	 else if (no_b > 51 .and. no_b < 55) then
            phase = rcfp_nine_half_ket_52_54(no_a-25,no_b-51)%phase
	    nom   = rcfp_nine_half_ket_52_54(no_a-25,no_b-51)%nom
	    denom = rcfp_nine_half_ket_52_54(no_a-25,no_b-51)%denom
 	 else if (no_b > 54 .and. no_b < 58) then
            phase = rcfp_nine_half_ket_55_57(no_a-25,no_b-54)%phase
	    nom   = rcfp_nine_half_ket_55_57(no_a-25,no_b-54)%nom
	    denom = rcfp_nine_half_ket_55_57(no_a-25,no_b-54)%denom
 	 else if (no_b > 57 .and. no_b < 61) then
            phase = rcfp_nine_half_ket_58_60(no_a-25,no_b-57)%phase
	    nom   = rcfp_nine_half_ket_58_60(no_a-25,no_b-57)%nom
	    denom = rcfp_nine_half_ket_58_60(no_a-25,no_b-57)%denom
 	 else if (no_b > 60 .and. no_b < 64) then
            phase = rcfp_nine_half_ket_61_63(no_a-25,no_b-60)%phase
	    nom   = rcfp_nine_half_ket_61_63(no_a-25,no_b-60)%nom
	    denom = rcfp_nine_half_ket_61_63(no_a-25,no_b-60)%denom
	 else if (rabs_use_stop) then 
            stop  "rcfp_get_coefficient(): program stop A."
	 end if
      case default
         stop  "rcfp_get_coefficient(): program stop B."
      end select
      !
      if (phase /= 0) then
         coeff = phase * sqrt(real(nom,kind=dp)/denom)
         if (no_bra >= no_ket  .and.                                       &
	     mod(terms_jj(no_bra)%Q - terms_jj(no_ket)%Q                 + &
	         terms_jj(no_bra)%subshellJ - terms_jj(no_ket)%subshellJ + &
	         terms_jj(no_bra)%j - 1, 4) /= 0) then
            coeff = - coeff
         end if
      end if
      !
   end function rcfp_get_coefficient
   !
   !
   function rcfp_get_reduced_W(no_bra,no_ket,k_q,k_j)              result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the reduced matrix element
   ! ( j QJ :: W^(k_q k_j) :: j Q'J') for subshells with j= 1/2, 3/2, 5/2,
   ! 7/2, and 9/2. 
   !
   ! Calls: 
   !-----------------------------------------------------------------------
      !
      integer, intent(in) :: no_bra, no_ket, k_q, k_j
      real(kind=dp)       :: coeff
      !
      integer             :: delta_J, delta_Q, no_run
      real(kind=dp)       :: phase
      !   
      coeff = zero
      !   
      if (terms_jj(no_bra)%j < 9) then
         coeff = rcfp_get_reduced_W_from_table(no_bra,no_ket,k_q,k_j)	
      else if (terms_jj(no_bra)%j == 9) then
         if (rcfp_min_even(no_bra)  /=  rcfp_min_even(no_ket)) return
         if (triangle(terms_jj(no_ket)%Q +1, 2*k_q +1,              &
	              terms_jj(no_bra)%Q +1) == 0)             return      
         if (triangle(terms_jj(no_bra)%subshellJ +1, 2*k_j +1,      &
                      terms_jj(no_ket)%subshellJ +1) == 0)     return
         do no_run = rcfp_min_odd(no_bra),rcfp_max_odd(no_bra)
	    delta_J = wigner_6j_triangle(terms_jj(no_bra)%j,         &
	                                 terms_jj(no_bra)%j,2*k_j,   &
		                         terms_jj(no_ket)%subshellJ, &
		                         terms_jj(no_bra)%subshellJ, & 
	                                 terms_jj(no_run)%subshellJ)
            if (delta_J /= 0 ) then
	       delta_Q = wigner_6j_triangle(1,1,2*k_q,terms_jj(no_ket)%Q, &
	                           terms_jj(no_bra)%Q,terms_jj(no_run)%Q)
	       if (delta_Q /= 0) then
	          coeff = coeff + rcfp_get_coefficient(no_bra,no_run)  * &
		                  rcfp_get_coefficient(no_run,no_ket)  * &
			  wigner_6j_symbol(terms_jj(no_bra)%j,           &
		                           terms_jj(no_bra)%j,2*k_j,     &
			                   terms_jj(no_ket)%subshellJ,   &
	                                   terms_jj(no_bra)%subshellJ,   &
				           terms_jj(no_run)%subshellJ,   &
					   .true.)                     * &
			  wigner_6j_symbol(1,1,2*k_q,terms_jj(no_ket)%Q, &
	                                   terms_jj(no_bra)%Q,           &
				           terms_jj(no_run)%Q,.true.)
	       end if
	    end if
	 end do
         !
	 if (abs(coeff) > eps10) then 
            coeff = coeff * sqrt((2 * k_q + one) * (2 * k_j + one))
	    phase = mod(terms_jj(no_bra)%Q + terms_jj(no_bra)%subshellJ + &
	                terms_jj(no_ket)%Q + terms_jj(no_ket)%subshellJ + &
		        2 * k_q + 2 * k_j, 4)
            if (rabs_use_stop         .and.  &
               (phase == 1  .or.  phase == 3)) then
               stop  "rcfp_get_reduced_W(): program stop A."
	    else if (phase == 2) then
	       coeff = - coeff
            end if
	 end if
      end if
      !
   end function rcfp_get_reduced_W
   !
   !
   function rcfp_get_reduced_W_from_table(no_bra,no_ket,k1,k2)     result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the reduced matrix element 
   ! ( j QJ :: W^(k_q k_j) :: j Q'J') for subshells with j= 1/2, 3/2, 5/2, 
   ! and 7/2 by looking up appropriate tables.
   !-----------------------------------------------------------------------
      !
      integer, intent(in) :: no_bra, no_ket, k1, k2
      real(kind=dp)       :: coeff
      !       
      coeff = zero
      !   
      if (rcfp_min_even(no_bra)  /=  rcfp_min_even(no_ket))     return
      if (triangle(terms_jj(no_ket)%Q +1, 2*k1 +1,                   &
                   terms_jj(no_bra)%Q +1) == 0)                 return      
      if (triangle(terms_jj(no_bra)%subshellJ +1, 2*k2 +1,           &
                   terms_jj(no_ket)%subshellJ +1) == 0)         return
      !
      if (k1 == 0  .and.  k2 == 0) then
         ! Determine the value of the reduced matrix element of W^(00)
         coeff = rcfp_get_reduced_W_all(W_00_one_half, W_00_three_half, &
	                                W_00_five_half,W_00_seven_half, &
				        no_bra,no_ket)
      else if (k1 == 1  .and.  k2 == 0) then
         ! Determine the value of the reduced matrix element of W^(10)
         coeff = rcfp_get_reduced_W_all(W_10_one_half, W_10_three_half, &
	                                W_10_five_half,W_10_seven_half, &
				        no_bra,no_ket)
      else if (k1 == 0  .and.  k2 == 1) then
         ! Determine the value of the reduced matrix element of W^(01)
         coeff = rcfp_get_reduced_W_all(W_01_one_half, W_01_three_half, &
	                                W_01_five_half,W_01_seven_half, &
				        no_bra,no_ket)
      else if (k1 == 1  .and.  k2 == 2) then
         ! Determine the value of the reduced matrix element of W^(12)
         select case (terms_jj(no_bra)%j)
         case (3)       
            ! j=3/2
            coeff = rcfp_get_reduced_W_three_half(W_12_three_half_odd,  &
	                                          W_12_three_half_even, &
					          no_bra,no_ket)
         case (5)       
            ! j=5/2
            coeff = rcfp_get_reduced_W_five_half (W_12_five_half_odd,   &
	                                          W_12_five_half_even,  &
					          no_bra,no_ket)
         case (7)
            ! j=7/2
            coeff = rcfp_get_reduced_W_seven_half(W_12_seven_half_odd,  &
	                                          W_12_seven_half_even, &
					          no_bra,no_ket)
         end select
      else if (k1 == 0  .and.  k2 == 3) then
         ! Determine the value of the reduced matrix element of W^(03)
         select case (terms_jj(no_bra)%j)
         case (3)
            ! j=3/2
            coeff = rcfp_get_reduced_W_three_half(W_03_three_half_odd,  &
	                                          W_03_three_half_even, &
				                  no_bra,no_ket)
         case (5)
            ! j=5/2
            coeff = rcfp_get_reduced_W_five_half (W_03_five_half_odd,   &
	                                          W_03_five_half_even,  &
					          no_bra,no_ket)
         case (7)
            ! j=7/2
            coeff = rcfp_get_reduced_W_seven_half(W_03_seven_half_odd,  &
	                                          W_03_seven_half_even, &
					          no_bra,no_ket)
         end select
      else if (k1 == 1  .and.  k2 == 4) then
         ! Determine the value of the reduced matrix element of W^(14)
         select case (terms_jj(no_bra)%j)
         case (5)
            ! j=5/2
            coeff = rcfp_get_reduced_W_five_half (W_14_five_half_odd,   &
	                                          W_14_five_half_even,  &
					          no_bra,no_ket)
         case (7)
            ! j=7/2
            coeff = rcfp_get_reduced_W_seven_half(W_14_seven_half_odd,  &
	                                          W_14_seven_half_even, &
					          no_bra,no_ket)
         end select
      else if (k1 == 0  .and.  k2 == 5) then
         ! Determine the value of the reduced matrix element of W^(05)
         select case (terms_jj(no_bra)%j)
         case (5)
            ! j=5/2
            coeff = rcfp_get_reduced_W_five_half (W_05_five_half_odd,   &
     	                                          W_05_five_half_even,  &
					          no_bra,no_ket)
         case (7)
            ! j=7/2
            coeff = rcfp_get_reduced_W_seven_half(W_05_seven_half_odd,  &
	                                          W_05_seven_half_even, &
					          no_bra,no_ket)
         end select
      else if (k1 == 1  .and.  k2 == 6) then	
         ! Determine the value of the reduced matrix element of W(05)
         select case (terms_jj(no_bra)%j)
         case (7)
            ! j=7/2
            coeff = rcfp_get_reduced_W_seven_half(W_16_seven_half_odd,  &
	                                          W_16_seven_half_even, &
					          no_bra,no_ket)
         end select
      else if (k1 == 0  .and.  k2 == 7) then	
         ! Determine the value of the reduced matrix element of W(05)
         select case (terms_jj(no_bra)%j)
         case (7)
            ! j=7/2
            coeff = rcfp_get_reduced_W_seven_half(W_07_seven_half_odd,  &
	                                          W_07_seven_half_even, &
					          no_bra,no_ket)
         end select
      end if
      !
   end function rcfp_get_reduced_W_from_table
   !
   !
   function rcfp_get_reduced_W_all(w1_2,w3_2,w5_2,w7_2,no_bra,no_ket)   &
                                                                   result(coeff)
   !-----------------------------------------------------------------------
   ! Auxiliarity routine for rcfp_get_reduced_W. Returns the value of the
   ! reduced matrix elements 
   ! ( j QJ :: W^(00) :: j Q'J'), ( j QJ :: W^(10) :: j Q'J') and 
   ! ( j QJ :: W^(01) :: j Q'J')  for subshells with j= 1/2, 3/2, 5/2, 7/2
   ! from tables.
   !-----------------------------------------------------------------------
      !
      type(reduced_coeff), dimension(1:2), intent(in)  :: w1_2
      type(reduced_coeff), dimension(1:3), intent(in)  :: w3_2
      type(reduced_coeff), dimension(1:6), intent(in)  :: w5_2
      type(reduced_coeff), dimension(1:14), intent(in) :: w7_2            
      integer, intent(in)                              :: no_bra, no_ket
      real(kind=dp)                                    :: coeff
      !       
      integer :: phase, nom, denom
      !
      coeff = zero
      !   
      if (no_bra /= no_ket) return
      !
      select case (terms_jj(no_bra)%j)
      case (1) 
         ! j=1/2      
         phase = w1_2(no_bra)%phase
	 nom   = w1_2(no_bra)%nom
	 denom = w1_2(no_bra)%denom
      case (3)
         ! j=3/2      
         phase = w3_2(no_bra-2)%phase
	 nom   = w3_2(no_bra-2)%nom
	 denom = w3_2(no_bra-2)%denom
      case (5)
         ! j=5/2	       
         phase = w5_2(no_bra-5)%phase
	 nom   = w5_2(no_bra-5)%nom
	 denom = w5_2(no_bra-5)%denom
      case (7)
         ! j=7/2
         phase = w7_2(no_bra-11)%phase
	 nom   = w7_2(no_bra-11)%nom
	 denom = w7_2(no_bra-11)%denom
      end select
      !
      if (phase /= 0) then
         coeff = phase * sqrt (real(nom,kind=dp)/denom)
      end if
      !
   end function rcfp_get_reduced_W_all
   !
   !
   function rcfp_get_reduced_W_three_half(w3_2_o,w3_2_e,no_bra,no_ket) &
                                                              result(coeff)   
   !-----------------------------------------------------------------------
   ! Auxiliarity routine for rcfp_get_reduced_W. Returns the value of the
   ! reduced matrix elements 
   ! ( j QJ :: W^(12) :: j Q'J') and ( j QJ :: W^(03) :: j Q'J') 
   ! for subshell j= 3/2 from tables.
   !-----------------------------------------------------------------------
      !
      type(reduced_coeff), dimension(1,1), intent(in)  :: w3_2_o
      type(reduced_coeff), dimension(2,2), intent(in)  :: w3_2_e
      integer, intent(in)                              :: no_bra, no_ket
      real(kind=dp)                                    :: coeff
      !
      integer :: phase, nom, denom
      !
      if (no_bra == 3) then
         phase = w3_2_o(no_bra-2, no_ket-2)%phase
	 nom   = w3_2_o(no_bra-2, no_ket-2)%nom
	 denom = w3_2_o(no_bra-2, no_ket-2)%denom
      else if (no_bra > 3 .and. no_bra < 6) then
         phase = w3_2_e(no_bra-3, no_ket-3)%phase
	 nom   = w3_2_e(no_bra-3, no_ket-3)%nom
	 denom = w3_2_e(no_bra-3, no_ket-3)%denom
      else if (rabs_use_stop) then
         stop  "rcfp_get_reduced_W_three_half(): program stop A."
      end if
      !
      if (phase /= 0) then
         coeff = phase * sqrt (real(nom,kind=dp)/denom)
      else
         coeff = zero
      end if
      !
   end function rcfp_get_reduced_W_three_half
   !
   !
   function rcfp_get_reduced_W_five_half(w5_2_o,w5_2_e,no_bra,no_ket) &
                                                              result(coeff)   
   !-----------------------------------------------------------------------
   ! Auxiliarity routine for rcfp_get_reduced_W. Returns the value of the
   ! reduced matrix elements 
   ! ( j QJ :: W^(12) :: j Q'J'),    ( j QJ :: W^(03) :: j Q'J'), 
   ! ( j QJ :: W^(14) :: j Q'J') and ( j QJ :: W^(05) :: j Q'J') 
   ! for subshell j= 5/2 from tables. 
   !-----------------------------------------------------------------------
      !
      type(reduced_coeff), dimension(3,3), intent(in)  :: w5_2_o
      type(reduced_coeff), dimension(3,3), intent(in)  :: w5_2_e
      integer, intent(in)                              :: no_bra, no_ket
      real(kind=dp)                                    :: coeff
      !       
      integer :: phase, nom, denom
      !
      if (no_bra > 5 .and. no_bra < 9) then
         phase = w5_2_o(no_bra-5, no_ket-5)%phase
	 nom   = w5_2_o(no_bra-5, no_ket-5)%nom
	 denom = w5_2_o(no_bra-5, no_ket-5)%denom
      else if (no_bra > 8 .and. no_bra < 12) then
         phase = w5_2_e(no_bra-8, no_ket-8)%phase
	 nom   = w5_2_e(no_bra-8, no_ket-8)%nom
	 denom = w5_2_e(no_bra-8, no_ket-8)%denom
      else if (rabs_use_stop) then
         stop  "rcfp_get_reduced_W_five_half(): program stop A."
      end if
      !
      if (phase /= 0) then
         coeff = phase * sqrt (real(nom,kind=dp)/denom)
      else
         coeff = zero
      end if
      !
   end function rcfp_get_reduced_W_five_half
   !
   !
   function rcfp_get_reduced_W_seven_half(w7_2_o,w7_2_e,no_bra,no_ket) &
                                                              result(coeff)
   !-----------------------------------------------------------------------
   ! Auxiliarity routine for rcfp_get_reduced_W. Returns the value of the
   ! reduced matrix elements 
   ! ( j QJ :: W^(12) :: j Q'J'),    ( j QJ :: W^(03) :: j Q'J'), 
   ! ( j QJ :: W^(14) :: j Q'J'),    ( j QJ :: W^(05) :: j Q'J'),
   ! ( j QJ :: W^(16) :: j Q'J') and ( j QJ :: W^(07) :: j Q'J') 
   ! for subshell j= 7/2 from tables. 
   !-----------------------------------------------------------------------
      !
      type(reduced_coeff), dimension(21), intent(in)  :: w7_2_o
      type(reduced_coeff), dimension(36), intent(in)  :: w7_2_e
      integer, intent(in)                             :: no_bra, no_ket
      real(kind=dp)                                   :: coeff
      !       
      integer :: phase, nom, denom, no_a, no_b, j_parameter, switch            
      integer, dimension(6), parameter ::                          &
         limits_1 = (/ 0, 5, 9, 12, 14, 15  /)
      integer, dimension(8), parameter ::                          &
         limits_2 = (/ 0, 7, 13, 18, 22, 25, 27, 28 /)
      !
      coeff = zero
      if (no_bra < 12 .and. no_bra > 25) return
      if (no_ket < 12 .and. no_ket > 25) return
      !   
      if (no_bra > no_ket) then
         no_a = no_ket;   no_b = no_bra
      else
         no_a = no_bra;   no_b = no_ket
      end if 
      !      
      if (no_bra > 17) then
         no_a        = no_a - 17
	 no_b        = no_b - 17	 
         j_parameter = limits_2(no_a) + no_b
	 switch      = 2
      else
         no_a        = no_a - 11
         no_b        = no_b - 11	 
         j_parameter = limits_1(no_a) + no_b
	 switch      = 1
      end if
      if (switch == 1) then
         phase = w7_2_o(j_parameter)%phase
	 nom   = w7_2_o(j_parameter)%nom
	 denom = w7_2_o(j_parameter)%denom
      else if (switch == 2 ) then 
         phase = w7_2_e(j_parameter)%phase
	 nom   = w7_2_e(j_parameter)%nom
	 denom = w7_2_e(j_parameter)%denom
      else if (rabs_use_stop) then 
         stop  "rcfp_get_reduced_W_seven_half(): program stop A."
      end if
      !
      if (phase /= 0) then
         coeff = phase * sqrt (real(nom,kind=dp)/denom)
         if (no_bra > no_ket .and.                                         &
             mod(terms_jj(no_ket)%subshellJ - terms_jj(no_bra)%subshellJ + &
	         terms_jj(no_ket)%Q - terms_jj(no_bra)%Q, 4) /= 0) then
	    coeff = - coeff
         end if
      else
         coeff = zero
      end if
      !
   end function rcfp_get_reduced_W_seven_half
   !
   !
   function rcfp_calculate_Wk_me(bra,ket,k_j,q_m1,q_m2)       result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the matrix element
   ! ( j^N QJ :: W^(k_j) :: j^N' Q'J') for subshells with j = 1/2, 3/2, 5/2,
   ! 7/2, 9/2
   ! and for j > 9/2 with N = 0, 1, 2. 
   !
   ! Calls: rcfp_Clebsch_Gordan_qusispin().
   !-----------------------------------------------------------------------
      !
      type (subshell_state), intent(in) :: bra, ket
      integer, intent(in)               :: k_j, q_m1, q_m2
      real(kind=dp)                     :: coeff
      !
      type (subshell_state) :: run
      integer               :: k_q, run_nu, run_subshellJ, min_run, max_run
      integer               :: run_i, delta_J, j, bra_subshellJ, ket_subshellJ
      !
      coeff = zero
      !
      if (bra%state < 64  .and.  ket%state < 64) then
         if (rcfp_min_even(bra%state)  /=  rcfp_min_even(ket%state)  .or.   &
             triangle(terms_jj(bra%state)%subshellJ +1,2 * k_j +1,          &
                      terms_jj(ket%state)%subshellJ + 1) == 0        .or.   &
             rcfp_Q_space_delta(bra) == 0                            .or.   &
             rcfp_Q_space_delta(ket) == 0                            .or.   &
             bra%nq - ket%nq - q_m1 - q_m2  /= 0) return
         !
         if (q_m1 == q_m2 ) then
            !
            ! cases    a * a    and    a^+ * a^+
            if (mod(k_j,2) /= 0) return
	    if (triangle(terms_jj(ket%state)%Q + 1,3,                       &
	                 terms_jj(bra%state)%Q + 1) == 0) return
            coeff = rcfp_Clebsch_Gordan_qusispin                            &
	                           (terms_jj(ket%state)%Q,ket%subshellMQ,   &
                                   2, q_m1 + q_m2,                          &
                                   terms_jj(bra%state)%Q,bra%subshellMQ)
            if (dabs(coeff) < eps10) return
            coeff = coeff * rcfp_get_reduced_W(bra%state,ket%state,1,k_j)
            if (dabs(coeff) < eps10) return
	    coeff = coeff / sqrt(terms_jj(bra%state)%Q + one)
         else
            !
            ! cases    a * a^+    and    a^+ * a
	    if (k_j  == 0) then
               if (bra%state /= ket%state) return
               if (q_m1  == 1) then
	          coeff = - bra%nq
	       else
	          coeff = terms_jj(bra%state)%j + one - bra%nq
	       end if
               coeff = coeff * sqrt((terms_jj(bra%state)%subshellJ + one)   &
	                          / (terms_jj(bra%state)%j + one))
	    else
               if (mod(k_j,2) == 0) then
	          k_q = 1
	       else
	          k_q = 0
	       end if	    
	       if (triangle(terms_jj(ket%state)%Q +1,2 * k_q+1,             &
	                    terms_jj(bra%state)%Q + 1) == 0) return
               coeff = rcfp_Clebsch_Gordan_qusispin                         &
	                              (terms_jj(ket%state)%Q,ket%subshellMQ,&
	                              2 * k_q, q_m1 + q_m2,                 &
                                      terms_jj(bra%state)%Q,bra%subshellMQ)
               if (dabs(coeff) < eps10) return
               coeff = coeff * rcfp_get_reduced_W(bra%state,ket%state,k_q,k_j)
               if (dabs(coeff) < eps10) return
	       coeff = coeff / sqrt(two * (terms_jj(bra%state)%Q + one))
               if (q_m1 == -1  .and.  mod(k_j,2) /= 0) then
	          coeff = - coeff
	       end if
	    end if
         end if
      else if (bra%state > 64  .and.  ket%state > 64) then
         j             =mod(bra%state/(1000*1000),1000)
         bra_subshellJ =mod(bra%state,1000);  ket_subshellJ =mod(ket%state,1000)
         if (rabs_use_stop) then
         if (ket%nq > 2) then
            stop  "rcfp_calculate_Wk_me(): program stop A."
	 end if
         if (bra%nq > 2) then
            stop  "rcfp_calculate_Wk_me(): program stop B."
	 end if
         if (q_m1 + q_m2 + ket%nq /= bra%nq) then
            stop  "rcfp_calculate_Wk_me(): program stop C."
	 end if
	 end if
         run%nq = ket%nq + q_m2
         select case (run%nq)
         case (0)
            min_run = 0;   max_run = 0
         case (1)
            min_run = j;   max_run = j
         case (2)
            min_run = 0;   max_run = 2*j - 2
         case default
            stop  "rcfp_calculate_Wk_me(): program stop D."
         end select
         do run_i = min_run, max_run, 4
            run_subshellJ = run_i
	    delta_J = wigner_6j_triangle(j,j,2*k_j,ket_subshellJ,bra_subshellJ,&
	                                                          run_subshellJ)
            if (delta_J /= 0) then
               select case (run%nq)
               case (0)
                  run_nu = 0
               case (1)
                  run_nu = 1
               case (2)
                  if (run_subshellJ == 0 ) then
                     run_nu = 0
                  else
                     run_nu = 2
                  end if
               case default
                  stop  "rcfp_calculate_Wk_me(): program stop E."
               end select
               run%subshellMQ = run%nq - (j + 1)/2
               run%state = ((j * 1000) + run_nu) * 1000 + run_subshellJ
               coeff = coeff + rcfp_calculate_a_me(bra,run,q_m1)*              &
                       rcfp_calculate_a_me(run,ket,q_m2)*                      &
	               wigner_6j_symbol(j,j,2*k_j,ket_subshellJ,bra_subshellJ, &
	                                                   run_subshellJ,.true.)
            end if
         end do
         coeff = coeff * sqrt(two*k_j + one)
         if(mod(bra_subshellJ + ket_subshellJ + 2 * k_j,4) /= 0) coeff = -coeff
      else if (rabs_use_stop) then
         stop  "rcfp_calculate_Wk_me(): program stop F."
      end if
      !
   end function rcfp_calculate_Wk_me
   !
   !
   function rcfp_calculate_Wk_times_Wk_0_me(bra,ket,k_j,                  &
                                     q_m1,q_m2,q_m3,q_m4)     result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the matrix element 
   ! ( j^N QJ :: [ W^(k_j) * W^(k_j)]^(0) :: j^N' QJ) 
   ! for subshells with j = 1/2, 3/2, 5/2, 7/2, 9/2 
   ! and for j > 9/2 with N = 0, 1, 2. 
   !
   ! Calls: rcfp_calculate_Wk_me().
   !-----------------------------------------------------------------------
      !
      type (subshell_state), intent(in) :: bra, ket
      integer, intent(in)               :: k_j, q_m1, q_m2, q_m3, q_m4
      real(kind=dp)                     :: coeff
      !
      integer               :: no_run, delta_J
      type (subshell_state) :: run
      real(kind=dp)         :: coeff1
      integer               :: run_nu, run_subshellJ, min_run, max_run
      integer               :: run_i, j, bra_subshellJ, ket_subshellJ
      !   
      coeff = zero
      !   
      if (bra%state < 64  .and.  ket%state < 64) then
         if(rcfp_Q_space_delta(bra)      ==                             0)return
         if(rcfp_Q_space_delta(ket)      ==                             0)return
         if(terms_jj(bra%state)%subshellJ/= terms_jj(ket%state)%subshellJ)return
         if(rcfp_min_even(bra%state)     /= rcfp_min_even(ket%state)     )return
         if(rcfp_max_even(bra%state)     /= rcfp_max_even(ket%state)     )return
         if(bra%nq - ket%nq - q_m1 - q_m2 - q_m3 - q_m4              /= 0)return
         !
         do no_run = rcfp_min_even(bra%state), rcfp_max_even(bra%state)
            run = subshell_state(no_run, bra%n, ket%nq +q_m3 +q_m4,     &
                           ket%nq +q_m3 +q_m4 -(terms_jj(no_run)%j +1)/2) 
            if (terms_jj(run%state)  %Q >=  abs(run%subshellMQ))then
	       delta_J = wigner_6j_triangle(2*k_j,2*k_j,0,        &
	                           terms_jj(ket%state)%subshellJ, &
		                   terms_jj(bra%state)%subshellJ, & 
	                           terms_jj(run%state)%subshellJ)
               if (delta_J /= 0) then
                  coeff1 = rcfp_calculate_Wk_me(bra,run,k_j,q_m1,q_m2) * &
	                   rcfp_calculate_Wk_me(run,ket,k_j,q_m3,q_m4)
	          if( mod(2*k_j - terms_jj(bra%state)%subshellJ +   &
                                  terms_jj(run%state)%subshellJ, 4) &
			          /= 0) coeff1 = - coeff1
                  coeff = coeff + coeff1		   
	       end if	 
	    end if
         end do
         coeff = coeff/sqrt((two*k_j +one)*(terms_jj(bra%state)%subshellJ +one))
      else if (bra%state > 64  .and.  ket%state > 64) then
         if (bra%state /= ket%state) return
         j             =mod(bra%state/(1000*1000),1000)
         bra_subshellJ =mod(bra%state,1000);  ket_subshellJ =mod(ket%state,1000)
         if (rabs_use_stop) then
         if (ket%nq > 2) then
            stop  "rcfp_calculate_Wk_times_Wk_0_me(): program stop A."
	 end if
         if (bra%nq > 2) then
            stop  "rcfp_calculate_Wk_times_Wk_0_me(): program stop B."
	 end if
         if (q_m1 + q_m2 + q_m3 + q_m4 + ket%nq /= bra%nq) then
            stop  "rcfp_calculate_Wk_times_Wk_0_me(): program stop C."
	 end if
	 end if
         !
         run%nq = ket%nq + q_m3 + q_m4
         select case (run%nq)
         case (0)
            min_run = 0;   max_run = 0
         case (1)
            min_run = j;   max_run = j
         case (2)
            min_run = 0;   max_run = 2*j - 2
         case default
            stop  "rcfp_calculate_Wk_times_Wk_0_me(): program stop D."
         end select
         do run_i = min_run, max_run, 4
            run_subshellJ = run_i
	    delta_J = wigner_6j_triangle(2*k_j,2*k_j,0,ket_subshellJ,          &
	                                            bra_subshellJ,run_subshellJ)
            if (delta_J /= 0) then
               select case (run%nq)
               case (0)
                  run_nu = 0
               case (1)
                  run_nu = 1
               case (2)
                  if (run_subshellJ == 0 ) then
                     run_nu = 0
                  else
                     run_nu = 2
                  end if
               case default
                  stop  "rcfp_calculate_Wk_times_Wk_0_me(): program stop E."
               end select
               run%subshellMQ = run%nq - (j + 1)/2
	       run%state = ((j * 1000) + run_nu) * 1000 + run_subshellJ
               coeff = coeff + rcfp_calculate_Wk_me(bra,run,k_j,q_m1,q_m2)*    &
                       rcfp_calculate_Wk_me(run,ket,k_j,q_m3,q_m4)*            &
	               wigner_6j_symbol(2*k_j,2*k_j,0,ket_subshellJ,           &
		                             bra_subshellJ,run_subshellJ,.true.)
            end if
         end do
         if(mod(bra_subshellJ + ket_subshellJ,4) /= 0) coeff = -coeff
      else if(rabs_use_stop) then
         stop  "rcfp_calculate_Wk_times_Wk_0_me(): program stop F."
      end if
      !       
   end function rcfp_calculate_Wk_times_Wk_0_me
   !
   !
   function rcfp_calculate_a_times_Wk_me(bra,ket,k_j1,kk_j2,              &
                                         q_m1,q_m2,q_m3)      result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the matrix element
   ! ( j^N QJ :: [ a^(q j)_m_q * W^(k_j1)]^(k_j2) :: j^N' QJ)
   ! for subshells with j= 1/2, 3/2, 5/2, 7/2, 9/2.
   !
   ! Calls: rcfp_Clebsch_Gordan_qusispin(), rcfp_calculate_Wk_me(), 
   !        wigner_6j_symbol().
   !-----------------------------------------------------------------------
      !
      type (subshell_state), intent(in) :: bra, ket
      integer, intent(in)               :: k_j1, kk_j2, q_m1, q_m2, q_m3
      real(kind=dp)                     :: coeff
      !
      type (subshell_state) :: run
      integer               :: no_run, delta_J
      real(kind=dp)         :: coeff1
      integer               :: run_nu, run_subshellJ, min_run, max_run
      integer               :: run_i, j, bra_subshellJ, ket_subshellJ
      !
      coeff = zero
      !
      if (bra%state < 64  .and.  ket%state < 64) then
         if (triangle(terms_jj(bra%state)%subshellJ +1, kk_j2+1,        &
                      terms_jj(ket%state)%subshellJ +1)      == 0) return
         if (rcfp_Q_space_delta(bra) ==                         0) return
         if (rcfp_Q_space_delta(ket) ==                         0) return
         if (rcfp_min_odd(bra%state) /=  rcfp_min_even(ket%state)) return
         if (rcfp_max_odd(bra%state) /=  rcfp_max_even(ket%state)) return
         if (bra%nq - ket%nq - q_m1 - q_m2 - q_m3            /= 0) return
         !
         do no_run = rcfp_min_odd(bra%state), rcfp_max_odd(bra%state)
            run = subshell_state(no_run, bra%n, ket%nq +q_m2 +q_m3, &
                       ket%nq +q_m2 +q_m3 -(terms_jj(no_run)%j +1)/2)
            if (terms_jj(run%state)%Q  >=  abs(run%subshellMQ)) then
	       delta_J = wigner_6j_triangle(terms_jj(bra%state)%j,2*k_j1,kk_j2,&
		                   terms_jj(ket%state)%subshellJ,              &
		                   terms_jj(bra%state)%subshellJ,              &
	                           terms_jj(run%state)%subshellJ)
               if (delta_J /= 0) then
                  coeff1 = rcfp_Clebsch_Gordan_qusispin                  &
		                   (terms_jj(run%state)%Q,run%subshellMQ,&
                                   1, q_m1,terms_jj(bra%state)%Q,bra%subshellMQ)
                  if (dabs(coeff1) > eps10) then
                     coeff1 = coeff1*rcfp_get_coefficient(bra%state,run%state)*&
                              rcfp_calculate_Wk_me(run,ket,k_j1,q_m2,q_m3)    *&
	                      wigner_6j_symbol(terms_jj(bra%state)%j,          &
			      2*k_j1,kk_j2,terms_jj(ket%state)%subshellJ,      &
		                           terms_jj(bra%state)%subshellJ,      &
	                                   terms_jj(run%state)%subshellJ,.true.)
                     coeff = coeff + coeff1 / sqrt(terms_jj(bra%state)%Q + one)
	          end if
	       end if
	    end if
         end do
         coeff = coeff * sqrt(kk_j2 + one)
         if (mod(kk_j2 + terms_jj(bra%state)%subshellJ +                       &
                 terms_jj(ket%state)%subshellJ + 2, 4) /= 0) coeff = - coeff
      else if (bra%state > 64  .and.  ket%state > 64) then
         j             =mod(bra%state/(1000*1000),1000)
         bra_subshellJ =mod(bra%state,1000);  ket_subshellJ =mod(ket%state,1000)
         if (rabs_use_stop) then
         if (ket%nq > 2) then
            stop  "rcfp_calculate_a_times_Wk_me(): program stop A."
	 end if
         if (bra%nq > 2) then
            stop  "rcfp_calculate_a_times_Wk_me(): program stop B."
	 end if
         if (q_m1 + q_m2 + q_m3  + ket%nq /= bra%nq) then
            stop  "rcfp_calculate_a_times_Wk_me(): program stop C."
	 end if
	 end if
         run%nq = ket%nq + q_m2 + q_m3
         select case (run%nq)
         case (0)
            min_run = 0;   max_run = 0
         case (1)
            min_run = j;   max_run = j
         case (2)
            min_run = 0;   max_run = 2*j - 2
         case default
            stop  "rcfp_calculate_a_times_Wk_me(): program stop D."
         end select
         do run_i = min_run, max_run, 4
            run_subshellJ = run_i
	    delta_J = wigner_6j_triangle(j,2*k_j1,kk_j2,ket_subshellJ,         &
	                                            bra_subshellJ,run_subshellJ)
            if (delta_J /= 0) then
               select case (run%nq)
               case (0)
                  run_nu = 0
               case (1)
                  run_nu = 1
               case (2)
                  if (run_subshellJ == 0 ) then
                     run_nu = 0
                  else
                     run_nu = 2
                  end if
               case default
                  stop  "rcfp_calculate_a_times_Wk_me(): program stop E."
               end select
               run%subshellMQ = run%nq - (j + 1)/2
	       run%state = ((j * 1000) + run_nu) * 1000 + run_subshellJ
               coeff = coeff + rcfp_calculate_a_me(bra,run,q_m1)*              &
                       rcfp_calculate_Wk_me(run,ket,k_j1,q_m2,q_m3)*           &
	               wigner_6j_symbol(j,2*k_j1,kk_j2,ket_subshellJ,          &
		                             bra_subshellJ,run_subshellJ,.true.)
            end if
         end do
         coeff = coeff * sqrt(kk_j2 + one)
         if(mod(bra_subshellJ + kk_j2 + ket_subshellJ,4) /= 0) coeff = -coeff
      else if (rabs_use_stop) then
         stop  "rcfp_calculate_a_times_Wk_me(): program stop F."
      end if
      !
   end function rcfp_calculate_a_times_Wk_me
   !
   !
   function rcfp_calculate_Wk_times_a_me(bra,ket,k_j1,kk_j2,q_m1,q_m2,q_m3) &
                                                                   result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the matrix element
   ! ( j^N QJ :: [ a^(q j)_m_q * W^(k_j1)]^(k_j2) :: j^N' QJ) 
   ! for subshells with j= 1/2, 3/2, 5/2, 7/2, 9/2. 
   !
   ! Calls: rcfp_Clebsch_Gordan_qusispin(), rcfp_get_coefficient(), 
   ! wigner_6j_symbol(), wigner_6j_triangle().
   !-----------------------------------------------------------------------
      !
      type (subshell_state), intent(in) :: bra, ket
      integer, intent(in)               :: k_j1, kk_j2, q_m1, q_m2, q_m3
      real(kind=dp)                     :: coeff
      !
      integer               :: no_run, delta_J
      type (subshell_state) :: run
      real(kind=dp)         :: coeff1
      integer               :: run_nu, run_subshellJ, min_run, max_run
      integer               :: run_i, j, bra_subshellJ, ket_subshellJ
      !   
      coeff = zero
      !   
      if (bra%state < 64  .and.  ket%state < 64) then
         if (triangle(terms_jj(bra%state)%subshellJ +1, kk_j2+1,       &
                   terms_jj(ket%state)%subshellJ +1)        == 0) return
         if (rcfp_Q_space_delta(bra)  ==                       0) return
         if (rcfp_Q_space_delta(ket)  ==                       0) return
         if (rcfp_min_even(bra%state) /= rcfp_min_odd(ket%state)) return
         if (rcfp_max_even(bra%state) /= rcfp_max_odd(ket%state)) return
         if (bra%nq -ket%nq - q_m1 - q_m2 - q_m3            /= 0) return
         !
         do no_run = rcfp_min_even(bra%state), rcfp_max_even(bra%state)
            run = subshell_state(no_run, bra%n, ket%nq + q_m3,               &
                                 ket%nq + q_m3 -(terms_jj(no_run)%j +1)/2)
             if (terms_jj(run%state)%Q >= abs(run%subshellMQ))then
	        delta_J = wigner_6j_triangle(2*k_j1,terms_jj(bra%state)%j,     &
	                                   kk_j2,terms_jj(ket%state)%subshellJ,&
		                           terms_jj(bra%state)%subshellJ,      &
	                                   terms_jj(run%state)%subshellJ)
               if (delta_J /= 0) then
                  coeff1 = rcfp_Clebsch_Gordan_qusispin                        &
		                         (terms_jj(ket%state)%Q,ket%subshellMQ,&
                          	          1, q_m3,                             &
                                          terms_jj(run%state)%Q,run%subshellMQ)	
                  if (dabs(coeff1) > eps10 ) then
                     coeff1 = coeff1*rcfp_get_coefficient(run%state,ket%state)*&
                              rcfp_calculate_Wk_me(bra,run,k_j1,q_m1,q_m2)    *&
	                      wigner_6j_symbol(2*k_j1,terms_jj(bra%state)%j,   &
	                                kk_j2,terms_jj(ket%state)%subshellJ,   &
		                              terms_jj(bra%state)%subshellJ,   &
	                                terms_jj(run%state)%subshellJ,.true.)
                     coeff = coeff + coeff1 / sqrt(terms_jj(run%state)%Q + one)
	          end if
	       end if
	    end if
         end do
         coeff = coeff * sqrt(kk_j2 + one)
         if (mod(kk_j2 + terms_jj(bra%state)%subshellJ +                    &
                 terms_jj(ket%state)%subshellJ + 2, 4)  /=  0) coeff = - coeff
      else if (bra%state > 64  .and.  ket%state > 64) then
         j             =mod(bra%state/(1000*1000),1000)
         bra_subshellJ =mod(bra%state,1000);  ket_subshellJ =mod(ket%state,1000)
         if (rabs_use_stop) then
         if (ket%nq > 2) then
            stop  "rcfp_calculate_Wk_times_a_me(): program stop A."
	 end if
         if (bra%nq > 2) then
            stop  "rcfp_calculate_Wk_times_a_me(): program stop B."
	 end if
         if (q_m1 + q_m2 + q_m3  + ket%nq /= bra%nq) then
            stop  "rcfp_calculate_Wk_times_a_me(): program stop C."
	 end if
	 end if
         !
         run%nq = ket%nq + q_m3
         select case (run%nq)
         case (0)
            min_run = 0;   max_run = 0
         case (1)
            min_run = j;   max_run = j
         case (2)
            min_run = 0;   max_run = 2*j - 2
         case default
            stop  "rcfp_calculate_Wk_times_a_me(): program stop D."
         end select
         do run_i = min_run, max_run, 4
            run_subshellJ = run_i
	    delta_J = wigner_6j_triangle(2*k_j1,j,kk_j2,ket_subshellJ,         &
	                                            bra_subshellJ,run_subshellJ)
            if (delta_J /= 0) then
               select case (run%nq)
               case (0)
                  run_nu = 0
               case (1)
                  run_nu = 1
               case (2)
                  if (run_subshellJ == 0 ) then
                     run_nu = 0
                  else
                     run_nu = 2
                  end if
               case default
                  stop  "rcfp_calculate_Wk_times_a_me(): program stop E."
               end select
               run%subshellMQ = run%nq - (j + 1)/2
	       run%state = ((j * 1000) + run_nu) * 1000 + run_subshellJ
               coeff = coeff + rcfp_calculate_Wk_me(bra,run,k_j1,q_m1,q_m2)*   &
                       rcfp_calculate_a_me(run,ket,q_m3)*                      &
	               wigner_6j_symbol(2*k_j1,j,kk_j2,ket_subshellJ,          &
		                             bra_subshellJ,run_subshellJ,.true.)
            end if
         end do
         coeff = coeff * sqrt(kk_j2 + one)
         if(mod(bra_subshellJ + kk_j2 + ket_subshellJ,4) /= 0) coeff = -coeff
      else if (rabs_use_stop) then
         stop  "rcfp_calculate_Wk_times_a_0_me(): program stop F."
      end if
      !
    end function rcfp_calculate_Wk_times_a_me
   !
   !
   function rcfp_calculate_a_me(bra,ket,q_m)                  result(coeff)
   !-----------------------------------------------------------------------
   ! Returns the value of the matrix elements
   ! (j^N QJ :: a^(q j)_m_q :: j^N' QJ)
   ! for subshells with j= 1/2, 3/2, 5/2, 7/2, 9/2
   ! and for j > 9/2 with N = 0, 1, 2.
   !-----------------------------------------------------------------------
      !
      type (subshell_state), intent(in) :: bra, ket
      integer, intent(in)               :: q_m
      real(kind=dp)                     :: coeff
      !   
      integer   :: j, bra_nu, ket_nu, bra_subshellJ, ket_subshellJ
      !
      coeff = zero
      !   
      if (bra%state < 64  .and.  ket%state < 64) then
         if (rcfp_min_odd(bra%state) /=  rcfp_min_even(ket%state)       ) return
         if (bra%nq - ket%nq - q_m                                  /= 0) return
         if (rcfp_Q_space_delta(bra)                                == 0) return
         if (rcfp_Q_space_delta(ket)                                == 0) return
         if (triangle(terms_jj(ket%state)%Q +1,2,                              &
                      terms_jj(bra%state)%Q +1)                     == 0) return
         if (triangle(terms_jj(bra%state)%subshellJ +1,terms_jj(bra%state)%j+1,&
                       terms_jj(ket%state)%subshellJ +1)            == 0) return
         ! 
         coeff = - rcfp_Clebsch_Gordan_qusispin                                &
	                          (terms_jj(ket%state)%Q, ket%subshellMQ,1,q_m,&
                                  terms_jj(bra%state)%Q, bra%subshellMQ)   *   &
                 rcfp_get_coefficient(bra%state,ket%state)                 /   &
                 sqrt(terms_jj(bra%state)%Q+one)
      else if (bra%state > 64  .and.  ket%state > 64) then
         bra_subshellJ = mod(bra%state,1000)
         if (q_m == -1) then
            j             = mod(bra%state/(1000*1000),1000)
            bra_nu        = mod(bra%state/1000,1000)
            ket_nu        = mod(ket%state/1000,1000)
            ket_subshellJ = mod(ket%state,1000)
            coeff = sqrt(ket%nq * (ket_subshellJ + one))
            if (mod(bra_subshellJ-ket_subshellJ-j+2*ket%nq,4)  &
	                                                  /= 0) coeff = - coeff
         else if (q_m == 1) then
            coeff = sqrt(bra%nq * (bra_subshellJ + one))
            if (mod(bra%nq,2) /= 0) coeff = - coeff
         end if
      else if (rabs_use_stop) then
         stop  "rcfp_calculate_a_me(): program stop A."
      end if
      !       
   end function rcfp_calculate_a_me
   !
   !
   function rcfp_Q_space_delta(sub_state)                     result(Delta)   
   !-----------------------------------------------------------------------
   ! Calculate trivial delta factors for Q space.
   !-----------------------------------------------------------------------
      !
      type (subshell_state), intent(in) :: sub_state   
      integer                           :: Delta
      !
      Delta = 0
      if (    terms_jj(sub_state%state)%Q < abs(sub_state%subshellMQ)   ) return
      if (mod(terms_jj(sub_state%state)%Q + sub_state%subshellMQ,2) /= 0) return
      Delta = 1
      !
   end function rcfp_Q_space_delta
   !
   !
   subroutine rcfp_input
   !-----------------------------------------------------------------------
   ! Reads in the input data and carries out the calculations of various
   ! (reduced) coefficients and reduced matrix elements. Calculations can
   ! currently be made for: 
   !    1) coefficients of fractional parentage, 
   !    2) reduced coefficient of fractional parentage, 
   !    3) reduced matrix element of W^{(k_q k_j)} operators, and 
   !    4) matrix element of the T^{(k)} unit operators 
   ! in jj coupling for subshells with j= 1/2, 3/2, 5/2, 7/2, and 9/2 . 
   ! The date input and output is done interactibely.
   !
   ! Calls:
   !-----------------------------------------------------------------------
      !
      integer      :: i, selection, rank_q, rank_j, rank_k
      integer      :: length_bra, length_ket, length_add, length_operator 
      integer      :: bra_nq, ket_nq, bra_j, ket_j, bra_subshellJ, ket_subshellJ
      integer      :: bra_nu, ket_nu, bra_nr, ket_nr, bra_MQ, ket_MQ
      integer      :: bra_no, ket_no
      real(kind=dp):: coeff
      logical      :: fail
      character(len=180) :: answer
      !
      ! Restart if another calculation is to be carried out
    1 continue
      !
      print *, "Select one item from the list for calculating:"
      print *, " "
      print *, "  1:  coefficients of fractional parentage,"
      print *, "  2:  reduced coefficients of fractional parentage,"
      print *, "  3:  completely reduced matrix elements of the operator" //&
                       " W^{k_q k_j},"
      print *, "  4:  reduced matrix elemenst of unit operator T^{(k)},"
      print *, "  q:  to quit the program."
      print *, " "
      read (*,"(a)") answer
      answer = adjustl(answer)
      if (answer(1:1) == "q"  .or. answer(1:1) == "Q")   return
      if (answer(1:1) /= "1"  .and.  answer(1:1) /= "2"  .and.  &
          answer(1:1) /= "3"  .and.  answer(1:1) /= "4") goto 1
      !
      selection = iachar(answer(1:1)) - iachar("1") + 1
      if (selection < 1 .or. selection > 4 .or. len_trim(answer) > 1) then
         print*, " error in input  !!! "
	 print*, " reenter "
	 go to 1
      end if	  
      !
    2 continue
      select case (answer(1:1))
      case ("1")
         print *, "Calculate a cfp coefficient" //&
                  " (j^N  nu J {| j^{N-1} nu' J',  j) :" 
         selection = 1     
      case ("2")
         print *, "Calculate a reduced cfp coefficient" //&
                  " (j   nu J ||| a^{(1/2 j)} ||| j  nu' J') :"      
         selection = 2     
      case ("3")
         print *, "Calculate a completely reduced matrix element" //&
                  " (j   nu J ||| W^{k_q k_j} ||| j  nu' J') :"
         selection = 3     
      case ("4")
         print *, "Calculate a reduced matrix element of the unit operator" //&
                  " (j^N   nu J || T^{k} || j^N  nu' J')"
         selection = 4     
      end select
      print *, " "
      !
      ! Now analyze the bra function
      read (*,"(a)",err=2) answer
      if (len_trim(answer) < 1) then
         goto 2
      else if (answer(1:1) == "b" .or. answer(1:1) == "B") then
         goto 1
      else if (answer(1:1) == "q" .or. answer(1:1) == "Q") then
         return
      end if
      !
      length_bra = len_trim(answer)
      call rcfp_analyse(answer(1:length_bra),bra_j,bra_nq,bra_nu,  &
                        bra_subshellJ,fail)
      if (fail) goto 2
      if (selection == 2 .or. selection == 3) then
	 if ( bra_nq /= -100) then
            print*, "No occupation number may be given for calculating"   //&
                    " reduced coefficients and matrix elements; reenter ..."
            goto 2
	 end if
      else if (selection == 1 .or. selection == 4) then
	 if (bra_nq == -100) then
	    print*, "A valid occupation number (^N) must be given for"  //&
                    " calculating cfp coefficients; reenter ..."
            goto 2
	  end if
      end if
      !
      ! Now distinguish degenerate subshell states for j=9/2 in the bra 
      ! function; if needed, an additional state identifier Nr. will be required
      if ( (bra_j == 9 .and. bra_nu == 4 .and. bra_subshellJ == 8)   .or.   &
	   (bra_j == 9 .and. bra_nu == 4 .and. bra_subshellJ == 12)) then
         !
       4 print*, "Enter an additional state identifier Nr = 1 or 2."
         answer(length_bra+1:length_bra+5)= " Nr= "	     
         write(unit=*,fmt="(a)",advance="no") answer(1:length_bra+5) 
         read (unit=*,fmt="(a)",err=4)        answer(length_bra+6:180) 
         length_add = len_trim(answer(length_bra+6:180))
         if (length_add == 0) goto 4
         bra_nr = get_integer_from_string(answer(                           &
	               length_bra+6:length_bra+5+length_add),fail) 
         if (fail) goto 4
         if (bra_nr > 2  .or.  bra_nr < 1) goto 4
         length_bra = length_bra + 5 + length_add
      else 
	 bra_nr = 0
      end if
      !
      ! Determine the internal index of the given subshell state.
      bra_no = rcfp_get_term_number                                            &
                                 (bra_j,bra_nq,bra_nu,bra_nr,bra_subshellJ,fail)
      if (fail) then
         print *, "Unable to recognize the subshell state with (j,nu,J) = (", &
                  bra_j,",",bra_nu,",",bra_subshellJ,"); reenter ..."
	 goto 2
      end if
      !
      ! Determine and set-up the operator
      select case (selection)
      case (1)
	 answer(length_bra+1:length_bra+4)  = " {| "
	 length_operator = 4
      case (2)
	 answer(length_bra+1:length_bra+21) = " ||| a^{(1/2 j)} ||| "
	 length_operator = 21
      case (3)
	 answer(length_bra+1:length_bra+9) = " ||| W^{ "
       5 write(unit=*,fmt="(a)",advance="no") answer(1:length_bra+9)
         read (unit=*,fmt="(a)",err=5)        answer(length_bra+10:180)
	 length_add = len_trim(answer(length_bra+10:180))
	 i = length_bra + 10
         !
       6 if (i > length_add+length_bra+9) then
	    print *, "Enter the rank k_q = 0 or 1."
            go to 5
	 end if
	 if (answer(i:i) == " ") then
	    i = i + 1
	    goto 6
	 end if
	 rank_q = get_integer_from_string(answer(i:i),fail)
	 if (fail  .or.  (rank_q /= 0  .and.  rank_q /= 1)) then
            print*, "The rank must be k_q = 0 or 1."
	    goto 5
         end if
	 i = i + 1
         !
       7 if (i > length_add+length_bra+10) then
	    print*, "Enter the ranks k_q = 0 or 1 and k_j. "
	    goto 5
	 end if
	 if (answer(i:i) == " ") then
	    i = i + 1
	    goto 7
	 end if	      
 	 rank_j = get_integer_from_string(answer(i:i),fail)
         if (fail .or. rank_j > bra_j) then
            print*, "The rank k_j must be must integer and k_j <= 2*j."
	    goto 5
         end if
	 answer(length_bra+11+length_add:length_bra+17+length_add) = " } ||| "
	 length_operator = 17 + length_add	  
      case (4)
	 answer(length_bra+1:length_bra+8)= " || T^{ "	  
         !	     
       8 write(unit=*,fmt="(a)",advance="no") answer(1:length_bra+8)
         read (unit=*,fmt="(a)",err=8)        answer(length_bra+9:180)
	 length_add = len_trim(answer(length_bra+9:180))
	 i = length_bra + 9
         !	     	    
       9 if (answer(i:i) == " ") then
	    i = i + 1
	    goto 9
	 end if
         !
	 rank_k = get_integer_from_string(answer(i:i),fail)
         if (fail .or. rank_k > bra_j) then
            print*, "The rank k must be must integer and k <= 2*j."
	    goto 8
         end if
	 answer(length_bra+9+length_add:length_bra+15+length_add) = " } || "
	 length_operator = 15 + length_add	  
      end select
      !
      answer(length_bra+length_operator+1:length_bra+length_operator+1) =     &
	                                                        achar(48+bra_j)
      answer(length_bra+length_operator+2:length_bra+length_operator+4) = "/2 "
      !
      if (selection == 1) then
         ket_nq = bra_nq -1
	 answer(length_bra+length_operator+5:length_bra+length_operator+5)= "^"
         answer(length_bra+length_operator+6:length_bra+length_operator+6)=   &
							       achar(48+ket_nq)
	 answer(length_bra+length_operator+7:length_bra+length_operator+7)= " "
	 length_operator = length_operator + 7
      else if (selection == 4) then
         ket_nq = bra_nq
	 answer(length_bra+length_operator+5:length_bra+length_operator+5)= "^"
         answer(length_bra+length_operator+6:length_bra+length_operator+6)=   &
							       achar(48+ket_nq)
	 answer(length_bra+length_operator+7:length_bra+length_operator+7)= " "
	 length_operator = length_operator + 7      
      else
         length_operator = length_operator + 4	 
      end if
      !
      ! Now analyze the bra function
   10 write(unit=*,fmt="(a)",advance="no") answer(1:length_bra+length_operator) 
      read (unit=*,fmt="(a)") answer(length_bra+length_operator+1:180)
      !
      if (selection == 1) then
         length_ket = len_trim(answer(length_bra+length_operator+1:180))
         if (1 > length_ket) then
	    print *, "Enter the quantum numbers  nu and J."
	    goto 10
         end if
	 length_add = index(answer(length_bra+length_operator+1:180),",")
         if (length_add == 0) then
	    print*, "There is missing the string ' ,' "
            goto 10
         end if
         answer(length_bra+length_operator+length_add:        &
	        length_bra+length_operator+length_add) = ")"
      else
         length_ket = len_trim(answer(length_bra+length_operator+1:180))
         if (1 > length_ket) then
	    print *, "Enter the quantum numbers  nu and J annd terminate" //&
                     " with )." 
	    goto 10
         end if
      end if
      call rcfp_analyse(answer(length_bra+length_operator+1:                &
	                       length_bra+length_operator+length_ket),      &
 	                       ket_j,ket_nq,ket_nu,ket_subshellJ,fail)      
      if (fail) goto 10
      !
      if (selection == 1) then
         answer(length_bra+length_operator+length_add:        &
	        length_bra+length_operator+length_add) = ","
         ket_nq = bra_nq - 1
      else if(selection == 4) then
         ket_nq = bra_nq
      end if
      ket_j = bra_j
      !
      ! Now distinguish degenerate subshell states for j=9/2 in the ket 
      ! function; if needed, an additional state identifier Nr. will be required
      if ( (ket_j == 9 .and. ket_nu == 4 .and. ket_subshellJ == 8)   .or.   &
           (ket_j == 9 .and. ket_nu == 4 .and. ket_subshellJ == 12)) then
         !
         ! Analys the j=9/2 subshell of ket function. If need 
         ! the code will ask about additional state identifier Nr. 
      11 answer(length_bra+length_operator+length_ket:                      &
	        length_bra+length_ket+length_operator+4)= " Nr= "
         print*, "Enter an additional state identifier Nr = 1 or 2."
	 write(unit=*,fmt="(a)",advance="no")                               &
	                    answer(1:length_bra+length_ket+length_operator+4) 
         read (unit=*,fmt="(a)",err=11)                                     &
	                  answer(length_bra+length_ket+length_operator+5:180)
         ! 
         if (selection == 1) then
	    length_add = index(answer(length_bra+length_ket+length_operator &
                               +5:180),",")
	    if (length_add == 0) then
	       print *, "Input must start with a parenthesis '(' and" //    &
	                " and by a comma ','." 
	       goto 11	
            end if       
         else
	    length_add = index(answer(length_bra+length_ket+                &
                                      length_operator+5:180),")")
	    if (length_add == 0) then
	       print*, "Iinput must start and end with parenthesis '( ...)'."
	       goto 11	
	    end if
         end if
	 ket_nr = get_integer_from_string(answer(                           &
	            length_bra+length_ket+length_operator+5:                &
	            length_bra+length_ket+length_operator+3+length_add),fail) 
	 if (fail) goto 11
	 if (ket_nr > 2 .or. ket_nr < 1) goto 11
  	 length_ket = length_ket + 4 + length_add
      else 
	 ket_nr = 0	   
      end if
      !
      ket_no = rcfp_get_term_number                                            &
                                 (ket_j,ket_nq,ket_nu,ket_nr,ket_subshellJ,fail)
      if (fail) then
         print *, "Unable to recognize the subshell state with (j,nu,J) = (", &
                  ket_j,",",ket_nu,",",ket_subshellJ,"); reenter ..."
         goto 10
      end if
      !
      ! Calculate the required reduced coefficient or matrix element
      select case (selection)
      case (1)
         answer(length_bra+length_ket+length_operator+1:                   &
	        length_bra+length_ket+length_operator+2) = "  "
         answer(length_bra+length_ket+length_operator+3:                   &
	        length_bra+length_ket+length_operator+3) = achar(48+bra_j)
         answer(length_bra+length_ket+length_operator+4:                   &
	        length_bra+length_ket+length_operator+7) = "/2 )"
  	 length_ket = length_ket + 7
         bra_MQ     = bra_nq-(bra_j +1)/2
	 ket_MQ     = ket_nq-(ket_j +1)/2
	 if (abs(bra_MQ) > terms_jj(bra_no)%Q .or.        &
	     mod(terms_jj(bra_no)%Q + bra_MQ,2) /= 0) then
	    coeff = zero
	 else if (abs(ket_MQ) > terms_jj(ket_no)%Q .or.        &
	          mod(terms_jj(ket_no)%Q + ket_MQ,2) /= 0) then	     
	    coeff = zero
	 else
	    coeff = rcfp_Clebsch_Gordan_qusispin                     &
	                           (terms_jj(ket_no)%Q,ket_MQ,       &
                                   1,1,terms_jj(bra_no)%Q, bra_MQ) * &
		    rcfp_get_coefficient(bra_no, ket_no)           / &
		    sqrt(bra_nq * (terms_jj(bra_no)%Q + one)       * &
		    (terms_jj(bra_no)%subshellJ + one))
            if(mod(bra_nq+1,2) /= 0) coeff = - coeff
	 end if
      case (2)
	 coeff = rcfp_get_coefficient(bra_no,ket_no)
      case (3)
	 coeff = rcfp_get_reduced_W(bra_no,ket_no,rank_q,rank_j)
      case (4)
         bra_MQ = bra_nq-(bra_j +1) / 2
	 ket_MQ = ket_nq-(ket_j +1) / 2
	 if (abs(bra_MQ) > terms_jj(bra_no)%Q .or.        &
	     mod(terms_jj(bra_no)%Q + bra_MQ,2) /= 0) then
	    coeff = zero
	 else if (abs(ket_MQ) > terms_jj(ket_no)%Q .or.        &
	          mod(terms_jj(ket_no)%Q + ket_MQ,2) /= 0) then	     
	    coeff = zero
	 else if (mod(rank_k,2) /= 0) then
	    if (bra_nq  /=  ket_nq) then
	       coeff = zero
	    else
	       coeff = -rcfp_get_reduced_W (bra_no, ket_no,0,rank_k)/ &
		        sqrt(two*(terms_jj(bra_no)%Q+one)*(two*rank_k+one))
            end if             
	 else if (rank_k == 0) then
	    if (bra_no /=  ket_no) then
               coeff = zero
	    else if (bra_nq  /=  ket_nq) then
	       coeff = zero
	    else
	       coeff = ket_nq * sqrt((terms_jj(bra_no)%subshellJ + one)/  &
                                     (terms_jj(bra_no)%j + one))
	    end if
	 else
	    coeff = - rcfp_Clebsch_Gordan_qusispin                        &
	                           (terms_jj(ket_no)%Q,ket_MQ,            &
                                   2,0,terms_jj(bra_no)%Q, bra_MQ)      * &
		    rcfp_get_reduced_W (bra_no, ket_no,1,rank_k)        / &
		    sqrt(two *(terms_jj(bra_no)%Q + one)*(two*rank_k +one))
	 end if
      end select
      if(abs(coeff) < eps10) coeff = zero
      write(unit=*,fmt="(a)",advance="no")                                    &
	                        answer(1:length_bra+length_operator+length_ket)
      write(*,"(a,es15.8)") " = ",coeff
      write(unit=*,fmt="(a)",advance="no") "Continue"
      read(*,*)
      go to 1
      !
   end subroutine rcfp_input
   !
   !
   subroutine rcfp_analyse(string,j,nq,nu,subshellJ,fail)
   !-----------------------------------------------------------------------
   ! Analyzes 'string' in order to obtain quantum numbers of symmetry adapted
   ! subshell <bra| or |ket> states. For <bra| states, it returns
   ! the values of angular momentum j, seniority nu, and subshell total 
   ! angular momentum subshellJ; for |ket> states, it returns seniority nu, 
   ! and subshell total angular momentum subshellJ. It returns fail,
   ! if one or more quantum numbers cannot be recognized.
   !
   ! Calls: get_dinteger_from_string().
   !-----------------------------------------------------------------------
      !
      character(len=*), intent(in) :: string
      logical, intent(out)         :: fail 
      integer, intent(out)         :: j, nq, nu, subshellJ
      integer                      :: i
      !
      character(len=80)            :: string_a
      !
      fail  = .false.
      nq    = -100 
      !
      if (index(string,"(") /= 0) then
         !
         ! Determine the angular momentum j value for the <bra| function
         i = index(string,"(")
         string_a = adjustl(string(i+1:))
         if (len_trim(string_a) < 1) then
            fail = .true.
	    print *, "Quantum numbers j, nu and J are missing; reenter ..."
            return
         end if
         j = get_dinteger_from_string(string_a(1:3),fail)
         if (fail) then
            print *, "Unable to decode the angular momentum j: '"    // &
                     string_a(1:4)//"', this must be 1/2, 3/2, 5/2 " // &
                     " or 9/2; reenter ..." 
	    return
         end if
         !
         ! Determine the occupation ^N for the <bra| function
         i        = index(string,"^")
         if (i /= 0) then
            string_a = adjustl(string(i+1:))
            nq       = get_integer_from_string(string_a(1:2),fail)
            if (fail) then
               print *, "Unable to decode the occupation number N: '" // &
                        string_a(i+11:i+2)//"', this must be an integer;"  // &
                        " reenter ..." 
	       return
            end if
            string_a = adjustl(string_a(3:))
	 else
            string_a = adjustl(string_a(4:))	    
         end if
      else if (index(string,")") /= 0) then
         i        = index(string,")")
         string_a = adjustl(string(1:i-1))
      else
	 fail = .true.
	 print*, "Input must either start with symbol '(' or " // &
	         " end with symbol ')'; reenter ..."
         return
      end if
      !
      ! Determine the seniority number nu
      nu = get_integer_from_string(string_a(1:2),fail)
      if (fail) then
         print *, "Unable to decode the seniority number nu: '" // &
                  string_a(1:2)//"', this must be an integer;"  // &
                  " reenter ..." 
         return
      end if
      string_a = adjustl(string_a(3:))
      !
      ! Determine the total subshell angular momentum subshellJ
      subshellJ = get_dinteger_from_string(string_a(1:),fail)
      if (fail) then
         print *, "Unable to decode the subshell angular momentum J: '" // &
                  string_a(1:)//"', this must be an integer;"          // &
                  " reenter ..." 
         return
      end if
      !
   end subroutine rcfp_analyse
   !
   !
   function rcfp_get_term_number(j,nq,nu,Nr,subshellJ,fail)      result(value)
   !-----------------------------------------------------------------------
   ! Returns the internal index for a subshell term which is given by 
   ! its angular momentum 2*j, the seniority quantum number as well as the
   ! total subshell angular momentum 2*J.
   !-----------------------------------------------------------------------
      !
      integer, intent(in)              :: j, nu, nq, Nr, subshellJ
      logical, optional, intent(out)   :: fail
      integer                          :: value
      !       
      integer, dimension(9), parameter :: no_min = (/1,0,3,0, 6,0,12,0,26/)
      integer, dimension(9), parameter :: no_max = (/2,0,5,0,11,0,25,0,63/)
      integer                          :: p_min, p_max, run
      !
      fail = .false.
      if (j <= 9) then
         p_min = no_min(j);   p_max = no_max(j)
         if (rabs_use_stop         .and.  &
            (p_min == 0  .or.  p_max == 0)) then
            stop  "rcfp_get_term_number(): program stop A."
         end if
      end if
      if (j < 9) then
         do run = p_min, p_max
	    if (nu == terms_jj(run)%nu ) then
               if (subshellJ == terms_jj(run)%subshellJ) then
	          value = run
                  return
               end if    
	    end if
         end do 
      else if (j == 9) then 
         if (Nr == -100) then
            if (nq < 3) then
            do run = p_min, p_max      
	       if (nu == terms_jj(run)%nu ) then
                  if (subshellJ == terms_jj(run)%subshellJ) then
	             value = run
	             return
	          end if    
	       end if
            end do
            else
            do run = p_min, p_max      
	       if (nu == terms_jj(run)%Nr ) then
                  if (subshellJ == terms_jj(run)%subshellJ) then
	             value = run
	             return
	          end if    
	       end if
            end do
            end if
	 else
            do run = p_min, p_max      
	       if (nu == terms_jj(run)%nu ) then
                  if (subshellJ == terms_jj(run)%subshellJ) then
                     if (Nr == terms_jj(run)%Nr) then		 
	                value = run
	                return
                     end if
	          end if    
	       end if
            end do
         end if
      else
         if (nq == 0) then
            if (subshellJ == 0) then
	       if (nu == 0) then
	          value = ((j * 1000) + nu) * 1000 + subshellJ
	          return
               end if
            end if
         else if (nq == 1) then
            if (subshellJ == j) then
	       if (nu == 1) then
	          value = ((j * 1000) + nu) * 1000 + subshellJ
	          return
               end if
            end if
         else if (nq == 2) then
            do run = 0, 2*j-2, 4
               if (subshellJ == run) then
                  if (subshellJ == 0) then
	             if (nu == 0) then
	                value = ((j * 1000) + nu) * 1000 + subshellJ
	                return
                     end if
                  else
	             if (nu == 2) then
	                value = ((j * 1000) + nu) * 1000 + subshellJ
	                return
                     end if
                  end if
               end if
	    end do
         end if
         fail = .true.
      end if
      !
      if (present(fail)) then
         fail = .true.
      else if (rabs_use_stop) then
         stop  "rcfp_get_term_number(): program stop B."
      end if
      !
   end function rcfp_get_term_number
   !
   !
   function rcfp_Clebsch_Gordan_qusispin(ja,ma,jb,mb,Jab,Mab)     result(CG)
   !-----------------------------------------------------------------------
   ! Calculates specific the Clebsch-Gordan coefficient 
   ! which need in quasispin formalism.
   !-----------------------------------------------------------------------
      !
      integer, intent(in)          :: ja, jb, Jab, ma, mb, Mab
      real(kind=dp)                :: CG
      ! real (kind=dp), dimension(2) :: factor
      !
      if (rabs_use_stop   .and.   mod(ja-jb+Mab,2) /= 0) then
         stop "rcfp_Clebsch_Gordan_qusispin(): program stop A."
      end if
      !
      ! Test the triangular condition and that for magnetic quantum numbers
      !
      CG = zero
      if (ma+mb /= Mab) then
	 return
      else if (triangle(ja+1,jb+1,Jab+1) == 0) then
	return
      else if (abs(ma)>ja .or. abs(mb)>jb .or. abs(Mab)>Jab) then
	return
      else if (mod(abs(ma)+ja,2) /= 0 .or. mod(abs(mb)+jb,2) /= 0 .or.    &
	       mod(abs(Mab)+Jab,2) /= 0) then
	return
      else if (rabs_use_stop                    .and. &
              (mod(ma+ja+ja,2)    /= mod(ja,2)  .or.  &
               mod(mb+jb+jb,2)    /= mod(jb,2)  .or.  &
               mod(Mab+Jab+Jab,2) /= mod(Jab,2)) ) then
         stop "rcfp_Clebsch_Gordan_qusispin(): program stop B."
      end if
      !
      ! Calculate
      select case (jb)
      case (0)
         if (ja /= Jab) return
         if (ma /= Mab) return
         CG = one
      case (1)
         if (ja + 1 == Jab) then
            CG = sqrt(half*(Jab+mb*Mab)/Jab)
         else if (ja - 1 == Jab) then
            CG = -mb*sqrt(half*(Jab-mb*Mab+two)/(Jab+two))
         end if
      case (2)
         if (mb == 0) then
            if (ja + 2 == Jab) then
               CG = sqrt((half*(Jab+Mab)*(Jab-Mab))/((Jab-one)*Jab))
            else if (ja == Jab) then
               CG = (half*Mab)/(half*sqrt(Jab*(Jab+two)))
            else if (ja - 2 == Jab) then
               CG = -sqrt(half*((Jab+Mab+two)*(Jab-Mab+two))/              &
	                       ((Jab+two)*(Jab+three)))
            end if
         else if (mb == 2  .or.  mb == -2) then
            if (ja + 2 == Jab) then
               CG = half*sqrt(((Jab+mb*Mab*half-two)*(Jab+mb*half*Mab))/       &
	                      ((Jab-one)*Jab))
            else if (ja == Jab) then
               CG =-mb*half*sqrt(half*((Jab-mb*Mab*half+two)*(Jab+mb*Mab*half))&
	                             /((Jab+two)*Jab))
            else if (ja - 2 == Jab) then
               CG = half*sqrt(((Jab-mb*Mab*half+two)*(Jab-mb*Mab*half+four))/  &
	                       ((Jab+two)*(Jab+three)))
            end if
         end if
      case default
         CG =  Clebsch_Gordan(ja,ma,jb,mb,Jab,Mab)
      end select
      !
   end function rcfp_Clebsch_Gordan_qusispin
   !
end module rabs_rcfp
