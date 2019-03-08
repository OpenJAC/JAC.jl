module rabs_anco
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module controls the computation of `angular coefficients' within the RATIP environment. At the present, angular 
! coefficients may be obtained for all scalar (and symmetric) two-particle interactions. They are sometimes called 
! 'pure angular coefficients' in order to distinguish them from earlier approaches in GRASP (and its derivatives) where 
! part of the `physical interaction' were included into the angular integrals. The 'pure' angular coefficients generally 
! appear in the following form:
!
!    <CSF_l | operator | CSF_r> = sum_{t} coeff_t (L,abcd) * X^L (abcd)
!
! where CSF is a standard (jj-coupled) configuration state function, t a summation index, and X^L(abcd) denotes the 
! effective interaction strength. This effective strength is the only part which depends on the (scalar) operator. 
! Thus, the `pure' coefficients can be applied to any interaction.
!
! The computation of the angular coefficients is based on the quasi-spin concept and the reduced coefficients of 
! fractional parentage (rcfp). For details, see G. Gaigalas et al., J. Phys. B30 (1997) 3747.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_rcfp
   use rabs_recoupling
   implicit none
   !
   public  :: anco_calculate_csf_block
                 ! Calculates the angular coefficients for a given CSF scheme.
   public  :: anco_calculate_csf_matrix
                 ! Calculates the angular coefficients for a sub-matrix of a given CSF scheme.
   public  :: anco_calculate_csf_matrix_1p
                 ! Calculates the angular coefficients for a one-particle operator A of given rank nu for a sub-matrix 
                 ! of a CSF scheme.
   public  :: anco_calculate_csf_pair
                 ! Calculates the angular coefficients for a given pair of CSF.
   public  :: anco_calculate_csf_pair_1p
                 ! Calculates the angular coefficients for a one-particle operator A of given rank nu for a given 
                 ! pair of a CSF.
   public  :: anco_case_1
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distribution  1  (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_2_to_5
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 2 - 5 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_6
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions  6 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_7_to_14
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 11 - 14 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_7_to_8
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 7 - 8 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_9_to_10
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 9 - 10 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_11_to_14
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 11 - 14 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_15_to_18
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 15 - 18 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_15_to_18_order
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 15 - 18 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_19_to_42
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 19 - 26 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_19_to_26
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 19 - 26 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_27_to_34
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 27 - 34 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_case_35_to_42
                 ! Calculates the angular coefficients for a given pair r, s of CSF for distributions 35 - 42 (of Table 1 in 
                 ! G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747).
   public  :: anco_codeing
                 ! Returns subshell states as appropriate to call routines from the module rabs_rcfp.
   public  :: anco_collect_input
                 ! Collect and proceeds all input for the xanco program.
   private :: anco_diagonal_angle_s
                 ! Calculates the angular coefficients of none-scalar one-particle operator for a given pair r, s of CSF
                 ! for diagonal matrix elements with respect to configurations.
   private :: anco_diff_occ_2
                 ! Calculates the angular coefficients for a given pair r, s of CSF for non-diagonal matrix elements with 
                 ! diff_occ = 2.
   public  :: anco_normal_form
                 ! Found the normal form of creation and annihilation operators.
   public  :: anco_normal_phase
                 ! Determinate the phase factor wich appear from permutation of operators of second quantization.
   private :: anco_one_particle_diag
                 ! Calculates the angular coefficients for a given pair r, s of CSF for diagonal matrix elements with 
                 ! respect to configurations for one-particle non scalar operator.
   private :: anco_one_particle_off
                 ! Calculates the angular coefficients for a given pair r, s of CSF for off diagonal matrix elements 
                 ! with respect to configurations for one-particle non scalar operator
   public  :: anco_one_particle_scalart
                 ! Calculates the angular coefficients for a given pair r, s for one-particle scalar operator.
   public  :: anco_open_vnu
                 ! Opens a .vnu V^nu (abcd) Angular Coefficient file for the ANCO program on stream 25.
   public  :: anco_print_results
                 ! Prints a summary of the calculation in a neat format to the .sum file.
   public  :: anco_print_summary
                 ! Appends a summary of the input data to the .sum file.
   !
   ! Define storage for the CSF basis to read in from a .csl file
   type(asf_basis), public :: asf
   !
   ! Define an internal structure to stores the 'pure' one- and two-particle
   ! coefficients
   type, bind(c) :: anco_T_coeff
      integer          :: nu    ! Should always be zero for scalar interaction.
      type(nkappa)     :: a, b
      real(kind=dp)    :: T
   end type anco_T_coeff
   !
   type, bind(c) :: anco_V_coeff
      integer          :: nu
      type(nkappa)     :: a, b, c, d
      real(kind=dp)    :: V
   end type anco_V_coeff
   !
   type :: anco_csf_pair
      integer          :: r, s
      integer          :: no_T_coeff, no_V_coeff
      type(anco_T_coeff), dimension(:), pointer :: T_coeff
      type(anco_V_coeff), dimension(:), pointer :: V_coeff
   end type anco_csf_pair
   !
   type(anco_T_coeff), dimension(1000), bind(c) :: anco_T_list
   type(anco_V_coeff), dimension(2000), bind(c) :: anco_V_list
   !
   integer :: number_of_pair_list       = 0,     &
              number_of_pair_list_alloc = 0,     &
              number_of_pair_list_max   = 10000, &
              anco_one_particle_rank
   !
   type(anco_csf_pair), dimension(:), allocatable   :: anco_pair_list
   !
   type, private  :: anco_temp
      real(kind=dp), dimension(50) :: recoupling
      real(kind=dp), dimension(50) :: operator
   end type anco_temp
   !
   ! Define global logical flags for the control of the ANCO program; the default values for some of these flags may 
   ! be overwritten interactively during input time.
   logical, public ::   anco_peel_shells_only             = .false.,  &
                        anco_one_particle                 = .true.,   &
                        anco_one_particle_nonscalar       = .true.,   &
                        anco_pure_one_particle            = .true.,   &
                        anco_pure_one_particle_nonscal    = .true.,   &
                        anco_two_particle                 = .true.,   &
                        anco_pure_two_particle            = .true.
   !
   real(kind=dp), private :: anco_total, anco_total_T, anco_total_Vnu
   real(kind=dp), private :: total_number_pair, numerator = zero
   real(kind=dp), private :: numerator_print              = 1000000.0
   !
contains
   !
   !
   subroutine anco_calculate_csf_block(csf_set)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for all rows (>= columns) and 
   ! columns of the given CSF scheme as defined in csf_set.
   !
   ! Calls: anco_calculate_csf_pair().
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(in) :: csf_set
      !
      character(len=4) :: o_a, o_b, o_c, o_d
      integer          :: i, no_T_coeff, no_V_coeff, r, s,                    &
                          ja, jb, jc, jd, la, lb, lc, ld, nu
      real(kind=dp)    :: xc
      !
      total_number_pair = csf_set%nocsf *(csf_set%nocsf + one)/two
      print *, "Calculate 'pure' angular coefficients for a total number of ",&
               total_number_pair," pairs of CSF ..."
      !
      do  s = 1,csf_set%nocsf
         do  r = s,csf_set%nocsf
            if (anco_one_particle .or. anco_two_particle) then
               if (.not. anco_one_particle_nonscalar) then
                 call anco_calculate_csf_pair(csf_set,r,s,no_T_coeff,no_V_coeff)
               end if
            else if (anco_one_particle_nonscalar) then
               if (anco_pure_one_particle)then
                 call anco_calculate_csf_pair_1p(anco_one_particle_rank,  &
                                                 csf_set,r,s,no_T_coeff)
               else
                 call anco_calculate_csf_pair_1p(anco_one_particle_rank,  &
                                                 csf_set,s,r,no_T_coeff)
               end if
            else
               stop "anco_calculate_csf_block(): program stop A."
            end if
            if (no_T_coeff > 0) then
               do  i = 1,no_T_coeff
                  anco_total   = anco_total + 1
                  anco_total_T = anco_total_T + one
                  ja  = angular_momentum_j(anco_T_list(i)%a%kappa)
                  jb  = angular_momentum_j(anco_T_list(i)%b%kappa)
                  nu  = anco_T_list(i)%nu
                  o_a = orbital_name(anco_T_list(i)%a%n,anco_T_list(i)%a%kappa)
                  o_b = orbital_name(anco_T_list(i)%b%n,anco_T_list(i)%b%kappa)
                  if (.not. anco_one_particle_nonscalar) then
                     if (anco_pure_one_particle)then
                        write(25,"(a,i3,a,i3,5a,es20.12)")                     &
                        "  pure one-particle [",r,",",s,"] (",o_a,",",o_b,")=",&
                         anco_T_list(i)%T
                     else
                        anco_T_list(i)%T = sqrt(ja + one) * anco_T_list(i)%T
                        write(25,"(a,i5,a,i5,5a,es20.12)")                     &
                        "  T^[",r,",",s,"] (",o_a,",",o_b,") =",anco_T_list(i)%T
                     end if
                  else if (anco_one_particle_nonscalar) then
                     if (anco_pure_one_particle)then
                        write(25,"(a,i3,a,i3,a,i2,5a,es20.12)")                &
                        "  pure one-particle [",r,",",s,"]  nu=",nu,           &
                        " (",o_a,",",o_b,") =",anco_T_list(i)%T
                     else
                        anco_T_list(i)%T = sqrt(ja+one)* anco_T_list(i)%T
                        write(25,"(a,i5,a,i5,a,i2,5a,es20.12)")                &
                        "  T^[",s,",",r,"]  nu=",nu,                           & 
                        " (",o_a,",",o_b,") =",anco_T_list(i)%T
                     end if
                  end if
               end do
            end if
            !
            if (no_V_coeff > 0) then
               do  i = 1,no_V_coeff
                  nu = anco_V_list(i)%nu
                  ja = angular_momentum_j(anco_V_list(i)%a%kappa)
                  jb = angular_momentum_j(anco_V_list(i)%b%kappa)
                  jc = angular_momentum_j(anco_V_list(i)%c%kappa)
                  jd = angular_momentum_j(anco_V_list(i)%d%kappa) 
                  !
                  la = angular_momentum_l(anco_V_list(i)%a%kappa)
                  lb = angular_momentum_l(anco_V_list(i)%b%kappa) 
                  lc = angular_momentum_l(anco_V_list(i)%c%kappa)
                  ld = angular_momentum_l(anco_V_list(i)%d%kappa) 
                  !
                  o_a = orbital_name(anco_V_list(i)%a%n,anco_V_list(i)%a%kappa)
                  o_b = orbital_name(anco_V_list(i)%b%n,anco_V_list(i)%b%kappa)
                  o_c = orbital_name(anco_V_list(i)%c%n,anco_V_list(i)%c%kappa)
                  o_d = orbital_name(anco_V_list(i)%d%n,anco_V_list(i)%d%kappa)
                  if (anco_pure_two_particle)then
                     anco_total = anco_total + 1
                     anco_total_Vnu = anco_total_Vnu + one
                     write(25,"(a,i1,a,i3,a,i3,a,4a,3a,a,es20.12)")            &
                     "  pure two-particle [( ",anco_V_list(i)%nu,")]_[",r,",", &
                     s,"] (",o_a,",",o_b,";",o_c,",",o_d,") =",anco_V_list(i)%V
                  else
                     if (triangle(ja+1,jc+1,nu+nu+1) *                         &
                        triangle(jb+1,jd+1,nu+nu+1)== 0  .or.                  &
                        mod(la+lc+nu,2) == 1 .or. mod(lb+ld+nu,2) == 1) then
                        cycle
                     end if               
                     xc= CL_reduced_me                                         &
                         (anco_V_list(i)%a%kappa,nu,anco_V_list(i)%c%kappa)*   &
                         CL_reduced_me                                         &
                         (anco_V_list(i)%b%kappa,nu,anco_V_list(i)%d%kappa)
                     if (mod(nu,2) == 1) then
                        xc = - xc!
                     end if
                     !
                     if (abs(xc) .gt. eps10) then
                        anco_total       = anco_total + 1
                        anco_total_Vnu   = anco_total_Vnu + one
                        anco_V_list(i)%V = xc * anco_V_list(i)%V
                        write(25,"(a,i2,a,i5,a,i5,a,4a,3a,a,es20.12)")         &
                        "  V^[(",anco_V_list(i)%nu,")]_[", r,",", s,"] (",     &
                        o_a,",", o_b,";", o_c,",", o_d,") =",anco_V_list(i)%V
                     end if
                  end if
               end do
            end if
           !!x end if
         end do
      end do
      !
      print *, " Total number of calculated coefficients =  ",nint(anco_total)
      print *, " "
      if (anco_pure_one_particle   .and.   anco_total_T > 0) then
         print *, " Total number of 'pure' T coefficients           =  ", &
                                                               anco_total_T
      else if (anco_total_T > 0) then
         print *, " Total number of GRASP92-like T coefficients     =  ", &
                                                               anco_total_T
      end if
      !
      if (anco_pure_two_particle   .and.   anco_total_Vnu > 0) then
         print *, " Total number of 'pure' V^L (abcd)  coefficients =  ", &
                                                             anco_total_Vnu
      else if (anco_total_Vnu > 0) then
         print *, " Total number of V^L_Coulomb (abcd) coefficients =  ", &
                                                             anco_total_Vnu
      end if
      !
   end subroutine anco_calculate_csf_block
   !
   !
   subroutine anco_calculate_csf_matrix(csf_set,row_low,row_up,col_low,col_up)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a sub-matrix with rows from
   ! row_low ... row_up and columns from col_low ... col_up of the
   ! CSF scheme as defined in csf_set.
   !
   ! Calls: anco_calculate_csf_pair().
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: row_low, row_up, col_low, col_up
      type(csf_basis), intent(in) :: csf_set
      !
      integer :: no_T_coeff, no_V_coeff, r, s, max_number
      !
      ! Modified for tests
      if (allocated(anco_pair_list)) then
         do  r = 1,number_of_pair_list
            deallocate( anco_pair_list(r)%V_coeff )
            deallocate( anco_pair_list(r)%T_coeff )
         end do
         deallocate( anco_pair_list )
      end if
      !
      max_number = (row_up - row_low + 1) * (col_up - col_low + 1)
      allocate( anco_pair_list(1:max_number))
      number_of_pair_list_alloc = max_number
      !
      number_of_pair_list = 0
      if (rabs_use_stop     .and.                                  &
         (row_low < 1       .or.   row_up > csf_set%nocsf   .or.   &
          col_low < 1       .or.   col_up > csf_set%nocsf   .or.   &
          row_low > row_up  .or.   col_low > col_up)) then
         stop "anco_calculate_csf_matrix(): program stop A."
      end if
      !
      do  s = col_low, col_up
         do  r = row_low, row_up
            call anco_calculate_csf_pair(csf_set,r,s,no_T_coeff,no_V_coeff)
            if (no_V_coeff > 0) then
               number_of_pair_list = number_of_pair_list + 1
               allocate(                                                    &
                  anco_pair_list(number_of_pair_list)%T_coeff(1:no_T_coeff),&
                  anco_pair_list(number_of_pair_list)%V_coeff(1:no_V_coeff) )
               anco_pair_list(number_of_pair_list)%r          = r
               anco_pair_list(number_of_pair_list)%s          = s
               anco_pair_list(number_of_pair_list)%no_T_coeff = no_T_coeff
               anco_pair_list(number_of_pair_list)%no_V_coeff = no_V_coeff
               anco_pair_list(number_of_pair_list)%T_coeff(1:no_T_coeff) =  &
                                                    anco_T_list(1:no_T_coeff)
               anco_pair_list(number_of_pair_list)%V_coeff(1:no_V_coeff) =  &
                                                    anco_V_list(1:no_V_coeff)

               !! write(65,"(3i5)")                                         &
               !!         r,s,anco_pair_list(number_of_pair_list)%no_V_coeff
               !! do i = 1,anco_pair_list(number_of_pair_list)%no_T_coeff
               !!    write(65,'(5i5,f17.7)')                                &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%nu,     &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%a,      &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%b,      &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%T
               !! end do
               !! do i = 1,anco_pair_list(number_of_pair_list)%no_V_coeff
               !!    write(65,'(9i5,f17.7)')                                &
               !!    anco_pair_list(number_of_pair_list)%V_coeff(i)%nu,     &
               !!    anco_pair_list(number_of_pair_list)%V_coeff(i)%a,      &
               !!    anco_pair_list(number_of_pair_list)%V_coeff(i)%b,      &
               !!    anco_pair_list(number_of_pair_list)%V_coeff(i)%c,      &
               !!    anco_pair_list(number_of_pair_list)%V_coeff(i)%d,      &
               !!    anco_pair_list(number_of_pair_list)%V_coeff(i)%V
               !! end do
            end if
         end do
      end do
      !
   end subroutine anco_calculate_csf_matrix
   !
   !
   subroutine anco_calculate_csf_matrix_1p(nu,csf_set,row_low,row_up,   &
                                                         col_low,col_up)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a one-particle operator A of
   ! rank nu for a sub-matrix with rows from row_low ... row_up and
   ! columns from col_low ... col_up of the CSF scheme as defined in
   ! csf_set.
   !
   ! Calls: anco_calculate_csf_pair_1p().
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: nu, row_low, row_up, col_low, col_up
      type(csf_basis), intent(in) :: csf_set
      !
      integer :: no_T_coeff, r, s, max_number
      !
      max_number = (row_up - row_low + 1) * (col_up - col_low + 1)
      allocate( anco_pair_list(1:max_number))
      !
      number_of_pair_list = 0
      if (rabs_use_stop     .and.                                 &
         (row_low < 1       .or.   row_up > csf_set%nocsf   .or.  &
          col_low < 1       .or.   col_up > csf_set%nocsf   .or.  &
          row_low > row_up  .or.   col_low > col_up)) then
         stop "anco_calculate_csf_matrix_1p(): program stop A."
      end if
      !
      do  s = col_low, col_up
         do  r = row_low, row_up
            call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
            if (no_T_coeff > 0) then
               number_of_pair_list = number_of_pair_list + 1
               allocate(                                                    &
               anco_pair_list(number_of_pair_list)%T_coeff(1:no_T_coeff))
               anco_pair_list(number_of_pair_list)%r          = r
               anco_pair_list(number_of_pair_list)%s          = s
               anco_pair_list(number_of_pair_list)%no_T_coeff = no_T_coeff
               anco_pair_list(number_of_pair_list)%T_coeff(1:no_T_coeff) =  &
                                                    anco_T_list(1:no_T_coeff)
               !! write(65,"(3i5)")                                         &
               !!         r,s,anco_pair_list(number_of_pair_list)%no_T_coeff
               !! do i = 1,anco_pair_list(number_of_pair_list)%no_T_coeff
               !!    write(65,'(5i5,f17.7)')                                &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%nu,     &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%a,      &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%b,      &
               !!    anco_pair_list(number_of_pair_list)%T_coeff(i)%T
               !! end do
            end if
         end do
      end do
      !
   end subroutine anco_calculate_csf_matrix_1p
   !
   !
   subroutine anco_calculate_csf_pair(csf_set,r,s,no_T_coeff,no_V_coeff)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF.
   !
   ! Calls:  anco_diagonal_angle_s, anco_diff_occ_2, anco_case_6,
   !         anco_case_15_to_18, anco_case_19_to_42
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: r, s
      type(csf_basis), intent(in) :: csf_set
      integer, intent(out)        :: no_T_coeff, no_V_coeff
      !
      integer :: count, diff, diff_occ, nshell
      integer :: creation_one, creation_two, annihilation_one, annihilation_two
      !
      no_T_coeff = 0;   no_V_coeff = 0
      !
      numerator = numerator + one
      if (modulo(numerator,numerator_print) == 0) &
         print *, " ANCO calculates the pair", numerator
      !
      if (csf_set%csf(r)%totalJ /= csf_set%csf(s)%totalJ) return
      diff_occ = 0;   count = 0
      creation_one     = 0;   creation_two     = 0
      annihilation_one = 0;   annihilation_two = 0
      do  nshell = csf_set%nwcore + 1,csf_set%nwshells
         diff = csf_set%csf(r)%occupation(nshell) -  &
                csf_set%csf(s)%occupation(nshell)
         if (iabs(diff) > 2) then
            return
         else if (diff /= 0) then
            diff_occ = diff_occ + iabs(diff);   count = count + 1
         end if
         !
         if (diff == 1) then
            if (creation_one == 0) then
                creation_one = nshell
            else if (creation_two == 0) then
                creation_two = nshell
            else
               return
            end if
         else if (diff == 2) then
            if (creation_one == 0) then
                creation_one = nshell
                creation_two = nshell
            else
               return
            end if
         else if (diff == -1) then
            if (annihilation_one == 0) then
                annihilation_one = nshell
            else if (annihilation_two == 0) then
                annihilation_two = nshell
            else
               return
            end if
         else if (diff == -2) then
            if (annihilation_one == 0) then
                annihilation_one = nshell
                annihilation_two = nshell
            else
               return
            end if
         end if
      end do
      !
      if (rabs_use_stop   .and.   mod(diff_occ,2) /= 0) then
         stop "anco_calculate_csf_pair(): program stop A."
      else
         select case(count)
         case (0)     
            call anco_diagonal_angle_s(csf_set%csf(r),csf_set%csf(s), &
                 csf_set%nwshells,csf_set%nwcore,csf_set%subshell,    &
                                                 no_T_coeff,no_V_coeff)
         case (2)
            if (diff_occ == 2) then
            call anco_diff_occ_2(csf_set%csf(r),csf_set%csf(s),       &
                 csf_set%nwshells,csf_set%nwcore,csf_set%subshell,    &
                   no_T_coeff,no_V_coeff,creation_one,annihilation_one)
            else if (diff_occ == 4 .and. anco_two_particle) then
               call anco_case_6(csf_set%csf(r),csf_set%csf(s),        &
                    csf_set%nwshells,csf_set%nwcore,csf_set%subshell, &
                              no_V_coeff,creation_one,annihilation_one)
            end if
         case (3)
            if (anco_two_particle) then
               call anco_case_15_to_18(csf_set%csf(r),csf_set%csf(s),        &
                 csf_set%nwshells,csf_set%nwcore,csf_set%subshell,no_V_coeff,&
                  creation_one,creation_two,annihilation_one,annihilation_two)
            end if
         case (4)
            if (anco_two_particle) then
               call anco_case_19_to_42(csf_set%csf(r),csf_set%csf(s),        &
                 csf_set%nwshells,csf_set%nwcore,csf_set%subshell,no_V_coeff,&
                  creation_one,creation_two,annihilation_one,annihilation_two)
            end if
         end select
      end if
      !
   end subroutine anco_calculate_csf_pair
   !
   !
   subroutine anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients of non-scalar one-particle
   ! operator for a given pair r, s of CSF.
   !
   ! Calls:  anco_diagonal_angle_s_1p
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: nu, r, s
      type(csf_basis), intent(in) :: csf_set
      integer, intent(out)        :: no_T_coeff
      !
      integer :: count, diff, diff_occ, nshell
      integer :: creation_one, annihilation_one
      !
      no_T_coeff = 0;
      !
      numerator = numerator + one
      if (modulo(numerator,numerator_print) == 0) &
         print *, " ANCO calculates the pair", numerator
      !
      if (triangle(csf_set%csf(r)%totalJ+1,2*nu+1, &
                   csf_set%csf(s)%totalJ+1) == 0)  return
      diff_occ = 0;   count = 0
      creation_one     = 0;   annihilation_one = 0;
      do  nshell = csf_set%nwcore + 1,csf_set%nwshells
         diff = csf_set%csf(r)%occupation(nshell) - &
                csf_set%csf(s)%occupation(nshell)
         if (iabs(diff) > 1) then
            return
         else if (diff /= 0) then
            diff_occ = diff_occ + iabs(diff);   count = count + 1
         end if
         !
         if (diff == 1) then
            creation_one = nshell
         else if (diff == -1) then
            annihilation_one = nshell
         end if
      end do
      !
      if (rabs_use_stop   .and.   mod(diff_occ,2) /= 0) then
         stop "anco_calculate_csf_pair_1p(): program stop A."
      else
         select case(count)
         case (0)     
            call anco_one_particle_diag(nu,csf_set%csf(r),csf_set%csf(s),   &
                 csf_set%nwshells,csf_set%nwcore,csf_set%subshell,no_T_coeff)
         case (2)
            if (diff_occ == 2) then
               call anco_one_particle_off(nu,csf_set%csf(r),csf_set%csf(s), &
               csf_set%nwshells,csf_set%nwcore,csf_set%subshell,            &
                                    no_T_coeff,creation_one,annihilation_one)
            end if
         end select
      end if
      !
   end subroutine anco_calculate_csf_pair_1p
   !
   !
   subroutine anco_case_1(csf_r,csf_s,nwshells,nwcore,subshell,no_T_coeff,  &
                                                           no_V_coeff,no_one)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distribution  1  (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 1.         Alpha   Alpha   Alpha   Alpha
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, rcfp_calculate_Wk_me,
   !         rcfp_calculate_Wk_times_Wk_0_me
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one
      integer, intent(inout)                 :: no_T_coeff, no_V_coeff
      !
      type(subshell_state) :: bra, ket
      real(kind=dp)        :: coeff, coeff_one
      integer              :: j_1, rank, rank_min
      !
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_one,no_one,no_one,no_one,&
                                            nwshells,nwcore,0)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra,ket)
      j_1 =  iabs(subshell(no_one)%kappa) * 2 - 1
      coeff_one = rcfp_calculate_Wk_me(bra,ket,0,1,-1)
      if (.not. anco_peel_shells_only .and. anco_one_particle) then
         if (abs(coeff_one) > eps10) then
            no_T_coeff                      = no_T_coeff + 1
            anco_T_list(no_T_coeff)%nu      = 0
            anco_T_list(no_T_coeff)%a%kappa = subshell(no_one)%kappa
            anco_T_list(no_T_coeff)%a%n     = subshell(no_one)%n
            anco_T_list(no_T_coeff)%b%kappa = subshell(no_one)%kappa
            anco_T_list(no_T_coeff)%b%n     = subshell(no_one)%n
            anco_T_list(no_T_coeff)%T       = -coeff_one                      &
                                            / sqrt(csf_r%subshellJ(no_one)+one)
         end if
      end if
      if (.not. anco_two_particle) return
      coeff_one = coeff_one / sqrt(j_1 + one)
      if (anco_peel_shells_only) then
         rank_min = 1
      else
         rank_min = 0
      end if
      do rank = rank_min, j_1
         coeff = coeff_one
         if (mod(j_1 + rank, 2) /= 0) coeff = -coeff
         coeff = half*((rcfp_calculate_Wk_times_Wk_0_me(bra,ket,rank,1,-1,1,-1)&
               /sqrt(two*rank+one)) - coeff) / sqrt(csf_r%subshellJ(no_one)+one)
         if (abs(coeff) > eps10) then
            no_V_coeff                      = no_V_coeff + 1
            anco_V_list(no_V_coeff)%nu      = rank
            anco_V_list(no_V_coeff)%a%kappa = subshell(no_one)%kappa
            anco_V_list(no_V_coeff)%a%n     = subshell(no_one)%n
            anco_V_list(no_V_coeff)%b%kappa = subshell(no_one)%kappa
            anco_V_list(no_V_coeff)%b%n     = subshell(no_one)%n
            anco_V_list(no_V_coeff)%c%kappa = subshell(no_one)%kappa
            anco_V_list(no_V_coeff)%c%n     = subshell(no_one)%n
            anco_V_list(no_V_coeff)%d%kappa = subshell(no_one)%kappa
            anco_V_list(no_V_coeff)%d%n     = subshell(no_one)%n
            anco_V_list(no_V_coeff)%V       = coeff
         end if
      end do
      !
   end subroutine anco_case_1
   !
   !
   subroutine anco_case_2_to_5(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                                                  no_one,no_two)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 2 - 5 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 2.         Alpha   Beta    Alpha   Beta
   ! 3.         Beta    Alpha   Beta    Alpha
   ! 4.         Alpha   Beta    Beta    Alpha
   ! 5.         Bate    Alpha   Alpha   Beta
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_2_shells,
   !         rcfp_calculate_Wk_me, wigner_6j_triangle, wigner_6j_symbol 
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one
      integer              :: j_1, j_2, rank, delta_J
      integer              :: i_max, i_sum, i_min_2, i_max_2
      !
      if (no_one == no_two) return
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_one,no_two,no_two,no_two, &
                                             nwshells,nwcore,1)) < eps10) return
      !
      ! cases 1212  + + - -     transform to 1122  + - + -
      !       2121                           1122
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      j_1 =  iabs(subshell(no_one)%kappa) * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa) * 2 - 1
      i_max = min(j_1,j_2)
      if (i_max + 1 > 20) stop "anco_case_2_to_5(): program stop A."
      do rank = 0, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_2_shells                &
                               (csf_r,csf_s,2*rank,no_one,no_two,.true.,delta_J)
         if (delta_J == 0) then
            temp%recoupling(rank + 1) = zero;   temp%operator(rank + 1) = zero
         else
            temp%operator(rank + 1)=                                           &
                               rcfp_calculate_Wk_me(bra_one,ket_one,rank,1,-1)*&
                               rcfp_calculate_Wk_me(bra_two,ket_two,rank,1,-1)
            if (abs(temp%operator(rank + 1)) < eps10) then
               temp%recoupling(rank + 1) = zero
            else
               temp%recoupling(rank + 1) = recoupling_matrix_2_shells          &
                              (csf_r,csf_s,2*rank,no_one,no_two,.false.,delta_J)
               coeff = temp%recoupling(rank + 1)                              *&
                                temp%operator(rank + 1) / sqrt(two * rank + one)
               if (abs(coeff) > eps10) then
                  if (anco_peel_shells_only) then
                     if (rank > 0) then
                        nocoeff                      = nocoeff + 1
                        anco_V_list(nocoeff)%nu      = rank
                        anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
                        anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
                        anco_V_list(nocoeff)%b%kappa = subshell(no_two)%kappa
                        anco_V_list(nocoeff)%b%n     = subshell(no_two)%n
                        anco_V_list(nocoeff)%c%kappa = subshell(no_one)%kappa
                        anco_V_list(nocoeff)%c%n     = subshell(no_one)%n
                        anco_V_list(nocoeff)%d%kappa = subshell(no_two)%kappa
                        anco_V_list(nocoeff)%d%n     = subshell(no_two)%n
                        anco_V_list(nocoeff)%V       = coeff
                     end if
                  else
                     nocoeff                      = nocoeff + 1
                     anco_V_list(nocoeff)%nu      = rank
                     anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
                     anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
                     anco_V_list(nocoeff)%b%kappa = subshell(no_two)%kappa
                     anco_V_list(nocoeff)%b%n     = subshell(no_two)%n
                     anco_V_list(nocoeff)%c%kappa = subshell(no_one)%kappa
                     anco_V_list(nocoeff)%c%n     = subshell(no_one)%n
                     anco_V_list(nocoeff)%d%kappa = subshell(no_two)%kappa
                     anco_V_list(nocoeff)%d%n     = subshell(no_two)%n
                     anco_V_list(nocoeff)%V       = coeff
                  end if
               end if
            end if
         end if
      end do
      !
      ! cases 1221  + + - -     transform to 1122  + - + -
      !       2112                           1122
      i_min_2 = iabs(j_1 - j_2)/2;   i_max_2 = (j_1 + j_2)/2
      do rank = i_min_2, i_max_2
         coeff = zero
         do i_sum = 0, i_max
            coeff_one = temp%recoupling(i_sum + 1) * temp%operator(i_sum + 1)
            if (abs(coeff_one) > eps10) then
               if (wigner_6j_triangle(j_1,j_2,2*rank,j_2,j_1,2*i_sum)/= 0) then
                  coeff = coeff + coeff_one * sqrt(two * i_sum + one)         *&
                         wigner_6j_symbol(j_1,j_2,2*rank,j_2,j_1,2*i_sum,.true.)
               end if
            end if
         end do
         if (abs(coeff) > eps10) then
            nocoeff                      = nocoeff + 1
            anco_V_list(nocoeff)%nu      = rank
            anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
            anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
            anco_V_list(nocoeff)%b%kappa = subshell(no_two)%kappa
            anco_V_list(nocoeff)%b%n     = subshell(no_two)%n
            anco_V_list(nocoeff)%c%kappa = subshell(no_two)%kappa
            anco_V_list(nocoeff)%c%n     = subshell(no_two)%n
            anco_V_list(nocoeff)%d%kappa = subshell(no_one)%kappa
            anco_V_list(nocoeff)%d%n     = subshell(no_one)%n
            anco_V_list(nocoeff)%V       = coeff
         end if
      end do
      !
   end subroutine anco_case_2_to_5
   !
   !
   subroutine anco_case_6(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff,        &
                                                                  no_one,no_two)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions  6 (of Table 1 in G. Gaigalas et al., 1997 
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 6.         Alpha   Alpha   Beta    Beta
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_2_shells,
   !         rcfp_calculate_Wk_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one
      integer              :: j_1, j_2, delta_J, rank
      integer              :: no_form_one, no_form_two
      integer              :: i_max, i_sum, i_min_2, i_max_2
      !
      if (nwshells == 1) return
      no_form_one = min(no_one,no_two)
      no_form_two = max(no_one,no_two)
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_form_one,no_form_two,     &
                     no_form_two,no_form_two,nwshells,nwcore,1)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      j_1 =  iabs(subshell(no_one)%kappa) * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa) * 2 - 1
      i_max = min(j_1,j_2)
      if (i_max +1 > 20) stop "anco_case_6(): program stop A."
      do rank = 0, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_2_shells                &
                     (csf_r,csf_s,2*rank,no_form_one,no_form_two,.true.,delta_J)
         if (delta_J == 0) then
            temp%recoupling(rank + 1) = zero;  temp%operator(rank + 1) = zero
         else
            temp%operator(rank + 1)=                                           &
                               rcfp_calculate_Wk_me(bra_one,ket_one,rank,1,1) *&
                               rcfp_calculate_Wk_me(bra_two,ket_two,rank,-1,-1)
            if (abs(temp%operator(rank + 1)) < eps10) then
               temp%recoupling(rank + 1) = zero
            else
               temp%recoupling(rank + 1) = recoupling_matrix_2_shells          &
                 (csf_r,csf_s,2*rank,no_form_one,no_form_two,.false.,delta_J)
            end if
         end if
      end do
      i_min_2 = iabs(j_1 - j_2)/2;            i_max_2 = (j_1 + j_2)/2
      do rank = i_min_2, i_max_2
         coeff = zero
         do i_sum = 0, i_max
            coeff_one = temp%recoupling(i_sum+1) * temp%operator(i_sum+1)
            if (abs(coeff_one) > eps10) then
               if (wigner_6j_triangle(j_1,j_2,2*rank,j_2,j_1,2*i_sum) /= 0) then
                  coeff_one = coeff_one * sqrt(two * i_sum + one)             *&
                         wigner_6j_symbol(j_1,j_2,2*rank,j_2,j_1,2*i_sum,.true.)
                  if (mod(j_1+j_2+2*rank+2*i_sum,4) /= 0) coeff_one = -coeff_one
                  coeff = coeff + coeff_one
               end if
            end if
         end do
         coeff = -half * coeff
         if (abs(coeff) > eps10) then
            nocoeff                      = nocoeff + 1
            anco_V_list(nocoeff)%nu      = rank
            anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
            anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
            anco_V_list(nocoeff)%b%kappa = subshell(no_one)%kappa
            anco_V_list(nocoeff)%b%n     = subshell(no_one)%n
            anco_V_list(nocoeff)%c%kappa = subshell(no_two)%kappa
            anco_V_list(nocoeff)%c%n     = subshell(no_two)%n
            anco_V_list(nocoeff)%d%kappa = subshell(no_two)%kappa
            anco_V_list(nocoeff)%d%n     = subshell(no_two)%n
            anco_V_list(nocoeff)%V       = coeff
         end if
      end do
      !
   end subroutine anco_case_6
   !
   !
   subroutine anco_case_7_to_14(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                                no_one,no_two,no_three,no_four)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 11 - 14 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   !  7.         Beta    Alpha   Alpha   Alpha
   !  8.         Alpha   Beta    Alpha   Alhpa
   !  9.         Beta    Beta    Beta    Alpha
   ! 10.         Beta    Beta    Alpha   Beta
   ! 11.         Beta    Gamma   Alpha   Gamma
   ! 12.         Gamma   Beta    Gamma   Alpha
   ! 13.         Gamma   Beta    Alhpa   Gamma
   ! 14.         Beta    Gamma   Gamma   Alpha
   !
   ! Calls:  anco_case_7_to_8, anco_case_9_to_10,  anco_case_11_to_14
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_three, no_four
      integer, intent(inout)                 :: nocoeff
      !
      if (nwshells < 2) return
      if (no_two == no_four) then
         if ((no_one == no_two) .or. (no_three == no_two)) then
            if (no_one == no_three) then
               stop "anco_case_7_to_14(): program stop A."
            else if (no_three == no_two)then
               call anco_case_7_to_8(csf_r,csf_s,nwshells,nwcore,subshell,    &
                        nocoeff,no_three,no_one,no_one,no_two,no_three,no_four)
            else
               call anco_case_9_to_10(csf_r,csf_s,nwshells,nwcore,subshell,   &
                        nocoeff,no_three,no_one,no_one,no_two,no_three,no_four)
            end if
         else
            call anco_case_11_to_14(csf_r,csf_s,nwshells,nwcore,subshell,     &
                 nocoeff,no_three,no_one,no_two,no_one,no_two,no_three,no_four)
         end if
      else if (no_one == no_three) then
         if ((no_two == no_one) .or. (no_four == no_one)) then
            if (no_two == no_four) then
               stop "anco_case_7_to_14(): program stop B."
            else if (no_four == no_one)then
               call anco_case_7_to_8(csf_r,csf_s,nwshells,nwcore,subshell,    &
                         nocoeff,no_four,no_two,no_one,no_two,no_three,no_four)
            else
               call anco_case_9_to_10(csf_r,csf_s,nwshells,nwcore,subshell,   &
                         nocoeff,no_four,no_two,no_one,no_two,no_three,no_four)
            end if
         else
            call anco_case_11_to_14(csf_r,csf_s,nwshells,nwcore,subshell,     &
                  nocoeff,no_four,no_two,no_one,no_one,no_two,no_three,no_four)
         end if      
      else if (no_one == no_four) then
         if ((no_two == no_one) .or. (no_three == no_one)) then
            if (no_two == no_three) then
               stop "anco_case_7_to_14(): program stop C."
            else if (no_three == no_four)then
               call anco_case_7_to_8(csf_r,csf_s,nwshells,nwcore,subshell,    &
                        nocoeff,no_three,no_two,no_one,no_two,no_three,no_four)
            else
               call anco_case_9_to_10(csf_r,csf_s,nwshells,nwcore,subshell,   &
                        nocoeff,no_three,no_two,no_one,no_two,no_three,no_four)
            end if
         else
            call anco_case_11_to_14(csf_r,csf_s,nwshells,nwcore,subshell,     &
                 nocoeff,no_three,no_two,no_one,no_one,no_two,no_four,no_three)
         end if
      else if (no_two == no_three) then
         if ((no_one == no_two) .or. (no_four == no_two)) then
            if (no_one == no_four) then
               stop "anco_case_7_to_14(): program stop D."
            else if (no_four == no_two)then
               call anco_case_7_to_8(csf_r,csf_s,nwshells,nwcore,subshell,    &
                         nocoeff,no_four,no_one,no_one,no_two,no_three,no_four)
            else
               call anco_case_9_to_10(csf_r,csf_s,nwshells,nwcore,subshell,   &
                         nocoeff,no_four,no_one,no_one,no_two,no_three,no_four)
            end if
         else
            call anco_case_11_to_14(csf_r,csf_s,nwshells,nwcore,subshell,     &
                  nocoeff,no_four,no_one,no_two,no_one,no_two,no_four,no_three)
         end if
      else
         stop "anco_case_7_to_14(): program stop E."
      end if
      !
   end subroutine anco_case_7_to_14
   !
   !
   subroutine anco_case_7_to_8(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff,  &
                          no_one,no_two,no_one_w,no_two_w,no_three_w,no_four_w)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 7 - 8 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 7.         Beta    Alpha   Alpha   Alpha
   ! 8.         Alpha   Beta    Alpha   Alhpa
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_2_shells,
   !         rcfp_calculate_a_times_Wk_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_one_w, no_two_w, no_three_w
      integer, intent(in)                    :: no_four_w
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(anco_temp)      :: temp
      real(kind=dp)        :: operator_a, coeff, coeff_one
      integer              :: j_1, j_2, l_1, l_2, rank, delta_J
      integer              :: i, no_form_one, no_form_two, occup
      integer              :: i_min, i_max, i_sum
      !
      if (nwshells == 1) return
      no_form_one = min(no_one,no_two);       no_form_two = max(no_one,no_two)
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_form_one,no_form_two,   &
                   no_form_two,no_form_two,nwshells,nwcore,1)) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa) * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa) * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa / iabs(subshell(no_one)%kappa))/2
      l_2 = (j_2 + subshell(no_two)%kappa / iabs(subshell(no_two)%kappa))/2
      temp%recoupling(1) = recoupling_matrix_2_shells(csf_r, csf_s,       &
                                j_2,no_form_one,no_form_two,.true.,delta_J)
      if (delta_J == 0) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      operator_a = rcfp_calculate_a_me(bra_two,ket_two,1)
      if (abs(operator_a) < eps10) return
      temp%recoupling(1) = recoupling_matrix_2_shells(csf_r,csf_s,        &
                               j_2,no_form_one,no_form_two,.false.,delta_J)
      i_min = iabs(j_1 - j_2)/2;     i_max = min(j_1 + j_1,j_1 + j_2)/2
      if (i_max +1 > 20) stop "anco_case_7_to_8(): program stop A."
      do i_sum = i_min, i_max
         temp%operator(i_sum + 1) = rcfp_calculate_a_times_Wk_me(bra_one, &
                                                 ket_one,i_sum,j_2,1,-1,-1)
      end do
      !
      do rank = i_min, i_max
         coeff = zero
         do i_sum = i_min, i_max
            if (mod(j_2 - i_sum + 1,2) == 0) then
               coeff_one = temp%operator(i_sum+1)
               if (abs(coeff_one) > eps10) then
                  if (wigner_6j_triangle(j_2,j_1,2*rank,j_1,j_1,2*i_sum)/=0)then
                     coeff_one = coeff_one*sqrt(two*i_sum+one)*wigner_6j_symbol&
                                         (j_2,j_1,2*rank,j_1,j_1,2*i_sum,.true.)
                     if (mod(j_1+rank+i_sum,2)/= 0) coeff_one=-coeff_one
                     coeff = coeff + coeff_one
                  end if
               end if
            end if
         end do
         coeff = -coeff * temp%recoupling(1) * operator_a
         if (abs(coeff) > eps10) then
            occup = 0
            do i = no_form_one, no_form_two-1
               occup = occup + csf_r%occupation(i)
            end do
            if (mod(occup,2) == 0) coeff = -coeff
            nocoeff                      = nocoeff + 1
            anco_V_list(nocoeff)%nu      = rank
            anco_V_list(nocoeff)%a%kappa = subshell(no_one_w)%kappa
            anco_V_list(nocoeff)%a%n     = subshell(no_one_w)%n
            anco_V_list(nocoeff)%b%kappa = subshell(no_two_w)%kappa
            anco_V_list(nocoeff)%b%n     = subshell(no_two_w)%n
            anco_V_list(nocoeff)%c%kappa = subshell(no_three_w)%kappa
            anco_V_list(nocoeff)%c%n     = subshell(no_three_w)%n
            anco_V_list(nocoeff)%d%kappa = subshell(no_four_w)%kappa
            anco_V_list(nocoeff)%d%n     = subshell(no_four_w)%n
            anco_V_list(nocoeff)%V       = coeff
         end if
      end do
      !
   end subroutine anco_case_7_to_8
   !
   !
   subroutine anco_case_9_to_10(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                          no_one,no_two,no_one_w,no_two_w,no_three_w,no_four_w)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 9 - 10 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   !  9.         Beta    Beta    Beta    Alpha
   ! 10.         Beta    Beta    Alpha   Beta
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_2_shells,
   !         rcfp_calculate_Wk_times_a_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_one_w, no_two_w, no_three_w
      integer, intent(in)                    :: no_four_w
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(anco_temp)      :: temp
      real(kind=dp)        :: operator_a, coeff, coeff_one
      integer              :: j_1, j_2, l_1, l_2, rank, delta_J
      integer              :: i, no_form_one, no_form_two, occup
      integer              :: i_min, i_max, i_sum
      !
      if (nwshells == 1) return
      no_form_one = min(no_one,no_two);       no_form_two = max(no_one,no_two)
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_form_one,no_form_two,   &
                   no_form_two,no_form_two,nwshells,nwcore,1)) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa) * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa) * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa / iabs(subshell(no_one)%kappa))/2
      l_2 = (j_2 + subshell(no_two)%kappa / iabs(subshell(no_two)%kappa))/2
      temp%recoupling(1) = recoupling_matrix_2_shells(csf_r,csf_s,        &
                                j_1,no_form_one,no_form_two,.true.,delta_J)
      if (delta_J == 0) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      operator_a = rcfp_calculate_a_me(bra_one,ket_one,-1)
      if (abs(operator_a) < eps10) return
      temp%recoupling(1) = recoupling_matrix_2_shells(csf_r,csf_s,        &
                               j_1,no_form_one,no_form_two,.false.,delta_J)
      i_min = iabs(j_1 - j_2)/2;     i_max = min(j_2 + j_2,j_1 + j_2)/2
      if (i_max +1 > 20) stop "anco_case_9_to_10(): program stop A."
      do i_sum = i_min, i_max
         temp%operator(i_sum + 1) = rcfp_calculate_Wk_times_a_me(bra_two, &
                                                  ket_two,i_sum,j_1,1,1,-1)
      end do
      !
      do rank = i_min, i_max
         coeff = zero
         do i_sum = i_min, i_max
            if (mod(j_1 - i_sum + 1,2) == 0) then
               coeff_one = temp%operator(i_sum+1)
               if (abs(coeff_one) > eps10) then
                  if (wigner_6j_triangle(j_2,j_1,2*rank,j_2,j_2,2*i_sum)/=0)then
                     coeff_one = coeff_one*sqrt(two*i_sum+one)*wigner_6j_symbol&
                                         (j_2,j_1,2*rank,j_2,j_2,2*i_sum,.true.)
                     if (mod(j_2+rank+i_sum,2)/= 0) coeff_one=-coeff_one
                     coeff = coeff + coeff_one
                  end if
               end if
            end if
         end do
         coeff = -coeff * temp%recoupling(1) * operator_a
         if (abs(coeff) > eps10) then
            occup = 0
            do i = no_form_one, no_form_two-1
               occup = occup + csf_r%occupation(i)
            end do
            if (mod(occup,2) == 0) coeff = -coeff
            nocoeff                      = nocoeff + 1
            anco_V_list(nocoeff)%nu      = rank
            anco_V_list(nocoeff)%a%kappa = subshell(no_one_w)%kappa
            anco_V_list(nocoeff)%a%n     = subshell(no_one_w)%n
            anco_V_list(nocoeff)%b%kappa = subshell(no_two_w)%kappa
            anco_V_list(nocoeff)%b%n     = subshell(no_two_w)%n
            anco_V_list(nocoeff)%c%kappa = subshell(no_three_w)%kappa
            anco_V_list(nocoeff)%c%n     = subshell(no_three_w)%n
            anco_V_list(nocoeff)%d%kappa = subshell(no_four_w)%kappa
            anco_V_list(nocoeff)%d%n     = subshell(no_four_w)%n
            anco_V_list(nocoeff)%V       = coeff
         end if
      end do
      !
   end subroutine anco_case_9_to_10
   !
   !
   subroutine anco_case_11_to_14(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                  no_one,no_two,no_three,no_one_w,no_two_w,no_three_w,no_four_w)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 11 - 14 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 11.         Beta    Gamma   Alpha   Gamma
   ! 12.         Gamma   Beta    Gamma   Alpha
   ! 13.         Gamma   Beta    Alhpa   Gamma
   ! 14.         Beta    Gamma   Gamma   Alpha
   !
   ! Calls:  anco_codeing, anco_normal_form, anco_normal_phase,
   !         recoupling_matrix_3_shells, rcfp_calculate_a_me, 
   !         rcfp_calculate_Wk_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_one_w, no_two_w, no_three_w
      integer, intent(in)                    :: no_four_w, no_three
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(subshell_state) :: bra_three, ket_three
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one, operator_a_a
      integer              :: j_1, j_2, j_3, l_1, l_2, l_3, rank, delta_J
      integer              :: i, no_form_one, no_form_two, no_form_three
      integer              :: phase, occup, i_min, i_max, i_sum, i_min1, i_max1
      !
      if (nwshells < 3) return
      call anco_normal_form(no_one,no_two,no_three,no_form_one,no_form_two, &
                                                               no_form_three)
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_form_one,no_form_three,&
                  no_form_two,no_form_two,nwshells,nwcore,2)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      call anco_codeing(no_three,subshell,csf_r,csf_s,bra_three,ket_three)
      operator_a_a = rcfp_calculate_a_me(bra_one,ket_one,-1) *              &
                     rcfp_calculate_a_me(bra_two,ket_two, 1)
      if (abs(operator_a_a) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa)   * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa)   * 2 - 1
      j_3 =  iabs(subshell(no_three)%kappa) * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa   / iabs(subshell(no_one)%kappa))  /2
      l_2 = (j_2 + subshell(no_two)%kappa   / iabs(subshell(no_two)%kappa))  /2
      l_3 = (j_3 + subshell(no_three)%kappa / iabs(subshell(no_three)%kappa))/2
      i_min = iabs(j_1 - j_2)/2;     i_max = min(j_1 + j_2,j_3 + j_3)/2
      if (i_max+1 > 20) stop "anco_case_11_to_14(): program stop A."
      do rank = i_min, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_3_shells(csf_r,csf_s,  &
                          j_2,j_1,2*rank,no_two,no_one,no_three,.true.,delta_J)
         if (delta_J == 0) then
            temp%recoupling(rank + 1) = zero;  temp%operator(rank + 1) = zero
         else
            temp%operator(rank+1)= rcfp_calculate_Wk_me(bra_three,ket_three,  &
                                                                     rank,1,-1)
            if (abs(temp%operator(rank + 1)) < eps10) then
               temp%recoupling(rank + 1) = zero
            else
               temp%recoupling(rank + 1) = recoupling_matrix_3_shells(csf_r,  &
                   csf_s,j_2,j_1,2*rank,no_two,no_one,no_three,.false.,delta_J)
            end if
         end if
      end do
      !
      phase = anco_normal_phase(no_two,no_one,no_three,no_three)
      occup = 0
      do i = min(no_one,no_two), max(no_one,no_two)-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) phase = -phase
      !
      ! cases 2313  + + - -     transform to 2133  + - + -
      !       3231                           2133
      do rank = i_min, i_max
         coeff = zero
         if (mod(2*rank,2) == 0) then
            coeff = operator_a_a * phase * temp%recoupling(rank+1) *          &
                    temp%operator(rank+1) / sqrt(two * rank + one)
         end if
         if (abs(coeff) > eps10) then
            nocoeff                      = nocoeff + 1
            anco_V_list(nocoeff)%nu      = rank
            anco_V_list(nocoeff)%a%kappa = subshell(no_one_w)%kappa
            anco_V_list(nocoeff)%a%n     = subshell(no_one_w)%n
            anco_V_list(nocoeff)%b%kappa = subshell(no_two_w)%kappa
            anco_V_list(nocoeff)%b%n     = subshell(no_two_w)%n
            anco_V_list(nocoeff)%c%kappa = subshell(no_three_w)%kappa
            anco_V_list(nocoeff)%c%n     = subshell(no_three_w)%n
            anco_V_list(nocoeff)%d%kappa = subshell(no_four_w)%kappa
            anco_V_list(nocoeff)%d%n     = subshell(no_four_w)%n
            anco_V_list(nocoeff)%V       = coeff
         end if
      end do
      !
      ! cases 3213  + + - -     transform to 2133  + - + -
      !       2331                           2133
      i_min1 = max(iabs(j_1 - j_3),iabs(j_2 - j_3))/2
      i_max1 = min(j_1 + j_3,j_2 + j_3)/2
      if (i_min1 <= i_max1) then
         do rank = i_min1, i_max1
            coeff = zero
            do i_sum = i_min, i_max
               coeff_one = operator_a_a * phase * temp%recoupling(i_sum+1) *  &
                                                         temp%operator(i_sum+1)
               if (abs(coeff_one) > eps10) then
                  if (wigner_6j_triangle(j_1,j_3,2*rank,j_3,j_2,2*i_sum)/=0)then
                     coeff_one = coeff_one * sqrt(two * i_sum + one) *     &
                     wigner_6j_symbol(j_1,j_3,2*rank,j_3,j_2,2*i_sum,.true.)
                     coeff = coeff + coeff_one
                  end if
               end if
            end do
            if (abs(coeff) > eps10) then
               nocoeff                      = nocoeff + 1
               anco_V_list(nocoeff)%nu      = rank
               anco_V_list(nocoeff)%a%kappa = subshell(no_one_w)%kappa
               anco_V_list(nocoeff)%a%n     = subshell(no_one_w)%n
               anco_V_list(nocoeff)%b%kappa = subshell(no_two_w)%kappa
               anco_V_list(nocoeff)%b%n     = subshell(no_two_w)%n
               anco_V_list(nocoeff)%c%kappa = subshell(no_four_w)%kappa
               anco_V_list(nocoeff)%c%n     = subshell(no_four_w)%n
               anco_V_list(nocoeff)%d%kappa = subshell(no_three_w)%kappa
               anco_V_list(nocoeff)%d%n     = subshell(no_three_w)%n
               anco_V_list(nocoeff)%V       = coeff
            end if
         end do
      end if
      !
   end subroutine anco_case_11_to_14
   !
   !
   subroutine anco_case_15_to_18(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                                 no_one,no_two,no_three,no_four)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 15 - 18 (of Table 1 in G. Gaigalas et al., 1997 
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 15.         Gamma   Gamma   Alhpa   Beta
   ! 16.         Gamma   Gamma   Beta    Alpha
   ! 17.         Alpha   Beta    Gamma   Gamma
   ! 18.         Beta    Alpha   Gamma   Gamma
   !
   ! Calls:  anco_case_15_to_18_order
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_three, no_four
      integer, intent(inout)                 :: nocoeff
      !
      if (nwshells < 3) return
      if (no_one == no_two) then
         call anco_case_15_to_18_order(csf_r,csf_s,nwshells,nwcore,subshell,  &
              nocoeff,no_three,no_four,no_one,no_one,no_two,no_three,no_four,1)
      else if (no_three == no_four) then
         call anco_case_15_to_18_order(csf_r,csf_s,nwshells,nwcore,subshell,  &
               nocoeff,no_one,no_two,no_three,no_one,no_two,no_three,no_four,2)
      else
         stop "anco_case_15_to_18(): program stop A."
      end if
      !
   end subroutine anco_case_15_to_18
   !
   !
   subroutine anco_case_15_to_18_order(csf_r,csf_s,nwshells,nwcore,subshell,  &
                   nocoeff,no_one,no_two,no_three,no_one_w,no_two_w,          &
                   no_three_w,no_four_w,irez)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 15 - 18 (of Table 1 in G. Gaigalas et al., 1997 
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 15.         Gamma   Gamma   Alhpa   Beta
   ! 16.         Gamma   Gamma   Beta    Alpha
   ! 17.         Alpha   Beta    Gamma   Gamma
   ! 18.         Beta    Alpha   Gamma   Gamma
   !
   ! Calls:  anco_codeing, anco_normal_form, anco_normal_phase,
   !         recoupling_matrix_check, recoupling_matrix_3_shells, 
   !         rcfp_calculate_a_me, rcfp_calculate_Wk_me
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_one_w, no_two_w, no_three_w
      integer, intent(in)                    :: no_four_w, no_three, irez
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(subshell_state) :: bra_three, ket_three
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one, operator_a_a
      integer              :: j_1, j_2, j_3, l_1, l_2, l_3, rank, delta_J
      integer              :: i, no_form_one, no_form_two, no_form_three
      integer              :: phase, occup, i_min, i_max, i_sum, i_min1, i_max1
      integer              :: qm1, qm2, qm3, qm4, phase_b
      !
      if (nwshells < 3) return
      call anco_normal_form(no_one,no_two,no_three,no_form_one,no_form_two,   &
                                                                 no_form_three)
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_form_one,no_form_three,  &
                    no_form_two,no_form_two,nwshells,nwcore,2)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      call anco_codeing(no_three,subshell,csf_r,csf_s,bra_three,ket_three)
      !
      if (irez == 1) then
         qm1 = -1; qm2 = -1; qm3 =  1; qm4 =  1
      else if( irez == 2) then
         qm1 =  1; qm2 =  1; qm3 = -1; qm4 = -1
      else
         stop "anco_case_15_to_18_order(): program stop A."
      end if
      !
      operator_a_a = rcfp_calculate_a_me(bra_one,ket_one,qm1) *               &
                     rcfp_calculate_a_me(bra_two,ket_two,qm2)
      if (abs(operator_a_a) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa)   * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa)   * 2 - 1
      j_3 =  iabs(subshell(no_three)%kappa) * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa   / iabs(subshell(no_one)%kappa))  /2
      l_2 = (j_2 + subshell(no_two)%kappa   / iabs(subshell(no_two)%kappa))  /2
      l_3 = (j_3 + subshell(no_three)%kappa / iabs(subshell(no_three)%kappa))/2
      i_min = iabs(j_1 - j_2)/2;     i_max = min(j_1 + j_2,j_3 + j_3)/2
      if (i_max+1 > 20) stop "anco_case_15_to_18_order(): program stop B."
      do rank = i_min, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_3_shells(csf_r,csf_s,  &
                          j_1,j_2,2*rank,no_one,no_two,no_three,.true.,delta_J)
         if (delta_J == 0) then
            temp%recoupling(rank + 1) = zero;  temp%operator(rank + 1) = zero
         else
            temp%operator(rank+1) = rcfp_calculate_Wk_me(bra_three,ket_three, &
                                                                  rank,qm3,qm4)
            if (abs(temp%operator(rank + 1)) < eps10) then
               temp%recoupling(rank + 1) = zero
            else
               temp%recoupling(rank + 1) = recoupling_matrix_3_shells(csf_r,  &
                   csf_s,j_1,j_2,2*rank,no_one,no_two,no_three,.false.,delta_J)
            end if
         end if
      end do
      !
      phase = anco_normal_phase(no_three,no_one,no_two,no_three)
      occup = 0
      do i = min(no_one,no_two), max(no_one,no_two)-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) phase = -phase
      !
      ! cases 3312  + + - -     transform to 1233  - - + +
      !       3321                           1233
      !                                                   (irez = 1)
      !
      ! cases 1233  + + - -     transform to 1233  + + - -
      !       2133                           1233
      !                                                   (irez = 2)
      i_min1 = max(iabs(j_1 - j_3),iabs(j_2 - j_3))/2
      i_max1 = min(j_1 + j_3,j_2 + j_3)/2
      if (i_min1 <= i_max1) then
         do rank = i_min1, i_max1
            coeff = zero
            do i_sum = i_min, i_max
               if(irez == 1) then
                  phase_b = j_2 - i_sum + 1
               else if (irez == 2) then
                  phase_b = j_1 - i_sum + 1
               end if               
               if (mod(phase_b,2) == 0) then
                  coeff_one = operator_a_a * phase * temp%recoupling(i_sum+1)* &
                                                          temp%operator(i_sum+1)
                  if (abs(coeff_one) > eps10) then
                     if (wigner_6j_triangle(j_3,j_1,2*rank,j_2,j_3,2*i_sum)/=0)&
                                                                            then
                        coeff_one = -coeff_one * sqrt(two * i_sum + one)     * &
                        wigner_6j_symbol(j_3,j_1,2*rank,j_2,j_3,2*i_sum,.true.)
                        if (irez == 1) then
                           if (mod(j_1 + j_3 + 2*i_sum + 2*rank,4) /= 0)       &
                                                            coeff_one=-coeff_one
                        else if (irez == 2) then
                           if (mod(j_2 + j_3 + 2*i_sum + 2*rank,4) /= 0)       &
                                                          coeff_one = -coeff_one
                        end if
                        coeff = coeff + coeff_one
                     end if
                  end if
               end if
            end do
            if (abs(coeff) > eps10) then
               nocoeff                      = nocoeff + 1
               anco_V_list(nocoeff)%nu      = rank
               anco_V_list(nocoeff)%a%kappa = subshell(no_one_w)%kappa
               anco_V_list(nocoeff)%a%n     = subshell(no_one_w)%n
               anco_V_list(nocoeff)%b%kappa = subshell(no_two_w)%kappa
               anco_V_list(nocoeff)%b%n     = subshell(no_two_w)%n
               anco_V_list(nocoeff)%c%kappa = subshell(no_three_w)%kappa
               anco_V_list(nocoeff)%c%n     = subshell(no_three_w)%n
               anco_V_list(nocoeff)%d%kappa = subshell(no_four_w)%kappa
               anco_V_list(nocoeff)%d%n     = subshell(no_four_w)%n
               anco_V_list(nocoeff)%V       = coeff
            end if
         end do
      end if
      !
   end subroutine anco_case_15_to_18_order
   !
   !
   subroutine anco_case_19_to_42(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                                 no_one,no_two,no_three,no_four)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 19 - 26 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 19.         Alpha   Beta    Gamma   Delta
   ! 20.         Beta    Alpha   Gamma   Delta
   ! 21.         Alpha   Beta    Delta   Gamma
   ! 22.         Beta    Alpha   Delta   Gamma
   ! 23.         Gamma   Delta   Alpha   Beta
   ! 24.         Gamma   Delta   Beta    Alpha
   ! 25.         Delta   Gamma   Alpha   Beta
   ! 26.         Delta   Gamma   Beta    Alpha
   ! 27.         Alpha   Gamma   Beta    Delta
   ! 28.         Alpha   Gamma   Delta   Beta
   ! 29.         Gamma   Alpha   Delta   Beta
   ! 30.         Gamma   Alpha   Beta    Delta
   ! 31.         Beta    Delta   Alpha   Gamma
   ! 32.         Delta   Beta    Gamma   Alpha
   ! 33.         Beta    Delta   Gamma   Alpha
   ! 34.         Delta   Beta    Alpha   Gamma
   ! 35.         Alhpa   Delta   Beta    Gamma
   ! 36.         Delta   Alhpa   Gamma   Beta
   ! 37.         Alpha   Delta   Gamma   Beta
   ! 38.         Delta   Alpha   Beta    Gamma
   ! 39.         Beta    Gamma   Alpha   Delta
   ! 40.         Gamma   Beta    Delta   Alhpa
   ! 41.         Beta    Gamma   Delta   Alhpa
   ! 42.         Gamma   Beta    Alpha   Beta
   !
   ! Calls:  anco_case_19_to_26, anco_case_27_to_34, anco_case_35_to_42
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_three, no_four
      integer, intent(inout)                 :: nocoeff
      !
      if (nwshells < 4) return
      if (no_two < no_three) then
         call anco_case_19_to_26(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                               no_one,no_two,no_three,no_four,1)
      else if ((no_one > no_four) .and. (no_two > no_four)) then
         call anco_case_19_to_26(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                               no_three,no_four,no_one,no_two,2)
      else if ((no_two > no_three) .and. (no_two < no_four) .and.              &
              (no_one < no_three)) then
         call anco_case_27_to_34(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                               no_one,no_three,no_two,no_four,1)
      else if ((no_two > no_three) .and. (no_two > no_four) .and.              &
              (no_one > no_three)) then
         call anco_case_27_to_34(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                               no_three,no_one,no_four,no_two,2)
      else if ((no_two > no_three) .and. (no_two > no_four) .and.              &
              (no_one < no_three)) then
         call anco_case_35_to_42(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                               no_one,no_three,no_four,no_two,1)
      else if ((no_two > no_three) .and. (no_two < no_four) .and.              &
              (no_one > no_three)) then
         call anco_case_35_to_42(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                               no_three,no_one,no_two,no_four,2)
      else
         stop "anco_case_19_to_42_order(): program stop A."
      end if
      !
   end subroutine anco_case_19_to_42
   !
   !
   subroutine anco_case_19_to_26(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                            no_one,no_two,no_three,no_four,irez)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 19 - 26 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 19.         Alpha   Beta    Gamma   Delta
   ! 20.         Beta    Alpha   Gamma   Delta
   ! 21.         Alpha   Beta    Delta   Gamma
   ! 22.         Beta    Alpha   Delta   Gamma
   ! 23.         Gamma   Delta   Alpha   Beta
   ! 24.         Gamma   Delta   Beta    Alpha
   ! 25.         Delta   Gamma   Alpha   Beta
   ! 26.         Delta   Gamma   Beta    Alpha
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_4_shells,
   !         rcfp_calculate_a_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_three, no_four, irez
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(subshell_state) :: bra_three, ket_three, bra_four, ket_four 
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one
      integer              :: i, j_1, j_2, j_3, j_4, l_1, l_2, l_3, l_4
      integer              :: rank, delta_J
      integer              :: phase, occup, i_min, i_max, i_sum, i_min1, i_max1
      integer              :: qm1, qm2, qm3, qm4, phase_b
      !
      if (nwshells < 4) return
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_one,no_four,no_three,     &
                                      no_two,nwshells,nwcore,3)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      call anco_codeing(no_three,subshell,csf_r,csf_s,bra_three,ket_three)
      call anco_codeing(no_four,subshell,csf_r,csf_s,bra_four,ket_four)
      if (irez == 1) then
         qm1 =  1; qm2 =  1; qm3 = -1; qm4 = -1
      else if (irez == 2) then
         qm1 = -1; qm2 = -1; qm3 =  1; qm4 =  1
      else
         stop "anco_case_19_to_26_order(): program stop A."
      end if
      temp%operator(1) = rcfp_calculate_a_me(bra_one,ket_one,qm1)             *&
                         rcfp_calculate_a_me(bra_two,ket_two,qm2)             *&
                         rcfp_calculate_a_me(bra_three,ket_three,qm3)         *&
                         rcfp_calculate_a_me(bra_four,ket_four,qm4)
      if (abs(temp%operator(1)) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa)   * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa)   * 2 - 1
      j_3 =  iabs(subshell(no_three)%kappa) * 2 - 1
      j_4 =  iabs(subshell(no_four)%kappa)  * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa   / iabs(subshell(no_one)%kappa))  /2
      l_2 = (j_2 + subshell(no_two)%kappa   / iabs(subshell(no_two)%kappa))  /2
      l_3 = (j_3 + subshell(no_three)%kappa / iabs(subshell(no_three)%kappa))/2
      l_4 = (j_4 + subshell(no_four)%kappa  / iabs(subshell(no_four)%kappa)) /2
      i_min = max(iabs(j_1 - j_2),iabs(j_3 - j_4))/2
      i_max = min(j_1 + j_2,j_3 + j_4)/2
      if (i_min > i_max) return
      if (i_max+1 > 20) stop "anco_case_19_to_26_order(): program stop B."

      do rank = i_min, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_4_shells(csf_r,csf_s,   &
                                     j_1,j_2,j_3,j_4,2*rank,no_one,no_two,     &
                                     no_three,no_four,.true.,delta_J)
         if(delta_J == 0) then
            temp%recoupling(rank + 1) = zero
         else
            temp%recoupling(rank + 1) = recoupling_matrix_4_shells(csf_r,csf_s,&
                                        j_1,j_2,j_3,j_4,2*rank,no_one,no_two,  &
                                        no_three,no_four,.false.,delta_J)
         end if
      end do
      phase = 1
      occup = 0
      do i = no_one, no_two-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) phase = -phase
      occup = 0
      do i = no_three, no_four-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) phase = -phase
      !
      ! cases 1234  + + - -     transform to 1234  + + - -
      !       2134                           1234
      !                                                   (irez = 1)
      !
      ! cases 3412  + + - -     transform to 1234  - - + +
      !       3421                           1234
      !                                                   (irez = 2)
      i_min1 = max(iabs(j_1 - j_3),iabs(j_2 - j_4))/2
      i_max1 = min(j_1 + j_3,j_2 + j_4)/2
      if (i_min1 <= i_max1) then
         do rank = i_min1, i_max1
            coeff = zero
            do i_sum = i_min, i_max
               if (irez == 1) then
                  phase_b = j_1 + j_4 + 2*i_sum
               else if (irez == 2) then
                  phase_b = j_2 + j_3 - 2*i_sum
               end if               
               if (mod(phase_b,2) == 0) then
                  coeff_one = phase * temp%recoupling(i_sum+1)* temp%operator(1)
                  if (abs(coeff_one) > eps10) then
                     if (wigner_6j_triangle(j_1,j_3,2*rank,j_4,j_2,2*i_sum)/=0)&
                                                                            then
                        coeff_one = -coeff_one * sqrt(two * i_sum + one)     * &
                        wigner_6j_symbol(j_1,j_3,2*rank,j_4,j_2,2*i_sum,.true.)
                        if (irez == 1) then
                           if (mod(j_2 + j_3 + 2*i_sum + 2*rank,4) /= 0)       &
                                                            coeff_one=-coeff_one
                        else if (irez == 2) then
                           if (mod(j_1 + j_4 - 2*i_sum - 2*rank,4) /= 0)       &
                                                          coeff_one = -coeff_one
                        end if
                        coeff = coeff + coeff_one
                     end if
                  end if
               end if
            end do
            if (abs(coeff) > eps10) then
               if (irez == 1) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%V       = coeff
               else if(irez == 2) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%V       = coeff
               end if
            end if
         end do
      end if
      !
      ! cases 1243  + + - -     transform to 1234  + + - -
      !       2134                           1234
      !                                                   (irez = 1)
      !
      ! cases 3421  + + - -     transform to 1234  - - + +
      !       4321                           1234
      !                                                   (irez = 2)
      i_min1 = max(iabs(j_1 - j_4),iabs(j_2 - j_3))/2
      i_max1 = min(j_1 + j_4,j_2 + j_3)/2
      if (i_min1 <= i_max1) then
         do rank = i_min1, i_max1
            coeff = zero
            do i_sum = i_min, i_max
               if (irez == 1) then
                  phase_b = j_1 - j_4
               else if (irez == 2) then
                  phase_b = j_3 - j_2
               end if               
               if (mod(phase_b,2) == 0) then
                  coeff_one = phase * temp%recoupling(i_sum+1)* temp%operator(1)
                  if (abs(coeff_one) > eps10) then
                     if (wigner_6j_triangle(j_1,j_4,2*rank,j_3,j_2,2*i_sum)/=0)&
                                                                            then
                        coeff_one = coeff_one * sqrt(two * i_sum + one)      * &
                        wigner_6j_symbol(j_1,j_4,2*rank,j_3,j_2,2*i_sum,.true.)
                        if (irez == 1) then
                           if (mod(j_2 + j_3 + 2*j_4 - 2*rank,4) /= 0)         &
                                                            coeff_one=-coeff_one
                        else if (irez == 2) then
                        if (mod(j_1+ 2*j_2 + j_4 + 4*i_sum + 2*rank,4) /= 0)&
                                                          coeff_one = -coeff_one
                        end if
                        coeff = coeff + coeff_one
                     end if
                  end if
               end if
            end do
            if (abs(coeff) > eps10) then
               if (irez == 1) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%V       = coeff
               else if(irez == 2) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%V       = coeff
               end if
            end if
         end do
      end if
      !
   end subroutine anco_case_19_to_26
   !
   !
   subroutine anco_case_27_to_34(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                            no_one,no_two,no_three,no_four,irez)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 27 - 34 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 27.         Alpha   Gamma   Beta    Delta
   ! 28.         Alpha   Gamma   Delta   Beta
   ! 29.         Gamma   Alpha   Delta   Beta
   ! 30.         Gamma   Alpha   Beta    Delta
   ! 31.         Beta    Delta   Alpha   Gamma
   ! 32.         Delta   Beta    Gamma   Alpha
   ! 33.         Beta    Delta   Gamma   Alpha
   ! 34.         Delta   Beta    Alpha   Gamma
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_4_shells,
   !         rcfp_calculate_a_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_three, no_four, irez
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(subshell_state) :: bra_three, ket_three, bra_four, ket_four 
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one
      integer              :: i, j_1, j_2, j_3, j_4, l_1, l_2, l_3, l_4
      integer              :: rank, delta_J
      integer              :: phase, occup, i_min, i_max, i_sum, i_min1, i_max1
      integer              :: qm1, qm2, qm3, qm4, phase_b
      !
      if (nwshells < 4) return
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_one,no_four,no_three,     &
                                      no_two,nwshells,nwcore,3)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      call anco_codeing(no_three,subshell,csf_r,csf_s,bra_three,ket_three)
      call anco_codeing(no_four,subshell,csf_r,csf_s,bra_four,ket_four)
      if (irez == 1) then
         qm1 =  1; qm2 = -1; qm3 =  1; qm4 = -1
      else if (irez == 2) then
         qm1 = -1; qm2 =  1; qm3 = -1; qm4 =  1
      else
         stop "anco_case_27_to_34_order(): program stop A."
      end if
      temp%operator(1) = rcfp_calculate_a_me(bra_one,ket_one,qm1)             *&
                         rcfp_calculate_a_me(bra_two,ket_two,qm2)             *&
                         rcfp_calculate_a_me(bra_three,ket_three,qm3)         *&
                         rcfp_calculate_a_me(bra_four,ket_four,qm4)
      if (abs(temp%operator(1)) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa)   * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa)   * 2 - 1
      j_3 =  iabs(subshell(no_three)%kappa) * 2 - 1
      j_4 =  iabs(subshell(no_four)%kappa)  * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa   / iabs(subshell(no_one)%kappa))  /2
      l_2 = (j_2 + subshell(no_two)%kappa   / iabs(subshell(no_two)%kappa))  /2
      l_3 = (j_3 + subshell(no_three)%kappa / iabs(subshell(no_three)%kappa))/2
      l_4 = (j_4 + subshell(no_four)%kappa  / iabs(subshell(no_four)%kappa)) /2
      i_min = max(iabs(j_1 - j_2),iabs(j_3 - j_4))/2
      i_max = min(j_1 + j_2,j_3 + j_4)/2
      if (i_min > i_max) return
      if (i_max+1 > 20) stop "anco_case_27_to_34_order(): program stop B."

      do rank = i_min, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_4_shells(csf_r,csf_s,   &
                                     j_1,j_2,j_3,j_4,2*rank,no_one,no_two,     &
                                     no_three,no_four,.true.,delta_J)
         if (delta_J == 0) then
            temp%recoupling(rank + 1) = zero
         else
            temp%recoupling(rank + 1) = recoupling_matrix_4_shells(csf_r,csf_s,&
                                        j_1,j_2,j_3,j_4,2*rank,no_one,no_two,  &
                                        no_three,no_four,.false.,delta_J)
         end if
      end do
      phase = 1
      occup = 0
      do i = no_one, no_two-1
         occup = occup + csf_r%occupation(i)
      end do
         if (mod(occup,2) == 0) phase = -phase
      occup = 0
      do i = no_three, no_four-1
         occup = occup + csf_r%occupation(i)
      end do
         if (mod(occup,2) == 0) phase = -phase
      !
      ! cases 1324  + + - -     transform to 1234  - + - +
      !       1342                           1234
      !                                                   (irez = 1)
      !
      ! cases 2413  + + - -     transform to 1234  + - + -
      !       4231                           1234
      !                                                   (irez = 2)
      do rank = i_min, i_max
         coeff=phase*temp%recoupling(rank+1)*temp%operator(1)/sqrt(two*rank+one)
         if (irez == 2) then
            if (mod(j_1+j_2+j_3+j_4+4*rank,4)/=0) coeff = -coeff
         end if
         if (abs(coeff) > eps10) then
            if (irez == 1) then
               nocoeff                      = nocoeff + 1
               anco_V_list(nocoeff)%nu      = rank
               anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
               anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
               anco_V_list(nocoeff)%b%kappa = subshell(no_three)%kappa
               anco_V_list(nocoeff)%b%n     = subshell(no_three)%n
               anco_V_list(nocoeff)%c%kappa = subshell(no_two)%kappa
               anco_V_list(nocoeff)%c%n     = subshell(no_two)%n
               anco_V_list(nocoeff)%d%kappa = subshell(no_four)%kappa
               anco_V_list(nocoeff)%d%n     = subshell(no_four)%n
               anco_V_list(nocoeff)%V       = coeff
            else if(irez == 2) then
               nocoeff                      = nocoeff + 1
               anco_V_list(nocoeff)%nu      = rank
               anco_V_list(nocoeff)%a%kappa = subshell(no_two)%kappa
               anco_V_list(nocoeff)%a%n     = subshell(no_two)%n
               anco_V_list(nocoeff)%b%kappa = subshell(no_four)%kappa
               anco_V_list(nocoeff)%b%n     = subshell(no_four)%n
               anco_V_list(nocoeff)%c%kappa = subshell(no_one)%kappa
               anco_V_list(nocoeff)%c%n     = subshell(no_one)%n
               anco_V_list(nocoeff)%d%kappa = subshell(no_three)%kappa
               anco_V_list(nocoeff)%d%n     = subshell(no_three)%n
               anco_V_list(nocoeff)%V       = coeff
            end if
         end if
      end do
      !
      ! cases 1342  + + - -     transform to 1234  + - + -
      !       3124                           1234
      !                                                   (irez = 1)
      !
      ! cases 2431  + + - -     transform to 1234  - + - +
      !       4213                           1234
      !                                                   (irez = 2)
      i_min1 = max(iabs(j_1 - j_4),iabs(j_2 - j_3))/2
      i_max1 = min(j_1 + j_4,j_2 + j_3)/2
      if (i_min1 <= i_max1) then
         do rank = i_min1, i_max1
            coeff = zero
            do i_sum = i_min, i_max
               if (irez == 1) then
                  phase_b = j_1 - j_3 + 2*i_sum
               else if (irez == 2) then
                  phase_b = j_4 - j_2 - 2*i_sum
               end if               
               if (mod(phase_b,2) == 0) then
                  coeff_one = phase * temp%recoupling(i_sum+1)* temp%operator(1)
                  if (abs(coeff_one) > eps10) then
                     if (wigner_6j_triangle(j_1,j_4,2*rank,j_3,j_2,2*i_sum)/=0)&
                                                                            then
                        coeff_one = -coeff_one * sqrt(two * i_sum + one)     * &
                        wigner_6j_symbol(j_1,j_4,2*rank,j_3,j_2,2*i_sum,.true.)
                        if (irez == 1) then
                           if (mod(2*j_3 - 4*i_sum - 4*rank,4) /= 0)           &
                                                            coeff_one=-coeff_one
                        else if (irez == 2) then
                           if (mod(j_1 + j_2 + j_3 - j_4 - 4*rank,4) /= 0)     &
                                                          coeff_one = -coeff_one
                        end if
                        coeff = coeff + coeff_one
                     end if
                  end if
               end if
            end do
            if (abs(coeff) > eps10) then
               if (irez == 1) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%V       = coeff
              else if(irez == 2) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%V       = coeff
               end if
            end if
         end do
      end if
      !
   end subroutine anco_case_27_to_34
   !
   !
   subroutine anco_case_35_to_42(csf_r,csf_s,nwshells,nwcore,subshell,nocoeff, &
                                            no_one,no_two,no_three,no_four,irez)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for distributions 35 - 42 (of Table 1 in G. Gaigalas et al., 1997
   ! J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):
   !
   ! 35.         Alhpa   Delta   Beta    Gamma
   ! 36.         Delta   Alhpa   Gamma   Beta
   ! 37.         Alpha   Delta   Gamma   Beta
   ! 38.         Delta   Alpha   Beta    Gamma
   ! 39.         Beta    Gamma   Alpha   Delta
   ! 40.         Gamma   Beta    Delta   Alhpa
   ! 41.         Beta    Gamma   Delta   Alhpa
   ! 42.         Gamma   Beta    Alpha   Beta
   !
   ! Calls:  anco_codeing, recoupling_matrix_check, recoupling_matrix_4_shells,
   !         rcfp_calculate_a_me, wigner_6j_triangle, wigner_6j_symbol
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore, no_one, no_two
      integer, intent(in)                    :: no_three, no_four, irez
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_one, ket_one, bra_two, ket_two
      type(subshell_state) :: bra_three, ket_three, bra_four, ket_four 
      type(anco_temp)      :: temp
      real(kind=dp)        :: coeff, coeff_one
      integer              :: i, j_1, j_2, j_3, j_4, l_1, l_2, l_3, l_4
      integer              :: rank, delta_J
      integer              :: phase, occup, i_min, i_max, i_sum, i_min1, i_max1
      integer              :: qm1, qm2, qm3, qm4, phase_b
      !
      if (nwshells < 4) return
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_one,no_four,no_three,     &
                                      no_two,nwshells,nwcore,3)) < eps10) return
      call anco_codeing(no_one,subshell,csf_r,csf_s,bra_one,ket_one)
      call anco_codeing(no_two,subshell,csf_r,csf_s,bra_two,ket_two)
      call anco_codeing(no_three,subshell,csf_r,csf_s,bra_three,ket_three)
      call anco_codeing(no_four,subshell,csf_r,csf_s,bra_four,ket_four)
      if (irez == 1) then
         qm1 =  1; qm2 = -1; qm3 = -1; qm4 =  1
      else if (irez == 2) then
         qm1 = -1; qm2 =  1; qm3 =  1; qm4 = -1
      else
         stop "anco_case_35_to_42_order(): program stop A."
      end if
      temp%operator(1) = rcfp_calculate_a_me(bra_one,ket_one,qm1)             *&
                         rcfp_calculate_a_me(bra_two,ket_two,qm2)             *&
                         rcfp_calculate_a_me(bra_three,ket_three,qm3)         *&
                         rcfp_calculate_a_me(bra_four,ket_four,qm4)
      if (abs(temp%operator(1)) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa)   * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa)   * 2 - 1
      j_3 =  iabs(subshell(no_three)%kappa) * 2 - 1
      j_4 =  iabs(subshell(no_four)%kappa)  * 2 - 1
      l_1 = (j_1 + subshell(no_one)%kappa   / iabs(subshell(no_one)%kappa))  /2
      l_2 = (j_2 + subshell(no_two)%kappa   / iabs(subshell(no_two)%kappa))  /2
      l_3 = (j_3 + subshell(no_three)%kappa / iabs(subshell(no_three)%kappa))/2
      l_4 = (j_4 + subshell(no_four)%kappa  / iabs(subshell(no_four)%kappa)) /2
      i_min = max(iabs(j_1 - j_2),iabs(j_3 - j_4))/2
      i_max = min(j_1 + j_2,j_3 + j_4)/2
      if (i_min > i_max) return
      if (i_max+1 > 20) stop "anco_case_35_to_42_order(): program stop B."

      do rank = i_min, i_max
         temp%recoupling(rank + 1) = recoupling_matrix_4_shells(csf_r,csf_s,   &
                                     j_1,j_2,j_3,j_4,2*rank,no_one,no_two,     &
                                     no_three,no_four,.true.,delta_J)
         if (delta_J == 0) then
            temp%recoupling(rank + 1) = zero
         else
            temp%recoupling(rank + 1) = recoupling_matrix_4_shells(csf_r,csf_s,&
                                        j_1,j_2,j_3,j_4,2*rank,no_one,no_two,  &
                                        no_three,no_four,.false.,delta_J)
         end if
      end do
      phase = 1
      occup = 0
      do i = no_one, no_two-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) phase = -phase
      occup = 0
      do i = no_three, no_four-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) phase = -phase
      !
      ! cases 1423  + + - -     transform to 1234  + - - +
      !       4132                           1234
      !                                                   (irez = 1)
      !
      ! cases 2314  + + - -     transform to 1234  - + + -
      !       3241                           1234
      !                                                   (irez = 2)
      do rank = i_min, i_max
         coeff=phase*temp%recoupling(rank+1)*temp%operator(1)/sqrt(two*rank+one)
         if (irez == 1) then
            if (mod(j_3+j_4-2*rank+2,4)/=0) coeff = -coeff
         else if (irez == 2) then
            if (mod(j_1+j_2-2*rank+2,4)/=0) coeff = -coeff
         end if
         if (abs(coeff) > eps10) then
            if (irez == 1) then
               nocoeff                      = nocoeff + 1
               anco_V_list(nocoeff)%nu      = rank
               anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
               anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
               anco_V_list(nocoeff)%b%kappa = subshell(no_four)%kappa
               anco_V_list(nocoeff)%b%n     = subshell(no_four)%n
               anco_V_list(nocoeff)%c%kappa = subshell(no_two)%kappa
               anco_V_list(nocoeff)%c%n     = subshell(no_two)%n
               anco_V_list(nocoeff)%d%kappa = subshell(no_three)%kappa
               anco_V_list(nocoeff)%d%n     = subshell(no_three)%n
               anco_V_list(nocoeff)%V       = coeff
            else if(irez == 2) then
               nocoeff                      = nocoeff + 1
               anco_V_list(nocoeff)%nu      = rank
               anco_V_list(nocoeff)%a%kappa = subshell(no_two)%kappa
               anco_V_list(nocoeff)%a%n     = subshell(no_two)%n
               anco_V_list(nocoeff)%b%kappa = subshell(no_three)%kappa
               anco_V_list(nocoeff)%b%n     = subshell(no_three)%n
               anco_V_list(nocoeff)%c%kappa = subshell(no_one)%kappa
               anco_V_list(nocoeff)%c%n     = subshell(no_one)%n
               anco_V_list(nocoeff)%d%kappa = subshell(no_four)%kappa
               anco_V_list(nocoeff)%d%n     = subshell(no_four)%n
               anco_V_list(nocoeff)%V       = coeff
            end if
         end if
      end do
      !
      ! cases 1432  + + - -     transform to 1234  + - - +
      !       4132                           1234
      !                                                   (irez = 1)
      !
      ! cases 2341  + + - -     transform to 1234  - + + -
      !       3214                           1234
      !                                                   (irez = 2)
      i_min1 = max(iabs(j_1 - j_3),iabs(j_2 - j_4))/2
      i_max1 = min(j_1 + j_3,j_2 + j_4)/2
      if (i_min1 <= i_max1) then
         do rank = i_min1, i_max1
            coeff = zero
            do i_sum = i_min, i_max
               if (irez == 1) then
                  phase_b = j_1 + j_4 + 2*i_sum
               else if (irez == 2) then
                  phase_b = j_2 + j_3 + 2*i_sum
               end if               
               if (mod(phase_b,2) == 0) then
                  coeff_one = phase * temp%recoupling(i_sum+1)* temp%operator(1)
                  if (abs(coeff_one) > eps10) then
                     if (wigner_6j_triangle(j_1,j_3,2*rank,j_4,j_2,2*i_sum)/=0)&
                                                                            then
                        coeff_one = coeff_one * sqrt(two * i_sum + one)       *&
                        wigner_6j_symbol(j_1,j_3,2*rank,j_4,j_2,2*i_sum,.true.)
                        if (irez == 1) then
                           if (mod(j_3 - j_4 - 4*rank + 2*i_sum,4) /= 0)       &
                                                            coeff_one=-coeff_one
                        else if (irez == 2) then
                           if (mod(j_1 - j_2 - 4*rank + 2*i_sum,4) /= 0)       &
                                                          coeff_one = -coeff_one
                        end if
                        coeff = coeff + coeff_one
                     end if
                  end if
               end if
            end do
            if (abs(coeff) > eps10) then
               if (irez == 1) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%V       = coeff
               else if(irez == 2) then
                  nocoeff                      = nocoeff + 1
                  anco_V_list(nocoeff)%nu      = rank
                  anco_V_list(nocoeff)%a%kappa = subshell(no_two)%kappa
                  anco_V_list(nocoeff)%a%n     = subshell(no_two)%n
                  anco_V_list(nocoeff)%b%kappa = subshell(no_three)%kappa
                  anco_V_list(nocoeff)%b%n     = subshell(no_three)%n
                  anco_V_list(nocoeff)%c%kappa = subshell(no_four)%kappa
                  anco_V_list(nocoeff)%c%n     = subshell(no_four)%n
                  anco_V_list(nocoeff)%d%kappa = subshell(no_one)%kappa
                  anco_V_list(nocoeff)%d%n     = subshell(no_one)%n
                  anco_V_list(nocoeff)%V       = coeff
               end if
            end if
         end do
      end if
      !
   end subroutine anco_case_35_to_42
   !
   !
   subroutine anco_codeing(number,subshell,csf_r,csf_s,bra,ket)
   !--------------------------------------------------------------------
   ! Returns subshell states from scf_r and scf_s as appropriate to call
   ! the module rabs_rcfp.
   !
   ! Calls:  rcfp_get_term_number
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s 
      type(subshell_state), intent(out)      :: bra, ket
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: number
      !
      integer :: j, nu_bra, nu_ket, Nr_bra, Nr_ket, subshellJ_bra, &
                 subshellJ_ket
      logical :: fail
      !
       j             = (iabs(subshell(number)%kappa)*2) -1
       nu_bra        = csf_r%seniority(number)  * 1
       subshellJ_bra = csf_r%subshellJ(number)  * 1;   Nr_bra = -100
       bra%nq        = csf_r%occupation(number) * 1;   fail   = .false.
       !
       bra%state=rcfp_get_term_number(j,bra%nq,nu_bra,Nr_bra,subshellJ_bra,fail)
       if (fail) then
          print *,                                                             &
               "Unable to recognize the bra subshell state with (j,nu,J) = (", &
                j,",",nu_bra,",",subshellJ_bra,")."
          stop "anco_codeing(): program stop A."
       end if
       !
       nu_ket        = csf_s%seniority(number)  * 1
       subshellJ_ket = csf_s%subshellJ(number)  * 1;   Nr_ket = -100
       ket%nq        = csf_s%occupation(number) * 1;   fail   = .false.
       ket%state=rcfp_get_term_number(j,ket%nq,nu_ket,Nr_ket,subshellJ_ket,fail)
       if (fail) then
          print *,                                                             &
               "Unable to recognize the ket subshell state with (j,nu,J) = (", &
                j,",",nu_bra,",",subshellJ_bra,")."
          stop "anco_codeing(): program stop B."
       end if
       bra%n = subshell(number)%n * 1;   bra%subshellMQ = bra%nq -(j+1)/2    
       ket%n = bra%n * 1;                ket%subshellMQ = ket%nq -(j+1)/2
       !
   end subroutine anco_codeing
   !
   !
   subroutine anco_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the ANCO program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: selection
      logical            :: yes
      character(len=180) :: answer
      !
      anco_total_Vnu = 0;  anco_total_T = 0
      !
      print *, "Generate only non-trivial angular coefficients which"//&
               " include (at least one) open shells ?"
      yes = get_yes_stream()
      if (yes) anco_peel_shells_only = .true.
      !
      print *, "Generate one-electron angular coefficients for"//&
               " scalar interactions ?"
      yes = get_yes_stream()
      if (.not.yes) then
         anco_one_particle = .false.
      else
         print *, " Generate GRASP92-like T coefficients for scalar"//&
                  " interactions ?"
         yes = get_yes_stream()
         if (yes) anco_pure_one_particle  = .false.        
      end if
      !
      if(.not. anco_one_particle) then
         print *, "Generate one-electron angular coefficients for"//&
                  " non - scalar interactions ?"
         yes = get_yes_stream()
         if (.not.yes) then
            anco_one_particle_nonscalar = .false.
         else
            print *, " Generate GRASP92-like d coefficients for scalar"//&
                     " interactions ?"
            yes = get_yes_stream()
            if (yes) anco_pure_one_particle  = .false.        
     1      print *, " Enter the rank of the tensor"
            read (*,"(a)") answer
            answer = adjustl(answer)
            if (answer(1:1) == "q"  .or. answer(1:1) == "Q")   return
            if (answer(1:1) /= "1"  .and.  answer(1:1) /= "2"  .and.  &
                answer(1:1) /= "3"  .and.  answer(1:1) /= "4"  .and.  &
                answer(1:1) /= "5"  .and.  answer(1:1) /= "6"  .and.  &
                answer(1:1) /= "7"  .and.  answer(1:1) /= "8"  .and.  &
                answer(1:1) /= "9"  .and.  answer(1:1) /= "10" .and.  &
                answer(1:1) /= "0") goto 1
         !
            selection = iachar(answer(1:1)) - iachar("1") + 1
            if (selection < 0 .or. selection > 10 .or. len_trim(answer) > 1) then
               print*, " error in input  !!! "
               print*, " reenter "
               go to 1
            end if
            anco_one_particle_rank = selection
         end if
      else
         anco_one_particle_nonscalar = .false.
      end if
      !
      print *, "Generate two-electron angular coefficients for"//&
               " scalar interactions ?"
      yes = get_yes_stream()
      if (.not.yes) then
         anco_two_particle = .false.
      else
         print *, " Generate GRASP92-like V^k coefficients for scalar"//&
                  " interactions ?"
         yes = get_yes_stream()
         if (yes) anco_pure_two_particle  = .false.        
      end if
      !
   end subroutine anco_collect_input
   !
   !
   subroutine anco_diagonal_angle_s(csf_r,csf_s,nwshells,nwcore,subshell,      &
                                                          no_T_coeff,no_V_coeff)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for diagonal matrix elements with respect to configurations.
   !
   ! Calls: anco_case_1, anco_case_2_to_5
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s   
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore 
      integer, intent(inout)                 :: no_T_coeff, no_V_coeff
      !
      integer :: first_shell, second_shell, first_shell_min, second_shell_min
      logical :: not_core_close
      ! 
      if (.not. anco_peel_shells_only) then
         first_shell_min = 0;   second_shell_min = 0
         not_core_close = .true.
      else if (anco_pure_one_particle .or. anco_pure_two_particle) then
         first_shell_min = nwcore;   second_shell_min = nwcore
         not_core_close = .false.
      else if ((.not. anco_pure_one_particle) .or.                             &
               (.not. anco_pure_two_particle)) then
         first_shell_min = nwcore;   second_shell_min = nwcore
         not_core_close = .false.
      else
         stop " anco_diagonal_angle_s(): program stop A."
      end if
      do second_shell = second_shell_min+1, nwshells
         if (csf_s%occupation(second_shell) /= 0 ) then
            if (csf_s%occupation(second_shell) <                               &
                                   iabs(subshell(second_shell)%kappa)*2 .or.   &
                                                            not_core_close) then
               call anco_case_1(csf_r,csf_s,nwshells,second_shell_min,subshell,&
                                             no_T_coeff,no_V_coeff,second_shell)
               if (anco_two_particle) then
                  if (first_shell_min+1 < second_shell) then
                     do first_shell = first_shell_min+1, second_shell-1
                        if (csf_s%occupation(first_shell) /= 0 ) then
                           if (csf_s%occupation(first_shell) <                 &
                                   iabs(subshell(first_shell)%kappa)*2  .or.   &
                                                            not_core_close) then
                              call anco_case_2_to_5(csf_r,csf_s,nwshells,      &
                              second_shell_min,subshell,no_V_coeff,first_shell,&
                                                                   second_shell)
                           end if
                        end if
                     end do
                  end if
               end if
            end if
         end if
      end do
      !
   end subroutine anco_diagonal_angle_s
   !
   !
   subroutine anco_diff_occ_2(csf_r,csf_s,nwshells,nwcore,subshell,no_T_coeff, &
                                               no_V_coeff,creation,annihilation)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s of CSF
   ! for non-diagonal matrix elements with diff_occ = 2.
   !
   ! Calls: anco_case_7_to_14
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s   
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore
      integer, intent(in)                    :: creation, annihilation 
      integer, intent(inout)                 :: no_T_coeff, no_V_coeff
      !
      integer :: first_shell
      integer :: creation_one, annihilation_one, creation_two, annihilation_two
      logical :: check
      !
      if (anco_one_particle) then
         call anco_one_particle_scalart(csf_r,csf_s,nwshells,nwcore,subshell,   &
                                               no_T_coeff,creation,annihilation)
      end if
      if (.not. anco_two_particle) return
      do first_shell = 1, nwshells
         check = .true.
         if (creation > first_shell) then
           creation_one     =  first_shell
           creation_two     =  creation
         else if (creation == first_shell) then
           creation_one     =  creation
           creation_two     =  creation
         else
           creation_one     =  creation
           creation_two     =  first_shell
         end if
         if (annihilation > first_shell) then
           annihilation_one =  first_shell
           annihilation_two =  annihilation
           if (csf_s%occupation(first_shell) == 0 ) check = .false.
         else if (annihilation == first_shell) then
           annihilation_one =  annihilation
           annihilation_two =  annihilation
           if (csf_s%occupation(first_shell) <= 1 ) check = .false.
         else
           annihilation_one =  annihilation
           annihilation_two =  first_shell
           if (csf_s%occupation(first_shell) == 0 ) check = .false.
         end if
         if (check) then
           call anco_case_7_to_14(csf_r,csf_s,nwshells,nwcore,subshell,        &
                no_V_coeff,creation_one,creation_two,annihilation_one,         &
                                                               annihilation_two)
         end if
      end do
      !
   end subroutine anco_diff_occ_2
   !
   !
   subroutine anco_normal_form(no_one,no_two,no_three,no_form_one,no_form_two, &
                                                                  no_form_three)
   !--------------------------------------------------------------------
   ! Found the normal form of creation and annihilation operators.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)     :: no_one, no_two, no_three
      integer, intent(out)    :: no_form_one, no_form_two, no_form_three
      !
      no_form_one   = no_one
      no_form_three = no_one
      if (no_form_one   > no_two)   no_form_one   = no_two
      if (no_form_three < no_two)   no_form_three = no_two
      if (no_form_one   > no_three) no_form_one   = no_three
      if (no_form_three < no_three) no_form_three = no_three
      if ((no_one   > no_form_one) .and. (no_one   < no_form_three))no_form_two&
                                                                      = no_one
      if ((no_two   > no_form_one) .and. (no_two   < no_form_three))no_form_two&
                                                                      = no_two
      if ((no_three > no_form_one) .and. (no_three < no_form_three))no_form_two&
                                                                      = no_three
      !
   end subroutine anco_normal_form
   !
   !
   function anco_normal_phase(no_one,no_two,no_three,no_four)      result(phase)
   !--------------------------------------------------------------------
   ! Determinate the phase factor wich appear from permutation
   ! of operators of second quantization
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: no_one, no_two, no_three, no_four
      integer             :: phase
      !
      phase = 1
      if (no_one   > no_two)   phase = -phase
      if (no_one   > no_three) phase = -phase
      if (no_one   > no_four)  phase = -phase
      if (no_two   > no_three) phase = -phase
      if (no_two   > no_four)  phase = -phase
      if (no_three > no_four)  phase = -phase
      !
   end function anco_normal_phase
   !
   !
   subroutine anco_one_particle_diag(nu,csf_r,csf_s,nwshells,nwcore,subshell,  &
                                                          no_T_coeff)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients of non-scalar one-particle
   ! operator for a given pair r, s of CSF for diagonal matrix elements
   ! with respect to configurations.
   !
   ! Calls: anco_codeing, rcfp_calculate_Wk_me
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s   
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nu, nwshells, nwcore 
      integer, intent(inout)                 :: no_T_coeff
      !
      type(subshell_state) :: bra, ket
      real(kind=dp)        :: coeff
      integer :: first_shell, first_shell_min, delta_J
      logical :: not_core_close
      ! 
      if (.not. anco_peel_shells_only) then
         first_shell_min = 0;
         not_core_close = .true.
      else if (anco_pure_one_particle) then
         first_shell_min = nwcore;
         not_core_close = .false.
      else if (.not. anco_pure_one_particle) then
         first_shell_min = nwcore;
         not_core_close = .false.
      else
         stop " anco_one_particle_diag(): program stop A."
      end if
      do first_shell = first_shell_min+1, nwshells
         if (csf_s%occupation(first_shell) /= 0 ) then
            if (csf_s%occupation(first_shell) <                                 &
                iabs(subshell(first_shell)%kappa)*2  .or.  not_core_close) then
               if (abs(recoupling_matrix_check_nonscal(csf_r,csf_s,2*nu,        &
                        first_shell,first_shell,nwshells,nwcore)) > eps10) then
                  coeff = recoupling_matrix_1p_shells                           &
                    (csf_r,csf_s,2*nu,first_shell,.true.,nwshells,nwcore,delta_J)
                  if (delta_j /= 0) then
                     call anco_codeing(first_shell,subshell,csf_r,csf_s,bra,ket)
                     coeff = rcfp_calculate_Wk_me(bra,ket,nu,1,-1)
                     if (abs(coeff) > eps10) then
                     coeff = coeff * recoupling_matrix_1p_shells(csf_r,csf_s, &
                              2*nu,first_shell,.false.,nwshells,nwcore,delta_J)
                     no_T_coeff                      =no_T_coeff + 1
                     anco_T_list(no_T_coeff)%nu      =nu
                     anco_T_list(no_T_coeff)%a%kappa =subshell(first_shell)%kappa
                     anco_T_list(no_T_coeff)%a%n     =subshell(first_shell)%n
                     anco_T_list(no_T_coeff)%b%kappa =subshell(first_shell)%kappa
                     anco_T_list(no_T_coeff)%b%n     =subshell(first_shell)%n
                     anco_T_list(no_T_coeff)%T       =-coeff/ sqrt((two*nu+one))
                     end if
                  end if
               end if
            end if
         end if
      end do
      !
   end subroutine anco_one_particle_diag
   !
   !
   subroutine anco_one_particle_off(nu,csf_r,csf_s,nwshells,nwcore,subshell,  &
                                              no_T_coeff,creation,annihilation)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients of non-scalar one-particle
   ! operator for a given pair r, s of CSF for off diagonal matrix elements
   ! with respect to configurations.
   !
   ! Calls: anco_codeing, rcfp_calculate_Wk_me
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s   
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nu, nwshells, nwcore 
      integer, intent(in)                    :: creation, annihilation
      integer, intent(inout)                 :: no_T_coeff
      !
      type(subshell_state) :: bra_creation, ket_creation
      type(subshell_state) :: bra_annihilation, ket_annihilation
      real(kind=dp)        :: coeff, wa 
      integer              :: j_1, j_2, i, no_one, no_two, delta_J, occup
      ! 
      if (creation == annihilation) return
      no_one = min(creation,annihilation)
      no_two = max(creation,annihilation)
      wa     = abs(recoupling_matrix_check_nonscal(csf_r,csf_s,2*nu,        &
                        no_one,no_two,nwshells,nwcore))
      !
      if (abs(recoupling_matrix_check_nonscal(csf_r,csf_s,2*nu,             &
                        no_one,no_two,nwshells,nwcore)) < eps10) return
      j_1 =  iabs(subshell(no_one)%kappa) * 2 - 1
      j_2 =  iabs(subshell(no_two)%kappa) * 2 - 1
      coeff = recoupling_matrix_2p_shells(csf_r,csf_s,j_1,j_2,2*nu,no_one,  &
                                       no_two,.true.,nwshells,nwcore,delta_J)
      !
      if (delta_J == 0) return
      call anco_codeing(creation,subshell,csf_r,csf_s,bra_creation,ket_creation)
      call anco_codeing(annihilation,subshell,csf_r,csf_s,bra_annihilation, &
                                                            ket_annihilation)
      coeff = rcfp_calculate_a_me(bra_creation,ket_creation,1) *            &
              rcfp_calculate_a_me(bra_annihilation,ket_annihilation,-1)
      coeff = coeff * recoupling_matrix_2p_shells(csf_r,csf_s,j_1,j_2,2*nu, &
                               no_one,no_two,.false.,nwshells,nwcore,delta_J)
      occup = 0
      do i = no_one, no_two-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) coeff = -coeff
      if (abs(coeff) > eps10) then
         no_T_coeff                      =no_T_coeff + 1
         anco_T_list(no_T_coeff)%nu      =nu
         anco_T_list(no_T_coeff)%a%kappa =subshell(creation)%kappa
         anco_T_list(no_T_coeff)%a%n     =subshell(creation)%n
         anco_T_list(no_T_coeff)%b%kappa =subshell(annihilation)%kappa
         anco_T_list(no_T_coeff)%b%n     =subshell(annihilation)%n
         anco_T_list(no_T_coeff)%T       =-coeff/ sqrt((two*nu+one))
      end if
      !
   end subroutine anco_one_particle_off
   !
   !
   subroutine anco_one_particle_scalart(csf_r,csf_s,nwshells,nwcore,subshell,  &
                                                  nocoeff,creation,annihilation)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given pair r, s
   ! for one--particle scalar operator
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in)          :: csf_r, csf_s
      type(nkappa), dimension(:), intent(in) :: subshell
      integer, intent(in)                    :: nwshells, nwcore
      integer, intent(in)                    :: creation, annihilation
      integer, intent(inout)                 :: nocoeff
      !
      type(subshell_state) :: bra_creation, ket_creation
      type(subshell_state) :: bra_annihilation, ket_annihilation
      real(kind=dp)        :: coeff, coeff_1
      integer              :: i, j, no_one, no_two, delta_j, occup
      !
      if (creation == annihilation) return
      if (subshell(creation)%kappa /= subshell(annihilation)%kappa) return
      no_one = min(creation,annihilation)
      no_two = max(creation,annihilation)
      if (abs(recoupling_matrix_check(csf_r,csf_s,no_one,no_two,no_two,no_two, &
                                             nwshells,nwcore,1)) < eps10) return
      j =  iabs(subshell(creation)%kappa) * 2 - 1
      coeff_1 = recoupling_matrix_2_shells                                     &
                                    (csf_r,csf_s,j,no_one,no_two,.true.,delta_J)
      if (delta_J == 0) return
      call anco_codeing(creation,subshell,csf_r,csf_s,bra_creation,ket_creation)
      call anco_codeing(annihilation,subshell,csf_r,csf_s,bra_annihilation,    &
                                                               ket_annihilation)
      coeff = rcfp_calculate_a_me(bra_creation,ket_creation,1)                *&
              rcfp_calculate_a_me(bra_annihilation,ket_annihilation,-1)
      coeff = -coeff * recoupling_matrix_2_shells                              &
                                   (csf_r,csf_s,j,no_one,no_two,.false.,delta_J)
      occup = 0
      do i = no_one, no_two-1
         occup = occup + csf_r%occupation(i)
      end do
      if (mod(occup,2) == 0) coeff = -coeff
      if (abs(coeff) > eps10) then
         nocoeff = nocoeff + 1
         anco_T_list(nocoeff)%nu      = 0
         anco_T_list(nocoeff)%a%kappa = subshell(creation)%kappa
         anco_T_list(nocoeff)%a%n     = subshell(creation)%n
         anco_T_list(nocoeff)%b%kappa = subshell(annihilation)%kappa
         anco_T_list(nocoeff)%b%n     = subshell(annihilation)%n
         anco_T_list(nocoeff)%T       = coeff
      end if
      !
   end subroutine anco_one_particle_scalart
   !
   !
   subroutine anco_open_vnu(csf_set,file_formatted)
   !--------------------------------------------------------------------
   ! Opens a .vnu  Angular Coefficient output File on stream 25. 
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(inout) :: csf_set
      logical, intent(in)            :: file_formatted
      !
      integer :: ierr
      character(len=256) :: anco_vnu_file
      !
    1 print *, "Enter a file name for the  anco.vnu  file:"
      read *,  anco_vnu_file
      if (file_formatted) then
         call file_open(25,anco_vnu_file,"formatted  ","new",ierr)
      else
         call file_open(25,anco_vnu_file,"unformatted","new",ierr)
      end if
      !
      if (ierr /= 0) goto 1
      !
      if (file_formatted) then
         write(25,*) "ANCO"
         write(25,*) csf_set%nocsf, csf_set%nwshells, csf_set%nwcore
      else
         write(25) "ANCO"
         write(25) csf_set%nocsf, csf_set%nwshells, csf_set%nwcore
      end if
      !
   end subroutine anco_open_vnu
   !
   !
   subroutine anco_print_results(stream)
   !--------------------------------------------------------------------
   ! Writes the ... in a neat format to the given Fortran stream.
   !
   ! Calls:  
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      print *, stream
      !
   end subroutine anco_print_results
   !
   !
   subroutine anco_print_summary(csf_set)
   !--------------------------------------------------------------------
   ! Appends a short summary about the input to the  anco.sum  file
   ! which is open on stream 24.
   !--------------------------------------------------------------------
      !
      implicit none
      type(csf_basis), intent(in) :: csf_set
      !
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of 
      ! the  cowf.sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(24,*) "ANCO run at "//ctime(1:2)//":"//ctime(3:4)//":"//   &
                  ctime(5:6)//" on "//month//                          &
                  " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
      ! Write out the basic dimensions of the electron cloud
      write(24,*)
      write(24,*) "There are ",csf_set%number_of_electrons, &
                  " electrons in the cloud"
      write(24,*) " in ",csf_set%nocsf," relativistic CSFs"
      write(24,*) " based on ",csf_set%nwshells," relativistic subshells."
      write(24,*) " Total number of pair is:",total_number_pair
      write(24,*)
      !
      if (anco_peel_shells_only)then
         write(24,*) "Generate only not trivial angular coefficients"
      else
         write(24,*) "Generate all posible angular coefficients"
      end if
      !
      if (anco_one_particle  .or.  anco_one_particle_nonscalar) then
         if (anco_pure_one_particle)then
            write(24,*) " there are",anco_total_T,      &
                        " pure one-particle angular coefficients"
         else
            write(24,*) " there are",anco_total_T,      &
                        " Grasp92-like T coefficients"
         end if
      end if
      !
      if (anco_two_particle) then
         if (anco_pure_two_particle)then
            write(24,*) " there are",anco_total_Vnu,    &
                        " pure two-particle angular coefficients"
         else
            write(24,*) " there are",anco_total_Vnu,    &
                        " Grasp92-like Vk coefficients"
         end if
      end if
      write(24,*) " the total number of coefficients is",anco_total
      !
   end subroutine anco_print_summary
   !
end module rabs_anco
