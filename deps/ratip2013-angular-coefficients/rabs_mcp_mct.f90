module rabs_mcp
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module provides an 'interface' to the routines of the MCP/MCBP components and the NJGRAH library as they are 
! implemented in GRASP92. These parts of the program are considered to be 'obsolet' for the further development of the 
! RABS package and will be replaced in the future by a more efficient treatment of 'angular integration'. At the present, 
! however, these routines for calculating angular coefficients are required to generate continuum orbital in the mean-
! field of some given symmetry-adapted atomic state or for calculating Auger matrix elements, for instance.
! In this module, a CSF basis for a pair of CSF (two_csf) is defined into which all required quantum numbers are defined. 
! Only the angular coefficients for one pair of CSF can be calculated by a call of the MCP routines. The quantum numbers 
! of these CSF are 'used' by the (independent, non-module) procedures in rabs_mcp_adaptation to set the proper COMMON 
! arrays of the MCP package. This somehow sophisticated looking procedure ensures however, that COMMON arrays and other
! features (which are obsolet in Fortran 90) need only to appear in a small number of procedures in rabs_mcp_adaptation.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   implicit none
   !
   public  :: mcbp_generate_coefficients
                 ! Generates all (non-zero) MCBP coefficients for a given 'configuration scheme'.
   public  :: mcp_generate_coefficients
                 ! Generates all (non-zero) MCP coefficients for a given 'configuration scheme'.
   public  :: mcp_initialize_two_csf
                 ! Initializes an intermediate CSF scheme for two CSF in order to tranfer proper quantum numbers to the
                 ! MCP and NJGRAF components.
   public  :: mct_generate_coefficients
                 ! Generates all (non-zero) MCT coefficients for a given 'configuration scheme'.
   !
   ! Define data structure to allow for an adaptation of the GRASP92 procedures of the MCP and NJGRAF components.
   integer, save         :: csf_a_No, csf_b_No
   logical, save         :: mcp_common_need_initialization = .true.
   type(csf_basis), save :: two_csf
   !
   ! Define a structure to handle with MCP or MCBP coefficients
   type :: mcp_coefficient
      integer           :: r, s
      integer(kind=i2b) :: a, b, c, d, mu, nu
      real(kind=dp)     :: V
   end type mcp_coefficient
   !
   ! Define a structure to handle with MCT coefficients
   type, bind(c) :: mct_coefficient
      integer           :: r, s
      integer(kind=i2b) :: a, b, nu
      real(kind=dp)     :: T
   end type mct_coefficient
   !
   integer :: number_of_mcp_coefficients,  number_of_mcp_coefficients_max
   integer :: number_of_mcbp_coefficients, number_of_mcbp_coefficients_max
   !
   integer, bind(c) :: number_of_mct_coefficients
   !
   integer, parameter :: number_of_mct_coefficients_max = 1000
   !
   type (mcp_coefficient), dimension(:), allocatable :: mcp_list, mcp_list_save
   type (mcp_coefficient), dimension(:), allocatable :: mcbp_list
   !
   type (mct_coefficient), dimension(number_of_mct_coefficients_max), bind(c) :: mct_list
   !
   ! Define debug variable to support a few GRASP92 debug features
   integer, save :: debug_ibug1, debug_ibug2, debug_ibug3, debug_ibug4, debug_ibug5, debug_ibug6
   !
   logical :: debug_mcp_coefficients = .false.
   logical :: debug_mct_coefficients = .false.
   !
contains
   !
   subroutine cfp_1
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------------------------------------------------------
      !
   end subroutine cfp_1
   !
   !
   subroutine mcbp_generate_coefficients(rl,ru,sl,su,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Generates the (non-zero) MCBP coefficients for all pairs of CSF in the range   r = rl,...,ru and  s = sl,...,su  
   ! of the given 'configuration scheme'. The procedure applies an adaption of the MCBP component from  GRASP92 to the 
   ! Fortran 90 standard. All coefficients for the full configuration scheme or for the selected part are stored in the 
   ! array mcbp_list() and can be 'used' later from any other module of the RABS package..
   ! 
   ! Calls: mcp_initialize_two_csf(), breid(), breit(), rkco() [from GRASP92/MCP].
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)         :: rl, ru, sl, su
      type(csf_basis), intent(in) :: csf_set
      !
      integer  :: incor, r, s
      external    breid, breit, rkco
      !
      ! Allocate fresh storage for the MCP coefficients
      if (allocated(mcbp_list))  deallocate( mcbp_list )
      number_of_mcbp_coefficients = 0;   number_of_mcbp_coefficients_max = 1000
      allocate( mcbp_list(1:number_of_mcbp_coefficients_max) )
      !
      if (rl < 1   .or.  ru > csf_set%nocsf   .or.  sl < 1   .or.  su > csf_set%nocsf)  then
         stop "mcbp_generate_coefficients(): program stop A."
      end if
      !
      debug_ibug1 = 0
      !
      ! Accumulate contributions from the two-electron Coulomb operator
      if (debug_mcp_coefficients) then
         write (99,*)
         write (99,*) "                                     kt"
         write (99,*) "                                    v  (abcd)"
         write (99,*) "  r  s   a    b    c    d   k  t     rs"
         debug_ibug1 = 1
      end if
      !
      do  r = rl,ru
         do  s = sl,su
            call mcp_initialize_two_csf(r,s,csf_set)
            !
            incor = 1
            call rkco (r,s,breit,breid,incor)
         end do
      end do
      !
   end subroutine mcbp_generate_coefficients
   !
   !
   subroutine mcp_generate_coefficients(rl,ru,sl,su,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Generates the (non-zero) MCP coefficients for all pairs of CSF in the range   r = rl,...,ru and  s = sl,...,su  of 
   ! the given 'configuration scheme'. The procedure applies an adaption of the MCP component from GRASP92 to the 
   ! Fortran 90 standard. All coefficients for the full configuration scheme or for the selected part are stored in the 
   ! array mcp_list() and can be 'used' later from any other module of the RABS package..
   ! 
   ! Calls: mcp_initialize_two_csf(), rkco() [from GRASP92/MCP].
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)         :: rl, ru, sl, su
      type(csf_basis), intent(in) :: csf_set
      !
      integer  :: incor, r, s
      external    cor, cord, rkco
      !
      ! Allocate fresh storage for the MCP coefficients
      if (allocated(mcp_list))  deallocate( mcp_list )
      number_of_mcp_coefficients = 0;   number_of_mcp_coefficients_max = 1000
      allocate( mcp_list(1:number_of_mcp_coefficients_max) )
      !
      if (rl < 1   .or.  ru > csf_set%nocsf   .or.  sl < 1   .or.  su > csf_set%nocsf)  then
         stop "mcp_generate_coefficients(): program stop A."
      end if
      !
      debug_ibug1 = 0
      !
      ! Accumulate contributions from the two-electron Coulomb operator
      if (debug_mcp_coefficients) then
         write (99,*)
         write (99,*) "                                    k"
         write (99,*) "                                   V  (abcd)"
         write (99,*) "  r  s   a    c    b    d   k       rs"
         debug_ibug1 = 1
      end if
      !
      do  r = rl,ru
         do  s = sl,su
            call mcp_initialize_two_csf(r,s,csf_set)
            !
            incor = 1
            call rkco(r,s,cor,cord,incor)
         end do
      end do
      !
   end subroutine mcp_generate_coefficients
   !
   !
   subroutine mcp_initialize_two_csf(r,s,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Initializes some storage of quantum numbers for the calculation MCP coefficients for the pair r, s of CSF.
   ! 
   ! Calls: mcp_initialize_grasp92().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)         :: r, s
      type(csf_basis), intent(in) :: csf_set
      !
      integer :: nw
      !
      call deallocate_csf_basis(two_csf)
      !
      csf_a_No = r;   csf_b_No = s
      nw = csf_set%nwshells
      two_csf%nwshells = csf_set%nwshells
      two_csf%nocsf    = 2
      allocate( two_csf%csf(1:2), two_csf%subshell(1:csf_set%nwshells) )
      two_csf%subshell(1:nw)%n     = csf_set%subshell(1:nw)%n
      two_csf%subshell(1:nw)%kappa = csf_set%subshell(1:nw)%kappa
      allocate( two_csf%csf(1)%occupation(1:csf_set%nwshells),   two_csf%csf(1)%seniority(1:csf_set%nwshells),   &
                two_csf%csf(1)%subshellJ(1:csf_set%nwshells),    two_csf%csf(1)%subshellX(1:csf_set%nwshells),   &
                two_csf%csf(2)%occupation(1:csf_set%nwshells),   two_csf%csf(2)%seniority(1:csf_set%nwshells),   &
                two_csf%csf(2)%subshellJ(1:csf_set%nwshells),    two_csf%csf(2)%subshellX(1:csf_set%nwshells) )
      two_csf%csf(1)%totalJ           = csf_set%csf(r)%totalJ
      two_csf%csf(1)%parity           = csf_set%csf(r)%parity
      two_csf%csf(1)%occupation(1:nw) = csf_set%csf(r)%occupation(1:nw)
      two_csf%csf(1)%seniority(1:nw)  = csf_set%csf(r)%seniority(1:nw)
      two_csf%csf(1)%subshellJ(1:nw)  = csf_set%csf(r)%subshellJ(1:nw)
      two_csf%csf(1)%subshellX(1:nw)  = csf_set%csf(r)%subshellX(1:nw)
      !
      two_csf%csf(2)%totalJ           = csf_set%csf(s)%totalJ
      two_csf%csf(2)%parity           = csf_set%csf(s)%parity
      two_csf%csf(2)%occupation(1:nw) = csf_set%csf(s)%occupation(1:nw)
      two_csf%csf(2)%seniority(1:nw)  = csf_set%csf(s)%seniority(1:nw)
      two_csf%csf(2)%subshellJ(1:nw)  = csf_set%csf(s)%subshellJ(1:nw)
      two_csf%csf(2)%subshellX(1:nw)  = csf_set%csf(s)%subshellX(1:nw)
      !
      call mcp_initialize_grasp92()
      !
   end subroutine mcp_initialize_two_csf
   !
   !
   subroutine mct_generate_coefficients(rl,ru,sl,su,csf_set,iopar,rank)
   !--------------------------------------------------------------------------------------------------------------------
   ! Generates the (non-zero) MCT coefficients for all pairs of CSF in the range   r = rl,...,ru and  s = sl,...,su  of 
   ! the given 'configuration scheme'. The procedure applies an adaption of the MCT component from GRASP92 to the 
   ! Fortran 90 standard. All coefficients for the full configuration scheme or for the selected part are stored in the 
   ! array mct_list() and can be 'used' later from any other module of the RABS package..
   ! 
   ! Calls: mcp_initialize_two_csf(), rkco() [from GRASP92/MCP].
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)           :: rl, ru, sl, su, iopar, rank
      type(csf_basis), intent(in)   :: csf_set
      !
      integer                       :: ia, ib, r, s, ncf, nw_grasp92, pntriq
      real(kind=dp), dimension(149) :: tshell
      !
      common/orb2/ ncf, nw_grasp92, pntriq
      !
      number_of_mct_coefficients = 0
      !
      if (rl < 1   .or.  ru > csf_set%nocsf   .or.  sl < 1   .or.  su > csf_set%nocsf)  then
         stop "mct_generate_coefficients(): program stop A."
      end if
      !
      debug_ibug1 = 0
      !
      ! Accumulate contributions from the two-electron Coulomb operator
      if (debug_mct_coefficients) then
         write (99,*)
         write (99,*) "                          k"
         write (99,*) "                         T  (ab)"
         write (99,*) "  r  s   a    b   k       rs"
         debug_ibug1 = 1
      end if
      !
      do  r = rl,ru
         do  s = sl,su
            call mcp_initialize_two_csf(r,s,csf_set)
            !
            call tnsrjj (rank,iopar,r,s,ia,ib,tshell)
            if (ia .ne. 0) then
               if (ia .eq. ib) then
                  do 2  ia = 1,nw_grasp92
                     if (abs (tshell(ia)) .gt.  eps10) then
                        number_of_mct_coefficients = number_of_mct_coefficients + 1
                        if(number_of_mct_coefficients > number_of_mct_coefficients_max) then
                           stop "mct_generate_coefficients(): program stop B."
                        end if
                        !
                        mct_list(number_of_mct_coefficients)%r  = r 
                        mct_list(number_of_mct_coefficients)%s  = s 
                        mct_list(number_of_mct_coefficients)%a  = ia
                        mct_list(number_of_mct_coefficients)%b  = ia
                        mct_list(number_of_mct_coefficients)%nu = rank
                        mct_list(number_of_mct_coefficients)%T  = tshell(ia)
                     endif
                2 continue
               else
                  if (abs (tshell(1)) .gt. eps10) then
                     number_of_mct_coefficients = number_of_mct_coefficients + 1
                     if(number_of_mct_coefficients > number_of_mct_coefficients_max) then
                        stop "mct_generate_coefficients(): program stop C."
                     end if
                     !
                     mct_list(number_of_mct_coefficients)%r  = r 
                     mct_list(number_of_mct_coefficients)%s  = s 
                     mct_list(number_of_mct_coefficients)%a  = ia
                     mct_list(number_of_mct_coefficients)%b  = ib
                     mct_list(number_of_mct_coefficients)%nu = rank
                     mct_list(number_of_mct_coefficients)%T  = tshell(1)
                  endif
               endif
            endif
         end do
      end do
      !
   end subroutine mct_generate_coefficients
   !
end module rabs_mcp
