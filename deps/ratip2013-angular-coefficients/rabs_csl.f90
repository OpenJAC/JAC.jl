module rabs_csl
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module processes and maintains all informations concerning the present Configuration Symmetry List and the 
! representation of a number of atomic states in a CSF basis. It enables reading of such a list from a given file and set 
! all necessary arrays. Several abstract data types are defined to "keep" all data of a single configuration state 
! function as well as of a CSF and ASF basis simply together. It also defines an 'atomic basis' in terms of a determinant 
! representation which might be useful in different applications.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_determinant
   use rabs_dirac_orbital
   use rabs_functions_string
   implicit none
   !
   public  :: add_csf_to_basis
                 ! Appends all CSF with given symmetries from an 'old' to a 'new' configuration scheme. 
   public  :: combine_configuration_scheme
                 ! Combines two configuration schemes into a third one due to given cases.
   public  :: copy_basis_to_basis
                 ! Copies an 'old' configuration scheme into a 'new' one.
   public  :: copy_csf_to_csf
                 ! Copies a CSF from an 'old' configuration scheme as CSF of a 'new' one.
   public  :: deallocate_asf_basis
                 ! Deallocates the storage of type(asf_basis) including the allocation of all substructures.
   public  :: deallocate_asf_det_basis
                 ! Deallocates the storage of type(asf_det_basis) including the allocation of all substructures.
   public  :: deallocate_csf_basis
                 ! Deallocates the storage of type(csf_basis) including the allocation of all substructures.
   public  :: get_subshell_seniority
                 ! Returns the seniority quantum number for a given antisymmetric subshell state.
   public  :: load_csl_from_file
                 ! Reads in an parses the configuration symmetry list.
   public  :: load_mix_file_grasp92
                 ! Loads a GRASP92 .mix file into an appropriate data structure.
   private :: parse_intermediate_X
                 ! Parses the intermediate and total angular momenta and parity for a given CSF.
   private :: parse_subshell_J_and_seniority
                 ! Parses the subshell angular momenta J_sub for a given CSF.
   private :: parse_subshell_labels
                 ! Parses the subshell labels and occupations for a given CSF.
   public  :: print_configuration_scheme
                 ! Prints all data about a given CSF scheme.
   public  :: readwrite_asf_det_basis
                 ! Reads or writes a type(asf_det_basis) data structure
		 ! from or to a file.
   public  :: set_configuration_scheme
                 ! Creates a new 'configuration scheme' by coupling an additional electron to previously defined set 
                 ! of CSF.
   !
   type, public :: cs_function
      integer(kind=i1b) :: totalJ
      character(len=1)  :: parity
      integer(kind=i1b), dimension(:), pointer :: occupation
      integer(kind=i1b), dimension(:), pointer :: seniority
      integer(kind=i1b), dimension(:), pointer :: subshellJ
      integer(kind=i1b), dimension(:), pointer :: subshellX
   end type cs_function
   !
   type, public :: csf_basis
      integer :: nocsf         ! Number of CSF in the basis.
      integer :: nwshells      ! Number of (relativistic) subshells.
      integer :: nwcore        ! Number of (closed) core subshells.
      integer :: number_of_electrons
      type(nkappa), dimension(:), pointer      :: subshell
      type(cs_function), dimension(:), pointer :: csf
   end type csf_basis
   !
   type, public :: as_function
      integer           :: level_No
      integer(kind=i1b) :: totalJ, totalM
      character(len=1)  :: parity
      real(kind=dp)     :: energy
      real(kind=dp), dimension(:), pointer :: eigenvector
   end type as_function
   !
   type, public :: as_hfs_function
      integer		:: level_No, levelJ_No
      integer(kind=i1b) :: totalI, totalJ, totalF, totalM   ! totalM is only used if explicit M-dependent
      character(len=1)  :: parity
      real(kind=dp)	:: energy
      real(kind=dp), dimension(:), pointer :: eigenvector
   end type as_hfs_function
   !
   type, public :: asf_basis
      integer         :: noasf           ! Number of considered ASF.
      real(kind=dp)   :: average_energy  ! Averaged energy of this set of ASF.
      type(as_function), dimension(:), pointer :: asf
      type(csf_basis) :: csf_set
   end type asf_basis
   !
   type, public :: asf_hfs_ijf_basis
      integer         :: noasf_ijf,xxx       ! Number of considered IJF-coupled ASF.
      integer         :: nocsf, nwshells, number_of_electrons, noasf
      real(kind=dp)   :: average_energy  ! Averaged energy of this set of ASF.
      type(as_hfs_function), dimension(:), pointer :: asf
   end type asf_hfs_ijf_basis
   !
   type, public :: asf_det_basis
      integer         :: noasf           ! Number of considered ASF.
      real(kind=dp)   :: average_energy  ! Averaged energy of this set of ASF.
      type(as_function), dimension(:), pointer :: asf
      type(det_basis) :: det_set
   end type asf_det_basis
   !
   type, public :: asf_nrdet_basis
      integer         :: noasf           ! Number of considered ASF.
      real(kind=dp)   :: average_energy  ! Averaged energy of this set of ASF.
      type(as_function), dimension(:), pointer :: asf
      type(nrdet_basis) :: nrdet_set
   end type asf_nrdet_basis
   !
contains
   !
   subroutine add_csf_to_basis(csf_old,csf_new,totalJ,parity,index)
   !--------------------------------------------------------------------------------------------------------------------
   ! Appends all CSF with symmetries totalJ and parity from the 'configuration scheme' csf_old to csf_new. It assumes 
   ! and checkes that both configuration schemes have the 'same set and order' of subshells and the same number of 
   ! electrons. The (optional) argument index(:) returns for the appended CSF their original index i in csf_old%csf(i).
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(csf_basis), intent(in)    :: csf_old
      type(csf_basis), intent(inout) :: csf_new
      integer, intent(in)            :: totalJ
      character(len=1), intent(in)   :: parity
      integer, dimension(:), optional, intent(out)  :: index
      !
      integer                                     :: i, ii, j, n, nw, nwp
      type(cs_function), dimension(csf_new%nocsf) :: csf_sav
      !
      ! Check that the set and order of subshells and the number of electrons
      ! is conform for both CSF schemes
      if (rabs_use_stop) then
      if (csf_old%nwshells > csf_new%nwshells) then
         print *, "csf_old%nwshells, csf_new%nwshells = ", csf_old%nwshells, csf_new%nwshells
         stop "add_csf_to_basis(): program stop A."
      else if (csf_old%number_of_electrons /= csf_new%number_of_electrons) then
         stop "add_csf_to_basis(): program stop B."
      end if
      !
      do  i = 1,csf_old%nwshells
         if (csf_old%subshell(i)%n     /= csf_new%subshell(i)%n    .or. &
	     csf_old%subshell(i)%kappa /= csf_new%subshell(i)%kappa) then
            stop "add_csf_to_basis(): program stop C."
         end if
      end do
      end if
      !
      ! Determine the number of additional CSF to be appended
      n = 0
      do  i = 1,csf_old%nocsf
         if (csf_old%csf(i)%totalJ == totalJ  .and.  csf_old%csf(i)%parity == parity) n = n + 1
      end do
      !
      if (rabs_use_stop   .and.  n == 0) then
         stop "add_csf_to_basis(): program stop D."
      end if
      !
      csf_new%nwcore = min(csf_new%nwcore, csf_old%nwcore)
      !
      nw = csf_new%nwshells
      do  i = 1,csf_new%nocsf
         csf_sav(i)%totalJ = csf_new%csf(i)%totalJ
         csf_sav(i)%parity = csf_new%csf(i)%parity
	 allocate( csf_sav(i)%occupation(1:nw), csf_sav(i)%seniority(1:nw), &
                   csf_sav(i)%subshellJ(1:nw),  csf_sav(i)%subshellX(1:nw) )
	 csf_sav(i)%occupation(1:nw) = csf_new%csf(i)%occupation(1:nw)	   
	 csf_sav(i)%seniority(1:nw)  = csf_new%csf(i)%seniority(1:nw)
	 csf_sav(i)%subshellJ(1:nw)  = csf_new%csf(i)%subshellJ(1:nw)
	 csf_sav(i)%subshellX(1:nw)  = csf_new%csf(i)%subshellX(1:nw)
	 !
	 deallocate( csf_new%csf(i)%occupation, csf_new%csf(i)%seniority, &
	             csf_new%csf(i)%subshellJ,  csf_new%csf(i)%subshellX )
      end do
      !
      deallocate( csf_new%csf );   allocate( csf_new%csf(1:csf_new%nocsf+n) )
      !
      do  i = 1,csf_new%nocsf
         csf_new%csf(i)%totalJ = csf_sav(i)%totalJ
         csf_new%csf(i)%parity = csf_sav(i)%parity
	 allocate( csf_new%csf(i)%occupation(1:nw), csf_new%csf(i)%seniority(1:nw),  &
	           csf_new%csf(i)%subshellJ(1:nw),  csf_new%csf(i)%subshellX(1:nw) )
	 csf_new%csf(i)%occupation(1:nw) = csf_sav(i)%occupation(1:nw)	   
	 csf_new%csf(i)%seniority(1:nw)  = csf_sav(i)%seniority(1:nw)
	 csf_new%csf(i)%subshellJ(1:nw)  = csf_sav(i)%subshellJ(1:nw)
	 csf_new%csf(i)%subshellX(1:nw)  = csf_sav(i)%subshellX(1:nw)
	 !
	 deallocate( csf_sav(i)%occupation, csf_sav(i)%seniority, csf_sav(i)%subshellJ, csf_sav(i)%subshellX )
      end do
      !
      ! Append additional CSF
      ii = csf_new%nocsf;  nwp = csf_old%nwshells
      do  i = 1,csf_old%nocsf
         if (csf_old%csf(i)%totalJ == totalJ  .and. csf_old%csf(i)%parity == parity) then
            ii = ii + 1
            csf_new%csf(ii)%totalJ = totalJ
            csf_new%csf(ii)%parity = parity
	    allocate( csf_new%csf(ii)%occupation(1:nw), csf_new%csf(ii)%seniority(1:nw),  &
	              csf_new%csf(ii)%subshellJ(1:nw),  csf_new%csf(ii)%subshellX(1:nw) )
	    csf_new%csf(ii)%occupation(1:nwp) = csf_old%csf(i)%occupation(1:nwp)	   
	    csf_new%csf(ii)%seniority(1:nwp)  = csf_old%csf(i)%seniority(1:nwp)
	    csf_new%csf(ii)%subshellJ(1:nwp)  = csf_old%csf(i)%subshellJ(1:nwp)
	    csf_new%csf(ii)%subshellX(1:nwp)  = csf_old%csf(i)%subshellX(1:nwp)
	    !
            do  j = nwp+1,nw   
	       csf_new%csf(ii)%occupation(j) = 0	   
	       csf_new%csf(ii)%seniority(j)  = 0 
	       csf_new%csf(ii)%subshellJ(j)  = 0 
	       csf_new%csf(ii)%subshellX(j)  = csf_new%csf(ii)%subshellX(nwp)
            end do 
	    if (present(index)) then
	       index(ii) = i
	    end if
	 end if
      end do
      !
      if (rabs_use_stop   .and.   ii /= csf_new%nocsf + n) then
         stop "add_csf_to_basis(): program stop E."
      else
         csf_new%nocsf = csf_new%nocsf + n
      end if
      !
   end subroutine add_csf_to_basis
   !
   !
   subroutine combine_configuration_scheme(key,csf_a,csf_b,csf_c)
   !--------------------------------------------------------------------------------------------------------------------
   ! Combines the two configuration schemes csf_a and csf_b into the new scheme csf_c for which the memory is allocated 
   ! inside of this procedure. All configuration schemes must represent the same number of electrons.
   ! Several cases are distinguished for 'combination':
   !
   ! key == "Case A":  
   !    (i)  Total number of shells is nw = csf_a%nwshells + 1
   !    (ii) For the shell list: subshell(1:nw-1 from csf_a) and
   !                             subshell(nw) is the last shell of csf_b
   !                     
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)   :: key
      type(csf_basis), intent(in)    :: csf_a, csf_b
      type(csf_basis), intent(inout) :: csf_c
      !
      integer :: i, ii, nwb, nw, nocsf
      !
      if (csf_a%number_of_electrons /= csf_b%number_of_electrons) then
         stop "combine_configuration_scheme(): program stop A."
      else if (csf_a%nwshells < csf_b%nwshells) then
         stop "combine_configuration_scheme(): program stop B."
      end if
      !
      do  i = 1,csf_b%nwshells-1
         if (csf_a%subshell(i)%n     /= csf_b%subshell(i)%n    .or. &
	     csf_a%subshell(i)%kappa /= csf_b%subshell(i)%kappa) then
            stop "combine_configuration_scheme(): program stop C."
         end if
      end do
      !
      if (key == "Case A") then
         !
         nw    = csf_a%nwshells + 1
         nocsf = csf_a%nocsf + csf_b%nocsf
         !
         csf_c%nwshells = nw
         csf_c%nocsf    = nocsf
         csf_c%nwcore   = min(csf_a%nwcore, csf_b%nwcore)
         csf_c%number_of_electrons = csf_a%number_of_electrons
         !
         allocate( csf_c%subshell(1:nw), csf_c%csf(1:nocsf) )
         csf_c%subshell(1:nw-1) = csf_a%subshell(1:nw-1)
         csf_c%subshell(nw)     = csf_b%subshell(csf_b%nwshells)
         !
         do  i = 1,csf_a%nocsf
            csf_c%csf(i)%totalJ = csf_a%csf(i)%totalJ
            csf_c%csf(i)%parity = csf_a%csf(i)%parity
	    allocate( csf_c%csf(i)%occupation(1:nw), csf_c%csf(i)%seniority(1:nw),  &
	              csf_c%csf(i)%subshellJ(1:nw),  csf_c%csf(i)%subshellX(1:nw) )
	    csf_c%csf(i)%occupation(1:nw-1) = csf_a%csf(i)%occupation(1:nw-1)      
	    csf_c%csf(i)%seniority(1:nw-1)  = csf_a%csf(i)%seniority(1:nw-1)
	    csf_c%csf(i)%subshellJ(1:nw-1)  = csf_a%csf(i)%subshellJ(1:nw-1)
	    csf_c%csf(i)%subshellX(1:nw-1)  = csf_a%csf(i)%subshellX(1:nw-1)
	    !
	    csf_c%csf(i)%occupation(nw) = 0           
	    csf_c%csf(i)%seniority(nw)  = 0 
	    csf_c%csf(i)%subshellJ(nw)  = 0 
	    csf_c%csf(i)%subshellX(nw)  = csf_c%csf(i)%subshellX(nw-1)
         end do
         !
         ii = csf_a%nocsf;   nwb = csf_b%nwshells
         do  i = 1,csf_b%nocsf
            ii = ii + 1
            csf_c%csf(ii)%totalJ = csf_b%csf(i)%totalJ
            csf_c%csf(ii)%parity = csf_b%csf(i)%parity
	    allocate( csf_c%csf(ii)%occupation(1:nw), csf_c%csf(ii)%seniority(1:nw),  &
	              csf_c%csf(ii)%subshellJ(1:nw),  csf_c%csf(ii)%subshellX(1:nw) )
	    csf_c%csf(ii)%occupation(1:nwb-1) = csf_b%csf(i)%occupation(1:nwb-1)
	    csf_c%csf(ii)%seniority(1:nwb-1)  = csf_b%csf(i)%seniority(1:nwb-1)
	    csf_c%csf(ii)%subshellJ(1:nwb-1)  = csf_b%csf(i)%subshellJ(1:nwb-1)
	    csf_c%csf(ii)%subshellX(1:nwb-1)  = csf_b%csf(i)%subshellX(1:nwb-1)
	    !
	    csf_c%csf(ii)%occupation(nwb:nw-1) = 0           
	    csf_c%csf(ii)%seniority(nwb:nw-1)  = 0 
	    csf_c%csf(ii)%subshellJ(nwb:nw-1)  = 0 
	    csf_c%csf(ii)%subshellX(nwb:nw-1)  = csf_b%csf(i)%subshellX(nwb-1)
            !
	    csf_c%csf(ii)%occupation(nw) = csf_b%csf(i)%occupation(nwb)  
	    csf_c%csf(ii)%seniority(nw)  = csf_b%csf(i)%seniority(nwb) 
	    csf_c%csf(ii)%subshellJ(nw)  = csf_b%csf(i)%subshellJ(nwb)
	    csf_c%csf(ii)%subshellX(nw)  = csf_b%csf(i)%subshellX(nwb)
         end do
      else
         stop "combine_configuration_scheme(): program stop G."
      end if
      !
      !
   end subroutine combine_configuration_scheme
   !
   !
   subroutine copy_basis_to_basis(csf_old,csf_new)
   !--------------------------------------------------------------------------------------------------------------------
   ! Copies and 'old' configuration scheme into a 'new' one. The memory for the 'new' scheme is allocated inside of this 
   ! procedure.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(csf_basis), intent(in)    :: csf_old
      type(csf_basis), intent(inout) :: csf_new
      !
      integer :: i, nw, nocsf
      !
      nw    = csf_old%nwshells
      nocsf = csf_old%nocsf
      !
      csf_new%nwshells = csf_old%nwshells
      csf_new%nocsf    = csf_old%nocsf
      csf_new%nwcore   = csf_old%nwcore
      csf_new%number_of_electrons = csf_old%number_of_electrons
      !
      allocate( csf_new%subshell(1:nw), csf_new%csf(1:nocsf) )
      !
      csf_new%subshell(1:nw) = csf_old%subshell(1:nw)
      !
      do  i = 1,nocsf
         csf_new%csf(i)%totalJ = csf_old%csf(i)%totalJ
         csf_new%csf(i)%parity = csf_old%csf(i)%parity
	 allocate( csf_new%csf(i)%occupation(1:nw), csf_new%csf(i)%seniority(1:nw),  &
	           csf_new%csf(i)%subshellJ(1:nw),  csf_new%csf(i)%subshellX(1:nw) )
	 csf_new%csf(i)%occupation(1:nw) = csf_old%csf(i)%occupation(1:nw)      
	 csf_new%csf(i)%seniority(1:nw)  = csf_old%csf(i)%seniority(1:nw)
	 csf_new%csf(i)%subshellJ(1:nw)  = csf_old%csf(i)%subshellJ(1:nw)
	 csf_new%csf(i)%subshellX(1:nw)  = csf_old%csf(i)%subshellX(1:nw)
	 !
      end do
      !
   end subroutine copy_basis_to_basis
   !
   !
   subroutine copy_csf_to_csf(csf_old,i_old,csf_new,i_new)
   !--------------------------------------------------------------------------------------------------------------------
   ! Copies a CSF from an 'old' configuration scheme as CSF of a 'new' one. It is assumed that all the necessary arrays 
   ! have been allocated properly before.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)            :: i_old, i_new
      type(csf_basis), intent(in)    :: csf_old
      type(csf_basis), intent(inout) :: csf_new
      !
      integer :: i, nw
      !
      nw    = csf_old%nwshells
      !
      csf_new%csf(i_new)%totalJ = csf_old%csf(i_old)%totalJ
      csf_new%csf(i_new)%parity = csf_old%csf(i_old)%parity
      csf_new%csf(i_new)%occupation(1:nw) = csf_old%csf(i_old)%occupation(1:nw)      
      csf_new%csf(i_new)%seniority(1:nw)  = csf_old%csf(i_old)%seniority(1:nw)
      csf_new%csf(i_new)%subshellJ(1:nw)  = csf_old%csf(i_old)%subshellJ(1:nw)
      csf_new%csf(i_new)%subshellX(1:nw)  = csf_old%csf(i_old)%subshellX(1:nw)
      !
      do  i = nw+1,csf_new%nwshells
         csf_new%csf(i_new)%occupation(i) = 0
         csf_new%csf(i_new)%seniority(i)  = 0
         csf_new%csf(i_new)%subshellJ(i)  = 0
         csf_new%csf(i_new)%subshellX(i)  = csf_new%csf(i_new)%subshellX(nw)
      end do
      !
   end subroutine copy_csf_to_csf
   !
   !
   subroutine deallocate_asf_basis(asf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Deallocates the storage of asf_set [type(asf_basis)] including the allocation of all underlying substructures.
   !
   ! Calls: deallocate_csf_basis().
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(asf_basis), intent(inout) :: asf_set
      !
      integer :: i
      !
      do  i = 1,asf_set%noasf
         !! deallocate( asf_set%asf(i)%eigenvector )
      end do
      deallocate( asf_set%asf )
      !
      call deallocate_csf_basis(asf_set%csf_set)
      !
   end subroutine deallocate_asf_basis
   !
   !
   subroutine deallocate_asf_det_basis(asf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Deallocates the storage of asf_set [type(asf_det_basis)] including the allocation of all underlying substructures.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(asf_det_basis), intent(inout) :: asf_set
      integer :: i
      !
      do  i = 1,asf_set%noasf
         deallocate( asf_set%asf(i)%eigenvector )
      end do
      deallocate( asf_set%asf )
      !
      deallocate( asf_set%det_set%orbital )
      do  i = 1,asf_set%det_set%nodmax
         deallocate( asf_set%det_set%determinant(i)%occupation )
      end do
      deallocate( asf_set%det_set%determinant )
      !
   end subroutine deallocate_asf_det_basis
   !
   !
   subroutine deallocate_csf_basis(csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Deallocates the storage of csf_set [type(csf_basis)] including the allocation of all underlying substructures.
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(csf_basis), intent(inout) :: csf_set
      integer :: i, ios
      !
      do  i = 1,csf_set%nocsf
         deallocate( csf_set%csf(i)%occupation,csf_set%csf(i)%seniority, csf_set%csf(i)%subshellJ, csf_set%csf(i)%subshellX,&
		     stat=ios)
      end do
      !
      deallocate( csf_set%csf, csf_set%subshell, stat=ios )
      !
      csf_set%nocsf               = 0
      csf_set%nwshells            = 0
      csf_set%nwcore              = 0
      csf_set%number_of_electrons = 0
      !
   end subroutine deallocate_csf_basis
   !
   !
   function get_subshell_seniority(kappa,occupation,J2sub)    result(nu)
   !--------------------------------------------------------------------------------------------------------------------
   ! Return the seniority quantum number of a given antisymmetric subshell state by looking up an appropriate table. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kappa, occupation, J2sub
      integer             :: nu
      !
      nu = -1
      !
      select case(kappa)
      case(-1, 1)
         select case(occupation)
         case(0,2)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1)
            if (J2sub /= 1) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-2, 2)
         select case(occupation)
         case(0, 4)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 3)
            if (J2sub /= 3) goto 1
            nu = 1
         case(2)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4) then
                nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-3, 3)
         select case(occupation)
         case(0, 6)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1,5)
            if (J2sub /= 5) goto 1
            nu = 1
         case(2, 4)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case(3)
            if (J2sub == 5) then
               nu = 1; goto 2
            else if (J2sub == 3   .or.   J2sub == 9) then
               nu = 3; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-4, 4)
         select case(occupation)
         case(0, 8)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 7)
            if (J2sub /= 7) goto 1
            nu = 1
         case(2, 6)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or. &
                     J2sub == 12) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case(3, 5)
            if (J2sub == 7) then
               nu = 1; goto 2
            else if (J2sub == 3   .or.   J2sub == 5   .or.   J2sub == 9  .or. &
                     J2sub == 11  .or.   J2sub == 15) then
               nu = 3; goto 2
            else
               goto 1
            end if
         case(4)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 12) then
               nu = 2; goto 2
            else if (J2sub == 10  .or.   J2sub == 16) then
               nu = 4; goto 2
            else if (J2sub == 4   .or.   J2sub == 8)  then
                print *, "For j = 7/2, subshell states are not uniquely"//&
                         " defined by their Jsub und occupation numbers."
                stop "get_subshell_seniority(): program stop A."
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-5, 5)
         select case(occupation)
         case(0, 10)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 9)
            if (J2sub /= 9) goto 1
            nu = 1
         case(2, 8)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or. &
                     J2sub == 12  .or.   J2sub == 16) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-6, 6)
         select case(occupation)
         case(0, 12)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 11)
            if (J2sub /= 11) goto 1
            nu = 1
         case(2, 10)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or.   J2sub == 12  .or.&
                     J2sub == 16  .or.   J2sub == 20) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-7, 7)
         select case(occupation)
         case(0, 14)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 13)
            if (J2sub /= 13) goto 1
            nu = 1
         case(2, 12)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or.   J2sub == 12  .or.&
                     J2sub == 16  .or.   J2sub == 20  .or.   J2sub == 24) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-8, 8)
         select case(occupation)
         case(0, 16)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 15)
            if (J2sub /= 15) goto 1
            nu = 1
         case(2, 14)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or.   J2sub == 12  .or.&
                     J2sub == 16  .or.   J2sub == 20  .or.   J2sub == 24  .or.&
		     J2sub == 28 ) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-9, 9)
         select case(occupation)
         case(0, 18)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 17)
            if (J2sub /= 17) goto 1
            nu = 1
         case(2, 16)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or.   J2sub == 12  .or.&
                     J2sub == 16  .or.   J2sub == 20  .or.   J2sub == 24  .or.&
		     J2sub == 28  .or.   J2sub == 32 ) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-10, 10)
         select case(occupation)
         case(0, 20)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 19)
            if (J2sub /= 19) goto 1
            nu = 1
         case(2, 18)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or.   J2sub == 12  .or.&
                     J2sub == 16  .or.   J2sub == 20  .or.   J2sub == 24  .or.&
		     J2sub == 28  .or.   J2sub == 32  .or.   J2sub == 36 ) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-11, 11)
         select case(occupation)
         case(0, 22)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 21)
            if (J2sub /= 21) goto 1
            nu = 1
         case(2, 20)
            if (J2sub == 0) then
               nu = 0; goto 2
            else if (J2sub == 4   .or.   J2sub == 8   .or.   J2sub == 12  .or.&
                     J2sub == 16  .or.   J2sub == 20  .or.   J2sub == 24  .or.&
		     J2sub == 28  .or.   J2sub == 32  .or.   J2sub == 36  .or.&
		     J2sub == 40 ) then
               nu = 2; goto 2
            else
               goto 1
            end if
         case default;  goto 1
         end select
      case(-12, 12)
         select case(occupation)
         case(0, 24)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 23)
            if (J2sub /= 23) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-13, 13)
         select case(occupation)
         case(0, 26)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 25)
            if (J2sub /= 25) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-14, 14)
         select case(occupation)
         case(0, 28)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 27)
            if (J2sub /= 27) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-15, 15)
         select case(occupation)
         case(0, 30)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 29)
            if (J2sub /= 29) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-16, 16)
         select case(occupation)
         case(0, 32)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 31)
            if (J2sub /= 31) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-17, 17)
         select case(occupation)
         case(0, 34)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 33)
            if (J2sub /= 33) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-18, 18)
         select case(occupation)
         case(0, 36)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 35)
            if (J2sub /= 35) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-19, 19)
         select case(occupation)
         case(0, 38)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 37)
            if (J2sub /= 37) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-20, 20)
         select case(occupation)
         case(0, 40)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 39)
            if (J2sub /= 39) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-21, 21)
         select case(occupation)
         case(0, 42)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 41)
            if (J2sub /= 41) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-22, 22)
         select case(occupation)
         case(0, 44)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 43)
            if (J2sub /= 43) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-23, 23)
         select case(occupation)
         case(0, 46)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 45)
            if (J2sub /= 45) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-24, 24)
         select case(occupation)
         case(0, 48)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 47)
            if (J2sub /= 47) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-25, 25)
         select case(occupation)
         case(0, 50)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 49)
            if (J2sub /= 49) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-26, 26)
         select case(occupation)
         case(0, 52)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 51)
            if (J2sub /= 51) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-27, 27)
         select case(occupation)
         case(0, 54)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 53)
            if (J2sub /= 53) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-28, 28)
         select case(occupation)
         case(0, 56)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 55)
            if (J2sub /= 55) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-29, 29)
         select case(occupation)
         case(0, 58)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 57)
            if (J2sub /= 57) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-30, 30)
         select case(occupation)
         case(0, 60)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 59)
            if (J2sub /= 59) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-31, 31)
         select case(occupation)
         case(0, 62)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 61)
            if (J2sub /= 61) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-32, 32)
         select case(occupation)
         case(0, 64)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 63)
            if (J2sub /= 63) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-33, 33)
         select case(occupation)
         case(0, 66)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 65)
            if (J2sub /= 65) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-34, 34)
         select case(occupation)
         case(0, 68)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 67)
            if (J2sub /= 67) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-35, 35)
         select case(occupation)
         case(0, 70)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 69)
            if (J2sub /= 69) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-36, 36)
         select case(occupation)
         case(0, 72)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 71)
            if (J2sub /= 71) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-37, 37)
         select case(occupation)
         case(0, 74)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 73)
            if (J2sub /= 73) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-38, 38)
         select case(occupation)
         case(0, 76)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 75)
            if (J2sub /= 75) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-39, 39)
         select case(occupation)
         case(0, 78)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 77)
            if (J2sub /= 77) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-40, 40)
         select case(occupation)
         case(0, 80)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 79)
            if (J2sub /= 79) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-41, 41)
         select case(occupation)
         case(0, 82)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 81)
            if (J2sub /= 81) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-42, 42)
         select case(occupation)
         case(0, 84)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 83)
            if (J2sub /= 83) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-43, 43)
         select case(occupation)
         case(0, 86)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 85)
            if (J2sub /= 85) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-44, 44)
         select case(occupation)
         case(0, 88)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 87)
            if (J2sub /= 87) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-45, 45)
         select case(occupation)
         case(0, 90)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 89)
            if (J2sub /= 89) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-46, 46)
         select case(occupation)
         case(0, 92)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 91)
            if (J2sub /= 91) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-47, 47)
         select case(occupation)
         case(0, 94)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 93)
            if (J2sub /= 93) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-48, 48)
         select case(occupation)
         case(0, 96)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 95)
            if (J2sub /= 95) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-49, 49)
         select case(occupation)
         case(0, 98)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 97)
            if (J2sub /= 97) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-50, 50)
         select case(occupation)
         case(0, 100)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 99)
            if (J2sub /= 99) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-51, 51)
         select case(occupation)
         case(0, 102)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 101)
            if (J2sub /= 101) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-52, 52)
         select case(occupation)
         case(0, 104)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 103)
            if (J2sub /= 103) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-53, 53)
         select case(occupation)
         case(0, 106)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 105)
            if (J2sub /= 105) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-54, 54)
         select case(occupation)
         case(0, 108)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 107)
            if (J2sub /= 107) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-55, 55)
         select case(occupation)
         case(0, 110)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 109)
            if (J2sub /= 109) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-56, 56)
         select case(occupation)
         case(0, 112)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 111)
            if (J2sub /= 111) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-57, 57)
         select case(occupation)
         case(0, 114)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 113)
            if (J2sub /= 113) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-58, 58)
         select case(occupation)
         case(0, 116)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 115)
            if (J2sub /= 115) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-59, 59)
         select case(occupation)
         case(0, 118)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 117)
            if (J2sub /= 117) goto 1
            nu = 1
         case default;  goto 1
         end select
      case(-60, 60)
         select case(occupation)
         case(0, 120)
            if (J2sub /= 0) goto 1
            nu = 0
         case(1, 119)
            if (J2sub /= 119) goto 1
            nu = 1
         case default;  goto 1
         end select
      case default
         print *, "No proper case; kappa = ",kappa
         goto 1
      end select
    2 if (nu == -1) then
         print *, "nu == -1"
         goto 1
      end if
      !
      return
      !
    1 print *, "Improper input for determining a subshell seniority qantum number or |kappa| > 11"
      print *, " kappa,occupation,J2sub = ",kappa,occupation,J2sub
      stop "get_subshell_seniority(): program stop B."
      !      
   end function get_subshell_seniority
   !
   subroutine load_csl_from_file(csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in an parses the configuration symmetry list. A number of checks are made to ensure correctness and consistency 
   ! of the .csl file.    
   !
   ! Calls: angular_momentum_j(), get_integer_from_string(), get_kappa_from_name(), parse_intermediate_X(),
   !        parse_subshell_labels(), parse_subshell_J_and_seniority().
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(csf_basis), intent(inout) :: csf_set
      !
      integer                             :: i, j, k, kappa, nwcore, nwshells, nocsf, pqn, subsh
      logical                             :: single_configuration
      character(len=1000)                 :: recorda, recordb, recordc, recordd
      type(nkappa), dimension(:), pointer :: subshell_collect
      !
      single_configuration = .false.
      !
      print *, "Loading configuration symmetry list file ..."
      !
      read(21,"(a)") recorda
      read(21,*)
      read(21,"(a,a,a)") recordb, recordc, recordd
      read(21,*)
      print *, "recordc = ",recordc
      print *, "recordd = ",recordd
      !
      allocate( subshell_collect(500) )
      !
      ! Get the list of core and peel subshells
      nwshells = 0;   nwcore = 0
      do  subsh = 1,1000
         if (len(trim(recorda)) == 0) goto 2
         pqn   = get_integer_from_string(recorda(1:3))
         kappa = get_kappa_from_name(recorda(4:5))
         nwcore = nwcore + 1;    nwshells = nwshells + 1
         subshell_collect(nwshells)%n     = pqn
         subshell_collect(nwshells)%kappa = kappa
         recorda = recorda(6:1000)
      end do
      stop "load_csl_from_file(): program stop A."
      !
    2 do  subsh = 1,1000
         if (len(trim(recordb)) == 0) goto 13
         pqn   = get_integer_from_string(recordb(1:3))
         kappa = get_kappa_from_name(recordb(4:5))
         nwshells = nwshells + 1
         subshell_collect(nwshells)%n     = pqn
         subshell_collect(nwshells)%kappa = kappa
         recordb = recordb(6:1000)
      end do
      !
   13 do  subsh = 1,1000
         if (len(trim(recordc)) == 0) goto 14
         pqn   = get_integer_from_string(recordc(1:3))
         kappa = get_kappa_from_name(recordc(4:5))
         nwshells = nwshells + 1
         subshell_collect(nwshells)%n     = pqn
         subshell_collect(nwshells)%kappa = kappa
         recordc = recordc(6:1000)
      end do
      !
   14 do  subsh = 1,1000
         if (len(trim(recordd)) == 0) goto 3
         pqn   = get_integer_from_string(recordd(1:3))
         kappa = get_kappa_from_name(recordd(4:5))
         nwshells = nwshells + 1
         subshell_collect(nwshells)%n     = pqn
         subshell_collect(nwshells)%kappa = kappa
         recordd = recordd(6:1000)
      end do
      stop "load_csl_from_file(): program stop B."
      !
      ! Ensure that the sets of core and peel subshell are disjoint
    3 do  j = nwcore+1,nwshells
         pqn   = subshell_collect(j)%n
         kappa = subshell_collect(j)%kappa
         do  i = 1,nwcore
            if (subshell_collect(i)%n     == pqn   .and.  subshell_collect(i)%kappa == kappa)  then
               print *, "load_csl_from_file(): The lists of core and peel "
               stop " subshells must form disjoint sets."
            end if
         end do
      end do
      !
      ! Define the orbital functions of the current configuration scheme
      csf_set%nwshells = nwshells
      csf_set%nwcore   = nwcore
      allocate( csf_set%subshell(1:nwshells) )
      do  j = 1,nwshells
         csf_set%subshell(j)%n     = subshell_collect(j)%n
         csf_set%subshell(j)%kappa = subshell_collect(j)%kappa
      end do
      deallocate( subshell_collect )
      !
      if (nwshells > 1) then
         print *, " There are ",nwshells," relativistic subshells;"
      else
         print *, " There is 1 relativistic subshell;"
      end if
      !
      ! Count the number of csl in the .csl file and allocate memory for them in the current configuration scheme
      nocsf = 0
      do  i = 1,1000000000
         read(21,*,end=4);   read(21,*,err=7);   read(21,*,err=7)
         nocsf = nocsf + 1
      end do
      !
    4 if (nocsf == 0) then
         ! Set values for a single closed-shell configuration which is not set up by GRASP92
         csf_set%nocsf = 1;   nocsf = 1   
	 single_configuration = .true.
      else
         csf_set%nocsf = nocsf
      end if
      !
      allocate( csf_set%csf(1:nocsf) )
      do  i = 1,nocsf
         allocate( csf_set%csf(i)%occupation(1:nwshells),  csf_set%csf(i)%seniority(1:nwshells),   &
                   csf_set%csf(i)%subshellJ(1:nwshells),   csf_set%csf(i)%subshellX(1:nwshells)  )
         !
         do  j = 1,nwcore
            csf_set%csf(i)%occupation(j) =  angular_momentum_j(csf_set%subshell(j)%kappa) + 1
            csf_set%csf(i)%seniority(j)  = 0
            csf_set%csf(i)%subshellJ(j)  = 0
            csf_set%csf(i)%subshellX(j)  = 0
         end do
         !
         do  j = nwcore+1,nwshells
            csf_set%csf(i)%occupation(j) = 0
            csf_set%csf(i)%seniority(j)  = 0
            csf_set%csf(i)%subshellJ(j)  = 0
            csf_set%csf(i)%subshellX(j)  = 0
         end do
      end do
      !
      if (single_configuration) then
         csf_set%number_of_electrons = sum(csf_set%csf(1)%occupation(1:csf_set%nwshells))
	 csf_set%csf(1)%totalJ = 0;   csf_set%csf(1)%parity = "+"
         return
      end if
      !
      ! Rewind the .csl file and parse all CSF in turn
      rewind 21;   
      read(21,*);   read(21,*);   read(21,*);   read(21,*);   read(21,*)
      !
      do  i = 1,nocsf
         read(21,"(a)",end=6) recorda
         read(21,"(a)")       recordb
         read(21,"(a)")       recordc
         call parse_subshell_labels(recorda,i,csf_set)
         call parse_subshell_J_and_seniority(recordb,i,csf_set)
         call parse_intermediate_X(recordc,i,csf_set)
         !
         ! Check if this CSF was already in the list
         if (i > 1) then
            do  k = 1,i-1
               do  j = 1,nwshells
                  if ( &
                   csf_set%csf(k)%occupation(j) /= csf_set%csf(i)%occupation(j) .or.&
                   csf_set%csf(k)%seniority(j)  /= csf_set%csf(i)%seniority(j)  .or.&
                   csf_set%csf(k)%subshellJ(j)  /= csf_set%csf(i)%subshellJ(j)  .or.&
                   csf_set%csf(k)%subshellX(j)  /= csf_set%csf(i)%subshellX(j))     then
                     goto 5
                  end if
               end do
               if (csf_set%csf(k)%totalJ /= csf_set%csf(i)%totalJ   .or. &
                   csf_set%csf(k)%parity /= csf_set%csf(i)%parity)  goto 5
            end do
            print *, "CSF i has been occured in the list before; i = ",i
            stop "load_csl_from_file(): program stop C."
         end if
    5    continue
      end do
      !
      ! Check a consistent number of electrons in all CSF
    6 csf_set%number_of_electrons = &
         sum(csf_set%csf(1)%occupation(1:csf_set%nwshells))
      do  i = 1,nocsf
         if (sum(csf_set%csf(i)%occupation(1:csf_set%nwshells)) /= csf_set%number_of_electrons) then
            print *, "ne, occu = ",csf_set%number_of_electrons, &
                     csf_set%csf(i)%occupation(1:csf_set%nwshells)
            stop "load_csl_from_file(): program stop D."
         end if
      end do
      !
      print *, " there are ",nocsf," relativistic CSFs;"
      print *, " ... load complete."
      return
      !
    7 stop "load_csl_from_file(): program stop E."
      !
   end subroutine load_csl_from_file
   !
   !
   subroutine load_mix_file_grasp92(asf_set,is_formatted,ierr)
   !--------------------------------------------------------------------------------------------------------------------
   ! Read in the eigenvectors of a number of ASF in a given ASF basis. It assumes and checks that the CSF basis of 
   ! asf_set has been read  before. For is_formatted = .true., is assumes a formatted .mix file and an (unformatted) 
   ! original GRASP92 .mix file otherwise. An error code ierr =/= 0 is returned if the routine was not successful.
   ! This file is always attached to stream 25.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(out) :: ierr
      logical, intent(in)  :: is_formatted
      type(asf_basis), intent(inout) :: asf_set
      !
      integer                              :: i, j, noasf, noasf_sofar, test_no_electrons, test_nocsf, test_nwshells
      integer, dimension(:), pointer       :: iatjpo, iaspar
      real(kind=dp), dimension(:), pointer :: eval
      !
      ierr = 0
      !
      ! Check consistency with the CSF basis
      if (is_formatted) then
         read (25,"(3i6)") test_no_electrons,test_nocsf,test_nwshells
      else
         read (25) test_no_electrons,test_nocsf,test_nwshells
      end if
      if (test_no_electrons /= asf_set%csf_set%number_of_electrons   .or.  test_nocsf /= asf_set%csf_set%nocsf  .or. &
          test_nwshells     /= asf_set%csf_set%nwshells) then
	 ierr = 1
         return
      end if
      !
      ! Allocate memory for the set of ASF and load the data from the .mix file
      print *, "Loading mixing coefficients file ..."
      if (is_formatted) then
         read (25,"(i6)") asf_set%noasf
      else
         read (25) asf_set%noasf
      end if
      allocate( asf_set%asf(1:asf_set%noasf) )
      do  i = 1,asf_set%noasf
         allocate( asf_set%asf(i)%eigenvector(1:asf_set%csf_set%nocsf) )
      end do
      allocate( iatjpo(1:asf_set%noasf), iaspar(1:asf_set%noasf), &
                eval(1:asf_set%noasf) ) 
      !
      if (is_formatted) then
         noasf = min(60,asf_set%noasf)
         read (25,"(60i6)")      (asf_set%asf(i)%level_No,i = 1,noasf)
         read (25,"(60(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, i = 1,noasf)
         read (25,"(61e26.19)")   asf_set%average_energy, (eval(i),i = 1,noasf)
         do  i = 1,noasf
            asf_set%asf(i)%energy = asf_set%average_energy + eval(i)
         end do
	 !
         do  j = 1,asf_set%csf_set%nocsf
            read (25,"(60e16.9)") (asf_set%asf(i)%eigenvector(j),i=1,noasf)
         end do
         !
         if (asf_set%noasf <= 60) goto 10
         !
         noasf_sofar = 0
       1 noasf_sofar = noasf_sofar + 60
         !
         noasf = min(noasf_sofar+60,asf_set%noasf)
         read (25,"(60i6)") (asf_set%asf(i)%level_No,i = noasf_sofar+1,noasf)
         read (25,"(60(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, i = noasf_sofar+1,noasf)
         !
         read (25,"(60e26.19)") (eval(i),i = noasf_sofar+1,noasf)
         do  i = noasf_sofar+1,noasf
            asf_set%asf(i)%energy = asf_set%average_energy + eval(i)
         end do
         !
         do  j = 1,asf_set%csf_set%nocsf
            read (25,"(60e16.9)") (asf_set%asf(i)%eigenvector(j), i=noasf_sofar+1,noasf)
         end do
         !
         if (asf_set%noasf > noasf_sofar+60) goto 1
         !
      else
         read (25) (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
         read (25) (iatjpo(i),iaspar(i),i = 1,asf_set%noasf)
         do  i = 1,asf_set%noasf
            asf_set%asf(i)%totalJ = iatjpo(i) - 1
            if (iaspar(i) == -1) then
               asf_set%asf(i)%parity = "-"
            else
               asf_set%asf(i)%parity = "+"
            end if
         end do
         read (25) asf_set%average_energy,(eval(i),i = 1,asf_set%noasf)
         do  i = 1,asf_set%noasf
            asf_set%asf(i)%energy = asf_set%average_energy + eval(i)
         end do
         !
         read (25) ((asf_set%asf(i)%eigenvector(j),j=1,asf_set%csf_set%nocsf), i=1,asf_set%noasf)
      end if
      !
      !
   10 continue
      !
      deallocate( iatjpo, iaspar, eval ) 
      print *, ' ... load complete;'
      close (25)
      !
   end subroutine load_mix_file_grasp92
   !
   !
   subroutine parse_intermediate_X(record,ncsf,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Parses the intermediate and total angular momenta and parity on the character string record. These are the angular 
   ! momenta of the CSF ncsf in the configuration scheme csf_set.
   !
   ! Calls: angular_momentum_j(), get_dinteger_from_string().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(inout) :: record
      integer, intent(in)             :: ncsf
      type(csf_basis), intent(inout)  :: csf_set
      !
      integer :: j
      logical :: open
      !
      ! Get total angular momentum and parity by looking 'backward' on the string
      j = scan(record,"+-",.true.)
      if (rabs_use_stop   .and.  j == 0) then
         stop "parse_intermediate_X(): program stop A."
      end if
      csf_set%csf(ncsf)%parity = record(j:j)
      csf_set%csf(ncsf)%totalJ = get_dinteger_from_string(record(j-4:j-1))
      !
      open = .false.
      do j = csf_set%nwcore+1,csf_set%nwshells
         !
         if (csf_set%csf(ncsf)%occupation(j) == 0) then
            ! For unoccupied shells, the intermediate X remains the same
            ! as from previous subshells
            if (j > 1) then
               csf_set%csf(ncsf)%subshellX(j) = csf_set%csf(ncsf)%subshellX(j-1)
            else
               csf_set%csf(ncsf)%subshellX(j) = 0   
            end if
            cycle
            !
         else if (csf_set%csf(ncsf)%occupation(j) < (angular_momentum_j(csf_set%subshell(j)%kappa)+1) .and. &
                  (open .eqv. .false.)) then
            ! For the 'first' open shell take X from the subshell J
            open = .true.
            csf_set%csf(ncsf)%subshellX(j) = csf_set%csf(ncsf)%subshellJ(j)
         else if (csf_set%csf(ncsf)%occupation(j) < &
                  (angular_momentum_j(csf_set%subshell(j)%kappa)+1) .and. &
                  (open .eqv. .true.)   .and.   record(10:12) == "   ") then
            ! There is no explicit X value for subshell-J = 0 even for
            ! an open subshell
            csf_set%csf(ncsf)%subshellX(j) =  csf_set%csf(ncsf)%subshellX(j-1)
         else if (csf_set%csf(ncsf)%occupation(j) < angular_momentum_j(csf_set%subshell(j)%kappa)+1) then
            ! For the second, ... open shell (with subshell-j /= 0) there must be a corresponding X value
            if (record(12:12) == "/") then
               csf_set%csf(ncsf)%subshellX(j) = get_dinteger_from_string(record(10:13))
            else if (record(11:11) == "/") then
               csf_set%csf(ncsf)%subshellX(j) = get_dinteger_from_string(record(9:12))
            else if (record(13:13) == "+"  .or.  record(13:13) == "-") then
               csf_set%csf(ncsf)%subshellX(j) = get_dinteger_from_string(record(10:12))
            else if (record(12:12) == "+"  .or.  record(12:12) == "-") then
               csf_set%csf(ncsf)%subshellX(j) = get_dinteger_from_string(record(10:11))
            else if (record(11:11) == "+"  .or.  record(11:11) == "-") then
               csf_set%csf(ncsf)%subshellX(j) = get_dinteger_from_string(record(10:10))
            else 
               csf_set%csf(ncsf)%subshellX(j) = get_dinteger_from_string(record(10:12))
            end if
         else if (csf_set%csf(ncsf)%occupation(j) == angular_momentum_j(csf_set%subshell(j)%kappa)+1) then
            ! For fully occupied shells take the X from the previous subshell
            if (j > 1) then
               csf_set%csf(ncsf)%subshellX(j) =  csf_set%csf(ncsf)%subshellX(j-1)
            else
               csf_set%csf(ncsf)%subshellX(j) = 0
            end if
         else if (rabs_use_stop) then 
            stop "parse_intermediate_X(): program stop B."
         end if
         record = record(10:256)
      end do
      csf_set%csf(ncsf)%totalJ = csf_set%csf(ncsf)%subshellX(csf_set%nwshells)
      !
   end subroutine parse_intermediate_X
   !
   !
   subroutine parse_subshell_J_and_seniority(record,ncsf,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Parses the subshell angular momenta J_sub for a given CSF ncsf on the character string record.
   !
   ! Calls: angular_momentum_j(), get_dinteger_from_string().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(inout) :: record
      integer, intent(in)             :: ncsf
      type(csf_basis), intent(inout)  :: csf_set
      !
      integer :: j, J2sub, occ, pos_semicolon
      !
      do  j = csf_set%nwcore+1,csf_set%nwshells
         if (csf_set%csf(ncsf)%occupation(j) == 0) cycle
         if (len_trim(record) == 0) return 
         if (csf_set%csf(ncsf)%occupation(j) == 0) then
            !
         else if (csf_set%csf(ncsf)%occupation(j) <  angular_momentum_j(csf_set%subshell(j)%kappa)+1) then
            csf_set%csf(ncsf)%subshellJ(j) = get_dinteger_from_string(record(6:9))
         end if
         !
         occ   = csf_set%csf(ncsf)%occupation(j)
         J2sub = csf_set%csf(ncsf)%subshellJ(j)
         if (angular_momentum_j(csf_set%subshell(j)%kappa) == 7   .and.  occ == 4  .and.  &
             (J2sub == 4  .or.  J2sub == 8) ) then
            pos_semicolon = scan(record(1:9),";")
            if (rabs_use_stop   .and.   pos_semicolon == 0) then
               stop "parse_subshell_J_and_seniority(): program stop A."
            end if
            csf_set%csf(ncsf)%seniority(j) = get_integer_from_string(record(pos_semicolon-2:pos_semicolon-1))
         else
            csf_set%csf(ncsf)%seniority(j) = get_subshell_seniority(csf_set%subshell(j)%kappa,occ,J2sub)
         end if
         record = record(10:256)
      end do
      !
   end subroutine parse_subshell_J_and_seniority
   !
   !
   subroutine parse_subshell_labels(record,ncsf,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Parses the subshell labels and the occupation of the (open) subshells as decribed on the character string record. 
   ! These are the labels and the occupation of the CSF ncsf in the configuration scheme csf_set.
   !
   ! Calls: get_integer_from_string(), get_kappa_from_name().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(inout) :: record
      integer, intent(in)             :: ncsf
      type(csf_basis), intent(inout)  :: csf_set
      !
      integer :: j
      !
      do  j = csf_set%nwcore+1,csf_set%nwshells
         if (len(trim(record)) == 0) exit
         if (csf_set%subshell(j)%n == get_integer_from_string(record(1:3)) .and. &
             csf_set%subshell(j)%kappa == get_kappa_from_name(record(4:5)) ) then
            csf_set%csf(ncsf)%occupation(j) = get_integer_from_string(" "//record(7:8))
         else
            cycle
         end if
         record = record(10:256)
      end do
      !
   end subroutine parse_subshell_labels
   !
   !
   subroutine print_configuration_scheme(stream,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Prints all information about the CSF scheme csf_set in a neat format on stream.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)            :: stream
      type(csf_basis), intent(in)    :: csf_set
      !
      integer :: i
      !
      write(stream,*) " "
      write(stream,*) "The current configuration scheme with ",csf_set%nocsf, " CSF is defined as follows:"
      write(stream,*) " "
      write(stream,*) " Number of (relativistic) subshells: ",csf_set%nwshells
      write(stream,*) " Number of (closed) core subshells:  ",csf_set%nwcore
      write(stream,*) " Number of electrons: ",csf_set%number_of_electrons
      write(stream,*) " "
      write(stream,*) " Subshells:"
      write(stream,1) (csf_set%subshell(i)%n,csf_set%subshell(i)%kappa, i=1,csf_set%nwshells)
      do  i = 1,csf_set%nocsf
         write(stream,*) " "
         write(stream,*) " CSF ",i,":   (",csf_set%csf(i)%totalJ,"/2", csf_set%csf(i)%parity,")"
         write(stream,*) csf_set%csf(i)%occupation(1:csf_set%nwshells)
         write(stream,*) csf_set%csf(i)%seniority(1:csf_set%nwshells)
         write(stream,*) csf_set%csf(i)%subshellJ(1:csf_set%nwshells)
         write(stream,*) csf_set%csf(i)%subshellX(1:csf_set%nwshells)
      end do
      write(stream,*) " "
    1 format(20(3x,i2,2x,i2))
      !
   end subroutine print_configuration_scheme
   !
   !
   subroutine readwrite_asf_det_basis(stream,read_from,asf)             
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads or writes a type(asf_det_basis) data structure from or to a file on stream.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                :: stream
      logical, intent(in)                :: read_from
      type(asf_det_basis), intent(inout) :: asf
      !
      integer       :: i, j
      !
      if (read_from) then
         read(stream,*)  asf%noasf, asf%average_energy
	 allocate( asf%asf(1:asf%noasf) )
	 call readwrite_det_basis(stream,read_from,asf%det_set)
      else
         write(stream,*) asf%noasf, asf%average_energy
	 call readwrite_det_basis(stream,read_from,asf%det_set)
      end if
      !
      do  i = 1,asf%noasf
         if (read_from) then
	    read(stream,*) asf%asf(i)%level_No, asf%asf(i)%totalJ, asf%asf(i)%parity, asf%asf(i)%energy
	    allocate( asf%asf(i)%eigenvector(1:asf%det_set%nodmax) )
            do  j = 1,asf%det_set%nodmax
               read(stream,*) asf%asf(i)%eigenvector(j)
            end do
	 else
	    write(stream,*) asf%asf(i)%level_No, asf%asf(i)%totalJ, asf%asf(i)%parity,   asf%asf(i)%energy
            do  j = 1,asf%det_set%nodmax
               write(stream,*) asf%asf(i)%eigenvector(j)
            end do
	 end if
      end do
      !
   end subroutine readwrite_asf_det_basis
   !
   !
   subroutine set_configuration_scheme(csf_old,csf_new,pqn,kappa, totalJ_old,parity_old,totalJ,parity,append,index)
   !--------------------------------------------------------------------------------------------------------------------
   ! Set up a (relativistic) 'configuration scheme' in arrays which are appropriate for calls to other procedures 
   ! including those for angular coefficients. It 'copies' all CSF with totalJ_old and parity_old from csf_old 
   ! [of type(csf_basis)] and extends these CSF by a (pqn,kappa)-subshell so that they become CSF with symmetries of
   ! totalJ and parity. If the (optional) argument append is present and append = .true. then these CSF are 'added' 
   ! to a previously established CSF basis; it is assumed that these CSF fit together by their number_of_electrons, 
   ! number of subshells, ... A test is made only after the 'new' CSF have been appended. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(csf_basis), intent(in)    :: csf_old
      type(csf_basis), intent(inout) :: csf_new
      integer, intent(in)            :: pqn, kappa, totalJ_old, totalJ
      character(len=1), intent(in)   :: parity_old, parity
      logical, optional, intent(in)  :: append
      integer, dimension(:), optional, intent(out)  :: index
      !
      integer         :: i, ipc, j, l, na, nocsf, nocsf_add, nw
      !
      ! First, check that a kappa-subshell is conform with the given symmetries
      if (parity_old == "+") then;   ipc = 1;   else;   ipc = -1; end if
      l = angular_momentum_l(kappa);   j = angular_momentum_j(kappa)
      if ((ipc * ((-1)**l) ==  1  .and.  parity == "-")  .or.  (ipc * ((-1)**l) == -1  .and.  parity == "+")  .or. &
           abs(totalJ_old - j) > totalJ  .or.  totalJ_old + j < totalJ)  then
         print *, "totalJ_old,parity_old,kappa,totalJ,parity = ",totalJ_old,parity_old,kappa,totalJ,parity
         stop "set_configuration_scheme(): program stop A."
      end if
      !
      ! Determine the number of 'additional' CSF
      nocsf_add = 0
      do  i = 1,csf_old%nocsf
         if (csf_old%csf(i)%totalJ == totalJ_old   .and.  csf_old%csf(i)%parity == parity_old)  then
            nocsf_add = nocsf_add + 1
         end if
      end do
      if (rabs_use_stop   .and.   nocsf_add == 0) then
         stop "set_configuration_scheme(): program stop B."
      end if
      !
      if (present(append)   .and.   append) then
         !
         ! Stop here if the 'additional' subshell already appears in the list
         ! of subshells in csf_new
         do  i = 1,csf_new%nwshells
            if (csf_new%subshell(i)%n == pqn  .and.  csf_new%subshell(i)%kappa == kappa) then
               stop "set_configuration_scheme(): program stop C."
            end if
         end do
         !
         ! Stop if subshell 1..csf_old%nwshells do not appear in the same order
         do  i = 1,csf_old%nwshells
            if (csf_new%subshell(i)%n     /= csf_old%subshell(i)%n  .or.   &
                csf_new%subshell(i)%kappa /= csf_old%subshell(i)%kappa) then
               stop "set_configuration_scheme(): program stop D."
            end if
         end do
         ! 
         if (csf_new%number_of_electrons /= csf_old%number_of_electrons+1) then
            stop "set_configuration_scheme(): program stop E."
         end if
         !
         stop "set_configuration_scheme(): program stop F."
         !         
      else
         nw = csf_old%nwshells
         csf_new%nocsf    = nocsf_add
         csf_new%nwshells = csf_old%nwshells + 1
         csf_new%nwcore   = csf_old%nwcore
         csf_new%number_of_electrons = csf_old%number_of_electrons + 1
         allocate( csf_new%subshell(1:csf_new%nwshells), csf_new%csf(1:nocsf_add) )
         csf_new%subshell(1:nw)             = csf_old%subshell(1:nw)
         csf_new%subshell(csf_new%nwshells) = nkappa(pqn,kappa)
         !
         na = 0
         do  i = 1,csf_old%nocsf
            if (csf_old%csf(i)%totalJ == totalJ_old   .and.  csf_old%csf(i)%parity == parity_old)  then
               na = na + 1
               csf_new%csf(na)%totalJ = totalJ
               csf_new%csf(na)%parity = parity
	       allocate( csf_new%csf(na)%occupation(1:nw+1), csf_new%csf(na)%seniority(1:nw+1),  &
	                 csf_new%csf(na)%subshellJ(1:nw+1),  csf_new%csf(na)%subshellX(1:nw+1) )
               csf_new%csf(na)%occupation(1:nw)= csf_old%csf(i)%occupation(1:nw)
               csf_new%csf(na)%seniority(1:nw) = csf_old%csf(i)%seniority(1:nw)
               csf_new%csf(na)%subshellJ(1:nw) = csf_old%csf(i)%subshellJ(1:nw)
               csf_new%csf(na)%subshellX(1:nw) = csf_old%csf(i)%subshellX(1:nw)
               !
               csf_new%csf(na)%occupation(nw+1)= 1
               csf_new%csf(na)%seniority(nw+1) = get_subshell_seniority(kappa,1,j)
               csf_new%csf(na)%subshellJ(nw+1) = j
               csf_new%csf(na)%subshellX(nw)   = csf_old%csf(i)%totalJ
               csf_new%csf(na)%subshellX(nw+1) = totalJ
	       if (present(index)) then
	          index(na) = i
	       end if
            end if
         end do
         if (rabs_use_stop   .and.   na /= nocsf_add) then
            stop "set_configuration_scheme(): program stop G."
         end if
      end if
      !
   end subroutine set_configuration_scheme
   !
end module rabs_csl
