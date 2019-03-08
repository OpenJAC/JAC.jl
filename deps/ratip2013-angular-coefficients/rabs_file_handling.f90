module rabs_file_handling
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module maintains most of the procedures to read in and to write out information to/from the (GRASP92) interface 
! files. Attempts has been made to formalize the inquiry of file names and the reading of the data.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_grasp92
   use rabs_io_dialog
   use rabs_nucleus
   implicit none
   private
   !
   public  :: file_get_csl_list
                 ! Opens a .csl file, checks it and loads the list of CSF into the internal arrays.
   public  :: file_get_eigenvectors_unformatd
                 ! Reads in the eigenvectors of a number of ASF from a standard GRASP-92 .mix mixing coefficient file.
   public  :: file_get_eigenvectors_formatted
                 ! Reads in the eigenvectors of a number of ASF from a formatted GRASP-92 .mix mixing coefficient file.
   public  :: file_get_isodat_grasp92    
                 ! Open, check, load data from and close the  .iso  file.
   public  :: file_get_mix
		 ! Opens a .mix GRASP92 mixing file and reads in the levels and mixing coefficients for the given states.
   public  :: file_get_rwf
                 ! Opens a .out file and reads in the radial orbitals for the bound states.
   public  :: file_get_rwf_det
                 ! Opens a .out file and reads in the radial orbitals for the bound states within a determinant basis.
   public  :: file_get_xpn
                 ! Opens a .xpn file and reads in the determinant expansion for some bound states.
   private :: file_load_xpn
		 ! Loads all necessary input from an .xpn file into an ASF determinant basis structure.
   public  :: file_read_isodat_grasp92  
                 ! Reads the nuclear parameters from a grasp92 .iso file.
   public  :: file_open            
                 ! Open a file on a given stream.
   public  :: file_open_formatted_stream            
                 ! Open a (new) formatted file on a given stream.
   public  :: file_open_formatted_stream_old            
                 ! Open a (existing) formatted file on a given stream.
   public  :: file_open_unformatted_stream            
                 ! Open a unformatted file on a given stream.
   public  :: file_write_mix
                 ! Writes a new .mix file.
   !
   ! Define global logical flags for the control of the file_handling module; the default values for these flags may 
   ! be overwritten interactively during input time
   logical, public :: file_use_formatted_rwf      = .false.
   !
contains
   !
   !
   subroutine file_get_csl_list(wstring,csf_set,csl_file)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a .csl file, checks it format and number of functions, and loads the list of CSF into the internal arrays. 
   ! This file is always attached to stream 21.
   !
   ! Calls: load_csl_from_file(), file_open().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)   :: wstring
      type(csf_basis), intent(inout) :: csf_set
      character(len=*), optional     :: csl_file
      !
      integer                        :: ierr, ios
      character(len=15)              :: record
      character(len=256)             :: buffer, csl_list_file
      !
    1 call print_input(wstring)
      buffer = read_input()
      read (buffer, "(a)") csl_list_file
      if (len(trim(csl_list_file)) == 0) goto 1
      !
      call file_open(21,csl_list_file,"formatted  ","old",ierr)
      if (ierr == 1) goto 1
      !
      ! Check the first record of the file; if not as expected, try again
      read (21,"(1a15)",iostat = ios) record
      if (ios /= 0   .or.   record(1:15) /= "Core subshells:") then
         print *, "ios, record(1:15) = ",ios, record(1:15)
         print *, "Not a configuration symmetry list file;"
         close (21)
         goto 1
      end if
      !
      ! Load data from the  .csl  file
      call load_csl_from_file(csf_set)
      !
      ! Close the  .csl  file
      close (21)
      !
      ! Return the name of the .csl file if required
      if (present(csl_file)) then
         csl_file = trim(csl_list_file)
      end if
      !
   end subroutine file_get_csl_list
   !
   !
   subroutine file_get_eigenvectors_unformatd(asf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Read in the eigenvectors of a number of ASF in a given CSF basis. Open, check, load data from and close the  
   ! standard .mix  file from GRASP-92. This file is always attached to stream 25.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=6)                     :: g92mix
      character(len=256)                   :: grasp_mix_file
      type(asf_basis), intent(inout)       :: asf_set
      !
      integer                              :: i, j, ios, ierr, test_no_electrons, test_nocsf, test_nwshells
      integer, dimension(:), pointer       :: iatjpo, iaspar
      real(kind=dp), dimension(:), pointer :: eval
      !
      ! Open the  .mix  file of the mixing coefficients from GRASP92.
    1 print *, "Enter the name of the GRASP92 mixing coefficient file:"
      read *,  grasp_mix_file
      call file_open(25,grasp_mix_file,"unformatted","old",ierr)
      !
      if (ierr /= 0) then
    2    print *, "Enter the name of the GRASP92 mixing coefficient file:"
         read *,  grasp_mix_file
         if (len_trim(grasp_mix_file) == 0) goto 2
         goto 1
      end if
      !
      ! Check the header of the file; if not as expected, try again
      read (25,iostat=ios) g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "g92mix = ",g92mix
         print *, "Not a GRASP92 mixing coefficients file;"
         close (25)
         goto 1
      end if
      !
      ! Check consistency with the CSF basis
      read (25) test_no_electrons,test_nocsf,test_nwshells
      if (test_no_electrons /= asf_set%csf_set%number_of_electrons   .or. test_nocsf  /= asf_set%csf_set%nocsf   .or. &
          test_nwshells     /= asf_set%csf_set%nwshells) then
         print *, "This mixing coefficients file is not appropriate to the CSF basis"
         print *, " as loaded from the configuration symmetry list file."
         close (25)
         goto 1
      endif
      !
      call save_input(grasp_mix_file,.false.)
      !
      ! Allocate memory for the set of ASF and load the data from the .mix file
      print *, "Loading mixing coefficients file ..."
      read (25) asf_set%noasf
      allocate( asf_set%asf(1:asf_set%noasf) )
      do  i = 1,asf_set%noasf
         allocate( asf_set%asf(i)%eigenvector(1:asf_set%csf_set%nocsf) )
      end do
      allocate( iatjpo(1:asf_set%noasf), iaspar(1:asf_set%noasf), &
                eval(1:asf_set%noasf) ) 
      !
      read (25) (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
      read (25) (iatjpo(i),iaspar(i),i = 1,asf_set%noasf)
      do  i = 1,asf_set%noasf
         asf_set%asf(i)%totalJ = iatjpo(i) - 1
         if (iaspar(i) == 1) then
            asf_set%asf(i)%parity = "+"
         else if (iaspar(i) == -1) then
            asf_set%asf(i)%parity = "-"
         else if (rabs_use_stop) then
            stop "set_eigenvectors_unformatted(): program stop A."
         end if
      end do
      !
      read (25) asf_set%average_energy,(eval(i),i = 1,asf_set%noasf)
      do  i = 1,asf_set%noasf
         asf_set%asf(i)%energy = asf_set%average_energy + eval(i)
      end do
      !
      read (25) ((asf_set%asf(i)%eigenvector(j),j=1,asf_set%csf_set%nocsf), i=1,asf_set%noasf)
      !
      deallocate( iatjpo, iaspar, eval ) 
      print *, ' ... load complete;'
      close (25)
      !
   end subroutine file_get_eigenvectors_unformatd
   !
   !
   subroutine file_get_eigenvectors_formatted(asf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Read in the eigenvectors of a number of ASF in a given CSF basis. Open, check, load data from and close the  
   ! formatted .mix  file a generated from an unformatted GRASP92 .mix file. This file is always attached to stream 25.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=6)                     :: g92mix
      character(len=256)                   :: grasp_mix_file
      type(asf_basis), intent(inout)       :: asf_set
      !
      integer                              :: i, j, ios, ierr, test_no_electrons, test_nocsf, test_nwshells
      integer, dimension(:), pointer       :: iatjpo, iaspar
      real(kind=dp), dimension(:), pointer :: eval
      !
      ! Open the  .mix  file of the mixing coefficients from GRASP92.
      print *, "Enter the name of the GRASP92 mixing coefficient file:"
      read *,  grasp_mix_file
    1 call file_open(25,grasp_mix_file,"formatted  ","old",ierr)
      !
      if (ierr /= 0) then
    2    print *, "Enter the name of the GRASP92 mixing coefficient file:"
         read *,  grasp_mix_file
         if (len_trim(grasp_mix_file) == 0) goto 2
         goto 1
      end if
      !
      ! Check the header of the file; if not as expected, try again
      read (25,"(a)",iostat=ios) g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 mixing coefficients file;"
         close (25)
         goto 1
      end if
      !
      ! Check consistency with the CSF basis
      read (25,"(3i6)") test_no_electrons,test_nocsf,test_nwshells
      if (test_no_electrons /= asf_set%csf_set%number_of_electrons   .or.  test_nocsf  /= asf_set%csf_set%nocsf   .or. &
          test_nwshells     /= asf_set%csf_set%nwshells) then
         print *, "This mixing coefficients file is not appropriate to the CSF basis"
         print *, " as loaded from the configuration symmetry list file."
         close (25)
         goto 1
      endif
      !
      call save_input(grasp_mix_file,.false.)
      !
      ! Allocate memory for the set of ASF and load the data from the .mix file
      print *, "Loading mixing coefficients file ..."
      read (25,"(i6)") asf_set%noasf
      allocate( asf_set%asf(1:asf_set%noasf) )
      do  i = 1,asf_set%noasf
         allocate( asf_set%asf(i)%eigenvector(1:asf_set%csf_set%nocsf) )
      end do
      allocate( iatjpo(1:asf_set%noasf), iaspar(1:asf_set%noasf), eval(1:asf_set%noasf) ) 
      !
      read (25,"(100i6)")      (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
      read (25,"(100(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, i = 1,asf_set%noasf)
      !
      read (25,"(101e16.9)") asf_set%average_energy, (eval(i),i = 1,asf_set%noasf)
      do  i = 1,asf_set%noasf
         asf_set%asf(i)%energy = asf_set%average_energy + eval(i)
      end do
      !
      do  j = 1,asf_set%csf_set%nocsf
         read (25,"(100e16.9)") (asf_set%asf(i)%eigenvector(j), i=1,asf_set%noasf)
      end do
      !
      deallocate( iatjpo, iaspar, eval ) 
      print *, ' ... load complete;'
      close (25)
      !
   end subroutine file_get_eigenvectors_formatted
   !
   !
   subroutine file_get_isodat_grasp92()
   !--------------------------------------------------------------------------------------------------------------------
   ! Open, check, load data from and close the  .iso  file. This file is always attached to stream  22.
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      integer            :: ierr, ios
      logical            :: found
      character(len=3)   :: file_status
      character(len=11)  :: file_format, file_default
      character(len=14)  :: record
      character(len=256) :: buffer, file_name
      !
      ! The  .iso  file is formatted; it must exist
      file_default = "grasp92.iso"
      file_format  = "formatted"
      file_status  = "old"
      !
      ! Look for  grasp92.iso
      inquire (file = file_default, exist = found)
      !
    1 if (found) then
         !
         ! File  grasp92.iso  exists; ascertain that it is to be used; if it is not to be used, determine another filename
         call print_input("file  grasp92.iso  found; enter another file name if this file")
         call print_input(" is not to be used as the isotope data file; <cr> otherwise:")
         buffer = read_input()
         read (buffer, "(a)") file_name
         if (len_trim (file_name) == 0) then
            file_name = file_default
         end if
      else
         !
         ! File  grasp92.iso  does not exist; determine the name of the .iso  file
         call print_input("Enter the name of the isotope data file:")
         buffer = read_input()
         read (buffer, "(a)") file_name
         if (len_trim (file_name) == 0) then
            found = .false.
            goto 1
         end if
      endif
      !
      call file_open (22,file_name,file_format,file_status,ierr)
      if (ierr .eq. 1) then
         found = .false.
         goto 1
      end if
      !
      ! Check the first record of the file; if not as expected, try again
      read (22,"(a)",iostat = ios) record
      if ((ios .ne. 0)  .or.  (record(1:14) .ne. 'Atomic number:')) then
         call print_input("not an isotope data file;")
         close (22)
         found = .false.
         goto 1
      endif
      !
      ! load data from the  .iso  file
      call file_read_isodat_grasp92()
      !
      ! close the  .iso  file
      close (22)
      !
   end subroutine file_get_isodat_grasp92
   !
   !
   subroutine file_get_mix(wstring,asf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a final-state .mix Mixing Coefficient File from GRASP92 and reads in all required levels and mixing coefficients.
   !
   ! Calls: load_mix_file(), file_open().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)   :: wstring
      type(asf_basis), intent(inout) :: asf_set
      !
      integer		             :: ierr, ios
      character(len=6)               :: g92mix
      character(len=256)             :: buffer, mix_file
      !
    2 call print_input(wstring)
      buffer = read_input()
      read (buffer, "(a)") mix_file
      if (len(trim(mix_file)) == 0) goto 2
      !
      call file_open(25,mix_file,"unformatted","old",ierr)
      if (ierr /= 0) goto 2
      !
      ! Check the header of the file; if not as expected for unformatted
      ! files, check formatted form
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
	 close (25)
	 goto 3
      else
	 call load_mix_file_grasp92(asf_set,.false.,ierr)
	 if (ierr /= 0) then
            call print_input("Not a proper .mix mixing coefficient file for the given .csl list; reenter ...")
	    close (25)
	    goto 2
	 end if
	 goto 4
      end if
      !
      ! Try formatted file format; check the header of the file; if not as expected for formatted files, try again
    3 call file_open(25,mix_file,"formatted  ","old",ierr)
      if (ierr /= 0) goto 2
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
	 print *, "ios, g92mix = ",ios, g92mix
	 print *, "Not a GRASP92 Mixing Coefficients File;"
	 close (25)
	 goto 2
      else
	 call load_mix_file_grasp92(asf_set,.true.,ierr)
	 if (ierr /= 0) then
	    close (25)
	    goto 2
	 end if
      end if
      !
      ! Close the  .mix  file
    4 close (25)
      !
   end subroutine file_get_mix
   !
   !
   subroutine file_get_rwf(wstring,asf_bound,wave_bound,is_allocated, notall_needed)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a bound-state .rwf Radial WaveFunction File from GRASP92 and reads in all required orbital wave functions. 
   !
   ! Calls: load_rwf_file_grasp92(), file_open().
    !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)         :: wstring
      type(asf_basis), intent(inout)       :: asf_bound
      type(grasp92_orbital), intent(inout) :: wave_bound
      logical, intent(in)                  :: is_allocated
      logical, intent(in), optional        :: notall_needed   
      !
      integer                              :: iorb, ierr, ios
      character(len=31)                    :: record
      character(len=256)                   :: buffer, file_rwf
      !
      print *, "aa"
      !
    1 call print_input(wstring)
      buffer = read_input()
      read (buffer, "(a)") file_rwf
      if (len(trim(file_rwf)) == 0) goto 1
      !
      call file_open(21,file_rwf,"unformatted","old",ierr)
      if (ierr == 1) goto 1
      !
      read (21,iostat = ios) record(1:6)
      if (record(1:6) /= 'G92RWF') then
         close (21)
         !
         call file_open(21,file_rwf,"formatted  ","old",ierr)
         if (ierr == 1) goto 1
         read (21,"(a31)",iostat = ios) record
         if (ios /= 0   .or.  record(1:31) /= "G92RWF (formatted file version)") then
            print *, "ios, record(1:31) = ",ios, record(1:31)
            print *, "Not a G92RWF Radial WaveFunction File;"
            close (21)
            goto 1
         end if
         file_use_formatted_rwf = .true.
      else
         file_use_formatted_rwf = .false.
      endif
      !
      if (.not.is_allocated) then
         !
         ! Allocate memory for the bound-state orbitals
         wave_bound%number_of_rwf = asf_bound%csf_set%nwshells
         allocate( wave_bound%rwf(1:wave_bound%number_of_rwf) )
         do  iorb = 1,asf_bound%csf_set%nwshells
            wave_bound%rwf(iorb)%orbital%n     = asf_bound%csf_set%subshell(iorb)%n
            wave_bound%rwf(iorb)%orbital%kappa = asf_bound%csf_set%subshell(iorb)%kappa
            wave_bound%rwf(iorb)%mtp           = 0
            wave_bound%rwf(iorb)%energy        = zero
            wave_bound%rwf(iorb)%gamma         = zero
            wave_bound%rwf(iorb)%pz            = zero
            wave_bound%rwf(iorb)%phase         = zero
         end do
      end if
      !
      ! Load data from the  .rwf  file
      if (present(notall_needed)  .and.  notall_needed) then
         call load_rwf_file_grasp92(wave_bound,file_use_formatted_rwf,ierr, notall_needed)
      else
         call load_rwf_file_grasp92(wave_bound,file_use_formatted_rwf,ierr)
      end if
      if (ierr /= 0) then
         goto 1
      end if
      !
      ! Close the  .rwf  file
      close (21)
      !
   end subroutine file_get_rwf
   !
   !
   subroutine file_get_rwf_det(wstring,wave_bound)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a .rwf Radial WaveFunction File from GRASP92 and reads in all required orbital wave functions. 
   !
   ! Calls: load_rwf_file_grasp92(), file_open().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)         :: wstring
      type(grasp92_orbital), intent(inout) :: wave_bound
      !
      integer		                   :: ierr, ios
      character(len=31)                    :: record
      character(len=256)                   :: buffer, file_rwf
      !
    1 call print_input(wstring)
      buffer = read_input()
      read (buffer, "(a)") file_rwf
      !!x print *, "Enter the name of the final-state Radial WaveFunction File:"
      !!x read (*,"(a)") file_rwf
      if (len(trim(file_rwf)) == 0) goto 1
      !
      call file_open(21,file_rwf,"unformatted","old",ierr)
      if (ierr == 1) goto 1
      !
      read (21,iostat = ios) record(1:6)
      if (record(1:6) /= 'G92RWF') then
	 close (21)
	 !
	 call file_open(21,file_rwf,"formatted  ","old",ierr)
	 if (ierr == 1) goto 1
	 read (21,"(a31)",iostat = ios) record
	 if (ios /= 0	.or.   &
	     record(1:31) /= "G92RWF (formatted file version)") then
	     print *, "ios, record(1:31) = ",ios, record(1:31)
	     print *, "Not a G92RWF Radial WaveFunction File;"
	     close (21)
	     goto 1
	 end if
	 file_use_formatted_rwf = .true.
      else
	 file_use_formatted_rwf = .false.
      endif
      !
      ! Load data from the  .rwf  file
      call load_rwf_file_grasp92(wave_bound,file_use_formatted_rwf,ierr)
      if (ierr /= 0) then
	 goto 1
      end if
      !
      call save_input(file_rwf,.false.)
      !
      ! Close the  .rwf  file
      close (21)
      !
   end subroutine file_get_rwf_det
   !
   !
   subroutine file_get_xpn(wstring,asf_load)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a .xpn file for some (bound) states and reads in all necessary data from this file. It also close the file 
   ! finally.
   !
   ! Calls: open_file(), reos_load_xpn_file().
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)         :: wstring
      type(asf_det_basis), intent(inout)   :: asf_load
      !
      integer		                   :: ierr, ios
      character(len=26)                    :: record
      character(len=256)                   :: buffer, file_xpn
      !
    1 call print_input(wstring)
      buffer = read_input()
      read (buffer, "(a)") file_xpn
      !!x print *, wstring
      !!x read (*,"(a)") file_xpn
      if (len(trim(file_xpn)) == 0) goto 1
      !
      call file_open(21,file_xpn,"formatted  ","old",ierr)
      if (ierr == 1) goto 1
      !
      ! Check the first record of the file; if not as expected, try again
      read (21,"(a26)",iostat = ios) record
      if (ios /= 0   .or.   record(1:26) /= "ASF-based CESD output file") then
	 print *, "ios, record(1:26) = ",ios, record(1:26)
	 print *, "Not a ASF-based CESD output file;"
	 close (21)
	 goto 1
      end if
      !
      ! Load data from the  .xpn  file
      call file_load_xpn(asf_load)
      !
      call save_input(file_xpn,.false.)
      !
      ! Close the  .xpn  file
      close (21)
      !
   end subroutine file_get_xpn
   !
   !
   subroutine file_load_xpn(asf_load)
   !--------------------------------------------------------------------------------------------------------------------
   ! Loads all necessary input from a CESD .xpn output file on stream 21 into a ASF determinant basis structure.
   !
   ! Calls: unpack_occupation_from_integer().
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(asf_det_basis), intent(inout) :: asf_load
      !
      integer    :: i, iorb,iorbtest, iblock, idet, idettest, imax, isu, iso, itest, nblock, nobit, noint, nn
      integer, dimension(:), allocatable :: occupation
      !
      read (21,*)
      read (21,*)
      read (21,*)
      read (21,*)
      read (21, "(i6)") asf_load%noasf
      allocate( asf_load%asf(1:asf_load%noasf) )
      !
      ! First read the orbital reference list and then the list of determinants in terms of their occupation numbers 
      ! with respect to the reference list
      !
      read (21, "(i6)") asf_load%det_set%norbital
      allocate( asf_load%det_set%orbital(1:asf_load%det_set%norbital) )
      !
      do  iorb = 1,asf_load%det_set%norbital,6
	 if (iorb + 5 > asf_load%det_set%norbital) then
	    imax = asf_load%det_set%norbital
	 else
	    imax = iorb + 5
	 end if
	 !
	 read (21, "(6(i4,1x,3i4,4x))") (iorbtest, asf_load%det_set%orbital(i)%n, asf_load%det_set%orbital(i)%kappa,      &
	                                 asf_load%det_set%orbital(i)%mm,i=iorb,imax)
	 if (rabs_use_stop   .and.   imax /= iorbtest) then
	    stop "reos_load_xpn_file(): program stop A."
	 end if
      end do
      !
      read (21,*)
      read (21,*)   ! List of determinants:
      read (21,*)   ! ---------------------
      read (21, "(i6)") asf_load%det_set%nod
      asf_load%det_set%nodmax = asf_load%det_set%nod
      !
      !
      ! Determine the number of bits in a standard integer, the number of
      ! integer cells for a single determinant, and allocate memory
      nobit = bit_size(itest)
      if (mod(asf_load%det_set%norbital,nobit) == 0 ) then
	 asf_load%det_set%noint = asf_load%det_set%norbital/nobit
      else
	 asf_load%det_set%noint = asf_load%det_set%norbital/nobit + 1
      end if
      !
      noint = asf_load%det_set%noint
      allocate( asf_load%det_set%determinant(1:asf_load%det_set%nod) )
      do  i = 1,asf_load%noasf
	 allocate( asf_load%asf(i)%eigenvector(1:asf_load%det_set%nod) )
      end do
      print *, "nobit, noint = ", nobit, noint
      !
      do  idet = 1,asf_load%det_set%nod
         print *, "idet = ", idet
	 allocate( asf_load%det_set%determinant(idet)%occupation(1:noint) )
	 if (nobit <= 32) then
	    read (21, "(i7,3x,i4,a1,2x,40(i12))") idettest, asf_load%det_set%determinant(idet)%totalM,  &
		                                            asf_load%det_set%determinant(idet)%parity,  &
	                                                   (asf_load%det_set%determinant(idet)%occupation(i),i=1,noint)
	 else
	    read (21,"(i7,3x,i4,a1,2x,40(i24))") idettest,  asf_load%det_set%determinant(idet)%totalM,  &
		                                            asf_load%det_set%determinant(idet)%parity,  &
	                                                   (asf_load%det_set%determinant(idet)%occupation(i),i=1,noint)
	 end if
	 if (rabs_use_stop   .and.   idet /= idettest) then
	    stop "reos_load_xpn_file(): program stop B."
	 end if
      end do
      !
      ! Determine number of electrons
      nn = noint * nobit
      allocate( occupation(nn) )
      call unpack_occupation_from_integer(asf_load%det_set%determinant(1),occupation,nn)
      asf_load%det_set%number_of_electrons = sum(occupation(1:nn))
      deallocate( occupation )
      !
      read (21,*)
      read (21,*)   ! CESD expansion of the ASF:
      read (21,*)   ! --------------------------
      read (21,*)
      !
      read (21, "(2x,e26.19)") asf_load%average_energy
      !
      ! Read the eigenvectors in 'blocks' of maximal 13 eigenvectors
      if (mod(asf_load%noasf,13) == 0) then
	 nblock = asf_load%noasf / 13
      else
	 nblock = asf_load%noasf / 13 + 1
      end if
      !
      do  iblock = 1,nblock
	 isu = (iblock-1)*13 + 1
	 iso = min( ((iblock-1)*13 + 13),asf_load%noasf )
	 read (21,*)
	 read (21,*) 
	 read (21, "(15(6x,i5,6x))")	(asf_load%asf(i)%level_No,i=isu,iso)
	 read (21, "(15(2x,e26.19))")   (asf_load%asf(i)%energy,  i=isu,iso)
	 read (21, "(15(6x,i5,6x))")	(asf_load%asf(i)%totalJ,   i=isu,iso)
	 read (21,*)   
	 read (21, "(15(6x,4x,a1,6x))") (asf_load%asf(i)%parity,i=isu,iso)
	 do  idet = 1,asf_load%det_set%nod
	    read (21, "(15(2x,es15.8))")  &
	       (asf_load%asf(i)%eigenvector(idet),i=isu,iso)
	 end do
      enddo
      !
   end subroutine file_load_xpn
   !
   !
   subroutine file_read_isodat_grasp92()
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads the data from the  grasp92 .iso  file. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      real(kind=dp) :: aparm, cparm
      !
      print *, 'loading isotope data file ...'
      !
      ! Read and echo pertinent information from  .iso  file; atomic number and nuclear geometry
      read (22,*) nuclear_charge
      read (22,*)
      read (22,*) atomic_mass
      read (22,*)
      read (22,*) aparm
      read (22,*)
      read (22,*) cparm
      !
      if (atomic_mass .ne. zero) then
         nuclear_model = "fermi"
         fermi_c_parameter = cparm * convert_fermi_to_bohr
         fermi_a_parameter = aparm * convert_fermi_to_bohr
      else
         nuclear_model = "point"
      endif
      !
      ! Nuclear mass
      read (22,*)
      read (22,*) nuclear_mass
      !
      ! Nuclear spin and moments
      read (22,*)
      read (22,*) nuclear_spin
      read (22,*)
      read (22,*) nuclear_dipole_moment
      read (22,*)
      read (22,*) nuclear_quadrupole_moment
      !
      print *, ' ... load complete;'
      !
   end subroutine file_read_isodat_grasp92
   !
   !
   subroutine file_open(nfile,file_name,file_format,file_status,ierr)
   !--------------------------------------------------------------------------------------------------------------------
   ! Issues open for file with unit number nfile, name  file_name, format  file_format, status file_status.  If this is 
   ! successful the head is positioned to the beginning of the file and ierr is 0; otherwise ierr is set to 1. This 
   ! subroutine was originally written by F A Parpia and adapted to rabs.
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer,            intent(in) :: nfile
      character(len=3),   intent(in) :: file_status 
      character(len=11),  intent(in) :: file_format
      character(len=256), intent(in) :: file_name
      !
      integer  :: loc, ios, ierr
      !
      open (nfile, file = file_name, form = file_format, status = file_status, iostat = ios)
      !
      if (ios == 0) then
         rewind (nfile)
         ierr = 0
      else
         loc = len_trim (file_name)
         print *, "file_open(): error opening file ",file_name(1:loc), " as ",file_status,";"
         ierr = 1
      end if
      !
   end subroutine file_open
   !
   !
   subroutine file_open_formatted_stream(stream,wstring,formatted_stream)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a formatted file on stream
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)          :: stream
      character(len=*), intent(in) :: wstring
      character(len=*), optional   :: formatted_stream
      !
      integer                      :: ierr
      character(len=256)           :: buffer, formatted_file
      !
    1 call print_input(wstring)
      buffer = read_input()
      read (buffer, *)  formatted_file
      call file_open(stream,formatted_file,"formatted  ","new",ierr)
      !
      if (ierr /= 0) goto 1
      !
      ! Return the name of the file if required
      if (present(formatted_stream)) then
         formatted_stream = trim(formatted_file)
      end if
      !
   end subroutine file_open_formatted_stream
   !
   !
   subroutine file_open_formatted_stream_old(stream,wstring,formatted_stream)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a formatted file on stream
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)          :: stream
      character(len=*), intent(in) :: wstring
      character(len=*), optional   :: formatted_stream
      !
      integer :: ierr
      character(len=256)           :: formatted_file
      !
    1 print *, wstring
      read *,  formatted_file
      call file_open(stream,formatted_file,"formatted  ","old",ierr)
      !
      if (ierr /= 0) goto 1
      !
      ! Return the name of the file if required
      if (present(formatted_stream)) then
         formatted_stream = trim(formatted_file)
      end if
      call save_input(formatted_file,.false.)
      !
   end subroutine file_open_formatted_stream_old
   !
   !
   subroutine file_open_unformatted_stream(stream,wstring,unformatted_stream)
   !--------------------------------------------------------------------------------------------------------------------
   ! Opens a unformatted file on stream
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)          :: stream
      character(len=*), intent(in) :: wstring
      character(len=*), optional   :: unformatted_stream
      !
      integer                      :: ierr
      character(len=256)           :: unformatted_file
      !
    1 print *, wstring
      read *,  unformatted_file
      call file_open(stream,unformatted_file,"unformatted","new",ierr)
      !
      if (ierr /= 0) goto 1
      !
      ! Return the name of the file if required
      if (present(unformatted_stream)) then
         unformatted_stream = trim(unformatted_file)
      end if
      !
      call save_input(unformatted_file,.false.)
      !
   end subroutine file_open_unformatted_stream
   !
   !
   subroutine file_write_mix(asf_set,levels,number_of_levels)
   !--------------------------------------------------------------------------------------------------------------------
   ! Writes out a .mix file.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)               :: number_of_levels
      integer, dimension(:), intent(in) :: levels
      type(asf_basis), intent(inout)    :: asf_set
      !
      integer                           :: i, j, ierr, ios, noasf_sofar, noasf
      character(len=6)                  :: g92mix = "G92MIX"
      character(len=256)                :: file_mix_new
      !
    1 print *, "Enter the name of the .mix  mixing coefficient file which has to be created:"       
      read (*,"(a)") file_mix_new
      call file_open(25,file_mix_new,"formatted  ","new",ierr)
      if (ierr /= 0  .or.  len_trim(file_mix_new) == 0) goto 1
      !
      ! Write out a file header
      write (25,"(a)") g92mix//" (formatted file version)."
      !
      print *, "Write formatted Mixing Coefficients File ..."
      print *, "   Total number of ASFs = ",number_of_levels
      !
      write (25,"(3i6)") asf_set%csf_set%number_of_electrons, asf_set%csf_set%nocsf, asf_set%csf_set%nwshells
      write (25,"(i6)")  number_of_levels    
     
      noasf = min(60,number_of_levels)
      write (25,"(60i6)")      (asf_set%asf(levels(i))%level_No,i=1,noasf)
      write (25,"(60(i5,a1))") (asf_set%asf(levels(i))%totalJ, asf_set%asf(levels(i))%parity, i=1,noasf)
      !
      write (25,"(61e26.19)") asf_set%average_energy, (asf_set%asf(levels(i))%energy, i=1,noasf)
      do j = 1, asf_set%csf_set%nocsf
         write (25,"(60e16.9)") (asf_set%asf(levels(i))%eigenvector(j), i=1,noasf)
      end do
      !  
      if (number_of_levels <= 60) goto 10
      !
      noasf_sofar = 0
    3 noasf_sofar = noasf_sofar + 60
      !
      noasf = min(noasf_sofar+60,number_of_levels)
      write (25,"(60i6)")      (asf_set%asf(levels(i))%level_No, i = noasf_sofar+1,noasf)
      write (25,"(60(i5,a1))") (asf_set%asf(levels(i))%totalJ, asf_set%asf(levels(i))%parity, i = noasf_sofar+1,noasf)
      !
      write (25,"(60e26.19)")  (asf_set%asf(levels(i))%energy, i = noasf_sofar+1,noasf)
      do  j = 1,asf_set%csf_set%nocsf
	 write (25,"(60e16.9)") (asf_set%asf(levels(i))%eigenvector(j), i=noasf_sofar+1,noasf)
      end do
      !
      if (number_of_levels > noasf_sofar+60) goto 3
      !
   10 print *, "... write out complete;"
      !
      call save_input(file_mix_new,.false.)
      close (25)
      !
   end subroutine file_write_mix
   !
end module rabs_file_handling
