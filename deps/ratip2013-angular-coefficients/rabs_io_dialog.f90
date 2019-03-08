module rabs_io_dialog
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
   !
   use iso_fortran_env   
   use rabs_constant
   use rabs_functions_math 
   use rabs_functions_string   
   use rabs_nucleus
   implicit none
   !
   public  :: get_yesno
                 ! Reads and interpretes 'y' or 'n' from the default input unit.
   public  :: io_energy_unit
                 ! Reads in and determines the energy unit for input/output. 
   public  :: io_grid_parameters
                 ! Reads in and determines the parameters for the radial grid.
   public  :: io_levels
                 ! Reads in an selected set of levels numbers.
   public  :: io_sharings
                 ! Reads in and determine a number of `sharings' between 0..1.
   public  :: io_time_unit
                 ! Reads in and determines the time unit for input/output. 
   public  :: io_transition_multipoles
                 ! Reads in the transition multipoles from a given record. 
   public  :: io_transition_pairs
                 ! Reads in and selects pairs of transitions.
   public  :: io_transition_triples
                 ! Reads in and selects triples of transitions.
   public  :: print_input
                 ! Prints some records to the standard stream (*) but keeps a `copy' for the input buffer.
   public  :: print_runtime_component
                 ! Prints the date, time and compiler (options) when the given component runs
   public  :: read_input
                 ! Reads a record from standard stream (*), stores it into the input buffer and returns it as string
                 ! of variable length for internal file reading.
   public  :: set_warning
                 ! Writes a warning to the warning buffer.
   public  :: write_input
                 ! Writes all input records to stream with some proper begin/end line.
   public  :: write_warnings
                 ! Writes all warnings records to stream with some proper begin/end line.
   !
   !
   ! Define some global data for the input/output of the RATIP code and for collecting warnings during its execution
   type, private :: io_buffer
      character(len=130), dimension(200) :: input
      character(len=130), dimension(100) :: warnings
   end type io_buffer
   !
   integer, private         :: no_of_input = 0,  no_of_warnings = 0
   type(io_buffer), private :: buffer
   !
   integer, dimension(200), public :: select_level_i, select_level_f, select_level_m
   !
contains
   !
   function get_yesno()   result(yes)
   !--------------------------------------------------------------------------------------------------------------------
   ! This function reads and interpretes 'y' or 'n' from the default input unit; it returns .true.  if 'y' is entered 
   ! and  .false.  if 'n' is entered.  
   !--------------------------------------------------------------------------------------------------------------------
      !
      logical          :: yes
      integer          :: ios
      character(len=1) :: response
      !
    2 read (unit=*,fmt="(a)",iostat=ios) response
      if (ios /= 0) then
         print *, "Expecting <y><cr> or <n><cr> ..."
         read (unit= *,fmt="(a)",iostat=ios) response
         if (ios /= 0) then
            print *, "Expecting <y><cr> or <n><cr> ..."
            read (unit=*,fmt="(a)",iostat=ios) response
         endif
      end if
      !
      if (response == "y") then
         yes = .true.
      else if (response == "n") then
         yes = .false.
      else
         print *, "Expecting <y><cr> or <n><cr> ..."
         goto 2
      end if
      !
      !!x call print_input(response)
      !!x !
      no_of_input               = no_of_input + 1
      buffer%input(no_of_input) = response
      !
   end function get_yesno
   !
   !
   subroutine io_energy_unit()
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in and determines the energy unit for input/output.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=256) :: buffer, record
      !
      ! Determine the units for the printout and further optional input
    1 call print_input("Which units are to be used to enter and to print the energies of the continuum orbitals ?")
      call print_input("    A       : Angstrom;")
      call print_input("    eV      : electron volts;")
      call print_input("    Hartree : Hartree atomic units;")
      call print_input("    Hz      : Hertz;")
      call print_input("    Kayser  : [cm**(-1)];")
      !!x print *, "Which units are to be used to enter and to print the energies of the continuum orbitals ?"
      !!x print *, "    A       : Angstrom;"
      !!x print *, "    eV      : electron volts;"
      !!x print *, "    Hartree : Hartree atomic units;"
      !!x print *, "    Hz      : Hertz;"
      !!x print *, "    Kayser  : [cm**(-1)];"
      buffer = read_input()
      read (buffer, "(a)") record
      if (len_trim(record) > 0) then
         record = adjustl(record)
         select case(record(1:7))
         case("A      ", "eV     ", "Hartree", "Hz     ", "Kayser ")
            energy_unit = record(1:7)      
         case default
            !!x print *, "Unable to decode the units '"//record(1:7)// "' to be used for the transition energies; reenter ..."
            call print_input("Unable to decode the units '"//record(1:7)// &
                             "' to be used for the transition energies; reenter ...")
            goto 1
         end select
      else 
         energy_unit = "eV     "
      end if
      energy_inverse = .false.
      !
      !!x call save_input(record,.false.)
      !!x !
      ! Determine conversion factor for calculating energies in this unit
      select case(energy_unit)
      case("A      ")
         energy_factor  = 1.0e8_dp / convert_au_to_kaysers
         energy_inverse = .true.
      case("eV     ")
         energy_factor  = convert_au_to_ev
      case("Hartree")
         energy_factor  = one
      case("Hz     ")
         energy_factor  = convert_au_to_per_sec
      case("Kayser ")
         energy_factor  = convert_au_to_kaysers
      end select
      !
   end subroutine io_energy_unit
   !
   !
   subroutine io_grid_parameters(keystring)
   !--------------------------------------------------------------------------------------------------------------------
   ! Determines of reads in the new parameters for specifying the radial grid.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in) :: keystring
      logical                      :: yes
      character(len=256)           :: buffer, record
      !
      select case(keystring)
      case("standard")
         !
         if (nuclear_model == "point") then
	    rnt_grasp92 = exp(-65.0_dp/16.0_dp) / nuclear_charge
	    h_grasp92   = half**4
	    n_grasp92   = 220
         else
	    rnt_grasp92 = 2.0e-6_dp
	    h_grasp92   = 5.0e-2_dp
	    n_grasp92   = 390
         endif
	 !
      case("modify")
	 !
         !!x print *, "The default radial grid parameters for this case are:"
         !!x print *, " rnt = ",rnt_grasp92,";"
         !!x print *, " h   = ",h_grasp92,  ";"
         !!x print *, " hp  = ",hp_grasp92, ";"
         !!x print *, " n   = ",n_grasp92,  ";"
         !!x print *, " revise these values ?"
         call print_input("The default radial grid parameters for this case are:")
         write(buffer, "(a, g0.3, a)" ) " rnt = ",rnt_grasp92,";"
         call print_input(buffer(1:20))
         write(buffer, "(a, g0.3, a)" ) " h   = ",h_grasp92,  ";"
         call print_input(buffer(1:20))
         write(buffer, "(a, g0.3, a)" ) " hp  = ",hp_grasp92, ";"
         call print_input(buffer(1:20))
         write(buffer, "(a, g0.3, a)" ) " n   = ",n_grasp92,  ";"
         call print_input(buffer(1:20))
         call print_input(" revise these values ?")
         !!x yes = get_yes_stream()
         yes = get_yesno()
         if (yes) then
            !!x print *, "Enter rnt:"
            !!x read *, rnt_grasp92
            call print_input("Enter rnt:")
            buffer = read_input()
            read(buffer, *) rnt_grasp92
	    !
	    !!x write(stream_input,"(es12.5)") rnt_grasp92
	    !!x backspace (stream_input)
	    !!x read (stream_input,*) record
            !!x call save_input(trim(record)//" :: rnt_grasp92",.false.)
	    !!x !
            !!x print *, "Enter h:"
            !!x read *, h_grasp92
            call print_input("Enter h:")
            buffer = read_input()
            read(buffer, *) h_grasp92
	    !
	    !!x write(stream_input,"(es12.5)") h_grasp92
	    !!x backspace (stream_input)
	    !!x read (stream_input,*) record
            !!x call save_input(trim(record)//" :: h_grasp92",.false.)
	    !!x !
            !!x print *, "enter hp:"
            !!x read *, hp_grasp92
            call print_input("Enter hp:")
            buffer = read_input()
            read(buffer, *) hp_grasp92
	    !
	    !!x write(stream_input,"(es12.5)") hp_grasp92
	    !!x backspace (stream_input)
	    !!x read (stream_input,*) record
            !!x call save_input(trim(record)//" :: hp_grasp92",.false.)
	    !!x !
            !!x print *, "enter n:"
            !!x read *, n_grasp92
            call print_input("Enter n:")
            buffer = read_input()
            read(buffer, *) n_grasp92
	    !
	    !!x write(stream_input,"(i6)") n_grasp92
	    !!x backspace (stream_input)
	    !!x read (stream_input,*) record
            !!x call save_input(trim(record)//"        :: n_grasp92",.false.)
	    !!x !
         end if
	 !
     case default
         stop "input_grid_parameters(): program stop A."
      end select
      !
   end subroutine io_grid_parameters
   !
   !
   subroutine io_levels(number_of_levels)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in some selected level numbers.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(inout) :: number_of_levels
      !
      integer            :: level_i, blank_position
      logical            :: yes
      character(len=256) :: buffer, record
      character(len=20 ) :: string
      !
      number_of_levels  = 0
      !
      !!x print *, "Select individual levels ?"
      call print_input("Select individual levels ?")
      yes = get_yes_stream()
      if (yes) then
    4    call print_input("Enter a list of levels numbers:")
         !!x print *, "Enter a list of levels numbers:"
         buffer = read_input()
         read (buffer, "(a)") record
         !!x call save_input(record,.false.)
	 record         = adjustl(record)
         !
    5    blank_position = scan(record," ")
	 string         = adjustl(record(1:blank_position-1))
	 number_of_levels  = number_of_levels + 1
	 select_level_i(number_of_levels) = get_integer_from_string(string)
	 record         = adjustl(record(blank_position+1:256))
         if (len_trim(record) > 0)  goto 5
      end if
      !
      if (number_of_levels == 0) then
         !!x print *, "No individual level number has been selected; the program continues the default path."
         call print_input("No individual level number has been selected; the program continues the default path.")
      end if
      !
   end subroutine io_levels
   !
   !
   subroutine io_sharings(string,number_of_sharings,sharing,weight,is_gl)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in and determines a number of sharings in the interval 0..1, for instance, for integration.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)               :: string             
      integer, intent(inout)                     :: number_of_sharings
      real(kind=dp), dimension(:), intent(inout) :: sharing,weight
      logical, intent(out)                       :: is_gl
      !
      integer                       :: i
      character(len=256)            :: buffer, record
      logical                       :: yes
      !!x real(kind=dp), dimension(500) :: wa
      !
    1 call print_input(string)
      call print_input(" sharings can be given explicitly or for Gauss-Legendre (GL) integration; use GL ?")
      yes = get_yes_stream()
      if (yes) then
         call print_input("Enter the order of the GL integration:") 
         buffer = read_input()
         read (buffer, *) number_of_sharings
         !
         call set_gauss_zeros(zero,one,number_of_sharings,sharing,weight)
         is_gl = .true.
      else
         call print_input("Enter the number of sharings:") 
         buffer = read_input()
         read (buffer, *) number_of_sharings
         call print_input("Enter the list of sharings:") 
         buffer = read_input()
         read (buffer, *) (sharing(i),i=1,number_of_sharings)
         is_gl = .false.
      end if
      !
   end subroutine io_sharings
   !
   !
   subroutine io_time_unit()
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in and determines the time unit for input/output.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=256) :: buffer, record
      !
      ! Determine the units for the printout and further optional input
    1 call print_input("Which units are to be used to enter and to print the time ?")
      call print_input("    sec     : seconds;")
      call print_input("    fs      : femtoseconds;")
      call print_input("    a.u.    : atomic units;")
      !!x print *, "Which units are to be used to enter and to print the time ?"
      !!x print *, "    sec     : seconds;"
      !!x print *, "    fs      : femtoseconds;"
      !!x print *, "    a.u.    : atomic units;"
      buffer = read_input()
      read (buffer, "(a)") record
      if (len_trim(record) > 0) then
         record = adjustl(record)
         select case(record(1:7))
         case("sec    ", "fs     ", "a.u.   ")
            time_unit = record(1:7)      
         case default
            !!x print *, "Unable to decode the units '"//record(1:7)// "' to be used for the time; reenter ..."
            call print_input("Unable to decode the units '"//record(1:7)// "' to be used for the time; reenter ...")
            goto 1
         end select
      else 
         time_unit = "sec    "
      end if
      !
      !!x call save_input(record,.false.)
      !!x !
      ! Determine conversion factor for calculating times in this unit
      select case(time_unit)
      case("sec    ")
         time_factor  = bohr_radius_in_cm * one_over_alpha / c_vacuum_in_cm_per_s
      case("fs     ")
         time_factor  = bohr_radius_in_cm * one_over_alpha / c_vacuum_in_cm_per_s * 1.0e15
      case("a.u.   ")
         time_factor  = one
      end select
      !
   end subroutine io_time_unit
   !
   !
   subroutine io_transition_multipoles(number_of_multipoles,multipole)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in the transition multipoles from a given record.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(inout)                        :: number_of_multipoles
      character(len=2), dimension(:), intent(inout) :: multipole
      !
      character(len=256) :: buffer, record
      !
    1 call print_input("Enter the transition multipoles, e.g.  E1 M2 ... :")
      !!x print *, "Enter the transition multipoles, e.g.  E1 M2 ... :"
      buffer = read_input()
      read (buffer, "(a)") record
      !!x call save_input(record,.false.)
      !
    2 if (len_trim(record) > 0) then
         record = adjustl(record)
         select case(record(1:2))
         case("E1", "E2", "E3", "E4", "E5", "M1", "M2", "M3", "M4", "M5")
            number_of_multipoles = number_of_multipoles + 1
            if (number_of_multipoles > 20) then
               stop "input_transition_multipoles(): program stop A."
            end if
            multipole(number_of_multipoles) = record(1:2)
            record(1:2) = "  "
            goto 2
         case default
            number_of_multipoles = 0
            call print_input("Unable to decode the transition multipole '"// record(1:2)//"'; reenter ...")
            !!x print *, "Unable to decode the transition multipole '"// record(1:2)//"'; reenter ..." 
            goto 1
         end select
      else if (number_of_multipoles == 0) then
         !!x print *, "At least one transition multipole must be specified; reenter ..."
         call print_input("At least one transition multipole must be specified; reenter ...")
         goto 1
      end if
      !
   end subroutine io_transition_multipoles
   !
   !
   subroutine io_transition_pairs(number_of_transitions)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in and selects pairs of transitions.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(inout) :: number_of_transitions
      !
      integer                :: level_i, level_f, score_position
      logical                :: yes
      character(len=256)     :: buffer, record
      character(len=20 )     :: string
      !
      !!x print *, "Select individual transitions ?"
      !!x yes = get_yes_stream()
      call print_input("Select individual transitions ?")
      yes = get_yesno()
      if (yes) then
    4    call print_input("Enter one pair |i> - |f> of level numbers, e.g. 2 - 3; 2 - 0; 0 - 17 ... ")
         call print_input(" (0 is here equivalent to all); <cr> if done.")
         !!x print *, "Enter one pair |i> - |f> of level numbers, e.g. 2 - 3; 2 - 0; 0 - 17 ... "
         !!x print *, " (0 is here equivalent to all); <cr> if done."
         buffer = read_input()
         read (buffer, "(a)") record
         !!x call save_input(record,.false.)
	 !
         if (len_trim(record) > 0) then
            score_position = scan(record,"-")
            if (score_position == 0) then
               !!x print *, "Unable to decode the transition; reenter ..."
               call print_input("Unable to decode the transition; reenter ...")
               goto 4
            else 
               string  = adjustl(record(1:score_position-1))
               level_i = get_integer_from_string(string)
               string  = adjustl(record(score_position+1:256))
               level_f = get_integer_from_string(string)
               !!x print *, "level_i, level_f = ",level_i, level_f
               if (level_i < 0   .or.   level_f < 0) then
                  !!x print *, "All level numbers must greater or equivalent 0; reenter ..."
                  call print_input("All level numbers must greater or equivalent 0; reenter ...")
                  goto 4
               elseif (level_i == 0   .and.   level_f == 0) then
                  !!x print *, "The pair of level numbers 0 - 0 is not allowed; reenter ...";   
                  call print_input("The pair of level numbers 0 - 0 is not allowed; reenter ...")
                  goto 4
               endif
               !
               number_of_transitions = number_of_transitions + 1
               if (number_of_transitions > 200) then
                  stop "input_transition_pairs(): program stop A."
               end if
               select_level_i(number_of_transitions) = level_i
               select_level_f(number_of_transitions) = level_f
               goto 4
            end if
         else if (number_of_transitions == 0) then
            !!x print *, "No pair of level numbers has been selected; the  program continues the default path."
            call print_input("No pair of level numbers has been selected; the  program continues the default path.")
         end if
      end if
      !
   end subroutine io_transition_pairs
   !
   !
   subroutine io_transition_triples(number_of_transitions)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads in and selects triples of transitions.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(inout) :: number_of_transitions
      !
      integer                :: level_i, level_m, level_f, score_position
      logical                :: yes
      character(len=256)     :: buffer, record
      character(len=20 )     :: string
      !
      !!x print *, "Select individual transitions |i> --> |m> --> |f> ?"
      !!x yes = get_yes_stream()
      call print_input("Select individual transitions |i> --> |m> --> |f> ?")
      yes = get_yesno()
      if (yes) then
    6    call print_input("Enter one triple |i> - |m> - |f> of level numbers, e.g. 2 - 3 - 1; 2 - 0 - 5; 0 - 17 - 3 ... ")
         call print_input(" (0 is here equivalent to all); <cr> if done.")
         !!x print *, "Enter one triple |i> - |m> - |f> of level numbers, e.g. 2 - 3 - 1; 2 - 0 - 5; 0 - 17 - 3 ... "
         !!x print *, " (0 is here equivalent to all); <cr> if done."
         buffer = read_input()
         read (buffer, "(a)") record
         if (len_trim(record) > 0) then
            score_position = scan(record,"-")
            if (score_position == 0) then
               call print_input("Unable to decode the transition; reenter ...")
               !!x print *, "Unable to decode the transition; reenter ..."
               goto 6
            else 
               string  = adjustl(record(1:score_position-1))
               level_i = get_integer_from_string(string)
               record  = adjustl(record(score_position+1:256))
               score_position = scan(record,"-")
               if (score_position == 0) then
                  call print_input("Unable to decode the transition; reenter ...")
                  !!x print *, "Unable to decode the transition; reenter ..."
                  goto 6
               else 
                  string  = adjustl(record(1:score_position-1))
                  level_m = get_integer_from_string(string)
                  string  = adjustl(record(score_position+1:256))
                  level_f = get_integer_from_string(string)
               end if
               !
               !!x print *, "level_i, level_m, level_f = ",level_i,level_m,level_f
               if (level_i < 0   .or.   level_m < 0   .or.   level_f < 0) then
                  !!x print *, "All level numbers must greater or equivalent 0; reenter ..."  
                  call print_input("All level numbers must greater or equivalent 0; reenter ...")
                  goto 6
               elseif (level_i == 0  .and. level_i == 0 .and. level_f == 0) then
                  !!x print *, "The pair of level numbers 0 - 0 - 0 is not allowed; reenter ..."  
                  call print_input("The pair of level numbers 0 - 0 - 0 is not allowed; reenter ...")
                  goto 6
               endif
               !
               number_of_transitions = number_of_transitions + 1
               if (number_of_transitions > 200) then
                  stop "dierec_collect_input(): program stop B."
               end if
               select_level_i(number_of_transitions) = level_i
               select_level_m(number_of_transitions) = level_m
               select_level_f(number_of_transitions) = level_f
               goto 6
            end if
         else if (number_of_transitions == 0) then
            !!x print *, "No pair of level numbers has been selected; the program continues the default path."
            call print_input("No pair of level numbers has been selected; the program continues the default path.")
         end if
      end if
      !
   end subroutine io_transition_triples
   !
   !
   subroutine print_input(string)                                           
   !--------------------------------------------------------------------------------------------------------------------
   ! Prints some records to the standard stream (*) but keeps a `copy' for the input buffer.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in) :: string
      !
      no_of_input               = no_of_input + 1
      buffer%input(no_of_input) = trim(string)
      print *, trim(string)
      !
   end subroutine print_input
   !
   !
   subroutine print_runtime_component(stream, progstring)                                           
   !--------------------------------------------------------------------------------------------------------------------
   ! Prints the components name and the date, time and compiler (options) for the present run.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)          :: stream
      character(len=*), intent(in) :: progstring
      !
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of the .sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(stream,*) progstring // " run at " // ctime(1:2) // ":" //ctime(3:4) // ":" // ctime(5:6) // &
                      " on " // month // " "// cdate(7:8) // " " // cdate(1:4) // "."
      write(stream,*) " "
      write(stream,*) "   Compiler version:  ",compiler_version()
      write(stream,*) "   Compiler options:  ",compiler_options()
      write(stream,*) " "
      !
   end subroutine print_runtime_component
   !
   !
   function read_input(string)                                                                              result(wa)
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads a record from standard stream (*) and stores it into the input buffer; it also returns the string with variable
   ! length to assign it to some local `buffer' for internal file reading.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in), optional :: string
      !
      character(len=250)  :: record
      character(len=250)  :: wa
      !
      read (*, "(a)", err=2) record
      !
      no_of_input               = no_of_input + 1
      buffer%input(no_of_input) = trim(record)
      wa                        = trim(record)
      goto 3
      !
    2 if (present(string)) then
         print *, string
         stop "read_input(): program stop A."
      else
         stop "read_input(): program stop A."
      end if
      !
    3 continue
      !
   end function read_input
   !
   !
   subroutine set_warning(string)                                           
   !--------------------------------------------------------------------------------------------------------------------
   ! Writes a warning to the warning buffer.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in) :: string
      !
      no_of_warnings                  = no_of_warnings + 1
      buffer%warnings(no_of_warnings) = trim(string)
      !
   end subroutine set_warning
   !
   !
   subroutine write_input(stream)                                           
   !--------------------------------------------------------------------------------------------------------------------
   ! Writes all input records to stream with some proper begin/end line.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      integer             :: i
      !
      write(stream, "(1x,a)") "The following input data were used for the computations:"
      write(stream, "(1x,a)") "--------------------------------------------------------"
      do i = 1,no_of_input
         write(stream, "(4x,a)") trim(buffer%input(i))
      end do
      write(stream, "(1x,a)") "--- end of input ---"
      !
   end subroutine write_input
   !
   !
   subroutine write_warnings(stream)                                           
   !--------------------------------------------------------------------------------------------------------------------
   ! Writes all warning records to stream with some proper begin/end line.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      integer             :: i
      !
      write(stream,*)         " "
      write(stream, "(1x,a)") "The following warnings were issued during the computations:"
      write(stream, "(1x,a)") "-----------------------------------------------------------"
      do i = 1,no_of_warnings
         write(stream, "(4x,a)") buffer%warnings(i)
      end do
      write(stream, "(1x,a)") "--- end of warnings ---"
      !
   end subroutine write_warnings
   !
end module rabs_io_dialog

