module rabs_functions_string
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module contains a number of string and related functions as they appear frequently in the computation of atomic 
! properties and structures.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   implicit none
   private
   !
   public  :: center_unit
                 ! Returns a string of given length where a given unit is centered in parenthesis.
   public  :: get_dinteger_from_string
                 ! Returns the double integer value as described by a given string.
   public  :: get_integer_from_string
                 ! Returns the integer value as described by a given string.
   public  :: get_month        
                 ! Returns a string (len = 3) for a given month.
   public  :: interprete_levels     
                 ! Interpretes a string of level numbers for further processing.
   public  :: interprete_orbitals    
                 ! Interpretes a string of orbital notations for further processing.
   public  :: parity_add
                 ! Returns "+" or "-" due to the summation of two parities.
   !
   !
   character(len=3), dimension(-70:70), private, parameter :: int_string =   &
    (/ "-70", "-69", "-68", "-67", "-66", "-65", "-64", "-63", "-62", "-61", &
       "-60", "-59", "-58", "-57", "-56", "-55", "-54", "-53", "-52", "-51", &
       "-50", "-49", "-48", "-47", "-46", "-45", "-44", "-43", "-42", "-41", &
       "-40", "-39", "-38", "-37", "-36", "-35", "-34", "-33", "-32", "-31", &
       "-30", "-29", "-28", "-27", "-26", "-25", "-24", "-23", "-22", "-21", &
       "-20", "-19", "-18", "-17", "-16", "-15", "-14", "-13", "-12", "-11", &
       "-10", " -9", " -8", " -7", " -6", " -5", " -4", " -3", " -2", " -1", &
       "  0",                                                                &
       "  1", "  2", "  3", "  4", "  5", "  6", "  7", "  8", "  9", " 10", &
       " 11", " 12", " 13", " 14", " 15", " 16", " 17", " 18", " 19", " 20", &
       " 21", " 22", " 23", " 24", " 25", " 26", " 27", " 28", " 29", " 30", &
       " 31", " 32", " 33", " 34", " 35", " 36", " 37", " 38", " 39", " 40", &
       " 41", " 42", " 43", " 44", " 45", " 46", " 47", " 48", " 49", " 50", &
       " 51", " 52", " 53", " 54", " 55", " 56", " 57", " 58", " 59", " 60", &
       " 61", " 62", " 63", " 64", " 65", " 66", " 67", " 68", " 69", " 70"  /)
   !
   character(len=5), dimension(0:160), private, parameter :: dint_string =   &
    (/ "    0", &
       "  1/2","    1","  3/2","    2","  5/2","    3","  7/2","    4",  &
       "  9/2","    5"," 11/2","    6"," 13/2","    7"," 15/2","    8",  & 
       " 17/2","    9"," 19/2","   10"," 21/2","   11"," 23/2","   12",  &
       " 25/2","   13"," 27/2","   14"," 29/2","   15"," 31/2","   16",  &
       " 33/2","   17"," 35/2","   18"," 37/2","   19"," 39/2","   20",  &
       " 41/2","   21"," 43/2","   22"," 45/2","   23"," 47/2","   24",  &
       " 49/2","   25"," 51/2","   26"," 53/2","   27"," 55/2","   28",  &
       " 57/2","   29"," 59/2","   30"," 61/2","   31"," 63/2","   32",  &
       " 65/2","   33"," 67/2","   34"," 69/2","   35"," 71/2","   36",  &
       " 73/2","   37"," 75/2","   38"," 77/2","   39"," 79/2","   40",  &
       " 81/2","   41"," 83/2","   42"," 85/2","   43"," 87/2","   44",  &
       " 89/2","   45"," 91/2","   46"," 93/2","   47"," 95/2","   48",  &
       " 97/2","   49"," 99/2","   50","101/2","   51","103/2","   52",  &
       "105/2","   53","107/2","   54","109/2","   55","111/2","   56",  &
       "113/2","   57","115/2","   58","117/2","   59","119/2","   60",  &
       "121/2","   61","123/2","   62","125/2","   63","127/2","   64",  &
       "129/2","   65","131/2","   66","133/2","   67","135/2","   68",  &
       "137/2","   69","139/2","   70","141/2","   71","143/2","   72",  &
       "145/2","   73","147/2","   74","149/2","   75","151/2","   76",  &
       "153/2","   77","155/2","   78","157/2","   79","159/2","   80"   /)
   !
   character(len=3), dimension(1:12), private, parameter :: month_string =  &
    (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)
   !
   !
contains
   !
   !
   function center_unit(string,length)                                  result(ws)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a string of given length where a given unit is centered in parenthesis.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in) :: string 
      integer, intent(in)          :: length
      character(len=30)            :: ws 
      !
      integer           :: lx, ly
      character(len=10) :: wa = "          "
      !
      lx = len_trim(string)
      ly = (length - lx - 2)/2
      !
      if (ly < 0) then
         stop "center_unit(): program stop A."
      end if
      !
      ws = wa(1:ly) // "(" // trim(string) // ")" // wa
      !
   end function center_unit
   !
   !
   function get_dinteger_from_string(string,fail)              result(L)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the double of the integer or half-integer value which is described by string.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)   :: string
      logical, optional, intent(out) :: fail 
      integer                        :: i, L
      character(len=5)               :: stringr
      !
      l = len(trim(string))
      stringr = "     "
      stringr(5-l+1:5) = trim(string)
      do  i = 0,160
         if (stringr == dint_string(i)) then
            L = i
            return
         end if
      end do
      !
      if (present(fail)) then
         fail = .true.
      else if (rabs_use_stop) then
         stop "get_dinteger_from_string(): program stop A."
      end if
      !
   end function get_dinteger_from_string 
   !
   !
   function get_integer_from_string(string,fail)               result(L)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the integer value which is described by string. The optional argument fail is set to .true. if string 
   ! cannot be interpreted correctly; if fail is not present, the procedure terminates.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)   :: string
      logical, optional, intent(out) :: fail 
      integer                        :: i, L
      logical                        :: minus
      !
      character(len=10) :: stringr
      !
      stringr = adjustl(string)
      stringr = adjustr(stringr)
      do  i = -70,70
         if (stringr(8:10) == int_string(i)) then
            L = i
            if (present(fail))  fail = .false.
            return
         end if
      end do
      !
      L     = 0
      minus = .false.
      do  i = 1,10
         if (stringr(i:i) == " "  .and.  minus) then
            stop "get_integer_from_string(): program stop A."
         else if (stringr(i:i) == " ") then
            cycle
         else if (stringr(i:i) == "-"  .and.  minus) then
            stop "get_integer_from_string(): program stop B."
         else if (stringr(i:i) == "-") then
            minus = .true.
         else if (stringr(i:i) == "0") then 
            L = 10*L + 0
         else if (stringr(i:i) == "1") then 
            L = 10*L + 1
         else if (stringr(i:i) == "2") then 
            L = 10*L + 2
         else if (stringr(i:i) == "3") then 
            L = 10*L + 3
         else if (stringr(i:i) == "4") then 
            L = 10*L + 4
         else if (stringr(i:i) == "5") then 
            L = 10*L + 5
         else if (stringr(i:i) == "6") then 
            L = 10*L + 6
         else if (stringr(i:i) == "7") then 
            L = 10*L + 7
         else if (stringr(i:i) == "8") then 
            L = 10*L + 8
         else if (stringr(i:i) == "9") then 
            L = 10*L + 9
         end if
      end do
      !
      if  (minus)  then
         L = -L
      end if
      !
      return
      !
      if (present(fail)) then
         fail = .true.
      else if (rabs_use_stop) then
         stop "get_integer_from_string(): program stop C."
      end if
      !
   end function get_integer_from_string 
   !
   !
   function get_month(string)                               result(mstr)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the double of the integer or half-integer value which is described by string.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=2), intent(in) :: string
      character(len=3)             :: mstr
      !
      select case(string)
      case("01");   mstr = month_string(1)
      case("02");   mstr = month_string(2)
      case("03");   mstr = month_string(3)
      case("04");   mstr = month_string(4)
      case("05");   mstr = month_string(5)
      case("06");   mstr = month_string(6)
      case("07");   mstr = month_string(7)
      case("08");   mstr = month_string(8)
      case("09");   mstr = month_string(9)
      case("10");   mstr = month_string(10)
      case("11");   mstr = month_string(11)
      case("12");   mstr = month_string(12)
      case default
         print *, "string = "//"'"//string//"'" 
         stop     "get_month: program stop A."
      end select
      !
   end function get_month
   !
   !
   subroutine interprete_levels(record,level,no_levels, max_level,fail)
   !--------------------------------------------------------------------------------------------------------------------
   ! Attempts to interprete the serial level numbers which are to be calculated in the present run from a given string 
   ! of numbers and 'intervals'. Levels can be given in the format:
   !                1 3 4  7 - 20  48  69 - 85
   ! The level numbers can be given in any order and 'overlapping' intervals are also allowed. The procedure returns 
   ! with fail = .true. if the level numbers cannot be interpreted properly (fail = .false. otherwise).
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)        :: record
      logical, intent(out)                :: fail
      integer, intent(out)                :: no_levels, max_level
      integer, dimension(:), pointer      :: level
      !
      logical, dimension(200)             :: low_to
      character(len=500)                  :: string
      integer                             :: a, i, lower, n
      integer, dimension(200)             :: low
      integer(kind=i1b), dimension(10000) :: eigpair
      !
      fail       = .true.
      !
      string = adjustl(record)
      !
      n = 0;   lower = 0
      low_to(:) = .false.
      !
    1 string = adjustl(string)
      if (string(1:1) == " ") then
         if (n == 0  .or.  low_to(n)) then
            return
         else
            goto 10
         end if
      else if (string(1:1) == "-") then
         if (n == 0)  return
         low_to(n)   = .true.
         string(1:1) = " "
         goto 1
      else if (string(1:1) == "0"   .or.   string(1:1) == "1" .or. string(1:1) == "2"   .or.   string(1:1) == "3" .or.  &
               string(1:1) == "4"   .or.   string(1:1) == "5" .or. string(1:1) == "6"   .or.   string(1:1) == "7" .or.  &
               string(1:1) == "8"   .or.   string(1:1) == "9") then
         a = get_integer_from_string(string(1:1))
         lower = 10*lower + a
         if (string(2:2) == " "   .or.   string(2:2) == "-") then
            n = n + 1
            low(n)      = lower
            lower       = 0
         end if
         string(1:1) = " "
         goto 1
      end  if
      !
      ! Determine no_levels and max_level
   10 eigpair(:) = 0
      do  i = 1,n
         if (low_to(i)) then
            if (low(i) <= low(i+1)) then;   eigpair(low(i):low(i+1)) = 1
            else;                           eigpair(low(i+1):low(i)) = 1
            end if
            !   
            if (low_to(i+1)) then;   return
            else;                    cycle
            end if
         else
            eigpair(low(i)) = 1
         end if
      end do
      !
      no_levels = 0;  max_level = 0
      do  i = 1,10000
         if (eigpair(i) == 1) then
            no_levels = no_levels + 1
            max_level = max( i, max_level )
         end if
      end do
      !
      allocate( level(1:no_levels) )     
      !
      no_levels = 0
      do  i = 1,10000
         if (eigpair(i) == 1) then
            no_levels = no_levels + 1
            level(no_levels) = i
         end if
      end do
      !
      fail = .false.
      !
   end subroutine interprete_levels
   !
   !
   subroutine interprete_orbitals(record,orbital,no_orbitals,max_orbitals,fail)
   !--------------------------------------------------------------------------------------------------------------------
   ! Attempts to interprete the orbital quantum numbers from a string of orbital notations including `*`notation for 
   ! all kappa's of a given n. The orbitals can be given in any order and may `overlap' with the `*` notation.
   ! The procedure returns with fail = .true. if the orbital notation cannot be interpreted properly (fail = .false. 
   ! otherwise).
   !
   ! ... max_orbitals -- maximal number of orbitals which can be returned in the record orbital.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=*), intent(in)      :: record
      integer, intent(in)               :: max_orbitals
      logical, intent(out)              :: fail
      integer, intent(out)              :: no_orbitals
      type(nkappa), dimension(:)        :: orbital
      !   
      integer             :: a, b, nc, pqn
      character(len=2)    :: wa
      character(len=120)  :: string
      !   
      fail        = .false.  
      no_orbitals = 0
      string      = adjustl(record)
      !
    1 string    = adjustl(string)
      !
      if (string(1:1) == " ") then
         return
      else if (string(1:1) == "0"   .or.   string(1:1) == "1" .or.  string(1:1) == "2"   .or.   string(1:1) == "3" .or.  &
               string(1:1) == "4"   .or.   string(1:1) == "5" .or.  string(1:1) == "6"   .or.   string(1:1) == "7" .or.  &
               string(1:1) == "8"   .or.   string(1:1) == "9") then
         a  = get_integer_from_string(string(1:1))
         nc = 1
         !
         if   (string(2:2) == "0"   .or.   string(2:2) == "1" .or.  string(2:2) == "2"   .or.   string(2:2) == "3" .or.  &
               string(2:2) == "4"   .or.   string(2:2) == "5" .or.  string(2:2) == "6"   .or.   string(2:2) == "7" .or.  &
               string(2:2) == "8"   .or.   string(2:2) == "9") then
            b = get_integer_from_string(string(2:2))
            pqn = a*10 + b
            nc  = 2
         else
            pqn = a
         end if
      else
         fail = .true.
      end if
      !
      wa = string(nc+1:nc+2)
      print *, "b: nc, wa, pqn, a, b = ",nc, wa, pqn, a, b
      if (wa == "* ") then
         select case(pqn)
         case (1:)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn, -1)
         end select
         !
         select case(pqn)
         case (2:)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn,  1)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn, -2)
         end select
         !
         select case(pqn)
         case (3:)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn,  2)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn, -3)
         end select
         !
         select case(pqn)
         case (4:)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn,  3)
            no_orbitals = no_orbitals + 1
            orbital(no_orbitals) = nkappa(pqn, -4)
         end select
         !
      else if (wa == "s ") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn, -1)
      else if (wa == "p-") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn,  1)
      else if (wa == "p ") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn, -2)
      else if (wa == "d-") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn,  2)
      else if (wa == "d ") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn, -3)
      else if (wa == "f-") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn,  3)
      else if (wa == "f ") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn, -4)
      else if (wa == "g-") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn,  4)
      else if (wa == "g ") then
         no_orbitals = no_orbitals + 1
         orbital(no_orbitals) = nkappa(pqn, -5)
      else
         fail = .true.
      end if
      !
      if (no_orbitals > max_orbitals - 1) then
         stop "interprete_orbitals(): program stop A."
      end if
      !
      string = string(nc+3:)
      goto 1
      !
   end subroutine interprete_orbitals
   !
   !
   function parity_add(pa,pb)                                 result(pc)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns "+" or "-" due to the sum of the two parities pa and pb.
   !--------------------------------------------------------------------------------------------------------------------
      !
      character(len=1), intent(in) :: pa, pb
      character(len=1)             :: pc
      !
      if (      (pa == "+"  .and.  pb == "+")  .or.  (pa == "-"  .and.  pb == "-") ) then 
         pc = "+"
      else if ( (pa == "+"  .and.  pb == "-")  .or.  (pa == "-"  .and.  pb == "+") ) then 
         pc = "-"
      else
         stop "parity_add(): program stop A."
      end if
      !
   end function parity_add
   !
end module rabs_functions_string
