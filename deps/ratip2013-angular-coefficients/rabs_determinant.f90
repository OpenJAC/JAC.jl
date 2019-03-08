module rabs_determinant
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module provides abstract data types and procedures to facilitate the work with determinants in atomic structure 
! calculations. The two main data types define a single 'determinant' and a 'determinant basis' for which the occupation 
! of the individual determinants are packed into standard integer arrays. 
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_dirac_orbital
   !?? use rabs_utilities
   implicit none
   !
   public  :: cofactor_1_of_overlap_matrix
                 ! Calculates under a logical mask the n x n minorants (cofactor D(k,l)) for a quadratic matrix of overlap-
                 ! integrals as given by two determinants.
   public  :: cofactor_2_of_overlap_matrix
                 ! Calculates under a logical mask the n x n x n x n minorants (cofactor D(k1,k2;l1,l2)) for a quadratic 
                 ! matrix of overlap-integrals as given by two determinants.
   public  :: distribute_m_in_ncells
                 ! Returns all possible distributions of M particle in N different cells.
   public  :: pack_occupation_in_integer
                 ! Packs the occupation numbers 0 and 1 of a determinant as 'bits' into an array of standard integers.
   public  :: pack_occupation_in_integer_nr
                 ! Packs the occupation numbers 0 and 1 of a nonrelativistic determinant as 'bits' into an array of 
                 ! standard integers.
   public  :: readwrite_det_basis
                 ! Reads or writes a type(det_basis) data structure from or to a file.
   public  :: unpack_occupation_from_integer
                 ! Unpacks the occupation numbers 0 and 1 of a determinant as 'bits' from an array of standard integers.
   public  :: unpack_occupation_from_int_nr
                 ! Unpacks the occupation numbers 0 and 1 of a nonrelativistic determinant as 'bits' from an array of 
                 ! standard integers.
   !
   type, public :: determinant
      integer(kind=i1b) :: totalM
      character(len=1)  :: parity
      integer, dimension(:), pointer :: occupation
   end type determinant
   !
   type, public :: det_basis
      integer :: nod, nodmax, norbital, number_of_electrons, noint
      type(nkappam), dimension(:), pointer     :: orbital
      type(determinant), dimension(:), pointer :: determinant
   end type det_basis
   !
   type, public :: nrdeterminant
      integer(kind=i1b) :: totalM, totalMs
      character(len=1)  :: parity
      integer, dimension(:), pointer :: occupation
   end type nrdeterminant
   !
   type, public :: nrdet_basis
      integer :: nod, nodmax, norbital, number_of_electrons, noint
      type(nlmms), dimension(:), pointer         :: orbital
      type(nrdeterminant), dimension(:), pointer :: determinant
   end type nrdet_basis
   !
contains
   !
   subroutine cofactor_1_of_overlap_matrix(norb,overlap,cofactors_1, cofactors_1_mask)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates under control of the logical cofactors_1_mask the norb x norb minorants (cofactor D(k,l)) for a quadratic 
   ! matrix overlap of overlap-integrals according to the two determinants orbital_f and orbital_i.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                              :: norb
      logical, dimension(norb,norb), intent(in)        :: cofactors_1_mask
      real(kind=dp), dimension(norb,norb), intent(in)  :: overlap
      real(kind=dp), dimension(norb,norb), intent(out) :: cofactors_1
      !
      integer       :: ki, kf
      real(kind=dp) :: minor
      real(kind=dp), dimension(norb-1,norb-1) :: matrix
      !
      ! Calculate the one-particle co-factors D(k/l)
      do  kf = 1,norb
         do  ki = 1,norb
            if (cofactors_1_mask(kf,ki)) then
               matrix(1:kf-1,1:ki-1)       = overlap(1:kf-1,1:ki-1)
               matrix(1:kf-1,ki:norb-1)    = overlap(1:kf-1,ki+1:norb)
               matrix(kf:norb-1,1:ki-1)    = overlap(kf+1:norb,1:ki-1)
               matrix(kf:norb-1,ki:norb-1) = overlap(kf+1:norb,ki+1:norb)
               !! lf = 0
               !! do  kff = 1,norb
               !!    if (kff /= kf) then
               !!       lf = lf + 1
               !!       li = 0
               !!       do  kii = 1,norb
               !!          if (kii /= ki) then
               !!             li = li + 1
               !!             matrix(lf,li) = overlap(kff,kii)
               !!             !!print *, "lf,li,matrix(lf,li) = ",lf,li,matrix(lf,li)
               !!          end if
               !!       end do
               !!    end if
               !! end do
               !
               ! Calculate the value of determinant
               call calculate_determinant(matrix,minor,norb-1)
	       !
               if (abs(minor) < eps20) then
                  cofactors_1(kf,ki) = zero
               else
                  cofactors_1(kf,ki) = minor
               end if
            else
               cofactors_1(kf,ki) = zero
            end if
         end do
      end do
      !
   end subroutine cofactor_1_of_overlap_matrix
   !
   !
   subroutine cofactor_2_of_overlap_matrix(norb,overlap,cofactors_2, cofactors_2_mask)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates under control of the logical cofactors_2_mask the norb x norb x norb x norb minorants (cofactor D(k,l)) 
   ! for a quadratic matrix of overlap-integrals according to the two determinants orbital_f and orbital_i.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                                        :: norb
      logical, dimension(norb,norb,norb,norb), intent(in)        :: cofactors_2_mask
      real(kind=dp), dimension(norb,norb), intent(in)            :: overlap
      real(kind=dp), dimension(norb,norb,norb,norb), intent(out) :: cofactors_2
      !
      integer       :: kr1, kr2, ks1, ks2, krr, krx, kss, ksx
      real(kind=dp) :: minor
      real(kind=dp), dimension(norb-2,norb-2) :: matrix
      !
      cofactors_2(:,:,:,:) = zero
      !
      ! Calculate the one-particle co-factors D(k1,k2/l1,l2)
      do  kr1 = 1,norb
      do  kr2 = 1,norb
         if (kr1 == kr2) cycle
	 do ks1 = 1,norb
	 do ks2 = 1,norb
            if (ks1 == ks2) cycle
            if (cofactors_2_mask(kr1,kr2,ks1,ks2)) then
	       krr = 0
	       do  krx = 1,norb
	          if (krx == kr1  .or.  krx == kr2) cycle
		  krr = krr + 1;   kss = 0
	          do  ksx = 1,norb
	             if (ksx == ks1  .or.  ksx == ks2) cycle
		     kss = kss + 1
		     matrix(krr,kss) = overlap(krx,ksx)
		  end do
	       end do
               !
               ! Calculate the value of determinant
               call calculate_determinant(matrix,minor,norb-2)
               if (abs(minor) < eps20) then
                  cofactors_2(kr1,kr2,ks1,ks2) = zero
               else
                  cofactors_2(kr1,kr2,ks1,ks2) = minor
               end if
            else
               cofactors_2(kr1,kr2,ks1,ks2) = zero
            end if
         end do
         end do
      end do
      end do
      !
   end subroutine cofactor_2_of_overlap_matrix
   !
   !
   subroutine distribute_m_in_ncells(m,n,d,dmax,distribution)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns all possible distributions of M particle in N different cells. In this version, the value of M must be m <= 7.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)  :: m, n, dmax
      integer, intent(out) :: d
      integer, dimension(1:n,1:dmax), intent(out) :: distribution
      !
      integer :: i1, i2, i3, i4, i5, i6, i7
      !
      do  d = 1,dmax;   do  i1 = 1,n
         distribution(i1,d) = 0
      end do;   end do
      !
      d = 0
      !
      select case(m)
      case(1)
         do  i1 = 1,n
            d = d + 1
            distribution(i1,d) = 1
         end do
      case(2)
         do  i1 = 1,n-1;   do  i2 = i1+1,n
            d = d + 1
            distribution(i1,d) = 1;   distribution(i2,d) = 1
         end do;   end do
      case(3)
         do  i1 = 1,n-2;   do  i2 = i1+1,n-1;   do  i3 = i2+1,n
            d = d + 1
            distribution(i1,d) = 1;   distribution(i2,d) = 1;   distribution(i3,d) = 1
         end do;   end do;   end do
      case(4)
         do  i1 = 1,n-3;   do  i2 = i1+1,n-2;   do  i3 = i2+1,n-1;  do  i4 = i3+1,n
            d = d + 1
            distribution(i1,d) = 1;   distribution(i2,d) = 1;  distribution(i3,d) = 1;   distribution(i4,d) = 1
         end do;   end do;   end do;   end do
      case(5)
         do  i1 = 1,n-4;     do  i2 = i1+1,n-3;   do  i3 = i2+1,n-2;   do  i4 = i3+1,n-1;  do  i5 = i4+1,n
            d = d + 1
            distribution(i1,d) = 1;   distribution(i2,d) = 1;  distribution(i3,d) = 1;   distribution(i4,d) = 1
            distribution(i5,d) = 1
         end do;   end do;   end do;   end do;   end do
      case(6)
         do  i1 = 1,n-5;     do  i2 = i1+1,n-4;   do  i3 = i2+1,n-3;   do  i4 = i3+1,n-2;  do  i5 = i4+1,n-1
         do  i6 = i5+1,n
            d = d + 1
            distribution(i1,d) = 1;   distribution(i2,d) = 1;  distribution(i3,d) = 1;   distribution(i4,d) = 1
            distribution(i5,d) = 1;   distribution(i6,d) = 1
         end do;   end do;   end do;   end do;   end do;   end do
      case(7)
         do  i1 = 1,n-6;     do  i2 = i1+1,n-5;   do  i3 = i2+1,n-4;   do  i4 = i3+1,n-3;  do  i5 = i4+1,n-2
         do  i6 = i5+1,n-1;  do  i7 = i6+1,n
            d = d + 1
            distribution(i1,d) = 1;   distribution(i2,d) = 1;   distribution(i3,d) = 1;   distribution(i4,d) = 1
            distribution(i5,d) = 1;   distribution(i6,d) = 1;   distribution(i7,d) = 1
         end do;   end do;   end do;   end do;   end do;   end do;   end do
      case default
         print *, "distribute_m_in_ncells(): program stop A."
         stop
      end select
      !
   end subroutine distribute_m_in_ncells
   !
   !
   subroutine pack_occupation_in_integer(determ,noint,occupation,norb)
   !--------------------------------------------------------------------------------------------------------------------
   ! Packs the occupation numbers 0 and 1 of a determinant as 'bits' into an array of standard integers. It assumes that 
   ! the array determinant%occupation(1:) has properly been allocated before.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: noint,norb
      integer, dimension(1:norb), intent(in) :: occupation
      type(determinant), intent(inout)       :: determ
      !
      integer :: i, int, nobit, pos
      !
      ! Clear up the integer storage of the determinant
      nobit = bit_size(determ%occupation(1))
      do  int = 1,noint;   do  pos = 0,nobit-1
         determ%occupation(int) = ibclr(determ%occupation(int),pos)
      end do;   end do
      !
      int   = 0
      do  i = 1,norb
         if (mod(i-1,nobit) == 0) then
            int = int + 1
         end if
         pos = mod(i,nobit) - 1
         if (pos == -1) pos = nobit -1
         !
         if (occupation(i) == 0) then
            determ%occupation(int) = ibclr(determ%occupation(int),pos)
         else if (occupation(i) == 1) then
            determ%occupation(int) = ibset(determ%occupation(int),pos)
         else if (rabs_use_stop) then
            print *, "pack_occupation_in_integer(): program stop A."
            stop
         end if
      end do
      !
   end subroutine pack_occupation_in_integer
   !
   !
   subroutine pack_occupation_in_integer_nr(nrdeterm,noint,occupation,norb)
   !--------------------------------------------------------------------------------------------------------------------
   ! Packs the occupation numbers 0 and 1 of a non-relativistic determinant as 'bits' into an array of standard integers. 
   ! It assumes that the array determinant%occupation(1:) has properly been allocated before.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: noint, norb
      integer, dimension(1:norb), intent(in) :: occupation
      type(nrdeterminant), intent(inout)     :: nrdeterm
      !
      integer :: i, int, nobit, pos
      !
      ! Clear up the integer storage of the determinant
      nobit = bit_size(nrdeterm%occupation(1))
      do  int = 1,noint;   do  pos = 0,nobit-1
         nrdeterm%occupation(int) = ibclr(nrdeterm%occupation(int),pos)
      end do;   end do
      !
      int   = 0
      do  i = 1,norb
         if (mod(i-1,nobit) == 0) then
            int = int + 1
         end if
         pos = mod(i,nobit) - 1
         if (pos == -1) pos = nobit -1
         !
         if (occupation(i) == 0) then
            nrdeterm%occupation(int) = ibclr(nrdeterm%occupation(int),pos)
         else if (occupation(i) == 1) then
            nrdeterm%occupation(int) = ibset(nrdeterm%occupation(int),pos)
         else if (rabs_use_stop) then
            print *, "pack_occupation_in_integer_nr(): program stop A."
            stop
         end if
      end do
      !
   end subroutine pack_occupation_in_integer_nr
   !
   !
   subroutine readwrite_det_basis(stream,read_from,det_set)             
   !--------------------------------------------------------------------------------------------------------------------
   ! Reads or writes a type(det_basis) data structure from or to a file on stream.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)            :: stream
      logical, intent(in)            :: read_from
      type(det_basis), intent(inout) :: det_set
      !
      integer :: i, j
      !
      if (read_from) then
         read(stream,*) det_set%nod, det_set%nodmax, det_set%norbital, det_set%number_of_electrons, det_set%noint
	 allocate( det_set%orbital(1:det_set%norbital), det_set%determinant(1:det_set%nodmax) )
      else
         write(stream,*) det_set%nod, det_set%nodmax, det_set%norbital, det_set%number_of_electrons, det_set%noint 
      end if
      !
      do  i = 1,det_set%norbital
         if (read_from) then
	    read(stream,*)  det_set%orbital(i)%n, det_set%orbital(i)%kappa, det_set%orbital(i)%mm
	 else
	    write(stream,*) det_set%orbital(i)%n, det_set%orbital(i)%kappa, det_set%orbital(i)%mm
	 end if
      end do
      !
      do  i = 1,det_set%nodmax
         if (read_from) then
	    read(stream,*) det_set%determinant(i)%totalM, det_set%determinant(i)%parity
            allocate( det_set%determinant(i)%occupation(1:det_set%noint) )
            read(stream,*) (det_set%determinant(i)%occupation(j), j=1,det_set%noint)
	 else
	    write(stream,*) det_set%determinant(i)%totalM, det_set%determinant(i)%parity
            write(stream,*) (det_set%determinant(i)%occupation(j), j=1,det_set%noint)
	 end if
      end do
      !
   end subroutine readwrite_det_basis
   !
   !
   subroutine unpack_occupation_from_integer(determ,occupation,norb)
   !--------------------------------------------------------------------------------------------------------------------
   ! Unpacks the occupation numbers 0 and 1 of the determinant determ which are given as 'bits' from an array of standard 
   ! integers.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                     :: norb
      integer, dimension(1:norb), intent(out) :: occupation
      type(determinant), intent(in)           :: determ
      !
      integer :: i, nobit, noint, pos
      !
      nobit = bit_size(determ%occupation(1))
      do  i = 1,norb
         if (mod(i,nobit) == 0) then
            noint = i/nobit
            pos   = nobit - 1
         else
            noint = i/nobit + 1
            pos   = mod(i,nobit) - 1
         end if
         !
         if (btest(determ%occupation(noint),pos)) then
            occupation(i) = 1
         else
            occupation(i) = 0
         end if
      end do
      !
   end subroutine unpack_occupation_from_integer
   !
   !
   subroutine unpack_occupation_from_int_nr(nrdeterm,occupation,norb)
   !--------------------------------------------------------------------------------------------------------------------
   ! Unpacks the occupation numbers 0 and 1 of the determinant determ which are given as 'bits' from an array of standard 
   ! integers.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)                     :: norb
      integer, dimension(1:norb), intent(out) :: occupation
      type(nrdeterminant), intent(in)         :: nrdeterm
      !
      integer :: i, nobit, noint, pos
      !
      nobit = bit_size(nrdeterm%occupation(1))
      do  i = 1,norb
         if (mod(i,nobit) == 0) then
            noint = i/nobit
            pos   = nobit - 1
         else
            noint = i/nobit + 1
            pos   = mod(i,nobit) - 1
         end if
         !
         if (btest(nrdeterm%occupation(noint),pos)) then
            occupation(i) = 1
         else
            occupation(i) = 0
         end if
      end do
      !
   end subroutine unpack_occupation_from_int_nr
   !
end module rabs_determinant
