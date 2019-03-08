module rabs_naglib
!
!-----------------------------------------------------------------------
! This module provides an easy and proper treatment of the NAG-library 
! procedures for the case that this library is NOT used 
! [cf.  rabs_use_naglib = .false.  in rabs_constant.f].
! Below, those procedures which occur in the RATIP package are provided
! as dummy procedure which stop the program if they are to be invoked.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   implicit none
   !
   private
   !
   public  :: nag_s14aaf 
   public  :: nag_s17def 
   !     
   contains
   !                        
   !                        
   function nag_s14aaf(x,ifail)                               result(wa)
   !--------------------------------------------------------------------
   ! 
   !--------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: x
      integer, intent(out)      :: ifail
      real(kind=dp)             :: wa
      !
      ifail= 0;  wa = x
      !
      print *, "*** ERROR: rabs_use_naglib = .true.  in rabs_constant.f ***"
      print *, "*** ERROR: rabs_use_naglib = .true.  in rabs_constant.f ***"
      print *, "*** ERROR: rabs_use_naglib = .true.  in rabs_constant.f ***"
      !
      stop "nag_s14aaf() in rabs_naglib: program stop A."
      !
   end function nag_s14aaf
   !
   !                        
   subroutine nag_s17def(nu,zz,k,string,cy,nz,ifail)   
   !--------------------------------------------------------------------
   ! 
   !--------------------------------------------------------------------
      !
      integer, intent(in)             :: k
      integer, intent(inout)          :: nz, ifail
      real(kind=dp), intent(inout)    :: nu
      complex(kind=dp), intent(inout) :: zz, cy
      character(len=*), intent(in)    :: string
      !
      ifail= 0; nz = k;  nu = zero;  zz = zero;  cy = zero;
      print *, "string = ",string
      !
      print *, "*** ERROR: rabs_use_naglib = .true.  in rabs_constant.f ***"
      print *, "*** ERROR: rabs_use_naglib = .true.  in rabs_constant.f ***"
      print *, "*** ERROR: rabs_use_naglib = .true.  in rabs_constant.f ***"
      !
      stop "nag_s17def() in rabs_naglib: program stop A."
      !
   end subroutine nag_s17def
   !
end module rabs_naglib
