!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This part of the program provides a set of functions which 'simulate' task from the GRASP92 package. These procedures 
! have the same names as before in GRASP92 (to enable a simple access of them) but get their data from structures as 
! defined in the module rabs_mcp. The procedures in this file become obsolet if the MCP and NJGRAF components are 
! 'replaced' by more efficient schemes in the future.
!
! Routines from this program part should not be invoked from any procedure of the varies RABS modules; they also do not 
! appear in the alphabetic list if procedures of the RABS package. In contrast, they 'use' data and procedures from RABS, 
! mainly rabs_constant, rabs_angular, and rabs_mcp.
!-----------------------------------------------------------------------------------------------------------------------
!
!
subroutine cfp(lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,coefp)
!-----------------------------------------------------------------------------------------------------------------------
! Selects a cfp coefficient from an appropriate table of fractional parentage  coefficients in jj-coupling. It has been 
! taken from the GRASP package and slightly adopted to the format of the RABS program.
!
! Input variables:
!
!   lock     : + or - (2*j + 1).
!   nel      : number of equivalent electrons in shell.
!   ijd/ijp  : total j of daughter/parent state.
!   ivd/ivp  : seniority of daughter/parent states.
!   iwd/iwp  : other quantum number (if needed).
!
! Output variable: 
!
!   coefp    : numerical result
!
! This control routine does not check the input variables for consistency, except the trivial case of j = 1/2. All  
! other checks are performed at a lower level. The package will return correct results for j = 3/2, 5/2, 7/2. Higher 
! values of j return a value 1.0 if NEL = 1 or 2; otherwise 0 with an error signal. 
!
! Calls: cfp_coefficient().
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   implicit none
   !
   integer, intent(in)        :: lock, nel, ijd, ivd, iwd, ijp, ivp, iwp
   real(kind=dp), intent(out) :: coefp
   !
   call cfp_coefficient(lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,coefp)
   !
end subroutine cfp
!
!
function ichop(isubsh,icsf)                                
!-----------------------------------------------------------------------------------------------------------------------
! ichop is -1 if subshell isubsh is empty in CSF  icsf,  +1 if the subshell is full, and 0 if it is open.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_mcp
   implicit none
   !
   integer, intent(in) :: isubsh, icsf
   integer             :: ichop, j
   !
   if (isubsh > two_csf%nwshells) then
      stop "ichop() in rabs_mcp_adaptation: program stop A."
   end if
   !
   j = angular_momentum_j(two_csf%subshell(isubsh)%kappa)
   if (icsf == csf_a_No) then
      if (two_csf%csf(1)%occupation(isubsh) == 0) then
         ichop = -1
      else if (two_csf%csf(1)%occupation(isubsh) == j + 1) then
         ichop =  1
      else
         ichop =  0
      end if
   else if (icsf == csf_b_No) then
      if (two_csf%csf(2)%occupation(isubsh) == 0) then
         ichop = -1
      else if (two_csf%csf(2)%occupation(isubsh) == j + 1) then
         ichop =  1
      else
         ichop =  0
      end if
   else
      stop "ichop() in rabs_mcp_adaptation: program stop B."
   end if
   !
end function ichop
!
!
function iq(isubsh,icsf)                                
!-----------------------------------------------------------------------------------------------------------------------
! iq  is the occupation of subshell isubsh in CSF  icsf.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_mcp
   implicit none
   !
   integer, intent(in) :: isubsh, icsf
   integer             :: iq
   !
   if (isubsh > two_csf%nwshells) then
      stop "iq() in rabs_mcp_adaptation: program stop A."
   else if (icsf == csf_a_No) then
      iq = two_csf%csf(1)%occupation(isubsh)
   else if (icsf == csf_b_No) then
      iq = two_csf%csf(2)%occupation(isubsh)
   else
      stop "iq() in rabs_mcp_adaptation: program stop B."
   end if
   !
end function iq
!
!
function ispar(icsf)                                
!-----------------------------------------------------------------------------------------------------------------------
! ispar  is the parity value of P (1 or -1) for CSF number icsf.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_mcp
   implicit none
   !
   integer, intent(in) :: icsf
   integer             :: ispar
   !
   if (icsf == csf_a_No) then
      if (two_csf%csf(1)%parity == "+") then
         ispar =  1
      else
         ispar = -1
      end if
   else if (icsf == csf_b_No) then
      if (two_csf%csf(2)%parity == "+") then
         ispar =  1
      else
         ispar = -1
      end if
   else
      stop "ispar() in rabs_mcp_adaptation: program stop A."
   end if
   !
end function ispar
!
!
function itjpo(icsf)                                
!-----------------------------------------------------------------------------------------------------------------------
! itjpo  is the value of 2J+1 for CSF number icsf.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_mcp
   implicit none
   !
   integer, intent(in) :: icsf
   integer             :: itjpo
   !
   if (icsf == csf_a_No) then
      itjpo = two_csf%csf(1)%totalJ+1
   else if (icsf == csf_b_No) then
      itjpo = two_csf%csf(2)%totalJ+1
   else
      stop "itjpo() in rabs_mcp_adaptation: program stop A."
   end if
   !
end function itjpo
!
!
function jcup(loc,icsf)                                
!-----------------------------------------------------------------------------------------------------------------------
! jcup is the 2J+1 value of the loc-th nontrivial intermediate angular momentum in CSF icsf.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_mcp
   implicit none
   !
   integer, intent(in) :: loc, icsf
   integer             :: i, jcup, locate, subshellX, ichop
   !
   locate = -1
   if (icsf == csf_a_No) then
      do  i = 1,two_csf%nwshells
         if (ichop(i,icsf) == 0) then
            locate    = locate + 1
            subshellX = two_csf%csf(1)%subshellX(i)
         end if
         !
         if (locate == loc) then
            jcup = subshellX + 1
            return
         end if
      end do
      print *, "subshellX(i) = ",two_csf%csf(1)%subshellX(1:two_csf%nwshells)
      print *, "loc, locate, subshellX = ",loc, locate, subshellX
      stop "jcup() in rabs_mcp_adaptation: program stop A."
   else if (icsf == csf_b_No) then
      do  i = 1,two_csf%nwshells
         !x if (ichop(i,icsf) == 0  .and.  jqs(3,i,icsf) /= 1) then
         if (ichop(i,icsf) == 0) then
            locate    = locate + 1
            subshellX = two_csf%csf(2)%subshellX(i)
         end if
         !
         if (locate == loc) then
            jcup = subshellX + 1
            return
         end if
      end do
      stop "jcup() in rabs_mcp_adaptation: program stop B."
   else
      stop "jcup() in rabs_mcp_adaptation: program stop C."
   end if
   !
end function jcup
!
!
function jqs(iwhich,isubsh,icsf)                                
!-----------------------------------------------------------------------------------------------------------------------
! jqs  is a subshell quantum number for subshell isubsh in configuration state function  icsf:  the seniority if  
! iwhich  is 1;  the quantum number w if  iwhich  is 2, and 2j+1 if  iwhich  is 3.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_mcp
   implicit none
   !
   integer, intent(in) :: iwhich, isubsh, icsf
   integer             :: jqs
   !
   if (isubsh > two_csf%nwshells) then
      stop "jqs() in rabs_mcp_adaptation: program stop A."
   else if (iwhich < 1   .or.   iwhich > 3) then
      print *, "iwhich = ",iwhich
      stop "jqs() in rabs_mcp_adaptation: program stop B."
   end if
   !
   if (icsf == csf_a_No) then
      select case(iwhich)
      case(1);   jqs = two_csf%csf(1)%seniority(isubsh)
      case(2);   jqs = 0
      case(3);   jqs = two_csf%csf(1)%subshellJ(isubsh) + 1
      end select
   else if (icsf == csf_b_No) then
      select case(iwhich)
      case(1);   jqs = two_csf%csf(2)%seniority(isubsh)
      case(2);   jqs = 0
      case(3);   jqs = two_csf%csf(2)%subshellJ(isubsh) + 1
      end select
   else
      stop "jqs() in rabs_mcp_adaptation: program stop C."
   end if
   !
end function jqs
!
!
subroutine mcp_initialize_grasp92()                                
!-----------------------------------------------------------------------------------------------------------------------
! Initializes all COMMON arrays from GRASP92 which are required to calculate angular coefficients by means of the MCP 
! and NJGRAF components.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_mcp
   implicit none
   !
   integer            :: i, icore, ncf, nw_grasp92, pntriq, np, nak, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6
   integer, parameter :: mfact=500
   real(kind=dp)      :: zero_92, half_92, tenth_92, one_92, two_92, three_92, ten_92, x, gam
   common/bcore/icore(149)
   common/cons/ zero_92,half_92,tenth_92,one_92,two_92,three_92,ten_92
   common/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6
   common/facts/gam(mfact)
   common/orb2/ ncf, nw_grasp92, pntriq
   common/orb4/ np(149), nak(149)
   !
   if (mcp_common_need_initialization) then
      zero_92 = 0.0_dp;   half_92 = 0.5_dp;   tenth_92 = 0.1_dp
      one_92  = 1.0_dp;   two_92  = 2.0_dp;   three_92 = 3.0_dp
      ten_92  = 10.0_dp
      !
      gam(1) = one_92;   gam(2) = one_92;   x = two_92
      do  i = 3,30
         gam(i) = gam(i-1)*x
         x      = x + one_92
      end do
      do  i = 1,30
         gam(i) = log(gam(i))
      end do
      x = 30_dp
      do  i = 31,mfact
         gam(i) = gam(i-1) + log(x)
         x = x + one_92
      end do
      !
      mcp_common_need_initialization = .false.
   end if
   !
   ibug1 = debug_ibug1;  ibug2 = debug_ibug2;  ibug3 = debug_ibug3
   ibug4 = debug_ibug4;  ibug5 = debug_ibug5;  ibug6 = debug_ibug6
   !
   !!x ibug1 = 1; ibug2 = 1; ibug3 = 1; ibug4 = 1; ibug5 = 1; ibug6 = 1
   ! 
   nw_grasp92 = two_csf%nwshells
   !
   if (nw_grasp92 > 149) then
      stop "mcp_initialize_grasp92() in rabs_mcp_adaptation: program stop A."
   end if
   !
   do  i = 1,nw_grasp92
      np(i)  = two_csf%subshell(i)%n
      nak(i) = two_csf%subshell(i)%kappa
   end do
   !
   do i = 1,nw_grasp92
      icore(i) = 0
      if (two_csf%csf(1)%occupation(1) == angular_momentum_j(nak(i)) + 1  .and. &
          two_csf%csf(2)%occupation(1) == angular_momentum_j(nak(i)) + 1)    then
         icore(i) = 1
      end if
   end do
   !
end subroutine mcp_initialize_grasp92
!
!
subroutine speak(ja,jb,ia1,ib1,ia2,ib2,k,x)                             
!-----------------------------------------------------------------------------------------------------------------------
! Output MCP coefficients and integral parameters to the vector mcp_list() which stores these coefficients for later use. 
! For IBUG1 = 1, it also prints these coefficients.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_dirac_orbital
   use rabs_mcp
   implicit none
   !
   integer, intent(in)       :: ja, jb, ia1, ib1, ia2, ib2, k
   real(kind=dp), intent(in) :: x
   !
   integer :: n
   !
   if (debug_ibug1 /= 0) then
      write (99,1) ja, jb, orbital_name(two_csf%subshell(ia1)%n,two_csf%subshell(ia1)%kappa),    &
                           orbital_name(two_csf%subshell(ib1)%n,two_csf%subshell(ib1)%kappa),    &
                           orbital_name(two_csf%subshell(ia2)%n,two_csf%subshell(ia2)%kappa),    &
                           orbital_name(two_csf%subshell(ib2)%n,two_csf%subshell(ib2)%kappa), k, x
    1 format(1x,2(1x,1i2), 4(1x,a4), 1x,i2,1x,es19.12)
   end if
   !
   ! Increment counter; ensure that arrays are of adequate size; reallocate if necessary
   number_of_mcp_coefficients = number_of_mcp_coefficients + 1
   if (number_of_mcp_coefficients > number_of_mcp_coefficients_max) then
      allocate( mcp_list_save(1:number_of_mcp_coefficients_max) )
      n = number_of_mcp_coefficients_max
      mcp_list_save(1:n) = mcp_list(1:n)
      deallocate( mcp_list )
      number_of_mcp_coefficients_max = number_of_mcp_coefficients_max + 1000
      allocate( mcp_list(1:number_of_mcp_coefficients_max) )
      mcp_list(1:n) = mcp_list_save(1:n)
      deallocate( mcp_list_save )
   end if
   !
   mcp_list(number_of_mcp_coefficients)%r  = ja
   mcp_list(number_of_mcp_coefficients)%s  = jb
   mcp_list(number_of_mcp_coefficients)%a  = ia1
   mcp_list(number_of_mcp_coefficients)%b  = ib1
   mcp_list(number_of_mcp_coefficients)%c  = ia2
   mcp_list(number_of_mcp_coefficients)%d  = ib2
   mcp_list(number_of_mcp_coefficients)%nu = k
   mcp_list(number_of_mcp_coefficients)%mu = 0
   mcp_list(number_of_mcp_coefficients)%V  = x
   !
end subroutine speak
!
!
subroutine talk (ja,jb,nu,ia,ib,ic,id,itype,x)                            
!-----------------------------------------------------------------------------------------------------------------------
! Output MCBP coefficients and integral parameters to the vector mcbp_list() which stores these coefficients for later 
! use. For IBUG1 = 1, it also prints these coefficients.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_dirac_orbital
   use rabs_mcp
   implicit none
   !
   integer, intent(in)       :: ja, jb, ia, ib, ic, id, itype, nu
   real(kind=dp), intent(in) :: x
   !
   integer :: n
   !
   if (debug_ibug1 /= 0) then
      write (99,1) ja, jb, orbital_name(two_csf%subshell(ia)%n,two_csf%subshell(ia)%kappa),            &
                           orbital_name(two_csf%subshell(ib)%n,two_csf%subshell(ib)%kappa),            &
                           orbital_name(two_csf%subshell(ic)%n,two_csf%subshell(ic)%kappa),            &
                           orbital_name(two_csf%subshell(id)%n,two_csf%subshell(id)%kappa), nu, itype, x
    1 format(1x,2(1x,1i2), 4(1x,a4), 1x,i2,1x,i2,1x,es19.12)
   end if
   !
   ! Increment counter; ensure that arrays are of adequate size; reallocate if necessary
   number_of_mcbp_coefficients = number_of_mcbp_coefficients + 1
   if (number_of_mcbp_coefficients > number_of_mcbp_coefficients_max) then
      allocate( mcp_list_save(1:number_of_mcbp_coefficients_max) )
      n = number_of_mcbp_coefficients_max
      mcp_list_save(1:n) = mcbp_list(1:n)
      deallocate( mcbp_list )
      number_of_mcbp_coefficients_max = number_of_mcbp_coefficients_max + 1000
      allocate( mcbp_list(1:number_of_mcbp_coefficients_max) )
      mcbp_list(1:n) = mcp_list_save(1:n)
      deallocate( mcp_list_save )
   end if
   !
   mcp_list(number_of_mcp_coefficients)%r  = ja
   mcp_list(number_of_mcp_coefficients)%s  = jb
   mcp_list(number_of_mcp_coefficients)%a  = ia
   mcp_list(number_of_mcp_coefficients)%b  = ib
   mcp_list(number_of_mcp_coefficients)%c  = ic
   mcp_list(number_of_mcp_coefficients)%d  = id
   mcp_list(number_of_mcp_coefficients)%nu = nu
   mcp_list(number_of_mcp_coefficients)%mu = itype
   mcp_list(number_of_mcp_coefficients)%V  = x
   !
end subroutine talk
!
