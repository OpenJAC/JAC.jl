module rabs_recoupling
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module provides a set of procedures which support the analytic integration over the angular coordinates for 
! matrix elements with symmetry-adapted subshell states. This includes a collection of routines which are concerned with 
! the calculation of recoupling matrices in jj--coupling (see in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. 
! Opt. Phys, Vol 30 3747 section 4).
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   use rabs_rcfp
   implicit none
   !
   private :: recoupling_C_1
                 ! Returns the value of the coefficient C_1 from the paper of G. Gaigalas et al., 1997 J. Phys. B: 
                 ! At. Mol. Opt. Phys, Vol 30 3747. The expression of C_1 is given in Eq. (15).
   private :: recoupling_C_2
                 ! Returns the value of the coefficient C_2 from the paper of G. Gaigalas et al., 1997 J. Phys. B: 
                 ! At. Mol. Opt. Phys, Vol 30 3747. The expression of C_2 is given in Eq. (16).
   private :: recoupling_C_3
                 ! Returns the value of the coefficient C_3 from the paper of G. Gaigalas et al., 1997 J. Phys. B: 
                 ! At. Mol. Opt. Phys, Vol 30 3747. The expression of C_3 is given in Eq. (17).
   private :: recoupling_C_4
                 ! Returns the value of the coefficient C_4 from the paper of G. Gaigalas et al., 1997 J. Phys. B: 
                 ! At. Mol. Opt. Phys, Vol 30 3756. The expression of C_4 is given in Eq. (21).
   private :: recoupling_C_5
                 ! Returns the value of the coefficient C_5 from the paper of G. Gaigalas et al., 1997 J. Phys. B: 
                 ! At. Mol. Opt. Phys, Vol 30 3747. The expression of C_5 is given in Eq. (23).
   public  :: recoupling_matrix_check
                 ! Checks the angular momentum selection rules for the recoupling coefficients.
   public  :: recoupling_matrix_check_nonscal
                 ! Checks the angular momentum selection rules for the recoupling coefficients of non scalar operator.
   public  :: recoupling_matrix_1p_shells
                 ! Calculates the recoupling matrix for a non scalar operator in the case of one interacting shells.
   public  :: recoupling_matrix_2_shells
                 ! Calculates the recoupling matrix for a scalar operator in the case of two interacting shells.
   public  :: recoupling_matrix_2p_shells
                 ! Calculates the recoupling matrix for a non scalar operator in the case of two interacting shells.
   public  :: recoupling_matrix_3_shells
                 ! Calculates the recoupling matrix for a scalar operator in the case of three interacting shells.
   private :: recoupling_matrix_3_ordered
                 ! Auxiliarity routine for the procedure recoupling_matrix_3_shells
   public  :: recoupling_matrix_4_shells
                 ! Calculates the recoupling matrix for a scalar operator in the case of four interacting shells.
contains
   !
   function recoupling_C_1(csf_r,csf_s,rank,shell,verify,delta_J)  result(coeff)
   !--------------------------------------------------------------------
   ! Returns the coefficient C_1 from the paper of G. Gaigalas et al., 
   ! 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (15). 
   ! The phase factor is given in the table 2 in this reference.
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, shell
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer :: J, csf_s_J, csf_s_X, csf_r_J, csf_r_X
      !
      coeff = zero
      csf_r_J = csf_r%subshellJ(shell);   csf_s_J = csf_s%subshellJ(shell)
      !
      select case(shell)
      case(1)
         csf_r_X = csf_r%subshellX(2);   csf_s_X = csf_s%subshellX(2) 
         J = csf_r%subshellJ(2)
      case(2)
         csf_r_X = csf_r%subshellX(2);   csf_s_X = csf_s%subshellX(2) 
         J = csf_r%subshellJ(1)
      case default
         csf_r_X = csf_r%subshellX(shell);   csf_s_X = csf_s%subshellX(shell) 
         J = csf_r%subshellX(shell-1) 
      end select
      !
      if (verify) then
         delta_J = wigner_6j_triangle(rank,csf_s_J,csf_r_J,J,csf_r_X,csf_s_X)
      else
         delta_J = 1
         coeff = wigner_6j_symbol(rank,csf_s_J,csf_r_J,J,csf_r_X,csf_s_X)
         coeff = coeff * sqrt((csf_r_J + one) * (csf_s_X + one))
         if (mod(J+csf_r_X+csf_s_J+rank ,4) /= 0) coeff = - coeff
         if (shell == 1) then
            if (mod(csf_r_J+csf_s_J+2*J-csf_r_X-csf_s_X, 4) /= 0) coeff = -coeff
         end if
      end if
      !
   end function recoupling_C_1  
   !
   !
   function recoupling_C_2(csf_r,csf_s,rank,shell_1,shell_2,verify,delta_J)    &
                                                                   result(coeff)
   !--------------------------------------------------------------------
   ! Returns the coefficient C_2 from the paper of G. Gaigalas et al., 
   ! 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (16). 
   ! The phase factor is given in the table 2 in this reference.
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, shell_1, shell_2
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer          :: run, csf_J, csf_r_X, csf_s_X, csf_r_X_1, csf_s_X_1
      real(kind=dp)    :: coeff_run
      !
      coeff = one
      !
      if (shell_1 == 1) then
         run = 3
      else
         run = shell_1 + 1
      end if
      !
      delta_J = 1
    1 if (run < shell_2) then
         csf_J     = csf_r%subshellJ(run)
         csf_r_X   = csf_r%subshellX(run-1);  csf_s_X   = csf_s%subshellX(run-1)
         csf_r_X_1 = csf_r%subshellX(run);    csf_s_X_1 = csf_s%subshellX(run)
         if (verify) then
            delta_J = wigner_6j_triangle(rank,csf_s_X,csf_r_X,csf_J,           &
                                                           csf_r_X_1,csf_s_X_1)
            if (delta_J == 0) return
         else
            coeff_run = wigner_6j_symbol(rank,csf_s_X,csf_r_X,csf_J,           &
                                                            csf_r_X_1,csf_s_X_1)
            coeff_run = coeff_run * sqrt((csf_r_X + one) * (csf_s_X_1 + one))
            if (mod(rank+csf_J+csf_r_X+csf_s_X_1,4) /= 0) coeff_run = -coeff_run
            coeff = coeff * coeff_run
         end if
         run = run + 1
         go to 1
      end if
      !
   end function recoupling_C_2
   !
   !
   function recoupling_C_3(csf_r,csf_s,rank,shell,nwshells,verify,delta_J)      &
                                                                  result(coeff)
   !--------------------------------------------------------------------
   ! Returns the coefficient C_3 from the paper of G. Gaigalas et al., 
   ! 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (17). 
   ! The phase factor is given in the table 3 in this reference.
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, shell, nwshells
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer :: T_r, T_s, T1_r, T1_s, JI
      !
      coeff = zero
      T1_r = csf_r%subshellX(nwshells);       T1_s = csf_s%subshellX(nwshells)
      if (shell == nwshells) then
         T_r = csf_r%subshellJ(nwshells);     T_s = csf_s%subshellJ(nwshells)
         JI  = csf_r%subshellX(nwshells-1)
      else
         T_r = csf_r%subshellX(nwshells-1);   T_s = csf_s%subshellX(nwshells-1)
         JI  = csf_r%subshellJ(nwshells)
      end if
      if (verify) then
         delta_J = wigner_6j_triangle(rank,T_s,T_r,JI,T1_r,T1_s)
      else
         coeff = wigner_6j_symbol(rank,T_s,T_r,JI,T1_r,T1_s)
         coeff = coeff * sqrt((T_r + one) * (T1_s + one))
         if (mod(rank+JI+T_s+T1_r,4) /= 0) coeff = - coeff
         delta_J = 1
         if (shell == nwshells) return
         if (mod(T_r+T_s-T1_s-T1_r+2*JI, 4) /= 0) coeff = -coeff
      end if
      !
   end function recoupling_C_3  
   !
   !
   function recoupling_C_4(csf_r,csf_s,rank_1,rank_2,rank,shell_1,shell_2,     &
                                              verify,delta_J)      result(coeff)
   !--------------------------------------------------------------------
   ! Returns the coefficient C_4 from the paper of G. Gaigalas et al., 
   ! 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (23). 
   ! The phase factor is given in the table 4 in this reference.
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, rank_1, rank_2, shell_1, shell_2
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer   :: csf_s_J_2, csf_r_J_2, csf_s_X_1, csf_r_X_1, csf_s_X, csf_r_X
      !
      coeff = zero
      !
      csf_r_J_2 = csf_r%subshellJ(shell_2); csf_s_J_2 = csf_s%subshellJ(shell_2)
      csf_r_X   = csf_r%subshellX(shell_2); csf_s_X   = csf_s%subshellX(shell_2)
      !
      if (shell_1 == 1  .and.  shell_2 == 2) then
         csf_r_X_1 = csf_r%subshellJ(shell_1)
         csf_s_X_1 = csf_s%subshellJ(shell_1)
      else
         csf_r_X_1 = csf_r%subshellX(shell_2-1)
         csf_s_X_1 = csf_s%subshellX(shell_2-1)
      end if
      !
      if (verify) then
         delta_J = wigner_9j_triangle(csf_s_X_1,rank_1,csf_r_X_1,              &
                                      csf_s_J_2,rank_2,csf_r_J_2,              &
                                      csf_s_X  ,rank  ,csf_r_X)
      else
         delta_J = 1
         coeff = wigner_9j_symbol(csf_s_X_1,rank_1,csf_r_X_1,                  &
                                  csf_s_J_2,rank_2,csf_r_J_2,                  &
                                  csf_s_X  ,rank  ,csf_r_X  ,.true.)
         coeff = coeff * sqrt((csf_r_X_1 + one)*(rank + one)*(csf_r_J_2 + one) &
                                                            *(csf_s_X   + one))
      end if
      !
   end function recoupling_C_4
   !
   !
   function recoupling_C_5(csf_r,csf_s,rank,shell_1,shell_2,verify,delta_J)    &
                                                                   result(coeff)
   !--------------------------------------------------------------------
   ! Returns the coefficient C_5 from the paper of G. Gaigalas et al., 
   ! 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (23). 
   ! The phase factor is given in the table 5 in this reference.
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, shell_1, shell_2
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer    :: csf_s_J, csf_r_J, csf_s_J_2, csf_r_J_2, csf_X_2
      !
      coeff = zero
      !
      csf_r_J_2 = csf_r%subshellJ(shell_2); csf_s_J_2 = csf_s%subshellJ(shell_2)
      csf_X_2   = csf_r%subshellX(shell_2)
      !
      if (shell_1 == 1  .and.  shell_2 == 2) then
         csf_r_J = csf_r%subshellJ(shell_1); csf_s_J = csf_s%subshellJ(shell_1)
      else
         csf_r_J = csf_r%subshellX(shell_2-1)
         csf_s_J = csf_s%subshellX(shell_2-1) 
      end if
      !
      if (verify) then
         delta_J = wigner_6j_triangle(rank,csf_s_J_2,csf_r_J_2,csf_X_2,        &
                                                               csf_r_J,csf_s_J)
      else
         delta_J = 1
         coeff = wigner_6j_symbol(rank,csf_s_J_2,csf_r_J_2,csf_X_2,            &
                                                               csf_r_J,csf_s_J)
         coeff = coeff * sqrt((csf_s_J_2 + one) * (csf_r_J + one))
         if (mod(csf_X_2+csf_s_J+csf_r_J_2+rank ,4) /= 0) coeff = - coeff
      end if
      !
   end function recoupling_C_5  
   !
   !
   function recoupling_matrix_check(csf_r,csf_s,shell_1,shell_2,shell_3,       &
                                   shell_4,nwshells,nwcore,switch) result(coeff)
   !--------------------------------------------------------------------
   ! Checks the angular momentum selection rules for the recoupling
   ! coefficients in the cases of one, two, three or four interacting
   ! shells.
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: shell_1, shell_2, shell_3, shell_4
      integer, intent(in)           :: nwshells, nwcore, switch
      real(kind=dp)                 :: coeff
      !
      integer     :: i, shell
      !
      coeff = one
      if (nwshells == 1) return
      !
      if (switch == 0) then
         do i = nwcore + 1, nwshells
            if (i > 1) then
               if (csf_r%subshellX(i) /= csf_s%subshellX(i)) coeff = zero
            end if
            if (csf_r%subshellJ(i) /= csf_s%subshellJ(i)) coeff = zero
            if (csf_r%seniority(i) /= csf_s%seniority(i)) coeff = zero
            if (coeff == zero) return
         end do
      else if (shell_1 == 1  .and.  shell_2 == 2) then
         do i = nwcore + 1, nwshells
            if (i > 1) then
               if (csf_r%subshellX(i) /= csf_s%subshellX(i)) coeff = zero
            end if
            if (i > 2) then 
               if (csf_r%subshellJ(i) /= csf_s%subshellJ(i)) coeff = zero
               if (csf_r%seniority(i) /= csf_s%seniority(i)) coeff = zero
            end if
            if (coeff == zero) return
         end do
      else
         if (shell_1 == 2) then
            shell = 1
         else
            shell = shell_1 
         end if
         do i = nwcore + 1, nwshells
            if (i < shell  .or.  i >= shell_2) then
               if(i > 1) then
                  if (csf_r%subshellX(i) /= csf_s%subshellX(i)) coeff = zero
               end if
            end if
            if(i == shell_1) then
            else if (i == shell_2) then
            else if (switch == 2 .and. i == shell_3) then
            else if (switch == 3 .and. i == shell_3) then
            else if (switch == 3 .and. i == shell_4) then
            else
               if (csf_r%subshellJ(i) /= csf_s%subshellJ(i)) coeff = zero
               if (csf_r%seniority(i) /= csf_s%seniority(i)) coeff = zero
            end if
            if (coeff == zero) return
            end do
      end if 
      !
   end function recoupling_matrix_check
   !
   !
   function recoupling_matrix_check_nonscal(csf_r,csf_s,rank,shell_1,       &
                                     shell_2,nwshells,nwcore) result(coeff)
   !--------------------------------------------------------------------
   ! Checks the angular momentum selection rules for the recoupling
   ! coefficients in the cases of one ore two interacting shells
   ! for non scalar operatro
   !--------------------------------------------------------------------
      !
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: shell_1, shell_2
      integer, intent(in)           :: nwshells, nwcore, rank
      real(kind=dp)                 :: coeff
      !
      integer     :: i
      !
      coeff = zero
      if (triangle(csf_r%totalJ+1,rank+1,csf_s%totalJ+1) == 0) return
      coeff = one
      if (nwshells == 1) return
      !
      do i = nwcore + 1, nwshells
         if (i /= shell_1) then
            if (i /= shell_2) then
               if (csf_r%subshellJ(i) /= csf_s%subshellJ(i)) coeff = zero
               if (csf_r%seniority(i) /= csf_s%seniority(i)) coeff = zero
            end if
         end if
         if (coeff == zero) return
      end do
      if (nwshells <= 2) return
      if (nwshells-shell_2 > 0) then
         do i=shell_2,nwshells
            if (triangle(csf_r%subshellX(i)+1,rank+1,csf_s%subshellX(i)+1) &
                == 0) coeff = zero
         end do
      end if
      if (coeff == zero) return
      if (shell_1 <= nwcore) return
      do i = nwcore, shell_1-1
         if (i == 0) cycle
         if (csf_r%subshellX(i) /= csf_s%subshellX(i)) coeff = zero
      end do
      !
   end function recoupling_matrix_check_nonscal
   !
   !
   function recoupling_matrix_1p_shells(csf_r,csf_s,rank,shell_1,verify,       &
                                       nwshells,nwcore,delta_J)   result(coeff)
   !--------------------------------------------------------------------
   ! Calculates the recoupling matrix for a scalar operator in the case
   ! of two interacting shells. See G. Gaigalas et al., 1997 J. Phys. B:
   ! At. Mol. Opt. Phys, Vol 30 3747, Eq. (14) (p. 3756). 
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, shell_1, nwshells, nwcore
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer                       :: run, run_correction
      !
      delta_J = 1
      coeff = one / sqrt((csf_r%subshellJ(shell_1) + one))
      if (nwshells == 1 ) return
      delta_J = 0
      if ((nwshells) /= 2) then
         coeff = coeff * recoupling_C_3                                       &
                         (csf_r,csf_s,rank,shell_1,nwshells,verify,delta_J)
         if (delta_J == 0) return 
         if (shell_1 == nwshells) return
         delta_J = 0
      end if
      coeff = coeff * recoupling_C_1(csf_r,csf_s,rank,shell_1,verify,delta_J)
      if (delta_J == 0) return
      run_correction = nwshells - shell_1
      if(shell_1 == 1) run_correction = run_correction - 1
      if(run_correction <= 1) return
      delta_J = 0
      coeff = coeff * recoupling_C_2                                          &
                            (csf_r,csf_s,rank,shell_1,nwshells,verify,delta_J)
      !
   end function recoupling_matrix_1p_shells
   !
   !
   function recoupling_matrix_2_shells(csf_r,csf_s,rank,shell_1,shell_2,verify,&
                                                        delta_J)   result(coeff)
   !--------------------------------------------------------------------
   ! Calculates the recoupling matrix for a scalar operator in the case
   ! of two interacting shells. See G. Gaigalas et al., 1997 J. Phys. B:
   ! At. Mol. Opt. Phys, Vol 30 3747, Eq. (22) (p. 3756). 
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank, shell_1, shell_2
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      delta_J = 0
      coeff = one / sqrt((csf_r%subshellJ(shell_1) + one) *                    &
                         (csf_r%subshellJ(shell_2) + one))
      if (rank /= 0 .or. verify) then
         coeff = coeff *                                                       &
                 recoupling_C_5(csf_r,csf_s,rank,shell_1,shell_2,verify,delta_J)
         if (delta_J == 0) return
         coeff = coeff * sqrt((csf_r%subshellJ(shell_2) + one) /               &
                              ((rank + one) * (csf_s%subshellJ(shell_2) + one)))
         if (shell_1 == 1 .and. shell_2 == 2) return
         coeff = coeff * recoupling_C_1(csf_r,csf_s,rank,shell_1,verify,delta_J)
         if (delta_J == 0) return
         if (shell_1 == 1) then
            if (shell_2 - 1 - shell_1 <= 1) return
         else
            if (shell_2 - shell_1     <= 1) return
         end if
         coeff = coeff *                                                       &
                 recoupling_C_2(csf_r,csf_s,rank,shell_1,shell_2,verify,delta_J)
      else
         delta_J = 1
      end if
      !
   end function recoupling_matrix_2_shells
   !
   !
   function recoupling_matrix_2p_shells(csf_r,csf_s,rank_1,rank_2,rank,shell_1,&
                         shell_2,verify,nwshells,nwcore,delta_J)   result(coeff)
   !--------------------------------------------------------------------
   ! Calculates the recoupling matrix for a non scalar operator in the case
   ! of two interacting shells. See G. Gaigalas et al., 1997 J. Phys. B:
   ! At. Mol. Opt. Phys, Vol 30 3747, Eq. (19) (p. 3756). 
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank_1, rank_2, rank, shell_1, shell_2
      integer, intent(in)           :: nwshells, nwcore
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      integer                       :: run, run_correction
      !
      delta_J = 1
      coeff = one / sqrt((csf_r%subshellJ(shell_1) + one)                      &
                        *(csf_r%subshellJ(shell_2) + one))
      if ((nwshells - shell_2) > 1) then
         delta_J = 0
         coeff = coeff * recoupling_C_2                                        &
                            (csf_r,csf_s,rank,shell_2,nwshells,verify,delta_J)
         if (delta_J == 0) return
         delta_J = 0
      end if
      if (shell_2 /= nwshells) then
         coeff = coeff * recoupling_C_3                                        &
                         (csf_r,csf_s,rank,shell_2,nwshells,verify,delta_J)
         if (delta_J == 0) return 
         delta_J = 0
      end if
      coeff = coeff * recoupling_C_4(csf_r,csf_s,rank_1,rank_2,rank,shell_1,   &
                                                         shell_2,verify,delta_J)
      if (delta_J == 0) return 
      if (shell_1 == 1 .and. shell_2 == 2) return
      delta_J = 0
      coeff = coeff * recoupling_C_1(csf_r,csf_s,rank_1,shell_1,verify,delta_J)
      if (delta_J == 0) return
      run_correction = shell_2 - shell_1
      if(shell_1 == 1) run_correction = run_correction - 1
      if(run_correction <= 1) return
      delta_J = 0
      coeff = coeff * recoupling_C_2                                          &
                            (csf_r,csf_s,rank_1,shell_1,shell_2,verify,delta_J)
      !
   end function recoupling_matrix_2p_shells
   !
   !
   function recoupling_matrix_3_shells(csf_r,csf_s,rank_1,rank_2,rank_3,       &
                          shell_1,shell_2,shell_3,verify,delta_J)  result(coeff)
   !--------------------------------------------------------------------
   ! Calculates the recoupling matrix for a scalar operator in the case
   ! of three interacting shells. See G. Gaigalas et al., 1997 J. Phys. B:
   ! At. Mol. Opt. Phys, Vol 30 3747, Eq. (26) (p. 3758).
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank_1, rank_2, rank_3
      integer, intent(in)           :: shell_1, shell_2, shell_3
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      if ((shell_3 > shell_1)  .and. (shell_3 > shell_2)) then
         if (shell_1 < shell_2) then
            coeff = recoupling_matrix_3_ordered(csf_r,csf_s,rank_1,rank_2,     &
                                  rank_3,shell_1,shell_2,shell_3,verify,delta_J)
         else if (shell_1 > shell_2) then
            coeff = recoupling_matrix_3_ordered(csf_r,csf_s,rank_2,rank_1,     &
                                  rank_3,shell_2,shell_1,shell_3,verify,delta_J)
            if (mod(rank_1 + rank_2 - rank_3,4) /= 0) coeff = -coeff
         else if (rabs_use_stop) then
            stop "recoupling_matrix_3_shells(): program stop A."
         end if
      else if ((shell_3 < shell_1)  .and. (shell_3 < shell_2)) then
         if (shell_1 < shell_2) then
            coeff = recoupling_matrix_3_ordered(csf_r,csf_s,rank_3,rank_1,     &
                                  rank_2,shell_3,shell_1,shell_2,verify,delta_J)
            if (mod(rank_3,2) /= 0) coeff = -coeff
         else if (shell_1 > shell_2) then
            coeff = recoupling_matrix_3_ordered(csf_r,csf_s,rank_3,rank_2,     &
                                  rank_1,shell_3,shell_2,shell_1,verify,delta_J)
            if (mod(rank_1 + rank_2 + rank_3,4) /= 0) coeff = -coeff
         else if (rabs_use_stop) then
            stop "recoupling_matrix_3_shells(): program stop B."
         end if
      else
         if (shell_1 < shell_2) then
            coeff = recoupling_matrix_3_ordered(csf_r,csf_s,rank_1,rank_3,     &
                                  rank_2,shell_1,shell_3,shell_2,verify,delta_J)
            if (mod(rank_1 - rank_2 - rank_3,4) /= 0) coeff = -coeff
         else if (shell_1 > shell_2) then
            coeff = recoupling_matrix_3_ordered(csf_r,csf_s,rank_2,rank_3,     &
                                  rank_1,shell_2,shell_3,shell_1,verify,delta_J)
            if (mod(rank_1,2) /= 0) coeff = -coeff
         else if (rabs_use_stop) then
            stop "recoupling_matrix_3_shell(): program stop C."
         end if
      end if
      !
   end function recoupling_matrix_3_shells
   !
   !
   function recoupling_matrix_3_ordered(csf_r,csf_s,rank_1,rank_2,rank,shell_1,&
                                 shell_2,shell_3,verify,delta_J)   result(coeff)
   !--------------------------------------------------------------------
   ! Auxiliary routine for recoupling_matrix_3_shells.
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank_1, rank_2, rank
      integer, intent(in)           :: shell_1, shell_2, shell_3
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      coeff = one / sqrt( (csf_r%subshellJ(shell_1) + one) *                   &
              (csf_r%subshellJ(shell_2) + one) * (rank + one) *                &
              (csf_s%subshellJ(shell_3) + one))
      if (shell_3 - shell_2 > 1) then
         coeff = coeff *                                                       &
                 recoupling_C_2(csf_r,csf_s,rank,shell_2,shell_3,verify,delta_J)
         if (delta_J == 0) return 
      end if
      coeff = coeff *                                                          &
                 recoupling_C_5(csf_r,csf_s,rank,shell_1,shell_3,verify,delta_J)
      if (delta_J == 0) return 
      coeff = coeff * recoupling_C_4(csf_r,csf_s,rank_1,rank_2,rank,shell_1,   &
                                                         shell_2,verify,delta_J)
      if (delta_J == 0) return 
      if (shell_1 == 1 .and. shell_2 == 2) return
      coeff = coeff * recoupling_C_1(csf_r,csf_s,rank_1,shell_1,verify,delta_J)
      if (delta_J == 0) return 
      if (shell_1 == 1) then
         if (shell_2 - 1 - shell_1 <= 1) return
      else
         if (shell_2 - shell_1     <= 1) return
      end if
      coeff = coeff *                                                          &
              recoupling_C_2(csf_r,csf_s,rank_1,shell_1,shell_2,verify,delta_J)
      !
   end function recoupling_matrix_3_ordered
   !
   !
   function recoupling_matrix_4_shells(csf_r,csf_s,rank_1,rank_2,rank_3,rank_4,&
            rank,shell_1,shell_2,shell_3,shell_4,verify,delta_J)   result(coeff)
   !--------------------------------------------------------------------
   ! Calculates the recoupling matrix for a scalar operator in the case
   ! of four interacting shells. See G. Gaigalas et al., 1997 J. Phys. B:
   ! At. Mol. Opt. Phys, Vol 30 3747, Eq. (33) (p. 3760).
   !--------------------------------------------------------------------
      !
      logical, intent(in)           :: verify
      type(cs_function), intent(in) :: csf_r, csf_s
      integer, intent(in)           :: rank_1, rank_2, rank_3, rank_4, rank
      integer, intent(in)           :: shell_1, shell_2, shell_3, shell_4
      integer, intent(out)          :: delta_J
      real(kind=dp)                 :: coeff
      !
      coeff = one / sqrt( (csf_r%subshellJ(shell_1) + one) *                   &
              (csf_r%subshellJ(shell_2) + one)*(csf_r%subshellJ(shell_3) + one)&
                            * (rank_4 + one)*(csf_s%subshellJ(shell_4) + one) )
      if (shell_3 - shell_2 > 1) then
         coeff = coeff *                                                       &
               recoupling_C_2(csf_r,csf_s,rank,shell_2,shell_3,verify,delta_J)
         if (delta_J == 0) return 
      end if 
      if (shell_4 - shell_3 > 1) then
         coeff=coeff *                                                         &
               recoupling_C_2(csf_r,csf_s,rank_4,shell_3,shell_4,verify,delta_J)
         if (delta_J == 0) return 
      end if
      coeff = coeff *                                                          &
              recoupling_C_5(csf_r,csf_s,rank_4,shell_1,shell_4,verify,delta_J)
      if (delta_J == 0) return 
      coeff = coeff * recoupling_C_4(csf_r,csf_s,rank_1,rank_2,rank,shell_1,   &
                                                         shell_2,verify,delta_J)
      if (delta_J == 0) return  
      coeff = coeff * recoupling_C_4(csf_r,csf_s,rank,rank_3,rank_4,shell_2,   &
                                                         shell_3,verify,delta_J)
      if (delta_J == 0) return  
      if (shell_1 == 1 .and. shell_2 == 2) return
      coeff = coeff * recoupling_C_1(csf_r,csf_s,rank_1,shell_1,verify,delta_J)
      if (delta_J == 0) return
      if (shell_1 == 1) then
         if (shell_2 - 1 - shell_1 <= 1) return
      else
         if (shell_2 - shell_1     <= 1) return
      end if
      coeff = coeff *                                                          &
              recoupling_C_2(csf_r,csf_s,rank_1,shell_1,shell_2,verify,delta_J)
      !
   end function recoupling_matrix_4_shells
   !
end module rabs_recoupling
