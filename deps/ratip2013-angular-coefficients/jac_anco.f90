module jac_anco

    use rabs_anco
    use rabs_mcp

    implicit none

    public :: jac_anco_load_csl
    ! Saves the CSL specified the arguments in csf_set. No consistency checks are performed.
    public :: jac_anco_calculate_csf_pair
    ! Calculate angular coefficients between the r-th and s-th CSF of csf_set
    public :: jac_anco_calculate_csf_pair_1p
    ! Calculate angular coefficients for a non-scalar one-particle operator between the r-th and s-th CSF of csf_set
    public :: jac_mct_generate_coefficients
    ! Calculate MCT coefficients between the rl..ru-th and sl..su-th CSFs of csf_set

    type(csf_basis)  :: csf_set

    integer, bind(c) :: no_anco_T_list = 0, &
                        no_anco_V_list = 0, &
                        no_mct_list    = 0

    contains

    subroutine jac_anco_load_csl(nocsf, nwshells, nwcore,               &
                                 number_of_electrons, subshell,         &
                                 totalJs, parities, occupations,        &
                                 seniorities, subshellJs, subshellXs)   &
                                                                  bind(c)
    !--------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------
        !
        implicit none
        integer, intent(in) :: nocsf, nwshells, nwcore, number_of_electrons
        type(nkappa), dimension(nwshells), intent(in)   :: subshell
        integer(kind=i1b), dimension(nocsf), intent(in) :: totalJs
        character, dimension(nocsf), intent(in)         :: parities
        integer(kind=i1b), dimension(nwshells,nocsf), intent(in) :: occupations, seniorities, subshellJs, subshellXs
        !
        integer :: i,j
        !
        call deallocate_csf_basis(csf_set)
        !
        csf_set%nocsf               = nocsf
        csf_set%nwshells            = nwshells
        csf_set%nwcore              = nwcore
        csf_set%number_of_electrons = number_of_electrons
        !
        allocate(csf_set%subshell(1:nwshells))
        csf_set%subshell = subshell
        !
        allocate(csf_set%csf(1:nocsf))
        csf_set%csf%totalJ = totalJs
        csf_set%csf%parity = parities
        do i = 1,nocsf
            allocate(csf_set%csf(i)%occupation(1:nwshells), csf_set%csf(i)%seniority(1:nwshells), &
                     csf_set%csf(i)%subshellJ(1:nwshells),  csf_set%csf(i)%subshellX(1:nwshells))
            csf_set%csf(i)%occupation = occupations(:,i)
            csf_set%csf(i)%seniority  = seniorities(:,i)
            csf_set%csf(i)%subshellJ  = subshellJs(:,i)
            csf_set%csf(i)%subshellX  = subshellXs(:,i)
        end do
        !
    end subroutine jac_anco_load_csl

    subroutine jac_anco_calculate_csf_pair(r, s)                  bind(c)
    !--------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------
        !
        implicit none
        integer, intent(in) :: r,s
        !
        call anco_calculate_csf_pair(csf_set, r, s, no_anco_T_list, no_anco_V_list)
        !
    end subroutine jac_anco_calculate_csf_pair

    subroutine jac_anco_calculate_csf_pair_1p(nu, r, s)           bind(c)
    !--------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------
        !
        implicit none
        integer, intent(in) :: nu,r,s
        !
        call anco_calculate_csf_pair_1p(nu, csf_set, r, s, no_anco_T_list)
        !
    end subroutine jac_anco_calculate_csf_pair_1p

    subroutine jac_mct_generate_coefficients(rl, ru, sl, su, iopar,     &
                                             rank)                bind(c)
    !--------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------
        !
        implicit none
        integer, intent(in) :: rl, ru, sl, su, iopar, rank
        !
        call mct_generate_coefficients(rl, ru, sl, su, csf_set, iopar, rank)
        !
    end subroutine jac_mct_generate_coefficients

end module jac_anco
