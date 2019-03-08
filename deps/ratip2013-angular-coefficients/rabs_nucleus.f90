module rabs_nucleus
!
!***** February 2012, wide *****
!-----------------------------------------------------------------------------------------------------------------------
! This module stores information about the nucleus and how it is modeled in the calculation. It also allows for to read 
! this information from a GRASP92 .iso file.
!-----------------------------------------------------------------------------------------------------------------------
   !
   use rabs_constant
   implicit none
   private
   !
   public  :: gamma_r               
                 ! Return the relativistic gamma factor for a given symmetry kappa.
   public  :: nuclear_charge_density     
                 ! Return the nuclear charge density rho(r) at r0.
   public  :: nuclear_charge_density_prime     
                 ! Return the derivative of the nuclear charge density d rho(r) / dr at r0.
   public  :: nuclear_extent        
                 ! Return a nuclear radius R_nuc which includes at least 99 percent of charge.
   public  :: nuclear_potential     
                 ! Return the nuclear potential at r0.
   public  :: nuclear_potential_es 
                 ! Utility subroutine for nuclear_potential.
   public  :: nuclear_radial_points
                 ! Determines the number of radial points inside of R_nuc.
   !
   real(kind=dp), public     :: nuclear_charge
   real(kind=dp), public     :: nuclear_mass
   real(kind=dp), public     :: nuclear_spin
   real(kind=dp), public     :: nuclear_dipole_moment
   real(kind=dp), public     :: nuclear_quadrupole_moment
   real(kind=dp), public     :: atomic_mass
   real(kind=dp), public     :: fermi_a_parameter
   real(kind=dp), public     :: fermi_c_parameter
   !
   real(kind=dp), public     :: nuclear_screening_lambda    = zero
   real(kind=dp), public     :: nuclear_ion_sphere_Nb       = zero
   real(kind=dp), public     :: nuclear_ion_sphere_R0       = zero
   real(kind=dp), public     :: nuclear_debye_sphere_R0     = zero
   !
   character(len=12), public :: nuclear_model
   !
contains
   !
   !
   function gamma_r(kappa)   result(gamma)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the relativistic gamma factor for a given symmetry kappa.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: kappa
      real(kind=dp)       :: gamma
      !
      gamma = kappa * kappa - nuclear_charge * nuclear_charge / (c * c)
      if (rabs_use_stop   .and.  gamma < zero) then
         stop "gamma_r(): program stop A."
      end if
      gamma = sqrt(gamma)
      !
   end function gamma_r
   !
   !
   function nuclear_charge_density(r0)                   result(density)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a nuclear charge density at a given radius r0.
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: r0
      real(kind=dp)             :: density, Rnuc
      !
      select case(nuclear_model)
      case("point")
         density = zero
      case("uniform")
         Rnuc = nuclear_extent()
	 if (r0 <= Rnuc) then
	    density = three / (two*two*pi*(Rnuc**3))
	 !!! else if (r0 < Rnuc + 0.0001_dp) then
	 !!!    density = three / (two*two*pi*(Rnuc**3)) &
	 !!!              *(one- (r0-Rnuc)/0.0001_dp)
	 else
	    density = zero
	 end if
      case("fermi")
         density = one /(exp( (r0-fermi_c_parameter)/fermi_a_parameter ) + one)
      case default
         stop "nuclear_charge_density(): program stop A."
      end select
      !
   end function nuclear_charge_density
   !
   !
   function nuclear_charge_density_prime(r0)              result(rho_prime)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a nuclear charge density at a given radius r0.
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: r0
      real(kind=dp)             :: rho_prime, wa
      !
      select case(nuclear_model)
      case("point")
         rho_prime = zero
      case("uniform")
         rho_prime = zero
      case("fermi")
         wa = exp( (r0-fermi_c_parameter)/fermi_a_parameter )
         rho_prime = - wa / ((one + wa)*(one + wa)*fermi_a_parameter)
      case default
         stop "nuclear_charge_density_prime(): program stop A."
      end select
      !
   end function nuclear_charge_density_prime
   !
   !
   function nuclear_extent()   result(Rnuc)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns a nuclear radius Rnuc which should include at least 99 percent of the total nuclear charge.
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp) :: Rnuc
      !
      select case(nuclear_model)
      case("point")
         Rnuc = zero
      case("uniform")
         Rnuc = 0.836_dp * (nuclear_mass)**(one/three) + 0.570_dp
         Rnuc = Rnuc   * convert_fermi_to_bohr
      case("fermi")
         Rnuc = 0.836_dp * (nuclear_mass)**(one/three) + 0.570_dp
         Rnuc = two * two * Rnuc   * convert_fermi_to_bohr
      case default
         stop "nuclear_extent(): program stop A."
      end select
      !
   end function nuclear_extent
   !
   !
   function nuclear_potential(r0)   result(nucpot)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the nuclear potential at point r0 in dependence of the nuclear model.
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      real(kind=dp), intent(in) :: r0
      real(kind=dp)             :: nucpot, Rnuc, a, abc, abc2, abc3, c, cba, tabc, thabc2, hpiac2, sabc3, dmsas, en, &
                                   rmc, rmcba, rbc, s2mcba, s3mcba, s2rcba, s3rcba, h3php, zz0, zbn
      !
      if (rabs_use_stop   .and.   r0 <= zero) then
         print *, "r0 = ",r0
         stop     "nuclear_potential(): program stop A."
      end if
      !
      select case(nuclear_model)
      case("point")
         !
         nucpot = - nuclear_charge / r0
      case("uniform")
         !
         Rnuc   = nuclear_extent()
         if( r0 <= Rnuc ) then
            nucpot =  - nuclear_charge/(two*Rnuc) * (three - r0*r0/(Rnuc*Rnuc))
         else
            nucpot = - nuclear_charge / r0
         end if            
      case("fermi")
         !
         Rnuc   = nuclear_extent()
         if( r0 <= two * two * Rnuc ) then
            c      = fermi_c_parameter
            a      = fermi_a_parameter
            abc    = a/c
            tabc   = two * abc
            abc2   = abc * abc
            thabc2 = three * abc2
            abc3   = abc2 * abc
            cba    = c / a
            hpiac2 = half * pi * pi * abc2
            h3php  = three / two + hpiac2
            call nuclear_potential_es(-cba,s2mcba,s3mcba)
            sabc3  = two * three * abc3
            dmsas  = -sabc3 * s3mcba
            en     = one + abc2 * pi * pi + dmsas
            zbn    = nuclear_charge / en
            !
            rmc    = r0 - c
            rmcba  = rmc / a
            rbc    = r0 / c
            if (rbc < one) then
               call nuclear_potential_es(rmcba,s2rcba,s3rcba)
               zz0 = zbn *(dmsas + sabc3 * s3rcba + rbc * (h3php - thabc2 * s2rcba - half * rbc * rbc))
            else
               call nuclear_potential_es(-rmcba,s2rcba,s3rcba)
               zz0 = nuclear_charge * (one + thabc2 * (rbc * s2rcba + tabc * s3rcba) / en)
            end if
            nucpot = - zz0 / r0
         else
            nucpot = - nuclear_charge / r0
         end if            
         !
      case("debye")
         !
         nucpot = - nuclear_charge / r0 * exp(-nuclear_screening_lambda * r0)
      case("debye-box")
         !
         if (r0  <= nuclear_debye_sphere_R0) then
            nucpot = - nuclear_charge / r0 * exp(-nuclear_screening_lambda * r0)
         else
            nucpot = zero
         end if
      case("ion-sphere")
         !
         if (r0  <= nuclear_ion_sphere_R0) then
            nucpot = - nuclear_charge / r0 + half * (nuclear_charge - nuclear_ion_sphere_Nb) / &
                       nuclear_ion_sphere_R0 * (three - r0*r0/(nuclear_ion_sphere_R0*nuclear_ion_sphere_R0) )
         else
            nucpot = zero
         end if
      case default
         stop "nuclear_potential(): program stop B."
      end select
      !
   end function nuclear_potential
   !
   !
   subroutine nuclear_potential_es(f,s2f,s3f) 
   !--------------------------------------------------------------------------------------------------------------------
   ! Evaluate the sum of the series     
   !
   !                       infinity      n              k 
   !              S  (f) =   Sum     (-1)  exp (n*f) / n  
   !               k        n = 0       
   !  
   ! for k = 2, 3 to machine precision. This is a utility subroutine for nuclear_potential(). This subroutine was 
   ! originally written by F A Parpia and adapted to the rabs package.
   !--------------------------------------------------------------------------------------------------------------------
      !
      implicit none
      real(kind=dp), intent(in)  :: f
      real(kind=dp), intent(out) :: s2f, s3f
      integer                    :: n
      real(kind=dp)              :: enf, fase, obn, s2last, term2, term3
      !
      n    = 0
      s2f  = zero
      s3f  = zero
      fase = one
      !
    1 n   = n+1
      obn   = one / n
      fase  = - fase
      enf   = exp (n*f)
      term2 = fase * enf * obn * obn
      term3 = term2* obn
      s2last= s2f
      s2f   = s2f + term2
      s3f   = s3f + term3
      if (abs(s2f) /= abs(s2last)) then
         goto 1
      end if
      !
   end subroutine nuclear_potential_es
   !
   !
   function nuclear_radial_points(r_grid)                  result(number)
   !--------------------------------------------------------------------------------------------------------------------
   ! Returns the number of radial points inside of the nuclear radius Rnuc as calculated in nuclear_extent().
   !--------------------------------------------------------------------------------------------------------------------
      !
      real(kind=dp), dimension(:), intent(in) :: r_grid
      !
      integer       :: i, number
      real(kind=dp) :: Rnuc
      !
      number = 0
      Rnuc   = nuclear_extent()
      !
      do  i = 1,1000
         if (r_grid(i) > Rnuc)  then
	    exit
	 end if
	 !
	 number = i
      end do
      !
   end function nuclear_radial_points
   !
end module rabs_nucleus
