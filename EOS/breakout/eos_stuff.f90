module actual_eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use eos_type_module
  use eos_data_module

  implicit none

  character (len=64) :: eos_name = "break"  
  
  double precision, save :: gamma_const

contains

  subroutine actual_eos_init

    use extern_probin_module, only: eos_gamma

    implicit none
 
    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       gamma_const = FIVE3RD
    end if

  end subroutine actual_eos_init



  subroutine actual_eos(input, state)

    use fundamental_constants_module, only: k_B, n_A

    implicit none

    integer,             intent(in   ) :: input
    type (eos_t_vector), intent(inout) :: state

    double precision, parameter :: R = k_B*n_A

    integer :: j
    double precision :: poverrho

    ! Calculate mu.

    ! xxxxxxx hack!
    ! The only difference between this EOS and gamma-law
    do j = 1, state % N
       state % mu = one / state % aux(j,2)
    end do

    select case (input)

    case (eos_input_rt)

       ! dens, temp and xmass are inputs
       do j = 1, state % N
          state % cv(j) = R / (state % mu(j) * (gamma_const-ONE)) 
          state % e(j) = state % cv(j) * state % T(j)
          state % p(j) = (gamma_const-ONE) * state % rho(j) * state % e(j)
          state % gam1(j) = gamma_const
       end do

    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       call bl_error('EOS: eos_input_rh is not supported in this EOS.')

    case (eos_input_tp)

       ! temp, pres, and xmass are inputs

       call bl_error('EOS: eos_input_tp is not supported in this EOS.')
       
    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       do j = 1, state % N
          poverrho = state % p(j) / state % rho(j)
          state % T(j) = poverrho * state % mu(j) * (ONE/R)
          state % e(j) = poverrho * (ONE/(gamma_const-ONE))
          state % gam1(j) = gamma_const
       end do

    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       do j = 1, state % N
          poverrho = (gamma_const - ONE) * state % e(j)

          state % p(j) = poverrho * state % rho(j)
          state % T(j) = poverrho * state % mu(j) * (ONE/R)
          state % gam1(j) = gamma_const
          
          ! sound speed
          state % cs(j) = sqrt(gamma_const * poverrho)

          state % dpdr_e(j) = poverrho
          state % dpde(j) = (gamma_const-ONE) * state % rho(j)

          ! Try to avoid the expensive log function.  Since we don't need entropy 
          ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
          state % s(j) = ONE  
       end do

    case (eos_input_ps)

       ! pressure entropy, and xmass are inputs

       call bl_error('EOS: eos_input_ps is not supported in this EOS.')
       
    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       call bl_error('EOS: eos_input_ph is not supported in this EOS.')
       
    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       ! This system is underconstrained.
       
       call bl_error('EOS: eos_input_th is not a valid input for the gamma law EOS.')

    case default
       
       call bl_error('EOS: invalid input.')
       
    end select
    
  end subroutine actual_eos

end module actual_eos_module
