module eos_module

  use bl_types
  use network

  implicit none

  integer, parameter :: NP = 1
  integer, parameter :: npts = 1

!--------------------------------------------

  real(kind=dp_t) :: ye_eos(NP)
  real(kind=dp_t) :: temp_eos(NP)
  real(kind=dp_t) :: den_eos(NP)
  real(kind=dp_t) :: e_eos(NP)
  real(kind=dp_t) :: p_eos(NP)
  real(kind=dp_t) :: h_eos(NP)
  real(kind=dp_t) :: cv_eos(NP)
  real(kind=dp_t) :: cp_eos(NP)
  real(kind=dp_t) :: xne_eos(NP)
  real(kind=dp_t) :: eta_eos(NP)
  real(kind=dp_t) :: pele_eos(NP)
  real(kind=dp_t) :: dpdt_eos(NP)
  real(kind=dp_t) :: dpdr_eos(NP)
  real(kind=dp_t) :: dedr_eos(NP)
  real(kind=dp_t) :: dedt_eos(NP)
  real(kind=dp_t) :: dedy_eos(NP)
  real(kind=dp_t) :: gam1_eos(NP)
  real(kind=dp_t) ::   cs_eos(NP)
  real(kind=dp_t) ::    s_eos(NP)
  real(kind=dp_t) :: dsdt_eos(NP)
  real(kind=dp_t) :: dsdr_eos(NP)
  real(kind=dp_t) :: dpdX_eos(NP,nspec)
  real(kind=dp_t) :: dhdX_eos(NP,nspec)

  integer, parameter :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter :: eos_input_re = 5   ! density, internal energy are inputs

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  integer        , save, private :: numel

  logical, save, private :: initialized = .false.

contains

!*************************************************************

  subroutine eos_init(small_temp, small_dens, gamma_in)

    implicit none

    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens
    real(kind=dp_t), intent(in), optional :: gamma_in

    if (present(small_temp)) then
       if (small_temp > 0.d0) then
          smallt = small_temp
       else
          ! smallt = 5.d6
          ! Since we're now in MeV, divide by 10^10
          smallt = 5.d-4
       end if
    else
       ! smallt = 5.d6
       ! Since we're now in MeV, divide by 10^10
       smallt = 5.d-4
    endif

    if (present(small_dens)) then
       if (small_dens > 0.d0) then
          smalld = small_dens
       else
          smalld = 1.d-5
       end if
    else
       smalld = 1.d-5
    end if

    call LOADMX()

    initialized = .true.

  end subroutine eos_init

!*************************************************************

  subroutine eos_get_small_temp(small_temp_out)

    real(kind=dp_t), intent(out) :: small_temp_out

    small_temp_out = smallt

  end subroutine eos_get_small_temp

!*************************************************************

  subroutine eos_get_small_dens(small_dens_out)

    real(kind=dp_t), intent(out) :: small_dens_out

    small_dens_out = smalld

  end subroutine eos_get_small_dens

!*************************************************************

  subroutine eos_given_ReX(G, P, C, T, dpdr, dpde, R, e, Y, pt_index)

! input/output variables
    real(kind=dp_t), intent(  out) :: G, P, C, dpdr, dpde
    real(kind=dp_t), intent(inout) :: T
    real(kind=dp_t), intent(in   ) :: R, e, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

    do_eos_diag = .false.
!   do_eos_diag = .true.

    temp_eos(1) = T
    den_eos(1)  = R
    e_eos(1)    = e

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_re, den_eos, temp_eos, ye_eos, &
             npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
             dedt_eos, dedy_eos, do_eos_diag)

    G = gam1_eos(1)
    P = p_eos(1)
    C = cs_eos(1)
    T = temp_eos(1)

! hack
    dpdr = 0.d0
    dpde = 0.d0
! end hack

  end subroutine eos_given_ReX

!*************************************************************

  subroutine eos_S_given_ReX(S, R, e, T, Y, pt_index)

! input/output variables
    real(kind=dp_t), intent(out) :: S
    real(kind=dp_t), intent(in)  :: R, e, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

    do_eos_diag = .false.
!   do_eos_diag = .true.

    temp_eos(1) = T
    den_eos(1)  = R
    e_eos(1)    = e

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_re, den_eos, temp_eos, ye_eos, &
             npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
             dedt_eos, dedy_eos, do_eos_diag)

    S = s_eos(1)

  end subroutine eos_S_given_ReX

!*************************************************************

  subroutine eos_given_RTX(e, P, R, T, Y, pt_index)

! input/output variables
    real(kind=dp_t), intent(out) :: e, P
    real(kind=dp_t), intent(in)  :: R, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    logical :: do_eos_diag

    do_eos_diag = .false.

    den_eos(1)  = R
    temp_eos(1) = T

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_rt, den_eos, temp_eos, ye_eos, &
             npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
             dedt_eos, dedy_eos, do_eos_diag)

    p = p_eos(1)
    e = e_eos(1)

  end subroutine eos_given_RTX


!*************************************************************

  subroutine eos_get_cv(cv, R, T, Y, pt_index)

! input/output variables
    real(kind=dp_t), intent(out) :: cv
    real(kind=dp_t), intent(in)  :: R, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    logical :: do_eos_diag

    do_eos_diag = .false.

    den_eos(1)  = R
    temp_eos(1) = T

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1) = Y(nspec+1)

    call eos(eos_input_rt, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         dedt_eos, dedy_eos, do_eos_diag)
    
    cv = dedt_eos(1)

  end subroutine eos_get_cv

!*************************************************************

  subroutine eos_get_dedx(dedT, dedY, R, T, Y, dT, dY, pt_index)

! input/output variables
    real(kind=dp_t), intent(out) :: dedT, dedY
    real(kind=dp_t), intent(in)  :: R, T, Y(:)
    real(kind=dp_t), optional, intent(in)  :: dT, dY
    integer, optional, intent(in   ) :: pt_index(:)

    logical :: do_eos_diag

    do_eos_diag = .false.

    den_eos(1)  = R
    temp_eos(1) = T

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1) = Y(nspec+1)

    call eos(eos_input_rt, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         dedt_eos, dedy_eos, do_eos_diag)
    
    dedT = dedt_eos(1)
    dedY = dedy_eos(1)

  end subroutine eos_get_dedx

end module eos_module
