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
  real(kind=dp_t) :: xn_eos(NP)
  real(kind=dp_t) :: xp_eos(NP)
  real(kind=dp_t) :: xne_eos(NP)
  real(kind=dp_t) :: eta_eos(NP)
  real(kind=dp_t) :: pele_eos(NP)
  real(kind=dp_t) :: dpdt_eos(NP)
  real(kind=dp_t) :: dpdr_eos(NP)
  real(kind=dp_t) :: dedr_eos(NP)
  real(kind=dp_t) :: dedt_eos(NP)
  real(kind=dp_t) :: gam1_eos(NP)
  real(kind=dp_t) ::   cs_eos(NP)
  real(kind=dp_t) ::    s_eos(NP)
  real(kind=dp_t) :: dsdt_eos(NP)
  real(kind=dp_t) :: dsdr_eos(NP)
  real(kind=dp_t) :: dpdX_eos(NP,nspec)
  real(kind=dp_t) :: dhdX_eos(NP,nspec)

  common /eos_common/ ye_eos,temp_eos,den_eos,e_eos,p_eos,h_eos,cv_eos,cp_eos
  common /eos_common/ xn_eos,xp_eos,xne_eos,eta_eos,pele_eos,dpdt_eos,dpdr_eos,dedr_eos
  common /eos_common/ dedt_eos,gam1_eos,cs_eos,s_eos,dsdt_eos,dsdr_eos,dpdX_eos,dhdX_eos
  SAVE /eos_common/
!$omp threadprivate(/eos_common/)

  integer, parameter :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter :: eos_input_re = 5   ! density, internal energy are inputs
  integer, parameter :: eos_input_cv = 6
  integer, parameter :: eos_input_emin = 9

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld

  logical, save, private :: initialized = .false.

! jbb added a constant to eint to have posistive internal energy
  real(kind=dp_t),parameter :: const_hack = 9.3d0*0.95655684d18

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

    call readtable("HShenEOS_rho220_temp180_ye65_version2.0_20111026_EOSmaker_svn9.h5")

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

  subroutine min_e_given_RX(R, Y, e, pt_index)

    ! input/output variables
    real(kind=dp_t)  , intent(  out) :: e
    real(kind=dp_t)  , intent(in   ) :: R, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

    do_eos_diag = .false.

    den_eos(1)  = R

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_emin, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

    e = e_eos(1)

  end subroutine min_e_given_RX

!*************************************************************

  subroutine eos_given_ReX(G, P, C, T, dpdr, dpde, R, e, Y, pt_index)

! input/output variables
    real(kind=dp_t)  , intent(  out) :: G, P, C, dpdr, dpde
    real(kind=dp_t)  , intent(inout) :: T
    real(kind=dp_t)  , intent(in   ) :: R, e, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

    do_eos_diag = .false.

    temp_eos(1) = T
    den_eos(1)  = R
    e_eos(1)    = e

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_re, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

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

!   input/output variables
    real(kind=dp_t), intent(in   ) :: R, e, Y(:)
    real(kind=dp_t), intent(inout) :: T
    real(kind=dp_t), intent(  out) :: S
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

!   do_eos_diag = .false.
    do_eos_diag = .false.

    temp_eos(1) = T
    den_eos(1)  = R
    e_eos(1)    = e

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_re, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

    S = s_eos(1)
    T = temp_eos(1)

  end subroutine eos_S_given_ReX

!*************************************************************

  subroutine eos_given_RTX(e, P, R, T, Y, pt_index)

! input/output variables
    real(kind=dp_t), intent(out) :: e, P
    real(kind=dp_t), intent(in)  :: R, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    logical :: do_eos_diag

!    do_eos_diag = .false.
    do_eos_diag = .false.

    den_eos(1)  = R
    temp_eos(1) = T

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_rt, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

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

    call eos(eos_input_cv, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

    cv = cv_eos(1)

  end subroutine eos_get_cv

!*************************************************************
!*************************************************************

! a wrapper for interacting with Ott's Shen EOS
  subroutine eos(input, dens, temp, ye_in, &
                 npoints, pres, eint, &
                 gam1, cs, entropy, cv, xn, xp, &
                 do_eos_diag, pt_index)

    use bl_error_module
    use eosmodule 

    implicit none

!   Arguments
    logical          , intent(in   ) :: do_eos_diag
    integer          , intent(in   ) :: input, npoints
    integer, optional, intent(in   ) :: pt_index(:)

    double precision , intent(in   ) ::  dens(npoints), ye_in(npoints)
    double precision , intent(inout) ::  temp(npoints), eint(npoints)
    double precision , intent(  out) ::  pres(npoints)
    double precision , intent(  out) ::  gam1(npoints),   cs(npoints), entropy(npoints)
    double precision , intent(  out) ::  cv(npoints), xn(npoints), xp(npoints)

!   Local variables
    integer :: k, keytemp, keyerr
    double precision :: xdpderho,xdpdrhoe,xmunu,eint_loc
    double precision, dimension(npoints) :: Ye_loc
    
    if (.not. initialized) call bl_error('EOS: not initialized')

    ! Check if density within bounds
    do k = 1, npoints
       if ( (dens(k) .lt. eos_rhomin) .or. (dens(k) .gt. eos_rhomax) ) then
          print *,'DENSITY OUT OF BOUNDS ',dens(k)
          print *,'LIMITS ARE            ',eos_rhomin,eos_rhomax
          print *,'FROM CALL WITH INPUT  ',input
          if (present(pt_index)) &
             print *,'AT POINT              ',pt_index(:)
          stop
       end if
    end do
    
    ! Enforce ye within bounds
    do k = 1, npoints
       Ye_loc(k) = ye_in(k)
       if ( (ye(k) .lt. eos_yemin) .or. (ye(k) .gt. eos_yemax) ) then
!          print *,'IN EOS_STUFF : input = ',input
!          print *,'YE   OUT OF BOUNDS   ',ye(k)
!          if (present(pt_index)) &
!               print *,'at point ',pt_index(:)
!          print *,'LIMITS ARE: Y1 Y2    ',eos_yemin,eos_yemax
!          print *,'SETTING TO LIMIT AND MOVING ON '
          ye(k) = max(min(ye_in(k),eos_yemax),eos_yemin)
!         stop
       end if
    end do

    if (input .eq. eos_input_emin) then

       keytemp = 1

       do k = 1, npoints

          call nuc_eos_short(dens(k),eos_tempmin*(1.d0+1.d-6),Ye_loc(k),eint(k),pres(k),entropy(k),&
               cs(k),cv(k),xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,precision)

          eint(k) = eint(k) + const_hack

       end do

    else if (input .eq. eos_input_rt) then

!---------------------------------------------------------------------------    
! dens, temp, and ye are inputs                                   
!---------------------------------------------------------------------------    

       keytemp = 1

       do k = 1, npoints

          if (temp(k) .lt. eos_tempmin .or. temp(k) .gt. eos_tempmax) then
             print *,'TEMP OUT OF BOUNDS WITH RTX CALL ', temp(k)
             print *,'LIMITS ARE: T1 T2 ',eos_tempmin,eos_tempmax
             stop
          end if

          call nuc_eos_short(dens(k),temp(k),Ye_loc(k),eint(k),pres(k),entropy(k),&
               cs(k),cv(k),xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,precision)

          gam1(k) = cs(k) * dens(k) / pres(k)
          cs(k) = sqrt(cs(k))
          eint(k) = eint(k) + const_hack
          
       enddo
 
    else if (input .eq. eos_input_re) then

!---------------------------------------------------------------------------
! input = 5: dens, energy, and ye are inputs
!---------------------------------------------------------------------------

       keytemp = 0
       
       do k = 1, npoints

          eint_loc = eint(k) - const_hack
          
          call nuc_eos_short(dens(k),temp(k),Ye_loc(k),eint_loc,pres(k),entropy(k),&
               cs(k),cv(k),xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,precision)
          
          gam1(k) = cs(k) * dens(k) / pres(k)
          cs(k) = sqrt(cs(k))
          
       enddo
       
    else

       call bl_error('EOS: conversion not implemented')

    end if

  end subroutine eos
!
end module eos_module
