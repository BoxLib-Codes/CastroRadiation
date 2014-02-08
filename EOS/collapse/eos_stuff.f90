module eos_module

  use bl_types
  use network
  use eos_type_module

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

  real(kind=dp_t), save :: table_Tmin, table_Tmax
  real(kind=dp_t), save :: table_Yemin, table_Yemax
  real(kind=dp_t), save :: table_rhomin, table_rhomax

  integer, parameter :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter :: eos_input_rs = 2   ! density, entropy are inputs
  integer, parameter :: eos_input_re = 5   ! density, internal energy are inputs
  integer, parameter :: eos_input_cv = 6
  integer, parameter :: eos_input_xnxp = 7
  integer, parameter :: eos_input_rp = 8
  integer, parameter :: eos_input_emin = 9

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  integer        , save, private :: numel

  logical, save, private :: initialized = .false.

! jbb added a constant to eint to have posistive internal energy
  real(kind=dp_t),parameter :: const_hack = 9.3d0  
  ! want to converge to the given energy (erg/g)
  ! energy (MeV) per baryon -> ergs/gm
  !   MeV   1.602176462e-6 erg   6.02214199e23
  !   --- x ------------------ x -------------
  !   bar          1 MeV            1 mole
  real(kind=dp_t),parameter :: enuc2ecgs = 0.95655684d18
  ! pressure (MeV/fm^3) -> dyn/cm^2
  real(kind=dp_t),parameter :: Pnuc2Pcgs = 1.60217733d33

  interface eos
     module procedure eos_old
     module procedure eos_new
  end interface eos

contains

!*************************************************************

  subroutine eos_init(small_temp, small_dens, gamma_in)

    use table_module, only : get_numel, collapse_init, t1, t2, r1, r2, y1, y2

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

    ! initialize collapse
    call collapse_init

    initialized = .true.

    call get_numel(numel)

    table_Tmin =  10.d0**t1
    table_Tmax =  10.d0**t2

    table_Yemin = y1
    table_Yemax = y2

    table_rhomin = 10.d0**r1
    table_rhomax = 10.d0**r2

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

  subroutine eos_given_RPX(R, P, T, e, Y, pt_index)

! input/output variables
    real(kind=dp_t)  , intent(  out) :: e
    real(kind=dp_t)  , intent(inout) :: T
    real(kind=dp_t)  , intent(in   ) :: R, P, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

    do_eos_diag = .false.

    temp_eos(1) = T
    den_eos(1)  = R
    p_eos(1)    = P

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_rp, den_eos, temp_eos, ye_eos, &
             npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
             cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

    e = e_eos(1)

  end subroutine eos_given_RPX

!*************************************************************

  subroutine eos_given_RSX(G, P, C, T, e, dpdr, dpde, R, S, Y, pt_index)

! input/output variables
    real(kind=dp_t), intent(  out) :: G, P, C, dpdr, dpde
    real(kind=dp_t), intent(inout) :: T, e
    real(kind=dp_t), intent(in   ) :: R, S, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    logical :: do_eos_diag

!   do_eos_diag = .false.
    do_eos_diag = .false.

    temp_eos(1) = T
    den_eos(1)  = R
    s_eos(1)    = S

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1)   = Y(nspec+1)

    call eos(eos_input_rs, den_eos, temp_eos, ye_eos, &
             npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
             cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

    G = gam1_eos(1)
    P = p_eos(1)
    C = cs_eos(1)
    T = temp_eos(1)
    e = e_eos(1)

! hack
    dpdr = 0.d0
    dpde = 0.d0
! end hack

  end subroutine eos_given_RSX

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

  subroutine eos_get_xnxp(xn, xp, R, T, Y, pt_index)

    real(kind=dp_t), intent(out) :: xn, xp
    real(kind=dp_t), intent(in)  :: R, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    logical :: do_eos_diag

    do_eos_diag = .false.

    den_eos(1) = R
    temp_eos(1) = T

    ! NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos(1) = Y(nspec+1)

    call eos(eos_input_xnxp, den_eos, temp_eos, ye_eos, &
         npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, &
         cv_eos, xn_eos, xp_eos, do_eos_diag, pt_index)

    xn = xn_eos(1)
    xp = xp_eos(1)

  end subroutine eos_get_xnxp


!*************************************************************

  !---------------------------------------------------------------------------
  ! new interface
  !---------------------------------------------------------------------------
  subroutine eos_new(input, eos_state, do_eos_diag, pt_index)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: eos_state
    logical,           intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

!    call eos_old(input, eos_state%rho, eos_state%T, &
    print *, "SHOULDN'T GET HERE"
    stop

  end subroutine eos_new

! a wrapper for interacting with Burrows eos
  subroutine eos_old(input, dens, temp, ye_in, &
                     npoints, pres, eint, &
                     gam1, cs, entropy, cv, xn, xp, &
                     do_eos_diag, pt_index)

    use bl_error_module
    use table_module, only : findthis, t1, t2, r1, r2, y1, y2

! dens     -- mass density (g/cc)
! temp     -- temperature (K)
! npoints  -- number of elements in i/o arrays
! ye       -- Ye
! pres     -- pressure (dyn/cm**2)
! enthalpy -- enthalpy (erg/g)
! eint     -- internal energy (erg/g)

    implicit none

!   Arguments
    logical          , intent(in   ) :: do_eos_diag
    integer          , intent(in   ) :: input, npoints
    integer, optional, intent(in   ) :: pt_index(:)

    double precision , intent(in   ) ::  dens(npoints), ye_in(npoints)
    double precision , intent(inout) ::  temp(npoints), eint(npoints)
    double precision , intent(  out) ::  pres(npoints)
    double precision , intent(  out) ::  gam1(npoints),   cs(npoints), entropy(npoints)

!   Local variables
    double precision :: error
    double precision :: energy_want(npoints)
    double precision :: entropy_want(npoints)
    double precision :: pressure_want(npoints)
    double precision :: tnew(npoints), dnew(npoints)
    double precision :: temp_row(npoints) 
    double precision :: cv(npoints)
    double precision :: ener1(npoints)
    double precision :: f(numel)
    double precision :: tempmev
    double precision :: xn(npoints)
    double precision :: xp(npoints)
    double precision :: ye(npoints)

    double precision, parameter :: ttol = 1.d-11
    double precision, parameter :: dtol = 1.d-11
    double precision, parameter :: stol = 1.d-11
    integer         , parameter :: max_newton = 100

    integer :: i, k, jy, jq, jr, icount
    integer :: nel
    integer :: iter, niter

    double precision :: min_temp, max_temp, temp_hold
    double precision :: temp0,depdt,dppdt,dt,dd
    double precision :: rho,to,tc,epo,epc,ppo,ppc

    if (.not. initialized) call bl_error('EOS: not initialized')

! table input:
!   - temp -- temperature (MeV)
!   - dens -- mass density (g/cc)
!   - ye   -- Ye
!
! table output:
!   - f(1)  -- energy (mev/baryon)
!   - f(2)  -- pressure (mev/fm^3)
!   - f(3)  -- entropy (per baryon per Boltzmann's constant)
!   - f(4)  -- specific heat at constant volume (ergs/MeV/gm)
!   - f(5)  -- mass fraction of neutrons
!   - f(6)  -- mass fraction of protons
!   - f(7)  -- mass fraction alpha particles
!   - f(8)  -- mass fraction heavy nuclei
!   - f(9)  -- not needed
!   - f(10) -- not needed
!   - f(11) -- not needed
!   - f(12) -- not needed
!   - f(13) -- not needed
!   - f(14) -- not needed
!   - f(15) -- not needed
!   - f(16) -- not needed [d pressure / d energy (??)]

    min_temp = 10.d0**t1
    max_temp = 10.d0**t2

    ! Check if density within bounds
    do k = 1, npoints
       if ( (dens(k) .lt. 10.d0**r1) .or. (dens(k) .gt. 10.d0**r2) ) then
          print *,'DENSITY OUT OF BOUNDS ',dens(k)
          print *,'LIMITS ARE            ',10.d0**r1,10.d0**r2
          print *,'FROM CALL WITH INPUT  ',input
          if (present(pt_index)) &
             print *,'AT POINT              ',pt_index(:)
          stop
       end if
    end do

    ! Enforce ye within bounds
    do k = 1, npoints
       ye(k) = ye_in(k)
       if ( (ye(k) .lt. y1) .or. (ye(k) .gt. y2) ) then
!          print *,'IN EOS_STUFF : input = ',input
!          print *,'YE   OUT OF BOUNDS   ',ye(k)
!          if (present(pt_index)) &
!               print *,'at point ',pt_index(:)
!          print *,'LIMITS ARE: Y1 Y2    ',y1,y2
!          print *,'SETTING TO LIMIT AND MOVING ON '

! Ye out of range can happen sometimes during problems run with radiation.
! This may be an indicator of insufficient resolution or missing physics,
! but for now suppress the per-point printing that was choking the output
! files and set Ye silently to the edge of the valid range:
          ye(k) = max(min(ye_in(k),y2),y1)
!         stop
       end if
    end do

    if (input .eq. eos_input_emin) then

       do k = 1, npoints
          call findthis(f,numel,min_temp+1.d-6*abs(min_temp),dens(k),ye(k),&
                        jy,jq,jr,input,pt_index)
          eint(k) = (f(1)+const_hack) * enuc2ecgs
       end do

    else if (input .eq. eos_input_rt) then

!---------------------------------------------------------------------------    
! input = 1: dens, temp, and ye are inputs                                   
!---------------------------------------------------------------------------    

       ! NOTE: we assume the temp is already in MeV
       ! tempmev = temp(k) / 1.160445d10

       do k = 1, npoints

          if (temp(k) .lt. min_temp .or. temp(k) .gt. max_temp) then
              print *,'TEMP OUT OF BOUNDS WITH RTX CALL ', temp(k)
              print *,'LIMITS ARE: T1 T2 ',min_temp, max_temp
              stop
          end if

          call findthis(f,numel,temp(k),dens(k),ye(k),jy,jq,jr,input,pt_index)

          ! energy (MeV) per baryon -> ergs/g
          eint(k) = (f(1)+const_hack) * enuc2ecgs

            cv(k) = f(4) ! * 1.160445d10

          ! pressure (MeV/fm^3) -> dyn/cm^2
          pres(k) = f(2) * Pnuc2Pcgs

          gam1(k) = f(12)
            cs(k) = sqrt(gam1(k)*pres(k)/dens(k))

          entropy(k) = f(3)

       enddo

    ! THIS IS BASED ON ADAM'S ITERATION SCHEME
    else if (input .eq. eos_input_re) then

!---------------------------------------------------------------------------
! input = 5: dens, energy, and ye are inputs
!---------------------------------------------------------------------------

       ! NOTE: we assume the temp is already in MeV
       do k = 1, npoints

          ! load initial guess
          temp_row(k) = temp(k)

          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)

          energy_want(k) = eint(k) / enuc2ecgs - const_hack

          if (do_eos_diag) print*,' '
          if (do_eos_diag) print*,'WANT e (erg/g) ',energy_want(k)

          icount = 0

          ! This is the initial guess
          temp0   = temp(k)

          call findthis(f,numel,max_temp,dens(k),ye(k),jy,jq,jr,input,pt_index)
          if (energy_want(k) .gt. f(1)) then
            print *,'eos_given_REX: energy too high :        ',energy_want(k)
            print *,'using rho temp ye ',dens(k), max_temp, ye(k)
            if (present(pt_index)) print *,'at point ',pt_index(:)
            print *,'max possible energy given this density  ',f(1)
            print *,'Setting temperature to max_temp ',max_temp
            temp_row(k) = max_temp
            go to 10
!           call bl_error('EOS: energy too high for table')
          end if

          call findthis(f,numel,min_temp+1.d-6*abs(min_temp),dens(k),ye(k),jy,jq,jr,input,pt_index)
          if (energy_want(k) .lt. f(1)) then
            print *,'eos_given_REX: energy too low :        ',energy_want(k)
            print *,'using rho temp ye ',dens(k), min_temp+1.d-6*abs(min_temp), ye(k)
            if (present(pt_index)) print *,'at point ',pt_index(:)
            print *,'min possible energy given this density  ',f(1)
            print *,'Setting temperature to min_temp ',min_temp*1.000001d0
            temp_row(k) = min_temp + 1.d-6*abs(min_temp)
            go to 10
!           call bl_error('EOS: energy too low for table')
          end if

          rho = dens(k)

          ! Make "old" value = 1.01 * incoming guess
          to = min(temp_row(k) * 1.01d0, max_temp)
          call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
          epo = f(1)

          ! Make "current" value = incoming guess
          tc  = temp_row(k)

          do i = 1,max_newton

             call findthis(f,numel,tc,rho,ye(k),jy,jq,jr,input,pt_index)
             epc = f(1)

             if (do_eos_diag) then
                print*,' '
                print*,'ITER ',i
                print*,'TEMP EP OLD ',to,epo
                print*,'TEMP EP NEW ',tc,epc
             end if

             ! Calculate the slope of energy vs temp
             depdt  = (epc-epo)/(tc-to)

             ! How far were we off to start with?
             dd = (energy_want(k)-epc)/depdt

             ! Add that much to the current guess.
             if (do_eos_diag) print *,'EPC - EPO ',epc, epo, epc-epo
             if (do_eos_diag) print *,'TC  -  TO ',tc,to,    tc-to
             if (do_eos_diag) print *,'DEPDT ',depdt

             if (do_eos_diag) print *,'ADDING DD      ',dd

             ! Reset "old" = "current".
             to = tc
             epo = epc

             ! Create new "current"
             tc = tc + dd
             temp_row(k) = tc

             ! Update the iteration counter
             icount = icount + 1

             ! If negative go back to the original guess and add 10%
!            if (temp_row(k).le.0.d0) then
             if (temp_row(k).le.min_temp) then
                 temp_row(k) = temp0 + 0.1d0*temp0
                 to = temp_row(k)
                 call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
                 epo = f(1)
                 tc = min(to * 1.01d0, max_temp)
             endif
!            if (temp_row(k).lt.min_temp) then
!                temp_row(k) = min_temp+1.e-6*abs(min_temp)
!                to = temp_row(k)
!                call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
!                epo = f(1)
!                if(epo.gt.energy_want(k))then
!                    print *,'energy at min_temp to high', min_temp
!                    print *, 'energy: ', epo, " want: ", energy_want(k)
!                    stop
!                endif
!                tc = min(to * 1.01d0, max_temp)
!            endif
             if (temp_row(k).gt.max_temp) then
                 temp_row(k) = max_temp*.99d0
                 to = temp_row(k)
!                temp_row(k) = temp0 - 0.1d0*temp0
!                to = temp_row(k)
                 call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
                 epo = f(1)
!                if(epo.lt.energy_want(k))then
!                    print *,'energy at max_temp to low', max_temp
!                    print *, 'energy: ', epo, " want: ", energy_want(k)
!                    stop
!                endif
                 tc = min(to * 0.99d0, max_temp)
             endif


             ! If the temperature isn't changing much we're done
             if (abs(dd/temp_row(k)).lt.ttol) goto 10

             if (temp_row(k) .lt. min_temp .or. temp_row(k) .gt. max_temp) then
                 print *,'TEMP OUT OF BOUNDS IN REX ITER ', temp(k),temp_row(k)
                 print *,'LIMITS ARE: T1 T2 ',min_temp, max_temp
                 stop
             end if

          enddo

 10       continue

          temp(k) = temp_row(k)

          ! If we iterated max_newton times and didn't solve:
          if (icount.eq.max_newton) then
             call bisection(energy_want(k),1,numel,temp0,rho,ye(k),input,pt_index)
             temp(k) = temp0
          end if

          ! pressure (MeV/fm^3) -> dyn/cm^2
          pres(k) = f(2) * Pnuc2Pcgs

          gam1(k) = f(12)

            cs(k) = sqrt(gam1(k)*pres(k)/rho)

          entropy(k) = f(3)

       end do ! k = 1,npoints

    ! THIS IS ALSO BASED ON ADAM'S ITERATION SCHEME
    else if (input .eq. eos_input_rp) then

!---------------------------------------------------------------------------
! input = 5: dens, pressure, and ye are inputs
!---------------------------------------------------------------------------

       ! NOTE: we assume the temp is already in MeV
       do k = 1, npoints

          ! load initial guess
          temp_row(k) = temp(k)

          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)

          pressure_want(k) = pres(k) / Pnuc2Pcgs

          if (do_eos_diag) print*,' '
          if (do_eos_diag) print*,'WANT SCALED PRESSURE ',pressure_want(k)

          icount = 0

          ! This is the initial guess
          temp0   = temp(k)

          rho = dens(k)

          ! Make "old" value = 1.01 * incoming guess
          to = min(temp_row(k) * 1.01d0, max_temp)
          call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
          ppo = f(2)

          ! Make "current" value = incoming guess
          tc  = temp_row(k)

          do i = 1,max_newton

             call findthis(f,numel,tc,rho,ye(k),jy,jq,jr,input,pt_index)
             ppc = f(2)

             if (do_eos_diag) then
                print*,' '
                print*,'ITER ',i
                print*,'TEMP PP OLD ',to,ppo
                print*,'TEMP PP NEW ',tc,ppc
             end if

             ! Calculate the slope of pressure vs temp
             dppdt  = (ppc-ppo)/(tc-to)

             ! How far were we off to start with?
             dd = (pressure_want(k)-ppc)/dppdt

             ! Add that much to the current guess.
             if (do_eos_diag) print *,'PPC - PPO ',ppc, ppo, ppc-ppo
             if (do_eos_diag) print *,'TC  -  TO ',tc,to,    tc-to
             if (do_eos_diag) print *,'DPPDT ',dppdt

             if (do_eos_diag) print *,'ADDING DD      ',dd

             ! Reset "old" = "current".
             to = tc
             ppo = ppc

             ! Create new "current"
             tc = tc + dd
             temp_row(k) = tc

             ! Update the iteration counter
             icount = icount + 1

             ! If negative go back to the original guess and add 10%
             if (temp_row(k).le.min_temp) then
                 temp_row(k) = temp0 + 0.1d0*temp0
                 to = temp_row(k)
                 call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
                 ppo = f(2)
                 tc = min(to * 1.01d0, max_temp)
             endif
             if (temp_row(k).gt.max_temp) then
                 temp_row(k) = max_temp*.99d0
                 to = temp_row(k)
                 call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
                 ppo = f(2)
                 tc = min(to * 0.99d0, max_temp)
             endif

             ! If the temperature isn't changing much we're done
             if (abs(dd/temp_row(k)).lt.ttol) goto 11

             if (temp_row(k) .lt. min_temp .or. temp_row(k) .gt. max_temp) then
                 print *,'TEMP OUT OF BOUNDS IN REX ITER ', temp(k),temp_row(k)
                 print *,'LIMITS ARE: T1 T2 ',min_temp, max_temp
                 stop
             end if

          enddo

 11       continue

          temp(k) = temp_row(k)

          ! If we iterated max_newton times and didn't solve:
          if (icount.eq.max_newton) then
             print *, 'no bisection routine for pressure solve in EOS'
             stop
          end if

          ! energy
          eint(k) = (f(1)+const_hack) * enuc2ecgs

       end do ! k = 1,npoints
 
    else if (input .eq. eos_input_rs) then

       ! NOTE: we assume the temp is already in MeV
       do k = 1, npoints

          ! load initial guess
          temp_row(k) = temp(k)

          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)

          entropy_want(k) = entropy(k)

          if (do_eos_diag) print*,' '
          if (do_eos_diag) print*,'WANT S  ',entropy_want(k)

          icount = 0

          ! This is the initial guess
          temp0   = temp(k)

          call findthis(f,numel,max_temp,dens(k),ye(k),jy,jq,jr,input,pt_index)
          if (entropy_want(k) .gt. f(3)) then
            print *,'eos_given_REX: entropy is too high given this rho and max_temp: ',entropy_want(k)
            print *,'maximum possible entropy given this density                     ',f(3)
            print *,'density and max_temp and Ye: ',dens(k), max_temp, ye(k)
            stop
          end if

          rho = dens(k)

          ! Make "old" value = 1.01 * incoming guess
          to = min(temp_row(k) * 1.01d0, max_temp)
          call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
          epo = f(3)

          ! Make "current" value = incoming guess
          tc  = temp_row(k)

          do i = 1,max_newton

             call findthis(f,numel,tc,rho,ye(k),jy,jq,jr,input,pt_index)
             epc = f(3)

             if (do_eos_diag) then
                print*,' '
                print*,'ITER ',i
                print*,'TEMP EP OLD ',to,epo
                print*,'TEMP EP NEW ',tc,epc
             end if

             ! Calculate the slope of energy vs temp
             depdt  = (epc-epo)/(tc-to)

             ! How far were we off to start with?
             dd = (entropy_want(k)-epc)/depdt

             ! Add that much to the current guess.
             if (do_eos_diag) print *,'EPC - EPO ',epc, epo, epc-epo
             if (do_eos_diag) print *,'TC  -  TO ',tc,to,    tc-to
             if (do_eos_diag) print *,'DEPDT ',depdt

             if (do_eos_diag) print *,'ADDING DD      ',dd

             ! Reset "old" = "current".
             to = tc
             epo = epc

             ! Create new "current"
             tc = tc + dd
             temp_row(k) = tc

             ! Update the iteration counter
             icount = icount + 1

             ! If negative go back to the original guess and add 10%
             if (temp_row(k).le.0.d0) then
                 temp_row(k) = temp0 + 0.1d0*temp0
                 to = temp_row(k)
                 call findthis(f,numel,to,rho,ye(k),jy,jq,jr,input,pt_index)
                 epo = f(3)
                 tc = min(to * 1.01d0, max_temp)
             endif

             ! If the temperature isn't changing much we're done
             if (abs(dd/temp_row(k)).lt.ttol) goto 210

             if (temp_row(k) .lt. min_temp .or. temp(k) .gt. max_temp) then
                 print *,'TEMP OUT OF BOUNDS IN REX ITER ', temp(k)
                 print *,'LIMITS ARE: T1 T2 ',min_temp, max_temp
                 stop
             end if

          enddo

210       continue

          temp(k) = temp_row(k)

          ! If we iterated max_newton times and didn't solve:
          if (icount.eq.max_newton) then
             call bisection(energy_want(k),1,numel,temp0,rho,ye(k),input,pt_index)
             call findthis(f,numel,temp0,rho,ye(k),jy,jq,jr,input,pt_index)
             temp(k) = temp0
          end if

          ! pressure (MeV/fm^3) -> dyn/cm^2
          pres(k) = f(2) * Pnuc2Pcgs

          gam1(k) = f(12)

            cs(k) = sqrt(gam1(k)*pres(k)/rho)

          eint(k) = (f(1)+const_hack) * enuc2ecgs

       end do ! k = 1,npoints

    else if (input .eq. eos_input_cv) then

       do k = 1, npoints

          call findthis(f,numel,temp(k),dens(k),ye(k),jy,jq,jr,input,pt_index)

          cv(k) = f(4)

       enddo

    else if (input .eq. eos_input_xnxp) then
          
          do k = 1, npoints
             call findthis(f,numel,temp(k),dens(k),ye(k),jy,jq,jr,input,pt_index)
             
             xn(k) = f(5)
             xp(k) = f(6)
          enddo

    else

       call bl_error('EOS: conversion not implemented')

    endif

  end subroutine eos_old

  subroutine bisection(target,nn,numel,temp,rho,ye,input,pt_index)

      use bl_error_module
      use table_module, only : findthis, t1, t2

      implicit none

      integer          , intent(in   ) :: nn, numel, input
      integer, optional, intent(in   ) :: pt_index(:)
      double precision , intent(in   ) :: target,rho,ye
      double precision , intent(inout) :: temp

      double precision :: ep,e,ep1,ep2,dt,depdt,dd
      double precision :: temp1,temp2,dtemp,tmid,temp00
      double precision :: f(numel),f1(numel),f2(numel),fmid(numel)
      double precision :: min_temp,max_temp,eps_ener
      integer          :: i,jy,jq,jr,mcount
      integer          :: input_bisect
      double precision, parameter :: eps_temp = 1.d-11
      integer, parameter :: max_bisec = 1000

      min_temp = 10.d0**t1
      max_temp = 10.d0**t2

      input_bisect = 100 + input

      ! Compute the energy of max_temp to be used in defining whether we're "close enough"
      call findthis(f1,numel,max_temp,rho,ye,jy,jq,jr,input_bisect+10000,pt_index)
      eps_ener = 1.d-12 * f1(nn)

      ! target is the target value for energy
      ! nn is the index in f, i.e. energy = f(1)
      ! temp is what we're solving for
      ! rho and ye are fixed values of density and Ye

      temp00 = temp
      mcount = 0

!     temp1 = 0.9d0 * temp
!     temp2 = 1.1d0 * temp
      temp1 = max(0.8d0 * temp,min_temp)
      temp2 = min(1.2d0 * temp,max_temp)

      call findthis(f,numel,temp,rho,ye,jy,jq,jr,input_bisect,pt_index)
      f(nn)=f(nn)-target

      call findthis(f1,numel,temp1,rho,ye,jy,jq,jr,input_bisect+1000,pt_index)
      f1(nn)=f1(nn)-target

      call findthis(f2,numel,temp2,rho,ye,jy,jq,jr,input_bisect+2000,pt_index)
      f2(nn)=f2(nn)-target

 2    continue

      if (f1(nn)*f2(nn).ge.0.d0) then

         mcount=mcount+1

         temp1 = max(0.8d0 * temp1,min_temp)
         temp2 = min(1.2d0 * temp2,max_temp)

         call findthis(f1,numel,temp1,rho,ye,jy,jq,jr,input_bisect+3000,pt_index)
         call findthis(f2,numel,temp2,rho,ye,jy,jq,jr,input_bisect+4000,pt_index)
 
         f1(nn)=f1(nn)-target
         f2(nn)=f2(nn)-target

         if (abs(f1(nn)) .le. eps_ener) then
            temp = temp1
            goto 10
         end if

         if (abs(f2(nn)) .le. eps_ener) then
            temp = temp2 
            goto 10
         end if

         if (mcount.le.max_bisec) then
            goto 2
         else
            write(6,*) 'BISECTION FAILED in eos_stuff'
            if (present(pt_index)) &
                  print *,'at point ',pt_index(:)
            write(6,*) 'MINTEMP / MAXTEMP ',min_temp,max_temp
            write(6,*) '  TEMP1 / TEMP2   ',temp1, temp2
            write(6,*) 'TARGET ENERGY     ',target
            call findthis(f1,numel,min_temp,rho,ye,jy,jq,jr,input_bisect+10000,pt_index)
            write(6,*) 'ENERGY OF MINTEMP ',f1(nn)
            call findthis(f1,numel,max_temp,rho,ye,jy,jq,jr,input_bisect+10000,pt_index)
            write(6,*) 'ENERGY OF MAXTEMP ',f1(nn)
            call bl_error('eos_given_ReX: bisection cant bracket energy')
         endif
      endif

      if (f1(nn).lt.0.d0)then
         temp=temp1
         dtemp=temp2-temp1
      else
         temp=temp2
         dtemp=temp1-temp2
      endif

      do i=1,max_bisec
         dtemp = dtemp*0.5d0
          tmid = temp+dtemp
         call findthis(fmid,numel,tmid,rho,ye,jy,jq,jr,input_bisect+5000,pt_index)
         fmid(nn) = fmid(nn) - target
         if (fmid(nn).le.0.d0) temp = tmid
         if (abs(dtemp).lt.eps_temp) goto 10
      enddo

      call bl_error('eos_given_ReX: bisection cant get to the energy we want after 200 iterations')

 10   continue

  end subroutine bisection

end module eos_module
