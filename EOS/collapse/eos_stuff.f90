module eos_module

  use bl_types
  use bl_error_module
  use network
  use eos_type_module
  use eos_data_module

  implicit none

  real(kind=dp_t), save :: table_Tmin, table_Tmax
  real(kind=dp_t), save :: table_Yemin, table_Yemax
  real(kind=dp_t), save :: table_rhomin, table_rhomax

  integer, save, private :: numel

! jbb added a constant to eint to have posistive internal energy
  real(kind=dp_t),parameter, private :: const_hack = 9.3d0  
  ! want to converge to the given energy (erg/g)
  ! energy (MeV) per baryon -> ergs/gm
  !   MeV   1.602176462e-6 erg   6.02214199e23
  !   --- x ------------------ x -------------
  !   bar          1 MeV            1 mole
  real(kind=dp_t),parameter, private :: enuc2ecgs = 0.95655684d18
  real(kind=dp_t),parameter, private :: ecgs2enuc = 1.d0/enuc2ecgs
  ! pressure (MeV/fm^3) -> dyn/cm^2
  real(kind=dp_t),parameter, private :: Pnuc2Pcgs = 1.60217733d33

contains

  ! EOS initialization routine
  subroutine eos_init(small_temp, small_dens)

    use table_module, only : get_numel, collapse_init, t1, t2, r1, r2, y1, y2

    implicit none

    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens

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


  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag_in, pt_index)

    use table_module, only : findthis

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state
    logical, optional, intent(in   ) :: do_eos_diag_in
    integer, optional, intent(in   ) :: pt_index(:)

    double precision, parameter :: ttol = 1.d-11
    double precision, parameter :: dtol = 1.d-11
    double precision, parameter :: stol = 1.d-11
    integer         , parameter :: max_newton = 100

    logical :: do_eos_diag
    integer :: jy, jq, jr, icount, i
    double precision :: rho, eint, ye, temp, energy_want, Tmin, to, epo, tc, epc
    double precision :: depdt, dd
    double precision :: f(numel)

    if (.not. initialized) call bl_error('EOS: not initialized')

    if (present(do_eos_diag_in)) then
       do_eos_diag = do_eos_diag_in
    else
       do_eos_diag = .false.
    end if

    ! Check if density within bounds
    if ( (state % rho .lt. table_rhomin) .or. (state % rho .gt. table_rhomax) ) then
       print *,'DENSITY OUT OF BOUNDS ',state % rho
       print *,'LIMITS ARE            ',table_rhomin, table_rhomax
       print *,'FROM CALL WITH INPUT  ',input
       if (present(pt_index)) &
            print *,'AT POINT              ',pt_index(:)
       stop
    end if

    ! Enforce ye within bounds
    ye = min(max(state % aux(1), table_Yemin), table_Yemax)

    select case (input)

    case (eos_input_rt)

       ! dens, temp, and ye are inputs
       ! NOTE: we assume the temp is already in MeV
       
       if (state % T .lt. table_Tmin .or. state % T .gt. table_Tmax) then
          print *,'TEMP OUT OF BOUNDS WITH RTX CALL ', state % T
          print *,'LIMITS ARE: T1 T2 ',table_Tmin, table_Tmax
          stop
       end if

       call findthis(f,numel,state%T,state%rho,ye,jy,jq,jr,input,pt_index)

       ! energy (MeV) per baryon -> ergs/g
       state % e = (f(1)+const_hack) * enuc2ecgs

       ! pressure (MeV/fm^3) -> dyn/cm^2
       state % p = f(2) * Pnuc2Pcgs

       state % s = f(3)
       
       state % cv = f(4)

       state % gam1 = f(12)

       state % cs = sqrt(state % gam1 * state % p / state % rho)

    case (eos_input_re)

       ! dens, energy, and ye are inputs
       ! T is initial guess
       rho  = state % rho
       eint = state % e
       temp = state % T

       if (do_eos_diag) print*,'T/D INIT ',temp, rho

       energy_want = eint * ecgs2enuc - const_hack

       if (do_eos_diag) print*,' '
       if (do_eos_diag) print*,'WANT e (erg/g) ',energy_want

       call findthis(f,numel,table_Tmax,rho,ye,jy,jq,jr,input,pt_index)
       if (energy_want .gt. f(1)) then
          print *,'eos: energy too high :        ',energy_want
          print *,'using rho temp ye ',rho,table_Tmax,ye
          if (present(pt_index)) print *,'at point ',pt_index(:)
          print *,'max possible energy given this density  ',f(1)
          print *,'Setting temperature to max_temp ',table_Tmax
          temp = table_Tmax
          go to 10
       end if

       Tmin = table_Tmin+1.d-6*abs(table_Tmin)

       call findthis(f,numel,Tmin,rho,ye,jy,jq,jr,input,pt_index)
       if (energy_want .lt. f(1)) then
          print *,'eos: energy too low :        ',energy_want
          print *,'using rho temp ye ',rho, Tmin, ye
          if (present(pt_index)) print *,'at point ',pt_index(:)
          print *,'min possible energy given this density  ',f(1)
          print *,'Setting temperature to min_temp ',Tmin
          temp = Tmin
          go to 10
       end if

       ! Make "old" value = 1.01 * incoming guess
       to = min(temp * 1.01d0, table_Tmax)
       call findthis(f,numel,to,rho,ye,jy,jq,jr,input,pt_index)
       epo = f(1)

       ! Make "current" value = incoming guess
       tc = temp

       icount = 0

       do i = 1,max_newton

          ! Update the iteration counter
          icount = icount + 1

          call findthis(f,numel,tc,rho,ye,jy,jq,jr,input,pt_index)
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
          dd = (energy_want-epc)/depdt

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
          temp = tc

          ! If negative go back to the original guess and add 10%
          if (temp.le.table_Tmin) then
             temp = state%T*1.1d0
             to = temp
             call findthis(f,numel,to,rho,ye,jy,jq,jr,input,pt_index)
             epo = f(1)
             tc = min(to * 1.01d0, table_Tmin)
          endif

          if (temp.gt.table_Tmax) then
             temp = table_Tmax*.99d0
             to = temp
             call findthis(f,numel,to,rho,ye,jy,jq,jr,input,pt_index)
             epo = f(1)
             tc = min(to * 0.99d0, table_Tmax)
          endif

          ! If the temperature isn't changing much we're done
          if (abs(dd/temp).lt.ttol) goto 10

       end do

10     continue

       ! If we iterated max_newton times and didn't solve:
       if (icount.eq.max_newton) then
          temp = state % T
          call bisection(energy_want,1,numel,temp,rho,ye,input,pt_index)
          call findthis(f,numel,temp,rho,ye,jy,jq,jr,input,pt_index)
       end if

       state % T = temp
       
       ! pressure (MeV/fm^3) -> dyn/cm^2
       state % p = f(2) * Pnuc2Pcgs

       state % s = f(3)
       
       state % cv = f(4)

       state % gam1 = f(12)

       state % cs = sqrt(state % gam1 * state % p / state % rho)

    case default

       call bl_error('EOS: invalid input.')

    end select

    state % dpdr   = 0.d0
    state % dpde   = 0.d0
    state % dpdr_e = 0.d0

  end subroutine eos


  subroutine bisection(tget,nn,numel,temp,rho,ye,input,pt_index)

      use bl_error_module
      use table_module, only : findthis

      implicit none

      integer          , intent(in   ) :: nn, numel, input
      integer, optional, intent(in   ) :: pt_index(:)
      double precision , intent(in   ) :: tget,rho,ye
      double precision , intent(inout) :: temp

      double precision :: ep,e,ep1,ep2,dt,depdt,dd
      double precision :: temp1,temp2,dtemp,tmid,temp00
      double precision :: f(numel),f1(numel),f2(numel),fmid(numel)
      double precision :: min_temp,max_temp,eps_ener
      integer          :: i,jy,jq,jr,mcount
      integer          :: input_bisect
      double precision, parameter :: eps_temp = 1.d-11
      integer, parameter :: max_bisec = 1000

      min_temp = table_Tmin
      max_temp = table_Tmax

      input_bisect = 100 + input

      ! Compute the energy of max_temp to be used in defining whether we're "close enough"
      call findthis(f1,numel,max_temp,rho,ye,jy,jq,jr,input_bisect+10000,pt_index)
      eps_ener = 1.d-12 * f1(nn)

      ! tget is the target value for energy
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
      f(nn)=f(nn)-tget

      call findthis(f1,numel,temp1,rho,ye,jy,jq,jr,input_bisect+1000,pt_index)
      f1(nn)=f1(nn)-tget

      call findthis(f2,numel,temp2,rho,ye,jy,jq,jr,input_bisect+2000,pt_index)
      f2(nn)=f2(nn)-tget

 2    continue

      if (f1(nn)*f2(nn).ge.0.d0) then

         mcount=mcount+1

         temp1 = max(0.8d0 * temp1,min_temp)
         temp2 = min(1.2d0 * temp2,max_temp)

         call findthis(f1,numel,temp1,rho,ye,jy,jq,jr,input_bisect+3000,pt_index)
         call findthis(f2,numel,temp2,rho,ye,jy,jq,jr,input_bisect+4000,pt_index)
 
         f1(nn)=f1(nn)-tget
         f2(nn)=f2(nn)-tget

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
            write(6,*) 'TARGET ENERGY     ',tget
            call findthis(f1,numel,min_temp,rho,ye,jy,jq,jr,input_bisect+10000,pt_index)
            write(6,*) 'ENERGY OF MINTEMP ',f1(nn)
            call findthis(f1,numel,max_temp,rho,ye,jy,jq,jr,input_bisect+10000,pt_index)
            write(6,*) 'ENERGY OF MAXTEMP ',f1(nn)
            call bl_error('eos: bisection cant bracket energy')
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
         fmid(nn) = fmid(nn) - tget
         if (fmid(nn).le.0.d0) temp = tmid
         if (abs(dtemp).lt.eps_temp) goto 10
      enddo

      call bl_error('eos: bisection cant get to the energy we want after many iterations')

 10   continue

  end subroutine bisection

end module eos_module
