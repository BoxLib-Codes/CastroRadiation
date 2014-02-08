module opacity_table_module

  integer         , save :: nromax, ntmax, nyemax, ngmax
  double precision, save :: romin, romax, tmin, tmaxi
  double precision, save :: yemin, yemax
  double precision, save ::    enuk(50)
  real            , save :: absopac(50,50,50,40,3)
  real            , save ::    emis(50,50,50,40,3)
  real            , save :: scaopac(50,50,50,40,3)
  real            , save ::  sdelta(50,50,50,40,3)

contains

  ! dummy subroutine
  subroutine get_opacities(kp, kr, rho, temp, rhoYe, nu, get_Planck_mean, get_Rosseland_mean)
    implicit none
    logical, intent(in) :: get_Planck_mean, get_Rosseland_mean
    double precision, intent(in) :: rho, temp, rhoYe, nu
    double precision, intent(out) :: kp, kr
    kp = 0.d0
    kr = 0.d0
  end subroutine get_opacities


  subroutine prep_opacity(g, inu, er, der)

    use rad_params_module, only : get_ispec, nugroup, dnugroup, Hz2MeV, etafactor

    implicit none

    integer, intent(in) :: g
    integer, intent(out) :: inu
    double precision, intent(out) :: er, der
  
    inu = get_ispec(g) + 1
    er = nugroup(g) * Hz2MeV
    der = dnugroup(g) * Hz2MeV * etafactor

  end subroutine prep_opacity

!     ============================================
!
!     opacity for a given temperature, density, and ye computed
!     by an interpolation of the precalculated opacity table
!
!     this is a simplified routine with all interpolations linear
!
!     input:   er - neutrino energy (mev)
!              temp   - temperature (mev)
!              yein  - electron fraction
!              rho (in state vector) - density     (g cm^-3)
!     output:  ab  - absorptive opacity  (1/cm)
!              sc  - scattering opacity  (1/cm)
!              delta - scattering anisotropy factor
!              eta   - emissivity
!
  subroutine get_opacity_emissivity( &
       ab, sc, delta, eta,           &
       rho, yein, temp, er, inu, comp_ab, comp_sc, comp_eta)

    implicit none

    integer, intent(in) :: inu 
    logical, intent(in) :: comp_ab, comp_sc, comp_eta
    double precision, intent(in) :: rho, yein, temp, er 
    double precision, intent(out) :: ab, sc, delta, eta
    
    integer jye, jr, jt, jf

    real*8 ye

    real*8 tl, rl, el
    real*8 frmin, frmax
    real*8 deltaye, deltar, deltat, deltaf
    real*8 r1i, r2i, dri
    real*8 t1i, t2i, dti
    real*8 ye1i, ye2i, dyei
    real*8 f1i, f2i, dfi
    
    real*8 opac1, opac2, opac3, opac4, opac5
    real*8 opac6, opac7, opac8, opac9, opac10
    real*8 opac11, opac12, opac13, opac14, opac15
    real*8 opac16, opac
    
    real*8 scat1, scat2, scat3, scat4, scat5
    real*8 scat6, scat7, scat8, scat9, scat10
    real*8 scat11, scat12, scat13, scat14, scat15
    real*8 scat16, scat
    
    real*8 delt1, delt2, delt3, delt4, delt5
    real*8 delt6, delt7, delt8, delt9, delt10
    real*8 delt11, delt12, delt13, delt14, delt15
    real*8 delt16, delt
    
    real*8 emis1, emis2, emis3, emis4, emis5
    real*8 emis6, emis7, emis8, emis9, emis10
    real*8 emis11, emis12, emis13, emis14, emis15
    real*8 emis16, emi
    
    ! LHH hack to go beyond high end of table
    real*8 yetrue

    if( (rho.gt.romax) .or. (rho.lt.romin) ) then
!       print*,'rho out of range ',romin,rho,romax,iflag
!       stop'done'
    endif
    
    if( (temp.gt.tmaxi) .or. (temp.lt.tmin) ) then
!       print*,'temp out of range ',tmin,temp,tmaxi,iflag
!       stop'done'
    endif

!           if( (yein.gt.yemax) .or. (yein.lt.yemin) ) then
!              print*,'ye out of range ',yemin,yein,yemax,iflag
!              stop'done'
!           endif

    if (yein .gt. 1.d0) then
       print *, 'Ye grossly out of range', yein
       stop
    endif

    tl = log(temp)
    rl = log(rho)
    el = log(er)
    
    ye = yein
    
    ! LHH hack to go beyond high end of table
    yetrue = ye

    ! hack to hold ye within range of table, see also getgderivsYe below
    ye = max(min(ye, yemax), yemin)
    tl = max(min(tl, log(tmaxi)),log(tmin))
    
    frmin = log(enuk(1))
    frmax = log(enuk(ngmax))
    deltaye = (ye-yemin)/(yemax-yemin)*float(nyemax-1)
    deltar=(rl-dlog(romin))/(dlog(romax/romin))*float(nromax-1)
    deltat=(tl-dlog(tmin))/(dlog(tmaxi/tmin))*float(ntmax-1)
    deltaf=(el-frmin)/(frmax-frmin)*float(ngmax-1)
    jye = 1 + idint(deltaye)
    jr  = 1 + idint(deltar)
    jt  = 1 + idint(deltat)
    jf  = 1 + idint(deltaf)
    if(jr.lt.1) jr = 1
    if(jr.gt.(nromax-1)) jr = nromax-1
    if(jt.lt.1) jt = 1
    if(jt.gt.(ntmax-1)) jt = ntmax-1
    if(jye.lt.1) jye = 1
    if(jye.gt.(nyemax-1)) jye = nyemax-1
    if(jf.lt.1) jf = 1
    if(jf.gt.(ngmax-1)) jf = ngmax-1
    
    r1i=dlog(romin) + (dlog(romax)-dlog(romin))*dfloat(jr-1)/ &
         dfloat(nromax-1)
    r2i=dlog(romin) + (dlog(romax)-dlog(romin))*dfloat(jr)/ &
         dfloat(nromax-1)
    dri=(rl-r1i)/(r2i-r1i)
    if(jr .eq. 1) dri = 0.d0
    
    t1i=dlog(tmin) + (dlog(tmaxi)-dlog(tmin))*dfloat(jt-1)/ &
         dfloat(ntmax-1)
    t2i=dlog(tmin) + (dlog(tmaxi)-dlog(tmin))*dfloat(jt)/ &
         dfloat(ntmax-1)
    dti=(tl-t1i)/(t2i-t1i)
    if(jt .eq. 1) dti = 0.d0

    ye1i=yemin + (yemax-yemin)*dfloat(jye-1)/dfloat(nyemax-1)
    ye2i=yemin + (yemax-yemin)*dfloat(jye)/dfloat(nyemax-1)
    dyei=(ye-ye1i)/(ye2i-ye1i)
    if(jye .eq. 1) dyei = 0.d0
    
    f1i=frmin + (frmax-frmin)*dfloat(jf-1)/dfloat(ngmax-1)
    f2i=frmin + (frmax-frmin)*dfloat(jf)/dfloat(ngmax-1)
    dfi=(el-f1i)/(f2i-f1i)
    if(jf .eq. 1) dfi = 0.d0
    
    ! LHH hack to go beyond high end of table
    dyei=(yetrue-ye1i)/(ye2i-ye1i)
    
    !
    ! absorption
    !
    ab = 0.d0

    if (comp_ab) then 
       opac1 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            absopac(jr,jt,jye,jf,inu)
       opac2 = dri*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            absopac(jr+1,jt,jye,jf,inu)
       opac3 = (1.d0-dri)*dti*(1.d0-dyei)*(1.d0-dfi)* &
            absopac(jr,jt+1,jye,jf,inu)
       opac4 = (1.d0-dri)*(1.d0-dti)*dyei*(1.d0-dfi)* &
            absopac(jr,jt,jye+1,jf,inu)
       opac5 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*dfi* &
            absopac(jr,jt,jye,jf+1,inu)
       opac6 = dri*dti*(1.d0-dyei)*(1.d0-dfi)* &
            absopac(jr+1,jt+1,jye,jf,inu)
       opac7 = dri*(1.d0-dti)*dyei*(1.d0-dfi)* &
            absopac(jr+1,jt,jye+1,jf,inu)
       opac8 = dri*(1.d0-dti)*(1.d0-dyei)*dfi* &
            absopac(jr+1,jt,jye,jf+1,inu)
       opac9 = (1.d0-dri)*dti*dyei*(1.d0-dfi)* &
            absopac(jr,jt+1,jye+1,jf,inu)
       opac10 = (1.d0-dri)*dti*(1.d0-dyei)*dfi* &
            absopac(jr,jt+1,jye,jf+1,inu)
       opac11 = (1.d0-dri)*(1.d0-dti)*dyei*dfi* &
            absopac(jr,jt,jye+1,jf+1,inu)
       opac12 = dri*dti*dyei*(1.d0-dfi)*absopac(jr+1,jt+1,jye+1,jf,inu)
       opac13 = dri*dti*(1.d0-dyei)*dfi*absopac(jr+1,jt+1,jye,jf+1,inu)
       opac14 = (1.d0-dri)*dti*dyei*dfi*absopac(jr,jt+1,jye+1,jf+1,inu)
       opac15 = dri*(1.d0-dti)*dyei*dfi*absopac(jr+1,jt,jye+1,jf+1,inu)
       opac16 = dri*dti*dyei*dfi*absopac(jr+1,jt+1,jye+1,jf+1,inu)
       
       opac = opac1 + opac2 + opac3 + opac4 + opac5 + opac6 + &
            opac7 + opac8
       opac = opac + opac9 + opac10 + opac11 + opac12 + opac13 + &
            opac14 + opac15 + opac16
       opac = exp(opac)
       
       ab=opac*rho
    end if

    !
    ! scattering
    !
    sc = 0.d0
    
    if (comp_sc) then

       scat1 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            scaopac(jr,jt,jye,jf,inu)
       scat2 = dri*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            scaopac(jr+1,jt,jye,jf,inu)
       scat3 = (1.d0-dri)*dti*(1.d0-dyei)*(1.d0-dfi)* &
            scaopac(jr,jt+1,jye,jf,inu)
       scat4 = (1.d0-dri)*(1.d0-dti)*dyei*(1.d0-dfi)* &
            scaopac(jr,jt,jye+1,jf,inu)
       scat5 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*dfi* &
            scaopac(jr,jt,jye,jf+1,inu)
       scat6 = dri*dti*(1.d0-dyei)*(1.d0-dfi)* &
            scaopac(jr+1,jt+1,jye,jf,inu)
       scat7 = dri*(1.d0-dti)*dyei*(1.d0-dfi)* &
            scaopac(jr+1,jt,jye+1,jf,inu)
       scat8 = dri*(1.d0-dti)*(1.d0-dyei)*dfi* &
            scaopac(jr+1,jt,jye,jf+1,inu)
       scat9 = (1.d0-dri)*dti*dyei*(1.d0-dfi)* &
            scaopac(jr,jt+1,jye+1,jf,inu)
       scat10 = (1.d0-dri)*dti*(1.d0-dyei)*dfi* &
            scaopac(jr,jt+1,jye,jf+1,inu)
       scat11 = (1.d0-dri)*(1.d0-dti)*dyei*dfi* &
            scaopac(jr,jt,jye+1,jf+1,inu)
       scat12 = dri*dti*dyei*(1.d0-dfi)* &
            scaopac(jr+1,jt+1,jye+1,jf,inu)
       scat13 = dri*dti*(1.d0-dyei)*dfi* &
            scaopac(jr+1,jt+1,jye,jf+1,inu)
       scat14 = (1.d0-dri)*dti*dyei*dfi* &
            scaopac(jr,jt+1,jye+1,jf+1,inu)
       scat15 = dri*(1.d0-dti)*dyei*dfi* &
            scaopac(jr+1,jt,jye+1,jf+1,inu)
       scat16 = dri*dti*dyei*dfi* &
            scaopac(jr+1,jt+1,jye+1,jf+1,inu)
       
       scat = scat1 + scat2 + scat3 + scat4 + scat5 + scat6 + &
            scat7 + scat8
       scat = scat + scat9 + scat10 + scat11 + scat12 + scat13 + &
            scat14 + scat15 + scat16
       scat = exp(scat)
       sc=scat*rho
    end if
       
    !
    ! anisotropy factor
    !
    delta = 0.d0
    
    if (comp_sc) then
       
       delt1 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            sdelta(jr,jt,jye,jf,inu)
       delt2 = dri*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            sdelta(jr+1,jt,jye,jf,inu)
       delt3 = (1.d0-dri)*dti*(1.d0-dyei)*(1.d0-dfi)* &
            sdelta(jr,jt+1,jye,jf,inu)
       delt4 = (1.d0-dri)*(1.d0-dti)*dyei*(1.d0-dfi)* &
            sdelta(jr,jt,jye+1,jf,inu)
       delt5 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*dfi* &
            sdelta(jr,jt,jye,jf+1,inu)
       delt6 = dri*dti*(1.d0-dyei)*(1.d0-dfi)* &
            sdelta(jr+1,jt+1,jye,jf,inu)
       delt7 = dri*(1.d0-dti)*dyei*(1.d0-dfi)* &
            sdelta(jr+1,jt,jye+1,jf,inu)
       delt8 = dri*(1.d0-dti)*(1.d0-dyei)*dfi* &
            sdelta(jr+1,jt,jye,jf+1,inu)
       delt9 = (1.d0-dri)*dti*dyei*(1.d0-dfi)* &
            sdelta(jr,jt+1,jye+1,jf,inu)
       delt10 = (1.d0-dri)*dti*(1.d0-dyei)*dfi* &
            sdelta(jr,jt+1,jye,jf+1,inu)
       delt11 = (1.d0-dri)*(1.d0-dti)*dyei*dfi* &
            sdelta(jr,jt,jye+1,jf+1,inu)
       delt12 = dri*dti*dyei*(1.d0-dfi)* &
            sdelta(jr+1,jt+1,jye+1,jf,inu)
       delt13 = dri*dti*(1.d0-dyei)*dfi* &
            sdelta(jr+1,jt+1,jye,jf+1,inu)
       delt14 = (1.d0-dri)*dti*dyei*dfi* &
            sdelta(jr,jt+1,jye+1,jf+1,inu)
       delt15 = dri*(1.d0-dti)*dyei*dfi* &
            sdelta(jr+1,jt,jye+1,jf+1,inu)
       delt16 = dri*dti*dyei*dfi* &
            sdelta(jr+1,jt+1,jye+1,jf+1,inu)
       
       delt = delt1 + delt2 + delt3 + delt4 + delt5 + delt6 + &
            delt7 + delt8
       delt = delt + delt9 + delt10 + delt11 + delt12 + delt13 + &
            delt14 + delt15 + delt16
       delta=delt
       
    endif
    
    !
    ! emissivity
    !
    eta = 0.d0    

    if (comp_eta) then
       emis1 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            emis(jr,jt,jye,jf,inu)
       emis2 = dri*(1.d0-dti)*(1.d0-dyei)*(1.d0-dfi)* &
            emis(jr+1,jt,jye,jf,inu)
       emis3 = (1.d0-dri)*dti*(1.d0-dyei)*(1.d0-dfi)* &
            emis(jr,jt+1,jye,jf,inu)
       emis4 = (1.d0-dri)*(1.d0-dti)*dyei*(1.d0-dfi)* &
            emis(jr,jt,jye+1,jf,inu)
       emis5 = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)*dfi* &
            emis(jr,jt,jye,jf+1,inu)
       emis6 = dri*dti*(1.d0-dyei)*(1.d0-dfi)* &
            emis(jr+1,jt+1,jye,jf,inu)
       emis7 = dri*(1.d0-dti)*dyei*(1.d0-dfi)* &
            emis(jr+1,jt,jye+1,jf,inu)
       emis8 = dri*(1.d0-dti)*(1.d0-dyei)*dfi* &
            emis(jr+1,jt,jye,jf+1,inu)
       emis9 = (1.d0-dri)*dti*dyei*(1.d0-dfi)* &
            emis(jr,jt+1,jye+1,jf,inu)
       emis10 = (1.d0-dri)*dti*(1.d0-dyei)*dfi* &
            emis(jr,jt+1,jye,jf+1,inu)
       emis11 = (1.d0-dri)*(1.d0-dti)*dyei*dfi* &
            emis(jr,jt,jye+1,jf+1,inu)
       emis12 = dri*dti*dyei*(1.d0-dfi)*emis(jr+1,jt+1,jye+1,jf,inu)
       emis13 = dri*dti*(1.d0-dyei)*dfi*emis(jr+1,jt+1,jye,jf+1,inu)
       emis14 = (1.d0-dri)*dti*dyei*dfi*emis(jr,jt+1,jye+1,jf+1,inu)
       emis15 = dri*(1.d0-dti)*dyei*dfi*emis(jr+1,jt,jye+1,jf+1,inu)
       emis16 = dri*dti*dyei*dfi*emis(jr+1,jt+1,jye+1,jf+1,inu)
       
       emi = emis1 + emis2 + emis3 + emis4 + emis5 + emis6 + &
            emis7 + emis8
       emi = emi + emis9 + emis10 + emis11 + emis12 + emis13 + &
            emis14 + emis15 + emis16
       eta = exp(emi)*rho

    end if

  end subroutine get_opacity_emissivity

end module opacity_table_module


  subroutine init_opacity_table(iverb)

    use opacity_table_module

    implicit none
    integer iverb
    
    integer i, j, k, ng
    integer irecl
    integer nu, nu_max
    
    character*80 rstfil
    
    open(3, file='opacity.50.50.30.40.shen.grid.den1.param', &
         status='old')
    read(3,*) ngmax
    read(3,*) romin, romax, nromax, tmin, tmaxi, ntmax, &
         yemin, yemax, nyemax, ngmax, (enuk(i), i=1,ngmax)
    
    if (iverb .gt. 0) then
       !         write(6,*) ngmax
       !         write(6,*) romin, romax, nromax, tmin, tmaxi, ntmax,
       !     +      yemin, yemax, nyemax, ngmax, (enuk(i), i=1,ngmax)
       
       write(6,*)'ngmax = ',ngmax
       write(6,*)'romin, romax, nromax = ',romin,romax,nromax
       write(6,*)'tmin, tmaxi, ntmax = ',tmin,tmaxi,ntmax
       write(6,*)'yemin, yemax, nyemax',yemin,yemax,nyemax
       write(6,*)'enuk ='
       write(6,*)(enuk(i), i=1,ngmax)
    endif
    
    rstfil = 'opacity.50.50.30.40.shen.grid.den1.bin'
    irecl = 4*nromax*ntmax*nyemax*ngmax + 10
    open(2,file=rstfil,form='unformatted',access='direct', &
         recl=irecl,status='old')
    nu_max = 3
    do nu = 1, nu_max
       read(2, rec = nu) ((((absopac(i,j,k,ng,nu), i = 1,nromax), &
            j=1,ntmax),k=1,nyemax), ng=1,ngmax)
    enddo
    do nu = 1, nu_max
       read(2, rec =  3 + nu) ((((scaopac(i,j,k,ng,nu), &
            i = 1,nromax), j=1,ntmax), k=1,nyemax), ng=1,ngmax)
    enddo
    do nu = 1, nu_max
       read(2, rec = 6 + nu) ((((emis(i,j,k,ng,nu), i = 1,nromax), &
            j=1,ntmax), k=1,nyemax), ng=1,ngmax)
    enddo
    do nu = 1, nu_max
       read(2, rec = 9 + nu) ((((sdelta(i,j,k,ng,nu), &
            i = 1,nromax), j=1,ntmax), k=1,nyemax), ng=1,ngmax)
    enddo
    
  end subroutine init_opacity_table

