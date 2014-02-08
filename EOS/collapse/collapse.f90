module table_module

  ! number of elements in AB table
  integer, parameter :: numel = 16

  ! table dimensions
  integer, parameter :: nt = 300
  integer, parameter :: nr = 300
  integer, parameter :: ny = 50

  ! Swesty Z Square; requires nt=300, nr=300, ny=50
!  double precision :: r1  =  1.d0  !  4.0d0
!  double precision :: r2  =  15.0d0
!  double precision :: t1  = -1.0d0
!  double precision :: t2  =  1.6d0 !0.8d0
!  double precision :: t12 = -1.0d0
!  double precision :: t22 =  1.6d0
!  double precision :: y1  =  0.05d0
!  double precision :: y2  =  0.513d0

  ! Shen 180 180 50

!  double precision :: r1  =  1.d0  !  4.0d0
!  double precision :: r2  =  15.0d0
!  double precision :: t1  = -1.0d0
!  double precision :: t2  =  1.6d0
!  double precision :: t12 = -1.0d0
!  double precision :: t22 =  1.6d0
!  double precision :: y1  =  0.05d0
!  double precision :: y2  =  0.513d0
 
  ! Shen 300 300 50, t lower 0.01 MeV

  double precision :: r1  =  1.d0  !  4.0d0
  double precision :: r2  =  15.0d0
  double precision :: t1  = -2.0d0
  double precision :: t2  =  1.6d0
  double precision :: t12 = -2.0d0
  double precision :: t22 =  1.6d0
  double precision :: y1  =  0.05d0
  double precision :: y2  =  0.513d0

  double precision,save :: table(numel,ny,nr,nt)

 contains

subroutine collapse_init

  use bl_error_module
  use parallel
! use eos_module

  implicit none

  integer i, jy, jr, jt, irec, irecl, tndx, tndx0
  double precision, allocatable :: table1d(:)

!--------------------------------------------
!  table(1,y,r,t)  :: energy per baryon      |
!  table(2,y,r,t)  :: pressure               |
!  table(3,y,r,t)  :: entropy per baryon     |
!  table(4,y,r,t)  :: cv                     |
!  table(5,y,r,t)  :: xn                     |
!  table(6,y,r,t)  :: xp                     |
!  table(7,y,r,t)  :: xa                     |
!  table(8,y,r,t)  :: xh                     |
!  table(9,y,r,t)  :: za                     |
!  table(10,y,r,t) :: aw                     |
!  table(11,y,r,t) :: muhat                  |
!  table(12,y,r,t) :: gamma!                 |
!  table(13,y,r,t) :: dhy                    |
!  table(14,y,r,t) :: zht                    |
!  table(15,y,r,t) :: mue                    |
!  table(16,y,r,t) :: dpde                   |
!--------------------------------------------

  irecl = 8 * ny * nr * nt
  allocate(table1d(numel*nt*nr*ny))

! open and read table on the io processor only
  if ( parallel_ioprocessor() ) then

     !  open(10,file='newtable30030050.swesty.220.ZZZ.key.dat', &
     !  open(10,file='newtable18018050.shen-grid.den1.dat', &

     ! go to line 1010 if there is a problem
     open(10,file='newtable30030050.shen-grid.den2.oct16.dat', &
         form='unformatted', access='direct', recl=irecl, &
         status='old', err=1010)
     goto 1011
1010 continue
     call bl_error('EOSINIT: Failed to open newtable')
1011 continue

     do irec = 1, numel
        read(10,rec=irec)(((table(irec,jy,jr,jt), jt = 1,nt),jr = 1, nr), jy = 1, ny)
     enddo

     i = 0 
     do irec = 1, numel
     do jt = 1,nt
     do jr = 1,nr
     do jy = 1,ny
        i = i+1
        table1d(i) = table(irec,jy,jr,jt)
     end do
     end do
     end do
     end do
  end if

  !
  ! We now must broadcast the info in "table" back to all CPUs.
  !
  call parallel_bcast_dv(table1d, parallel_IOProcessorNode())

  i = 0 
  do irec = 1, numel
  do jt = 1,nt
  do jr = 1,nr
  do jy = 1,ny
     i = i+1
     table(irec,jy,jr,jt) = table1d(i)
  end do
  end do
  end do
  end do

  deallocate(table1d)

! do jr = 1, nr
!    do tndx0 = 1, 5
!       tndx = 60 * tndx0
!       write(tndx,50) r1+jr*(r2-r1)/300.d0, table(1,tndx,jr,40), table(2,tndx,jr,40), table(3,tndx,jr,40), table(12,tndx,jr,40)
!    enddo
! enddo
!50      FORMAT(99999(D15.6,X))

  if ( parallel_ioprocessor() ) &
     write(6,*) 'successfully read ',numel,' records'

  close(unit=10)

  if ( parallel_ioprocessor() ) then
     write(6,*)
     write(6,*)'finished reading eos table'
     write(6,*)'table(numel,1,1,1) = ',table(numel,1,1,1)
     write(6,*)'table(1,ny,1,1) = ',table(1,ny,1,1)
     write(6,*)'table(1,1,nr,1) = ',table(1,1,nr,1)
     write(6,*)'table(1,1,1,nt) = ',table(1,1,1,nt)
     write(6,*)
     write(6,*)
  end if

end subroutine collapse_init

subroutine get_numel(num_elements)

  integer, intent(out) :: num_elements
  num_elements = numel

end subroutine get_numel
 

!----------------------------------------------------------------------

subroutine findthis (f,nelem,temp,rho,ye,jy_out,jq_out,jr_out,input,pt_index)

! use eos_module

  implicit none

  double precision, intent(inout) :: f(:)
  integer         , intent(in   ) :: nelem,input
  integer,optional, intent(in   ) :: pt_index(:)
  double precision, intent(in   ) :: temp,rho,ye
  integer         , intent(  out) :: jy_out,jq_out,jr_out

  double precision :: rl,tl,rfrac,yfrac
  double precision :: q,p,delty,yl0,yl1,yl2,told
  double precision :: ql,alpha,beta,delta,q0,q1,q2
  double precision :: rl0,rl1,rl2,tl0,tl1,tl2,t
  integer          :: jy,jq,jr,nn
  double precision :: coeff(10)
  double precision :: pq,pp,qq
  double precision :: dy10,dy20,dy21,yefac0,yefac1,yefac2,yefac3

  rl    = dlog10(rho)
  tl    = dlog10(temp)

! if ( (rl .lt. r1) .or. (rl .gt. r2) ) then
!    print *,'LOG(DENSITY) OUT OF BOUNDS ',rl
!    print *,'LIMITS ARE: R1 R2     ',r1,r2
!    print *,'FROM CALL WITH INPUT  ',input
!    if (present(pt_index)) &
!       print *,'AT POINT              ',pt_index(:)
!    stop
! end if

  if ( (tl .lt. t1) .or. (tl .gt. t2) ) then
     print *,'TEMP OUT OF BOUNDS   ',tl
     print *,'LIMITS ARE: T1 T2    ',t1,t2
     print *,'FROM CALL WITH INPUT ',input
     if (present(pt_index)) &
        print *,'AT POINT              ',pt_index(:)
     stop
  end if

! if ( (ye .lt. y1) .or. (ye .gt. y2) ) then
!    print *,'YE   OUT OF BOUNDS   ',ye
!    print *,'LIMITS ARE: Y1 Y2    ',y1,y2
!    print *,'FROM CALL WITH INPUT ',input
!    if (present(pt_index)) &
!       print *,'AT POINT              ',pt_index(:)
!    stop
! end if

  ! Checking limits on rho
  rfrac = (rl-r1)/(r2-r1)
  delta = dble(nr-1) * rfrac
  jr = 1 + int(delta)
  jr = max(2,min(jr,nr-1))
  p  = delta - dble(jr-1)

  if (p .lt. 0.d0) p = 0.d0
  if (p .gt. 1.d0) p = 1.d0

! if (p .lt. 0.d0 .or. p .gt. 1.d0) then
!    print *,'P OUT OF BOUNDS ',p
!    print *,'LOG(DEN) LOG(TEMP) ',rl,tl
!    print *,'FROM CALL WITH INPUT ',input
!    stop
! end if

  ! Checking limits on temp (given rho)
  alpha = t1    +            (t12-t1)  * rfrac
  beta  = t2-t1 + ((t22-t12)-( t2-t1)) * rfrac
  ql    = (tl - alpha)/beta
  jq    = 1 + idint(dble(nt-1)*ql)
  jq    = max(2,min(jq,nt-1))
  q     = dble(nt-1)*ql - dble(jq-1)

  if (q .lt. 0.d0) q = 0.d0
  if (q .gt. 1.d0) q = 1.d0

! if (q .lt. 0.d0 .or. q .gt. 1.d0) then
!    print *,'Q OUT OF BOUNDS ',q
!    print *,'LOG(DEN) LOG(TEMP) ',rl,tl
!    print *,'FROM CALL WITH INPUT ',input
!    print *,'    ********      '
!    print *,'ALPHA, BETA, QL   ',alpha, beta, ql
!    print *,'JQ BEFORE MIN/MAX ', 1 + idint(dble(nt-1)*ql)
!    print *,'    ********      '
!    stop
! end if
  
  ! Set the Y-related indices
  yfrac = (ye-y1)/(y2-y1)
  delty = dble(ny-1) * yfrac
  jy = 1 + int(delty)
  jy = max(2,min(jy,ny-1))
  yl0 = y1+(y2-y1)*dble(jy-1-1)/dble(ny-1)
  yl1 = y1+(y2-y1)*dble(jy-1  )/dble(ny-1)
  yl2 = y1+(y2-y1)*dble(jy+1-1)/dble(ny-1)


  !Calculate coefficients that will be common among the 
  pq = p*q
  pp = p*p
  qq = q*q
  
  coeff(1) = 1.D0-P-Q+PQ
  coeff(2) = P-PQ
  coeff(3) = Q-PQ
  coeff(4) = PQ
  coeff(5) = 0.5D0*(QQ-Q)
  coeff(6) = 0.5D0*(PP-P)
  coeff(7) = 1.D0+PQ-PP-QQ
  coeff(8) = 0.5D0*(PP-2.D0*PQ+P)
  coeff(9) = 0.5D0*(QQ-2.D0*PQ+Q)
  coeff(10) = coeff(4)

  !coeff(1) = (1.d0-p)*(1.d0-q)
  !coeff(2) = p*(1.d0-q)
  !coeff(3) = (1.d0-p)*q
  !coeff(4) = p*q
  !coeff(5) = 0.5d0*q*(q-1.d0)
  !coeff(6) = 0.5d0*p*(p-1.d0)
  !coeff(7) = (1.d0+p*q-p*p-q*q)
  !coeff(8) = 0.5d0*p*(p-2.d0*q+1.d0)
  !coeff(9) = 0.5d0*q*(q-2.d0*p+1.d0)
  !coeff(10) = coeff(4)
  !factors for quadratic interpolation
  dy10 = 1.d0/(yl1-yl0)
  dy20 = 1.d0/(yl2-yl0)
  dy21 = 1.d0/(yl2-yl1)
  !yefac0 = (ye-yl1)*dy21
  yefac1 = (ye-yl1)*(ye-yl2)*dy10*dy20
  yefac2 = -(ye-yl0)*(ye-yl2)*dy10*dy21
  yefac3 = (ye-yl0)*(ye-yl1)*dy20*dy21


  do nn= 1,nelem

   if (nn.eq.20) then

     q0 = coeff(1)*table(nn,jy-1,jr  ,jq  ) &
         + coeff(2)*table(nn,jy-1,jr+1,jq  ) &
         + coeff(3)*table(nn,jy-1,jr  ,jq+1) &
         + coeff(4)*table(nn,jy-1,jr+1,jq+1)

     q1 = coeff(1)*table(nn,jy,jr  ,jq  ) &
         +coeff(2)*table(nn,jy,jr+1,jq  ) &
         +coeff(3)*table(nn,jy,jr  ,jq+1) &
         +coeff(4)*table(nn,jy,jr+1,jq+1)

     q2 = coeff(1)*table(nn,jy+1,jr  ,jq  ) &
         +coeff(2)*table(nn,jy+1,jr+1,jq  ) &
         +coeff(3)*table(nn,jy+1,jr  ,jq+1) &
         +coeff(4)*table(nn,jy+1,jr+1,jq+1)

!    f(nn)  = linter(ye,yl1,q1,yl2,q2)
     f(nn) = q0*yefac1+ &
             q1*yefac2+ &
             q2*yefac3

     !f(nn)  = qinter(ye,yl0,q0,yl1,q1,yl2,q2)

   else if (nn .ne. 20) then

     q0 = coeff(5)*table(nn,jy-1,jr  ,jq-1) &
          + coeff(6)*table(nn,jy-1,jr-1,jq  ) &
          + coeff(7)*table(nn,jy-1,jr  ,jq  ) &
          + coeff(8)*table(nn,jy-1,jr+1,jq  ) &
          + coeff(9)*table(nn,jy-1,jr  ,jq+1) &
          + coeff(10)*table(nn,jy-1,jr+1,jq+1)
     q1 = coeff(5)*table(nn,jy,jr,jq-1) &
          + coeff(6)*table(nn,jy,jr-1,jq) &
          + coeff(7)*table(nn,jy,jr,jq) &
          + coeff(8)*table(nn,jy,jr+1,jq) &
          + coeff(9)*table(nn,jy,jr,jq+1) &
          + coeff(10)*table(nn,jy,jr+1,jq+1)
     q2 = coeff(5)*table(nn,jy+1,jr,jq-1) &
          + coeff(6)*table(nn,jy+1,jr-1,jq) &
          + coeff(7)*table(nn,jy+1,jr,jq) &
          + coeff(8)*table(nn,jy+1,jr+1,jq) &
          + coeff(9)*table(nn,jy+1,jr,jq+1) &
          + coeff(10)*table(nn,jy+1,jr+1,jq+1)

     f(nn) = q0*yefac1+ &
             q1*yefac2+ &
             q2*yefac3
     !f(nn)  = qinter(ye,yl0,q0,yl1,q1,yl2,q2)
     if (nn.eq.5.or.nn.eq.6.or.nn.eq.7.or.nn.eq.8) &
          f(nn)  = min(max(f(nn),0.d0),1.d0)

   end if

  end do

  jy_out = jy
  jq_out = jq
  jr_out = jr

contains

  double precision function qinter(x,x1,y1,x2,y2,x3,y3)
    implicit none
    double precision, intent(in) ::  x,x1,y1,x2,y2,x3,y3
    qinter = y1*(x-x2)*(x-x3)/(x1-x2)/(x1-x3)+ &
         y2*(x-x1)*(x-x3)/(x2-x1)/(x2-x3)+ &
         y3*(x-x1)*(x-x2)/(x3-x1)/(x3-x2)
  end function qinter
  
  double precision function linter(x,x0,y0,x1,y1)
    implicit none
    double precision, intent(in) :: x0,y0,x1,y1,x
    linter = y0+(y1-y0)/(x1-x0)*(x-x0)
  end function linter

end subroutine findthis

end module table_module
