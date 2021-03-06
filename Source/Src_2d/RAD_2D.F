
#include "LO_BCTYPES.H"

#define dims(a) a l0, a l1, a h0, a h1
#define dimdec(a) dims(a)
#define dimv(a) a l0:a h0,a l1:a h1

#define tiny 1.d-50
c#define big  1.d+50

#define BIGKR 1.d25
c#define BIGKP 1.d0

#define REAL_T real*8

#include "FLD_limiter.F"

      subroutine setlayout(Den_in, Xmom_in, Eden_in, Eint_in, Temp_in, 
     @     FirstSpec, FirstAux, NUM_STATE)
      implicit none
      INCLUDE 'StateLayout.H'
      integer Den_in, Xmom_in, Eden_in, Eint_in, Temp_in
      integer FirstSpec, FirstAux, NUM_STATE
      DEN = Den_in
      XMOM = Xmom_in
      YMOM = Xmom_in+1
      ZMOM = Xmom_in+2
      EDEN = Eden_in
      EINT = Eint_in
      ITEMP = Temp_in
      IFS = FirstSpec
      IFX = FirstAux
      SVSIZE = NUM_STATE
      end

      subroutine multrs(d,
     @  dims(dbox),
     @  dims(reg),
     @  r, s)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      real*8 d(dimv(dbox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            d(i,j) = d(i,j) * r(i) * s(j)
         enddo
      enddo
      end

      subroutine sphc(r, s,
     @  dims(reg), dx)
      implicit none
      integer dimdec(reg)
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 dx(2)
      real*8 h1, h2, d1, d2
      integer i, j
      h1 = 0.5d0 * dx(1)
      h2 = 0.5d0 * dx(2)
      d1 = 1.d0 / (3.d0 * dx(1))
      d2 = 1.d0 / dx(2)
      do i = regl0, regh0
         r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
      enddo
      do j = regl1, regh1
         s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
      enddo
      end

      subroutine sphe(r, s, n,
     @  dims(reg), dx)
      implicit none
      integer dimdec(reg)
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      integer n
      real*8 dx(2)
      real*8 h1, h2, d1, d2
      integer i, j
      if (n .eq. 0) then
         do i = regl0, regh0
            r(i) = r(i)**2
         enddo
         h2 = 0.5d0 * dx(2)
         d2 = 1.d0 / dx(2)
         do j = regl1, regh1
            s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
         enddo
      else
         h1 = 0.5d0 * dx(1)
         d1 = 1.d0 / (3.d0 * dx(1))
         do i = regl0, regh0
            r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
         enddo
         do j = regl1, regh1
            s(j) = sin(s(j))
         enddo
      endif
      end

      subroutine lacoef(a,
     @  dims(abox),
     @  dims(reg),
     @  fkp, eta, etainv, r, s, c, dt, theta)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 fkp(dimv(abox))
      real*8 eta(dimv(abox))
      real*8 etainv(dimv(abox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 c, dt, theta
      integer i, j
      real*8 dtm
      dtm = 1.d0 / dt
      do j = regl1, regh1
         do i = regl0, regh0
            a(i,j) = r(i) * s(j) *
     @               (fkp(i,j) * etainv(i,j) * c + dtm) /
     @               (1.d0 - (1.d0 - theta) * eta(i,j))
         enddo
      enddo
      end

c arithmetic average, geometrically correct(?) but underestimates surface flux

#define KAVG0(a,b) (0.5d0 * (a + b + tiny))

c#define KAVG(a,b,d) KAVG0(a,b)

c harmonic average, overestimates surface flux

#define KAVG1(a,b) ((2.d0 * a * b) / (a + b + tiny) + tiny)

c#define KAVG(a,b,d) KAVG1(a,b)

c chooses arithmetic where optically thin, harmonic where optically thick,
c surface flux approximation at a thick/thin boundary

#define KAVG2(a,b,d) min(KAVG0(a,b), max(KAVG1(a,b), 4.d0 / (3.d0 * d)))

c#define KAVG(a,b,d) KAVG2(a,b,d)

#define KAVG(a,b,d) kavg(a,b,d,-1)

      real*8 function kavg(a, b, d, iopt)
      implicit none
      real*8 a, b, d
      integer iopt
      integer, save :: opt=100
      if (iopt .ge. 0) then
         opt = iopt
         if (opt .gt. 2) then
            print *, "Fortran KAVG: invalid averaging option"
         endif
         return
      endif
      if (opt .eq. 0) then
         kavg = KAVG0(a,b)
      else if (opt .eq. 1) then
         kavg = KAVG1(a,b)
      else
         kavg = KAVG2(a,b,d)
      endif
      end

      subroutine bclim(b,
     @  lambda, dims(bbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  r, s, c, dx)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n
      real*8 b(dimv(bbox))
      real*8 lambda(dimv(bbox))
      real*8 kappar(dimv(kbox))
      real*8 r(regl0:regh0+1)
      real*8 s(regl1:regh1+1)
      real*8 c, dx(2)
      real*8 kavg
      integer i, j
      real*8 kap
      if (n .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0 + 1
               kap = KAVG(kappar(i-1,j), kappar(i,j), dx(1))
               b(i,j) = r(i) * s(j) * c * lambda(i,j) / kap
            enddo
         enddo
      else
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               kap = KAVG(kappar(i,j-1), kappar(i,j), dx(2))
               b(i,j) = r(i) * s(j) * c * lambda(i,j) / kap
            enddo
         enddo
      endif
      end

      subroutine flxlim(lambda,
     @  dims(rbox),
     @  dims(reg), limiter)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer limiter
      real*8 lambda(dimv(rbox)), FLDlambda
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            lambda(i,j) = FLDlambda(lambda(i,j),limiter)
         enddo
      enddo
      end

      subroutine eddfac(efact,
     @  dims(rbox),
     @  dims(reg), limiter, n)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer n, limiter
      real*8 efact(dimv(rbox))
      integer i, j
      real*8 r, lambda, FLDlambda
      if (n .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0 + 1
               r = efact(i,j)
               lambda = FLDlambda(r,limiter)
               efact(i,j) = lambda + (lambda * r)**2
            enddo
         enddo
      else
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               r = efact(i,j)
               lambda = FLDlambda(r,limiter)
               efact(i,j) = lambda + (lambda * r)**2
            enddo
         enddo
      endif
      end

      subroutine scgrd1(r,
     @  dims(rbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  er, dx)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n
      real*8 r(dimv(rbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 dx(2)
      real*8 kavg
      integer i, j
      real*8 kap
      if (n .eq. 0) then
         do j = regl1, regh1
c x derivatives, gradient assembly:
            do i = regl0, regh0 + 1
               r(i,j) = abs(er(i,j) - er(i-1,j)) / dx(1)
            enddo
            i = regl0
            if (er(i-1,j) .eq. -1.d0) then
               r(i,j) = abs(er(i+1,j) - er(i,j)) / dx(1)
            endif
            i = regh0 + 1
            if (er(i,j) .eq. -1.d0) then
               r(i,j) = abs(er(i-1,j) - er(i-2,j)) / dx(1)
            endif
c construct coefficients:
            do i = regl0, regh0 + 1
               kap = KAVG(kappar(i-1,j), kappar(i,j), dx(1))
               r(i,j) = r(i,j) /
     @                  (kap * max(er(i-1,j), er(i,j), tiny))
            enddo
         enddo
      else
c y derivatives, gradient assembly:
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               r(i,j) = abs(er(i,j) - er(i,j-1)) / dx(2)
            enddo
         enddo
         do i = regl0, regh0
            j = regl1
            if (er(i,j-1) .eq. -1.d0) then
               r(i,j) = abs(er(i,j+1) - er(i,j)) / dx(2)
            endif
            j = regh1 + 1
            if (er(i,j) .eq. -1.d0) then
               r(i,j) = abs(er(i,j-1) - er(i,j-2)) / dx(2)
            endif
         enddo
c construct coefficients:
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               kap = KAVG(kappar(i,j-1), kappar(i,j), dx(2))
               r(i,j) = r(i,j) /
     @                  (kap * max(er(i,j-1), er(i,j), tiny))
            enddo
         enddo
      endif
      end

      subroutine scgrd2(r,
     @  dims(rbox),
     @  dims(reg),
     @  n, kappar, dims(kbox), er,
     @  dims(dbox), d, dx)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(dbox)
      integer n
      real*8 r(dimv(rbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 d(dimv(dbox))
      real*8 dx(2)
      real*8 kavg
      integer i, j
      real*8 kap
      if (n .eq. 0) then
c y derivatives:
         do j = regl1, regh1
            do i = regl0 - 1, regh0 + 1
               d(i,j) = er(i,j+1) - er(i,j-1)
            enddo
         enddo
c check y derivatives at regl1 - 1 and regh1 + 1:
         do i = regl0 - 1, regh0 + 1
            j = regl1
            if (er(i,j-1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
            endif
            j = regh1
            if (er(i,j+1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
            endif
         enddo
         do j = regl1, regh1
c check y derivatives at regl0 - 1 and regh0 + 1:
c (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
            i = regl0 - 1
            if (er(i,j-1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
            else if (er(i,j+1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
            endif
            i = regh0 + 1
            if (er(i,j-1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
            else if (er(i,j+1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
            endif
c x derivatives, gradient assembly:
            do i = regl0, regh0 + 1
               r(i,j) = ((er(i,j) - er(i-1,j)) / dx(1)) ** 2 +
     @                  ((d(i-1,j) + d(i,j)) /
     @                   (4.d0 * dx(2))) ** 2
            enddo
            i = regl0
            if (er(i-1,j) .eq. -1.d0) then
               r(i,j) = ((er(i+1,j) - er(i,j)) / dx(1)) ** 2 +
     @                  (d(i,j) /
     @                   (2.d0 * dx(2))) ** 2
            endif
            i = regh0 + 1
            if (er(i,j) .eq. -1.d0) then
               r(i,j) = ((er(i-1,j) - er(i-2,j)) / dx(1)) ** 2 +
     @                  (d(i-1,j) /
     @                   (2.d0 * dx(2))) ** 2
            endif
c construct scaled gradient:
            do i = regl0, regh0 + 1
               kap = KAVG(kappar(i-1,j), kappar(i,j), dx(1))
               r(i,j) = sqrt(r(i,j)) /
     @                  (kap * max(er(i-1,j), er(i,j), tiny))
            enddo
         enddo
      else
c x derivatives:
         do j = regl1 - 1, regh1 + 1
            do i = regl0, regh0
               d(i,j) = er(i+1,j) - er(i-1,j)
            enddo
c check x derivatives at regl0 - 1 and regh0 + 1:
            i = regl0
            if (er(i-1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
            endif
            i = regh0
            if (er(i+1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
            endif
         enddo
         do i = regl0, regh0
c check x derivatives at regl1 - 1 and regh1 + 1:
c (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
            j = regl1 - 1
            if (er(i-1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
            else if (er(i+1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
            endif
            j = regh1 + 1
            if (er(i-1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
            else if (er(i+1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
            endif
         enddo
c y derivatives, gradient assembly:
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               r(i,j) = ((er(i,j) - er(i,j-1)) / dx(2)) ** 2 +
     @                  ((d(i,j-1) + d(i,j)) /
     @                   (4.d0 * dx(1))) ** 2
            enddo
         enddo
         do i = regl0, regh0
            j = regl1
            if (er(i,j-1) .eq. -1.d0) then
               r(i,j) = ((er(i,j+1) - er(i,j)) / dx(2)) ** 2 +
     @                  (d(i,j) /
     @                   (2.d0 * dx(1))) ** 2
            endif
            j = regh1 + 1
            if (er(i,j) .eq. -1.d0) then
               r(i,j) = ((er(i,j-1) - er(i,j-2)) / dx(2)) ** 2 +
     @                  (d(i,j-1) /
     @                   (2.d0 * dx(1))) ** 2
            endif
         enddo
c construct coefficients:
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               kap = KAVG(kappar(i,j-1), kappar(i,j), dx(2))
               r(i,j) = sqrt(r(i,j)) /
     @                  (kap * max(er(i,j-1), er(i,j), tiny))
            enddo
         enddo
      endif
      end

      subroutine scgrd3(r,
     @  dims(rbox),
     @  dims(reg),
     @  n, kappar, dims(kbox), er,
     @  dims(dbox), d, dx)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(dbox)
      integer n
      real*8 r(dimv(rbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 d(dimv(dbox))
      real*8 dx(2)
      real*8 kavg
      integer i, j
      real*8 kap
      if (n .eq. 0) then
c y derivatives:
         do j = regl1, regh1
            do i = regl0 - 1, regh0 + 1
               d(i,j) = er(i,j+1) - er(i,j-1)
            enddo
         enddo
c check y derivatives at regl1 - 1 and regh1 + 1:
         do i = regl0 - 1, regh0 + 1
            j = regl1
            if (er(i,j-1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
            endif
            j = regh1
            if (er(i,j+1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
            endif
         enddo
         do j = regl1, regh1
c check y derivatives at regl0 - 1 and regh0 + 1:
c (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
            i = regl0 - 1
            if (er(i,j-1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
            else if (er(i,j+1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
            endif
            i = regh0 + 1
            if (er(i,j-1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
            else if (er(i,j+1) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
            endif
c x derivatives, gradient assembly:
            do i = regl0, regh0 + 1
               r(i,j) = ((er(i,j) - er(i-1,j)) / dx(1)) ** 2 +
     @                  ((d(i-1,j) + d(i,j)) /
     @                   (4.d0 * dx(2))) ** 2
            enddo
            i = regl0
            if (er(i-1,j) .eq. -1.d0) then
               r(i,j) = ((er(i+1,j) - er(i,j)) / dx(1)) ** 2 +
     @                  (d(i,j) /
     @                   (2.d0 * dx(2))) ** 2
            endif
            i = regh0 + 1
            if (er(i,j) .eq. -1.d0) then
               r(i,j) = ((er(i-1,j) - er(i-2,j)) / dx(1)) ** 2 +
     @                  (d(i-1,j) /
     @                   (2.d0 * dx(2))) ** 2
            endif
c construct scaled gradient:
            do i = regl0, regh0 + 1
               kap = KAVG(kappar(i-1,j), kappar(i,j), dx(1))
               r(i,j) = sqrt(r(i,j)) /
     @                  (kap * max(er(i-1,j), er(i,j),
     @                             er(i-1,j-1), er(i,j-1),
     @                             er(i-1,j+1), er(i,j+1), tiny))
            enddo
         enddo
      else
c x derivatives:
         do j = regl1 - 1, regh1 + 1
            do i = regl0, regh0
               d(i,j) = er(i+1,j) - er(i-1,j)
            enddo
c check x derivatives at regl0 - 1 and regh0 + 1:
            i = regl0
            if (er(i-1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
            endif
            i = regh0
            if (er(i+1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
            endif
         enddo
         do i = regl0, regh0
c check x derivatives at regl1 - 1 and regh1 + 1:
c (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
            j = regl1 - 1
            if (er(i-1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
            else if (er(i+1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
            endif
            j = regh1 + 1
            if (er(i-1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
            else if (er(i+1,j) .eq. -1.d0) then
               d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
            endif
         enddo
c y derivatives, gradient assembly:
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               r(i,j) = ((er(i,j) - er(i,j-1)) / dx(2)) ** 2 +
     @                  ((d(i,j-1) + d(i,j)) /
     @                   (4.d0 * dx(1))) ** 2
            enddo
         enddo
         do i = regl0, regh0
            j = regl1
            if (er(i,j-1) .eq. -1.d0) then
               r(i,j) = ((er(i,j+1) - er(i,j)) / dx(2)) ** 2 +
     @                  (d(i,j) /
     @                   (2.d0 * dx(1))) ** 2
            endif
            j = regh1 + 1
            if (er(i,j) .eq. -1.d0) then
               r(i,j) = ((er(i,j-1) - er(i,j-2)) / dx(2)) ** 2 +
     @                  (d(i,j-1) /
     @                   (2.d0 * dx(1))) ** 2
            endif
         enddo
c construct coefficients:
         do j = regl1, regh1 + 1
            do i = regl0, regh0
               kap = KAVG(kappar(i,j-1), kappar(i,j), dx(2))
               r(i,j) = sqrt(r(i,j)) /
     @                  (kap * max(er(i,j-1), er(i,j),
     @                             er(i-1,j-1), er(i-1,j),
     @                             er(i+1,j-1), er(i+1,j), tiny))
            enddo
         enddo
      endif
      end

      subroutine lrhs(rhs,
     @  dims(rbox),
     @  dims(reg),
     @  temp, fkp, eta, etainv, frhoem, frhoes, dfo,
     @  ero, dims(ebox), edot,
     @  r, s, dt, sigma, c, theta)
      implicit none
      integer dimdec(rbox)
      integer dimdec(ebox)
      integer dimdec(reg)
      real*8 rhs(dimv(rbox))
      real*8 temp(dimv(rbox))
      real*8 fkp(dimv(rbox))
      real*8 eta(dimv(rbox))
      real*8 etainv(dimv(rbox))
      real*8 frhoem(dimv(rbox))
      real*8 frhoes(dimv(rbox))
      real*8 dfo(dimv(rbox))
      real*8 ero(dimv(ebox))
      real*8 edot(dimv(rbox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 dt, sigma, c, theta
      integer i, j
      real*8 dtm, ek, bs, es, ekt
      dtm = 1.d0 / dt
      do j = regl1, regh1
         do i = regl0, regh0
            ek = fkp(i,j) * eta(i,j)
            bs = etainv(i,j) *
     @           4.d0 * sigma * fkp(i,j) * temp(i,j)**4
            es = eta(i,j) * (frhoem(i,j) - frhoes(i,j))
            ekt = (1.d0 - theta) * eta(i,j)
            rhs(i,j) = (rhs(i,j) + r(i) * s(j) *
     @                    (bs + dtm * (ero(i,j) + es) +
     @                     ek * c * edot(i,j) -
     @                     ekt * dfo(i,j))) /
     @                   (1.d0 - ekt)
         enddo
      enddo
      end

      subroutine anatw2(test,
     @  dims(reg),
     @  temp, p, xf, Tc, dx, xlo, lo)
      implicit none
      integer dimdec(reg)
      real*8 test(dimv(reg), 0:1)
      real*8 temp(dimv(reg))
      real*8 p, xf, Tc, dx(2), xlo(2)
      integer lo(2)
      integer i, j
      real*8 x, y, r2
      do j = regl1, regh1
         y = xlo(2) + dx(2) * ((j-lo(2)) + 0.5d0)
         do i = regl0, regh0
            x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5d0)
            r2 = x*x + y*y
            test(i,j,0) = Tc * max((1.d0-r2/xf**2), 0.d0)**(1.d0/p)
            test(i,j,1) = temp(i,j) - test(i,j,0)
         enddo
      enddo
      end

      subroutine cfrhoe( dims(reg),
     @  frhoe,
     @  dims(fb),
     @  state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(fb)
      integer dimdec(sb)
      real*8 frhoe(dimv(fb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
c      real*8 kin
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
c            kin = 0.5d0 * (state(i,j,XMOM)   ** 2 +
c     @                     state(i,j,XMOM+1) ** 2) /
c     @                    state(i,j,DEN)
c            frhoe(i,j) = state(i,j,EDEN) - kin
            frhoe(i,j) = state(i,j,EINT)
         enddo
      enddo
      end

c temp contains frhoe on input:

      subroutine gtemp(dims(reg),
     @  temp, dims(tb),
     @  const, em, en, 
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 temp(dimv(tb))
      real*8 const(0:1), em(0:1), en(0:1)
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 alpha, teff, ex, frhoal
      integer i, j
      if (en(0) .ge. 1.d0) then
         print *, "Bad exponent for cv calculation"
         stop
      endif
      ex = 1.d0 / (1.d0 - en(0))
      do j = regl1, regh1
         do i = regl0, regh0
            if (em(0) .eq. 0.d0) then
               alpha = const(0)
            else
               alpha = const(0) * state(i,j,DEN) ** em(0)
            endif
            frhoal = state(i,j,DEN) * alpha + tiny
            if (en(0) .eq. 0.d0) then
               temp(i,j) = temp(i,j) / frhoal
            else
               teff = max(temp(i,j), tiny)
               temp(i,j) = ((1.d0 - en(0)) * teff / frhoal) ** ex
            endif
         enddo
      enddo
      end

c temp contains temp on input:

      subroutine gcv(dims(reg),
     @  cv,dims(cbox),
     @  temp, dims(tbox),
     @  const, em, en, tf,
     @  state, dims(sbox))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(cbox)
      integer dimdec(tbox)
      integer dimdec(sbox)
      real*8 cv(dimv(cbox))
      real*8 temp(dimv(tbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 state(dimv(sbox), 0:SVSIZE-1)
      real*8 alpha, teff, frhoal
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            if (em(0) .eq. 0.d0) then
               alpha = const(0)
            else
               alpha = const(0) * state(i,j,DEN) ** em(0)
            endif
            frhoal = state(i,j,DEN) * alpha + tiny
            if (en(0) .eq. 0.d0) then
               cv(i,j) = alpha
            else
               teff = max(temp(i,j), tiny)
               teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
               cv(i,j) = alpha * teff ** (-en(0))
            endif
         enddo
      enddo
      end

c exch contains temp on input:

      subroutine cexch( dims(reg),
     @  exch, dims(xbox),
     @  er  , dims(ebox),
     @  fkp , dims(kbox),
     @  sigma, c)
      implicit none
      integer dimdec(reg)
      integer dimdec(xbox)
      integer dimdec(ebox)
      integer dimdec(kbox)
      real*8 exch(dimv(xbox))
      real*8 er  (dimv(ebox))
      real*8 fkp (dimv(kbox))
      real*8 sigma, c
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            exch(i,j) = fkp(i,j) *
     @                  (4.d0 * sigma * exch(i,j)**4
     @                    - c * er(i,j))
         enddo
      enddo
      end

      subroutine ceta2(dims(reg),
     @  eta, etainv, dims(etab),
     @  frho, dims(sb),
     @  temp, dims(tb),
     @  cv, dims(cb),
     @  fkp, dims(fb),
     @  er, dims(ebox),
     @  dtemp, dtime, sigma, c, underr, lagpla)
      implicit none
      integer dimdec(reg)
      integer dimdec(etab)
      integer dimdec(sb)
      integer dimdec(tb)
      integer dimdec(cb)
      integer dimdec(fb)
      integer dimdec(ebox)
      real*8 eta(dimv(etab))
      real*8 etainv(dimv(etab))
      real*8 frho(dimv(sb))
      real*8 temp(dimv(tb))
      real*8 cv(dimv(cb))
      real*8 fkp(dimv(fb))
      real*8 er(dimv(ebox))
      real*8 dtemp, dtime, sigma, c, underr
      integer lagpla
      real*8 d, frc, fac0, fac1, fac2
      integer i, j
      fac1 = 16.d0 * sigma * dtime
      if (lagpla .eq. 0) then
         fac0 = 0.25d0 * fac1 / dtemp
         fac2 = dtime * c / dtemp
      endif
      do j = regl1, regh1
         do i = regl0, regh0
            if (lagpla .ne. 0) then
c assume eta and fkp are the same
               d = fac1 * fkp(i,j) * temp(i,j) ** 3
            else
               d = fac0 * (eta(i,j) * (temp(i,j) + dtemp) ** 4 -
     @                     fkp(i,j) * (temp(i,j)        ) ** 4) -
     @             fac2 * (eta(i,j) - fkp(i,j)) * er(i,j)
c alternate form, sometimes worse, sometimes better:
c               d = fac1 * fkp(i,j) * temp(i,j) ** 3 +
c     @             fac0 * (eta(i,j) - fkp(i,j)) * temp(i,j) ** 4 -
c     @             fac2 * (eta(i,j) - fkp(i,j)) * er(i,j)
c analytic derivatives for specific test problem:
c               d = (1.2d+6 * sigma * temp(i,j) ** 2 +
c     @              1.d+5 * c * er(i,j) * (temp(i,j) + tiny) ** (-2)) * dtime
c another alternate form (much worse):
c               d = fac1 * fkp(i,j) * (temp(i,j) + dtemp) ** 3 +
c     @             fac0 * (eta(i,j) - fkp(i,j)) * (temp(i,j) + dtemp) ** 4 -
c     @             fac2 * (eta(i,j) - fkp(i,j)) * er(i,j)
            endif
	    frc = frho(i,j) * cv(i,j) + tiny
            eta(i,j) = d / (d + frc)
            etainv(i,j) = underr * frc / (d + frc)
            eta(i,j) = 1.d0 - etainv(i,j)
c            eta(i,j) = 1.d0 - underr * (1.d0 - eta(i,j))
         enddo
      enddo
      end

      subroutine ceup( dims(reg), relres, absres, 
     @  frhoes, dims(grd),
     @  frhoem, eta, etainv, dfo, dfn, exch,
     @  dt, theta)
      implicit none
      integer dimdec(reg)
      integer dimdec(grd)
      real*8 frhoes(dimv(grd))
      real*8 frhoem(dimv(grd))
      real*8 eta(dimv(grd))
      real*8 etainv(dimv(grd))
      real*8 dfo(dimv(grd))
      real*8 dfn(dimv(grd))
      real*8 exch(dimv(grd))
      real*8 dt, theta, relres, absres
      real*8 tmp, chg, tot
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            chg = 0.d0
            tot = 0.d0
            tmp = eta(i,j) * frhoes(i,j) +
     @        etainv(i,j) *
     @        (frhoem(i,j) -
     @         dt * ((1.d0 - theta) *
     @               (dfo(i,j) - dfn(i,j)) +
     @               exch(i,j)))
            chg = abs(tmp - frhoes(i,j))
            tot = abs(frhoes(i,j))
            frhoes(i,j) = tmp
            absres = max(absres, chg)
            relres = max(relres, chg / (tot + tiny))
         enddo
      enddo
      end

      subroutine ceupdterm( dims(reg), relres, absres, 
     @  frhoes, dims(grd),
     @  frhoem, eta, etainv, dfo, dfn, exch, dterm,
     @  dt, theta)
      implicit none
      integer dimdec(reg)
      integer dimdec(grd)
      real*8 frhoes(dimv(grd))
      real*8 frhoem(dimv(grd))
      real*8 eta(dimv(grd))
      real*8 etainv(dimv(grd))
      real*8 dfo(dimv(grd))
      real*8 dfn(dimv(grd))
      real*8 exch(dimv(grd))
      real*8 dterm(dimv(grd))
      real*8 dt, theta, relres, absres
      real*8 tmp, chg, tot
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            chg = 0.d0
            tot = 0.d0
            tmp = eta(i,j) * frhoes(i,j) +
     @        etainv(i,j) *
     @        (frhoem(i,j) -
     @         dt * ((1.d0 - theta) *
     @               (dfo(i,j) - dfn(i,j)) +
     @               exch(i,j)))
     @        + dt * dterm(i,j)
            chg = abs(tmp - frhoes(i,j))
            tot = abs(frhoes(i,j))
            frhoes(i,j) = tmp
            absres = max(absres, chg)
            relres = max(relres, chg / (tot + tiny))
         enddo
      enddo
      end

c nonconservative form based on delta B
      subroutine nceup(dims(reg), relres, absres, 
     @  frhoes, dims(grd),
     @  frhoem, eta, etainv,
     @  er, dims(ebox),
     @  dfo, dfn, temp, fkp, cv,
     @  state, dims(sb),
     @  sigma, c, dt, theta)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(grd)
      integer dimdec(sb)
      integer dimdec(ebox)
      real*8 frhoes(dimv(grd))
      real*8 frhoem(dimv(grd))
      real*8 eta(dimv(grd))
      real*8 etainv(dimv(grd))
      real*8 er(dimv(ebox))
      real*8 dfo(dimv(grd))
      real*8 dfn(dimv(grd))
      real*8 temp(dimv(grd))
      real*8 fkp(dimv(grd))
      real*8 cv(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 sigma, c, dt, theta, relres, absres
      real*8 tmp, chg, tot, exch, b, db, dbdt, frhocv
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            chg = 0.d0
            tot = 0.d0
            frhocv = state(i,j,DEN) * cv(i,j)
            dbdt = 16.d0 * sigma * temp(i,j)**3
            b = 4.d0 * sigma * temp(i,j)**4
            exch = fkp(i,j) * (b - c * er(i,j))
            tmp = eta(i,j) * frhoes(i,j) + etainv(i,j) *
     @            (frhoem(i,j) -
     @             dt * ((1.d0 - theta) *
     @                   (dfo(i,j) - dfn(i,j)) +
     @                   exch))
#if 1
            if (frhocv .gt. tiny .and. tmp .gt. frhoes(i,j)) then
               db = (tmp - frhoes(i,j)) * dbdt / frhocv
               if (b + db .le. 0.d0) then
                  print *, i, j, b, db, b+db
               endif
               tmp = ((b + db) / (4.d0 * sigma))**0.25d0
               tmp = frhoes(i,j) + frhocv * (tmp - temp(i,j))
            endif
#endif
            chg = abs(tmp - frhoes(i,j))
            tot = abs(frhoes(i,j))
            frhoes(i,j) = tmp
            absres = max(absres, chg)
            relres = max(relres, chg / (tot + tiny))
         enddo
      enddo
      end

      subroutine cetot(dims(reg),
     @  state, dims(sb),
     @  frhoe, dims(fb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      integer dimdec(fb)
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 frhoe(dimv(fb))
      real*8 kin
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
c            kin = 0.5d0 * (state(i,j,XMOM)   ** 2 +
c     @                     state(i,j,XMOM+1) ** 2) /
c     @                    state(i,j,DEN)
            kin = state(i,j,EDEN) - state(i,j,EINT)
            state(i,j,EINT) = frhoe(i,j)
            state(i,j,EDEN) = frhoe(i,j) + kin
         enddo
      enddo
      end

      subroutine fkpn( dims(reg),
     @  fkp, dims(fb),
     @  const, em, en,
     @  ep, nu, tf, 
     @  temp, dims(tb),
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(fb)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 fkp(dimv(fb))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 temp(dimv(tb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 teff
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            teff = max(temp(i,j), tiny)
            teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
            fkp(i,j) = const(0) *
     @                   (state(i,j,DEN) ** em(0)) *
     @                   (teff ** (-en(0))) *
     @                   (nu ** (ep(0)))
         enddo
      enddo
      end

      subroutine rosse1( dims(reg),
     @  kappar, dims(kbox),
     @  const, em, en,
     @  ep, nu,
     @  tf, kfloor, 
     @  temp, dims(tb),
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 kappar(dimv(kbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 temp(dimv(tb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 kfloor
      real*8 kf, teff
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            teff = max(temp(i,j), tiny)
            teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
            kf = const(0) *
     @             (state(i,j,DEN) ** em(0)) *
     @             (teff ** (-en(0))) *
     @             (nu ** (ep(0)))
            kappar(i,j) = max(kf, kfloor)
         enddo
      enddo
      end

      subroutine rosse1s( dims(reg),
     @  kappar, dims(kbox),
     @  const, em, en,
     @  ep, 
     @  sconst, sem, sen,
     @  sep, 
     @  nu,
     @  tf, kfloor, 
     @  temp, dims(tb),
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 kappar(dimv(kbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 sconst(0:1), sem(0:1), sen(0:1), sep(0:1)
      real*8 temp(dimv(tb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 kfloor
      real*8 kf, teff, sct
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            teff = max(temp(i,j), tiny)
            teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
            kf = const(0) *
     @             (state(i,j,DEN) ** em(0)) *
     @             (teff ** (-en(0))) *
     @             (nu ** (ep(0)))
            sct = sconst(0) *
     @             (state(i,j,DEN) ** sem(0)) *
     @             (teff ** (-sen(0))) *
     @             (nu ** (sep(0)))
            kappar(i,j) = max(kf+sct, kfloor)
         enddo
      enddo
      end

      subroutine nfloor(dest,
     @  dims(dbox),
     @  dims(reg),
     @  nflr, flr, nvar)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      integer nvar, nflr
      real*8 dest(dimv(dbox), 0:nvar-1)
      real*8 flr
      integer i, j, n
      nflr = 0
      do n = 0, nvar-1
         do j = regl1, regh1
            do i = regl0, regh0
               if (dest(i,j,n) .lt. flr) then
                  dest(i,j,n) = flr
                  nflr = nflr + 1
               endif
            enddo
         enddo
      enddo
      end

c *********************************
c ** BEGIN MGFLD routines        **
c *********************************

      subroutine lacoefmgfld(a,
     +     dims(abox),
     +     dims(reg),
     +     kappa,
     +     dims(kbox),
     +     r, s,
     +     dt, c)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      integer dimdec(kbox)

      real*8 a(dimv(abox))
      real*8 kappa(dimv(kbox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 dt, c

      integer i, j

      do j = regl1, regh1
         do i = regl0, regh0

            a(i,j) = c*kappa(i,j) + 1.d0/dt
            a(i,j) = r(i) * s(j) * a(i,j)

         enddo
      enddo
      end

c *********************************
c ** END MGFLD routines          **
c *********************************

      subroutine rfface(fine,
     @     dims(fbox),
     @     crse,
     @     dims(cbox),
     @     idim, irat)
      implicit none
      integer dimdec(fbox)
      integer dimdec(cbox)
      real*8 fine(dimv(fbox))
      real*8 crse(dimv(cbox))
      integer idim, irat(0:1)
      integer i, j
      if (idim .eq. 0) then
         do j = fboxl1, fboxh1
            fine(fboxl0,j) = crse(cboxl0, j/irat(1)) / irat(1)
         enddo
      else
         do i = fboxl0, fboxh0
            fine(i,fboxl1) = crse(i/irat(0), cboxl1) / irat(0)
         enddo
      endif
      end


      subroutine bextrp(f, fboxl0, fboxl1, fboxh0, fboxh1,
     @      regl0, regl1, regh0, regh1)

      implicit none
      integer fboxl0, fboxl1, fboxh0, fboxh1
      integer  regl0,  regl1,  regh0,  regh1
      real*8 f(fboxl0:fboxh0,fboxl1:fboxh1)
      integer i, j

c     i direction first:
      do j = regl1, regh1
         i = regl0
         f(i-1,j) = 2.d0 * f(i,j) - f(i+1,j)
         i = regh0
         f(i+1,j) = 2.d0 * f(i,j) - f(i-1,j)
      enddo

c     j direction second, including corners:
      do i = regl0 - 1, regh0 + 1
         j = regl1
         f(i,j-1) = 2.d0 * f(i,j) - f(i,j+1)
         j = regh1
         f(i,j+1) = 2.d0 * f(i,j) - f(i,j-1)
      enddo
      
c  corner results are the same whichever direction we extrapolate first
      end subroutine bextrp


      subroutine lbcoefna(bcoef, 
     @      bcgrp, bboxl0, bboxl1, bboxh0, bboxh1, 
     @      regl0, regl1, regh0, regh1, 
     @      spec, sboxl0, sboxl1, sboxh0, sboxh1, 
     @      idim)

      implicit none
      integer idim
      integer  regl0,  regl1,  regh0,  regh1
      integer bboxl0, bboxl1, bboxh0, bboxh1
      integer sboxl0, sboxl1, sboxh0, sboxh1
      real*8 bcoef(bboxl0:bboxh0,bboxl1:bboxh1)
      real*8 bcgrp(bboxl0:bboxh0,bboxl1:bboxh1)
      real*8 spec(sboxl0:sboxh0,sboxl1:sboxh1)
      integer i, j
      if (idim .eq. 0) then
         do j = regl1, regh1
         do i = regl0, regh0
            bcoef(i,j) = bcoef(i,j) 
     @          + 0.5d0 * (spec(i-1,j) + spec(i,j)) * bcgrp(i,j)
         enddo
         enddo
      else
         do j = regl1, regh1
         do i = regl0, regh0
            bcoef(i,j) = bcoef(i,j)
     @          + 0.5d0 * (spec(i,j-1) + spec(i,j)) * bcgrp(i,j)
         enddo
         enddo
      endif

      end subroutine lbcoefna


      subroutine ljupna(jnew, jboxl0, jboxl1, jboxh0, jboxh1, 
     @      regl0, regl1, regh0, regh1, 
     @      spec, sboxl0, sboxl1, sboxh0, sboxh1, 
     @      accel, aboxl0, aboxl1, aboxh0, aboxh1, 
     @      nTotal)

      implicit none
      integer nTotal
      integer  regl0,  regl1,  regh0, regh1
      integer jboxl0, jboxl1, jboxh0, jboxh1
      integer sboxl0, sboxl1, sboxh0, sboxh1
      integer aboxl0, aboxl1, aboxh0, aboxh1
      real*8 jnew(jboxl0:jboxh0,jboxl1:jboxh1,0:nTotal-1)
      real*8 spec(sboxl0:sboxh0,sboxl1:sboxh1,0:nTotal-1)
      real*8 accel(aboxl0:aboxh0,aboxl1:aboxh1)

      integer i, j, n
      do n = 0, nTotal - 1
      do j = regl1, regh1
      do i = regl0, regh0
         jnew(i,j,n) = jnew(i,j,n) + spec(i,j,n) * accel(i,j)
      enddo
      enddo
      enddo

      end subroutine ljupna
