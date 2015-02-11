! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the Laplacian.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_laplac_error(tag,tagl1,tagl2,tagh1,tagh2, &
                                 set,clear, &
                                 var,varl1,varl2,varh1,varh2, &
                                 lo,hi,nd,domlo,domhi, &
                                 delta,xlo,problo,time,level)
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagh1,tagh2
      integer          :: varl1,varl2,varh1,varh2
      integer          :: lo(2), hi(2), domlo(2), domhi(2)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2)
      double precision :: var(varl1:varh1,varl2:varh2)
      double precision :: delta(2), xlo(2), problo(2), time
      integer          :: i,j

      double precision ::  delu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)
      double precision :: delua(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)
      double precision :: delu2(4), delu3(4), delu4(4)
      double precision :: num, denom, error

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8
      double precision, parameter :: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,1) =     var(i+1,j) -      var(i-1,j)
         delua(i,j,1) = abs(var(i+1,j)) + abs(var(i-1,j))
      end do
      end do

      ! d/dy
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,2) =     var(i,j+1) -      var(i,j-1)
         delua(i,j,2) = abs(var(i,j+1)) + abs(var(i,j-1))
      end do
      end do

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         ! d/dxdx
         delu2(1) =     delu(i+1,j,1)  -     delu(i-1,j,1)
         delu3(1) = abs(delu(i+1,j,1)) + abs(delu(i-1,j,1))
         delu4(1) =    delua(i+1,j,1)  +    delua(i-1,j,1)
                                                         
         ! d/dydx                                        
         delu2(2) =     delu(i,j+1,1)  -     delu(i,j-1,1)
         delu3(2) = abs(delu(i,j+1,1)) + abs(delu(i,j-1,1))
         delu4(2) =    delua(i,j+1,1)  +    delua(i,j-1,1)
                                                         
         ! d/dxdy                                        
         delu2(3) =     delu(i+1,j,2)  -     delu(i-1,j,2)
         delu3(3) = abs(delu(i+1,j,2)) + abs(delu(i-1,j,2))
         delu4(3) =    delua(i+1,j,2)  +    delua(i-1,j,2)
                                                         
         ! d/dydy                                        
         delu2(4) =     delu(i,j+1,2)  -     delu(i,j-1,2)
         delu3(4) = abs(delu(i,j+1,2)) + abs(delu(i,j-1,2))
         delu4(4) =    delua(i,j+1,2)  +    delua(i,j-1,2)

         ! compute the error
         num   =  delu2(1)**2 + delu2(2)**2 + delu2(3)**2 + delu2(4)**2
         denom = (delu3(1) + (epsil*delu4(1)+1.d-99))**2 + &
                 (delu3(2) + (epsil*delu4(2)+1.d-99))**2 + &
                 (delu3(3) + (epsil*delu4(3)+1.d-99))**2 + &
                 (delu3(4) + (epsil*delu4(4)+1.d-99))**2

         error = sqrt(num/denom)

         if (error .gt. ctore) tag(i,j)=set

      end do
      end do

      end subroutine ca_laplac_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_denerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             den,denl1,denl2,denh1,denh2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use tagging_params_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer denl1,denl2,denh1,denh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision den(denl1:denh1,denl2:denh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (den(i,j,1) .ge. denerr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(den(i+1,j,1) - den(i,j,1))
               ay = ABS(den(i,j+1,1) - den(i,j,1))
               ax = MAX(ax,ABS(den(i,j,1) - den(i-1,j,1)))
               ay = MAX(ay,ABS(den(i,j,1) - den(i,j-1,1)))
               if ( MAX(ax,ay) .ge. dengrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the temperature
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => temperature array
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine ca_temperror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             temp,templ1,templ2,temph1,temph2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use tagging_params_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer templ1,templ2,temph1,temph2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision temp(templ1:temph1,templ2:temph2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (temp(i,j,1) .ge. temperr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(temp(i+1,j,1) - temp(i,j,1))
               ay = ABS(temp(i,j+1,1) - temp(i,j,1))
               ax = MAX(ax,ABS(temp(i,j,1) - temp(i-1,j,1)))
               ay = MAX(ay,ABS(temp(i,j,1) - temp(i,j-1,1)))
               if ( MAX(ax,ay) .ge. tempgrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the pressure
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: press     => pressure array
! ::: np        => number of components in press array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_presserror(tag,tagl1,tagl2,tagh1,tagh2, &
                               set,clear, &
                               press,pressl1,pressl2, &
                                     pressh1,pressh2, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)

      use tagging_params_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagh1,tagh2
      integer pressl1,pressl2,pressh1,pressh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision press(pressl1:pressh1,pressl2:pressh2,np)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (press(i,j,1) .ge. presserr) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(press(i+1,j,1) - press(i,j,1))
                  ay = ABS(press(i,j+1,1) - press(i,j,1))
                  ax = MAX(ax,ABS(press(i,j,1) - press(i-1,j,1)))
                  ay = MAX(ay,ABS(press(i,j,1) - press(i,j-1,1)))
                  if ( MAX(ax,ay) .ge. pressgrad) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the velocity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vel       => velocity array
! ::: nv        => number of components in vel array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_velerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             vel,vell1,vell2,velh1,velh2, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use tagging_params_module
      implicit none

      integer set, clear, nv, level
      integer tagl1,tagl2,tagh1,tagh2
      integer vell1,vell2,velh1,velh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision vel(vell1:velh1,vell2:velh2,nv)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(vel(i+1,j,1) - vel(i,j,1))
               ay = ABS(vel(i,j+1,1) - vel(i,j,1))
               ax = MAX(ax,ABS(vel(i,j,1) - vel(i-1,j,1)))
               ay = MAX(ay,ABS(vel(i,j,1) - vel(i,j-1,1)))
               if ( MAX(ax,ay) .ge. velgrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the radiation
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: rad       => radiation array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_raderror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             rad,radl1,radl2, &
                                 radh1,radh2, &
                             lo,hi,nr,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use tagging_params_module
      implicit none

      integer set, clear, nr, level
      integer tagl1,tagl2,tagh1,tagh2
      integer radl1,radl2,radh1,radh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision rad(radl1:radh1,radl2:radh2,nr)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high radiation
      if (level .lt. max_raderr_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (rad(i,j,1) .ge. raderr) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

!     Tag on regions of high radiation gradient
      if (level .lt. max_radgrad_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(rad(i+1,j,1) - rad(i,j,1))
                  ay = ABS(rad(i,j+1,1) - rad(i,j,1))
                  ax = MAX(ax,ABS(rad(i,j,1) - rad(i-1,j,1)))
                  ay = MAX(ay,ABS(rad(i,j,1) - rad(i,j-1,1)))
                  if ( MAX(ax,ay) .ge. radgrad) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the entropy
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: ent       => entropy array
! ::: nd        => number of components in ent array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_enterror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             ent,entl1,entl2,enth1,enth2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use tagging_params_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer entl1,entl2,enth1,enth2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision ent(entl1:enth1,entl2:enth2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high entropy
      if (level .lt. max_enterr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (ent(i,j,1) .ge. enterr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high entropy gradient
      if (level .lt. max_entgrad_lev .and. time > 1.d-50) then
         do j = lo(2), hi(2)
            if (j.eq.domhi(2) .or. j.eq.domlo(2)) then
               cycle
            end if
            do i = lo(1), hi(1)
               if (i.eq.domhi(1)) then
                  cycle
               end if
               ax = ABS(ent(i+1,j,1) - ent(i,j,1))
               ay = ABS(ent(i,j+1,1) - ent(i,j,1))
               ax = MAX(ax,ABS(ent(i,j,1) - ent(i-1,j,1)))
               ay = MAX(ay,ABS(ent(i,j,1) - ent(i,j-1,1)))
               if ( MAX(ax,ay) .ge. entgrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on Ye
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: Ye        => Ye array
! ::: nd        => number of components in ent array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_yeerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             ye,yel1,yel2,yeh1,yeh2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer yel1,yel2,yeh1,yeh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision ye(yel1:yeh1,yel2:yeh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of low Ye
      if (level .lt. max_yeerr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (ye(i,j,1) .lt. yeerr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high Ye gradient
      if (level .lt. max_yegrad_lev) then
         do j = lo(2), hi(2)
            if (j.eq.domhi(2) .or. j.eq.domlo(2)) then
               cycle
            end if
            do i = lo(1), hi(1)
               if (i.eq.domhi(1)) then
                  cycle
               end if
               ax = ABS(ye(i+1,j,1) - ye(i,j,1))
               ay = ABS(ye(i,j+1,1) - ye(i,j,1))
               ax = MAX(ax,ABS(ye(i,j,1) - ye(i-1,j,1)))
               ay = MAX(ay,ABS(ye(i,j,1) - ye(i,j-1,1)))
               if ( MAX(ax,ay) .ge. yegrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on mass
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: grav      => grav array
! ::: nd        => number of components in ent array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_graverror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             grav,gravl1,gravl2,gravh1,gravh2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      use prob_params_module, only : center
      use fundamental_constants_module, only : Gconst
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagh1,tagl2,tagh2
      integer gravl1,gravh1,gravl2,gravh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision grav(gravl1:gravh1,gravl2:gravh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision, parameter :: Msun = 1.98892d33 

      double precision mass, x, y, r2
      integer i, j

      if (level .lt. max_masserr_lev) then
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            x = xlo(1) + (i-lo(1)+0.5d0)*delta(1) - center(1)
            y = xlo(2) + (j-lo(2)+0.5d0)*delta(2) - center(2)
            r2 = x**2 + y**2
            mass = grav(i,j,1) * r2 / (Gconst * Msun) ! grav is mag. grav
            if (mass .lt. masserr .and. mass .gt. 0.d0) then
               tag(i,j) = set
            endif
         enddo
         enddo
      endif
      
      end

