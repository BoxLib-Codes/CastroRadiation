       subroutine eos(input, den_eos, temp_eos, ye_eos, 
     1     npts, p_eos, e_eos, gam1_eos, cs_eos, s_eos, 
     2     dedt_eos, dedy_eos, do_eos_diag)

       implicit none
       double precision den_eos, temp_eos, ye_eos
       double precision p_eos, e_eos, gam1_eos, cs_eos, s_eos 
       double precision dedt_eos, dedy_eos

       double precision  eint
       double precision const_hack 
c       PARAMETER(const_hack=9.3d0)
       PARAMETER(const_hack=0.0d0)

       double precision LNLOW, TLOW, YELOW, YEHI
       PARAMETER(LNLOW=-8.92D0, TLOW=0.05d0, YELOW=0.03d0, YEHI=0.55d0)
       
       double precision logbry

       integer iflg, eflg, fflg, sf
       double precision ipvar(4), told, inpye, rho, xpr, ppr, inpt
       integer jsfirst

       common /faster/ ipvar,ppr
       common /faster/ jsfirst
       
       data jsfirst / 1 /

       integer input, npts

       logical do_eos_diag

       include 'eos_m4c.inc'

       fflg = 0

       if (jsfirst.eq.1) then

          ipvar(2) =  0.155D0
          ipvar(3) = -15.d0
          ipvar(4) = -10.d0
          
          jsfirst = 0
c     write(6,*)"starting eos"
          
       endif

       rho = den_eos / 1.67262158d15
       ppr = rho*ye_eos
       told = temp_eos

       logbry = dlog10(rho)
       if (logbry .le. LNLOW) then
          write(6,*)" eos failure"
          write(6,*) input
          write(6,*)" density", den_eos, rho
          write(6,*)"Ye", ye_eos
          write(6,*)"temp", temp_eos

          write(6,*) "Density is too low."
          write(6,*) "log baryon density is ", logbry
          write(6,*) "The lower bound for log baryon density is ", LNLOW
          call flush()
          stop
       endif

       if (ye_eos .le. YELOW .or. ye_eos .ge. YEHI) then
          write(6,*)" eos failure"
          write(6,*) input
          write(6,*)" density", den_eos, rho
          write(6,*)"Ye", ye_eos
          write(6,*)"temp", temp_eos

          write(6,*) "Ye is out of bound."
          write(6,*) "The bounds for Ye are ", YELOW, YEHI
          call flush()
          stop
       endif

       if (input .eq. 1 .or. input .eq. 6) then

          iflg = 1
          
          ipvar(1) = temp_eos
          
          if (temp_eos .le. TLOW) then
             write(6,*)" eos failure"
             write(6,*) input
             write(6,*)" density", den_eos, rho
             write(6,*)"Ye", ye_eos
             write(6,*)"temp", temp_eos
             
             write(6,*) "T is too low."
             write(6,*) "The lower bound for T is ", TLOW
             call flush()
             stop
          endif

          call INVEOS(ipvar,temp_eos,ye_eos,rho,iflg,eflg,fflg,sf,xpr,ppr)
          
          p_eos = PTOT*1.60217733d33
          e_eos = (UTOT+const_hack)*0.95655684d18
          gam1_eos = GAM_S
          cs_eos = sqrt(gam1_eos*p_eos/den_eos)
          s_eos = STOT
          dedt_eos = DUDT * 0.95655684d18
          dedy_eos = DUDY * 0.95655684d18
          
       else if (input .eq. 5) then
          
          iflg = 2
          
          eint =  e_eos  / .95655684d18 - const_hack
          ipvar(1) = eint
          
          call INVEOS(ipvar,temp_eos,ye_eos,rho,iflg,eflg,fflg,sf,xpr,ppr)
          
          p_eos = PTOT*1.60217733d33
          gam1_eos = GAM_S
          
          if(gam1_eos.le.0.d0)then
             told = 5.d0
             call INVEOS(ipvar,told,ye_eos,rho,iflg,eflg,fflg,sf,xpr,ppr)
             p_eos = PTOT*1.60217733d33
             gam1_eos = GAM_S
             if(gam1_eos .le.0d0)then
                write(6,*)"gamma negative ", gam1_eos
                write(6,*)ipvar(1),ipvar(2),ipvar(3),ipvar(4)
                write(6,*)told, ye_eos, rho
                stop
             endif
          endif
          
          cs_eos = sqrt(gam1_eos*p_eos/den_eos)
          s_eos = STOT
          dedt_eos = DUDT * 0.95655684d18
          dedy_eos = DUDY * 0.95655684d18
          
       else
          
          stop "invalid eos input"
          
       endif
       
       if(sf.ne.1)then
          
          write(6,*)" eos failure"
          write(6,*) input
          write(6,*)" density", den_eos, rho
          write(6,*)"Ye", ye_eos
          write(6,*)"temp", temp_eos
          write(6,*)"eint", e_eos
          stop

       endif

       return
       end
