C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C     Here is a sample program that calls the EOS.  You can also use
C     this code to test the EOS at a separate point.
C
C
C***********************************************************************
C
      PROGRAM PNTTST
C
      IMPLICIT NONE
C
      INTEGER IFLG, EFLG, FFLG, SF
      DOUBLE PRECISION IPVAR(4), TOLD, INPYE, RHO, XPR, PPR, INPT
C
C                   The following include files are only needed if you
C                   want to examine the EOS outputs.  Otherwise they
C                   can remain commented out.
      INCLUDE 'eos_m4c.inc'
C      INCLUDE 'el_eos.inc'
C      INCLUDE 'maxwel.inc'
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C                   Initialize the boundary and fermi integral tables
      CALL LOADMX()
C
 10   CONTINUE
C
      WRITE(*,'(T2,A,$)') ' Enter T, rho, Ye (in nuclear units): '
      READ(*,*) INPT, RHO, INPYE


      RHO = RHO / 1.67262158d15
C
C                            Now set the EOS inputs:
C
C                   Set input flag to indicate temperature as input
      IFLG = 1
C                   Set "forcing" flag to zero
      FFLG = 0
C                   Set temperature
      IPVAR(1) = INPT
C                   Set initial guess at nuclear density
      IPVAR(2) = 0.155D0
C                   Set initial guess at proton eta
      IPVAR(3) = -15.d0
C                   Set initial guess at neutron eta
      IPVAR(4) = -10.0d0
C                   Set initial guess at exterior proton fraction
      PPR = RHO*INPYE
C
C
C
C                   Now call the EOS
      CALL INVEOS(IPVAR,TOLD,INPYE,RHO,IFLG,EFLG,FFLG,SF,XPR,PPR)
C
C                   Check to make sure call was successfull
      
      IF(SF.NE.1) THEN
        WRITE(*,*) ' '
        WRITE(*,*) ' EOS failed !!!!!'
        WRITE(*,*) ' '
      ELSE
        WRITE(*,*) ' '
        WRITE(*,*) ' EOS call was successfull! '
        WRITE(*,*) ' Pressure =  ',PTOT* 1.60217733d33
        WRITE(*,*) ' Internal energy = ',UTOT *   0.95655684d18
        WRITE(*,*) ' Internal energy  shifted = ',(UTOT +9.3d0)
     1           * 0.95655684d18
        WRITE(*,*) ' Entropy = ',STOT
        WRITE(*,*) ' Free energy = ',FTOT
        WRITE(*,*) ' Gamma c = ',GAM_S
        WRITE(*,*) ' Gamma e = ',1.d0+PTOT* 1.60217733d33/
     1          (RHO*1.67262158d15*(UTOT +9.3d0) * 0.95655684d18)
        WRITE(*,*) ' '
      ENDIF
C
      GOTO 10
C
      STOP
C
      END










