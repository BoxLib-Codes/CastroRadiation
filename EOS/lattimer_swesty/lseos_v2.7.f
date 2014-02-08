C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C
C              LATTIMER & SWESTY EOS (RELEASE # 2.7  9/1/95)
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         INVEOS
C    MODULE:       INVEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         5/23/90
C                  Bug fixed on (5/24/90) (affected only performance
C                  of code NOT the results!)
C
C
C    CALL LINE:    CALL INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,XPREV)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  T_OLD = INITIAL GUESS AT THE TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C                  IFLAG = 1 --> INPVAR IS TEMPERATURE
C                          2 --> INPVAR IS INTERNAL ENERGY
C                          3 --> INPVAR IS ENTROPY (NOT IMPLEM)
C
C    OUTPUTS       EOSFLG = 1 --> "NO NUCLEI" EOS
C                           2 --> GENERAL EOS
C                           3 --> BULK EOS FOR DENSITIES ABOVE NUCLEAR
C                  XPREV = UNUSED
C                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY
C
C
C
C 
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,
     1                  FORFLG,SF,XPREV,P_PREV)
C
C
C

C
      IMPLICIT NONE
C
      INCLUDE 'eos_m4c.inc'
C
C                         Local variables
C
      DOUBLE PRECISION INP_V, INP_VO, INP_VN, UFTN, DUFTN, DT
      DOUBLE PRECISION T_OLD, T_NEW, T_TEMP, T_LB, T_UB, PERDIF
      INTEGER LOOP, SF, NEW_F
C
      RSFLAG = 1
C                         Input is the temperature; call the EOS
C                         normally and then return
      IF(IFLAG.EQ.1) THEN
        CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,
     1 XPREV,P_PREV)
        T_OLD = INPVAR(1)
        RETURN
      ENDIF
C
C
C                         The input variable must be the internal
C                         energy so calc the internal energy for
C                         the initial guess at the temperature
        INP_V = INPVAR(1)
C
        T_LB = 0.15
        T_UB = 50.0
C
        INPVAR(1) = T_OLD
        CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,
     1 XPREV,P_PREV)
CCC      CALL EOS_M1D(T_OLD,YE,BRYDNS,1,EOSFLG,XPREV,P_PREV)
C
C                         Save the value of the internal energy
      IF(IFLAG.EQ.2) THEN
        INP_VO = UTOT
      ELSEIF(IFLAG.EQ.3) THEN
        INP_VO = STOT
      ENDIF
C
C
C                         Tweak the initial guess slightly so as to
C                         get a new value of the internal energy
C
      T_NEW = 1.1*T_OLD
C
      NEW_F = 1
C
      DO 20 LOOP=1,50,1
C
        INPVAR(1) = T_NEW
        CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,
     1 XPREV,P_PREV)
CCC        CALL EOS_M1D(T_NEW,YE,BRYDNS,1,EOSFLG,XPREV,P_PREV)
C
C
        IF(SF.NE.1.AND.NEW_F.EQ.1) THEN
cc          WRITE(*,*) 'INVEOS: EOS fatally failed at try:'
cc          WRITE(*,*) T_NEW,BRYDNS,YE
          T_NEW = T_NEW-0.25
          T_LB = DMIN1(T_LB,T_NEW-1.0D-1)
          GOTO 20
        ELSEIF(SF.NE.1) THEN
           DT = 0.5*DT
           T_NEW = T_NEW+DT
        ELSE
C
          NEW_F = 0
C
C                         Save this value of the internal energy too
          IF(IFLAG.EQ.2) THEN
            INP_VN = UTOT
          ELSEIF(IFLAG.EQ.3) THEN
            INP_VN = STOT
          ENDIF
C
          IF(INP_VN.LT.INP_V) THEN
            T_LB = T_NEW
c            write(*,*) 'l @ ',t_new,inp_vn,inp_v
          ELSEIF(INP_VN.GT.INP_V) THEN
c            write(*,*) 'u @ ',t_new,inp_vn,inp_v
            T_UB = T_NEW
          ENDIF
        ENDIF
C
        UFTN = INP_VN-INP_V
C
        IF(LOOP.LT.20) THEN
C                         This is the function to be zeroed by the
C                         Newton-Raphson iteration
C
C                         Numerical derivative of the above function
C                         w.r.t. the temperature
CC          DUFTN = ((INP_VN-INP_VO)/(T_NEW-T_OLD))+1.0D-15
C
C                         Analytic derivatives
          IF(IFLAG.EQ.2) THEN
            DUFTN = DUDT
          ELSEIF(IFLAG.EQ.3) THEN
            DUFTN = DSDT
          ENDIF
C
C                         Estimated correction to temperature
          DT = UFTN/DUFTN
C
C
 10       CONTINUE
C                         Temporarily store the new temperature
          T_TEMP = T_NEW-DT
C
C                         Is the new temp within a valid range?
          IF((T_TEMP.GT.T_LB).AND.(T_TEMP.LT.T_UB)) THEN
            T_OLD = T_NEW
            INP_VO = INP_VN
            T_NEW = T_TEMP
          ELSE
C                         If not cut the step size in half & try again
            DT = 0.5*DT
            IF(T_TEMP.EQ.T_NEW) THEN
              T_OLD = T_NEW
              T_NEW = 0.5*(T_LB+T_UB)
              DT = T_NEW-T_OLD
            ELSE
              GOTO 10
            ENDIF
          ENDIF
C
        ELSE
C
          T_OLD = T_NEW
          T_NEW = 0.5*(T_LB+T_UB)
          DT = T_NEW-T_OLD
C
        ENDIF
C
C                         If relative change in T is less than 1.0e-5
C                         then quit out of loop
        IF(ABS(DT/T_NEW).LT.1.0D-5) GOTO 30
C
C                End of the Do loop
 20   CONTINUE
C
C                Didn't meet convergence criterion
      SF = 0
      WRITE(*,*) ' INVERSION OF EOS FAILED TO CONVERGE',DT,T_NEW
C
C                Met the convergence criterion!!!
 30   CONTINUE
C
c                This stuff is commented out for speed reasons; it
c                virtually never gets tripped anyway
cc      INPVAR(1) = T_OLD
cc      CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,XPREV,P_PREV)
cc      IF(IFLAG.EQ.2) THEN
cc        PERDIF = INP_V-UTOT
cc      ELSE
cc        PERDIF = INP_V-STOT
cc      ENDIF
cc      IF(ABS(PERDIF).GT.1.0D-4) THEN
cc        WRITE(*,*) 'INVEOS: FAILURE',INP_V,STOT
cc        write(*,*) uftn,dt,loop
cc        WRITE(*,*) 'TRY:',T_NEW,BRYDNS,YE
cc        SF = 0
cc        RETURN
cc      ENDIF
C 
C                Return this value for T
      INPVAR(1) = INP_V
      T_OLD = T_NEW
C
C                Time to call it quits!
 999  RETURN
C
      END
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       EOS_M4C
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         3/3/92  Model 4C modifications completed
C                  12/15/90 Modified from model 4A to include the
C                  phase boundary cutoffs and Maxwell construction
C                  boundaries.
C                  7/13/90 Modified from model 1-d to include Maxwell
C                  construction
C                  5/25/90  MODEL 1D
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C
C    CALL LINE:    CALL EOS_M4C(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,
C                  XPREV,P_PREV)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C                  IFLAG = 1 --> INPVAR IS TEMPERATURE
C                          2 --> INPVAR IS INTERNAL ENERGY (NOT IMPLEM)
C                          3 --> INPVAR IS ENTROPY (NOT IMPLEMENTED)
C                          (IFLAG=1 is now assumed at this level)
C                  FFLAG = "FORCING FLAG"  0 --> NO FORCING
C                                          1 --> FORCE A PARTICULAR
C                                                SCHEME TO BE USED
C
C
C    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B
C                           2 --> GENERAL EOS
C                           3 --> BULK EOS (includes alpha's)
C                  XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON
C                          FIRST CALL)
C                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY (MUST BE
C                          SUPPLIED ON FIRST CALL)
C
C
C
C 
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EOS_M4C(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,SSFLAG,
     1                   XPREV,P_PREV)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION OUTVAR(4)
C
C
C                       This include file contains all variable
C                       declarations.  NOTE:: no implicit typing
C                       scheme is followed in this code; if you
C                       have any doubt as to a variables type CHECK
C                       IT!!!!.  Also note that ALL variables are
C                       declared explicitly.
C
      INCLUDE 'eos_m4c.inc'
C
C
C                         Set the "switch" flag to zero
      SWTFLG = 0
C
C                         Set T equal to the input variable (the entropy
C                         and internal energy options should go through
C                         INVEOS untill further notice)
      T = INPVAR(1)
C
C
C                         If the "forcing" flag is set then skip
C                         the EOS determination logic and go straight
C                         to the EOS determined by EOSFLG
      IF(FFLAG.EQ.1) THEN
        GOTO 10
      ELSE
C                         Otherwise let the EOS logic module determine
C                         the correct EOS to use
        CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)
      ENDIF
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                        Try NUCEOS first and if not successfull
C                        then try bulk EOS
 10   CONTINUE
      IF(EOSFLG.EQ.1) THEN
C
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C                    If the nuclear EOS failed and the reset flag is set
C                    then reset the initial guesses and try again
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)
          OUTVAR(1) = INPVAR(1)
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    Make a last ditch effort at convergence
          IF(SSFLAG.NE.1) THEN
            OUTVAR(2) = 0.155
            OUTVAR(3) = -15.0
            OUTVAR(4) = -20.0
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
          ENDIF
C
        ENDIF
C
C
C
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN
C                    Set EOS flag to full scheme
          EOSFLG = 2
C
C                    Else if fraction of nuclei is less than the minimum
C                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE
          IF(FFLAG.NE.1) THEN
C
            CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG)
C
            IF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN
              EOSFLG = 1
              WRITE(*,*) 'A2 failed at try = ',T,BRYDNS,YE
              GOTO 999
            ENDIF
C
C                    Set nuclei to bulk EOS
            EOSFLG = 3
C                    Save value of proton fraction
            P_PREV = YE*BRYDNS
C
            GOTO 999
C
          ELSE
            IF(NF_FLG.EQ.1) 
     1          WRITE(*,*) 'NUC failed at t,rho = ',t,brydns
            GOTO 999
          ENDIF
        ENDIF
C
      ENDIF
C
C
C          End of NUCEOS--BULK EOS calculations
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
C
C
C
C
C
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                            CALCULATE FULL EOS (INCLUDING NUCLEI)
      IF(EOSFLG.EQ.2) THEN
C
C                    Call the nuclear EOS
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    If the nuclear EOS failed and the reset flag is set
C                    then reset the initial guesses and try again
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN
cccc          WRITE(*,*) ' EOS_M4C:: r.i.gs.'
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)
          OUTVAR(1) = INPVAR(1)
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    Make a last ditch effort at convergence
          IF(SSFLAG.NE.1) THEN
            OUTVAR(2) = 0.155
            OUTVAR(3) = -15.0
            OUTVAR(4) = -20.0
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
          ENDIF
C
C
C
          IF(SSFLAG.NE.1) THEN
cccc            WRITE(*,*) '     r.i.gs. failure @ try: ',inpvar
            GOTO 999
          ELSE
            INPVAR(2) = OUTVAR(2)
            INPVAR(3) = OUTVAR(3)
            INPVAR(4) = OUTVAR(4)
          ENDIF
C                    Otherwise quit and return
        ELSEIF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN
          GOTO 999
        ENDIF
C
C
C
C
C                    If fraction of heavies is greater than the minimum
C                    parameter, then this EOS is OK
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN
C                    Set EOS flag to full scheme
          EOSFLG = 2
C
C                    Else if fraction of nuclei is less than the minimum
C                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE
C                    If the forcing flag is not set
          IF(FFLAG.NE.1) THEN
C                    Set nuclei to no nuclei EOS
            EOSFLG = 3
C                    Set flag to indicate switch is being made
            SWTFLG = 1
C
            WRITE(*,*) ' NUCEOS failed at try =',t,brydns,ye
            WRITE(*,*) ' where it shouldnt have; Bulk EOS was used'
            WRITE(*,*) ' IV = ',INPVAR
            WRITE(*,*) ' '
C
C                    Branch to bulk EOS
            GOTO 50
C
C                    Otherwise since forcing flag is set then declare
C                    a failure and return
          ELSE
C                      If the failure message flag is set then announce
C                      the failure
            IF(NF_FLG.EQ.1) 
     1          WRITE(*,*) 'NUC failed at t,r = ',t,brydns
            GOTO 999
          ENDIF
        ENDIF
C
      ENDIF
C                              END OF FULL EOS CALULATIONS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
C
C
C
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                              CALCULATE BULK EOS
 50   CONTINUE
      IF(EOSFLG.EQ.3) THEN
C
        CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG)
C
        IF((SSFLAG.EQ.0).AND.(FFLAG.EQ.1).AND.(NF_FLG.EQ.1)) THEN
          WRITE(*,*) 'A1 failed at t,rho = ',t,brydns
          GOTO 999
        ENDIF
C                           If this EOS was used as a result of the
C                           nuclear EOS failing then set the
C                           success flag to indicate a warning
        IF(SWTFLG.EQ.1) THEN
          SSFLAG = 2
        ENDIF
C
C                           Save the value of the proton fraction
        P_PREV = YE*BRYDNS
C
        GOTO 999
C
      ENDIF
C                END OF BULK EOS CALCULATIONS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                              CALCULATE VIA MAXWELL CONSTRUCTION
      IF(EOSFLG.EQ.4) THEN
C
        CALL MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C                 Save the value of the proton fraction
        P_PREV = YE*BRYDNS
C
C                 If Maxwell EOS failed then announce the failure
        IF(SSFLAG.NE.1) THEN
          WRITE(*,*) ' MAXWEL failed at try = '
          WRITE(*,*) T,BRYDNS,YE
        ENDIF
C
          GOTO 999
C
      ENDIF
C                END OF MAXWELL CONSTRUCTION CALCULATIONS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
 999  RETURN
C
      END
C
Cnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnu
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         NUCEOS.FOR
C
C***********************************************************************
C
C    MODULE:       NUCEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/13/90 Modified from model 1-d
C
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C    CALL LINE:    CALL NUCEOS(INPVAR,YE,BRYDNS,X_PREV,SSFLAG)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON
C                  FIRST CALL)
C                  SSFLAG = SUCCESS FLAG 0 --> FAILURE
C                                        1 --> SUCCESS
C
C
C 
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
      IMPLICIT NONE
C
C
      INCLUDE 'eos_m4c.inc'
      INCLUDE 'el_eos.inc'
C
C
C                       Function type declarations
C
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFI, FHALFO
      double precision fhalf
C
      DOUBLE PRECISION ZNG, ZPG
      INTEGER TCFLAG, ftflag
C
      INTEGER KKI,LLI
      DOUBLE PRECISION RESULT(5), R_CHECK(5)
      double precision a_tmp(5,5)
      DOUBLE PRECISION NI_MIN
C
      integer cflag, schflg
      double precision dtst1, dtst2
      double precision break, dnsi, dtmp8
      double precision dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6,dtmp7
cc      double precision tbsph, tbph, tbnh, tbspl, tbpl, tbnl
cc      double precision dbspdx, dbpdx, dbndx, dbspdu, dbpdu, dbndu
cc      double precision tsgl, tsgh, thl, thh, dsgdx, dhfdx, ds2dx,dzdx
cc      double precision dpt1dx, dpt2dx
c
      INCLUDE 'force.inc'
C
C                         Set the scheme flag to zero
      SCHFLG = 0
C
C
 5    CONTINUE
C
C
C
C
C                         Set T equal to the input variable (the entropy
C                         and internal energy options are not implemented
C                         in this version)
      T = INPVAR(1)
      NSUBI = INPVAR(2)
      ETA_PO = INPVAR(3)
      ETA_NO = INPVAR(4)
C
C
C                         Calc the quantum concentration of nucleons
      NQ = 2.36D-4*T**1.5 
C
C                         Calc the Fermi integral coefficent
      UQ = 20.721
C
      MQ = (T/UQ)**1.5
C
      KQ = ((T/UQ)**2.5)/(2.0*PI**2)
C
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2))
C
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ)
C
      IF(ETA_PO.GE.ETAMAX) ETA_PO = ETAMAX-0.1
      IF(ETA_NO.GE.ETAMAX) ETA_NO = ETAMAX-0.1
      NI_MIN = DMAX1(4.5D-2,BRYDNS)
      IF(NSUBI.LT.NI_MIN) NSUBI = NI_MIN+1.0D-3
C
      TCFLAG = 0
C
      cflag = 0
C
      NEWFLG = 1
C
C                    Start Newton-Raphson iteration here
C
C
      DO 30 I=1,MAXIT,1
C
        IT_NUM = I
C                       Set the "Negative" flag
        NGFLAG = 0
C
C
C
        NNOUT = MQ*F_1_2(ETA_NO)/(2.0*PI**2)
        NPOUT = MQ*F_1_2(ETA_PO)/(2.0*PI**2)
C
        NOUT = NPOUT+NNOUT
C
c20        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
        VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1    CC*(1.0+DD)*NOUT**DD+DELTAM)
        VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
        F32_NO = F_3_2(ETA_NO)
C
        F32_PO = F_3_2(ETA_PO)
C
c20        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
        MUN_O = T*ETA_NO+VNOUT
C
        MUP_O = T*ETA_PO+VPOUT
C
        MUALFA = 2.0*MUN_O+2.0*MUP_O+BALPHA-BPROUT*V_ALFA
C
        IF(ABS(MUALFA/T).LT.30.0) THEN
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T)
        ELSEIF((MUALFA/T).LT.-30.0) THEN
          ALFDNS = 0.0
        ELSE
          ALFDNS = 8.0*NQ*DEXP(3.0D1)
        ENDIF
C
C
C                   These statements take out the alfas if the
C                   alpha particle enable flag is not set
        IF(ALFLAG.NE.1) THEN
          ALFDNS = 0.0
          MUALFA = -300.0
        ENDIF
C
C
C
        EXALFA = 1.0-ALFDNS*V_ALFA
C
C
        BPRALF = ALFDNS*T
C
c---------------------------------------------------
C
C
C             Calculate fraction of space occupied by nuclei
        U_NUC = (BRYDNS-EXALFA*NOUT-4.0*ALFDNS)/
     1        (NSUBI-EXALFA*NOUT-4.0*ALFDNS)
C
C
C            Is volume occupied by nuclei within acceptable limits?
cc        IF((U_NUC.LT.0.0).OR.((U_NUC-1.0).GT.-1.0E-20)) THEN
cc        IF((U_NUC.LT.0.0).OR.(U_NUC.GT.0.996)) THEN
        IF((U_NUC.LT.1.0d-17).OR.(U_NUC.GT.0.996)) THEN
          NGFLAG = 1
          GOTO 29
        ENDIF
C
C
C            Volume exclusion factor due to nuclei
        EXCLU = 1.0-U_NUC
C
C
C            If calculated nucleon and alfa densities are larger
C            than the baryon density then reduce the eta's
        IF((EXCLU*EXALFA*NOUT+EXCLU*4.0*ALFDNS).GT.BRYDNS) THEN
          NGFLAG = 1
          GOTO 29
        ENDIF
C
C
C            Calculate the internal (inside nuclei) proton fraction
C
        X = (BRYDNS*YE-(1.0-U_NUC)*(EXALFA*NPOUT+2.0*ALFDNS))/
     1    (U_NUC*NSUBI)
        COMPX = 1.0-X
C
C
C            Is X within reasonable (but not necessarily correct)
C            limits? (YE may not be the lower bound on X !!!)
cccc        X_MIN = DMAX1(1.0D-2,(YE-0.05))
        X_MIN = DMAX1(1.0D-2,(0.8*YE))
cc        x_min = 0.95*ye
        IF((X.LT.X_MIN).OR.(X.GT.0.6)) THEN
          NGFLAG = 1
          GOTO 29
        ENDIF
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                     Calculate critical temperature & its X derivative
        TSC_12 = 87.76*((COMP/375.0)**0.5)*((0.155/NSUBS)**OVR3)
c
cccdebug      tsc_12 = 1.0d8
c
        TSUBC = TSC_12*X*COMPX
        DTCDX = TSC_12*(1.0-2.0*X)
        DTCDXX = -2.0*TSC_12
        H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4
C
cc        tsubc = tsc_12*0.25
cc        dtcdx = 0.0
cc        dtcdxx = 0.0
C
C
CC        TSUBC = 80.0*X*COMPX
CC        DTCDX = 80.0*(1.0-2.0*X)
CC        DTCDXX = -160.0
C
C                     If the X is such that T is greater than the
C                     critical temperature then fix NSUBI so that
C                     it lies in the bounds of acceptable parameter
C                     space
        ftflag = 0
        IF(((T.GT.TSUBC).OR.(H.LE.0.0)).AND.(SCHFLG.EQ.0)) THEN
C                       If this is an initial guess, then lower
C                       NSUBI untill we get a good X
          IF(NEWFLG.EQ.1) THEN
cc        write(*,*) ' nuc exceeded Tc'
cc        write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc
            ZNG = 2.0*(PI**2)*BRYDNS*(1.0-0.1*YE)/(1.5*MQ)
            ZPG = 2.0*(PI**2)*BRYDNS*0.1*YE/(1.5*MQ)
            IF(TCFLAG.NE.1) THEN
              TCFLAG = 1
              ETA_PO = FINV12(ZPG)-0.0
              ETA_NO = FINV12(ZNG)-0.0
            ELSE
              ETA_PO = ETA_PO-2.0/T
              ETA_NO = ETA_NO-2.0/T
              NSUBI = DMAX1(0.9*NSUBI,5.1D-2)
            ENDIF
            IF(DBFLAG.EQ.1) THEN
              WRITE(*,2000) '1',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
            ENDIF
            GOTO 30
          ELSE
C                       Otherwise go back and cut the stepsize in
C                       half since it was obviously too big
            NGFLAG = 1
            GOTO 29
          ENDIF
        ELSEIF(((T.GT.TSUBC).OR.(H.LE.0.0)).AND.(SCHFLG.EQ.1)) THEN
          ftflag = 1
          tsubc = 80.0*(0.25+0.5*ye)*(0.75-0.25*ye)
C
        ENDIF
C
C
        R_0 = (0.75/(PI*NSUBS))**OVR3
        Q = (384.0*PI*(R_0**2)*SIG_S/SYM_S)-16.0
C
C                        Calculate surface functions of the internal
C                        (nuclear) proton fraction, X
        SIGMA = 1.0/(Q+1.0/(X**3)+1.0/(COMPX**3))
        OVRX4 = (1.0/X**4)-(1.0/COMPX**4)
        DSIGDX = 3.0*(SIGMA**2)*OVRX4
        SIGSGP = DSIGDX/SIGMA
        SIGSG2 = 18.0*(SIGMA**2)*OVRX4**2-12.0*SIGMA*((1.0/X**5)+
     1  (1.0/COMPX**5))
C
C                        If T is less than critical temp then
        IF(T.LT.TSUBC) THEN
C                        Calculate the surface energy temperature factor
C                        and its X and T derivatives
          H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4
          HPRIM = -4.0*T/(TSUBC**2)+4.0*((T/TSUBC)**3)/TSUBC
          HPPRIM = -4.0/(TSUBC**2)+12.0*(T**2)/(TSUBC**4)
          DHDX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDX
          DHDXX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDXX+
     1    4.0*(-3.0*T**2/TSUBC**4+5.0*T**4/TSUBC**6)*(DTCDX**2)
          HX = DHDX/H
          DHDTDX = 8.0*(T/TSUBC**3-2.0*(T**3)/TSUBC**5)*DTCDX
C
C
C                        X independent version of TZERO
c          TZERO = 0.25*TSC_12
c          DTZDX = 0.0
c          DTZDXX = 0.0
C                        X dependent version of TZERO
c          TZERO = TSUBC
c          DTZDX = DTCDX
c          DTZDXX = DTCDXX
C
C
C
C                        Coulomb liquid correction factors and their
C                        derivatives
c          W = 1-(T/TZERO)**2
c          DWDX = 2.0*(T**2)*DTZDX/(TZERO**3)
c          DWDT = -2.0*T/(TZERO**2)
c          DWDTDX = 4.0*T*DTZDX/(TZERO**3)
c          DWDXDX = 2.0*(T**2)*
c     1    (DTZDXX/(TZERO**3)-3.0*(DTZDX**2)/(TZERO**4))
c          DWDTDT = -2.0/(TZERO**2)
C
          w = 1.0
          dwdt = 0.0
          dwdx = 0.0
          dwdtdx = 0.0
          dwdxdx = 0.0
          dwdtdt = 0.0
C
C
C
C                        Calc lattice factor & derivatives & products
C
          EXCLU = 1.0-U_NUC
          COMPU = 1.0-U_NUC
C
          DU = DMAX1(1.0D-15, (1.0-1.5*W*U_NUC**OVR3+0.5*U_NUC))
          DMU = DMAX1(1.0D-15,(1.0-1.5*W*(1.0-U_NUC+1.0E-20)**OVR3+
     1    0.5*(1.0-U_NUC)))
C
          DUP = -0.5*W*U_NUC**M2OVR3+0.5
          DMUP =-0.5*W*(1.0-U_NUC+1.0E-20)**M2OVR3+0.5
          DUPP = OVR3*W*((U_NUC+1.0D-20)**M5OVR3)
          DMUPP = OVR3*W*((1.0-U_NUC)+1.0E-20)**M5OVR3
C
C                Derivatives w.r.t. T
C
          DUT = -1.5*DWDT*U_NUC**OVR3
          DMUT = -1.5*DWDT*(1.0-U_NUC+1.0E-20)**OVR3
          DUPT = -0.5*DWDT*U_NUC**M2OVR3
          DMUPT = -0.5*DWDT*(1.0-U_NUC+1.0E-20)**M2OVR3
C
C                Derivatives w.r.t. X
C
          DUX = -1.5*DWDX*U_NUC**OVR3
          DMUX = -1.5*DWDX*(1.0-U_NUC+1.0E-20)**OVR3
          DUPX = -0.5*DWDX*U_NUC**M2OVR3
          DMUPX = -0.5*DWDX*(1.0-U_NUC+1.0E-20)**M2OVR3
C
C                Second derivatives w.r.t. X
C
          DUXX = -1.5*DWDXDX*U_NUC**OVR3
          DMUXX = -1.5*DWDXDX*(1.0-U_NUC+1.0E-20)**OVR3
C
C                Second derivatives w.r.t. T
C
          DUTT = -1.5*DWDTDT*U_NUC**OVR3
          DMUTT = -1.5*DWDTDT*(1.0-U_NUC+1.0E-20)**OVR3
C
C                Second derivatives w.r.t. X & T
C
          DUXT = -1.5*DWDTDX*U_NUC**OVR3
          DMUXT = -1.5*DWDTDX*(1.0-U_NUC+1.0E-20)**OVR3
C
C
          TMP1 = (U_NUC**2)+(COMPU**2)+0.6*(U_NUC*COMPU)**2
          TMP1P = 4.0*U_NUC-2.0+
     1    2.0*0.6*(U_NUC*COMPU**2-COMPU*U_NUC**2)
          TMP1PP = 4.0+2.0*0.6*(COMPU**2-4.0*U_NUC*COMPU+U_NUC**2)
C
          TMP2 = COMPU*(DU**OVR3)
          TMP2P = -1.0*DU**OVR3+OVR3*COMPU*(DU**M2OVR3)*DUP
          TMP2PP = -OVR23*(DU**M2OVR3)*DUP-OVR29*COMPU*
     1    (DU**M5OVR3)*DUP**2+OVR3*COMPU*(DU**M2OVR3)*DUPP
C
          TMP2T = OVR3*COMPU*(DU**M2OVR3)*DUT
          TMP2X = OVR3*COMPU*(DU**M2OVR3)*DUX
          TMP2XX = OVR3*COMPU*(DU**M2OVR3)*DUXX+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUX**2)
          TMP2TT = OVR3*COMPU*(DU**M2OVR3)*DUTT+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUT**2)
          TMP2XT = OVR3*COMPU*(DU**M2OVR3)*DUXT+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUX*DUT
          TMP2PT = -OVR3*(DU**M2OVR3)*DUT+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUT+
     2        OVR3*COMPU*(DU**M2OVR3)*DUPT
          TMP2PX = -OVR3*(DU**M2OVR3)*DUX+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUX+
     2        OVR3*COMPU*(DU**M2OVR3)*DUPX
C
C
C
          TMP3 = U_NUC*(DMU**OVR3)
          TMP3P = (DMU**OVR3)-OVR3*U_NUC*(DMU**M2OVR3)*DMUP
          TMP3PP = -OVR23*(DMU**M2OVR3)*DMUP-OVR29*U_NUC*
     1    (DMU**M5OVR3)*(DMUP**2)+OVR3*U_NUC*(DMU**M2OVR3)*DMUPP
C
          TMP3T = OVR3*U_NUC*(DMU**M2OVR3)*DMUT
          TMP3X = OVR3*U_NUC*(DMU**M2OVR3)*DMUX
          TMP3XX = OVR3*U_NUC*(DMU**M2OVR3)*DMUXX+
     1        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUX**2)
          TMP3TT = OVR3*U_NUC*(DMU**M2OVR3)*DMUTT+
     1        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUT**2)
          TMP3XT = OVR3*U_NUC*(DMU**M2OVR3)*DMUXT+
     1        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*DMUX*DMUT
          TMP3PT = OVR3*(DMU**M2OVR3)*DMUT-OVR3*M2OVR3*U_NUC*
     1      (DMU**M5OVR3)*DMUP*DMUT-OVR3*U_NUC*(DMU**M2OVR3)*DMUPT

          TMP3PX = OVR3*(DMU**M2OVR3)*DMUX-OVR3*M2OVR3*U_NUC*
     1      (DMU**M5OVR3)*DMUP*DMUX-OVR3*U_NUC*(DMU**M2OVR3)*DMUPX
C
C
C                 Combination D function
C
          SCRDU = U_NUC*COMPU*(TMP2+TMP3)/TMP1
          SCRDUT = U_NUC*COMPU*(TMP2T+TMP3T)/TMP1
          SCRDUX = U_NUC*COMPU*(TMP2X+TMP3X)/TMP1
          SCRDXX = U_NUC*COMPU*(TMP2XX+TMP3XX)/TMP1
          SCRDTT = U_NUC*COMPU*(TMP2TT+TMP3TT)/TMP1
          SCRDXT = U_NUC*COMPU*(TMP2XT+TMP3XT)/TMP1
C
          SCRD = SCRDU/U_NUC
          SCRDT = SCRDUT/U_NUC
          SCRDX = SCRDUX/U_NUC
C
          SCRD2 = SCRDU/COMPU
          SCRD2T = SCRDUT/COMPU
          SCRD2X = SCRDUX/COMPU
C
          SCRDUP = SCRD-SCRD2+U_NUC*COMPU*
     1    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)
C
          SCRDPT = SCRDT-SCRD2T+U_NUC*COMPU*
     1    ((TMP2PT+TMP3PT)/TMP1-(TMP2T+TMP3T)*TMP1P/TMP1**2)
C
          SCRDPX = SCRDX-SCRD2X+U_NUC*COMPU*
     1    ((TMP2PX+TMP3PX)/TMP1-(TMP2X+TMP3X)*TMP1P/TMP1**2)
C
          SCRDPP = (SCRDUP-SCRD)/U_NUC-(SCRD2+SCRDUP)/COMPU+
     1    (1.0-2.0*U_NUC)*
     2    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)+U_NUC*COMPU*
     3    ((TMP2PP+TMP3PP)/TMP1-2.0*(TMP2P+TMP3P)*TMP1P/TMP1**2-
     4    (TMP2+TMP3)*TMP1PP/TMP1**2+
     5    2.0*(TMP2+TMP3)*(TMP1P**2)/TMP1**3)
C
C
c
c           bubble D function
cbub          scrdu = (1.0-u_nuc)*dmu**ovr3
cbub          scrd = scrdu/u_nuc
cbub          scrd2 = dmu**ovr3
cbub          scrdup = -1.0*dmu**ovr3-
cbub     1    ovr3*(1.0-u_nuc)*dmup*dmu**m2ovr3
cbub          scrdpp = ovr23*dmup*dmu**m2ovr3-ovr29*(1.0-u_nuc)*
cbub     1    dmu**m5ovr3*dmup**2+ovr3*(1.0-u_nuc)*dmu**m2ovr3*dmupp
c
c         
c           nuclei D function
cnuc          scrdu = u_nuc*du**ovr3
cnuc          scrd = du**ovr3
cnuc          scrd2 = scrdu/(1.0-u_nuc)
cnuc          scrdup = du**ovr3+ovr3*u_nuc*dup*du**m2ovr3
cnuc          scrdpp = ovr23*dup*du**m2ovr3-ovr29*u_nuc*
cnuc     1    (du**m5ovr3)*(dup**2)+ovr3*u_nuc*(du**m2ovr3)*dupp
c
c
C
          ZETA_0 = CSSCAL*6.035204*(SIG_S*(16.0+Q))**OVR23
C
C                        Surface energy coefficent
          ZETA = ZETA_0*(H*SIGMA*X*NSUBI)**OVR23
C
C                        Derivative of Zeta w.r.t. X
          DZDT = OVR23*ZETA*HPRIM/H
C
C                        Derivative of Zeta w.r.t. X
          DZDX = OVR23*ZETA*(DHDX/H+SIGSGP+1.0/X)
C
C                        Derivative of Zeta w.r.t. NSUBI
          DZDNI = OVR23*ZETA/NSUBI
C
C
C
C                        Nuclear radius
          RSUBN = 9.0*H*SIGMA*SIG_S*(16.0D0+Q)*U_NUC*(1.0-U_NUC)/
     1    (2.0*ZETA*SCRDU)
C
C                        Nuclear volume
          VSUBN = 4.0*PI*(RSUBN**3)/3.0
C
C                        Atomic number
          A = NSUBI*VSUBN
C
C                        Now calc surface, Coulomb free energies
C
          FSUBSC = ZETA*SCRDU/BRYDNS
          FSUBS = OVR23*ZETA*SCRDU/BRYDNS
          FSUBC = OVR3*ZETA*SCRDU/BRYDNS
C
C
C
C                   Translational chemical potential
          MUSUBT = TRSCAL*
     1        T*DLOG((1.0-U_NUC)*(U_NUC*NSUBI)/(NQ*AZERO**2.5))
C
C                   Derivative of trans. chem. potential w.r.t. T
          DMUTDT = TRSCAL*(MUSUBT/T-1.5)
C
C                   Translational free energy per baryon
          FTRANS = TRSCAL*H*(MUSUBT-T)/AZERO
C
C                            if T is above the critical temperature
        ELSE
          A = 0.0
          RSUBN = 0.0
          VSUBN = 0.0
          FSUBS = 0.0
          FSUBC = 0.0
          FTRANS = 0.0
        ENDIF
C                            Calc ratio of NSUBI to NSUBS
        NRATIO = NSUBI/NSUBS
C
C
c20        VNI = 2.0*AA*NSUBI+4.0*BB*X*NSUBI+CC*(1.0+DD)*NSUBI**DD
        VNI = PVN(X*NSUBI,(1.0-X)*NSUBI)
C
c20        VPI = 2.0*AA*NSUBI+4.0*BB*(1.0-X)*NSUBI+
c20     1    CC*(1.0+DD)*NSUBI**DD+DELTAM
        VPI = PVP(X*NSUBI,(1.0-X)*NSUBI)
C
c---------------------------------------------------
C
        ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ
C
        ZPI = 2.0*(PI**2)*NSUBI*X/MQ
C
        ETA_NI = FINV12(ZNI)
C
        ETA_PI = FINV12(ZPI)
C
        MUN_I = T*ETA_NI+VNI
C
        MUP_I = T*ETA_PI+VPI
C
        F32_NI = F_3_2(ETA_NI)
C
        F32_PI = F_3_2(ETA_PI)
C
c20        PSUBI = LQ*(F32_NI+F32_PI)+
c20     1    (NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+DD*CC*NSUBI**(1.0+DD)
        PSUBI = LQ*(F32_NI+F32_PI)+PV_PR(X*NSUBI,(1.0-X)*NSUBI)
C
C
        BN = OVR23*ZETA*SCRD*(SIGSGP+HX+1.5*SCRDUX/SCRDU)*X/NSUBI-
     1  TRSCAL*(1.0-U_NUC)*(MUSUBT*(H-X*DHDX)/AZERO+X*DHDX*T/AZERO)
C
        BP = -OVR23*ZETA*SCRD*
     1 ((SIGSGP+HX+1.5*SCRDUX/SCRDU)*COMPX+1.0/X)/NSUBI-
     1 TRSCAL*(1.0-U_NUC)*
     2 (MUSUBT*(H+DHDX*COMPX)/AZERO-DHDX*T*COMPX/AZERO)
C
        BSUBP = ZETA*SCRDUP-OVR23*ZETA*SCRD-
     1        TRSCAL*U_NUC*NSUBI*H*MUSUBT/AZERO
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
cc        GPI = 2.0*FHALFI(ETA_PI)
cc        GPI = 2.0*FHALFI(2.0*(pi**2)*x*nsubi/mq)
cc        GNI = 2.0*FHALFI(ETA_NI)
cc        GNI = 2.0*FHALFI(2.0*(pi**2)*(1.0-x)*nsubi/mq)
c
cc        GPO = 2.0*FHALFO(ETA_PO)
cc        GNO = 2.0*FHALFO(ETA_NO)
C
c
        GPO = 2.0*FHALF(ETA_PO)
        GNO = 2.0*FHALF(ETA_NO)
        GPI = 2.0*FHALF(ETA_PI)
        GNI = 2.0*FHALF(ETA_NI)
C
C                  Derivatives of inside potentials
C
c20        DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
        DVPIDP = DPVPDP(X*NSUBI,(1.0-X)*NSUBI)
c20        DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
        DVPIDN = DPVPDN(X*NSUBI,(1.0-X)*NSUBI)
c20        DVNIDP = DVPIDN
        DVNIDP = DPVNDP(X*NSUBI,(1.0-X)*NSUBI)
c20        DVNIDN = DVPIDP
        DVNIDN = DPVNDN(X*NSUBI,(1.0-X)*NSUBI)
C
C                  Derivatives of outside potentials
C
c20        DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
        DVPODP = EIFLAG*DPVPDP(NPOUT,NNOUT)
c20        DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)))
        DVPODN = EIFLAG*DPVPDN(NPOUT,NNOUT)
c20        DVNODP = DVPODN
        DVNODP = EIFLAG*DPVNDP(NPOUT,NNOUT)
c20        DVNODN = DVPODP
        DVNODN = EIFLAG*DPVNDN(NPOUT,NNOUT)
C
C                  Derivatives of inside K.E. densities
C
        MSSCON = 3.0*MASSN/((HBAR*C)**2)
        DTPIDP = MSSCON*T*GPI
        DTPIDN = 0.0
        DTNIDP = 0.0
        DTNIDN = MSSCON*T*GNI
C
C                  Derivatives of outside K.E. densities
C
        DTPODP = MSSCON*T*GPO
        DTPODN = 0.0
        DTNODP = 0.0
        DTNODN = MSSCON*T*GNO
C
C
C                  Derivatives of inside chem. potentials
C
        DMPIDP = T*GPI/(X*NSUBI)+DVPIDP
        DMPIDN = DVPIDN
        DMNIDP = DVNIDP
        DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN
C
C                  Derivatives of outside chem. potentials
C
        DMPODP = T+DVPODP*NPOUT/GPO
        DMPODN = DVPODN*NNOUT/GNO
        DMNODP = DVNODP*NPOUT/GPO
        DMNODN = T+DVNODN*NNOUT/GNO
C
C                  Derivatives of inside pressure
C
        DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP
        DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN
C
C                  Derivatives of outside pressure
C
        DPODP = NPOUT*DMPODP+NNOUT*DMNODP
        DPODN = NPOUT*DMPODN+NNOUT*DMNODN
C
C                  Derivatives of alpha pressure
C
        DPADP = ALFDNS*
     1  ( (2.0-NPOUT*V_ALFA)*DMPODP+(2.0-NNOUT*V_ALFA)*DMNODP )
        DPADN = ALFDNS*
     1  ( (2.0-NPOUT*V_ALFA)*DMPODN+(2.0-NNOUT*V_ALFA)*DMNODN )
C
C
        N1 = NSUBI-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS
        N2 = NSUBI*X-EXALFA*NPOUT-2.0*ALFDNS
C
C                  Derivatives of U
C
        DUDPO = -EXCLU*(EXALFA*NPOUT/GPO+
     1           (4.0-NOUT*V_ALFA)*DPADP/T)/N1
        DUDNO = -EXCLU*(EXALFA*NNOUT/GNO+
     1           (4.0-NOUT*V_ALFA)*DPADN/T)/N1
        DUDNI = -U_NUC/N1
C
C                  Derivatives of X
C
        DXDPO = -(N2*DUDPO+EXCLU*(EXALFA*NPOUT/GPO+
     1           (2.0-NPOUT*V_ALFA)*DPADP/T))/(U_NUC*NSUBI)
        DXDNO = -(N2*DUDNO+EXCLU*(2.0-NPOUT*V_ALFA)*DPADN/T)/
     1           (U_NUC*NSUBI)
        DXDNI = (N2-X*N1)/(NSUBI*N1)
C
C                  Derivatives of B's w.r.t. NSUBI
C
        DB1DNI = TRSCAL*( -U_NUC*H*(MUSUBT+T)/AZERO )+
     1      OVR23*ZETA*(SCRDUP-OVR23*SCRD)/NSUBI
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
C
        DB2DNI = -2.0*ZETA*SCRD*TMP4/(9.0*NSUBI**2)-
     1  TRSCAL*( (COMPU*T/(AZERO*NSUBI))*(H+COMPX*DHDX) )
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)
        DB3DNI = -2.0*ZETA*SCRD*X*TMP4/(9.0*NSUBI**2)-
     1          TRSCAL*( ((COMPU*T)/(AZERO*NSUBI))*(H-X*DHDX) )
C
c
c
C                  Derivatives of B's w.r.t. X
C
        DB1DX = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*(SIGSGP+DHDX/H+1.0/X)+
     1  OVR23*ZETA*(SCRDPX-OVR23*SCRDX)-
     2  TRSCAL*( U_NUC*NSUBI*DHDX*MUSUBT/AZERO )
C
C
C
C
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
C
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+(X**(-2))+(X-1.0)*
     1  (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+1.5*SCRDXX/SCRDU-
     2  1.5*(SCRDUX/SCRDU)**2)
C
C
        DB2DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+
     1      OVR23*ZETA*SCRD*TMP5/NSUBI-TRSCAL*( 
     2      COMPU*(DHDX*MUSUBT+(DHDXX*(1.0-X)-DHDX)*(MUSUBT-T))/AZERO)
C
C
C
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X
C
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+X*
     1         (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+
     2       1.5*SCRDXX/SCRDU-1.5*(SCRDUX/SCRDU)**2)
C
        DB3DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+
     1      OVR23*ZETA*SCRD*TMP5/NSUBI-
     2      TRSCAL*( COMPU*(DHDX*T-X*DHDXX*(MUSUBT-T))/AZERO )
C
C
C
C                  Derivatives of B's w.r.t. U_NUC
C
        DB1DU = ZETA*(SCRDPP-OVR23*SCRDUP/U_NUC+OVR23*SCRD/U_NUC)-
     1  TRSCAL*( NSUBI*H*(MUSUBT+T*(1.0-2.0*U_NUC)/(1.0-U_NUC))/AZERO )
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
        TMP5 = (X-1.0)*1.5*(SCRDPX/SCRDU-SCRDUX*SCRDUP/SCRDU**2)
        DB2DU = (OVR23*ZETA*SCRD/NSUBI)*TMP4*(SCRDUP/SCRDU-1.0/U_NUC)+
     1    OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+
     1    TRSCAL*( (H*MUSUBT+DHDX*COMPX*(MUSUBT-T))/AZERO-
     2    (T*(1.0-2.0*U_NUC)/U_NUC)*(H+DHDX*COMPX)/AZERO )
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X
        TMP5 = X*1.5*(SCRDPX/SCRDU-SCRDUP*SCRDUX/SCRDU**2)
        DB3DU = OVR23*ZETA*SCRD*TMP4*(U_NUC*SCRDUP/SCRDU-1.0)/
     1 (U_NUC*NSUBI)+OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+
     2  TRSCAL*( (H*MUSUBT-X*DHDX*(MUSUBT-T))/AZERO-
     3 T*(1.0-2.0*U_NUC)*(H-X*DHDX)/(AZERO*U_NUC) )
C
C
C                      A1 derivatives
C
        DA1ID1 = X*DPIDP+(1.0-X)*DPIDN+NSUBI*(DPIDP-DPIDN)*DXDNI
        DA1ID2 = NSUBI*(DPIDP-DPIDN)*DXDPO
        DA1ID3 = NSUBI*(DPIDP-DPIDN)*DXDNO
C
        DA1OD1 = 0.0
        DA1OD2 = DPODP+DPADP
        DA1OD3 = DPODN+DPADN
C
        DB1D1 = DB1DNI+DB1DX*DXDNI+DB1DU*DUDNI
        DB1D2 = DB1DX*DXDPO+DB1DU*DUDPO
        DB1D3 = DB1DX*DXDNO+DB1DU*DUDNO
C
        DA1D1 = DA1ID1-DB1D1-DA1OD1
        DA1D2 = DA1ID2-DB1D2-DA1OD2
        DA1D3 = DA1ID3-DB1D3-DA1OD3
C
C                      A3 derivatives
C
        DA3ID1 = X*DMNIDP+(1.0-X)*DMNIDN+NSUBI*(DMNIDP-DMNIDN)*DXDNI
        DA3ID2 = NSUBI*(DMNIDP-DMNIDN)*DXDPO
        DA3ID3 = NSUBI*(DMNIDP-DMNIDN)*DXDNO
C
        DA3OD1 = 0.0
        DA3OD2 = DMNODP
        DA3OD3 = DMNODN
C
        DB3D1 = DB3DNI+DB3DX*DXDNI+DB3DU*DUDNI
        DB3D2 = DB3DX*DXDPO+DB3DU*DUDPO
        DB3D3 = DB3DX*DXDNO+DB3DU*DUDNO
C
        DA3D1 = DA3ID1-DB3D1-DA3OD1
        DA3D2 = DA3ID2-DB3D2-DA3OD2
        DA3D3 = DA3ID3-DB3D3-DA3OD3
C
C                      A2 derivatives
C
        DA2ID1 = X*DMPIDP+(1.0-X)*DMPIDN+NSUBI*(DMPIDP-DMPIDN)*DXDNI
        DA2ID2 = NSUBI*(DMPIDP-DMPIDN)*DXDPO
        DA2ID3 = NSUBI*(DMPIDP-DMPIDN)*DXDNO
C
        DA2OD1 = 0.0
        DA2OD2 = DMPODP
        DA2OD3 = DMPODN
C
        DB2D1 = DB2DNI+DB2DX*DXDNI+DB2DU*DUDNI
        DB2D2 = DB2DX*DXDPO+DB2DU*DUDPO
        DB2D3 = DB2DX*DXDNO+DB2DU*DUDNO
C
        DA2D1 = DA2ID1-DB2D1-DA2OD1
        DA2D2 = DA2ID2-DB2D2-DA2OD2
        DA2D3 = DA2ID3-DB2D3-DA2OD3
C
C
C                      Eta derivatives
C
        DNDETN = NNOUT/GNO
        DPDETP = NPOUT/GPO
C
        DA1DN = DA1D1
        DA1ETP = DA1D2
        DA1ETN = DA1D3
C
        DA2DN = DA2D1
        DA2ETP = DA2D2
        DA2ETN = DA2D3
C
        DA3DN = DA3D1
        DA3ETP = DA3D2
        DA3ETN = DA3D3
C
C
C
C
        A1 = PSUBI-BSUBP-BPROUT-BPRALF
        A2 = MUP_I-BP-MUP_O
        A3 = MUN_I-BN-MUN_O
C
C
C                          Unset the "new" flag
        NEWFLG = 0
C
        DETERM = DA1DN*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)-
     1           DA1ETP*(DA2DN*DA3ETN-DA2ETN*DA3DN)+
     2           DA1ETN*(DA2DN*DA3ETP-DA2ETP*DA3DN)
C
        DNSUBI = -1.0*(A1*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)+
     1           A2*(DA3ETP*DA1ETN-DA1ETP*DA3ETN)+
     2           A3*(DA1ETP*DA2ETN-DA1ETN*DA2ETP))/DETERM
C
C
        DETAP = -1.0*(A1*(DA2ETN*DA3DN-DA2DN*DA3ETN)+
     1          A2*(DA1DN*DA3ETN-DA1ETN*DA3DN)+
     2          A3*(DA1ETN*DA2DN-DA1DN*DA2ETN))/DETERM
C
C
        DETAN = -1.0*(A1*(DA2DN*DA3ETP-DA2ETP*DA3DN)+
     1          A2*(DA1ETP*DA3DN-DA1DN*DA3ETP)+
     2          A3*(DA1DN*DA2ETP-DA1ETP*DA2DN))/DETERM
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C
C                        Check the step size in NSUBI
        IF(ABS(DNSUBI/NSUBI).GT.0.04) THEN
          DNSUBI = 0.04*DNSUBI*NSUBI/ABS(DNSUBI)
        ENDIF
 26     CONTINUE
        NSUBIN = NSUBI+DNSUBI
        IF((NSUBIN.LT.DMAX1(4.5D-2,BRYDNS)).OR.(NSUBIN.GT.0.25)) THEN
          DNSUBI = 0.5*DNSUBI
          GOTO 26
        ENDIF
C
C                        Check the step size in ETA_PO
        IF(ABS(DETAP).GT.4.0) THEN
          DETAP = 4.0*DETAP/ABS(DETAP)
        ENDIF
 27     CONTINUE
        NETAP = ETA_PO+DETAP
        IF((NETAP.LT.-5000.0).OR.(NETAP.GT.ETAMAX)) THEN
          DETAP = 0.5*DETAP
          GOTO 27
        ENDIF
C
C                        Check the step size in ETA_NO
        IF(ABS(DETAN).GT.4.0) THEN
          DETAN = 4.0*DETAN/ABS(DETAN)
        ENDIF
 28     CONTINUE
        NETAN = ETA_NO+DETAN
        IF((NETAN.LT.-5000.0).OR.(NETAN.GT.ETAMAX)) THEN
          DETAN = 0.5*DETAN
          GOTO 28
        ENDIF
C
C
C                        Update the variables
ccc        if(i.lt.30) write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc
 1205   format(i3,1p9e21.14)
c
        NSUBI = NSUBI+DNSUBI
        ETA_PO = ETA_PO+DETAP
        ETA_NO = ETA_NO+DETAN
C
C
C
C                        If the required tolarences have been met
C                        break out of the loop
        IF((ABS(DNSUBI).LT.NSIACC).AND.(ABS(DETAP).LT.PRTACC)
     1    .AND.(ABS(DETAN).LT.NUTACC) ) THEN
          GOTO 40
        ELSE
      IF(DBFLAG.EQ.1) THEN
        WRITE(*,2000) '2',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
      ENDIF
          GOTO 30
        ENDIF
C
C
 29     CONTINUE
        IF(NEWFLG.NE.1) THEN
          cflag = cflag+1
          DNSUBI = 0.5*DNSUBI
          NSUBI = NSUBI-DNSUBI
          DETAP = 0.5*DETAP
          ETA_PO = ETA_PO-DETAP
          DETAN = 0.5*DETAN
          ETA_NO = ETA_NO-DETAN
          IF(DBFLAG.EQ.1) THEN
            WRITE(*,2000) '3',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
          ENDIF
          GOTO 30
        ELSE
          NSUBI = NSUBS
cc          ETA_PO = ETA_PO-0.5/T
cc          ETA_NO = ETA_NO-0.5/T
          ETA_PO = ETA_PO-2.0/T
          ETA_NO = ETA_NO-2.0/T
        ENDIF
C
C
      IF(DBFLAG.EQ.1) THEN
        WRITE(*,2000) '4',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
      ENDIF
 2000   FORMAT(t2,a,1x,i3,1x,f8.5,3(1X,G13.5))
C
C
 30   CONTINUE
C
C            If scheme 1 has failed try scheme 2
      if(schflg.eq.0) then
        schflg = 1
        goto 5
      endif
c
c
      SSFLAG = 0
      GOTO 999
C
C                    Branch label to break out of DO 30 iteration
 40   CONTINUE
C
C
C                    The following logic determines whether this was
C                    the correct scheme to use, and if not then which
C                    one should be used
C
      if(ftflag.ne.0) then
        ssflag = 4
        goto 999
      endif
C
C                    If calculated critical temperature is less than T,
C                    then switch to the scheme with no nuclei
      IF(T.GE.TSUBC) THEN
C                    Set flag to indicate FAILURE
        SSFLAG = 0
        GOTO 999
      ENDIF
C
C
C                    If fraction of nuclei present is zero and no switch
C                    has been made then switch to the no nuclei scheme
      IF(U_NUC.LE.0.0) THEN
C                    Set flag to indicate FAILURE
        SSFLAG = 0
        GOTO 999
      ELSEIF(U_NUC.GT.1.0) THEN
C                    Set flag to indicate FAILURE
        SSFLAG = 0
        GOTO 999
      ELSE
C                    Set flag to indicate success
        SSFLAG = 1
      ENDIF
C
C
C
C                    If eqns aren't really zeroed then fail
C
C
      IF( (ABS(A1).GT.1.0D-5).OR.(ABS(A2).GT.1.0D-5).OR.
     1    (ABS(A3).GT.1.0D-5) ) THEN
        SSFLAG = 0
cc        WRITE(*,*) ' NUCEOS: False convg; A = ',A1,A2,A3
        GOTO 999
      ENDIF
C
C
C
C
      IF(NSUBI.LT.0.05) THEN
        WRITE(*,*) 'NUCEOS:: <<WARNING>> NSUBI GETTING CLOSE TO LB'
      ENDIF
C
C
C
      ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ
C
      ZPI = 2.0*(PI**2)*NSUBI*X/MQ
C
      ETA_NI = FINV12(ZNI)
C
      ETA_PI = FINV12(ZPI)
C
      MUN_I = T*ETA_NI+VNI
C
      MUP_I = T*ETA_PI+VPI
C
      F32_NI = F_3_2(ETA_NI)
C
      F32_PI = F_3_2(ETA_PI)
C
      EXCLU = 1.0-U_NUC
      EXALFA = 1.0-ALFDNS*V_ALFA
C
C
C
C                    Calculate particle fractions
C
      XALFA = 4.0*EXCLU*ALFDNS/BRYDNS
      XNUT = NNOUT*EXCLU*EXALFA/BRYDNS
      XPROT = NPOUT*EXCLU*EXALFA/BRYDNS
      XH = 1.0-XPROT-XNUT-XALFA
      XHCHK = U_NUC*NSUBI/BRYDNS
C
      IF((XH.LT.HEAVCT).OR.(XHCHK.LT.HEAVCT)) THEN
C                    Set flag to indicate switch is being made
        SSFLAG = 0
cc        write(*,*) ' xh,xhchk = ',xh,xhchk
        GOTO 999
      ENDIF
C
      IF((XALFA.LT.0.0).OR.(XH.LT.0.0).OR.
     1   (XNUT.LT.0.0).OR.(XPROT.LT.0.0)) THEN   
        SSFLAG = 0
        write(*,*) ' Xs hnpa = ',xh,xnut,xprot,xalfa
        GOTO 999
      ENDIF
C
C
C
C
C                    Baryons
C
C
      MUPROT = MUP_O
      MUN = MUN_O
      MUHAT = MUN-MUPROT
C
C 
      IF(ABS((XH-XHCHK)/XHCHK).GT.1.0D-4) THEN
        SSFLAG = 0
        GOTO 999
CCC        WRITE(*,*) ' INCONSISTENCEY IN XH AT',T,BRYDNS,YE,XH,XHCHK
      ENDIF
C   
      NUCDNS = BRYDNS*XH
C
      TAU_PO = KQ*F32_PO
      TAU_PI = KQ*F32_PI
C
      TAU_NO = KQ*F32_NO
      TAU_NI = KQ*F32_NI
C
      IF(NOUT.GT.0.0) XOUT = NPOUT/NOUT
C
C
C                    Calculate internal energy of outside nucleons,
C                    alpha particles, and nuclei (per baryon)
C
c20      BUOUT = (EXCLU*EXALFA/BRYDNS)*( UQ*(TAU_PO+TAU_NO)+EIFLAG*
c20     1    ( (NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+
c20     2    CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )
      BUOUT = (EXCLU*EXALFA/BRYDNS)*( 
     1    UQ*(TAU_PO+TAU_NO)+EIFLAG*PV_E(NPOUT,NNOUT) )
C
c20      BUNUC = XH*( ( UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*
c20     1 (AA+4.0*BB*X*(1.0-X))+CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM )/
c20     2 NSUBI)+FSUBSC*(1.0-T*(SCRDUT/SCRDU+OVR23*HPRIM/H))+
c20     3 TRSCAL*
c20     4 (1.0-U_NUC)*XH*(FTRANS*(1.0-T*HPRIM/H)-H*(MUSUBT-2.5*T)/AZERO)
      BUNUC = XH*( (UQ*(TAU_PI+TAU_NI)+PV_E(X*NSUBI,(1.0-X)*NSUBI))/
     2 NSUBI)+FSUBSC*(1.0-T*(SCRDUT/SCRDU+OVR23*HPRIM/H))+
     3 TRSCAL*
     4 (1.0-U_NUC)*XH*(FTRANS*(1.0-T*HPRIM/H)-H*(MUSUBT-2.5*T)/AZERO)
C
C
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA)
C
      BU = BUOUT+BUALFA+BUNUC
C
C
      BSOUT = (EXCLU*EXALFA/BRYDNS)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)-
     1 NNOUT*ETA_NO-NPOUT*ETA_PO )
C
C
C                    Calculate entropy of alpha particles (per baryon)
      BSALFA = -0.25*XALFA*(MUALFA/T-2.5)
C
C
      BSNUC = XH*( (5.0*UQ/(3.0*T))*(TAU_NI+TAU_PI)-
     1 NSUBI*(1.0-X)*ETA_NI-NSUBI*X*ETA_PI )/NSUBI-
     2 FSUBSC*(SCRDUT/SCRDU+OVR23*HPRIM/H)-
     3 XH*TRSCAL*(1.0-U_NUC)*
     4 ((FTRANS*HPRIM/H)+H*(MUSUBT/T-2.5)/AZERO)
C
C                    Calculate total baryon entropy (per baryon)
      BS = BSOUT+BSNUC+BSALFA
C
C                    Calculate free energy of outside nucleons (per baryon)
      BFOUT = BUOUT-T*BSOUT
C
C                    Calculate free energy of alpha particles (per baryon)
      BFALFA = BUALFA-T*BSALFA
C
C                    Calculate free energy of nuclei (per baryon)
      BFNUC = BUNUC-T*BSNUC
C
C                    Calculate total baryon free energy (per baryon)
      BFTOT = BFOUT+BFNUC+BFALFA
C
C                    Calculate pressure due to nuclei
      BPRNUC = -ZETA*(SCRDU-U_NUC*SCRDUP)+
     1 TRSCAL*U_NUC*NSUBI*H*((1.0-U_NUC)*T-U_NUC*MUSUBT)/AZERO
C
C
C                    Calculate total baryon pressure
      BPRESS = BPROUT+BPRALF+BPRNUC
C
C
C                    Leptons & Photons
C
      CALL EL_EOS(T,YE,BRYDNS)
C
C
C
C                    Total pressure and eng/ent per baryon
C
      FBARY = BFTOT+FSUBE
      PBARY = BPRESS+EPRESS
      MUBARY = YE*MUPROT+(1.0-YE)*MUN
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN
C
      FTOT = BFTOT+FSUBE+PF
      UTOT = BU+EU+PU
      STOT = BS+ES+PS
      PTOT = BPRESS+EPRESS+PPRESS
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------
C                Derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C                 ------------------------------------
C                 !      Derivatives of exterior     !
C                 !      quantities                  !
C                 !      (w.r.t. Temp. and ETA's)    !
C                 !                                  !
C                 ------------------------------------
C
C
C                  Derivatives of exterior potentials
C                  w.r.t. particle densities
c20      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
      DVPODP = EIFLAG*DPVPDP(NPOUT,NNOUT)
c20      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)))
      DVPODN = EIFLAG*DPVPDN(NPOUT,NNOUT)
c20      DVNODP = DVPODN
      DVNODP = EIFLAG*DPVNDP(NPOUT,NNOUT)
c20      DVNODN = DVPODP
      DVNODN = EIFLAG*DPVNDN(NPOUT,NNOUT)
C
C
C                  Derviatives of exterior chem. pot. w.r.t. ETA's
C                  (at fixed T)
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C                  Derivatives of pressure potential w.r.t.
C                  particle densities
c20      DV_DPO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DPO = EIFLAG*DPVRDP(NPOUT,NNOUT)
c20      DV_DNO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DNO = EIFLAG*DPVRDN(NPOUT,NNOUT)
C
C                  Derivatives of pressure potential w.r.t. ETA's
C                  (at fixed T)
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
C                  Derivatives of outside pressure w.r.t. ETA's
C                  (at fixed T)
      DPODEP = NPOUT*T+DV_DEP
      DPODEN = NNOUT*T+DV_DEN
C
C                  Derivatives of alpha density w.r.t. ETA's
C                  (at fixed T)
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T
C
C                  Derivatives of alpha pressure w.r.t. ETA's
C                  (at fixed T)
      DPADEP = T*DNADEP
      DPADEN = T*DNADEN
C
C                  Derivatives of particle densities w.r.t. T
C                  (at fixed ETA's)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
C                  Derivatives of exterior chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C                  Derivative of pressure potential w.r.t. T
C                  (at fixed ETA's)
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
C                  (at fixed ETA's)
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C                  Derivative of alpha chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
C                  (at fixed ETA's)
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C                  Derivative of alpha particle pressure w.r.t. T
C                  (at fixed ETA's)
      DPADT = ALFDNS+T*DNADT
C
C
C                 ------------------------------------
C                 !      Derivatives of interior     !
C                 !      quantities                  !
C                 !      (w.r.t. Temp. and density)  !
C                 !                                  !
C                 ------------------------------------
C
C
C                   Derivatives of kinetic energy densities w.r.t. T
C                   (holding the number densities (X & NSUBI) fixed)
      DTPIDT =2.5*TAU_PI/T-2.25*X*NSUBI*GPI/UQ
      DTNIDT =2.5*TAU_NI/T-2.25*(1.0-X)*NSUBI*GNI/UQ
C
C                   Derivatives of pressures w.r.t. T
C                   (holding the number densities (X & NSUBI) fixed)
      DPIDT = OVR23*UQ*(DTPIDT+DTNIDT)
C
C                   Derivatives of interior chem. pot. w.r.t. T
C                   (holding the number densities (X & NSUBI) fixed)
      DMPIDT = ETA_PI-1.5*GPI
      DMNIDT = ETA_NI-1.5*GNI
C
C
C                  Derivatives of inside potentials w.r.t.
C                  interior proton and neutron densities
C                  (at fixed T)
c20      DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
      DVPIDP = DPVPDP(X*NSUBI,(1.0-X)*NSUBI)
c20      DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
      DVPIDN = DPVPDN(X*NSUBI,(1.0-X)*NSUBI)
c20      DVNIDP = DVPIDN
      DVNIDP = DPVNDP(X*NSUBI,(1.0-X)*NSUBI)
c20      DVNIDN = DVPIDP
      DVNIDN = DPVNDN(X*NSUBI,(1.0-X)*NSUBI)
C
C
C                   Derivatives of interior chemical potentials
C                   w.r.t. interior neutron and proton densities
C                  (at fixed T)
      DMPIDP = T*GPI/(X*NSUBI)+DVPIDP
      DMPIDN = DVPIDN
      DMNIDP = DVNIDP
      DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN
C
C                   Derivatives of interior pressure
C                   w.r.t. interior neutron and proton densities
C                  (at fixed T)
      DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP
      DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN
C
C
C
C
C                 ------------------------------------
C                 !      Derivatives of "B" terms    !
C                 !      from the chemical and       !
C                 !      pressure equilibrium        !
C                 !      equations                   !
C                 !                                  !
C                 !      (w.r.t. Temperature )       !
C                 !                                  !
C                 ------------------------------------
C
C
C             Derivative of term from pressure equilibrium eqn.
C
      DB1DT = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*HPRIM/H+
     1    ZETA*(SCRDPT-OVR23*SCRDT)-
     2    TRSCAL*U_NUC*NSUBI*(HPRIM*MUSUBT+H*DMUTDT)/AZERO
C
C
C             Derivative of term from proton equilibrium eqn.
C
      TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+
     1 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2
C
      DB2DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*TMP4+
     1    OVR23*ZETA*SCRDT*TMP4/NSUBI+
     2    OVR23*(ZETA*SCRD/NSUBI)*(X-1.0)*TMP5-
     3    TRSCAL*EXCLU*(DMUTDT*(H+DHDX*(1.0-X))+MUSUBT*
     4    (HPRIM+DHDTDX*(1.0-X))-DHDX*(1.0-X)-T*DHDX*(1.0-X))/AZERO
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C             Derivative of term from neutron equilibrium eqn.
C
      TMP4 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+
     1 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2    
      DB3DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*X*TMP4+
     1        OVR23*(ZETA*SCRDT/NSUBI)*X*TMP4+
     2        OVR23*(ZETA*SCRD/NSUBI)*X*TMP5-
     3        TRSCAL*EXCLU*(HPRIM*MUSUBT+H*DMUTDT-X*DHDTDX*(MUSUBT-T)-
     4        X*DHDX*(DMUTDT-1.0))/AZERO
C
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C                 ------------------------------------
C                 !      Derivatives of constraint   !
C                 !      and equilibrium equations   !
C                 !      with respect to the five    !
C                 !      compositional variables     !
C                 !      (U,x,n_i,eta_po,eta_no)     !
C                 !      and the three independent   !
C                 !      variables                   !
C                 !      (Baryon density, T, and Ye) !
C                 !                                  !
C                 ------------------------------------
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 1 (Baryon conservation)
C
      DFDOM(1,1) = NOUT*EXALFA+4.0*ALFDNS-NSUBI
C
      DFDOM(1,2) = 0.0
C
      DFDOM(1,3) = -U_NUC
C
      DFDOM(1,4) = -EXCLU*EXALFA*NPOUT/GPO+
     1             V_ALFA*DNADEP*EXCLU*NOUT-4.0*EXCLU*DNADEP
C
      DFDOM(1,5) = -EXCLU*EXALFA*NNOUT/GNO+
     1             V_ALFA*DNADEN*EXCLU*NOUT-4.0*EXCLU*DNADEN
C
C
C
      DFDL_1(1) = -1.0
C
      DFDL_2(1) = EXCLU*EXALFA*(DNPODT+DNNODT)-EXCLU*V_ALFA*NOUT*DNADT+
     1     4.0*EXCLU*DNADT
C            
      DFDL_3(1) = 0.0     
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 2 (Charge conservation)
C
      DFDOM(2,1) = EXALFA*NPOUT+2.0*ALFDNS-X*NSUBI
C
      DFDOM(2,2) = -U_NUC*NSUBI
C
      DFDOM(2,3) = -X*U_NUC
C
      DFDOM(2,4) = -EXCLU*EXALFA*NPOUT/GPO+
     1     V_ALFA*EXCLU*NPOUT*DNADEP-2.0*EXCLU*DNADEP
C
      DFDOM(2,5) = V_ALFA*EXCLU*NPOUT*DNADEN-2.0*EXCLU*DNADEN             
C
C
C
      DFDL_1(2) = -1.0*YE
C
      DFDL_2(2) = EXCLU*EXALFA*DNPODT-V_ALFA*EXCLU*NPOUT*DNADT+
     1     2.0*EXCLU*DNADT
C            
      DFDL_3(2) = -1.0*BRYDNS
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 3 (Proton chemical equilibrium)
C
      DFDOM(3,1) = -DB2DU
C
      DFDOM(3,2) = NSUBI*(DMPIDP-DMPIDN)-DB2DX
C
      DFDOM(3,3) = (1.0-X)*DMPIDN+X*DMPIDP-DB2DNI
C
      DFDOM(3,4) = -DMPDEP
C
      DFDOM(3,5) = -DMPDEN
C
      DFDL_1(3) = 0.0
      DFDL_2(3) = -1.0*(DMPIDT-DMPODT-DB2DT)
      DFDL_3(3) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 4 (Neutron chemical equilibrium)
C
      DFDOM(4,1) = -DB3DU
C
      DFDOM(4,2) = NSUBI*(DMNIDP-DMNIDN)-DB3DX
C
      DFDOM(4,3) = (1.0-X)*DMNIDN+X*DMNIDP-DB3DNI
C
      DFDOM(4,4) = -DMNDEP
C
      DFDOM(4,5) = -DMNDEN
C
      DFDL_1(4) = 0.0
      DFDL_2(4) = -1.0*(DMNIDT-DMNODT-DB3DT)
      DFDL_3(4) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 5 (Pressure equilibrium)
C
      DFDOM(5,1) = -DB1DU
C
      DFDOM(5,2) = NSUBI*(DPIDP-DPIDN)-DB1DX
C
      DFDOM(5,3) = (1.0-X)*DPIDN+X*DPIDP-DB1DNI
      ncomp = dfdom(5,3)
C
      DFDOM(5,4) = -DPODEP-DPADEP
C
      DFDOM(5,5) = -DPODEN-DPADEN
C
      DFDL_1(5) = 0.0
      DFDL_2(5) = -1.0*(DPIDT-DPODT-DPADT-DB1DT)
      DFDL_3(5) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                    LU decompose the DFDOM matrix
      CALL MATLUD(DFDOM,DFDMLU,5,IPVT)
C
C                    Now solve the LU decomposed linear system
C                    to get the density derivatives
      CALL MLUSLV(DFDMLU,RESULT,DFDL_1,IPVT,5)
C
      DU_DN = RESULT(1)
      DX_DN = RESULT(2)
      DNI_DN = RESULT(3)
      DEP_DN = RESULT(4)
      DEN_DN = RESULT(5)
C
C
C                    Now solve the LU decomposed linear system
C                    to get the Temperature derivatives
      CALL MLUSLV(DFDMLU,RESULT,DFDL_2,IPVT,5)
C
      DU_DT = RESULT(1)
      DX_DT = RESULT(2)
      DNI_DT = RESULT(3)
      DEP_DT = RESULT(4)
      DEN_DT = RESULT(5)
C
C                    Now solve the LU decomposed linear system
C                    to get the Ye derivatives
      CALL MLUSLV(DFDMLU,RESULT,DFDL_3,IPVT,5)
C
      DU_DY = RESULT(1)
      DX_DY = RESULT(2)
      DNI_DY = RESULT(3)
      DEP_DY = RESULT(4)
      DEN_DY = RESULT(5)
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 ------------------------------------
C                 !      Derivatives of finite size  !
C                 !      terms in the internal       !
C                 !      energy and entropy          !
C                 !      densities w.r.t. to U,X,n_i !
C                 !      and T.  These are used in   !
C                 !      calculating the derivatives !
C                 !      w.r.t. the independant vars !
C                 !      (Baryon density, T, and Ye) !
C                 !                                  !
C                 ------------------------------------
C
C                        Free energy Surface & Coulomb terms
C                                  (Densities)
C
      F_SC = ZETA*SCRDU
C
      DFSCDU = ZETA*SCRDUP
C
      DFSCDX = ZETA*SCRDUX+SCRDU*DZDX
C
      DFSCDN = SCRDU*DZDNI
C
      DFSCDT = ZETA*SCRDUT+SCRDU*DZDT
C
C
C                        Free energy translational terms
C                                  (Densities)
      FTR = U_NUC*EXCLU*NSUBI*FTRANS
C
      DFTRDT = FTR*(HPRIM/H+1.0/T)-
     1    1.5*TRSCAL*U_NUC*EXCLU*NSUBI*H/AZERO
C
      DFTRDX = FTR*DHDX/H
C
      DFTRDU = FTR/U_NUC-FTR/EXCLU+
     1    TRSCAL*NSUBI*H*(1.0-2.0*U_NUC)/AZERO
C
      DFTRDN = FTR/NSUBI+TRSCAL*U_NUC*EXCLU*H*T/AZERO
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C                        Internal energy Surface & Coulomb terms
C                                  (Densities)
C
      TMP4 = 1.0-T*SCRDUT/SCRDU-OVR23*T*HPRIM/H
C
      E_SC = F_SC*TMP4
C
      DESCDU = DFSCDU*TMP4+
     1    F_SC*(T*SCRDUT*SCRDUP/SCRDU**2-T*SCRDPT/SCRDU)
C
      DESCDX = DFSCDX*TMP4+
     1    F_SC*(T*SCRDUT*SCRDUX/SCRDU**2-T*SCRDXT/SCRDU+
     2    OVR23*T*HPRIM*DHDX/H**2-OVR23*T*DHDTDX/H)
C
      DESCDN = DFSCDN*TMP4
C
      DESCDT = DFSCDT*TMP4+F_SC*
     1   (T*(SCRDUT**2)/SCRDU**2-SCRDUT/SCRDU-T*SCRDTT/SCRDU+
     2    OVR23*T*(HPRIM**2)/H**2-OVR23*HPRIM/H-OVR23*T*HPPRIM/H)
C
C                        Internal energy translational terms
C                                  (Densities)
C
      TMP4 = 1.5*H*T/AZERO-T*HPRIM*(MUSUBT-T)/AZERO
C
      E_TR = TRSCAL*EXCLU*BRYDNS*XH*TMP4
C
      DETRDU = TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4-
     1    NSUBI*(T**2)*HPRIM*(1.0-2.0*U_NUC)/AZERO)
C
      DETRDX = TRSCAL*BRYDNS*XH*EXCLU*
     1    (1.5*T*DHDX/AZERO-T*(MUSUBT-T)*DHDTDX/AZERO)
C
      DETRDN = TRSCAL*(U_NUC*EXCLU*TMP4-
     1    BRYDNS*XH*EXCLU*(T**2)*HPRIM/(NSUBI*AZERO))
C
      DETRDT = TRSCAL*BRYDNS*XH*EXCLU*
     1    (1.5*(H+T*HPRIM)/AZERO-(HPRIM+T*HPPRIM)*(MUSUBT-T)/AZERO-
     2    T*HPRIM*(MUSUBT/T-2.5)/AZERO )
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C                        Entropy Surface & Coulomb terms
C                                  (Densities)
C
      S_SC = (E_SC-F_SC)/T
C
      DSSCDU = (DESCDU-DFSCDU)/T
C
      DSSCDX = (DESCDX-DFSCDX)/T
C
      DSSCDN = (DESCDN-DFSCDN)/T
C
      DSSCDT = (DESCDT-DFSCDT)/T-(E_SC-F_SC)/T**2
C
C                        Entropy translational terms
C                                  (Densities)
C
      TMP4 = MUSUBT*(HPRIM+H/T)/AZERO-(T*HPRIM+2.5*H)/AZERO
C
      S_TR = -TRSCAL*BRYDNS*XH*EXCLU*TMP4
C
      DSTRDU = -TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4+
     1    NSUBI*T*(1.0-2.0*U_NUC)*(HPRIM+H/T)/AZERO)
C
      DSTRDX = -TRSCAL*BRYDNS*XH*EXCLU*
     1    (MUSUBT*(DHDTDX+DHDX/T)/AZERO-
     2    (T*DHDTDX+2.5*DHDX)/AZERO)
C
      DSTRDN = -TRSCAL*
     1    (U_NUC*EXCLU*TMP4+U_NUC*EXCLU*T*(HPRIM+H/T)/AZERO)
C
      DSTRDT = -(BRYDNS*XH*EXCLU*((MUSUBT/T-1.5)*(HPRIM+H/T)/AZERO+
     1    MUSUBT*(HPPRIM+HPRIM/T-H/T**2)/AZERO-
     2    (3.5*HPRIM+T*HPPRIM)/AZERO ))*TRSCAL
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of interior bulk !
C                 !      terms in the internal        !
C                 !      energy and entropy           !
C                 !      densities w.r.t. to U,X,n_i  !
C                 !      and T.  These are used in    !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
C
      S_NUC =(OVR53*UQ/T)*(TAU_NI+TAU_PI)-
     1    NSUBI*((1.0-X)*ETA_NI+X*ETA_PI)
C
c20      E_NUC = UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+
c20     1    CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM
      E_NUC = UQ*(TAU_PI+TAU_NI)+PV_E(X*NSUBI,(1.0-X)*NSUBI)
C
C
C                    Interior particle densties
      NPI = X*NSUBI
      NNI = (1.0-X)*NSUBI
C
      DTPIDT = 2.5*TAU_PI/T-2.25*NPI*GPI/UQ
      DTNIDT = 2.5*TAU_NI/T-2.25*NNI*GNI/UQ
C
C               Derivative of interior entropy density w.r.t. T
      DSIDT = UQ*(DTPIDT+DTNIDT)/T
C
C               Derivative of interior internal energy density w.r.t. T
      DEIDT = T*DSIDT
C
C
C
C
C                    Derivatives of eta's w.r.t. X and NSUBI
      DETPDX = GPI/X
      DETNDX = -GNI/(1.0-X)
      DETPDN = GPI/NSUBI
      DETNDN = GNI/NSUBI
C
C                    Derivatives of Tau's w.r.t. X and NSUBI
      DTPIDX = 1.5*T*NPI*DETPDX/UQ
      DTNIDX = 1.5*T*NNI*DETNDX/UQ
      DTPDNI = 1.5*T*NPI*DETPDN/UQ
      DTNDNI = 1.5*T*NNI*DETNDN/UQ
C
C
C
C           Derivative of interior entropy density w.r.t. X
      DSIDX = OVR53*UQ*(DTPIDX+DTNIDX)/T-NSUBI*(ETA_PI-ETA_NI)-
     1    NSUBI*((1.0-X)*DETNDX+X*DETPDX)
C
C           Derivative of interior internal energy density w.r.t. X
c20      DEIDX = UQ*(DTPIDX+DTNIDX)+
c20     1    (NSUBI**2)*4.0*BB*(1.0-2.0*X)+NSUBI*DELTAM
      DEIDX = UQ*(DTPIDX+DTNIDX)+DPVEDX(NSUBI,X)
C
C
C           Derivative of interior entropy density w.r.t. NSUBI
      DSIDN = OVR53*UQ*(DTPDNI+DTNDNI)/T-((1.0-X)*ETA_NI+X*ETA_PI)-
     1    NSUBI*((1.0-X)*DETNDN+X*DETPDN)
C
C
C           Derivative of interior internal energy density w.r.t. NSUBI
c20      DEIDN = UQ*(DTPDNI+DTNDNI)+2.0*NSUBI*(AA+4.0*BB*X*(1.0-X))+
c20     1    CC*(1.0+DD)*(NSUBI**DD)+X*DELTAM
      DEIDN = UQ*(DTPDNI+DTNDNI)+DPVEDN(NSUBI,X)
C
C
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of exterior bulk !
C                 !      nucleon internal energy &    !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO
C
c20      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*
c20     1((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*PV_E(NPOUT,NNOUT)
C
C                   Derivative of exterior entropy density w.r.t. T
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-
     1     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T
C
      DEODT = T*DSODT
C
C                    Derivatives of exterior particle densities w.r.t.
C                    Temperature (ETA's fixed)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C
      DNPDEP = NPOUT/GPO
      DNNDEN = NNOUT/GNO
C
      DTPDEP = 1.5*T*NPOUT/UQ
      DTNDEN = 1.5*T*NNOUT/UQ
C
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN
C
C
C                    Exterior particle potentials
c20      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD )
      VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
c20      VPOUT = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)
      VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN
C
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of alpha         !
C                 !      particle internal energy &   !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
C
C
      S_ALFA = ALFDNS*(2.5-MUALFA/T)
C
C
      E_ALFA = ALFDNS*(1.5*T-BALPHA)
C
C                  Derivative of pressure potential w.r.t. T
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2
C
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS
C
C
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN
C
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN
C
      DNADEP = ALFDNS*DMADEP/T
      DNADEN = ALFDNS*DMADEN/T
C
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T
C
      DEADEP = DNADEP*(1.5*T-BALPHA)
      DEADEN = DNADEN*(1.5*T-BALPHA)
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      S_DENS = U_NUC*S_NUC+EXCLU*EXALFA*S_OUT+EXCLU*S_ALFA+S_SC+S_TR
C
      E_DENS = U_NUC*E_NUC+EXCLU*EXALFA*E_OUT+EXCLU*E_ALFA+E_SC+E_TR
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !      Temperature Derivatives     !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT
C
C
      DBSDT = (DU_DT*S_NUC-
     1    DU_DT*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DT*S_OUT
     2    -DU_DT*S_ALFA+
     3    U_NUC*(DSIDT+DSIDX*DX_DT+DSIDN*DNI_DT)+
     4    EXCLU*EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+
     5    EXCLU*(DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT)+
     6    DSSCDT+DSSCDU*DU_DT+DSSCDX*DX_DT+DSSCDN*DNI_DT+
     7    DSTRDT+DSTRDU*DU_DT+DSTRDX*DX_DT+DSTRDN*DNI_DT)/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDT = T*DBSDT
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+
     1    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDT = BRYDNS*(DBMUDT-DBFDT)
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !       Density Derivatives        !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN
C
C
      DBSDN = (DU_DN*S_NUC-
     1    DU_DN*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DN*S_OUT-DU_DN*S_ALFA+
     2    U_NUC*(DSIDX*DX_DN+DSIDN*DNI_DN)+
     3    EXCLU*EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+
     4    EXCLU*(DSADEP*DEP_DN+DSADEN*DEN_DN)+
     5    DSSCDU*DU_DN+DSSCDX*DX_DN+DSSCDN*DNI_DN+
     6    DSTRDU*DU_DN+DSTRDX*DX_DN+DSTRDN*DNI_DN)/BRYDNS-
     7    S_DENS/BRYDNS**2
C
C
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDN = (DU_DN*E_NUC-
     1    DU_DN*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DN*E_OUT-DU_DN*E_ALFA+
     2    U_NUC*(DEIDX*DX_DN+DEIDN*DNI_DN)+
     3    EXCLU*EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+
     4    EXCLU*(DEADEP*DEP_DN+DEADEN*DEN_DN)+
     5    DESCDU*DU_DN+DESCDX*DX_DN+DESCDN*DNI_DN+
     6    DETRDU*DU_DN+DETRDX*DX_DN+DETRDN*DNI_DN)/BRYDNS-
     7    E_DENS/BRYDNS**2
C
C
C
C
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDN = DBUDN-T*DBSDN
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+
     1    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !         Ye Derivatives           !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
C
C
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY
C
C
      DBSDY = (DU_DY*S_NUC-
     1    DU_DY*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DY*S_OUT-DU_DY*S_ALFA+
     2    U_NUC*(DSIDX*DX_DY+DSIDN*DNI_DY)+
     3    EXCLU*EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+
     4    EXCLU*(DSADEP*DEP_DY+DSADEN*DEN_DY)+
     5    DSSCDU*DU_DY+DSSCDX*DX_DY+DSSCDN*DNI_DY+
     6    DSTRDU*DU_DY+DSTRDX*DX_DY+DSTRDN*DNI_DY)/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDY = (DU_DY*E_NUC-
     1    DU_DY*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DY*E_OUT-DU_DY*E_ALFA+
     2    U_NUC*(DEIDX*DX_DY+DEIDN*DNI_DY)+
     3    EXCLU*EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+
     4    EXCLU*(DEADEP*DEP_DY+DEADEN*DEN_DY)+
     5    DESCDU*DU_DY+DESCDX*DX_DY+DESCDN*DNI_DY+
     6    DETRDU*DU_DY+DETRDX*DX_DY+DETRDN*DNI_DY)/BRYDNS
C
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDY = DBUDY-T*DBSDY
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+
     1    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDY = BRYDNS*(DBMUDY-DBFDY)
C
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C-----------------------------------------------------------------------
C                End of derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C                  Total derivatives
C                  (Baryons+Electrons+Photons)
C
      DUDT = DBUDT+DEUDT+DPUDT
      DUDN = DBUDN+DEUDN+DPUDN
      DUDY = DBUDY+DEUDY+DPUDY
C
C
      DSDT = DBSDT+DESDT+DPSDT
      DSDN = DBSDN+DESDN+DPSDN
      DSDY = DBSDY+DESDY+DPSDY
C
C
      DPDT = DBPDT+DEPDT+DPPDT
      DPDN = DBPDN+DEPDN+DPPDN
      DPDY = DBPDY+DEPDY+DPPDY
C
C
      DMUDT = DBMUDT+YE*DEMUDT
      DMUDN = DBMUDN+YE*DEMUDN
      DMUDY = DBMUDY+YE*DEMUDY
C
C                Calculate the adiabatic index
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT)
C
C
C                Set the value of XPREV to X for use the next 
C                time through
C
      XPREV = X
C
C                Save the value of the proton density to be used
C                by the "no nuclei" scheme on the next call
      P_PREV = NPOUT
C
C
C                Return the three internal compositional variables
      INPVAR(2) = NSUBI
      INPVAR(3) = ETA_PO
      INPVAR(4) = ETA_NO
C
C
C
C                Rejoice for this routine is finished!!!!!!!
 999  RETURN
C
C
      END
Cnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnu
Calfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalf
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         ALFEOS.FOR
C
C***********************************************************************
C
C    MODULE:       ALFEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         8/30/90 Modified from model 4-A
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU
C                            FSWESTY@SBAST3.SUNYSB.EDU
C
C
C    CALL LINE:    CALL ALFEOS(INPVAR,YE,BRYDNS)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:
C
C 
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG)
C
      IMPLICIT NONE
C
C
      INCLUDE 'eos_m4c.inc'
      INCLUDE 'el_eos.inc'
C
C                       "ZERO" flag
      INTEGER ZFLAG
C
C                       "Negative" flag
      INTEGER NFLAG
C
C                       Function type declarations
C
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFO, FHALF
C
C                       Include the nucleon-nucleon interaction
C                       statement function definitions
      INCLUDE 'force.inc'
C
C
C
C
C
C                         Unset the zero flag
      ZFLAG = 0
C
C                         Ratio of baryon density to saturation density
      Y = BRYDNS/NSUBS
C
C
C                         Set T equal to the input variable (the entropy
C                         and internal energy options are not implemented
C                         in this version)
      T = INPVAR(1)
C
C
C                         Calc the quantum concentration of nucleons
      NQ = 2.36D-4*T**1.5 
C
C                         Calc the Fermi integral coefficent
      UQ = 20.721
C
      MQ = (T/UQ)**1.5
C
      KQ = ((T/UQ)**2.5)/(2.0*PI**2)
C
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2))
C
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ)
C
C
C
C                              Set the proton density to its old value
      NPOUT = P_PREV
C
      IF(BRYDNS.GT.(0.98*2.0/(YE*V_ALFA))) THEN
        NPOUT = YE*BRYDNS
        NNOUT = (1.0-YE)*BRYDNS
        NOUT = BRYDNS
C
C
c20        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
        VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1    CC*(1.0+DD)*NOUT**DD+DELTAM)
        VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
        ZNO = 2.0*(PI**2)*NNOUT/MQ
C
        ZPO = 2.0*(PI**2)*NPOUT/MQ
C
        ETA_NO = FINV12(ZNO)
C
        ETA_PO = FINV12(ZPO)
C
        F32_NO = F_3_2(ETA_NO)
        F32_PO = F_3_2(ETA_PO)
C
        TAU_NO = KQ*F32_NO
        TAU_PO = KQ*F32_PO
C
C
c20        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
C
        MUN_O = T*ETA_NO+VNOUT
        MUN = MUN_O
C
        MUP_O = T*ETA_PO+VPOUT
        MUPROT = MUP_O
C
C                              Calculate diff. of chem. potentials
        MUHAT = MUN-MUPROT
C
C
C                              Calculate the alpha particle
C                              chemical potential
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA
C
        ALFDNS = 0.0
C
        EXALFA = 1.0-ALFDNS*V_ALFA
C
      ELSE
C
C                              Calculate the neutron density
        NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+
     1            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/
     2            (2.0-BRYDNS*YE*V_ALFA)
C
C                              Calculate density of outside nucleons
        NOUT = NPOUT+NNOUT
C
c20        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
        VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1    CC*(1.0+DD)*NOUT**DD+DELTAM)
        VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
        ZNO = 2.0*(PI**2)*NNOUT/MQ
C
        ZPO = 2.0*(PI**2)*NPOUT/MQ
C
        ETA_NO = FINV12(ZNO)
C
        ETA_PO = FINV12(ZPO)
C
        F32_NO = F_3_2(ETA_NO)
        F32_PO = F_3_2(ETA_PO)
C
        TAU_NO = KQ*F32_NO
        TAU_PO = KQ*F32_PO
C
C
c20        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
C
        MUN_O = T*ETA_NO+VNOUT
        MUN = MUN_O
C
        MUP_O = T*ETA_PO+VPOUT
        MUPROT = MUP_O
C
C                              Calculate diff. of chem. potentials
        MUHAT = MUN-MUPROT
C
C
C                              Calculate the alpha particle
C                              chemical potential
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA
C
C                              Calculate density of alpha particles
C
        IF(ABS(MUALFA/T).LT.30.0) THEN
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T)
        ELSEIF((MUALFA/T).LT.-30.0) THEN
          ALFDNS = 0.0
        ELSE
          ALFDNS = 8.0*NQ*DEXP(3.0D1)
        ENDIF
C
C
        EXALFA = 1.0-ALFDNS*V_ALFA
C
C                              Calculate "non-zeroness" of baryon
C                              conservation equation and save the
C                              value to be used in the finite
C                              difference approximation of DGDPRT
        GOLD = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS
        PRTOLD = NPOUT
C
C                              Take a small step to get derivative
        NPOUT = NPOUT+0.001*BRYDNS
C
        DO 11 I=1,30,1
C
C                              Unset the negative flag
          NFLAG = 0
C
 14       CONTINUE
C                              Calculate the neutron density
          NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+
     1            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/
     2            (2.0-BRYDNS*YE*V_ALFA)
C
          IF((NNOUT.LT.0.0).AND.(I.EQ.1)) THEN
            NPOUT = PRTOLD-0.5*DPRT
          ELSEIF((NNOUT.LT.0.0).AND.(I.EQ.1).AND.(NFLAG.NE.1)) THEN
            NPOUT = 0.99*P_PREV
            NFLAG = 1
          ELSEIF((NNOUT.LT.0.0).AND.(I.EQ.1).AND.(NFLAG.EQ.1)) THEN
            SSFLAG = 0
            GOTO 999
          ENDIF
C                              Calculate density of outside nucleons
          NOUT = NPOUT+NNOUT
C
c20          VNOUT = EIFLAG*
c20     1      (2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
          VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20          VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1      CC*(1.0+DD)*NOUT**DD+DELTAM)
          VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
          ZNO = 2.0*(PI**2)*NNOUT/MQ
C
          ZPO = 2.0*(PI**2)*NPOUT/MQ
C
          ETA_NO = FINV12(ZNO)
C
          ETA_PO = FINV12(ZPO)
C
          F32_NO = F_3_2(ETA_NO)
C
          F32_PO = F_3_2(ETA_PO)
C
          TAU_NO = KQ*F32_NO
          TAU_PO = KQ*F32_PO
C
c20          BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*
c20     1      (AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
          BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
          MUN_O = T*ETA_NO+VNOUT
          MUN = MUN_O
C
          MUP_O = T*ETA_PO+VPOUT
          MUPROT = MUP_O
C
C                              Calc difference of potentials
          MUHAT = MUN-MUPROT
C
C                              Calc alpha particle chemical potentials
          MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA
C
C                              Calc alpha particle density
C
          IF(ABS(MUALFA/T).LT.30.0) THEN
            ALFDNS = 8.0*NQ*DEXP(MUALFA/T)
          ELSEIF((MUALFA/T).LT.-30.0) THEN
            ALFDNS = 0.0
          ELSE
            ALFDNS = 8.0*NQ*DEXP(3.0D1)
          ENDIF
C
C
          EXALFA = 1.0-ALFDNS*V_ALFA
C
C                              Calc "non-zeroness" of baryon cons. eqn.
          G = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS
C
C                              Calculate derivative of baryon conservation
C                              equation w.r.t. proton density by finite
C                              diference approximation
          DGDPRT = (G-GOLD)/(NPOUT-PRTOLD)
C
C                              If rate of change is near zero
C                              and zero flag is not set
          IF((ABS(DGDPRT).LT.1.0D-25).AND.(ZFLAG.EQ.0)) THEN
C                              Tweak the step size
            NPOUT = PRTOLD-0.5*DPRT
C                              and set the zero flag
            ZFLAG = 1
C                              and go back and try again
            GOTO 14
C                              If failure occurs again
          ELSEIF((ABS(DGDPRT).LT.1.0D-25).AND.(ZFLAG.EQ.1)) THEN
C                              declare an EOS failure
            SSFLAG = 0
C                              and return
            GOTO 999
          ENDIF
C
C                              Calculate new Newton-Raphson step
          DPRT = G/DGDPRT
C
C                              Save old value of proton density & G
          PRTOLD = NPOUT
          GOLD = G
C
C
 13       CONTINUE
C
C                              Potential "new" value of proton density
          PRTNEW = NPOUT-DPRT
C
C                              If new proton density is less than the
C                              baryon density and greater than zero 
C                              then update the proton density
          IF(PRTNEW*(BRYDNS-PRTNEW).GT.0.0) THEN
            NPOUT = NPOUT-DPRT
C                              Else cut the step size in half and try again
          ELSE
            DPRT = DPRT*0.5
            GOTO 13
          ENDIF
C
C                              If step size is small enough break out of
C                              the DO 11 loop, otherwise continue
          IF(ABS(DPRT/NPOUT).LT.10E-11) GOTO 12
 11     CONTINUE
C
c      write(*,*) 'A failed to converge; switching to F' ! take out later
        SSFLAG = 0
        GOTO 999
C
C
 12     CONTINUE
C
      ENDIF
C                              Set the success flag
      SSFLAG = 1
C
C
C                              Calc outside nucleon density
      NOUT = NNOUT+NPOUT
C
C                              Calc outside nucleon fraction
      XOUT = NPOUT/NOUT
C
C                              Calculate particle fractions
      XALFA = 4.0*ALFDNS/BRYDNS
      XPROT = EXALFA*NPOUT/BRYDNS
      XNUT = EXALFA*NNOUT/BRYDNS
      XH = 0.0
C
C                              Baryons
C
      F32_NO = F_3_2(ETA_NO)
C
      F32_PO = F_3_2(ETA_PO)
C
      TAU_PO = KQ*F32_PO
C
      TAU_NO = KQ*F32_NO
C
C
C
C
C
C
C                    Calculate internal energy of outside nucleons
c20      BUOUT = (XNUT+XPROT)*( UQ*(TAU_PO+TAU_NO)+
c20     1    EIFLAG*((NOUT**2)*AA+
c20     2   4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )/NOUT
      BUOUT = (XNUT+XPROT)*( UQ*(TAU_PO+TAU_NO)+
     1    EIFLAG*PV_E(NPOUT,NNOUT) )/NOUT
C
C
C                                Calc alfa particle internal energy
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA)
C
C                                Set nuclei internal energy to zero
      BUNUC = 0.0
C                                Calculate total baryon internal energy
C                                (per baryon)
      BU = BUOUT+BUALFA+BUNUC
C
C
C                                Calc entropy of outside nucleons
      BSOUT = (XNUT+XPROT)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)-
     1   NNOUT*ETA_NO-NPOUT*ETA_PO )/NOUT
C
C                                Calc alpha particle entropy
      BSALFA = 0.25*XALFA*(2.5-MUALFA/T)
C
C                                Set nuclei entropy to zero
      BSNUC = 0.0
C
C                                Calc total baryon entropy (per baryon)
      BS = BSOUT+BSALFA+BSNUC
C
C
C
C                                Calc outside free energy
      BFOUT = BUOUT-T*BSOUT
C                                Calc alpha particle free energy
      BFALFA = BUALFA-T*BSALFA
C                                Set nuclei free energy to zero
      BFNUC = BUNUC-T*BSNUC
C                                Calc total baryon free energy (per nucleon)
      BFTOT = BFOUT+BFALFA+BFNUC
C
C
C
C
C
C                                Calc outside pressure
c20      BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)))
      BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
C                                Calc alfa particle pressure
      BPRALF = ALFDNS*T
C
C                                Set nuclei pressure to zero
      BPRNUC = 0.0
C
C                                Calc total baryon pressure
      BPRESS = BPROUT+BPRALF+BPRNUC
C
C
C
C
C
C
C
C
C                           Leptons & Photons
      CALL EL_EOS(T,YE,BRYDNS)
C
C
C
C                           Total pressure and eng/ent per baryon
C
      FBARY = BFTOT+FSUBE
      PBARY = BPRESS+EPRESS
      MUBARY = YE*MUPROT+(1.0-YE)*MUN
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN
C
      FTOT = BFTOT+FSUBE+PF
      UTOT = BU+EU+PU
      STOT = BS+ES+PS
      PTOT = BPRESS+EPRESS+PPRESS
C
C
C
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------
C                Derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C
C
cc      GPO = 2.0*FHALFO(ETA_PO)
cc      GNO = 2.0*FHALFO(ETA_NO)
C
C
      GPO = 2.0*FHALF(ETA_PO)
      GNO = 2.0*FHALF(ETA_NO)
C
C
C                 ------------------------------------
C                 !      Derivatives of exterior     !
C                 !      quantities                  !
C                 !      (w.r.t. Temp. and ETA's)    !
C                 !                                  !
C                 ------------------------------------
C
C
C                  Derivatives of exterior potentials
C                  w.r.t. particle densities
c20      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
      DVPODP = EIFLAG*DPVPDP(NPOUT,NNOUT)
c20      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
      DVPODN = EIFLAG*DPVPDN(NPOUT,NNOUT)
c20      DVNODP = DVPODN
      DVNODP = EIFLAG*DPVNDP(NPOUT,NNOUT)
c20      DVNODN = DVPODP
      DVNODN = EIFLAG*DPVNDN(NPOUT,NNOUT)
C
C
C                  Derviatives of exterior chem. pot. w.r.t. ETA's
C                  (at fixed T)
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C                  Derivatives of pressure potential w.r.t.
C                  particle densities
c20      DV_DPO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DPO = EIFLAG*DPVRDP(NPOUT,NNOUT)
c20      DV_DNO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DNO = EIFLAG*DPVRDN(NPOUT,NNOUT)
C
C                  Derivatives of pressure potential w.r.t. ETA's
C                  (at fixed T)
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
C                  Derivatives of outside pressure w.r.t. ETA's
C                  (at fixed T)
      DPODEP = NPOUT*T+DV_DEP
      DPODEN = NNOUT*T+DV_DEN
C
C                  Derivatives of alpha density w.r.t. ETA's
C                  (at fixed T)
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T
C
C
C                  Derivatives of particle densities w.r.t. T
C                  (at fixed ETA's)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
C                  Derivatives of exterior chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C                  Derivative of pressure potential w.r.t. T
C                  (at fixed ETA's)
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
C                  (at fixed ETA's)
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C                  Derivative of alpha chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
C                  (at fixed ETA's)
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C                  Derivative of alpha particle pressure w.r.t. T
C                  (at fixed ETA's)
      DPADT = ALFDNS+T*DNADT
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C                 ------------------------------------
C                 !      Derivatives of constraint   !
C                 !      and equilibrium equations   !
C                 !      with respect to the five    !
C                 !      compositional variables     !
C                 !      (U,x,n_i,eta_po,eta_no)     !
C                 !      and the three independent   !
C                 !      variables                   !
C                 !      (Baryon density, T, and Ye) !
C                 !                                  !
C                 ------------------------------------
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 1 (Baryon conservation)
C
C
C
      DG1DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NOUT-4.0)*DNADEP
C
      DG1DO2 = -EXALFA*NNOUT/GNO+(V_ALFA*NOUT-4.0)*DNADEN
C
C
      DG1DL1 = 1.0
C
      DG1DL2 = -EXALFA*(DNNODT+DNPODT)+(V_ALFA*NOUT-4.0)*DNADT
C
      DG1DL3 = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 2 (Charge conservation)
C
C
      DG2DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NPOUT-2.0)*DNADEP
C
      DG2DO2 = (V_ALFA*NPOUT-2.0)*DNADEN
C
C
      DG2DL1 = YE
C
      DG2DL2 = -EXALFA*DNPODT+(V_ALFA*NPOUT-2.0)*DNADT
C            
      DG2DL3 = BRYDNS
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      DET_GT = DG1DO1*DG2DO2-DG1DO2*DG2DO1
C
C
      DEP_DN = (DG1DO2*DG2DL1-DG2DO2*DG1DL1)/DET_GT
      DEN_DN = (DG2DO1*DG1DL1-DG1DO1*DG2DL1)/DET_GT
C
C
      DEP_DT = (DG1DO2*DG2DL2-DG2DO2*DG1DL2)/DET_GT
      DEN_DT = (DG2DO1*DG1DL2-DG1DO1*DG2DL2)/DET_GT
C
C
      DEP_DY = (DG1DO2*DG2DL3-DG2DO2*DG1DL3)/DET_GT
      DEN_DY = (DG2DO1*DG1DL3-DG1DO1*DG2DL3)/DET_GT
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of exterior bulk !
C                 !      nucleon internal energy &    !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO
C
c20      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*
c20     1((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*PV_E(NPOUT,NNOUT)
C
C                   Derivative of exterior entropy density w.r.t. T
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-
     1     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T
C
      DEODT = T*DSODT
C
C                    Derivatives of exterior particle densities w.r.t.
C                    Temperature (ETA's fixed)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C
      DNPDEP = NPOUT/GPO
      DNNDEN = NNOUT/GNO
C
      DTPDEP = 1.5*T*NPOUT/UQ
      DTNDEN = 1.5*T*NNOUT/UQ
C
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN
C
C
C                    Exterior particle potentials
c20      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD )
      VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
c20      VPOUT = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)
      VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN
C
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of alpha         !
C                 !      particle internal energy &   !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
C
C
      S_ALFA = ALFDNS*(2.5-MUALFA/T)
C
C
      E_ALFA = ALFDNS*(1.5*T-BALPHA)
C
C                  Derivative of pressure potential w.r.t. T
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2
C
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS
C
C
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN
C
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN
C
      DNADEP = ALFDNS*DMADEP/T
      DNADEN = ALFDNS*DMADEN/T
C
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T
C
      DEADEP = DNADEP*(1.5*T-BALPHA)
      DEADEN = DNADEN*(1.5*T-BALPHA)
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      S_DENS = EXALFA*S_OUT+S_ALFA
C
      E_DENS = EXALFA*E_OUT+E_ALFA
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !      Temperature Derivatives     !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT
C
C
      DBSDT = (-V_ALFA*DNA_DT*S_OUT+
     1    EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+
     2    (DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT) )/BRYDNS
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDT = T*DBSDT
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+
     1    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDT = BRYDNS*(DBMUDT-DBFDT)
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !       Density Derivatives        !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN
C
C
      DBSDN = (-V_ALFA*DNA_DN*S_OUT+
     1    EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+
     2   (DSADEP*DEP_DN+DSADEN*DEN_DN) )/BRYDNS-S_DENS/BRYDNS**2
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBUDN = (-V_ALFA*DNA_DN*E_OUT+
     1    EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+
     2   (DEADEP*DEP_DN+DEADEN*DEN_DN) )/BRYDNS-E_DENS/BRYDNS**2
C
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDN = DBUDN-T*DBSDN
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+
     1    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !         Ye Derivatives           !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
C
C
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY
C
C
      DBSDY = (-V_ALFA*DNA_DY*S_OUT+
     1    EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+
     2   (DSADEP*DEP_DY+DSADEN*DEN_DY) )/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBUDY = (-V_ALFA*DNA_DY*E_OUT+
     1    EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+
     2   (DEADEP*DEP_DY+DEADEN*DEN_DY) )/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDY = DBUDY-T*DBSDY
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+
     1    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDY = BRYDNS*(DBMUDY-DBFDY)
C
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C-----------------------------------------------------------------------
C                End of derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C                  Total derivatives
C                  (Baryons+Electrons+Photons)
C
      DUDT = DBUDT+DEUDT+DPUDT
      DUDN = DBUDN+DEUDN+DPUDN
      DUDY = DBUDY+DEUDY+DPUDY
C
C
      DSDT = DBSDT+DESDT+DPSDT
      DSDN = DBSDN+DESDN+DPSDN
      DSDY = DBSDY+DESDY+DPSDY
C
C
      DPDT = DBPDT+DEPDT+DPPDT
      DPDN = DBPDN+DEPDN+DPPDN
      DPDY = DBPDY+DEPDY+DPPDY
C
C
      DMUDT = DBMUDT+YE*DEMUDT
      DMUDN = DBMUDN+YE*DEMUDN
      DMUDY = DBMUDY+YE*DEMUDY
C
C
C                  Adiabatic index
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT)
C
C
      INPVAR(2) = NSUBS
      INPVAR(3) = ETA_PO
      INPVAR(4) = ETA_NO
C
C
C                           Approximate the nuclear density
      NSUBI = NSUBS
C
C                           Use 0.45 as the nuclear proton fraction
      X = 0.45
      A = 4.0D0
C
C                           Save the proton number density for use
C                           as the initial guess on next call
      P_PREV = NPOUT
C
C
 999  RETURN
C
      END
Calfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalf
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         MAXWEL.FOR
C
C***********************************************************************
C
C    MODULE:       MAXWEL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/13/90
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU
C
C
C    CALL LINE:    CALL MAXWEL(INPVAR,YE,BRYDNS)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:
C
C
C
C 
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION OUTVAR(4)
      DOUBLE PRECISION DPN_DT, DPN_DN, DPN_DY
      DOUBLE PRECISION DSN_DT, DSN_DN, DSN_DY
      DOUBLE PRECISION DSB_DT, DSB_DN, DSB_DY
      DOUBLE PRECISION DMU_DT, DMU_DN, DMU_DY
      DOUBLE PRECISION DPHADT, DPHADY, DELDNS
      DOUBLE PRECISION N_XH, N_XA, N_XN, N_XP, B_XA, B_XN, B_XP
C
C
      INCLUDE 'eos_m4c.inc'
      INCLUDE 'el_eos.inc'
      INCLUDE 'maxwel.inc'
C
C
C                   Set the temperature
      T = INPVAR(1)
C
C
C                   Calculate and save chem. pot. and thermodynamic
C                   quantaties from low end of two phase region
      CALL NUCEOS(INPVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG)
C
C
C
C                    If the nuclear EOS failed and the reset flag is set
C                    then reset the initial guesses and try again
      IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN
        CALL RESET(INPVAR,YE,LOWDNS,OUTVAR)
        OUTVAR(1) = INPVAR(1)
        CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    Make a last ditch effort at convergence
        IF(SSFLAG.NE.1) THEN
          OUTVAR(2) = 0.155
          OUTVAR(3) = -15.0
          OUTVAR(4) = -20.0
          CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG)
        ELSE
          INPVAR(2) = OUTVAR(2)
          INPVAR(3) = OUTVAR(3)
          INPVAR(4) = OUTVAR(4)
        ENDIF
C
      ENDIF
C
C
C
C
      PRLOW = PTOT-PPRESS
      S_LOW = STOT-PS
      F_LOW = FTOT-PF
      MUTLOW = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE)
      MUELOW = MUSUBE
      MUHLOW = MUHAT
C
      DPN_DT = DPDT
      DPN_DN = DPDN
      DPN_DY = DPDY
C
      DMU_DT = DMUDT
      DMU_DN = DMUDN
      DMU_DY = DMUDY
C
      DSN_DT = DSDT-DPSDT
      DSN_DN = DSDN
      DSN_DY = DSDY
C
      N_XH = XH
      N_XA = XALFA
      N_XP = XPROT
      N_XN = XNUT
C
C
      IF(SSFLAG.NE.1) THEN
        WRITE(*,*) 'MAXWEL:  Nuclear EOS failed at try:'
        WRITE(*,*) T,LOWDNS,YE
        WRITE(*,*) INPVAR
        GOTO 999
      ENDIF
C                   Calculate and save chem. pot. and thermodynamic
C                   quantaties from high end of two phase region
      CALL ALFEOS(INPVAR,YE,HIDNS,P_PREV,SSFLAG)
C
      PRHI = PTOT-PPRESS
      S_HI = STOT-PS
      F_HI = FTOT-PF
      MUTHI = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE)
      MUEHI = MUSUBE
      MUHHI = MUHAT
C
C
      DSB_DT = DSDT-DPSDT
      DSB_DN = DSDN
      DSB_DY = DSDY
C
C
      B_XA = XALFA
      B_XP = XPROT
      B_XN = XNUT
C
C
      IF(SSFLAG.NE.1) THEN
        WRITE(*,*) 'MAXWEL:  Alfa EOS failed at try:'
        WRITE(*,*) T,HIDNS,YE
        WRITE(*,*) INPVAR
        GOTO 999
      ENDIF
C
C                   Calculate "average" chem. pot. and pressure
C                   in order to avoid numerical problems
      MUTILD = (MUTLOW+MUTHI)/2.0
      PRTILD = (PRLOW+PRHI)/2.0
C
C                   Calculate phase fraction
      PHASEF = (BRYDNS-LOWDNS)/(HIDNS-LOWDNS)
C
C
C                   Electron number density
      NSUBE = BRYDNS*YE
C
C                   Call electron EOS to determine the
C                   electron chemical potential
      CALL EL_EOS(T,YE,BRYDNS)
C
C
      MUHAT = MUSUBE+(1.0-PHASEF)*(MUHLOW-MUELOW)+PHASEF*(MUHHI-MUEHI)
C
      MUN = MUTILD+YE*(MUHAT-MUSUBE)
C
      MUPROT = MUN-MUHAT
C
C                   Calculate thermodynamic quantities
C
      STOT = ((1.0-PHASEF)*S_LOW*LOWDNS+PHASEF*S_HI*HIDNS)/BRYDNS+PS
C
      FTOT = (LOWDNS*F_LOW+MUTILD*(BRYDNS-LOWDNS))/BRYDNS+PF
C
      UTOT = FTOT+T*STOT+PU
C
      PTOT = PRTILD+PPRESS
C
C
      XH = (1.0-PHASEF)*N_XH
      XALFA = (1.0-PHASEF)*N_XA
      XNUT = (1.0-PHASEF)*N_XN
      XPROT = (1.0-PHASEF)*N_XP
      XALFA2 = PHASEF*B_XA
      XNUT2 = PHASEF*B_XN
      XPROT2 = PHASEF*B_XP
C
C
C
C
      DELDNS = HIDNS-LOWDNS
C
C
      DPHADT = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DT-
     1    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DT
C
      DPDT = DPN_DT+DPN_DN*DNL_DT
      DMUDT = DMU_DT+DMU_DN*DNL_DT
      DSDT = (1.0-PHASEF)*LOWDNS*(DSN_DT+DSN_DN*DNL_DT)/BRYDNS+
     2 (1.0-PHASEF)*S_LOW*DNL_DT/BRYDNS-LOWDNS*S_LOW*DPHADT/BRYDNS+
     3    (DPHADT*S_HI*HIDNS+PHASEF*DNH_DT*S_HI+
     4    PHASEF*HIDNS*(DSB_DT+DSB_DN*DNH_DT))/BRYDNS+DPSDT
      DUDT = DMUDT-DPDT/BRYDNS+STOT+T*DSDT
C 
C
      DPDN = 0.0
      DMUDN = 0.0
      DSDN = -DPDT/BRYDNS**2
      DUDN = (LOWDNS*(MUTILD-FTOT)/BRYDNS**2)+T*DSDN
C
C
      DPHADY = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DY-
     1    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DY
C
      DPDY = DPN_DY+DPN_DN*DNL_DY
      DMUDY = DMU_DY+DMU_DN*DNL_DY
      DSDY = (1.0-PHASEF)*LOWDNS*(DSN_DY+DSN_DN*DNL_DY)/BRYDNS+
     2 (1.0-PHASEF)*S_LOW*DNL_DY/BRYDNS-LOWDNS*S_LOW*DPHADY/BRYDNS+
     3    (DPHADY*S_HI*HIDNS+PHASEF*DNH_DY*S_HI+
     4    PHASEF*HIDNS*(DSB_DY+DSB_DN*DNH_DY))/BRYDNS
      DUDY = DMUDY-DPDY/BRYDNS+T*DSDY
C
C
C
C
C             Adiabatic index
C             (Note that the first term vanishes in this expression)
      GAM_S = T*(DPDT**2)/(BRYDNS*PTOT*DUDT)
C
C
 999  RETURN
C
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EOSLOG.FOR
C
C***********************************************************************
C
C    MODULE:       EOSLOG
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         12/15/90
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C
C    CALL LINE:    CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)
C
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C
C
C    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B
C                           2 --> GENERAL EOS
C                           3 --> BULK EOS (includes alpha's)
C
C
C
C 
C    INCLUDE FILES:  EOS_M4C.INC, MAXWEL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)
C
C
      IMPLICIT NONE
C
C                       This include file contains all variable
C                       declarations.  NOTE:: no implicit typing
C                       scheme is followed in this code; if you
C                       have any doubt as to a variables type CHECK
C                       IT!!!!.  Also note that ALL variables are
C                       declared explicitly.
C
      INCLUDE 'eos_m4c.inc'
      INCLUDE 'maxwel.inc'
C
      DOUBLE PRECISION NLOW, NHI, N_CUT, TEMP_1, TEMP_2, T_BNDY
C
      DOUBLE PRECISION LMM, LMP, LPM, LPP
      DOUBLE PRECISION DNDY1, DNDY2
C
C
C
 10   CONTINUE
C
C
C                         Set T equal to the input variable (any calls
C                         with entropy or internal energy should go
C                         the the EOS_M4B subroutine)
C
      T = INPVAR(1)
C
C-----------------------------------------------------------------------
C         code to figure out the boundaries from the tables
C-----------------------------------------------------------------------
C
C
C
C
      IF(YE.GT.Y_HI) THEN
C                         Ye is too large for EOS
C
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time'
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_HI,' instead'
        YE = Y_HI-1.0D-6
        GOTO 10
C
      ELSEIF(YE.GE.Y_LOW) THEN
C                         Calculate high and low boundary densities
C                         for the Maxwell construction
C
C----------------------------------------------------------
C           Calc Ye index
C----------------------------------------------------------
C
        YFRAC = (YE-Y_LOW)/(Y_HI-Y_LOW)
        J_MXWL = INT(YFRAC*(NUMYE-1))+1
        DELT_Y = (Y_HI-Y_LOW)/DBLE(NUMYE-1)
C
        YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
        YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
C
C
        IF((YE.GE.YMINUS).AND.(YE.LE.YPLUS)) THEN
          J_BD = J_MXWL
          J_BNDY = J_MXWL
        ELSEIF(YE.GT.YPLUS) THEN
          J_MXWL = J_MXWL+1
          J_BD = J_MXWL
          J_BNDY = J_MXWL
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
        ELSE
          J_MXWL = J_MXWL-1
          J_BD = J_MXWL
          J_BNDY = J_MXWL
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
        ENDIF
C
C
        IF(J_MXWL.GT.(NUMYE-1)) THEN
          J_MXWL = NUMYE-1
          J_BD = J_MXWL
          J_BNDY = J_MXWL
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
        ENDIF
C
C
        YINTRP = (YE-YMINUS)/(YPLUS-YMINUS)
C
C----------------------------------------------------------
C           Calc T index
C----------------------------------------------------------
C
C
        TFRAC = (T-T_LOW)/(T_HI-T_LOW)
        I_MXWL = INT(TFRAC*(NUMTMP-1))+1
        DELT_T = (T_HI-T_LOW)/DBLE(NUMTMP-1)
C
        TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
        TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
C
C
        IF((T.GT.TMINUS).AND.(T.LE.TPLUS)) THEN
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ELSEIF(T.GT.TPLUS) THEN   
          I_MXWL = I_MXWL+1
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ELSE
          I_MXWL = I_MXWL-1
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ENDIF
C
C
        IF(I_MXWL.GT.(NUMTMP-1)) THEN
          I_MXWL = NUMTMP-1
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ENDIF
C
C
        TINTRP = (T-TMINUS)/(TPLUS-TMINUS)
C
C
C
C
C                Find the temperature and density at the top of the
C                Maxwel construction
C
CC      T_MXWL = YINTRP*(T_H(J_MXWL+1)-T_H(J_MXWL))+T_H(J_MXWL)
CC      D_MXWL = YINTRP*(D_H(J_MXWL+1)-D_H(J_MXWL))+D_H(J_MXWL)
        T_MXWL = DMIN1(T_H(J_MXWL+1),T_H(J_MXWL))
        IF(T_H(J_MXWL+1).GT.T_H(J_MXWL)) THEN
          D_MXWL = D_H(J_MXWL)
        ELSE
          D_MXWL = D_H(J_MXWL+1)
        ENDIF
C
C
C
C--------------------------------------------------------------------
C            Interpolate to get Maxwell construction densities
C--------------------------------------------------------------------
C
C
C
        DNS_1 = YINTRP*(BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))+
     1               BRYLOW(I_MXWL,J_MXWL)
        DNS_2 = YINTRP*
     1        (BRYLOW(I_MXWL+1,J_MXWL+1)-BRYLOW(I_MXWL+1,J_MXWL))+
     2               BRYLOW(I_MXWL+1,J_MXWL)
C
        LOWDNS = TINTRP*(DNS_2-DNS_1)+DNS_1
C
C                Derivative of lower density w.r.t. T
        DNL_DT = (DNS_2-DNS_1)/DELT_T
C
        DNDY1 = (BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))/DELT_Y
        DNDY2 = (BRYLOW(I_MXWL+1,J_MXWL+1)-
     1      BRYLOW(I_MXWL+1,J_MXWL))/DELT_Y
        DNL_DY = TINTRP*(DNDY2-DNDY1)+DNDY1
C
C
C
C
        IF(YE.GT.Y_CUT) THEN
C
          DNS_1 = YINTRP*
     1        (BRYHI(I_MXWL,J_MXWL+1)-BRYHI(I_MXWL,J_MXWL))+
     2        BRYHI(I_MXWL,J_MXWL)
          DNS_2 = YINTRP*
     1        (BRYHI(I_MXWL+1,J_MXWL+1)-BRYHI(I_MXWL+1,J_MXWL))+
     2               BRYHI(I_MXWL+1,J_MXWL)
C
          HIDNS = TINTRP*(DNS_2-DNS_1)+DNS_1
C
C                Derivative of higher density w.r.t. T
          DNH_DT = (DNS_2-DNS_1)/DELT_T
C
C
        DNDY1 = (BRYHI(I_MXWL,J_MXWL+1)-
     1      BRYHI(I_MXWL,J_MXWL))/DELT_Y
        DNDY2 = (BRYHI(I_MXWL+1,J_MXWL+1)-
     1      BRYHI(I_MXWL+1,J_MXWL))/DELT_Y
        DNH_DY = TINTRP*(DNDY2-DNDY1)+DNDY1
C
C
        ELSE
          HIDNS = LOWDNS
        ENDIF
C
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
C                       Ye is too low
      ELSE
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time'
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_LOW,' instead'
        YE = Y_LOW+1.0D-6
        GOTO 10
      ENDIF
C
C
C
C
C
      DLTLN1 = (LNCUT-LNLOW)/DBLE(NUMLOW-1)
      DLTLN2 = (LNHI-LNCUT)/DBLE(NUMHI-1)
C
C
      NLOW = 10.0**LNLOW
      NHI = 10.0**LNHI
      N_CUT = 10.0**LNCUT
      LOGBRY = DLOG10(BRYDNS)
      LOGBCH = LOGBRY
C
C
C
C
C----------------------------------------------------------
C           Calc T index
C----------------------------------------------------------
C
C
      IF(LOGBRY.GE.LNHI) THEN
        I_BD = NBPNTS
        I_BNDY = NBPNTS
        T_BNDY = YINTRP*
     1           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+
     2            LBOUND(I_BNDY,J_BNDY)
        TCHK_B = 1.01*T_BNDY
        TCHK_N = 0.95*T_BNDY
        GOTO 70
      ELSEIF((LOGBRY.LT.LNHI).AND.(LOGBRY.GT.LNCUT)) THEN
C
        I_BD = INT((LOGBRY-LNCUT)/DLTLN2)+NUMLOW
        LNMINS = LNCUT+DBLE(I_BD-NUMLOW)*DLTLN2
        LNPLUS = LNCUT+DBLE(I_BD-NUMLOW+1)*DLTLN2
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN
          I_BNDY = I_BD
        ELSEIF(LOGBCH.GT.LNPLUS) THEN
          I_BD = I_BD+1
          I_BNDY = I_BD
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2
        ELSE
          I_BD = I_BD-1
          I_BNDY = I_BD
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2
        ENDIF
C
      ELSEIF((LOGBRY.LE.LNCUT).AND.(LOGBRY.GT.LNLOW)) THEN
C
        I_BD = INT((LOGBRY-LNLOW)/DLTLN1)+1
        LNMINS = LNLOW+DBLE(I_BD-1)*DLTLN1
        LNPLUS = LNLOW+DBLE(I_BD)*DLTLN1
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN
          I_BNDY = I_BD
        ELSEIF(LOGBCH.GT.LNPLUS) THEN
          I_BD = I_BD+1
          I_BNDY = I_BD
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1
        ELSE
          I_BD = I_BD-1
          I_BNDY = I_BD
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1
        ENDIF
C
      ENDIF
C
      IF(I_BNDY.GT.(NBPNTS-1)) THEN
        I_BD = NBPNTS-1
        I_BNDY = I_BD
        LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2
        LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2
      ENDIF
C
C
C
      LMM = LBOUND(I_BNDY,J_BNDY)
      LPM = LBOUND(I_BNDY+1,J_BNDY)
      LMP = LBOUND(I_BNDY,J_BNDY+1)
      LPP = LBOUND(I_BNDY+1,J_BNDY+1)
C
      LNFRAC = (LOGBCH-LNMINS)/(LNPLUS-LNMINS)
C
C                Interpolate in Ye first
C
      TEMP_1 = YINTRP*
     1           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+
     2            LBOUND(I_BNDY,J_BNDY)
      TEMP_2 = YINTRP*
     1        (LBOUND(I_BNDY+1,J_BNDY+1)-LBOUND(I_BNDY+1,J_BNDY))+
     2               LBOUND(I_BNDY+1,J_BNDY)
C
C                Interpolate in density between the two Ye
C                interpolated values
C
      T_BNDY = LNFRAC*(TEMP_2-TEMP_1)+TEMP_1
C
      TCHK_B = 1.01*T_BNDY
      TCHK_N = 0.95*T_BNDY
C
      IF((LMM.GE.LPM).OR.(LMP.GT.LPP)) THEN
        TCHK_N = DMAX1(0.0D0,DMIN1(0.95*TCHK_N,T_BNDY-3.0))
      ENDIF
C
 70   CONTINUE
C
C----------------------------------------------------------
C----------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
C               EOS Logic
C-----------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C                     If T is below the maximum for maxwel construction
      IF(T.LT.T_MXWL) THEN
C                       If rho is greater than the upper max. con.
C                       density the use the bulk EOS
        IF(BRYDNS.GT.HIDNS) THEN
          EOSFLG = 3
C                       Else if rho is greater than the lower max. con.
C                       density then
        ELSEIF(BRYDNS.GT.LOWDNS) THEN
C                         If Ye is large enough to have a signifigant
C                         max con then use the maxwell con. EOS
          IF(YE.GT.Y_CUT) THEN
            EOSFLG = 4
C                         Otherwise use the bulk EOS
          ELSE
            EOSFLG = 3
          ENDIF
C
C                       If density is greater than the minimum
C                       Maxwell con. density, then we know that we are
C                       in the Nuclear EOS density
        ELSEIF(BRYDNS.GT.D_MXWL) THEN
          EOSFLG = 2
C
C
C                       Otherwise check the Boundary table
        ELSE
C
C                         If T is well below the phase boundary curve
C                         then use the nuclear EOS
          IF(T.LT.TCHK_N) THEN
            EOSFLG = 2
C                         Otherwise if T is near the boundary, first
C                         try the nuclear EOS and if not successfull
C                         then use the bulk EOS
          ELSEIF(T.LT.TCHK_B) THEN
            EOSFLG = 1
          ELSE
C                         Otherwise T is well above the boundary so
C                         use the bulk EOS
            EOSFLG = 3
          ENDIF
        ENDIF
C
C                     Otherwise T is above the maximum for a maxwell
C                     construction
      ELSE
C                       If density is greater than that at the top of
C                       the maxwell construction then use the bulk EOS
        IF(BRYDNS.GT.D_MXWL) THEN
          EOSFLG = 3
C
C                       Otherwise density is below the maxwell con.
        ELSE
C
C                         If T is well below the phase boundary curve
C                         then use the nuclear EOS
          IF(T.LT.TCHK_N) THEN
            EOSFLG = 2
C
C                         Otherwise if T is near the phase boundary
C                         curve then try the nuclear EOS and if not
C                         successfull then use the bulk EOS
          ELSEIF(T.LT.TCHK_B) THEN
            EOSFLG = 1
C
C                         Otherwise T is well above the phase boundary
C                         curve so use the bulk EOS
          ELSE
            EOSFLG = 3
          ENDIF
        ENDIF
      ENDIF  
C
C
C-----------------------------------------------------------------------
C                         Done with EOS logic so return EOSFLG
C-----------------------------------------------------------------------
C
 999  RETURN
C
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         RESET.FOR
C
C***********************************************************************
C
C    MODULE:       RESET
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         12/21/90
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C
C    CALL LINE:    CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)
C
C
C    INPUTS:       INPVAR = TEMP, NSUBI, ETA_PO, ETA_NO
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C
C
C    OUTPUTS:      OUTVAR = ARRAY OF LENGTH 4 CONTAINING RESET VALUES
C                  FOR THE INITIAL GUESSES
C
C
C
C 
C    INCLUDE FILES: NONE
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE RESET(INPVAR,YE,BRYDNS,OUTVAR)
C
C
      IMPLICIT NONE
C
C
C                      Subroutine parameters
C
      DOUBLE PRECISION INPVAR(4), OUTVAR(4), YE, BRYDNS
C
C
C                      Local variables
C
      DOUBLE PRECISION ZPG, ZNG, ETA_PG, ETA_NG, PI, UQ, MQ, T, EFRAC
C
C                      Functions
C
      DOUBLE PRECISION FINV12
C
C-----------------------------------------------------------------------
C
      T = INPVAR(1)
C
C
      PI = 3.1415927
      UQ = 20.721
C
      MQ = (T/UQ)**1.5
C
C
      EFRAC = 0.5*YE
C
      ZNG = 2.0*(PI**2)*BRYDNS*(1.0-EFRAC)/MQ
C
      ZPG = 2.0*(PI**2)*BRYDNS*EFRAC/MQ
C
      ETA_NG = FINV12(ZNG)
C
      ETA_PG = FINV12(ZPG)
C
      OUTVAR(1) = INPVAR(1)
      OUTVAR(2) = INPVAR(2)
      OUTVAR(3) = ETA_PG
      OUTVAR(4) = ETA_NG
C
C
C-----------------------------------------------------------------------
C
 999  RETURN
C
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         ALOADMX.FOR
C    MODULE:       LOADMX
C    TYPE:         LOADMX
C
C    PURPOSE:      LOAD THE LOOK-UP TABLE FOR THE MAXWELL CONSTRUCTION
C
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/16/90
C
C    CALL LINE:    CALL LOADMX
C
C    INPUTS:       N/A
C
C    OUTPUTS       N/A
C
C    SUBROUTINE CALLS: EOS_M4C
C 
C    INCLUDE FILES:  EOS_M4C.INC, MAXWEL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE LOADMX()
C
      IMPLICIT NONE
C
C
      INCLUDE 'eos_m4c.inc'
      INCLUDE 'maxwel.inc'
C
      INTEGER NTMP, NYE, NYE2, NUM_BP
      INTEGER LUN1, LUN2, KK, KMIN
      PARAMETER(LUN1=54,LUN2=55)
C
C
      INTEGER FNML1, FNML2
      CHARACTER*60 FNAME1, FNAME2
C
      DOUBLE PRECISION N_SM, SYMM_M, COMP_M, BINDEM, SYM_SM, SIG_SM
      DOUBLE PRECISION N_SB, SYMM_B, COMP_B, BINDEB, SYM_SB, SIG_SB
      DOUBLE PRECISION MSCDD3, BSCDD3
C
      INCLUDE 'force.inc'
C
C
C
cc      CALL GETFNM(FNAME1,FNML1,'Enter ASCII Maxwell fname:',26)
      FNAME1 = 'maxwel.atb'
      FNML1=10
C
cc      CALL GETFNM(FNAME2,FNML2,'Enter ASCII boundary fname:',27)
      FNAME2 = 'bound.atb'
      FNML2=9
C
C
C
C-----------------------------------------------------------------------
C        Read the file Maxwell construction data file
C-----------------------------------------------------------------------
C
C
C
      OPEN(UNIT=LUN1,FILE=FNAME1(1:FNML1),STATUS='OLD')
C
C
C
C
C
C
      READ(LUN1,*) N_SM, SYMM_M
      READ(LUN1,*) COMP_M,BINDEM
      READ(LUN1,*) SYM_SM, SIG_SM
C
C
C
C
      READ(LUN1,*) NTMP,NYE
      READ(LUN1,*) T_LOW,T_HI
      READ(LUN1,*) Y_LOW,Y_HI
C
C
C
      IF((NTMP.NE.NUMTMP).OR.(NYE.NE.NUMYE)) THEN
        WRITE(*,*) 'LOADMX:  MXWL TABLE IS INCOMPATIBLE WITH ARRAYS'
        STOP
      ENDIF
C
C
      DO 101 J=1,NUMYE,1
        DO 100 I=1,NUMTMP,3
          KMIN = MIN0(I+2,NUMTMP)
          READ(LUN1,*) (BRYLOW(KK,J),KK=I,KMIN,1)
 100    CONTINUE
 101  CONTINUE
C
C
      DO 103 J=1,NUMYE,1
        DO 102 I=1,NUMTMP,3
          KMIN = MIN0(I+2,NUMTMP)
          READ(LUN1,*) (BRYHI(KK,J),KK=I,KMIN,1)
 102    CONTINUE
 103  CONTINUE
C
C
C
      DO 104 I=1,NUMYE,3
        KMIN = MIN0(I+2,NUMYE)
        READ(LUN1,*) (T_H(KK),KK=I,KMIN,1)
 104  CONTINUE
C
C
      DO 105 I=1,NUMYE,3
        KMIN = MIN0(I+2,NUMYE)
        READ(LUN1,*) (D_H(KK),KK=I,KMIN,1)
 105  CONTINUE
C
      READ(LUN1,*) YCUT
      READ(LUN1,*) MSCDD3
C
C
      CLOSE(UNIT=LUN1,STATUS='KEEP')
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  MAXWELL CON. TABLE IS INITIALIZED>>'
      WRITE(*,*)
C
C
C
C-----------------------------------------------------------------------
C        Read the file Boundary data file
C-----------------------------------------------------------------------
C
C
C
      OPEN(UNIT=LUN2,FILE=FNAME2(1:FNML2),STATUS='OLD')
C
C
C
C
      READ(LUN2,*) N_SB,SYMM_B
      READ(LUN2,*) COMP_B,BINDEB
      READ(LUN2,*) SYM_SB,SIG_SB
C
C
C
C
C
C
      READ(LUN2,*) NUM_BP,NYE2
      READ(LUN2,*) LNL,LNH,LNC
      READ(LUN2,*) Y_LOW2,Y_HI2
C
C
      IF((NBPNTS.NE.NUM_BP).OR.(NYE2.NE.NUMYE)) THEN
        WRITE(*,*) 'LOADMX:  BNDY TABLE IS INCOMPATIBLE WITH ARRAYS'
        STOP
      ENDIF
C
      IF(ABS(LNL-LNLOW).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  LOWER END OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
C
      IF(ABS(LNH-LNHI).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  UPPER END OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
C
      IF(ABS(LNC-LNCUT).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  MID CUT OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
      IF(ABS(Y_LOW-Y_LOW2).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  LOWER YE LIMITS ARE INCONSIST.'
        STOP
      ENDIF
C
      IF(ABS(Y_HI-Y_HI2).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  UPPER YE LIMITS ARE INCONSIST.'
        STOP
      ENDIF
C
C
      DO 201 J=1,NUMYE,1
        DO 200 I=1,NBPNTS,3
          KMIN = MIN0(I+2,NBPNTS)
          READ(LUN2,*) (LBOUND(KK,J),KK=I,KMIN,1)
 200    CONTINUE
 201  CONTINUE
C
C
      DO 203 J=1,NUMYE,1
        DO 202 I=1,NBPNTS,3
          KMIN = MIN0(I+2,NBPNTS)
          READ(LUN2,*) (UBOUND(KK,J),KK=I,KMIN,1)
 202    CONTINUE
 203  CONTINUE
C
      READ(LUN2,*) BSCDD3
C
      IF(ABS(MSCDD3-BSCDD3).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  SCRDD3 VALUES ARE INCONSIST.'
        STOP
      ENDIF
C
C
C
      CLOSE(UNIT=LUN2,STATUS='KEEP')
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  BOUNDARY TABLE IS INITIALIZED>>'
      WRITE(*,*)
C
C
C
C-----------------------------------------------------------------------
C                  All arrays are now loaded so return
C-----------------------------------------------------------------------
C
      SCRDD3 = BSCDD3
      N_S = N_SM
      NSUBS = N_SM
      SYMM = SYMM_M
      COMP = COMP_M
      BIND_E = BINDEM
      SYM_S = SYM_SM
      SIG_S = SIG_SM
C
c20      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23
c20      DD = (COMP+2.0*SKYRMC)/(3.0*SKYRMC+9.0*BIND_E)
c20      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S
c20      AA = (OVR23*SKYRMC-DD*(SKYRMC+BIND_E))/(N_S*(DD-1.0))-BB
c20      CC = (COMP+2.0*SKYRMC)/(9.0*DD*(DD-1.0)*N_S**DD)
      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23
      DD = (COMP+2.0*SKYRMC+SCRDD3*(COMP-4.0*SKYRMC-18.0*BIND_E))/
     1    ((1.0-SCRDD3)*(3.0*SKYRMC+9.0*BIND_E))
      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S
      CC = ((OVR3*SKYRMC+BIND_E)*(1.0+SCRDD3)**2)/((N_S**DD)*(DD-1.0))
      AA = ((OVR23*SKYRMC-DD*(SKYRMC+BIND_E)-SCRDD3*(OVR3*SKYRMC+BIND_E)
     1    )/(N_S*(DD-1.0)) )-BB
      DD3 = SCRDD3/(N_S**(DD-1.0))
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  SKYRME PARAMETERS FOR THIS RUN ARE:>>'
      WRITE(*,*) 'ABCD: ',AA,BB,CC,DD,SCRDD3
      WRITE(*,*) ' Satur. density, symmetry engy, & compression mod.:'
      WRITE(*,*) N_SM, SYMM_M, COMP_M
      WRITE(*,*) N_SB, SYMM_B, COMP_B
      WRITE(*,*) ' Binding engy, surf. symm. engy, & surface tension:'
      WRITE(*,*) BINDEM,SYM_SM,SIG_SM
      WRITE(*,*) BINDEB,SYM_SB,SIG_SB
C
      WRITE(*,*)
C
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX: FERMI INTEGRAL TABLES ARE INITIALIZED>>'
      WRITE(*,*)
C
C
 999  RETURN
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       GETFNM.FOR
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         8/5/89
C
C    PURPOSE:      OBTAINS A FILE NAME FROM USER
C
C    CALL LINE:    CALL GETFNM(FNAME,FNML,PROMPT,PROMTL)
C
C    INPUTS:       PROMPT = STRING TO POMPT USER WITH (C*60)
C                  PROMTL = LENGTH OF PROMPT STRING (I)
C
C    OUTPUTS:      FNAME = FILE NAME (C*60)
C                  FNML = FILE NAME LENGTH (I)
C*************************************************************************
C
      SUBROUTINE GETFNM(FNAME,FNML,PROMPT,PROMTL)
C
      IMPLICIT NONE
C
      INTEGER FNML, PROMTL
      CHARACTER*60 FNAME, PROMPT
C
C                       Local variables
C
      INTEGER STDOUT, STDIN, I
      DATA STDOUT/6/, STDIN/5/
C
C                       Prompt user for file name
C
      WRITE(STDOUT,'(T2,A,$)') PROMPT(1:PROMTL)
      READ(STDIN,'(A)') FNAME
C
C                        Figure out input file name length
      DO 10 I=1,20,1
        IF(FNAME(I:I).EQ.' ') GOTO 20
 10   CONTINUE
C
 20   CONTINUE
      FNML = I-1
C
C
 999  RETURN
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_1_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    CALL LINE:    F_1_2(Y)      (1/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       1/2th Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_1_2(Y)
      IMPLICIT NONE
C                     
      INTEGER NLOW,NHIGH,N
      PARAMETER(N=201)
      DOUBLE PRECISION Y, A(7), ETA(N), F12(N), F12A(N), TH
      DOUBLE PRECISION F0, F1, X2
C
C                Retain values between calls
      SAVE NLOW,NHIGH
C
      DATA A/6.16850274D0,1.77568655D0,6.92965606D0,
     $   .176776695D0,6.41500299D-02,.4D0,1.32934039D0/
      DATA TH,NLOW,NHIGH/0.33333333333D0,1,N/
C
C
C                    Cubic Spline data
      DATA ETA(  1),F12(  1),F12A(  1)
     $/-1.00000000E+01, 4.02339983E-05, 0.000000000000E+00/
      DATA ETA(  2),F12(  2),F12A(  2)
     $/-9.80000000E+00, 4.91417786E-05, 5.971803206329E-05/
      DATA ETA(  3),F12(  3),F12A(  3)
     $/-9.60000000E+00, 6.00216308E-05, 5.693865674682E-05/
      DATA ETA(  4),F12(  4),F12A(  4)
     $/-9.40000000E+00, 7.33101788E-05, 7.383171094946E-05/
      DATA ETA(  5),F12(  5),F12A(  5)
     $/-9.20000000E+00, 8.95406574E-05, 8.902408945532E-05/
      DATA ETA(  6),F12(  6),F12A(  6)
     $/-9.00000000E+00, 1.09364317E-04, 1.090490812293E-04/
      DATA ETA(  7),F12(  7),F12A(  7)
     $/-8.80000000E+00, 1.33576705E-04, 1.330888456275E-04/
      DATA ETA(  8),F12(  8),F12A(  8)
     $/-8.60000000E+00, 1.63148990E-04, 1.625800862606E-04/
      DATA ETA(  9),F12(  9),F12A(  9)
     $/-8.40000000E+00, 1.99267717E-04, 1.985571093302E-04/
      DATA ETA( 10),F12( 10),F12A( 10)
     $/-8.20000000E+00, 2.43381819E-04, 2.424977264186E-04/
      DATA ETA( 11),F12( 11),F12A( 11)
     $/-8.00000000E+00, 2.97260769E-04, 2.961791849954E-04/
      DATA ETA( 12),F12( 12),F12A( 12)
     $/-7.80000000E+00, 3.63065723E-04, 3.616861335996E-04/
      DATA ETA( 13),F12( 13),F12A( 13)
     $/-7.60000000E+00, 4.43435136E-04, 4.417451306063E-04/
      DATA ETA( 14),F12( 14),F12A( 14)
     $/-7.40000000E+00, 5.41591857E-04, 5.394295439752E-04/
      DATA ETA( 15),F12( 15),F12A( 15)
     $/-7.20000000E+00, 6.61470011E-04, 6.587516434928E-04/
      DATA ETA( 16),F12( 16),F12A( 16)
     $/-7.00000000E+00, 8.07873978E-04, 8.044358320536E-04/
      DATA ETA( 17),F12( 17),F12A( 17)
     $/-6.80000000E+00, 9.86669358E-04, 9.822169782926E-04/
      DATA ETA( 18),F12( 18),F12A( 18)
     $/-6.60000000E+00, 1.20501561E-03, 1.199327054776E-03/
      DATA ETA( 19),F12( 19),F12A( 19)
     $/-6.40000000E+00, 1.47165314E-03, 1.464166502604E-03/
      DATA ETA( 20),F12( 20),F12A( 20)
     $/-6.20000000E+00, 1.79724747E-03, 1.787526934808E-03/
      DATA ETA( 21),F12( 21),F12A( 21)
     $/-6.00000000E+00, 2.19481561E-03, 2.181797258165E-03/
      DATA ETA( 22),F12( 22),F12A( 22)
     $/-5.80000000E+00, 2.68023420E-03, 2.662851532532E-03/
      DATA ETA( 23),F12( 23),F12A( 23)
     $/-5.60000000E+00, 3.27287085E-03, 3.249505611708E-03/
      DATA ETA( 24),F12( 24),F12A( 24)
     $/-5.40000000E+00, 3.99634114E-03, 3.964172020638E-03/
      DATA ETA( 25),F12( 25),F12A( 25)
     $/-5.20000000E+00, 4.87942214E-03, 4.835412805738E-03/
      DATA ETA( 26),F12( 26),F12A( 26)
     $/-5.00000000E+00, 5.95717982E-03, 5.895678756409E-03/
      DATA ETA( 27),F12( 27),F12A( 27)
     $/-4.80000000E+00, 7.27229903E-03, 7.186101668628E-03/
      DATA ETA( 28),F12( 28),F12A( 28)
     $/-4.60000000E+00, 8.87672161E-03, 8.755420069081E-03/
      DATA ETA( 29),F12( 29),F12A( 29)
     $/-4.40000000E+00, 1.08335980E-02, 1.066028955505E-02/
      DATA ETA( 30),F12( 30),F12A( 30)
     $/-4.20000000E+00, 1.32195998E-02, 1.297223321072E-02/
      DATA ETA( 31),F12( 31),F12A( 31)
     $/-4.00000000E+00, 1.61277414E-02, 1.577174760206E-02/
      DATA ETA( 32),F12( 32),F12A( 32)
     $/-3.80000000E+00, 1.96706583E-02, 1.915707138105E-02/
      DATA ETA( 33),F12( 33),F12A( 33)
     $/-3.60000000E+00, 2.39845128E-02, 2.324060687373E-02/
      DATA ETA( 34),F12( 34),F12A( 34)
     $/-3.40000000E+00, 2.92335232E-02, 2.815388612401E-02/
      DATA ETA( 35),F12( 35),F12A( 35)
     $/-3.20000000E+00, 3.56152142E-02, 3.404593863024E-02/
      DATA ETA( 36),F12( 36),F12A( 36)
     $/-3.00000000E+00, 4.33663810E-02, 4.108372935503E-02/
      DATA ETA( 37),F12( 37),F12A( 37)
     $/-2.80000000E+00, 5.27697736E-02, 4.945301394966E-02/
      DATA ETA( 38),F12( 38),F12A( 38)
     $/-2.60000000E+00, 6.41614366E-02, 5.934477484635E-02/
      DATA ETA( 39),F12( 39),F12A( 39)
     $/-2.40000000E+00, 7.79383797E-02, 7.095990166492E-02/
      DATA ETA( 40),F12( 40),F12A( 40)
     $/-2.20000000E+00, 9.45664588E-02, 8.448601849395E-02/
      DATA ETA( 41),F12( 41),F12A( 41)
     $/-2.00000000E+00, 1.14587830E-01, 1.000898393593E-01/
      DATA ETA( 42),F12( 42),F12A( 42)
     $/-1.80000000E+00, 1.38627354E-01, 1.178775440690E-01/
      DATA ETA( 43),F12( 43),F12A( 43)
     $/-1.60000000E+00, 1.67396817E-01, 1.378908343646E-01/
      DATA ETA( 44),F12( 44),F12A( 44)
     $/-1.40000000E+00, 2.01696221E-01, 1.600502684727E-01/
      DATA ETA( 45),F12( 45),F12A( 45)
     $/-1.20000000E+00, 2.42410529E-01, 1.841436917447E-01/
      DATA ETA( 46),F12( 46),F12A( 46)
     $/-1.00000000E+00, 2.90500917E-01, 2.097869645484E-01/
      DATA ETA( 47),F12( 47),F12A( 47)
     $/-8.00000000E-01, 3.46989460E-01, 2.364317000617E-01/
      DATA ETA( 48),F12( 48),F12A( 48)
     $/-6.00000000E-01, 4.12937023E-01, 2.633392352050E-01/
      DATA ETA( 49),F12( 49),F12A( 49)
     $/-4.00000000E-01, 4.89414580E-01, 2.897104591184E-01/
      DATA ETA( 50),F12( 50),F12A( 50)
     $/-2.00000000E-01, 5.77470496E-01, 3.145727783213E-01/
      DATA ETA( 51),F12( 51),F12A( 51)
     $/ 0.00000000E+00, 6.78093925E-01, 3.371253775963E-01/
      DATA ETA( 52),F12( 52),F12A( 52)
     $/ 2.00000000E-01, 7.92181447E-01, 3.565396612936E-01/
      DATA ETA( 53),F12( 53),F12A( 53)
     $/ 4.00000000E-01, 9.20506015E-01, 3.722728772295E-01/
      DATA ETA( 54),F12( 54),F12A( 54)
     $/ 6.00000000E-01, 1.06369475E+00, 3.839938797886E-01/
      DATA ETA( 55),F12( 55),F12A( 55)
     $/ 8.00000000E-01, 1.22221592E+00, 3.916168536161E-01/
      DATA ETA( 56),F12( 56),F12A( 56)
     $/ 1.00000000E+00, 1.39637545E+00, 3.952927057472E-01/
      DATA ETA( 57),F12( 57),F12A( 57)
     $/ 1.20000000E+00, 1.58632329E+00, 3.954588233952E-01/
      DATA ETA( 58),F12( 58),F12A( 58)
     $/ 1.40000000E+00, 1.79206851E+00, 3.924790006721E-01/
      DATA ETA( 59),F12( 59),F12A( 59)
     $/ 1.60000000E+00, 2.01349622E+00, 3.869986739163E-01/
      DATA ETA( 60),F12( 60),F12A( 60)
     $/ 1.80000000E+00, 2.25039083E+00, 3.795613036627E-01/
      DATA ETA( 61),F12( 61),F12A( 61)
     $/ 2.00000000E+00, 2.50245792E+00, 3.706281114328E-01/
      DATA ETA( 62),F12( 62),F12A( 62)
     $/ 2.20000000E+00, 2.76934439E+00, 3.608332506061E-01/
      DATA ETA( 63),F12( 63),F12A( 63)
     $/ 2.40000000E+00, 3.05065972E+00, 3.503678861428E-01/
      DATA ETA( 64),F12( 64),F12A( 64)
     $/ 2.60000000E+00, 3.34598833E+00, 3.396872048225E-01/
      DATA ETA( 65),F12( 65),F12A( 65)
     $/ 2.80000000E+00, 3.65490490E+00, 3.290772945673E-01/
      DATA ETA( 66),F12( 66),F12A( 66)
     $/ 3.00000000E+00, 3.97698528E+00, 3.185751169082E-01/
      DATA ETA( 67),F12( 67),F12A( 67)
     $/ 3.20000000E+00, 4.31181109E+00, 3.084367378000E-01/
      DATA ETA( 68),F12( 68),F12A( 68)
     $/ 3.40000000E+00, 4.65897715E+00, 2.987154318921E-01/
      DATA ETA( 69),F12( 69),F12A( 69)
     $/ 3.60000000E+00, 5.01809514E+00, 2.894910346313E-01/
      DATA ETA( 70),F12( 70),F12A( 70)
     $/ 3.80000000E+00, 5.38879550E+00, 2.806759295830E-01/
      DATA ETA( 71),F12( 71),F12A( 71)
     $/ 4.00000000E+00, 5.77072680E+00, 2.724462470369E-01/
      DATA ETA( 72),F12( 72),F12A( 72)
     $/ 4.20000000E+00, 6.16355908E+00, 2.646860822693E-01/
      DATA ETA( 73),F12( 73),F12A( 73)
     $/ 4.40000000E+00, 6.56698239E+00, 2.574639238859E-01/
      DATA ETA( 74),F12( 74),F12A( 74)
     $/ 4.60000000E+00, 6.98070586E+00, 2.504822221872E-01/
      DATA ETA( 75),F12( 75),F12A( 75)
     $/ 4.80000000E+00, 7.40445435E+00, 2.443601873650E-01/
      DATA ETA( 76),F12( 76),F12A( 76)
     $/ 5.00000000E+00, 7.83797658E+00, 2.381380283526E-01/
      DATA ETA( 77),F12( 77),F12A( 77)
     $/ 5.20000000E+00, 8.28102899E+00, 2.326146992243E-01/
      DATA ETA( 78),F12( 78),F12A( 78)
     $/ 5.40000000E+00, 8.73338914E+00, 2.275641747505E-01/
      DATA ETA( 79),F12( 79),F12A( 79)
     $/ 5.60000000E+00, 9.19485021E+00, 2.222666017741E-01/
      DATA ETA( 80),F12( 80),F12A( 80)
     $/ 5.80000000E+00, 9.66520906E+00, 2.180364181530E-01/
      DATA ETA( 81),F12( 81),F12A( 81)
     $/ 6.00000000E+00, 1.01442864E+01, 2.133612256141E-01/
      DATA ETA( 82),F12( 82),F12A( 82)
     $/ 6.20000000E+00, 1.06319029E+01, 2.093926793908E-01/
      DATA ETA( 83),F12( 83),F12A( 83)
     $/ 6.40000000E+00, 1.11278965E+01, 2.056330568230E-01/
      DATA ETA( 84),F12( 84),F12A( 84)
     $/ 6.60000000E+00, 1.16321146E+01, 2.017500933175E-01/
      DATA ETA( 85),F12( 85),F12A( 85)
     $/ 6.80000000E+00, 1.21444066E+01, 1.984515699067E-01/
      DATA ETA( 86),F12( 86),F12A( 86)
     $/ 7.00000000E+00, 1.26646369E+01, 1.951886270556E-01/
      DATA ETA( 87),F12( 87),F12A( 87)
     $/ 7.20000000E+00, 1.31926749E+01, 1.919489218711E-01/
      DATA ETA( 88),F12( 88),F12A( 88)
     $/ 7.40000000E+00, 1.37283938E+01, 1.891506854598E-01/
      DATA ETA( 89),F12( 89),F12A( 89)
     $/ 7.60000000E+00, 1.42716773E+01, 1.861383362905E-01/
      DATA ETA( 90),F12( 90),F12A( 90)
     $/ 7.80000000E+00, 1.48224099E+01, 1.836609693776E-01/
      DATA ETA( 91),F12( 91),F12A( 91)
     $/ 8.00000000E+00, 1.53804867E+01, 1.808477861992E-01/
      DATA ETA( 92),F12( 92),F12A( 92)
     $/ 8.20000000E+00, 1.59458020E+01, 1.787228858260E-01/
      DATA ETA( 93),F12( 93),F12A( 93)
     $/ 8.40000000E+00, 1.65182614E+01, 1.758756704959E-01/
      DATA ETA( 94),F12( 94),F12A( 94)
     $/ 8.60000000E+00, 1.70977635E+01, 1.741794321910E-01/
      DATA ETA( 95),F12( 95),F12A( 95)
     $/ 8.80000000E+00, 1.76842275E+01, 1.716916007395E-01/
      DATA ETA( 96),F12( 96),F12A( 96)
     $/ 9.00000000E+00, 1.82775617E+01, 1.695841648522E-01/
      DATA ETA( 97),F12( 97),F12A( 97)
     $/ 9.20000000E+00, 1.88776803E+01, 1.676317398519E-01/
      DATA ETA( 98),F12( 98),F12A( 98)
     $/ 9.40000000E+00, 1.94845052E+01, 1.658338757391E-01/
      DATA ETA( 99),F12( 99),F12A( 99)
     $/ 9.60000000E+00, 2.00979619E+01, 1.638027571928E-01/
      DATA ETA(100),F12(100),F12A(100)
     $/ 9.80000000E+00, 2.07179742E+01, 1.622950954886E-01/
      DATA ETA(101),F12(101),F12A(101)
     $/ 1.00000000E+01, 2.13444734E+01, 1.600518608534E-01/
      DATA ETA(102),F12(102),F12A(102)
     $/ 1.02000000E+01, 2.19773812E+01, 1.587874610986E-01/
      DATA ETA(103),F12(103),F12A(103)
     $/ 1.04000000E+01, 2.26166368E+01, 1.569682947511E-01/
      DATA ETA(104),F12(104),F12A(104)
     $/ 1.06000000E+01, 2.32621732E+01, 1.554593598982E-01/
      DATA ETA(105),F12(105),F12A(105)
     $/ 1.08000000E+01, 2.39139295E+01, 1.541792656549E-01/
      DATA ETA(106),F12(106),F12A(106)
     $/ 1.10000000E+01, 2.45718484E+01, 1.522135774835E-01/
      DATA ETA(107),F12(107),F12A(107)
     $/ 1.12000000E+01, 2.52358613E+01, 1.510664244108E-01/
      DATA ETA(108),F12(108),F12A(108)
     $/ 1.14000000E+01, 2.59059148E+01, 1.496107248727E-01/
      DATA ETA(109),F12(109),F12A(109)
     $/ 1.16000000E+01, 2.65819535E+01, 1.482706760991E-01/
      DATA ETA(110),F12(110),F12A(110)
     $/ 1.18000000E+01, 2.72639241E+01, 1.470915707300E-01/
      DATA ETA(111),F12(111),F12A(111)
     $/ 1.20000000E+01, 2.79517770E+01, 1.457080409816E-01/
      DATA ETA(112),F12(112),F12A(112)
     $/ 1.22000000E+01, 2.86454568E+01, 1.441112653439E-01/
      DATA ETA(113),F12(113),F12A(113)
     $/ 1.24000000E+01, 2.93449082E+01, 1.435868976411E-01/
      DATA ETA(114),F12(114),F12A(114)
     $/ 1.26000000E+01, 3.00500951E+01, 1.418661440931E-01/
      DATA ETA(115),F12(115),F12A(115)
     $/ 1.28000000E+01, 3.07609620E+01, 1.409485259852E-01/
      DATA ETA(116),F12(116),F12A(116)
     $/ 1.30000000E+01, 3.14774652E+01, 1.397847519674E-01/
      DATA ETA(117),F12(117),F12A(117)
     $/ 1.32000000E+01, 3.21995592E+01, 1.385324661454E-01/
      DATA ETA(118),F12(118),F12A(118)
     $/ 1.34000000E+01, 3.29271975E+01, 1.377303834490E-01/
      DATA ETA(119),F12(119),F12A(119)
     $/ 1.36000000E+01, 3.36603403E+01, 1.362210000606E-01/
      DATA ETA(120),F12(120),F12A(120)
     $/ 1.38000000E+01, 3.43989420E+01, 1.362206163067E-01/
      DATA ETA(121),F12(121),F12A(121)
     $/ 1.40000000E+01, 3.51429758E+01, 1.337115347148E-01/
      DATA ETA(122),F12(122),F12A(122)
     $/ 1.42000000E+01, 3.58923769E+01, 1.340282448335E-01/
      DATA ETA(123),F12(123),F12A(123)
     $/ 1.44000000E+01, 3.66471262E+01, 1.324054859504E-01/
      DATA ETA(124),F12(124),F12A(124)
     $/ 1.46000000E+01, 3.74071779E+01, 1.317098113657E-01/
      DATA ETA(125),F12(125),F12A(125)
     $/ 1.48000000E+01, 3.81724978E+01, 1.309852685862E-01/
      DATA ETA(126),F12(126),F12A(126)
     $/ 1.50000000E+01, 3.89430513E+01, 1.293891142900E-01/
      DATA ETA(127),F12(127),F12A(127)
     $/ 1.52000000E+01, 3.97187891E+01, 1.291032742542E-01/
      DATA ETA(128),F12(128),F12A(128)
     $/ 1.54000000E+01, 4.04996882E+01, 1.283927886917E-01/
      DATA ETA(129),F12(129),F12A(129)
     $/ 1.56000000E+01, 4.12857180E+01, 1.269305709802E-01/
      DATA ETA(130),F12(130),F12A(130)
     $/ 1.58000000E+01, 4.20768328E+01, 1.266349273854E-01/
      DATA ETA(131),F12(131),F12A(131)
     $/ 1.60000000E+01, 4.28730059E+01, 1.252747194805E-01/
      DATA ETA(132),F12(132),F12A(132)
     $/ 1.62000000E+01, 4.36741991E+01, 1.252811946920E-01/
      DATA ETA(133),F12(133),F12A(133)
     $/ 1.64000000E+01, 4.44803934E+01, 1.237655017505E-01/
      DATA ETA(134),F12(134),F12A(134)
     $/ 1.66000000E+01, 4.52915468E+01, 1.235217983048E-01/
      DATA ETA(135),F12(135),F12A(135)
     $/ 1.68000000E+01, 4.61076365E+01, 1.225923050320E-01/
      DATA ETA(136),F12(136),F12A(136)
     $/ 1.70000000E+01, 4.69286280E+01, 1.213789815681E-01/
      DATA ETA(137),F12(137),F12A(137)
     $/ 1.72000000E+01, 4.77544832E+01, 1.214467686948E-01/
      DATA ETA(138),F12(138),F12A(138)
     $/ 1.74000000E+01, 4.85851870E+01, 1.201239436519E-01/
      DATA ETA(139),F12(139),F12A(139)
     $/ 1.76000000E+01, 4.94207048E+01, 1.201574566968E-01/
      DATA ETA(140),F12(140),F12A(140)
     $/ 1.78000000E+01, 5.02610178E+01, 1.185262295614E-01/
      DATA ETA(141),F12(141),F12A(141)
     $/ 1.80000000E+01, 5.11060801E+01, 1.181326250585E-01/
      DATA ETA(142),F12(142),F12A(142)
     $/ 1.82000000E+01, 5.19558725E+01, 1.184582702039E-01/
      DATA ETA(143),F12(143),F12A(143)
     $/ 1.84000000E+01, 5.28103838E+01, 1.158692941259E-01/
      DATA ETA(144),F12(144),F12A(144)
     $/ 1.86000000E+01, 5.36695566E+01, 1.172895532908E-01/
      DATA ETA(145),F12(145),F12A(145)
     $/ 1.88000000E+01, 5.45334024E+01, 1.159224927128E-01/
      DATA ETA(146),F12(146),F12A(146)
     $/ 1.90000000E+01, 5.54018793E+01, 1.136854758577E-01/
      DATA ETA(147),F12(147),F12A(147)
     $/ 1.92000000E+01, 5.62749338E+01, 1.159756038573E-01/
      DATA ETA(148),F12(148),F12A(148)
     $/ 1.94000000E+01, 5.71525927E+01, 1.130721087129E-01/
      DATA ETA(149),F12(149),F12A(149)
     $/ 1.96000000E+01, 5.80347986E+01, 1.137859612890E-01/
      DATA ETA(150),F12(150),F12A(150)
     $/ 1.98000000E+01, 5.89215441E+01, 1.127240461327E-01/
      DATA ETA(151),F12(151),F12A(151)
     $/ 2.00000000E+01, 5.98127985E+01, 1.116528541808E-01/
      DATA ETA(152),F12(152),F12A(152)
     $/ 2.02000000E+01, 6.07085314E+01, 1.124395371436E-01/
      DATA ETA(153),F12(153),F12A(153)
     $/ 2.04000000E+01, 6.16087389E+01, 1.097789972451E-01/
      DATA ETA(154),F12(154),F12A(154)
     $/ 2.06000000E+01, 6.25133677E+01, 1.116394738743E-01/
      DATA ETA(155),F12(155),F12A(155)
     $/ 2.08000000E+01, 6.34224367E+01, 1.096931072587E-01/
      DATA ETA(156),F12(156),F12A(156)
     $/ 2.10000000E+01, 6.43359013E+01, 1.089280970924E-01/
      DATA ETA(157),F12(157),F12A(157)
     $/ 2.12000000E+01, 6.52537327E+01, 1.096145043705E-01/
      DATA ETA(158),F12(158),F12A(158)
     $/ 2.14000000E+01, 6.61759281E+01, 1.072138854249E-01/
      DATA ETA(159),F12(159),F12A(159)
     $/ 2.16000000E+01, 6.71024418E+01, 1.092749539287E-01/
      DATA ETA(160),F12(160),F12A(160)
     $/ 2.18000000E+01, 6.80332966E+01, 1.068512988636E-01/
      DATA ETA(161),F12(161),F12A(161)
     $/ 2.20000000E+01, 6.89684391E+01, 1.064748506159E-01/
      DATA ETA(162),F12(162),F12A(162)
     $/ 2.22000000E+01, 6.99078465E+01, 1.069842986745E-01/
      DATA ETA(163),F12(163),F12A(163)
     $/ 2.24000000E+01, 7.08515186E+01, 1.052929546852E-01/
      DATA ETA(164),F12(164),F12A(164)
     $/ 2.26000000E+01, 7.17994175E+01, 1.058638825818E-01/
      DATA ETA(165),F12(165),F12A(165)
     $/ 2.28000000E+01, 7.27515431E+01, 1.052565149906E-01/
      DATA ETA(166),F12(166),F12A(166)
     $/ 2.30000000E+01, 7.37078724E+01, 1.036650574545E-01/
      DATA ETA(167),F12(167),F12A(167)
     $/ 2.32000000E+01, 7.46683674E+01, 1.049382551915E-01/
      DATA ETA(168),F12(168),F12A(168)
     $/ 2.34000000E+01, 7.56330357E+01, 1.025769217803E-01/
      DATA ETA(169),F12(169),F12A(169)
     $/ 2.36000000E+01, 7.66018314E+01, 1.038640576851E-01/
      DATA ETA(170),F12(170),F12A(170)
     $/ 2.38000000E+01, 7.75747700E+01, 1.034018474827E-01/
      DATA ETA(171),F12(171),F12A(171)
     $/ 2.40000000E+01, 7.85518284E+01, 1.004985523838E-01/
      DATA ETA(172),F12(172),F12A(172)
     $/ 2.42000000E+01, 7.95329456E+01, 1.034239429815E-01/
      DATA ETA(173),F12(173),F12A(173)
     $/ 2.44000000E+01, 8.05181599E+01, 1.003706756894E-01/
      DATA ETA(174),F12(174),F12A(174)
     $/ 2.46000000E+01, 8.15074177E+01, 1.016183542594E-01/
      DATA ETA(175),F12(175),F12A(175)
     $/ 2.48000000E+01, 8.25007267E+01, 1.008359072745E-01/
      DATA ETA(176),F12(176),F12A(176)
     $/ 2.50000000E+01, 8.34980640E+01, 9.928301664272E-02/
      DATA ETA(177),F12(177),F12A(177)
     $/ 2.52000000E+01, 8.44993916E+01, 1.005770261543E-01/
      DATA ETA(178),F12(178),F12A(178)
     $/ 2.54000000E+01, 8.55047245E+01, 9.920387874048E-02/
      DATA ETA(179),F12(179),F12A(179)
     $/ 2.56000000E+01, 8.65140324E+01, 9.885745888275E-02/
      DATA ETA(180),F12(180),F12A(180)
     $/ 2.58000000E+01, 8.75272999E+01, 9.930628572935E-02/
      DATA ETA(181),F12(181),F12A(181)
     $/ 2.60000000E+01, 8.85445194E+01, 9.671739819958E-02/
      DATA ETA(182),F12(182),F12A(182)
     $/ 2.62000000E+01, 8.95656452E+01, 9.976912147403E-02/
      DATA ETA(183),F12(183),F12A(183)
     $/ 2.64000000E+01, 9.05907154E+01, 9.586611590140E-02/
      DATA ETA(184),F12(184),F12A(184)
     $/ 2.66000000E+01, 9.16196613E+01, 9.812141492108E-02/
      DATA ETA(185),F12(185),F12A(185)
     $/ 2.68000000E+01, 9.26525059E+01, 9.645322441392E-02/
      DATA ETA(186),F12(186),F12A(186)
     $/ 2.70000000E+01, 9.36892185E+01, 9.626568742494E-02/
      DATA ETA(187),F12(187),F12A(187)
     $/ 2.72000000E+01, 9.47297840E+01, 9.641902588607E-02/
      DATA ETA(188),F12(188),F12A(188)
     $/ 2.74000000E+01, 9.57741947E+01, 9.483820903090E-02/
      DATA ETA(189),F12(189),F12A(189)
     $/ 2.76000000E+01, 9.68224201E+01, 9.643313798854E-02/
      DATA ETA(190),F12(190),F12A(190)
     $/ 2.78000000E+01, 9.78744755E+01, 9.392923901678E-02/
      DATA ETA(191),F12(191),F12A(191)
     $/ 2.80000000E+01, 9.89303150E+01, 9.546490594284E-02/
      DATA ETA(192),F12(192),F12A(192)
     $/ 2.82000000E+01, 9.99899540E+01, 9.413613721391E-02/
      DATA ETA(193),F12(193),F12A(193)
     $/ 2.84000000E+01, 1.01053362E+02, 9.334054520234E-02/
      DATA ETA(194),F12(194),F12A(194)
     $/ 2.86000000E+01, 1.02120516E+02, 9.440168197065E-02/
      DATA ETA(195),F12(195),F12A(195)
     $/ 2.88000000E+01, 1.03191431E+02, 9.320272691961E-02/
      DATA ETA(196),F12(196),F12A(196)
     $/ 2.90000000E+01, 1.04266077E+02, 9.243741035120E-02/
      DATA ETA(197),F12(197),F12A(197)
     $/ 2.92000000E+01, 1.05344431E+02, 9.324763167608E-02/
      DATA ETA(198),F12(198),F12A(198)
     $/ 2.94000000E+01, 1.06426508E+02, 9.302206294354E-02/
      DATA ETA(199),F12(199),F12A(199)
     $/ 2.96000000E+01, 1.07512262E+02, 8.621411654839E-02/
      DATA ETA(200),F12(200),F12A(200)
     $/ 2.98000000E+01, 1.08601709E+02, 1.160714708634E-01/
      DATA ETA(201),F12(201),F12A(201)
     $/ 3.00000000E+01, 1.09694826E+02, 0.000000000000E+00/

C
C
C          If Y is between -10 and 30 then use spline table
      IF((Y.LE.30.0D0).AND.(Y.GE.-10.0D0)) THEN
        CALL INTRP(Y,F1,ETA,F12,F12A,N,NLOW,NHIGH)
C          Else if Y is greater than 30 use degenerate approximation
      ELSEIF(Y.GT.30.0D0) THEN
        X2=Y**(-2)
        F1=A(6)*Y*SQRT(Y)*TH*(5.0D0+(A(1)+(3*A(2)+7.0D0*X2*A(3))*X2)*X2)  
C          Else if Y is less than -10 use nondegenerate approximation
      ELSE
        F0=DEXP(Y)   
        F1=A(7)*TH*F0*(2.0D0-(4.0D0*A(4)-(6.0D0*A(5)-0.25D0*F0)*F0)*F0) 
      ENDIF
C
      F_1_2=F1  
 999  RETURN
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_3_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    CALL LINE:    F_3_2(Y)      (3/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       3/2th Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_3_2(Y)
      IMPLICIT NONE 
      INTEGER N, NLOW, NHIGH
      PARAMETER(N=201)
      DOUBLE PRECISION Y, A(7), ETA(N), F32(N), F32A(N), X2, F0, F1
C
C                Retain values between calls
      SAVE NLOW,NHIGH
C
      DATA A/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0,
     $ 6.41500299D-02,.4D0,1.32934039D0/
      DATA NLOW, NHIGH/1,N/
C
C
C                    Cubic Spline data
      DATA ETA(  1),F32(  1),F32A(  1)
     $/-1.00000000E+01, 6.03514791E-05, 0.000000000000E+00/
      DATA ETA(  2),F32(  2),F32A(  2)
     $/-9.80000000E+00, 7.37133833E-05, 8.958056137096E-05/
      DATA ETA(  3),F32(  3),F32A(  3)
     $/-9.60000000E+00, 9.00335163E-05, 8.541207451613E-05/
      DATA ETA(  4),F32(  4),F32A(  4)
     $/-9.40000000E+00, 1.09966868E-04, 1.107539455645E-04/
      DATA ETA(  5),F32(  5),F32A(  5)
     $/-9.20000000E+00, 1.34313381E-04, 1.335463382257E-04/
      DATA ETA(  6),F32(  6),F32A(  6)
     $/-9.00000000E+00, 1.64050060E-04, 1.635856015326E-04/
      DATA ETA(  7),F32(  7),F32A(  7)
     $/-8.80000000E+00, 2.00370374E-04, 1.996565056440E-04/
      DATA ETA(  8),F32(  8),F32A(  8)
     $/-8.60000000E+00, 2.44731440E-04, 2.439011758913E-04/
      DATA ETA(  9),F32(  9),F32A(  9)
     $/-8.40000000E+00, 2.98913459E-04, 2.978817407908E-04/
      DATA ETA( 10),F32( 10),F32A( 10)
     $/-8.20000000E+00, 3.65090447E-04, 3.638172109455E-04/
      DATA ETA( 11),F32( 11),F32A( 11)
     $/-8.00000000E+00, 4.45917576E-04, 4.443705654272E-04/
      DATA ETA( 12),F32( 12),F32A( 12)
     $/-7.80000000E+00, 5.44637980E-04, 5.426917773457E-04/
      DATA ETA( 13),F32( 13),F32A( 13)
     $/-7.60000000E+00, 6.65211541E-04, 6.628358751901E-04/
      DATA ETA( 14),F32( 14),F32A( 14)
     $/-7.40000000E+00, 8.12475410E-04, 8.095109218942E-04/
      DATA ETA( 15),F32( 15),F32A( 15)
     $/-7.20000000E+00, 9.92335991E-04, 9.886272372330E-04/
      DATA ETA( 16),F32( 16),F32A( 16)
     $/-7.00000000E+00, 1.21200623E-03, 1.207428829174E-03/
      DATA ETA( 17),F32( 17),F32A( 17)
     $/-6.80000000E+00, 1.48029524E-03, 1.474473096071E-03/
      DATA ETA( 18),F32( 18),F32A( 18)
     $/-6.60000000E+00, 1.80795780E-03, 1.800711286541E-03/
      DATA ETA( 19),F32( 19),F32A( 19)
     $/-6.40000000E+00, 2.20812770E-03, 2.198782757766E-03/
      DATA ETA( 20),F32( 20),F32A( 20)
     $/-6.20000000E+00, 2.69683736E-03, 2.685121682395E-03/
      DATA ETA( 21),F32( 21),F32A( 21)
     $/-6.00000000E+00, 3.29366449E-03, 3.278351012652E-03/
      DATA ETA( 22),F32( 22),F32A( 22)
     $/-5.80000000E+00, 4.02250013E-03, 4.002750766996E-03/
      DATA ETA( 23),F32( 23),F32A( 23)
     $/-5.60000000E+00, 4.91251063E-03, 4.886874919365E-03/
      DATA ETA( 24),F32( 24),F32A( 24)
     $/-5.40000000E+00, 5.99928957E-03, 5.965015555545E-03/
      DATA ETA( 25),F32( 25),F32A( 25)
     $/-5.20000000E+00, 7.32625567E-03, 7.281136858453E-03/
      DATA ETA( 26),F32( 26),F32A( 26)
     $/-5.00000000E+00, 8.94638640E-03, 8.885131510642E-03/
      DATA ETA( 27),F32( 27),F32A( 27)
     $/-4.80000000E+00, 1.09242697E-02, 1.084122259898E-02/
      DATA ETA( 28),F32( 28),F32A( 28)
     $/-4.60000000E+00, 1.33386545E-02, 1.322520309344E-02/
      DATA ETA( 29),F32( 29),F32A( 29)
     $/-4.40000000E+00, 1.62855038E-02, 1.612764002725E-02/
      DATA ETA( 30),F32( 30),F32A( 30)
     $/-4.20000000E+00, 1.98816718E-02, 1.966204179756E-02/
      DATA ETA( 31),F32( 31),F32A( 31)
     $/-4.00000000E+00, 2.42694099E-02, 2.395970778252E-02/
      DATA ETA( 32),F32( 32),F32A( 32)
     $/-3.80000000E+00, 2.96217115E-02, 2.918365207237E-02/
      DATA ETA( 33),F32( 33),F32A( 33)
     $/-3.60000000E+00, 3.61488024E-02, 3.552407892800E-02/
      DATA ETA( 34),F32( 34),F32A( 34)
     $/-3.40000000E+00, 4.41058287E-02, 4.321034221562E-02/
      DATA ETA( 35),F32( 35),F32A( 35)
     $/-3.20000000E+00, 5.38020553E-02, 5.251459720952E-02/
      DATA ETA( 36),F32( 36),F32A( 36)
     $/-3.00000000E+00, 6.56117518E-02, 6.375175394630E-02/
      DATA ETA( 37),F32( 37),F32A( 37)
     $/-2.80000000E+00, 7.99869169E-02, 7.729867700529E-02/
      DATA ETA( 38),F32( 38),F32A( 38)
     $/-2.60000000E+00, 9.74722300E-02, 9.357573803256E-02/
      DATA ETA( 39),F32( 39),F32A( 39)
     $/-2.40000000E+00, 1.18722076E-01, 1.130783058645E-01/
      DATA ETA( 40),F32( 40),F32A( 40)
     $/-2.20000000E+00, 1.44520123E-01, 1.363411885096E-01/
      DATA ETA( 41),F32( 41),F32A( 41)
     $/-2.00000000E+00, 1.75800983E-01, 1.639788900970E-01/
      DATA ETA( 42),F32( 42),F32A( 42)
     $/-1.80000000E+00, 2.13674326E-01, 1.966157011023E-01/
      DATA ETA( 43),F32( 43),F32A( 43)
     $/-1.60000000E+00, 2.59450115E-01, 2.349252054937E-01/
      DATA ETA( 44),F32( 44),F32A( 44)
     $/-1.40000000E+00, 3.14665116E-01, 2.795652769229E-01/
      DATA ETA( 45),F32( 45),F32A( 45)
     $/-1.20000000E+00, 3.81109037E-01, 3.311516868147E-01/
      DATA ETA( 46),F32( 46),F32A( 46)
     $/-1.00000000E+00, 4.60848816E-01, 3.902066758184E-01/
      DATA ETA( 47),F32( 47),F32A( 47)
     $/-8.00000000E-01, 5.56249276E-01, 4.571237599118E-01/
      DATA ETA( 48),F32( 48),F32A( 48)
     $/-6.00000000E-01, 6.69988349E-01, 5.320902345345E-01/
      DATA ETA( 49),F32( 49),F32A( 49)
     $/-4.00000000E-01, 8.05064514E-01, 6.150791019500E-01/
      DATA ETA( 50),F32( 50),F32A( 50)
     $/-2.00000000E-01, 9.64795128E-01, 7.057607076653E-01/
      DATA ETA( 51),F32( 51),F32A( 51)
     $/ 0.00000000E+00, 1.15280381E+00, 8.035882673888E-01/
      DATA ETA( 52),F32( 52),F32A( 52)
     $/ 2.00000000E-01, 1.37299815E+00, 9.077349227796E-01/
      DATA ETA( 53),F32( 53),F32A( 53)
     $/ 4.00000000E-01, 1.62953690E+00, 1.017133541493E+00/
      DATA ETA( 54),F32( 54),F32A( 54)
     $/ 6.00000000E-01, 1.92678872E+00, 1.130691411249E+00/
      DATA ETA( 55),F32( 55),F32A( 55)
     $/ 8.00000000E-01, 2.26928741E+00, 1.247131313511E+00/
      DATA ETA( 56),F32( 56),F32A( 56)
     $/ 1.00000000E+00, 2.66168267E+00, 1.365268834706E+00/
      DATA ETA( 57),F32( 57),F32A( 57)
     $/ 1.20000000E+00, 3.10869223E+00, 1.483938347666E+00/
      DATA ETA( 58),F32( 58),F32A( 58)
     $/ 1.40000000E+00, 3.61505681E+00, 1.602230774630E+00/
      DATA ETA( 59),F32( 59),F32A( 59)
     $/ 1.60000000E+00, 4.18550170E+00, 1.719185053814E+00/
      DATA ETA( 60),F32( 60),F32A( 60)
     $/ 1.80000000E+00, 4.82470143E+00, 1.834255010113E+00/
      DATA ETA( 61),F32( 61),F32A( 61)
     $/ 2.00000000E+00, 5.53725398E+00, 1.946717905735E+00/
      DATA ETA( 62),F32( 62),F32A( 62)
     $/ 2.20000000E+00, 6.32765782E+00, 2.056566866947E+00/
      DATA ETA( 63),F32( 63),F32A( 63)
     $/ 2.40000000E+00, 7.20030272E+00, 2.163173626476E+00/
      DATA ETA( 64),F32( 64),F32A( 64)
     $/ 2.60000000E+00, 8.15945459E+00, 2.266784127147E+00/
      DATA ETA( 65),F32( 65),F32A( 65)
     $/ 2.80000000E+00, 9.20925546E+00, 2.367039864934E+00/
      DATA ETA( 66),F32( 66),F32A( 66)
     $/ 3.00000000E+00, 1.03537161E+01, 2.464021913115E+00/
      DATA ETA( 67),F32( 67),F32A( 67)
     $/ 3.20000000E+00, 1.15967200E+01, 2.558361482607E+00/
      DATA ETA( 68),F32( 68),F32A( 68)
     $/ 3.40000000E+00, 1.29420350E+01, 2.649197156458E+00/
      DATA ETA( 69),F32( 69),F32A( 69)
     $/ 3.60000000E+00, 1.43933013E+01, 2.737544891559E+00/
      DATA ETA( 70),F32( 70),F32A( 70)
     $/ 3.80000000E+00, 1.59540503E+01, 2.823028277308E+00/
      DATA ETA( 71),F32( 71),F32A( 71)
     $/ 4.00000000E+00, 1.76277032E+01, 2.905926999207E+00/
      DATA ETA( 72),F32( 72),F32A( 72)
     $/ 4.20000000E+00, 1.94175764E+01, 2.986308725864E+00/
      DATA ETA( 73),F32( 73),F32A( 73)
     $/ 4.40000000E+00, 2.13268934E+01, 3.065408097336E+00/
      DATA ETA( 74),F32( 74),F32A( 74)
     $/ 4.60000000E+00, 2.33587976E+01, 3.140138884795E+00/
      DATA ETA( 75),F32( 75),F32A( 75)
     $/ 4.80000000E+00, 2.55163160E+01, 3.216166363482E+00/
      DATA ETA( 76),F32( 76),F32A( 76)
     $/ 5.00000000E+00, 2.78024469E+01, 3.287070661277E+00/
      DATA ETA( 77),F32( 77),F32A( 77)
     $/ 5.20000000E+00, 3.02200609E+01, 3.358015991410E+00/
      DATA ETA( 78),F32( 78),F32A( 78)
     $/ 5.40000000E+00, 3.27719889E+01, 3.427965373083E+00/
      DATA ETA( 79),F32( 79),F32A( 79)
     $/ 5.60000000E+00, 3.54610072E+01, 3.493667516261E+00/
      DATA ETA( 80),F32( 80),F32A( 80)
     $/ 5.80000000E+00, 3.82897883E+01, 3.561784561870E+00/
      DATA ETA( 81),F32( 81),F32A( 81)
     $/ 6.00000000E+00, 4.12610064E+01, 3.624744236258E+00/
      DATA ETA( 82),F32( 82),F32A( 82)
     $/ 6.20000000E+00, 4.43772212E+01, 3.688743493097E+00/
      DATA ETA( 83),F32( 83),F32A( 83)
     $/ 6.40000000E+00, 4.76409770E+01, 3.751431791354E+00/
      DATA ETA( 84),F32( 84),F32A( 84)
     $/ 6.60000000E+00, 5.10547763E+01, 3.812054341488E+00/
      DATA ETA( 85),F32( 85),F32A( 85)
     $/ 6.80000000E+00, 5.46210566E+01, 3.872500842692E+00/
      DATA ETA( 86),F32( 86),F32A( 86)
     $/ 7.00000000E+00, 5.83422213E+01, 3.930602287744E+00/
      DATA ETA( 87),F32( 87),F32A( 87)
     $/ 7.20000000E+00, 6.22206164E+01, 3.989650006333E+00/
      DATA ETA( 88),F32( 88),F32A( 88)
     $/ 7.40000000E+00, 6.62585851E+01, 4.046837686926E+00/
      DATA ETA( 89),F32( 89),F32A( 89)
     $/ 7.60000000E+00, 7.04584142E+01, 4.102059245965E+00/
      DATA ETA( 90),F32( 90),F32A( 90)
     $/ 7.80000000E+00, 7.48223363E+01, 4.158875329212E+00/
      DATA ETA( 91),F32( 91),F32A( 91)
     $/ 8.00000000E+00, 7.93525945E+01, 4.212854437187E+00/
      DATA ETA( 92),F32( 92),F32A( 92)
     $/ 8.20000000E+00, 8.40513631E+01, 4.266266922047E+00/
      DATA ETA( 93),F32( 93),F32A( 93)
     $/ 8.40000000E+00, 8.89207860E+01, 4.320222874619E+00/
      DATA ETA( 94),F32( 94),F32A( 94)
     $/ 8.60000000E+00, 9.39630071E+01, 4.372571579482E+00/
      DATA ETA( 95),F32( 95),F32A( 95)
     $/ 8.80000000E+00, 9.91801320E+01, 4.425060807447E+00/
      DATA ETA( 96),F32( 96),F32A( 96)
     $/ 9.00000000E+00, 1.04574244E+02, 4.475250190736E+00/
      DATA ETA( 97),F32( 97),F32A( 97)
     $/ 9.20000000E+00, 1.10147364E+02, 4.525138429613E+00/
      DATA ETA( 98),F32( 98),F32A( 98)
     $/ 9.40000000E+00, 1.15901507E+02, 4.577646090805E+00/
      DATA ETA( 99),F32( 99),F32A( 99)
     $/ 9.60000000E+00, 1.21838717E+02, 4.624327207175E+00/
      DATA ETA(100),F32(100),F32A(100)
     $/ 9.80000000E+00, 1.27960940E+02, 4.676995080485E+00/
      DATA ETA(101),F32(101),F32A(101)
     $/ 1.00000000E+01, 1.34270176E+02, 4.719642470896E+00/
      DATA ETA(102),F32(102),F32A(102)
     $/ 1.02000000E+01, 1.40768269E+02, 4.772985035935E+00/
      DATA ETA(103),F32(103),F32A(103)
     $/ 1.04000000E+01, 1.47457218E+02, 4.816817385355E+00/
      DATA ETA(104),F32(104),F32A(104)
     $/ 1.06000000E+01, 1.54338871E+02, 4.865345422653E+00/
      DATA ETA(105),F32(105),F32A(105)
     $/ 1.08000000E+01, 1.61415135E+02, 4.913450924021E+00/
      DATA ETA(106),F32(106),F32A(106)
     $/ 1.10000000E+01, 1.68687886E+02, 4.953900881277E+00/
      DATA ETA(107),F32(107),F32A(107)
     $/ 1.12000000E+01, 1.76158863E+02, 5.004845550873E+00/
      DATA ETA(108),F32(108),F32A(108)
     $/ 1.14000000E+01, 1.83829975E+02, 5.046966915220E+00/
      DATA ETA(109),F32(109),F32A(109)
     $/ 1.16000000E+01, 1.91702992E+02, 5.093036788259E+00/
      DATA ETA(110),F32(110),F32A(110)
     $/ 1.18000000E+01, 1.99779728E+02, 5.138735931733E+00/
      DATA ETA(111),F32(111),F32A(111)
     $/ 1.20000000E+01, 2.08061970E+02, 5.177919484819E+00/
      DATA ETA(112),F32(112),F32A(112)
     $/ 1.22000000E+01, 2.16551396E+02, 5.227186128993E+00/
      DATA ETA(113),F32(113),F32A(113)
     $/ 1.24000000E+01, 2.25249821E+02, 5.263185999195E+00/
      DATA ETA(114),F32(114),F32A(114)
     $/ 1.26000000E+01, 2.34158879E+02, 5.315019874243E+00/
      DATA ETA(115),F32(115),F32A(115)
     $/ 1.28000000E+01, 2.43280430E+02, 5.350684503817E+00/
      DATA ETA(116),F32(116),F32A(116)
     $/ 1.30000000E+01, 2.52616062E+02, 5.394392110504E+00/
      DATA ETA(117),F32(117),F32A(117)
     $/ 1.32000000E+01, 2.62167458E+02, 5.436347054167E+00/
      DATA ETA(118),F32(118),F32A(118)
     $/ 1.34000000E+01, 2.71936318E+02, 5.479819672817E+00/
      DATA ETA(119),F32(119),F32A(119)
     $/ 1.36000000E+01, 2.81924325E+02, 5.516424254575E+00/
      DATA ETA(120),F32(120),F32A(120)
     $/ 1.38000000E+01, 2.92133065E+02, 5.564433308870E+00/
      DATA ETA(121),F32(121),F32A(121)
     $/ 1.40000000E+01, 3.02564278E+02, 5.596792509966E+00/
      DATA ETA(122),F32(122),F32A(122)
     $/ 1.42000000E+01, 3.13219429E+02, 5.639096651265E+00/
      DATA ETA(123),F32(123),F32A(123)
     $/ 1.44000000E+01, 3.24100167E+02, 5.684870884964E+00/
      DATA ETA(124),F32(124),F32A(124)
     $/ 1.46000000E+01, 3.35208199E+02, 5.715519808889E+00/
      DATA ETA(125),F32(125),F32A(125)
     $/ 1.48000000E+01, 3.46544991E+02, 5.767049879467E+00/
      DATA ETA(126),F32(126),F32A(126)
     $/ 1.50000000E+01, 3.58112282E+02, 5.791130673258E+00/
      DATA ETA(127),F32(127),F32A(127)
     $/ 1.52000000E+01, 3.69911385E+02, 5.840227427499E+00/
      DATA ETA(128),F32(128),F32A(128)
     $/ 1.54000000E+01, 3.81944008E+02, 5.875959616731E+00/
      DATA ETA(129),F32(129),F32A(129)
     $/ 1.56000000E+01, 3.94211678E+02, 5.912984105596E+00/
      DATA ETA(130),F32(130),F32A(130)
     $/ 1.58000000E+01, 4.06715920E+02, 5.957903960865E+00/
      DATA ETA(131),F32(131),F32A(131)
     $/ 1.60000000E+01, 4.19458352E+02, 5.983900050964E+00/
      DATA ETA(132),F32(132),F32A(132)
     $/ 1.62000000E+01, 4.32440285E+02, 6.031645835280E+00/
      DATA ETA(133),F32(133),F32A(133)
     $/ 1.64000000E+01, 4.45663338E+02, 6.057516607911E+00/
      DATA ETA(134),F32(134),F32A(134)
     $/ 1.66000000E+01, 4.59128884E+02, 6.112237733043E+00/
      DATA ETA(135),F32(135),F32A(135)
     $/ 1.68000000E+01, 4.72838723E+02, 6.137482459949E+00/
      DATA ETA(136),F32(136),F32A(136)
     $/ 1.70000000E+01, 4.86794106E+02, 6.169432427160E+00/
      DATA ETA(137),F32(137),F32A(137)
     $/ 1.72000000E+01, 5.00996376E+02, 6.217837831413E+00/
      DATA ETA(138),F32(138),F32A(138)
     $/ 1.74000000E+01, 5.15447221E+02, 6.245466247192E+00/
      DATA ETA(139),F32(139),F32A(139)
     $/ 1.76000000E+01, 5.30147965E+02, 6.285147179775E+00/
      DATA ETA(140),F32(140),F32A(140)
     $/ 1.78000000E+01, 5.45100114E+02, 6.324695033745E+00/
      DATA ETA(141),F32(141),F32A(141)
     $/ 1.80000000E+01, 5.60305131E+02, 6.346272685250E+00/
      DATA ETA(142),F32(142),F32A(142)
     $/ 1.82000000E+01, 5.75764237E+02, 6.403564225258E+00/
      DATA ETA(143),F32(143),F32A(143)
     $/ 1.84000000E+01, 5.91479142E+02, 6.409320413722E+00/
      DATA ETA(144),F32(144),F32A(144)
     $/ 1.86000000E+01, 6.07450822E+02, 6.475404119805E+00/
      DATA ETA(145),F32(145),F32A(145)
     $/ 1.88000000E+01, 6.23681230E+02, 6.498263107104E+00/
      DATA ETA(146),F32(146),F32A(146)
     $/ 1.90000000E+01, 6.40171525E+02, 6.514593451770E+00/
      DATA ETA(147),F32(147),F32A(147)
     $/ 1.92000000E+01, 6.56922746E+02, 6.582263085821E+00/
      DATA ETA(148),F32(148),F32A(148)
     $/ 1.94000000E+01, 6.73936845E+02, 6.588054204946E+00/
      DATA ETA(149),F32(149),F32A(149)
     $/ 1.96000000E+01, 6.91214799E+02, 6.643770094354E+00/
      DATA ETA(150),F32(150),F32A(150)
     $/ 1.98000000E+01, 7.08758256E+02, 6.662315417678E+00/
      DATA ETA(151),F32(151),F32A(151)
     $/ 2.00000000E+01, 7.26568315E+02, 6.697268234941E+00/
      DATA ETA(152),F32(152),F32A(152)
     $/ 2.02000000E+01, 7.44646317E+02, 6.740061642552E+00/
      DATA ETA(153),F32(153),F32A(153)
     $/ 2.04000000E+01, 7.62993791E+02, 6.763285194860E+00/
      DATA ETA(154),F32(154),F32A(154)
     $/ 2.06000000E+01, 7.81611955E+02, 6.810297577950E+00/
      DATA ETA(155),F32(155),F32A(155)
     $/ 2.08000000E+01, 8.00502335E+02, 6.827924493401E+00/
      DATA ETA(156),F32(156),F32A(156)
     $/ 2.10000000E+01, 8.19665971E+02, 6.866404448437E+00/
      DATA ETA(157),F32(157),F32A(157)
     $/ 2.12000000E+01, 8.39104264E+02, 6.905007712841E+00/
      DATA ETA(158),F32(158),F32A(158)
     $/ 2.14000000E+01, 8.58818620E+02, 6.923014700220E+00/
      DATA ETA(159),F32(159),F32A(159)
     $/ 2.16000000E+01, 8.78810136E+02, 6.976933486221E+00/
      DATA ETA(160),F32(160),F32A(160)
     $/ 2.18000000E+01, 8.99080522E+02, 6.999751354928E+00/
      DATA ETA(161),F32(161),F32A(161)
     $/ 2.20000000E+01, 9.19630815E+02, 7.010111094081E+00/
      DATA ETA(162),F32(162),F32A(162)
     $/ 2.22000000E+01, 9.40461930E+02, 7.083104268755E+00/
      DATA ETA(163),F32(163),F32A(163)
     $/ 2.24000000E+01, 9.61575822E+02, 7.074021830882E+00/
      DATA ETA(164),F32(164),F32A(164)
     $/ 2.26000000E+01, 9.82973161E+02, 7.137858407662E+00/
      DATA ETA(165),F32(165),F32A(165)
     $/ 2.28000000E+01, 1.00465578E+03, 7.166544538534E+00/
      DATA ETA(166),F32(166),F32A(166)
     $/ 2.30000000E+01, 1.02662479E+03, 7.154613438199E+00/
      DATA ETA(167),F32(167),F32A(167)
     $/ 2.32000000E+01, 1.04888077E+03, 7.260501708652E+00/
      DATA ETA(168),F32(168),F32A(168)
     $/ 2.34000000E+01, 1.07142618E+03, 7.217879727200E+00/
      DATA ETA(169),F32(169),F32A(169)
     $/ 2.36000000E+01, 1.09426114E+03, 7.300479382516E+00/
      DATA ETA(170),F32(170),F32A(170)
     $/ 2.38000000E+01, 1.11738761E+03, 7.306702742765E+00/
      DATA ETA(171),F32(171),F32A(171)
     $/ 2.40000000E+01, 1.14080655E+03, 7.343209646437E+00/
      DATA ETA(172),F32(172),F32A(172)
     $/ 2.42000000E+01, 1.16451920E+03, 7.376958671475E+00/
      DATA ETA(173),F32(173),F32A(173)
     $/ 2.44000000E+01, 1.18852677E+03, 7.386955667671E+00/
      DATA ETA(174),F32(174),F32A(174)
     $/ 2.46000000E+01, 1.21283023E+03, 7.458718657792E+00/
      DATA ETA(175),F32(175),F32A(175)
     $/ 2.48000000E+01, 1.23743155E+03, 7.457169701212E+00/
      DATA ETA(176),F32(176),F32A(176)
     $/ 2.50000000E+01, 1.26233133E+03, 7.481602537356E+00/
      DATA ETA(177),F32(177),F32A(177)
     $/ 2.52000000E+01, 1.28753067E+03, 7.550420149389E+00/
      DATA ETA(178),F32(178),F32A(178)
     $/ 2.54000000E+01, 1.31303140E+03, 7.525216865068E+00/
      DATA ETA(179),F32(179),F32A(179)
     $/ 2.56000000E+01, 1.33883389E+03, 7.612712390278E+00/
      DATA ETA(180),F32(180),F32A(180)
     $/ 2.58000000E+01, 1.36494022E+03, 7.599933573862E+00/
      DATA ETA(181),F32(181),F32A(181)
     $/ 2.60000000E+01, 1.39135086E+03, 7.634053314309E+00/
      DATA ETA(182),F32(182),F32A(182)
     $/ 2.62000000E+01, 1.41806705E+03, 7.696353168910E+00/
      DATA ETA(183),F32(183),F32A(183)
     $/ 2.64000000E+01, 1.44509061E+03, 7.686034010016E+00/
      DATA ETA(184),F32(184),F32A(184)
     $/ 2.66000000E+01, 1.47242203E+03, 7.738510790971E+00/
      DATA ETA(185),F32(185),F32A(185)
     $/ 2.68000000E+01, 1.50006278E+03, 7.759422826178E+00/
      DATA ETA(186),F32(186),F32A(186)
     $/ 2.70000000E+01, 1.52801395E+03, 7.786797904321E+00/
      DATA ETA(187),F32(187),F32A(187)
     $/ 2.72000000E+01, 1.55627664E+03, 7.821385556534E+00/
      DATA ETA(188),F32(188),F32A(188)
     $/ 2.74000000E+01, 1.58485208E+03, 7.840159869515E+00/
      DATA ETA(189),F32(189),F32A(189)
     $/ 2.76000000E+01, 1.61374137E+03, 7.895474965363E+00/
      DATA ETA(190),F32(190),F32A(190)
     $/ 2.78000000E+01, 1.64294608E+03, 7.890940269096E+00/
      DATA ETA(191),F32(191),F32A(191)
     $/ 2.80000000E+01, 1.67246671E+03, 7.928763958255E+00/
      DATA ETA(192),F32(192),F32A(192)
     $/ 2.82000000E+01, 1.70230460E+03, 7.983003897891E+00/
      DATA ETA(193),F32(193),F32A(193)
     $/ 2.84000000E+01, 1.73246121E+03, 7.947220450177E+00/
      DATA ETA(194),F32(194),F32A(194)
     $/ 2.86000000E+01, 1.76293668E+03, 8.057114301331E+00/
      DATA ETA(195),F32(195),F32A(195)
     $/ 2.88000000E+01, 1.79373355E+03, 8.034322344553E+00/
      DATA ETA(196),F32(196),F32A(196)
     $/ 2.90000000E+01, 1.82485221E+03, 8.074096320488E+00/
      DATA ETA(197),F32(197),F32A(197)
     $/ 2.92000000E+01, 1.85629361E+03, 8.080292373472E+00/
      DATA ETA(198),F32(198),F32A(198)
     $/ 2.94000000E+01, 1.88805936E+03, 8.257234185640E+00/
      DATA ETA(199),F32(199),F32A(199)
     $/ 2.96000000E+01, 1.92014981E+03, 7.595770883867E+00/
      DATA ETA(200),F32(200),F32A(200)
     $/ 2.98000000E+01, 1.95256705E+03, 1.037818227902E+01/
      DATA ETA(201),F32(201),F32A(201)
     $/ 3.00000000E+01, 1.98531168E+03, 0.000000000000E+00/
C
C          If Y is between -10 and 30 then use spline table
      IF((Y.LE.30.0D0).AND.(Y.GE.-10.0D0)) THEN
        CALL INTRP(Y,F1,ETA,F32,F32A,N,NLOW,NHIGH)
C          Else if Y is greater than 30 use degenerate approximation
      ELSEIF(Y.GT.30.0D0) THEN
        X2=Y**(-2)
        F1=A(6)*SQRT(Y)*(1.0D0/X2+A(1)-(A(2)+X2*A(3))*X2)  
      ELSE
C          Else if Y is less than -10 use nondegenerate approximation
        F0=DEXP(Y)   
        F1=A(7)*F0*(1.0D0-(A(4)-(A(5)-0.03125D0*F0)*F0)*F0) 
      ENDIF
C
      F_3_2=F1  
 999  RETURN
      END   
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FINV12
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    CALL LINE:    FINV12(Y)      (Inverse of the 1/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       Inverse of Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FINV12(Y)
      IMPLICIT NONE
      INTEGER N, NLOW, NHIGH
      PARAMETER(N=201)
      DOUBLE PRECISION Y, AI(8),F12(N),FIA(N),ETA(N), X2, X4, F1
C
C                Retain values between calls
      SAVE NLOW,NHIGH
C
      DATA AI/-.822467032D0,-1.21761363D0,-9.16138616D0,
     $ 0.398942281D0,.0732748216D0,-1.310707D0,1.12837917D0,  
     $ 8.2810645D-3/ 
      DATA NLOW,NHIGH/1,N/
C
C
C                    Cubic Spline data
      DATA ETA(  1),F12(  1),FIA(  1)
     $/-1.00000000E+01, 4.02339983E-05, 0.000000000000E+00/
      DATA ETA(  2),F12(  2),FIA(  2)
     $/-9.80000000E+00, 4.91417786E-05,-5.518729757333E+08/
      DATA ETA(  3),F12(  3),FIA(  3)
     $/-9.60000000E+00, 6.00216308E-05,-2.369113017381E+08/
      DATA ETA(  4),F12(  4),FIA(  4)
     $/-9.40000000E+00, 7.33101788E-05,-1.908760407552E+08/
      DATA ETA(  5),F12(  5),FIA(  5)
     $/-9.20000000E+00, 8.95406574E-05,-1.202176248026E+08/
      DATA ETA(  6),F12(  6),FIA(  6)
     $/-9.00000000E+00, 1.09364317E-04,-8.245454770317E+07/
      DATA ETA(  7),F12(  7),FIA(  7)
     $/-8.80000000E+00, 1.33576705E-04,-5.481585297269E+07/
      DATA ETA(  8),F12(  8),FIA(  8)
     $/-8.60000000E+00, 1.63148990E-04,-3.685695721234E+07/
      DATA ETA(  9),F12(  9),FIA(  9)
     $/-8.40000000E+00, 1.99267717E-04,-2.467974532472E+07/
      DATA ETA( 10),F12( 10),FIA( 10)
     $/-8.20000000E+00, 2.43381819E-04,-1.655022961155E+07/
      DATA ETA( 11),F12( 11),FIA( 11)
     $/-8.00000000E+00, 2.97260769E-04,-1.109333129137E+07/
      DATA ETA( 12),F12( 12),FIA( 12)
     $/-7.80000000E+00, 3.63065723E-04,-7.436285795690E+06/
      DATA ETA( 13),F12( 13),FIA( 13)
     $/-7.60000000E+00, 4.43435136E-04,-4.985364210851E+06/
      DATA ETA( 14),F12( 14),FIA( 14)
     $/-7.40000000E+00, 5.41591857E-04,-3.341805327658E+06/
      DATA ETA( 15),F12( 15),FIA( 15)
     $/-7.20000000E+00, 6.61470011E-04,-2.240407265052E+06/
      DATA ETA( 16),F12( 16),FIA( 16)
     $/-7.00000000E+00, 8.07873978E-04,-1.501972226436E+06/
      DATA ETA( 17),F12( 17),FIA( 17)
     $/-6.80000000E+00, 9.86669358E-04,-1.006911415258E+06/
      DATA ETA( 18),F12( 18),FIA( 18)
     $/-6.60000000E+00, 1.20501561E-03,-6.751030308338E+05/
      DATA ETA( 19),F12( 19),FIA( 19)
     $/-6.40000000E+00, 1.47165314E-03,-4.526103593217E+05/
      DATA ETA( 20),F12( 20),FIA( 20)
     $/-6.20000000E+00, 1.79724747E-03,-3.034927371978E+05/
      DATA ETA( 21),F12( 21),FIA( 21)
     $/-6.00000000E+00, 2.19481561E-03,-2.034893180359E+05/
      DATA ETA( 22),F12( 22),FIA( 22)
     $/-5.80000000E+00, 2.68023420E-03,-1.364626916017E+05/
      DATA ETA( 23),F12( 23),FIA( 23)
     $/-5.60000000E+00, 3.27287085E-03,-9.151916808167E+04/
      DATA ETA( 24),F12( 24),FIA( 24)
     $/-5.40000000E+00, 3.99634114E-03,-6.137988251686E+04/
      DATA ETA( 25),F12( 25),FIA( 25)
     $/-5.20000000E+00, 4.87942214E-03,-4.117642510051E+04/
      DATA ETA( 26),F12( 26),FIA( 26)
     $/-5.00000000E+00, 5.95717982E-03,-2.762388797406E+04/
      DATA ETA( 27),F12( 27),FIA( 27)
     $/-4.80000000E+00, 7.27229903E-03,-1.853723945326E+04/
      DATA ETA( 28),F12( 28),FIA( 28)
     $/-4.60000000E+00, 8.87672161E-03,-1.244247237626E+04/
      DATA ETA( 29),F12( 29),FIA( 29)
     $/-4.40000000E+00, 1.08335980E-02,-8.353189544372E+03/
      DATA ETA( 30),F12( 30),FIA( 30)
     $/-4.20000000E+00, 1.32195998E-02,-5.610484545822E+03/
      DATA ETA( 31),F12( 31),FIA( 31)
     $/-4.00000000E+00, 1.61277414E-02,-3.769624980571E+03/
      DATA ETA( 32),F12( 32),FIA( 32)
     $/-3.80000000E+00, 1.96706583E-02,-2.534200684298E+03/
      DATA ETA( 33),F12( 33),FIA( 33)
     $/-3.60000000E+00, 2.39845128E-02,-1.704676561770E+03/
      DATA ETA( 34),F12( 34),FIA( 34)
     $/-3.40000000E+00, 2.92335232E-02,-1.147572737681E+03/
      DATA ETA( 35),F12( 35),FIA( 35)
     $/-3.20000000E+00, 3.56152142E-02,-7.732385076806E+02/
      DATA ETA( 36),F12( 36),FIA( 36)
     $/-3.00000000E+00, 4.33663810E-02,-5.215910937804E+02/
      DATA ETA( 37),F12( 37),FIA( 37)
     $/-2.80000000E+00, 5.27697736E-02,-3.523253054057E+02/
      DATA ETA( 38),F12( 38),FIA( 38)
     $/-2.60000000E+00, 6.41614366E-02,-2.383626995336E+02/
      DATA ETA( 39),F12( 39),FIA( 39)
     $/-2.40000000E+00, 7.79383797E-02,-1.615785832573E+02/
      DATA ETA( 40),F12( 40),FIA( 40)
     $/-2.20000000E+00, 9.45664588E-02,-1.097815116306E+02/
      DATA ETA( 41),F12( 41),FIA( 41)
     $/-2.00000000E+00, 1.14587830E-01,-7.479631187407E+01/
      DATA ETA( 42),F12( 42),FIA( 42)
     $/-1.80000000E+00, 1.38627354E-01,-5.112413248282E+01/
      DATA ETA( 43),F12( 43),FIA( 43)
     $/-1.60000000E+00, 1.67396817E-01,-3.507905015366E+01/
      DATA ETA( 44),F12( 44),FIA( 44)
     $/-1.40000000E+00, 2.01696221E-01,-2.417700470772E+01/
      DATA ETA( 45),F12( 45),FIA( 45)
     $/-1.20000000E+00, 2.42410529E-01,-1.674984464907E+01/
      DATA ETA( 46),F12( 46),FIA( 46)
     $/-1.00000000E+00, 2.90500917E-01,-1.167337092503E+01/
      DATA ETA( 47),F12( 47),FIA( 47)
     $/-8.00000000E-01, 3.46989460E-01,-8.190719927811E+00/
      DATA ETA( 48),F12( 48),FIA( 48)
     $/-6.00000000E-01, 4.12937023E-01,-5.790648969028E+00/
      DATA ETA( 49),F12( 49),FIA( 49)
     $/-4.00000000E-01, 4.89414580E-01,-4.128944296030E+00/
      DATA ETA( 50),F12( 50),FIA( 50)
     $/-2.00000000E-01, 5.77470496E-01,-2.971061447648E+00/
      DATA ETA( 51),F12( 51),FIA( 51)
     $/ 0.00000000E+00, 6.78093925E-01,-2.159721811273E+00/
      DATA ETA( 52),F12( 52),FIA( 52)
     $/ 2.00000000E-01, 7.92181447E-01,-1.586687451373E+00/
      DATA ETA( 53),F12( 53),FIA( 53)
     $/ 4.00000000E-01, 9.20506015E-01,-1.178969413999E+00/
      DATA ETA( 54),F12( 54),FIA( 54)
     $/ 6.00000000E-01, 1.06369475E+00,-8.863670774363E-01/
      DATA ETA( 55),F12( 55),FIA( 55)
     $/ 8.00000000E-01, 1.22221592E+00,-6.744468494356E-01/
      DATA ETA( 56),F12( 56),FIA( 56)
     $/ 1.00000000E+00, 1.39637545E+00,-5.194863803958E-01/
      DATA ETA( 57),F12( 57),FIA( 57)
     $/ 1.20000000E+00, 1.58632329E+00,-4.051201736579E-01/
      DATA ETA( 58),F12( 58),FIA( 58)
     $/ 1.40000000E+00, 1.79206851E+00,-3.197436337379E-01/
      DATA ETA( 59),F12( 59),FIA( 59)
     $/ 1.60000000E+00, 2.01349622E+00,-2.554204570394E-01/
      DATA ETA( 60),F12( 60),FIA( 60)
     $/ 1.80000000E+00, 2.25039083E+00,-2.064307210907E-01/
      DATA ETA( 61),F12( 61),FIA( 61)
     $/ 2.00000000E+00, 2.50245792E+00,-1.687058961462E-01/
      DATA ETA( 62),F12( 62),FIA( 62)
     $/ 2.20000000E+00, 2.76934439E+00,-1.394151646274E-01/
      DATA ETA( 63),F12( 63),FIA( 63)
     $/ 2.40000000E+00, 3.05065972E+00,-1.163727819558E-01/
      DATA ETA( 64),F12( 64),FIA( 64)
     $/ 2.60000000E+00, 3.34598833E+00,-9.810926234532E-02/
      DATA ETA( 65),F12( 65),FIA( 65)
     $/ 2.80000000E+00, 3.65490490E+00,-8.349866180907E-02/
      DATA ETA( 66),F12( 66),FIA( 66)
     $/ 3.00000000E+00, 3.97698528E+00,-7.167023124865E-02/
      DATA ETA( 67),F12( 67),FIA( 67)
     $/ 3.20000000E+00, 4.31181109E+00,-6.203344907539E-02/
      DATA ETA( 68),F12( 68),FIA( 68)
     $/ 3.40000000E+00, 4.65897715E+00,-5.410768842926E-02/
      DATA ETA( 69),F12( 69),FIA( 69)
     $/ 3.60000000E+00, 5.01809514E+00,-4.753934908271E-02/
      DATA ETA( 70),F12( 70),FIA( 70)
     $/ 3.80000000E+00, 5.38879550E+00,-4.203693557621E-02/
      DATA ETA( 71),F12( 71),FIA( 71)
     $/ 4.00000000E+00, 5.77072680E+00,-3.741510484228E-02/
      DATA ETA( 72),F12( 72),FIA( 72)
     $/ 4.20000000E+00, 6.16355908E+00,-3.349161263478E-02/
      DATA ETA( 73),F12( 73),FIA( 73)
     $/ 4.40000000E+00, 6.56698239E+00,-3.014726514264E-02/
      DATA ETA( 74),F12( 74),FIA( 74)
     $/ 4.60000000E+00, 6.98070586E+00,-2.725049073768E-02/
      DATA ETA( 75),F12( 75),FIA( 75)
     $/ 4.80000000E+00, 7.40445435E+00,-2.478810817325E-02/
      DATA ETA( 76),F12( 76),FIA( 76)
     $/ 5.00000000E+00, 7.83797658E+00,-2.259805758188E-02/
      DATA ETA( 77),F12( 77),FIA( 77)
     $/ 5.20000000E+00, 8.28102899E+00,-2.071312823782E-02/
      DATA ETA( 78),F12( 78),FIA( 78)
     $/ 5.40000000E+00, 8.73338914E+00,-1.906423734281E-02/
      DATA ETA( 79),F12( 79),FIA( 79)
     $/ 5.60000000E+00, 9.19485021E+00,-1.756402815298E-02/
      DATA ETA( 80),F12( 80),FIA( 80)
     $/ 5.80000000E+00, 9.66520906E+00,-1.628975776444E-02/
      DATA ETA( 81),F12( 81),FIA( 81)
     $/ 6.00000000E+00, 1.01442864E+01,-1.510245658209E-02/
      DATA ETA( 82),F12( 82),FIA( 82)
     $/ 6.20000000E+00, 1.06319029E+01,-1.407142912242E-02/
      DATA ETA( 83),F12( 83),FIA( 83)
     $/ 6.40000000E+00, 1.11278965E+01,-1.314229513765E-02/
      DATA ETA( 84),F12( 84),FIA( 84)
     $/ 6.60000000E+00, 1.16321146E+01,-1.228449605398E-02/
      DATA ETA( 85),F12( 85),FIA( 85)
     $/ 6.80000000E+00, 1.21444066E+01,-1.153091235501E-02/
      DATA ETA( 86),F12( 86),FIA( 86)
     $/ 7.00000000E+00, 1.26646369E+01,-1.083803827012E-02/
      DATA ETA( 87),F12( 87),FIA( 87)
     $/ 7.20000000E+00, 1.31926749E+01,-1.019988677668E-02/
      DATA ETA( 88),F12( 88),FIA( 88)
     $/ 7.40000000E+00, 1.37283938E+01,-9.631384626898E-03/
      DATA ETA( 89),F12( 89),FIA( 89)
     $/ 7.60000000E+00, 1.42716773E+01,-9.093442710361E-03/
      DATA ETA( 90),F12( 90),FIA( 90)
     $/ 7.80000000E+00, 1.48224099E+01,-8.618263823991E-03/
      DATA ETA( 91),F12( 91),FIA( 91)
     $/ 8.00000000E+00, 1.53804867E+01,-8.160345480975E-03/
      DATA ETA( 92),F12( 92),FIA( 92)
     $/ 8.20000000E+00, 1.59458020E+01,-7.762531673437E-03/
      DATA ETA( 93),F12( 93),FIA( 93)
     $/ 8.40000000E+00, 1.65182614E+01,-7.360350480125E-03/
      DATA ETA( 94),F12( 94),FIA( 94)
     $/ 8.60000000E+00, 1.70977635E+01,-7.030118584822E-03/
      DATA ETA( 95),F12( 95),FIA( 95)
     $/ 8.80000000E+00, 1.76842275E+01,-6.688634687327E-03/
      DATA ETA( 96),F12( 96),FIA( 96)
     $/ 9.00000000E+00, 1.82775617E+01,-6.382668439251E-03/
      DATA ETA( 97),F12( 97),FIA( 97)
     $/ 9.20000000E+00, 1.88776803E+01,-6.100108206762E-03/
      DATA ETA( 98),F12( 98),FIA( 98)
     $/ 9.40000000E+00, 1.94845052E+01,-5.838946705110E-03/
      DATA ETA( 99),F12( 99),FIA( 99)
     $/ 9.60000000E+00, 2.00979619E+01,-5.584582666939E-03/
      DATA ETA(100),F12(100),FIA(100)
     $/ 9.80000000E+00, 2.07179742E+01,-5.361259022162E-03/
      DATA ETA(101),F12(101),FIA(101)
     $/ 1.00000000E+01, 2.13444734E+01,-5.126492837303E-03/
      DATA ETA(102),F12(102),FIA(102)
     $/ 1.02000000E+01, 2.19773812E+01,-4.934700434242E-03/
      DATA ETA(103),F12(103),FIA(103)
     $/ 1.04000000E+01, 2.26166368E+01,-4.735625230311E-03/
      DATA ETA(104),F12(104),FIA(104)
     $/ 1.06000000E+01, 2.32621732E+01,-4.556040567237E-03/
      DATA ETA(105),F12(105),FIA(105)
     $/ 1.08000000E+01, 2.39139295E+01,-4.391444865609E-03/
      DATA ETA(106),F12(106),FIA(106)
     $/ 1.10000000E+01, 2.45718484E+01,-4.216034760485E-03/
      DATA ETA(107),F12(107),FIA(107)
     $/ 1.12000000E+01, 2.52358613E+01,-4.071259985244E-03/
      DATA ETA(108),F12(108),FIA(108)
     $/ 1.14000000E+01, 2.59059148E+01,-3.924867157171E-03/
      DATA ETA(109),F12(109),FIA(109)
     $/ 1.16000000E+01, 2.65819535E+01,-3.788293366519E-03/
      DATA ETA(110),F12(110),FIA(110)
     $/ 1.18000000E+01, 2.72639241E+01,-3.661781142180E-03/
      DATA ETA(111),F12(111),FIA(111)
     $/ 1.20000000E+01, 2.79517770E+01,-3.535793689118E-03/
      DATA ETA(112),F12(112),FIA(112)
     $/ 1.22000000E+01, 2.86454568E+01,-3.410614036319E-03/
      DATA ETA(113),F12(113),FIA(113)
     $/ 1.24000000E+01, 2.93449082E+01,-3.315348054643E-03/
      DATA ETA(114),F12(114),FIA(114)
     $/ 1.26000000E+01, 3.00500951E+01,-3.196961694843E-03/
      DATA ETA(115),F12(115),FIA(115)
     $/ 1.28000000E+01, 3.07609620E+01,-3.101468403870E-03/
      DATA ETA(116),F12(116),FIA(116)
     $/ 1.30000000E+01, 3.14774652E+01,-3.004323966790E-03/
      DATA ETA(117),F12(117),FIA(117)
     $/ 1.32000000E+01, 3.21995592E+01,-2.909424474011E-03/
      DATA ETA(118),F12(118),FIA(118)
     $/ 1.34000000E+01, 3.29271975E+01,-2.827365948037E-03/
      DATA ETA(119),F12(119),FIA(119)
     $/ 1.36000000E+01, 3.36603403E+01,-2.734508928049E-03/
      DATA ETA(120),F12(120),FIA(120)
     $/ 1.38000000E+01, 3.43989420E+01,-2.674531887656E-03/
      DATA ETA(121),F12(121),FIA(121)
     $/ 1.40000000E+01, 3.51429758E+01,-2.568768499103E-03/
      DATA ETA(122),F12(122),FIA(122)
     $/ 1.42000000E+01, 3.58923769E+01,-2.520360850761E-03/
      DATA ETA(123),F12(123),FIA(123)
     $/ 1.44000000E+01, 3.66471262E+01,-2.437521643781E-03/
      DATA ETA(124),F12(124),FIA(124)
     $/ 1.46000000E+01, 3.74071779E+01,-2.374823121522E-03/
      DATA ETA(125),F12(125),FIA(125)
     $/ 1.48000000E+01, 3.81724978E+01,-2.313510297832E-03/
      DATA ETA(126),F12(126),FIA(126)
     $/ 1.50000000E+01, 3.89430513E+01,-2.239491881527E-03/
      DATA ETA(127),F12(127),FIA(127)
     $/ 1.52000000E+01, 3.97187891E+01,-2.190399198766E-03/
      DATA ETA(128),F12(128),FIA(128)
     $/ 1.54000000E+01, 4.04996882E+01,-2.135558627693E-03/
      DATA ETA(129),F12(129),FIA(129)
     $/ 1.56000000E+01, 4.12857180E+01,-2.070549224394E-03/
      DATA ETA(130),F12(130),FIA(130)
     $/ 1.58000000E+01, 4.20768328E+01,-2.026371946726E-03/
      DATA ETA(131),F12(131),FIA(131)
     $/ 1.60000000E+01, 4.28730059E+01,-1.966932919251E-03/
      DATA ETA(132),F12(132),FIA(132)
     $/ 1.62000000E+01, 4.36741991E+01,-1.930467763422E-03/
      DATA ETA(133),F12(133),FIA(133)
     $/ 1.64000000E+01, 4.44803934E+01,-1.872060355250E-03/
      DATA ETA(134),F12(134),FIA(134)
     $/ 1.66000000E+01, 4.52915468E+01,-1.834576476007E-03/
      DATA ETA(135),F12(135),FIA(135)
     $/ 1.68000000E+01, 4.61076365E+01,-1.788034097328E-03/
      DATA ETA(136),F12(136),FIA(136)
     $/ 1.70000000E+01, 4.69286280E+01,-1.739134736556E-03/
      DATA ETA(137),F12(137),FIA(137)
     $/ 1.72000000E+01, 4.77544832E+01,-1.709653214469E-03/
      DATA ETA(138),F12(138),FIA(138)
     $/ 1.74000000E+01, 4.85851870E+01,-1.661786300207E-03/
      DATA ETA(139),F12(139),FIA(139)
     $/ 1.76000000E+01, 4.94207048E+01,-1.633815786441E-03/
      DATA ETA(140),F12(140),FIA(140)
     $/ 1.78000000E+01, 5.02610178E+01,-1.584346715242E-03/
      DATA ETA(141),F12(141),FIA(141)
     $/ 1.80000000E+01, 5.11060801E+01,-1.552916203350E-03/
      DATA ETA(142),F12(142),FIA(142)
     $/ 1.82000000E+01, 5.19558725E+01,-1.531255731001E-03/
      DATA ETA(143),F12(143),FIA(143)
     $/ 1.84000000E+01, 5.28103838E+01,-1.473403678003E-03/
      DATA ETA(144),F12(144),FIA(144)
     $/ 1.86000000E+01, 5.36695566E+01,-1.467474616122E-03/
      DATA ETA(145),F12(145),FIA(145)
     $/ 1.88000000E+01, 5.45334024E+01,-1.426838647689E-03/
      DATA ETA(146),F12(146),FIA(146)
     $/ 1.90000000E+01, 5.54018793E+01,-1.377498261708E-03/
      DATA ETA(147),F12(147),FIA(147)
     $/ 1.92000000E+01, 5.62749338E+01,-1.383183816947E-03/
      DATA ETA(148),F12(148),FIA(148)
     $/ 1.94000000E+01, 5.71525927E+01,-1.327524740849E-03/
      DATA ETA(149),F12(149),FIA(149)
     $/ 1.96000000E+01, 5.80347986E+01,-1.315576157266E-03/
      DATA ETA(150),F12(150),FIA(150)
     $/ 1.98000000E+01, 5.89215441E+01,-1.283351521326E-03/
      DATA ETA(151),F12(151),FIA(151)
     $/ 2.00000000E+01, 5.98127985E+01,-1.252219378154E-03/
      DATA ETA(152),F12(152),FIA(152)
     $/ 2.02000000E+01, 6.07085314E+01,-1.242160880550E-03/
      DATA ETA(153),F12(153),FIA(153)
     $/ 2.04000000E+01, 6.16087389E+01,-1.194983711474E-03/
      DATA ETA(154),F12(154),FIA(154)
     $/ 2.06000000E+01, 6.25133677E+01,-1.197567559790E-03/
      DATA ETA(155),F12(155),FIA(155)
     $/ 2.08000000E+01, 6.34224367E+01,-1.159493320310E-03/
      DATA ETA(156),F12(156),FIA(156)
     $/ 2.10000000E+01, 6.43359013E+01,-1.135131594020E-03/
      DATA ETA(157),F12(157),FIA(157)
     $/ 2.12000000E+01, 6.52537327E+01,-1.125982877246E-03/
      DATA ETA(158),F12(158),FIA(158)
     $/ 2.14000000E+01, 6.61759281E+01,-1.085955718327E-03/
      DATA ETA(159),F12(159),FIA(159)
     $/ 2.16000000E+01, 6.71024418E+01,-1.091435662477E-03/
      DATA ETA(160),F12(160),FIA(160)
     $/ 2.18000000E+01, 6.80332966E+01,-1.052359718121E-03/
      DATA ETA(161),F12(161),FIA(161)
     $/ 2.20000000E+01, 6.89684391E+01,-1.034522675630E-03/
      DATA ETA(162),F12(162),FIA(162)
     $/ 2.22000000E+01, 6.99078465E+01,-1.025329275552E-03/
      DATA ETA(163),F12(163),FIA(163)
     $/ 2.24000000E+01, 7.08515186E+01,-9.955965775667E-04/
      DATA ETA(164),F12(164),FIA(164)
     $/ 2.26000000E+01, 7.17994175E+01,-9.877545065490E-04/
      DATA ETA(165),F12(165),FIA(165)
     $/ 2.28000000E+01, 7.27515431E+01,-9.690227539283E-04/
      DATA ETA(166),F12(166),FIA(166)
     $/ 2.30000000E+01, 7.37078724E+01,-9.420188557201E-04/
      DATA ETA(167),F12(167),FIA(167)
     $/ 2.32000000E+01, 7.46683674E+01,-9.412102489942E-04/
      DATA ETA(168),F12(168),FIA(168)
     $/ 2.34000000E+01, 7.56330357E+01,-9.082185476753E-04/
      DATA ETA(169),F12(169),FIA(169)
     $/ 2.36000000E+01, 7.66018314E+01,-9.080175988176E-04/
      DATA ETA(170),F12(170),FIA(170)
     $/ 2.38000000E+01, 7.75747700E+01,-8.923514704106E-04/
      DATA ETA(171),F12(171),FIA(171)
     $/ 2.40000000E+01, 7.85518284E+01,-8.566171195270E-04/
      DATA ETA(172),F12(172),FIA(172)
     $/ 2.42000000E+01, 7.95329456E+01,-8.706095710928E-04/
      DATA ETA(173),F12(173),FIA(173)
     $/ 2.44000000E+01, 8.05181599E+01,-8.344314949249E-04/
      DATA ETA(174),F12(174),FIA(174)
     $/ 2.46000000E+01, 8.15074177E+01,-8.346049903250E-04/
      DATA ETA(175),F12(175),FIA(175)
     $/ 2.48000000E+01, 8.25007267E+01,-8.180174375162E-04/
      DATA ETA(176),F12(176),FIA(176)
     $/ 2.50000000E+01, 8.34980640E+01,-7.958462623517E-04/
      DATA ETA(177),F12(177),FIA(177)
     $/ 2.52000000E+01, 8.44993916E+01,-7.966215802492E-04/
      DATA ETA(178),F12(178),FIA(178)
     $/ 2.54000000E+01, 8.55047245E+01,-7.763841457450E-04/
      DATA ETA(179),F12(179),FIA(179)
     $/ 2.56000000E+01, 8.65140324E+01,-7.646832835317E-04/
      DATA ETA(180),F12(180),FIA(180)
     $/ 2.58000000E+01, 8.75272999E+01,-7.591181446306E-04/
      DATA ETA(181),F12(181),FIA(181)
     $/ 2.60000000E+01, 8.85445194E+01,-7.309023009314E-04/
      DATA ETA(182),F12(182),FIA(182)
     $/ 2.62000000E+01, 8.95656452E+01,-7.452741062596E-04/
      DATA ETA(183),F12(183),FIA(183)
     $/ 2.64000000E+01, 9.05907154E+01,-7.079420670479E-04/
      DATA ETA(184),F12(184),FIA(184)
     $/ 2.66000000E+01, 9.16196613E+01,-7.165112531131E-04/
      DATA ETA(185),F12(185),FIA(185)
     $/ 2.68000000E+01, 9.26525059E+01,-6.963409034097E-04/
      DATA ETA(186),F12(186),FIA(186)
     $/ 2.70000000E+01, 9.36892185E+01,-6.873393148553E-04/
      DATA ETA(187),F12(187),FIA(187)
     $/ 2.72000000E+01, 9.47297840E+01,-6.807802206753E-04/
      DATA ETA(188),F12(188),FIA(188)
     $/ 2.74000000E+01, 9.57741947E+01,-6.623413987856E-04/
      DATA ETA(189),F12(189),FIA(189)
     $/ 2.76000000E+01, 9.68224201E+01,-6.661242492164E-04/
      DATA ETA(190),F12(190),FIA(190)
     $/ 2.78000000E+01, 9.78744755E+01,-6.418176402524E-04/
      DATA ETA(191),F12(191),FIA(191)
     $/ 2.80000000E+01, 9.89303150E+01,-6.453620029091E-04/
      DATA ETA(192),F12(192),FIA(192)
     $/ 2.82000000E+01, 9.99899540E+01,-6.295277664925E-04/
      DATA ETA(193),F12(193),FIA(193)
     $/ 2.84000000E+01, 1.01053362E+02,-6.176912881028E-04/
      DATA ETA(194),F12(194),FIA(194)
     $/ 2.86000000E+01, 1.02120516E+02,-6.181417564908E-04/
      DATA ETA(195),F12(195),FIA(195)
     $/ 2.88000000E+01, 1.03191431E+02,-6.038839429712E-04/
      DATA ETA(196),F12(196),FIA(196)
     $/ 2.90000000E+01, 1.04266077E+02,-5.927623476398E-04/
      DATA ETA(197),F12(197),FIA(197)
     $/ 2.92000000E+01, 1.05344431E+02,-5.918987997347E-04/
      DATA ETA(198),F12(198),FIA(198)
     $/ 2.94000000E+01, 1.06426508E+02,-5.840707347615E-04/
      DATA ETA(199),F12(199),FIA(199)
     $/ 2.96000000E+01, 1.07512262E+02,-5.368048493408E-04/
      DATA ETA(200),F12(200),FIA(200)
     $/ 2.98000000E+01, 1.08601709E+02,-7.132081093603E-04/
      DATA ETA(201),F12(201),FIA(201)
     $/ 3.00000000E+01, 1.09694826E+02, 0.000000000000E+00/
C
C
C          If Y is between 4.0234D-5 and 109.695 then use spline table
      IF((Y.LE.109.695D0).AND.(Y.GE.4.0234D-5)) THEN
        CALL INTRP(Y,F1,F12,ETA,FIA,N,NLOW,NHIGH)
C          Else if Y is greater than 109.695 use degen. approximation
      ELSEIF(Y.GT.109.695D0) THEN
        X2=(1.5D0*Y)**(0.666666667D0)  
        X4=1.0D0/(X2*X2) 
        F1=X2*(1.0D0+(AI(1)+(AI(2)+AI(3)*X4)*X4)*X4)
C          Else if Y is less than 4.0234D-5 use nondegen. approximation
      ELSE
        F1=LOG(AI(7)*MAX(Y,1.0D-20)*(1.D0+(AI(4)+(AI(5)+AI(8)*Y)*Y)*Y))
      ENDIF
C
      FINV12=F1  
 999  RETURN
      END   
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FHALF
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    CALL LINE:    FHALF(Y)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       Ratio of 1/2th Fermi Integral to the -1/2th Fermi
C                  Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FHALF(Y)
      IMPLICIT NONE
C
      INTEGER N, NLOW, NHIGH
      PARAMETER(N=201)
      DOUBLE PRECISION Y, A(7), ETA(N), FR(N), FRA(N), TH, X2, F0, F1
C
C
C                Retain values between calls
      SAVE NLOW,NHIGH
C
      DATA A/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0,
     $ 6.41500299D-02,.4D0,1.32934039D0/
      DATA TH,NLOW,NHIGH/0.3333333333333D0,1,N/
C
C
C                  Cubic spline data
      DATA ETA(  1),FR(  1),FRA(  1)
     $/-1.00000000E+01, 5.00008047E-01, 0.000000000000E+00/
      DATA ETA(  2),FR(  2),FRA(  2)
     $/-9.80000000E+00, 5.00009809E-01, 1.160669214176E-05/
      DATA ETA(  3),FR(  3),FRA(  3)
     $/-9.60000000E+00, 5.00011969E-01, 1.327318909372E-05/
      DATA ETA(  4),FR(  4),FRA(  4)
     $/-9.40000000E+00, 5.00014639E-01, 1.170462092124E-05/
      DATA ETA(  5),FR(  5),FRA(  5)
     $/-9.20000000E+00, 5.00017836E-01, 1.896244399075E-05/
      DATA ETA(  6),FR(  6),FRA(  6)
     $/-9.00000000E+00, 5.00021791E-01, 2.619339170395E-05/
      DATA ETA(  7),FR(  7),FRA(  7)
     $/-8.80000000E+00, 5.00026692E-01, 1.825697414678E-05/
      DATA ETA(  8),FR(  8),FRA(  8)
     $/-8.60000000E+00, 5.00032512E-01, 3.857146281817E-05/
      DATA ETA(  9),FR(  9),FRA(  9)
     $/-8.40000000E+00, 5.00039730E-01, 3.706621670383E-05/
      DATA ETA( 10),FR( 10),FRA( 10)
     $/-8.20000000E+00, 5.00048524E-01, 4.966876577432E-05/
      DATA ETA( 11),FR( 11),FRA( 11)
     $/-8.00000000E+00, 5.00059289E-01, 5.986611168059E-05/
      DATA ETA( 12),FR( 12),FRA( 12)
     $/-7.80000000E+00, 5.00072436E-01, 6.815575901197E-05/
      DATA ETA( 13),FR( 13),FRA( 13)
     $/-7.60000000E+00, 5.00088423E-01, 9.353156612777E-05/
      DATA ETA( 14),FR( 14),FRA( 14)
     $/-7.40000000E+00, 5.00108035E-01, 1.015077511200E-04/
      DATA ETA( 15),FR( 15),FRA( 15)
     $/-7.20000000E+00, 5.00131901E-01, 1.384040424797E-04/
      DATA ETA( 16),FR( 16),FRA( 16)
     $/-7.00000000E+00, 5.00161156E-01, 1.532660352027E-04/
      DATA ETA( 17),FR( 17),FRA( 17)
     $/-6.80000000E+00, 5.00196763E-01, 2.014816386790E-04/
      DATA ETA( 18),FR( 18),FRA( 18)
     $/-6.60000000E+00, 5.00240334E-01, 2.353001371116E-04/
      DATA ETA( 19),FR( 19),FRA( 19)
     $/-6.40000000E+00, 5.00293492E-01, 2.953124554528E-04/
      DATA ETA( 20),FR( 20),FRA( 20)
     $/-6.20000000E+00, 5.00358457E-01, 3.546103132566E-04/
      DATA ETA( 21),FR( 21),FRA( 21)
     $/-6.00000000E+00, 5.00437760E-01, 4.368757029944E-04/
      DATA ETA( 22),FR( 22),FRA( 22)
     $/-5.80000000E+00, 5.00534601E-01, 5.285417635068E-04/
      DATA ETA( 23),FR( 23),FRA( 23)
     $/-5.60000000E+00, 5.00652775E-01, 6.489897161003E-04/
      DATA ETA( 24),FR( 24),FRA( 24)
     $/-5.40000000E+00, 5.00797055E-01, 7.913279951832E-04/
      DATA ETA( 25),FR( 25),FRA( 25)
     $/-5.20000000E+00, 5.00973167E-01, 9.605790184029E-04/
      DATA ETA( 26),FR( 26),FRA( 26)
     $/-5.00000000E+00, 5.01188039E-01, 1.180342230219E-03/
      DATA ETA( 27),FR( 27),FRA( 27)
     $/-4.80000000E+00, 5.01450297E-01, 1.425959520518E-03/
      DATA ETA( 28),FR( 28),FRA( 28)
     $/-4.60000000E+00, 5.01770117E-01, 1.750084178541E-03/
      DATA ETA( 29),FR( 29),FRA( 29)
     $/-4.40000000E+00, 5.02160252E-01, 2.120996973156E-03/
      DATA ETA( 30),FR( 30),FRA( 30)
     $/-4.20000000E+00, 5.02635834E-01, 2.582978496851E-03/
      DATA ETA( 31),FR( 31),FRA( 31)
     $/-4.00000000E+00, 5.03215376E-01, 3.140940050287E-03/
      DATA ETA( 32),FR( 32),FRA( 32)
     $/-3.80000000E+00, 5.03921295E-01, 3.809961012388E-03/
      DATA ETA( 33),FR( 33),FRA( 33)
     $/-3.60000000E+00, 5.04780575E-01, 4.623412578668E-03/
      DATA ETA( 34),FR( 34),FRA( 34)
     $/-3.40000000E+00, 5.05825854E-01, 5.596028708118E-03/
      DATA ETA( 35),FR( 35),FRA( 35)
     $/-3.20000000E+00, 5.07096250E-01, 6.760222613033E-03/
      DATA ETA( 36),FR( 36),FRA( 36)
     $/-3.00000000E+00, 5.08638542E-01, 8.147420929689E-03/
      DATA ETA( 37),FR( 37),FRA( 37)
     $/-2.80000000E+00, 5.10508521E-01, 9.803125317419E-03/
      DATA ETA( 38),FR( 38),FRA( 38)
     $/-2.60000000E+00, 5.12772505E-01, 1.174090923418E-02/
      DATA ETA( 39),FR( 39),FRA( 39)
     $/-2.40000000E+00, 5.15508351E-01, 1.401245248421E-02/
      DATA ETA( 40),FR( 40),FRA( 40)
     $/-2.20000000E+00, 5.18807170E-01, 1.665511284146E-02/
      DATA ETA( 41),FR( 41),FRA( 41)
     $/-2.00000000E+00, 5.22774696E-01, 1.967341586786E-02/
      DATA ETA( 42),FR( 42),FRA( 42)
     $/-1.80000000E+00, 5.27531885E-01, 2.310044392523E-02/
      DATA ETA( 43),FR( 43),FRA( 43)
     $/-1.60000000E+00, 5.33215712E-01, 2.692057340338E-02/
      DATA ETA( 44),FR( 44),FRA( 44)
     $/-1.40000000E+00, 5.39978788E-01, 3.110469072503E-02/
      DATA ETA( 45),FR( 45),FRA( 45)
     $/-1.20000000E+00, 5.47988069E-01, 3.559136204354E-02/
      DATA ETA( 46),FR( 46),FRA( 46)
     $/-1.00000000E+00, 5.57422393E-01, 4.028629703731E-02/
      DATA ETA( 47),FR( 47),FRA( 47)
     $/-8.00000000E-01, 5.68468657E-01, 4.505445242759E-02/
      DATA ETA( 48),FR( 48),FRA( 48)
     $/-6.00000000E-01, 5.81316396E-01, 4.971717382570E-02/
      DATA ETA( 49),FR( 49),FRA( 49)
     $/-4.00000000E-01, 5.96151066E-01, 5.411642817285E-02/
      DATA ETA( 50),FR( 50),FRA( 50)
     $/-2.00000000E-01, 6.13146980E-01, 5.800386974609E-02/
      DATA ETA( 51),FR( 51),FRA( 51)
     $/ 0.00000000E+00, 6.32458839E-01, 6.125961642856E-02/
      DATA ETA( 52),FR( 52),FRA( 52)
     $/ 2.00000000E-01, 6.54215519E-01, 6.368105824286E-02/
      DATA ETA( 53),FR( 53),FRA( 53)
     $/ 4.00000000E-01, 6.78513328E-01, 6.518526072255E-02/
      DATA ETA( 54),FR( 54),FRA( 54)
     $/ 6.00000000E-01, 7.05412292E-01, 6.575130141701E-02/
      DATA ETA( 55),FR( 55),FRA( 55)
     $/ 8.00000000E-01, 7.34934891E-01, 6.535476854033E-02/
      DATA ETA( 56),FR( 56),FRA( 56)
     $/ 1.00000000E+00, 7.67065889E-01, 6.408943410837E-02/
      DATA ETA( 57),FR( 57),FRA( 57)
     $/ 1.20000000E+00, 8.01755749E-01, 6.211687075523E-02/
      DATA ETA( 58),FR( 58),FRA( 58)
     $/ 1.40000000E+00, 8.38925803E-01, 5.947211293743E-02/
      DATA ETA( 59),FR( 59),FRA( 59)
     $/ 1.60000000E+00, 8.78471892E-01, 5.639992416069E-02/
      DATA ETA( 60),FR( 60),FRA( 60)
     $/ 1.80000000E+00, 9.20271987E-01, 5.302910851018E-02/
      DATA ETA( 61),FR( 61),FRA( 61)
     $/ 2.00000000E+00, 9.64191670E-01, 4.942180766085E-02/
      DATA ETA( 62),FR( 62),FRA( 62)
     $/ 2.20000000E+00, 1.01008812E+00, 4.579899792645E-02/
      DATA ETA( 63),FR( 63),FRA( 63)
     $/ 2.40000000E+00, 1.05781651E+00, 4.217259322424E-02/
      DATA ETA( 64),FR( 64),FRA( 64)
     $/ 2.60000000E+00, 1.10723239E+00, 3.863512430722E-02/
      DATA ETA( 65),FR( 65),FRA( 65)
     $/ 2.80000000E+00, 1.15819499E+00, 3.529307034029E-02/
      DATA ETA( 66),FR( 66),FRA( 66)
     $/ 3.00000000E+00, 1.21057018E+00, 3.208270748942E-02/
      DATA ETA( 67),FR( 67),FRA( 67)
     $/ 3.20000000E+00, 1.26423024E+00, 2.910594156410E-02/
      DATA ETA( 68),FR( 68),FRA( 68)
     $/ 3.40000000E+00, 1.31905603E+00, 2.635426039667E-02/
      DATA ETA( 69),FR( 69),FRA( 69)
     $/ 3.60000000E+00, 1.37493737E+00, 2.380801215631E-02/
      DATA ETA( 70),FR( 70),FRA( 70)
     $/ 3.80000000E+00, 1.43177253E+00, 2.148777375912E-02/
      DATA ETA( 71),FR( 71),FRA( 71)
     $/ 4.00000000E+00, 1.48946852E+00, 1.936381862584E-02/
      DATA ETA( 72),FR( 72),FRA( 72)
     $/ 4.20000000E+00, 1.54794050E+00, 1.745753296313E-02/
      DATA ETA( 73),FR( 73),FRA( 73)
     $/ 4.40000000E+00, 1.60711200E+00, 1.573191622765E-02/
      DATA ETA( 74),FR( 74),FRA( 74)
     $/ 4.60000000E+00, 1.66691366E+00, 1.413981428697E-02/
      DATA ETA( 75),FR( 75),FRA( 75)
     $/ 4.80000000E+00, 1.72728245E+00, 1.277958760511E-02/
      DATA ETA( 76),FR( 76),FRA( 76)
     $/ 5.00000000E+00, 1.78816286E+00, 1.148276396969E-02/
      DATA ETA( 77),FR( 77),FRA( 77)
     $/ 5.20000000E+00, 1.84950378E+00, 1.036724436250E-02/
      DATA ETA( 78),FR( 78),FRA( 78)
     $/ 5.40000000E+00, 1.91126027E+00, 9.382681356530E-03/
      DATA ETA( 79),FR( 79),FRA( 79)
     $/ 5.60000000E+00, 1.97339221E+00, 8.419828946720E-03/
      DATA ETA( 80),FR( 80),FRA( 80)
     $/ 5.80000000E+00, 2.03586235E+00, 7.668390275724E-03/
      DATA ETA( 81),FR( 81),FRA( 81)
     $/ 6.00000000E+00, 2.09863908E+00, 6.894960215311E-03/
      DATA ETA( 82),FR( 82),FRA( 82)
     $/ 6.20000000E+00, 2.16169262E+00, 6.274318363239E-03/
      DATA ETA( 83),FR( 83),FRA( 83)
     $/ 6.40000000E+00, 2.22499741E+00, 5.694183220240E-03/
      DATA ETA( 84),FR( 84),FRA( 84)
     $/ 6.60000000E+00, 2.28853030E+00, 5.163412244369E-03/
      DATA ETA( 85),FR( 85),FRA( 85)
     $/ 6.80000000E+00, 2.35227011E+00, 4.691055087970E-03/
      DATA ETA( 86),FR( 86),FRA( 86)
     $/ 7.00000000E+00, 2.41619827E+00, 4.323617268449E-03/
      DATA ETA( 87),FR( 87),FRA( 87)
     $/ 7.20000000E+00, 2.48029883E+00, 3.875965713102E-03/
      DATA ETA( 88),FR( 88),FRA( 88)
     $/ 7.40000000E+00, 2.54455567E+00, 3.614743308902E-03/
      DATA ETA( 89),FR( 89),FRA( 89)
     $/ 7.60000000E+00, 2.60895650E+00, 3.262132860466E-03/
      DATA ETA( 90),FR( 90),FRA( 90)
     $/ 7.80000000E+00, 2.67348849E+00, 3.011856988495E-03/
      DATA ETA( 91),FR( 91),FRA( 91)
     $/ 8.00000000E+00, 2.73814091E+00, 2.753550229065E-03/
      DATA ETA( 92),FR( 92),FRA( 92)
     $/ 8.20000000E+00, 2.80290404E+00, 2.582380523235E-03/
      DATA ETA( 93),FR( 93),FRA( 93)
     $/ 8.40000000E+00, 2.86776968E+00, 2.292083501014E-03/
      DATA ETA( 94),FR( 94),FRA( 94)
     $/ 8.60000000E+00, 2.93272843E+00, 2.215304435771E-03/
      DATA ETA( 95),FR( 95),FRA( 95)
     $/ 8.80000000E+00, 2.99777480E+00, 1.991168746407E-03/
      DATA ETA( 96),FR( 96),FRA( 96)
     $/ 9.00000000E+00, 3.06290146E+00, 1.862794325246E-03/
      DATA ETA( 97),FR( 97),FRA( 97)
     $/ 9.20000000E+00, 3.12810241E+00, 1.701422191793E-03/
      DATA ETA( 98),FR( 98),FRA( 98)
     $/ 9.40000000E+00, 3.19337213E+00, 1.646713595709E-03/
      DATA ETA( 99),FR( 99),FRA( 99)
     $/ 9.60000000E+00, 3.25870691E+00, 1.470861367159E-03/
      DATA ETA(100),FR(100),FRA(100)
     $/ 9.80000000E+00, 3.32410133E+00, 1.415936859359E-03/
      DATA ETA(101),FR(101),FRA(101)
     $/ 1.00000000E+01, 3.38955181E+00, 1.273721641867E-03/
      DATA ETA(102),FR(102),FRA(102)
     $/ 1.02000000E+01, 3.45505393E+00, 1.235581482298E-03/
      DATA ETA(103),FR(103),FRA(103)
     $/ 1.04000000E+01, 3.52060505E+00, 1.134274717011E-03/
      DATA ETA(104),FR(104),FRA(104)
     $/ 1.06000000E+01, 3.58620184E+00, 1.077069548628E-03/
      DATA ETA(105),FR(105),FRA(105)
     $/ 1.08000000E+01, 3.65184173E+00, 1.022991960613E-03/
      DATA ETA(106),FR(106),FRA(106)
     $/ 1.10000000E+01, 3.71752230E+00, 9.341272710622E-04/
      DATA ETA(107),FR(107),FRA(107)
     $/ 1.12000000E+01, 3.78324041E+00, 8.692110436287E-04/
      DATA ETA(108),FR(108),FRA(108)
     $/ 1.14000000E+01, 3.84899382E+00, 8.850650648729E-04/
      DATA ETA(109),FR(109),FRA(109)
     $/ 1.16000000E+01, 3.91478159E+00, 7.459549148300E-04/
      DATA ETA(110),FR(110),FRA(110)
     $/ 1.18000000E+01, 3.98060041E+00, 7.864194213587E-04/
      DATA ETA(111),FR(111),FRA(111)
     $/ 1.20000000E+01, 4.04644996E+00, 7.192672054582E-04/
      DATA ETA(112),FR(112),FRA(112)
     $/ 1.22000000E+01, 4.11232811E+00, 6.253942766733E-04/
      DATA ETA(113),FR(113),FRA(113)
     $/ 1.24000000E+01, 4.17823227E+00, 6.811530337942E-04/
      DATA ETA(114),FR(114),FRA(114)
     $/ 1.26000000E+01, 4.24416274E+00, 5.969604189414E-04/
      DATA ETA(115),FR(115),FRA(115)
     $/ 1.28000000E+01, 4.31011748E+00, 5.709059055339E-04/
      DATA ETA(116),FR(116),FRA(116)
     $/ 1.30000000E+01, 4.37609516E+00, 5.615397884449E-04/
      DATA ETA(117),FR(117),FRA(117)
     $/ 1.32000000E+01, 4.44209499E+00, 5.044814122883E-04/
      DATA ETA(118),FR(118),FRA(118)
     $/ 1.34000000E+01, 4.50811538E+00, 5.045695992625E-04/
      DATA ETA(119),FR(119),FRA(119)
     $/ 1.36000000E+01, 4.57415563E+00, 4.552324352849E-04/
      DATA ETA(120),FR(120),FRA(120)
     $/ 1.38000000E+01, 4.64021477E+00, 5.087014064024E-04/
      DATA ETA(121),FR(121),FRA(121)
     $/ 1.40000000E+01, 4.70629289E+00, 3.581955061927E-04/
      DATA ETA(122),FR(122),FRA(122)
     $/ 1.42000000E+01, 4.77238707E+00, 4.664301602782E-04/
      DATA ETA(123),FR(123),FRA(123)
     $/ 1.44000000E+01, 4.83849842E+00, 3.514279312654E-04/
      DATA ETA(124),FR(124),FRA(124)
     $/ 1.46000000E+01, 4.90462510E+00, 4.274941163720E-04/
      DATA ETA(125),FR(125),FRA(125)
     $/ 1.48000000E+01, 4.97076776E+00, 3.365813609877E-04/
      DATA ETA(126),FR(126),FRA(126)
     $/ 1.50000000E+01, 5.03692462E+00, 3.549707958916E-04/
      DATA ETA(127),FR(127),FRA(127)
     $/ 1.52000000E+01, 5.10309513E+00, 2.910214326093E-04/
      DATA ETA(128),FR(128),FRA(128)
     $/ 1.54000000E+01, 5.16927849E+00, 4.096124892769E-04/
      DATA ETA(129),FR(129),FRA(129)
     $/ 1.56000000E+01, 5.23547613E+00, 2.118142438439E-04/
      DATA ETA(130),FR(130),FRA(130)
     $/ 1.58000000E+01, 5.30168451E+00, 3.539701323030E-04/
      DATA ETA(131),FR(131),FRA(131)
     $/ 1.60000000E+01, 5.36790514E+00, 2.099303736157E-04/
      DATA ETA(132),FR(132),FRA(132)
     $/ 1.62000000E+01, 5.43413621E+00, 3.726817539688E-04/
      DATA ETA(133),FR(133),FRA(133)
     $/ 1.64000000E+01, 5.50037998E+00, 2.035344437886E-04/
      DATA ETA(134),FR(134),FRA(134)
     $/ 1.66000000E+01, 5.56663346E+00, 2.697503364877E-04/
      DATA ETA(135),FR(135),FRA(135)
     $/ 1.68000000E+01, 5.63289722E+00, 2.603197709474E-04/
      DATA ETA(136),FR(136),FRA(136)
     $/ 1.70000000E+01, 5.69917089E+00, 1.752530698713E-04/
      DATA ETA(137),FR(137),FRA(137)
     $/ 1.72000000E+01, 5.76545287E+00, 2.848424287755E-04/
      DATA ETA(138),FR(138),FRA(138)
     $/ 1.74000000E+01, 5.83174481E+00, 1.787223011498E-04/
      DATA ETA(139),FR(139),FRA(139)
     $/ 1.76000000E+01, 5.89804520E+00, 2.684866184639E-04/
      DATA ETA(140),FR(140),FRA(140)
     $/ 1.78000000E+01, 5.96435493E+00, 1.488933760423E-04/
      DATA ETA(141),FR(141),FRA(141)
     $/ 1.80000000E+01, 6.03067153E+00, 1.660321432524E-04/
      DATA ETA(142),FR(142),FRA(142)
     $/ 1.82000000E+01, 6.09699556E+00, 3.007372019870E-04/
      DATA ETA(143),FR(143),FRA(143)
     $/ 1.84000000E+01, 6.16332905E+00, 5.039185382223E-05/
      DATA ETA(144),FR(144),FRA(144)
     $/ 1.86000000E+01, 6.22966780E+00, 2.880169444751E-04/
      DATA ETA(145),FR(145),FRA(145)
     $/ 1.88000000E+01, 6.29601555E+00, 1.459654809555E-04/
      DATA ETA(146),FR(146),FRA(146)
     $/ 1.90000000E+01, 6.36236965E+00, 8.022559364040E-05/
      DATA ETA(147),FR(147),FRA(147)
     $/ 1.92000000E+01, 6.42872883E+00, 2.967073127773E-04/
      DATA ETA(148),FR(148),FRA(148)
     $/ 1.94000000E+01, 6.49509677E+00, 4.547040400385E-05/
      DATA ETA(149),FR(149),FRA(149)
     $/ 1.96000000E+01, 6.56146941E+00, 2.275888933824E-04/
      DATA ETA(150),FR(150),FRA(150)
     $/ 1.98000000E+01, 6.62784906E+00, 9.594149535815E-05/
      DATA ETA(151),FR(151),FRA(151)
     $/ 2.00000000E+01, 6.69423382E+00, 1.539878947448E-04/
      DATA ETA(152),FR(152),FRA(152)
     $/ 2.02000000E+01, 6.76062461E+00, 1.936619445003E-04/
      DATA ETA(153),FR(153),FRA(153)
     $/ 2.04000000E+01, 6.82702164E+00, 6.161813072681E-06/
      DATA ETA(154),FR(154),FRA(154)
     $/ 2.06000000E+01, 6.89342166E+00, 2.302557839252E-04/
      DATA ETA(155),FR(155),FRA(155)
     $/ 2.08000000E+01, 6.95982860E+00, 1.118128460033E-04/
      DATA ETA(156),FR(156),FRA(156)
     $/ 2.10000000E+01, 7.02624073E+00, 1.013932536890E-04/
      DATA ETA(157),FR(157),FRA(157)
     $/ 2.12000000E+01, 7.09265730E+00, 1.470875659220E-04/
      DATA ETA(158),FR(158),FRA(158)
     $/ 2.14000000E+01, 7.15907881E+00, 5.141533931037E-05/
      DATA ETA(159),FR(159),FRA(159)
     $/ 2.16000000E+01, 7.22550396E+00, 1.950632132465E-04/
      DATA ETA(160),FR(160),FRA(160)
     $/ 2.18000000E+01, 7.29193509E+00, 6.331854420200E-05/
      DATA ETA(161),FR(161),FRA(161)
     $/ 2.20000000E+01, 7.35836991E+00, 1.054146921341E-04/
      DATA ETA(162),FR(162),FRA(162)
     $/ 2.22000000E+01, 7.42480867E+00, 1.072633188645E-04/
      DATA ETA(163),FR(163),FRA(163)
     $/ 2.24000000E+01, 7.49125168E+00, 1.027449618964E-04/
      DATA ETA(164),FR(164),FRA(164)
     $/ 2.26000000E+01, 7.55769894E+00, 1.183006063983E-04/
      DATA ETA(165),FR(165),FRA(165)
     $/ 2.28000000E+01, 7.62415059E+00, 8.253924684004E-05/
      DATA ETA(166),FR(166),FRA(166)
     $/ 2.30000000E+01, 7.69060579E+00, 8.451551828560E-05/
      DATA ETA(167),FR(167),FRA(167)
     $/ 2.32000000E+01, 7.75706437E+00, 8.685419149568E-05/
      DATA ETA(168),FR(168),FRA(168)
     $/ 2.34000000E+01, 7.82352627E+00, 6.485775084105E-05/
      DATA ETA(169),FR(169),FRA(169)
     $/ 2.36000000E+01, 7.88999116E+00, 1.043033206627E-04/
      DATA ETA(170),FR(170),FRA(170)
     $/ 2.38000000E+01, 7.95646040E+00, 1.690790286878E-04/
      DATA ETA(171),FR(171),FRA(171)
     $/ 2.40000000E+01, 8.02293437E+00,-7.229290158689E-05/
      DATA ETA(172),FR(172),FRA(172)
     $/ 2.42000000E+01, 8.08940876E+00, 1.853851230301E-04/
      DATA ETA(173),FR(173),FRA(173)
     $/ 2.44000000E+01, 8.15588780E+00, 2.682857763816E-05/
      DATA ETA(174),FR(174),FRA(174)
     $/ 2.46000000E+01, 8.22236941E+00, 9.315857978144E-05/
      DATA ETA(175),FR(175),FRA(175)
     $/ 2.48000000E+01, 8.28885454E+00, 1.279132160594E-04/
      DATA ETA(176),FR(176),FRA(176)
     $/ 2.50000000E+01, 8.35534335E+00,-5.151623010689E-05/
      DATA ETA(177),FR(177),FRA(177)
     $/ 2.52000000E+01, 8.42183289E+00, 1.870374120185E-04/
      DATA ETA(178),FR(178),FRA(178)
     $/ 2.54000000E+01, 8.48832695E+00,-1.961965051994E-05/
      DATA ETA(179),FR(179),FRA(179)
     $/ 2.56000000E+01, 8.55482254E+00, 1.229752535397E-04/
      DATA ETA(180),FR(180),FRA(180)
     $/ 2.58000000E+01, 8.62132190E+00, 9.217292731240E-05/
      DATA ETA(181),FR(181),FRA(181)
     $/ 2.60000000E+01, 8.68782416E+00,-5.687752100593E-05/
      DATA ETA(182),FR(182),FRA(182)
     $/ 2.62000000E+01, 8.75432688E+00, 2.031850501482E-04/
      DATA ETA(183),FR(183),FRA(183)
     $/ 2.64000000E+01, 8.82083439E+00,-3.579723391714E-05/
      DATA ETA(184),FR(184),FRA(184)
     $/ 2.66000000E+01, 8.88734279E+00, 7.422832091008E-05/
      DATA ETA(185),FR(185),FRA(185)
     $/ 2.68000000E+01, 8.95385391E+00, 1.447580583760E-04/
      DATA ETA(186),FR(186),FRA(186)
     $/ 2.70000000E+01, 9.02036919E+00,-2.735083889056E-05/
      DATA ETA(187),FR(187),FRA(187)
     $/ 2.72000000E+01, 9.08688502E+00, 4.542808910866E-05/
      DATA ETA(188),FR(188),FRA(188)
     $/ 2.74000000E+01, 9.15340256E+00, 1.035480068597E-04/
      DATA ETA(189),FR(189),FRA(189)
     $/ 2.76000000E+01, 9.21992316E+00,-6.679452567842E-07/
      DATA ETA(190),FR(190),FRA(190)
     $/ 2.78000000E+01, 9.28644514E+00, 1.046883666256E-04/
      DATA ETA(191),FR(191),FRA(191)
     $/ 2.80000000E+01, 9.35297006E+00, 2.438337814088E-05/
      DATA ETA(192),FR(192),FRA(192)
     $/ 2.82000000E+01, 9.41949693E+00, 8.980690836273E-05/
      DATA ETA(193),FR(193),FRA(193)
     $/ 2.84000000E+01, 9.48602609E+00,-3.931880271592E-05/
      DATA ETA(194),FR(194),FRA(194)
     $/ 2.86000000E+01, 9.55255558E+00, 1.158160177556E-04/
      DATA ETA(195),FR(195),FRA(195)
     $/ 2.88000000E+01, 9.61908818E+00, 4.256278618984E-05/
      DATA ETA(196),FR(196),FRA(196)
     $/ 2.90000000E+01, 9.68562279E+00, 1.531962767634E-05/
      DATA ETA(197),FR(197),FRA(197)
     $/ 2.92000000E+01, 9.75215876E+00, 1.018906708333E-04/
      DATA ETA(198),FR(198),FRA(198)
     $/ 2.94000000E+01, 9.81869734E+00,-3.387023032781E-05/
      DATA ETA(199),FR(199),FRA(199)
     $/ 2.96000000E+01, 9.88523639E+00, 1.063057625796E-04/
      DATA ETA(200),FR(200),FRA(200)
     $/ 2.98000000E+01, 9.95177816E+00, 1.564710994435E-05/
      DATA ETA(201),FR(201),FRA(201)
     $/ 3.00000000E+01, 1.00183211E+01, 0.000000000000E+00/
C
C          If Y is between -10 and 30 then use the spline tables
      IF((Y.GE.-10.0D0).AND.(Y.LE.30.0D0)) THEN
        CALL INTRP(Y,F1,ETA,FR,FRA,N,NLOW,NHIGH)
C          Else if Y is greater than 30 use degenerate approximation
      ELSEIF(Y.GT.30.0D0) THEN
        X2=Y**(-2)
        F1=Y*TH*(1.0D0+(0.2D0*A(1)+(0.6D0*A(2)+1.4D0*X2*A(3))*X2)*X2)
     >  /(1.0D0-(0.2D0*TH*A(1)-(A(2)-4.2D0*X2*A(3))*X2)*X2)  
      ELSE
C          Else if Y is less than -10 use nondegenerate approximation
        F0=EXP(Y)   
        F1=(1.0D0-(2.0D0*A(4)-(3*A(5)-0.125D0*F0)*F0)*F0)/
     >  (2.0D0-(8.0D0*A(4)-(18.0D0*A(5)-F0)*F0)*F0)
      ENDIF
C
      FHALF=F1  
 999  RETURN
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
      SUBROUTINE INTRP(X,Y,XTABLE,YTABLE,SC,NPTS,NLOW,NHIGH)
C
      IMPLICIT NONE
C
      INTEGER NPTS
      DOUBLE PRECISION X, Y, XTABLE(NPTS), YTABLE(NPTS), SC(NPTS)
C
C
      DOUBLE PRECISION C1, C2, S1, S2, DELTX, SIXTH
      PARAMETER(SIXTH=1.0D0/6.0D0)
      INTEGER NLOW, NHIGH, NMID
C
      NHIGH = MAX(MIN(NPTS-1,NHIGH),2)
      NLOW = MAX(MIN(NPTS-1,NLOW),2)
C
C             If old values of bracketing pointers work then use them...
      IF((NHIGH-NLOW).EQ.1) THEN
        IF((XTABLE(NHIGH).GT.X).AND.(XTABLE(NLOW).LE.X)) THEN
          GOTO 20
C                 Try the next lower pair
        ELSEIF((XTABLE(NHIGH-1).GT.X).AND.(XTABLE(NLOW-1).LE.X)) THEN
          NHIGH = MAX(NHIGH-1,1)
          NLOW = MAX(NLOW-1,1)
          GOTO 20
C                 Try the next higher pair
        ELSEIF((XTABLE(NHIGH+1).GT.X).AND.(XTABLE(NLOW+1).LE.X)) THEN
          NHIGH = MIN(NHIGH+1,NPTS)
          NLOW = MIN(NLOW+1,NPTS)
          GOTO 20
        ENDIF          
      ENDIF
C
C       Otherwise find the ordinates that bracket the X value by
C       bisecting the table...
C
C                         Set the high & low pointers to their extrema
      NLOW = 1
      NHIGH = NPTS
C
 10   NMID = (NHIGH+NLOW)/2
      IF(X.GT.XTABLE(NMID)) THEN
        NLOW = NMID
      ELSE
        NHIGH = NMID
      ENDIF
      IF((NHIGH-NLOW).NE.1) GOTO 10 
C
C
C                  Now that the bracketing ordinates have been found,
C                  interpolate to get the Y value
C
 20   CONTINUE
      IF((NHIGH.GT.NPTS).OR.(NHIGH.LT.1).OR.
     *   (NLOW.GE.NHIGH).OR.(NLOW.LT.1) ) THEN
        WRITE(*,*) ' INTRP: NH,NL = ',NHIGH,NLOW,NPTS
      ENDIF
C                  Distance between nearest X entries
      DELTX = XTABLE(NHIGH)-XTABLE(NLOW)
C
C                  Now evaluate the polynomial...
C
C                        Calculate linear interp. coefficients
      C1 = (XTABLE(NHIGH)-X)/DELTX
      C2 = (X-XTABLE(NLOW))/DELTX
C
C                        Calculate spline interp. coefficients
      S1 = (C1**3)-C1
      S2 = (C2**3)-C2
C
C                        Calculate spline interp. value
      Y = C1*YTABLE(NLOW)+C2*YTABLE(NHIGH)+
     1  (S1*SC(NLOW)+S2*SC(NHIGH))*(DELTX**2)*SIXTH
C
C                        If interpolated value of Y is not monotonic
C                        then linearly interpolate
      IF((Y.LT.YTABLE(NLOW)).OR.(Y.GT.YTABLE(NHIGH))) THEN
        Y = C1*YTABLE(NLOW)+C2*YTABLE(NHIGH)
      ENDIF
C
      RETURN
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C                  This file contains the complete electron EOS
C                  package el_eos_pak.f.  Also needed is the include
C                  file el_eos.inc.  This complete set of routines
C                  is included in the LS EOS code as of version 2.7
C
C                  The routines contained in this package are:
C
C                  EL_EOS, EL_REL, EL_IMT, IMTRULE, GL16, & F2
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_EOS.FOR
C
C***********************************************************************
C
C    MODULE:       EL_EOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, 
C                  Laboratory for Computational Astrophysics
C                  Dept. of Astronomy & NCSA
C                  University of Illinois at Urbana-Champaign
C
C    EMAIL:        dswesty@ncsa.uiuc.edu
C
C    DATE:         (original) 2/12/91
C                  (v2.7)     9/15/95
C
C    PURPOSE:      The electron and photon equation of state
C
C
C    CALL LINE:    CALL EL_EOS(T,YE,BRYDNS)
C
C    INPUTS:       T = TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      NONE
C
C
C 
C    INCLUDE FILES:  EL_EOS.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EL_EOS(T,YE,BRYDNS)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION T, YE, BRYDNS
C
      INCLUDE 'el_eos.inc'
C
C
C
C                           Plancks constant & speed of light
      DOUBLE PRECISION HBAR, C
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23)
C
C                           Pi and 1/3
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23
      PARAMETER(PI=3.1415927, PI2=9.8696044)
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667)
C
C
      IF( ((BRYDNS*YE).GT.1.0D-4).OR.(T.GT.5.0D0) ) THEN
C                      Assume that the electrons are relativistic
        CALL EL_REL(T,YE,BRYDNS)
      ELSE
C                      Otherwise call the Gauss-Laguerre version
        CALL EL_IMT(T,YE,BRYDNS)
C                      If the stuff is really degenerate call the
C                      relativistic version anyway
        IF((MUSUBE/T).GT.2.0D3) THEN
          CALL EL_REL(T,YE,BRYDNS)
        ENDIF
      ENDIF
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
 999  RETURN
C
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_REL.FOR
C
C***********************************************************************
C
C    MODULE:       EL_REL
C
C    TYPE:         SUBROUTINE
C
C    AUTHOR:       F. DOUGLAS SWESTY
C                  Laboratory for Computational Astrophysics
C                  Dept. of Astronomy & NCSA
C                  University of Illinois at Urbana-Champaign
C
C    EMAIL:        dswesty@ncsa.uiuc.edu
C
C    DATE:         2/12/91
C
C
C    PURPOSE:      The relativistic electron and photon EOS
C
C
C    CALL LINE:    CALL EL_REL(T,YE,BRYDNS)
C
C    INPUTS:       T = TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      NONE
C
C
C 
C    INCLUDE FILES:  EL_REL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EL_REL(T,YE,BRYDNS)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION T, YE, BRYDNS
C
      INCLUDE 'el_eos.inc'
C
C
C
C                           Plancks constant & speed of light
      DOUBLE PRECISION HBAR, C
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23)
C
C                           Pi and 1/3
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23
      PARAMETER(PI=3.1415927, PI2=9.8696044)
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667)
C
C                           2nd Fermi integral
      DOUBLE PRECISION F_2
C
C                           Positron degeneracy parameter
      DOUBLE PRECISION ELPETA
C
C
C                    Leptons
C
C                    Electron number density
      NSUBE = BRYDNS*YE
C
C                    Coefficants for chemical potential
C                    and thermodynamics quantities
      QSUBE = 1.0/( 3.0*(PI**2)*((HBAR*C)**3) )
C
      ACOEF = 0.5*NSUBE/QSUBE
C
      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0)**0.5
C
      DBDT = (PI**6)*(T**5)/(9.0*BCOEF)
C
      CCOEF = (ACOEF+BCOEF)**OVR3
C
C
C                    Electron chemical potential
      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF
C
C                    Positron degeneracy parameter
      ELPETA = -MUSUBE/T
C
C                    Positron number density
      NEPLUS =  3.0*QSUBE*(T**3)*F_2(ELPETA)
C
C
C                    Electron pressure for rel. case
      EPRESS = 0.25*QSUBE*(MUSUBE**4+2.0*(PI*T*MUSUBE)**2+
     1 7.0*((PI*T)**4)/15.0)
C
C
C                    Electron internal energy per baryon
      EU = 0.75*QSUBE*(MUSUBE**4+2.0*(PI*MUSUBE*T)**2+
     1 7.0*((PI*T)**4)/15.0)/BRYDNS
C
C
C                    Electron free energy per baryon
      FSUBE = ((MUSUBE*NSUBE)-EPRESS)/BRYDNS
C
C                    Electron entropy per baryon
      ES = QSUBE*(((PI*MUSUBE)**2)*T+7.0*(PI**4)*(T**3)/
     1 15.0)/BRYDNS
C
C                    Photons
C
C                    Photon pressure
      PPRESS = (PI**2)*(T**4)/(45.0*((HBAR*C)**3))
C                    Photon entropy per baryon
      PS = 4.0*PPRESS/(T*BRYDNS)
C
C                    Photon internal energy per baryon
      PU = 3.0*PPRESS/BRYDNS
C
C                    Photon free energy per baryon
      PF = PU-T*PS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C                    Derivatives of chem. potential w.r.t. T,
C                    BRYDNS, YE
C
      DEMUDT = DBDT/(3.0*CCOEF**2)-OVR23*(PI**2)*T/CCOEF+
     1         DBDT*((PI*T)**2)/(9.0*CCOEF**4)
C
      DEMUDN = (YE*PI2*(HBAR*C)**3)/(MUSUBE**2+OVR3*PI2*T**2)
C
      DEMUDY = BRYDNS*DEMUDN/YE
C
C
C                    Derivatives of pressure w.r.t. BRYDNS,YE,T
C
      DEPDN = BRYDNS*YE*DEMUDN
C
      DEPDY = BRYDNS*DEPDN/YE
C
      DEPDT = BRYDNS*(ES+YE*DEMUDT)
C
C
C                    Derivatives of entropy w.r.t. T,BRYDNS,YE
C
      DESDT = ES/T+OVR23*(7.0*PI2*(T**2)/15.0+MUSUBE*T*DEMUDT)/
     1        (BRYDNS*(HBAR*C)**3)
C
      DESDN = -1.0*DEPDT/(BRYDNS**2)
C
      DESDY = 2.0*T*QSUBE*PI2*MUSUBE*DEMUDY/BRYDNS
C
C
C                    Derivatives of internal energy w.r.t.
C                    T,BRYDNS,YE
      DEUDT = T*DESDT
C
      DEUDN = (YE*(MUSUBE-T*DEMUDT)-EU)/BRYDNS
C
      DEUDY = 3.0*QSUBE*((MUSUBE**3)+PI2*(T**2)*MUSUBE)*
     1        DEMUDY/BRYDNS
C
C
C                               Photons
C
C                    Derivatives of photon pressure
      DPPDN = 0.0
      DPPDT = BRYDNS*PS
      DPPDY = 0.0
C
C                    Derivatives of photon entropy
      DPSDN = -PS/BRYDNS
      DPSDT = 3.0*PS/T
      DPSDY = 0.0
C
C                    Derivatives of internal energy
      DPUDN = -0.75*T*PS/BRYDNS
      DPUDT = 3.0*PS
      DPUDY = 0.0
C
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
 999  RETURN
C
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_IMT.FOR
C
C***********************************************************************
C
C    MODULE:       EL_IMT
C
C    TYPE:         SUBROUTINE
C
C    AUTHOR:       F. DOUGLAS SWESTY
C                  Laboratory for Computational Astrophysics
C                  Dept. of Astronomy & NCSA
C                  University of Illinois at Urbana-Champaign
C
C    EMAIL:        dswesty@ncsa.uiuc.edu
C
C    DATE:         9/15/95
C
C    PURPOSE:      The elctron and photon equation of state
C
C
C    CALL LINE:    CALL EL_IMT(T,YE,BRYDNS)
C
C    INPUTS:       T = TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      NONE
C
C
C 
C    INCLUDE FILES:  EL_EOS.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EL_IMT(T,YE,BRYDNS)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION T, YE, BRYDNS
C
      INCLUDE 'el_eos.inc'
C
C
C
C                           Plancks constant & speed of light
      DOUBLE PRECISION HBAR, C
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23)
      DOUBLE PRECISION MASS_E
      PARAMETER(MASS_E=0.51100D0)
C
C                           Pi and 1/3
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23
      PARAMETER(PI=3.1415927, PI2=9.8696044)
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667)
C
C                           2nd Fermi integral
      DOUBLE PRECISION F_2
C                           Positron degeneracy parameter
      DOUBLE PRECISION ELPETA
C                           Multiplicative factors in the G-L
C                           quadrature
      DOUBLE PRECISION TFAC, DTFDT, EFAC, PFAC
C                           Electron eta and beta
      DOUBLE PRECISION ETA_E, BETA, DBTDT, DETDT, DETDNE, D_ETA
C                           Number density and it's derivatives
      DOUBLE PRECISION NE_CHK, NEN, DNENDB, DNENDE
      DOUBLE PRECISION NEM, DNEMDB, DNEMDE
      DOUBLE PRECISION NEP, DNEPDB, DNEPDE
C                           Energy density & derivatives
      DOUBLE PRECISION E_ENG, DEDB, DEDE
      DOUBLE PRECISION EM_ENG, DEMDB, DEMDE
      DOUBLE PRECISION EP_ENG, DEPDB, DEPDE
C                           Pressure & derivatives
      DOUBLE PRECISION E_PR, DPDB, DPDE
      DOUBLE PRECISION PM, DPMDB, DPMDE
      DOUBLE PRECISION PP, DPPDB, DPPDE
C
C                           Loop variables
      INTEGER J, MAXITR
      PARAMETER(MAXITR=30)
C
C                           Convergence tolerance for N-R iteration
      DOUBLE PRECISION EPSIL
      PARAMETER(EPSIL=1.0D-10)
C
      DOUBLE PRECISION X_LO, X_UP, ETA
C
C-----------------------------------------------------------------------
C
C
C
C
C                    Electron number density
      NSUBE = BRYDNS*YE
C
C                    Coefficants for chemical potential
C                    and thermodynamics quantities
      QSUBE = 1.0D0/( 3.0D0*(PI**2)*((HBAR*C)**3) )
C
      ACOEF = 0.5D0*NSUBE/QSUBE
C
      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0D0)**0.5D0
C
      DBDT = (PI**6)*(T**5)/(9.0D0*BCOEF)
C
      CCOEF = (ACOEF+BCOEF)**OVR3
C
C
C                    Electron chemical potential
      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF
C
C
C
C                            Initial guess at electron degeneracy 
C                            parameter
      ETA_E = MUSUBE/T
C
      BETA = MASS_E/T
      DBTDT = -BETA/T
C                            Multiplicative factor in Fermi integrals
      TFAC = (T**3)/((PI**2)*(HBAR*C)**3)
      DTFDT = 3.0D0*TFAC/T
      PFAC = TFAC*T*OVR3
      EFAC = TFAC*T/BRYDNS
C
C                            Loop until the N-R iteration for the
C                            chemical potential converges
      DO 20 J=1,MAXITR,1
C
C                            Zero the electron number density integral
        NEM = 0.0D0
        DNEMDE = 0.0D0
        DNEMDB = 0.0D0
        NEP = 0.0D0
        DNEPDE = 0.0D0
        DNEPDB = 0.0D0
C
C                            Zero the electron energy density integral
        EM_ENG = 0.0D0
        DEMDE = 0.0D0
        DEMDB = 0.0D0
        EP_ENG = 0.0D0
        DEPDE = 0.0D0
        DEPDB = 0.0D0
C
C                            Zero the electron pressure integral
        PM = 0.0D0
        DPMDE = 0.0D0
        DPMDB = 0.0D0
        PP = 0.0D0
        DPPDE = 0.0D0
        DPPDB = 0.0D0
C
C                            Zero the positron number density
        NEPLUS = 0.0D0
C
C       ----------------------------------------------------------------
C       |  Do the quadrature for electrons first                       |
C       ----------------------------------------------------------------
        ETA = ETA_E
        IF(ETA.GT.5.0D-1) THEN
C       |                    Do the IMT Quadrature in two parts:       |
C       |                    (0,ETA) & (ETA,ETA+60)                    |
          X_LO = 0.0D0
          X_UP = ETA
          CALL IMTRULE(ETA,BETA,X_LO,X_UP,NEM,DNEMDB,DNEMDE,
     *                 EM_ENG,DEMDB,DEMDE,PM,DPMDB,DPMDE)
          X_LO = ETA
          X_UP = ETA+60.0D0
          CALL IMTRULE(ETA,BETA,X_LO,X_UP,NEM,DNEMDB,DNEMDE,
     *                 EM_ENG,DEMDB,DEMDE,PM,DPMDB,DPMDE)
C       |                                                              |
        ELSE                                                           !
C       ----------------------------------------------------------------
C       |                    Do the IMT Quadrature in 1 part           |
C       ----------------------------------------------------------------
          X_LO = 0.0D0
          X_UP = 60.0D0
          CALL IMTRULE(ETA,BETA,X_LO,X_UP,NEM,DNEMDB,DNEMDE,
     *                 EM_ENG,DEMDB,DEMDE,PM,DPMDB,DPMDE)
cc          CALL GL16(ETA,BETA,X_LO,X_UP,NEM,DNEMDB,DNEMDE,
cc     *                 EM_ENG,DEMDB,DEMDE,PM,DPMDB,DPMDE)
C       |                                                              |
        ENDIF
C       ----------------------------------------------------------------
C       |               Now do the quadrature for positrons            |
C       ----------------------------------------------------------------
        ETA = -ETA_E
        IF(ETA.GT.5.0D-1) THEN
C       |                    Do the IMT Quadrature in two parts:       |
C       |                    (0,ETA) & (ETA,ETA+60)                    |
C       |                                                              |
          X_LO = 0.0D0
          X_UP = ETA
          CALL IMTRULE(ETA,BETA,X_LO,X_UP,NEP,DNEPDB,DNEPDE,
     *                 EP_ENG,DEPDB,DEPDE,PP,DPPDB,DPPDE)
          X_LO = ETA
          X_UP = ETA+60.0D0
          CALL IMTRULE(ETA,BETA,X_LO,X_UP,NEP,DNEPDB,DNEPDE,
     *                 EP_ENG,DEPDB,DEPDE,PP,DPPDB,DPPDE)
C       |                                                              |
        ELSE
          X_LO = 0.0D0
          X_UP = 60.0D0
          CALL IMTRULE(ETA,BETA,X_LO,X_UP,NEP,DNEPDB,DNEPDE,
     *                 EP_ENG,DEPDB,DEPDE,PP,DPPDB,DPPDE)
cc          CALL GL16(ETA,BETA,X_LO,X_UP,NEP,DNEPDB,DNEPDE,
cc     *                 EP_ENG,DEPDB,DEPDE,PP,DPPDB,DPPDE)
        ENDIF
C       |                                                              |
C       |                                                              |
C       |                                                              |
C       |                                                              |
C       ----------------------------------------------------------------
C
C                    Now sum the electron & positron contributions
        NEN = NEM-NEP
        DNENDE = DNEMDE+DNEPDE
        DNENDB = DNEMDB-DNEPDB
        E_ENG = EM_ENG+EP_ENG
        DEDE = DEMDE-DEPDE
        DEDB = DEMDB+DEPDB
        E_PR = PM+PP
        DPDE = DPMDE-DPPDE
        DPDB = DPMDB+DPPDB
        NEPLUS = NEP
C
C                            Multiply by the temperature factor to get
C                            number density
        NE_CHK = TFAC*NEN
C                            Multiply by the temperature factor to get
C                            number density
        NEPLUS = TFAC*NEPLUS
C
C                            Calculate the new change in ETA
        D_ETA = -(NE_CHK-NSUBE)/(DNENDE*TFAC)
C
C                   If we've met the convergence criterion...
        IF(DABS(D_ETA).LT.(EPSIL*ETA_E)) THEN
C                          Then break out the chemical potential
C                          loop
          GOTO 30
        ELSE
C                          Otherwise update the chemical potential
          ETA_E = ETA_E+D_ETA
        ENDIF
C
 20   CONTINUE
C
C                                      If we reached this point the
C                                      N-R iteration didn't converge
      WRITE(*,*) ' EL_EOS: N-R iteration didnt converge!'
      WRITE(*,*) ' EL_EOS: TRY = ',T,BRYDNS,YE
      STOP
C
C               If we reached this point then the N-R iteration has
C               converged.
 30   CONTINUE
C
C                     Is the result consistent with electron number
C                     density?
      IF(DABS((NE_CHK-NSUBE)/NSUBE).GT.1.0D-5) THEN
        WRITE(*,*) ' EL_EOS: Gaussian quadrature converged badly!'
        STOP
      ENDIF
C
C             Calculate thermodynamic quantities...
C
C                          Electron chemical potential
      MUSUBE = T*ETA_E
C                          Electron pressure
      EPRESS = PFAC*E_PR
C                          Electron internal energy per baryon
      EU = EFAC*E_ENG
C                          Electron free energy per baryon
      FSUBE = YE*MUSUBE-EPRESS/BRYDNS
C                          Electron entropy per baryon
      ES = (EU-FSUBE)/T
C
C                          Derivative of the electron eta w.r.t. T
      DETDT = -(DTFDT*NEN+TFAC*DNENDB*DBTDT)/(TFAC*DNENDE)
C                          Derivative of the electron eta w.r.t. NSUBE
      DETDNE = 1.0D0/(TFAC*DNENDE)
C
C                    Derivatives of chem. potential w.r.t. T,
C                    BRYDNS, YE
      DEMUDT = T*DETDT+MUSUBE/T
      DEMUDN = YE*T*DETDNE
      DEMUDY = BRYDNS*T*DETDNE
C
C
C                    Derivatives of pressure w.r.t. BRYDNS,YE,T
      DEPDN = YE*PFAC*DPDE*DETDNE
      DEPDY = BRYDNS*PFAC*DPDE*DETDNE
      DEPDT = (4.0D0*PFAC/T)*E_PR+PFAC*(DPDE*DETDT+DPDB*DBTDT)
C
C                    Derivatives of internal energy w.r.t.
C                    T,BRYDNS,YE
      DEUDT = (4.0D0*EFAC/T)*E_ENG+EFAC*(DEDB*DBTDT+DEDE*DETDT)
      DEUDN = (-1.0D0*EFAC/BRYDNS)*E_ENG+YE*EFAC*DEDE*DETDNE
      DEUDY = BRYDNS*EFAC*DEDE*DETDNE
C
C                    Derivatives of entropy w.r.t. T,BRYDNS,YE
      DESDT = -ES/T+(DEUDT-YE*DEMUDT+DEPDT/BRYDNS)/T
      DESDN = (DEUDN-YE*DEMUDN+DEPDN/BRYDNS-EPRESS/(BRYDNS**2))/T
      DESDY = (DEUDY-MUSUBE-YE*DEMUDY+DEPDY/BRYDNS)/T
C
C
C
C                    Photon pressure
      PPRESS = (PI**2)*(T**4)/(45.0*((HBAR*C)**3))
C                    Photon entropy per baryon
      PS = 4.0*PPRESS/(T*BRYDNS)
C
C                    Photon internal energy per baryon
      PU = 3.0*PPRESS/BRYDNS
C
C                    Photon free energy per baryon
      PF = PU-T*PS
C
C
C                    Derivatives of photon pressure
      DPPDN = 0.0
      DPPDT = BRYDNS*PS
      DPPDY = 0.0
C
C                    Derivatives of photon entropy
      DPSDN = -PS/BRYDNS
      DPSDT = 3.0*PS/T
      DPSDY = 0.0
C
C                    Derivatives of internal energy
      DPUDN = -0.75*T*PS/BRYDNS
      DPUDT = 3.0*PS
      DPUDY = 0.0
C
 999  RETURN
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       IMTRULE
C
C    TYPE:         SUBROUTINE
C
C    AUTHOR:       F. DOUGLAS SWESTY
C                  Laboratory for Computational Astrophysics
C                  Dept. of Astronomy & NCSA
C                  University of Illinois at Urbana-Champaign
C
C    EMAIL:        dswesty@ncsa.uiuc.edu
C
C    DATE:         8/25/95
C
C
C    PURPOSE:      Do the generalized Fermi integrals via the IMT rule
C
C
C    CALL LINE:          SUBROUTINE IMTRULE(ETA,BETA,X_LO,X_UP,
C                       * N,DNDB,DNDE,E,DEDB,DEDN,P,DPDB,DPDE)
C
C    INPUTS:       
C
C    OUTPUTS:      NONE
C
C
C 
C    INCLUDE FILES:  NONE
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE IMTRULE(ETA,BETA,X_LO,X_UP,
     *                   N,DNDB,DNDE,E,DEDB,DEDE,P,DPDB,DPDE)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ETA, BETA, X_LO, X_UP
      DOUBLE PRECISION N, DNDB, DNDE
      DOUBLE PRECISION E, DEDB, DEDE
      DOUBLE PRECISION P, DPDB, DPDE
C
      INTEGER I
      DOUBLE PRECISION FE, DFEDB, DFEDE
      DOUBLE PRECISION XROOT, XN1, DXN1DB
      DOUBLE PRECISION XN2, DXN2DB, XN3, DXN3DB
C
      INTEGER N_IMT
      PARAMETER(N_IMT=64)
      DOUBLE PRECISION W_IMT(1:N_IMT-1), X_IMT(1:N_IMT-1)
C
C                           Integration weights & limits
      DOUBLE PRECISION WGHT, X
C
C-----------------------------------------------------------------------
C         63 pt IMT quadrature abscissae & weights
C-----------------------------------------------------------------------
C
      DATA X_IMT(  1),W_IMT(  1 ) /
     *    1.9570031734585200D-30 ,  8.2607320509999998D-27 /
      DATA X_IMT(  2),W_IMT(  2 ) /
     *    5.9130760189713709D-16 ,  6.4169110630000000D-13 /
      DATA X_IMT(  3),W_IMT(  3 ) /
     *    5.4675241449266301D-11 ,  2.7067693377000001D-08 /
      DATA X_IMT(  4),W_IMT(  4 ) /
     *    1.9306978596659198D-08 ,  5.5092726333000005D-06 /
      DATA X_IMT(  5),W_IMT(  5 ) /
     *    7.1039233799069106D-07 ,  1.3273484743000000D-04 /
      DATA X_IMT(  6),W_IMT(  6 ) /
     *    8.2971031730510406D-06 ,  1.0999107621000000D-03 /
      DATA X_IMT(  7),W_IMT(  7 ) /
     *    4.9827396171457499D-05 ,  4.9514468819000000D-03 /
      DATA X_IMT(  8),W_IMT(  8 ) /
     *    1.9629223097988500D-04 ,  1.5218120421000000D-02 /
      DATA X_IMT(  9),W_IMT(  9 ) /
     *    5.8152639478483402D-04 ,  3.6255696685999997D-02 /
      DATA X_IMT( 10),W_IMT( 10 ) /
     *    1.4073702294191501D-03 ,  7.2251913311999996D-02 /
      DATA X_IMT( 11),W_IMT( 11 ) /
     *    2.9342961722595702D-03 ,  1.2642067111999999D-01 /
      DATA X_IMT( 12),W_IMT( 12 ) /
     *    5.4624206134241907D-03 ,  2.0058566844999998D-01 /
      DATA X_IMT( 13),W_IMT( 13 ) /
     *    9.3088894298794592D-03 ,  2.9511499967999999D-01 /
      DATA X_IMT( 14),W_IMT( 14 ) /
     *    1.4786180477344599D-02 ,  4.0908191706999997D-01 /
      DATA X_IMT( 15),W_IMT( 15 ) /
     *    2.2183887527440597D-02 ,  5.4053256120000004D-01 /
      DATA X_IMT( 16),W_IMT( 16 ) /
     *    3.1754957727637700D-02 ,  6.8677770086000001D-01 /
      DATA X_IMT( 17),W_IMT( 17 ) /
     *    4.3706349739020198D-02 ,  8.4466152386999993D-01 /
      DATA X_IMT( 18),W_IMT( 18 ) /
     *    5.8193563828888392D-02 ,  1.0107865713999999D+00 /
      DATA X_IMT( 19),W_IMT( 19 ) /
     *    7.5318304334592898D-02 ,  1.1816897903000001D+00 /
      DATA X_IMT( 20),W_IMT( 20 ) /
     *    9.5128529563062905D-02 ,  1.3539730362000000D+00 /
      DATA X_IMT( 21),W_IMT( 21 ) /
     *    1.1762022895325901D-01 ,  1.5243949451000001D+00 /
      DATA X_IMT( 22),W_IMT( 22 ) /
     *    1.4274038487381499D-01 ,  1.6899319837000000D+00 /
      DATA X_IMT( 23),W_IMT( 23 ) /
     *    1.7039069590573200D-01 ,  1.8478160150000000D+00 /
      DATA X_IMT( 24),W_IMT( 24 ) /
     *    2.0043174541744399D-01 ,  1.9955546823999999D+00 /
      DATA X_IMT( 25),W_IMT( 25 ) /
     *    2.3268738852436499D-01 ,  2.1309397372999999D+00 /
      DATA X_IMT( 26),W_IMT( 26 ) /
     *    2.6694920180823400D-01 ,  2.2520473296999999D+00 /
      DATA X_IMT( 27),W_IMT( 27 ) /
     *    3.0298089532632300D-01 ,  2.3572333323999999D+00 /
      DATA X_IMT( 28),W_IMT( 28 ) /
     *    3.4052262810720502D-01 ,  2.4451259898000002D+00 /
      DATA X_IMT( 29),W_IMT( 29 ) /
     *    3.7929519918613303D-01 ,  2.5146175750999999D+00 /
      DATA X_IMT( 30),W_IMT( 30 ) /
     *    4.1900410866587701D-01 ,  2.5648562654000000D+00 /
      DATA X_IMT( 31),W_IMT( 31 ) /
     *    4.5934349927577800D-01 ,  2.5952390874000000D+00 /
      DATA X_IMT( 32),W_IMT( 32 ) /
     *    5.0000000000000000D-01 ,  2.6054065145000003D+00 /
      DATA X_IMT( 33),W_IMT( 33 ) /
     *    5.4065650072422200D-01 ,  2.5952390874000000D+00 /
      DATA X_IMT( 34),W_IMT( 34 ) /
     *    5.8099589133412299D-01 ,  2.5648562654000000D+00 /
      DATA X_IMT( 35),W_IMT( 35 ) /
     *    6.2070480081386692D-01 ,  2.5146175750999999D+00 /
      DATA X_IMT( 36),W_IMT( 36 ) /
     *    6.5947737189279498D-01 ,  2.4451259898000002D+00 /
      DATA X_IMT( 37),W_IMT( 37 ) /
     *    6.9701910467367700D-01 ,  2.3572333323999999D+00 /
      DATA X_IMT( 38),W_IMT( 38 ) /
     *    7.3305079819176600D-01 ,  2.2520473296999999D+00 /
      DATA X_IMT( 39),W_IMT( 39 ) /
     *    7.6731261147563501D-01 ,  2.1309397372999999D+00 /
      DATA X_IMT( 40),W_IMT( 40 ) /
     *    7.9956825458255598D-01 ,  1.9955546823999999D+00 /
      DATA X_IMT( 41),W_IMT( 41 ) /
     *    8.2960930409426803D-01 ,  1.8478160150000000D+00 /
      DATA X_IMT( 42),W_IMT( 42 ) /
     *    8.5725961512618498D-01 ,  1.6899319837000000D+00 /
      DATA X_IMT( 43),W_IMT( 43 ) /
     *    8.8237977104674103D-01 ,  1.5243949451000001D+00 /
      DATA X_IMT( 44),W_IMT( 44 ) /
     *    9.0487147043693705D-01 ,  1.3539730362000000D+00 /
      DATA X_IMT( 45),W_IMT( 45 ) /
     *    9.2468169566540714D-01 ,  1.1816897903000001D+00 /
      DATA X_IMT( 46),W_IMT( 46 ) /
     *    9.4180643617111159D-01 ,  1.0107865713999999D+00 /
      DATA X_IMT( 47),W_IMT( 47 ) /
     *    9.5629365026097979D-01 ,  8.4466152386999993D-01 /
      DATA X_IMT( 48),W_IMT( 48 ) /
     *    9.6824504227236230D-01 ,  6.8677770086000001D-01 /
      DATA X_IMT( 49),W_IMT( 49 ) /
     *    9.7781611247255940D-01 ,  5.4053256120000004D-01 /
      DATA X_IMT( 50),W_IMT( 50 ) /
     *    9.8521381952265541D-01 ,  4.0908191706999997D-01 /
      DATA X_IMT( 51),W_IMT( 51 ) /
     *    9.9069111057012049D-01 ,  2.9511499967999999D-01 /
      DATA X_IMT( 52),W_IMT( 52 ) /
     *    9.9453757938657583D-01 ,  2.0058566844999998D-01 /
      DATA X_IMT( 53),W_IMT( 53 ) /
     *    9.9706570382774040D-01 ,  1.2642067111999999D-01 /
      DATA X_IMT( 54),W_IMT( 54 ) /
     *    9.9859262977058083D-01 ,  7.2251913311999996D-02 /
      DATA X_IMT( 55),W_IMT( 55 ) /
     *    9.9941847360521519D-01 ,  3.6255696685999997D-02 /
      DATA X_IMT( 56),W_IMT( 56 ) /
     *    9.9980370776902017D-01 ,  1.5218120421000000D-02 /
      DATA X_IMT( 57),W_IMT( 57 ) /
     *    9.9995017260382857D-01 ,  4.9514468819000000D-03 /
      DATA X_IMT( 58),W_IMT( 58 ) /
     *    9.9999170289682693D-01 ,  1.0999107621000000D-03 /
      DATA X_IMT( 59),W_IMT( 59 ) /
     *    9.9999928960766205D-01 ,  1.3273484743000000D-04 /
      DATA X_IMT( 60),W_IMT( 60 ) /
     *    9.9999998069302143D-01 ,  5.5092726333000005D-06 /
      DATA X_IMT( 61),W_IMT( 61 ) /
     *    9.9999999994532474D-01 ,  2.7067693377000001D-08 /
      DATA X_IMT( 62),W_IMT( 62 ) /
     *    9.9999999999999944D-01 ,  6.4169110630000000D-13 /
      DATA X_IMT( 63),W_IMT( 63 ) /
     *    1.0000000000000000D+00 ,  8.2607320509999998D-27 /
C
C-----------------------------------------------------------------------
C
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
      DO I=1,N_IMT-1,1                                         !   |
        X = (X_UP-X_LO)*X_IMT(I)+X_LO                          !   |
        WGHT = W_IMT(I)*(X_UP-X_LO)/(DBLE(N_IMT))        !   |
C       |                    Electron occupation numbers       !   |
        FE = 1.0D0/(DEXP(X+BETA-ETA)+1.0D0)                         !   |
C       |                    Derivatives w.r.t. beta           !   |
        DFEDB = -FE*(1.0D0-FE)                                 !   |
C       |                    Derivatives w.r.t. Eta            !   |
        DFEDE = FE*(1.0D0-FE)                                  !   |
C       |                    Number density & derivatives      !   |
        XROOT = DSQRT(X*(X+2.0D0*BETA))                        !   |
        XN1 = XROOT*(X+BETA)                                   !   |
        DXN1DB = XROOT+X*(X+BETA)/XROOT                        !   |
C       |                    Sum for electron # integral       !   |
        N = N+WGHT*XN1*FE                                      !   |
        DNDB = DNDB+WGHT*(DXN1DB*FE+XN1*DFEDB)                 !   |
        DNDE = DNDE+WGHT*XN1*DFEDE                             !   |
C       |                    Energy integral & derivatives     !   |
        XN2 = XROOT*(X+BETA)**2                                !   |
        DXN2DB = 2.0D0*XROOT*(X+BETA)+X*((X+BETA)**2)/XROOT    !   |
        E = E+WGHT*XN2*FE                                      !   |
        DEDB = DEDB+WGHT*(DXN2DB*FE+XN2*DFEDB)                 !   |
        DEDE = DEDE+WGHT*XN2*DFEDE                             !   |
C       |                    Pressure integral & derivatives   !   |
        XN3 = XROOT**3                                         !   |
        DXN3DB = 3.0D0*XROOT*X                                 !   |
        P = P+WGHT*XN3*FE                                      !   |
        DPDB = DPDB+WGHT*(DXN3DB*FE+XN3*DFEDB)                 !   |
        DPDE = DPDE+WGHT*XN3*DFEDE                             !   |
      ENDDO                                                    !   |
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
C
 999  RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       GL16
C
C    TYPE:         SUBROUTINE
C
C    AUTHOR:       F. DOUGLAS SWESTY
C                  Laboratory for Computational Astrophysics
C                  Dept. of Astronomy & NCSA
C                  University of Illinois at Urbana-Champaign
C
C    EMAIL:        dswesty@ncsa.uiuc.edu
C
C    DATE:         8/25/95
C
C    PURPOSE:      Do the generalized Fermi integrals via 16 point
C                  Gauss--Laguerre quadrature
C
C
C    CALL LINE:          SUBROUTINE GL16(ETA,BETA,X_LO,X_UP,
C                       * N,DNDB,DNDE,E,DEDB,DEDN,P,DPDB,DPDE)
C
C    INPUTS:       
C
C    OUTPUTS:      NONE
C
C
C 
C    INCLUDE FILES:  NONE
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE GL16(ETA,BETA,X_LO,X_UP,
     *                   N,DNDB,DNDE,E,DEDB,DEDE,P,DPDB,DPDE)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ETA, BETA, X_LO, X_UP
      DOUBLE PRECISION N, DNDB, DNDE
      DOUBLE PRECISION E, DEDB, DEDE
      DOUBLE PRECISION P, DPDB, DPDE
C
      INTEGER I
      DOUBLE PRECISION FE, DFEDB, DFEDE
      DOUBLE PRECISION XROOT, XN1, DXN1DB
      DOUBLE PRECISION XN2, DXN2DB, XN3, DXN3DB
C
      INTEGER N_IMT
      PARAMETER(N_IMT=64)
      DOUBLE PRECISION W_IMT(1:N_IMT-1), X_IMT(1:N_IMT-1)
C
C                           Integration weights & limits
      DOUBLE PRECISION WGHT, X
C

      INTEGER N_GL16
      PARAMETER(N_GL16=16)
C
      DOUBLE PRECISION W_GL16(N_GL16), X_GL16(N_GL16)
C
      DATA X_GL16  /.087649410479D00
     *,.46269632892D00,.11410577748D01,.21292836451D01,.34370866339D01
     *,.50780186145D01,.70703385350D01,.94383143364D01,.12214223369D02
     *,.15441527369D02,.19180156857D02,.23515905694D02,.28578729743D02
     *,.34583398702D02,.41940452648D02,.51701160340D02/
      DATA W_GL16 /.22503631486D00
     *,.52583605276D00,.83196139169D00,.11460992410D01,.14717513170D01
     *,.18131346874D01,.21755175197D01,.25657627502D01,.29932150864D01
     *,.34712344831D01,.40200440864D01,.46725166077D01,.54874206580D01
     *,.65853612333D01,.82763579844D01,.11824277552D02/
C
C
C-----------------------------------------------------------------------
C
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
      DO I=1,N_GL16,1                                          !   |
        X = X_GL16(I)                                          !   |
        WGHT = W_GL16(I)                                       !   |
C       |                    Electron occupation numbers       !   |
        FE = 1.0D0/(DEXP(X+BETA-ETA)+1.0D0)                    !   |
C       |                    Derivatives w.r.t. beta           !   |
        DFEDB = -FE*(1.0D0-FE)                                 !   |
C       |                    Derivatives w.r.t. Eta            !   |
        DFEDE = FE*(1.0D0-FE)                                  !   |
C       |                    Number density & derivatives      !   |
        XROOT = DSQRT(X*(X+2.0D0*BETA))                        !   |
        XN1 = XROOT*(X+BETA)                                   !   |
        DXN1DB = XROOT+X*(X+BETA)/XROOT                        !   |
C       |                    Sum for electron # integral       !   |
        N = N+WGHT*XN1*FE                                      !   |
        DNDB = DNDB+WGHT*(DXN1DB*FE+XN1*DFEDB)                 !   |
        DNDE = DNDE+WGHT*XN1*DFEDE                             !   |
C       |                    Energy integral & derivatives     !   |
        XN2 = XROOT*(X+BETA)**2                                !   |
        DXN2DB = 2.0D0*XROOT*(X+BETA)+X*((X+BETA)**2)/XROOT    !   |
        E = E+WGHT*XN2*FE                                      !   |
        DEDB = DEDB+WGHT*(DXN2DB*FE+XN2*DFEDB)                 !   |
        DEDE = DEDE+WGHT*XN2*DFEDE                             !   |
C       |                    Pressure integral & derivatives   !   |
        XN3 = XROOT**3                                         !   |
        DXN3DB = 3.0D0*XROOT*X                                 !   |
        P = P+WGHT*XN3*FE                                      !   |
        DPDB = DPDB+WGHT*(DXN3DB*FE+XN3*DFEDB)                 !   |
        DPDE = DPDE+WGHT*XN3*DFEDE                             !   |
      ENDDO                                                    !   |
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
C
 999  RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY
C                  Laboratory for Computational Astrophysics
C                  Dept. of Astronomy & NCSA
C                  University of Illinois at Urbana-Champaign
C
C    EMAIL:        dswesty@ncsa.uiuc.edu
C    DATE:         12/16/91
C
C                  Version 2: 1/24/93
C
C    CALL LINE:    F_2(Y)      (2nd Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       2nd Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_2(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y, YEXP
      IF(Y.GT.3.0D0) THEN
        WRITE(*,*) ' F_2(Y) FAILS FOR Y .GT. 3; Y =',Y
        STOP
      ENDIF
C
C                       Note: This approximation is based on the
C                       Bludman & Van Riper approximation (see
C                       Ap. J. Vol. 212 page 866-867 (1977))
C                       equation (3.6)
C
      IF(Y.LT.-1.0D0) THEN
        YEXP = EXP(Y)
        F_2 = 2.0*YEXP*(1.0-0.125*YEXP+0.037037*(YEXP**2))
      ELSE
        F_2 = 1.803D0+1.645D0*Y+0.6931D0*(Y**2)+0.1666667D0*(Y**3)+
     1        2.0833333D-2*(Y**4)-3.4722D-4*(Y**6)
      ENDIF
 999  RETURN
      END   
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       MATLUD
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         6/16/94
C
C    PURPOSE:      LU decomposes an NxN matrix
C
C    NOTE:         This subroutine employs Crout's algorithm with
C                  implicit pivoting as described in "Numerical
C                  Recipes", by Press, Flannery, Teukolsky, and
C                  Vetterling, (pub. Cambridge Univ. Press) first
C                  edition, whose authors desrve the credit for
C                  this algorithm.
C
C    CALL LINE:    CALL MATLUD(A,LU,N,IPVT)
C
C    INPUTS:       A = NxN array to be LU decomposed  (D)
C                  N = Size of A (I)
C
C    OUTPUTS:      LU = Array containing LU decomposition of A (D)
C                  IPVT = Vector of pivot indices (I)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE MATLUD(A,LU,N,IPVT)
C
      IMPLICIT NONE
C
C                 Parameters
      INTEGER N
      INTEGER IPVT(N)
      DOUBLE PRECISION A(N,N), LU(N,N)
C
C                 Local variables
C
C                 Loop variables
      INTEGER I, J, L
C                 A small floating point value to prevent division
C                 by zero
      DOUBLE PRECISION SMALLV
      PARAMETER(SMALLV=1.0D-20)
C
C                 Value & original row index or largest L value
      DOUBLE PRECISION E_MAX, TSTVAL
      INTEGER ROWPVT
C
C                 A scratch variable to use when swapping rows
      DOUBLE PRECISION SCRVAR
C
C                 An array for the largest value of each row of A
      INTEGER NMAX
      PARAMETER(NMAX=100)
      DOUBLE PRECISION MAXVAL(NMAX)
C-----------------------------------------------------------------------
C
C                 Find the maximum absolute value in each row of A
      DO 20 I=1,N,1
        MAXVAL(I) = -1.0D0
        DO 10 J=1,N,1
          MAXVAl(I)  = DMAX1(SMALLV,DABS(A(I,J)),MAXVAL(I))
          LU(I,J) = A(I,J) 
10     CONTINUE
 20   CONTINUE
C
C                 Now employ Crout's algorithm with implicit pivoting
      DO 90 J=1,N,1
C
C                 Calculate column J, U matrix elements
        DO 40 I=1,J-1,1
          DO 30 L=1,I-1,1
            LU(I,J) = LU(I,J)-LU(I,L)*LU(L,J)
 30       CONTINUE
 40     CONTINUE
C
C                 Calculate column J, L matrix elements.  The element is
C                 scaled by the largest element of the original matrix.
C                 Also, the column from the diagonal down is searched
C                 for the largest element.
        E_MAX = -1.0D0
        DO 60 I=J,N,1
          DO 50 L=1,J-1,1
            LU(I,J) = LU(I,J)-LU(I,L)*LU(L,J)
 50       CONTINUE
          TSTVAL = DABS( LU(I,J) )/MAXVAL(I)
          IF(TSTVAL.GT.E_MAX) THEN
            E_MAX=TSTVAL
            ROWPVT = I
          ENDIF
 60     CONTINUE
C                    Keep track of which row was pivoted into row J
       IPVT(J) = ROWPVT
C
C
C
C                 If the original diagonal element wasn't the
C                 largest, then swap row J with row ROWPVT
        IF(ROWPVT.NE.J) THEN
          DO 70 L=1,N,1
             SCRVAR = LU(ROWPVT,L)
             LU(ROWPVT,L) = LU(J,L)
             LU(J,L) = SCRVAR
 70       CONTINUE
C                    Set the MAXVAL for row ROWPVT to that of the one
C                    that was just swapped in. 
          MAXVAL(ROWPVT) = MAXVAL(J)
C
        ENDIF
C
C                  Now divide the rest of the L column by the pivot
C                  element
        IF(DABS(LU(J,J)).GT.SMALLV )  THEN
          SCRVAR = LU(J,J)
        ELSE
          SCRVAR = SIGN(SMALLV,LU(J,J))
        ENDIF
        DO 80 I=J+1,N,1
          LU(I,J) = LU(I,J)/SCRVAR
 80     CONTINUE
C
C
 90   CONTINUE
C
C
 999  RETURN
C
      END
C
C***********************************************************************
C
C    MODULE:       MLUSLV
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         6/16/94
C
C    PURPOSE:      LU decomposes an NxN matrix
C
C    NOTE:         This subroutine employs a forward & back substitution
C                  algorithm with implicit pivoting as described in 
C                  "Numerical Recipes", by Press, Flannery, Teukolsky, &
C                  Vetterling, (pub. Cambridge Univ. Press) first
C                  edition, whose authors desrve the credit for
C                  this algorithm.
C
C    CALL LINE:    CALL MLUSLV(LU,X,B,IPVT,N)
C
C    INPUTS:       A = NxN LU decomposed array (D)
C                  B  = RHS vector (D)
C                  IPVT = Vector of pivots (I)
C                  N = Size of A (I)
C
C    OUTPUTS:      X = Solution vector(D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE MLUSLV(LU,X,B,IPVT,N)
C
      IMPLICIT NONE
C
C                 Parameters
      INTEGER N
      INTEGER IPVT(N)
      DOUBLE PRECISION LU(N,N), X(N), B(N)
C
C                 Local variables
C
C                 Loop variables
      INTEGER I, J, L
C
C                 A scratch variable
      DOUBLE PRECISION SCRV
C
C-----------------------------------------------------------------------
C
C                 Copy the RHS into X so we don't destroy B by
C                 unscrambling the pivots
      DO 10 I=1,N,1
        X(I) = B(I)
 10   CONTINUE
C
C                 Do the forward substitution
      DO 30 I=1,N,1
        L = IPVT(I)
        SCRV = X(L)
        X(L) = X(I)
        X(I) = SCRV
        DO 20 J=1,I-1,1
          X(I) = X(I)-LU(I,J)*X(J)
 20     CONTINUE
 30   CONTINUE
C
C
C                 Do the backward substitution
      DO 50 I=N,1,-1
        DO 40 J=I+1,N,1
          X(I) = X(I)-LU(I,J)*X(J)
 40     CONTINUE
        X(I) = X(I)/LU(I,I)
 50   CONTINUE
C
C
 999  RETURN
C
      END
C
C***********************************************************************
