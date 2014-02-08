C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       MAXWEL.INC
C    TYPE:         INCLUDE FILE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         6/13/90
C
C
C    CALL LINE:    INCLUDE 'MAXWEL.INC/LIST'
C
C*************************************************************************
C
C                            Number of T points in each boundary
      INTEGER NUMTMP
      PARAMETER (NUMTMP=201)
cc      PARAMETER (NUMTMP=21)
C
C                            Number of Ye points in each boundary
      INTEGER NUMYE
cc      PARAMETER(NUMYE=49)
cc      PARAMETER(NUMYE=34)
      PARAMETER(NUMYE=49)
C
C                            Number of points in alpha-nuclei
C                            boundary
      INTEGER NBPNTS, NUMLOW, NUMHI
      PARAMETER(NBPNTS=101,NUMLOW=51,NUMHI=51)
C
C
C                            Baryon densities
      DOUBLE PRECISION BRYLOW(NUMTMP,NUMYE), BRYHI(NUMTMP,NUMYE)
      DOUBLE PRECISION LOWDNS, HIDNS, DNS_1, DNS_2
      DOUBLE PRECISION DNL_DT, DNH_DT, DNL_DY, DNH_DY
C
C
      DOUBLE PRECISION LNLOW, LNHI, LNCUT, LNMINS, LNPLUS
      DOUBLE PRECISION LOGBRY, LOGBCH, DLTLN1, DLTLN2, LNFRAC
      PARAMETER(LNLOW=-8.92D0,LNHI=-0.92D0,LNCUT=-2.92D0)
C
C
      DOUBLE PRECISION YLOW,YHI,Y_CUT
cccc 180 MeV parms      PARAMETER(YLOW=0.03,YHI=0.51,Y_CUT=0.175)
cccc 220 MeV parms      PARAMETER(YLOW=0.03,YHI=0.51,Y_CUT=0.135)
cccc 375 MeV parms      PARAMETER(YLOW=0.03,YHI=0.51,Y_CUT=0.155)
cccc Pearson et al. parms PARAMETER(YLOW=0.03,YHI=0.51,Y_CUT=0.125)
      PARAMETER(YLOW=0.03D0,YHI=0.51D0,Y_CUT=0.155D0)
C
C
      DOUBLE PRECISION TCHK_B, TCHK_N, T_MXWL, D_MXWL
C
C
C                            Highest temperatures where Coexistence
C                            of bulk and nuclear phases occurs
      DOUBLE PRECISION T_H(NUMYE), D_H(NUMYE)
C
C
C                            Arrays containing number of boundary
C                            points
      DOUBLE PRECISION LBOUND(NBPNTS,NUMYE), UBOUND(NBPNTS,NUMYE)
C
C                            Minimum density code will work at
C
      DOUBLE PRECISION MINDNS
      PARAMETER(MINDNS=1.0D-10)
C
C
C
      DOUBLE PRECISION NEWDNS
      DOUBLE PRECISION T_LOW, T_HI, TFRAC
      DOUBLE PRECISION Y_LOW, Y_LOW2, Y_HI, Y_HI2, YFRAC, YCH, YCUT
      DOUBLE PRECISION LNL, LNH, LNC
C
      DOUBLE PRECISION YMINUS, YPLUS, YINTRP, TMINUS, TPLUS, TINTRP
      DOUBLE PRECISION DELT_Y, DELT_T
      INTEGER I_MXWL, J_MXWL
      INTEGER I_BD, I_BNDY, J_BD, J_BNDY
C
C                            Total chemical potential
      DOUBLE PRECISION MUTILD, MUTLOW, MUTHI
C
C                            Total pressure
      DOUBLE PRECISION PRTILD, PRLOW, PRHI
C
C
C                            Muhat
      DOUBLE PRECISION MUHLOW, MUHHI
C
C                            Electron chemical potential
      DOUBLE PRECISION MUELOW, MUEHI
C
C                            Total entropy per baryon
      DOUBLE PRECISION S_LOW, S_HI
C
C                            Total free energy density
      DOUBLE PRECISION F_LOW, F_HI
C
C
C                            Phase fraction
      DOUBLE PRECISION PHASEF
C
C
      COMMON /MAXCOM/ T_LOW,T_HI,Y_LOW,Y_HI,
     1       BRYLOW,BRYHI,LBOUND,UBOUND,T_H,D_H,
     2       LOWDNS,HIDNS,TCHK_N,TCHK_B,T_MXWL,D_MXWL,
     3       DNL_DT,DNH_DT,DNL_DY,DNH_DY
C
C























