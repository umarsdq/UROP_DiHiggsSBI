      SUBROUTINE MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=4)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=2)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=82, NLOOPGROUPS=19, NCTAMPS=41)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=123)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=13,NLOOPWAVEFUNCS=142)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=3, NSQUAREDSO=3, NAMPSO=4)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL MP_SXXXXX(P(0,4),+1*IC(4),W(1,4))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),GC_65,AMP(1))
      CALL MP_VVS12_3(W(1,1),W(1,2),GC_336,MDL_MH,MDL_WH,W(1,5))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_SSS6_0(W(1,3),W(1,4),W(1,5),GC_322,AMP(2))
      CALL MP_SSS6_1(W(1,3),W(1,4),GC_420,MDL_MH,MDL_WH,W(1,6))
C     Counter-term amplitude(s) for loop diagram number 3
      CALL MP_VVS11_0(W(1,1),W(1,2),W(1,6),R2GC_631_1814,AMPL(1,1))
C     At this point, all CT amps needed for (NP=4 QCD=2 QED=3), i.e.
C      of split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000
C     Counter-term amplitude(s) for loop diagram number 5
      CALL MP_VVSS14_0(W(1,1),W(1,2),W(1,4),W(1,3),R2GC_618_1801
     $ ,AMPL(1,2))
      CALL MP_SSS6_1(W(1,3),W(1,4),GC_322,MDL_MH,MDL_WH,W(1,7))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL MP_VVS11_0(W(1,1),W(1,2),W(1,7),R2GC_631_1814,AMPL(1,3))
C     At this point, all CT amps needed for (NP=2 QCD=2 QED=4), i.e.
C      of split order ID=2, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.2) GOTO 2000
C     Counter-term amplitude(s) for loop diagram number 13
      CALL MP_VVSS20_0(W(1,1),W(1,2),W(1,3),W(1,4),R2GC_569_1756
     $ ,AMPL(1,4))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253_1EPS
     $ ,AMPL(2,5))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253_1EPS
     $ ,AMPL(2,6))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253_1EPS
     $ ,AMPL(2,7))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253_1EPS
     $ ,AMPL(2,8))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253_1EPS
     $ ,AMPL(2,9))
      CALL MP_VVSS18_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1833_208_1EPS
     $ ,AMPL(2,10))
      CALL MP_VVSS19_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1753_108_1EPS
     $ ,AMPL(2,11))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253
     $ ,AMPL(1,12))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253
     $ ,AMPL(1,13))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253
     $ ,AMPL(1,14))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253
     $ ,AMPL(1,15))
      CALL MP_VVSS15_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_253
     $ ,AMPL(1,16))
      CALL MP_VVSS15_18_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1861_254
     $ ,UVGC_1833_208,AMPL(1,17))
      CALL MP_VVSS19_0(W(1,1),W(1,2),W(1,3),W(1,4),UVGC_1753_108
     $ ,AMPL(1,18))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_VVS16_0(W(1,1),W(1,2),W(1,7),R2GC_570_1757,AMPL(1,19))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351_1EPS,AMPL(2
     $ ,20))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351_1EPS,AMPL(2
     $ ,21))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351_1EPS,AMPL(2
     $ ,22))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351_1EPS,AMPL(2
     $ ,23))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351_1EPS,AMPL(2
     $ ,24))
      CALL MP_VVS13_0(W(1,1),W(1,2),W(1,7),UVGC_1794_165_1EPS,AMPL(2
     $ ,25))
      CALL MP_VVS15_0(W(1,1),W(1,2),W(1,7),UVGC_1755_110_1EPS,AMPL(2
     $ ,26))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351,AMPL(1,27))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351,AMPL(1,28))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351,AMPL(1,29))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351,AMPL(1,30))
      CALL MP_VVS12_0(W(1,1),W(1,2),W(1,7),UVGC_1918_351,AMPL(1,31))
      CALL MP_VVS12_13_0(W(1,1),W(1,2),W(1,7),UVGC_1918_352
     $ ,UVGC_1794_165,AMPL(1,32))
      CALL MP_VVS15_0(W(1,1),W(1,2),W(1,7),UVGC_1755_110,AMPL(1,33))
      CALL MP_VVS12P0_1(W(1,1),W(1,3),GC_336,ZERO,ZERO,W(1,8))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL MP_VVS11_0(W(1,2),W(1,8),W(1,4),R2GC_631_1814,AMPL(1,34))
      CALL MP_VVS12P0_1(W(1,1),W(1,4),GC_336,ZERO,ZERO,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL MP_VVS11_0(W(1,2),W(1,9),W(1,3),R2GC_631_1814,AMPL(1,35))
      CALL MP_VVS12P0_1(W(1,2),W(1,3),GC_336,ZERO,ZERO,W(1,10))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL MP_VVS11_0(W(1,1),W(1,10),W(1,4),R2GC_631_1814,AMPL(1,36))
      CALL MP_VVS12P0_1(W(1,2),W(1,4),GC_336,ZERO,ZERO,W(1,11))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL MP_VVS11_0(W(1,1),W(1,11),W(1,3),R2GC_631_1814,AMPL(1,37))
C     Counter-term amplitude(s) for loop diagram number 29
      CALL MP_VVSS11_14_16_0(W(1,1),W(1,2),W(1,3),W(1,4),R2GC_630_1813
     $ ,R2GC_1297_467,R2GC_629_1812,AMPL(1,38))
C     Counter-term amplitude(s) for loop diagram number 38
      CALL MP_VVS11_14_0(W(1,1),W(1,2),W(1,7),R2GC_1298_468
     $ ,R2GC_595_1779,AMPL(1,39))
      CALL MP_SSS10_1(W(1,3),W(1,4),GC_422,MDL_MH,MDL_WH,W(1,12))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL MP_VVS11_0(W(1,1),W(1,2),W(1,12),R2GC_631_1814,AMPL(1,40))
      CALL MP_SSS6_1(W(1,3),W(1,4),GC_431,MDL_MH,MDL_WH,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 40
      CALL MP_VVS11_0(W(1,1),W(1,2),W(1,13),R2GC_631_1814,AMPL(1,41))
C     At this point, all CT amps needed for (NP=4 QCD=2 QED=4), i.e.
C      of split order ID=3, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.3) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

