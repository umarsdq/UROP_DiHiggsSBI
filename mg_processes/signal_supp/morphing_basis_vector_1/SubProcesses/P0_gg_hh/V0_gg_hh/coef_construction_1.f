      SUBROUTINE COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
C     
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
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=3, NSQUAREDSO=3, NAMPSO=4)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 3
      CALL FFV87L1_2(PL(0,0),W(1,1),GC_11,MDL_MT,MDL_WT,PL(0,1),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL FFV87L1_2(PL(0,1),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,2),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL FFS37L1_2(PL(0,2),W(1,6),GC_546,MDL_MT,MDL_WT,PL(0,3),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,3))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,3),3,4,1,1,1,42,H)
C     Coefficient construction for loop diagram with ID 4
      CALL FFV87L2_1(PL(0,0),W(1,1),GC_11,MDL_MT,MDL_WT,PL(0,4),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,4))
      CALL FFV87L2_1(PL(0,4),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,5),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,5))
      CALL FFS37L2_1(PL(0,5),W(1,6),GC_546,MDL_MT,MDL_WT,PL(0,6),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,6))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,6),3,4,1,1,1,43,H)
C     At this point, all loop coefficients needed for (NP=4 QCD=2
C      QED=3), i.e. of split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000
C     Coefficient construction for loop diagram with ID 5
      CALL FFS37L1_2(PL(0,2),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,7),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,7))
      CALL FFS37L1_2(PL(0,7),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,8),COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,8),4,4,2,1,1,44,H)
C     Coefficient construction for loop diagram with ID 6
      CALL FFS37L1_2(PL(0,2),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,9),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,9))
      CALL FFS37L1_2(PL(0,9),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,10)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,10))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,10),4,4,3,1,1,45,H)
C     Coefficient construction for loop diagram with ID 7
      CALL FFS37L1_2(PL(0,2),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,11)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,11))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,11),3,4,1,1,1,46,H)
C     Coefficient construction for loop diagram with ID 8
      CALL FFS37L2_1(PL(0,5),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,12)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,12))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,12),3,4,1,1,1,47,H)
C     Coefficient construction for loop diagram with ID 9
      CALL FFS37L2_1(PL(0,5),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,13)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,13))
      CALL FFS37L2_1(PL(0,13),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,14)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,14))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,14),4,4,2,1,1,48,H)
C     Coefficient construction for loop diagram with ID 10
      CALL FFS37L1_2(PL(0,1),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,15)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,15))
      CALL FFV87L1_2(PL(0,15),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,16)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,16))
      CALL FFS37L1_2(PL(0,16),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,17)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,16),4,COEFS,4,4,WL(1,0,1,17))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,17),4,4,4,1,1,49,H)
C     Coefficient construction for loop diagram with ID 11
      CALL FFS37L2_1(PL(0,5),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,18)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,18))
      CALL FFS37L2_1(PL(0,18),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,19)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,18),4,COEFS,4,4,WL(1,0,1,19))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,19),4,4,3,1,1,50,H)
C     Coefficient construction for loop diagram with ID 12
      CALL FFS37L2_1(PL(0,4),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,20)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,20))
      CALL FFV87L2_1(PL(0,20),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,21)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,21))
      CALL FFS37L2_1(PL(0,21),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,22)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,22))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,22),4,4,4,1,1,51,H)
C     At this point, all loop coefficients needed for (NP=2 QCD=2
C      QED=4), i.e. of split order ID=2, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.2) GOTO 4000
C     Coefficient construction for loop diagram with ID 13
      CALL VVV8L2P0_1(PL(0,0),W(1,1),GC_10,ZERO,ZERO,PL(0,23),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,23))
      CALL VVV8L2P0_1(PL(0,23),W(1,2),GC_10,ZERO,ZERO,PL(0,24),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,24))
      CALL VVSS15L2P0_1(PL(0,24),W(1,3),W(1,4),GC_65,ZERO,ZERO,PL(0,25)
     $ ,COEFS)
      CALL UPDATE_WL_2_2(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,25))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,25),4,4,5,1,1,52,H)
C     Coefficient construction for loop diagram with ID 14
      CALL VVS12L2P0_1(PL(0,24),W(1,7),GC_336,ZERO,ZERO,PL(0,26),COEFS)
      CALL UPDATE_WL_2_2(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,26))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,26),4,4,5,1,1,53,H)
C     Coefficient construction for loop diagram with ID 15
      CALL VVVS17L2P0_1(PL(0,23),W(1,2),W(1,7),GC_357,ZERO,ZERO,PL(0
     $ ,27),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,27))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,27),2,4,6,2,1,54,H)
C     Coefficient construction for loop diagram with ID 16
      CALL VVVSS25L2P0_1(PL(0,23),W(1,2),W(1,3),W(1,4),GC_153,ZERO
     $ ,ZERO,PL(0,28),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,28))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,28),2,4,6,2,1,55,H)
C     Coefficient construction for loop diagram with ID 17
      CALL VVV8L2P0_1(PL(0,0),W(1,2),GC_10,ZERO,ZERO,PL(0,29),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,29))
      CALL VVVS17L2P0_1(PL(0,29),W(1,1),W(1,7),GC_357,ZERO,ZERO,PL(0
     $ ,30),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,29),4,COEFS,4,4,WL(1,0,1,30))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,30),2,4,7,2,1,56,H)
C     Coefficient construction for loop diagram with ID 18
      CALL VVVV16L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,31)
     $ ,COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,31))
      CALL VVVV19L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,32)
     $ ,COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,32))
      CALL VVVV20L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,33)
     $ ,COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,33))
      CALL VVS12L2P0_1(PL(0,31),W(1,7),GC_336,ZERO,ZERO,PL(0,34),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,31),4,COEFS,4,4,WL(1,0,1,34))
      CALL VVS12L2P0_1(PL(0,32),W(1,7),GC_336,ZERO,ZERO,PL(0,35),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,32),4,COEFS,4,4,WL(1,0,1,35))
      CALL VVS12L2P0_1(PL(0,33),W(1,7),GC_336,ZERO,ZERO,PL(0,36),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,33),4,COEFS,4,4,WL(1,0,1,36))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,34),2,4,8,2,1,57,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,35),2,4,8,2,1,58,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,36),2,4,8,2,1,59,H)
C     Coefficient construction for loop diagram with ID 19
      CALL VVSS15L2P0_1(PL(0,31),W(1,3),W(1,4),GC_65,ZERO,ZERO,PL(0,37)
     $ ,COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,31),4,COEFS,4,4,WL(1,0,1,37))
      CALL VVSS15L2P0_1(PL(0,32),W(1,3),W(1,4),GC_65,ZERO,ZERO,PL(0,38)
     $ ,COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,32),4,COEFS,4,4,WL(1,0,1,38))
      CALL VVSS15L2P0_1(PL(0,33),W(1,3),W(1,4),GC_65,ZERO,ZERO,PL(0,39)
     $ ,COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,33),4,COEFS,4,4,WL(1,0,1,39))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,37),2,4,8,2,1,60,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,38),2,4,8,2,1,61,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,39),2,4,8,2,1,62,H)
C     Coefficient construction for loop diagram with ID 20
      CALL VVVSS25L2P0_1(PL(0,29),W(1,1),W(1,3),W(1,4),GC_153,ZERO
     $ ,ZERO,PL(0,40),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,29),4,COEFS,4,4,WL(1,0,1,40))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,40),2,4,7,2,1,63,H)
C     Coefficient construction for loop diagram with ID 21
      CALL FFV87L2_1(PL(0,0),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,41),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,41))
      CALL FFS37L2_1(PL(0,41),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,42)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,41),4,COEFS,4,4,WL(1,0,1,42))
      CALL FFV87L2_1(PL(0,42),W(1,8),GC_11,MDL_MT,MDL_WT,PL(0,43)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,42),4,COEFS,4,4,WL(1,0,1,43))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,43),3,4,9,1,1,64,H)
C     Coefficient construction for loop diagram with ID 22
      CALL FFV87L1_2(PL(0,0),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,44),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,44))
      CALL FFS37L1_2(PL(0,44),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,45)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,44),4,COEFS,4,4,WL(1,0,1,45))
      CALL FFV87L1_2(PL(0,45),W(1,8),GC_11,MDL_MT,MDL_WT,PL(0,46)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,45),4,COEFS,4,4,WL(1,0,1,46))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,46),3,4,9,1,1,65,H)
C     Coefficient construction for loop diagram with ID 23
      CALL FFS37L2_1(PL(0,41),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,47)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,41),4,COEFS,4,4,WL(1,0,1,47))
      CALL FFV87L2_1(PL(0,47),W(1,9),GC_11,MDL_MT,MDL_WT,PL(0,48)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,47),4,COEFS,4,4,WL(1,0,1,48))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,48),3,4,10,1,1,66,H)
C     Coefficient construction for loop diagram with ID 24
      CALL FFS37L1_2(PL(0,44),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,49)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,44),4,COEFS,4,4,WL(1,0,1,49))
      CALL FFV87L1_2(PL(0,49),W(1,9),GC_11,MDL_MT,MDL_WT,PL(0,50)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,49),4,COEFS,4,4,WL(1,0,1,50))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,50),3,4,10,1,1,67,H)
C     Coefficient construction for loop diagram with ID 25
      CALL FFS37L2_1(PL(0,4),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,51)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,51))
      CALL FFV87L2_1(PL(0,51),W(1,10),GC_11,MDL_MT,MDL_WT,PL(0,52)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,51),4,COEFS,4,4,WL(1,0,1,52))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,52),3,4,11,1,1,68,H)
C     Coefficient construction for loop diagram with ID 26
      CALL FFS37L1_2(PL(0,1),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,53)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,53))
      CALL FFV87L1_2(PL(0,53),W(1,10),GC_11,MDL_MT,MDL_WT,PL(0,54)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,53),4,COEFS,4,4,WL(1,0,1,54))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,54),3,4,11,1,1,69,H)
C     Coefficient construction for loop diagram with ID 27
      CALL FFV87L2_1(PL(0,20),W(1,11),GC_11,MDL_MT,MDL_WT,PL(0,55)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,55))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,55),3,4,12,1,1,70,H)
C     Coefficient construction for loop diagram with ID 28
      CALL FFV87L1_2(PL(0,15),W(1,11),GC_11,MDL_MT,MDL_WT,PL(0,56)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,56))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,56),3,4,12,1,1,71,H)
C     Coefficient construction for loop diagram with ID 29
      CALL FFSS26L1_2(PL(0,2),W(1,3),W(1,4),GC_1822,MDL_MT,MDL_WT,PL(0
     $ ,57),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,57))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,57),3,4,1,1,1,72,H)
C     Coefficient construction for loop diagram with ID 30
      CALL FFS37L1_2(PL(0,7),W(1,3),GC_1842,MDL_MT,MDL_WT,PL(0,58)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,58))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,58),4,4,2,1,1,73,H)
C     Coefficient construction for loop diagram with ID 31
      CALL FFS37L1_2(PL(0,2),W(1,4),GC_1842,MDL_MT,MDL_WT,PL(0,59)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,59))
      CALL FFS37L1_2(PL(0,59),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,60)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,59),4,COEFS,4,4,WL(1,0,1,60))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,60),4,4,2,1,1,74,H)
C     Coefficient construction for loop diagram with ID 32
      CALL FFS37L1_2(PL(0,9),W(1,4),GC_1842,MDL_MT,MDL_WT,PL(0,61)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,61))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,61),4,4,3,1,1,75,H)
C     Coefficient construction for loop diagram with ID 33
      CALL FFS37L1_2(PL(0,2),W(1,3),GC_1842,MDL_MT,MDL_WT,PL(0,62)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,62))
      CALL FFS37L1_2(PL(0,62),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,63)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,62),4,COEFS,4,4,WL(1,0,1,63))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,63),4,4,3,1,1,76,H)
C     Coefficient construction for loop diagram with ID 34
      CALL FFV110L1_2(PL(0,1),W(1,2),GC_358,MDL_MT,MDL_WT,PL(0,64)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,64))
      CALL FFS37L1_2(PL(0,64),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,65)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,64),4,COEFS,4,4,WL(1,0,1,65))
      CALL FFS37L1_2(PL(0,65),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,66)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,65),4,COEFS,4,4,WL(1,0,1,66))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,66),4,4,2,1,1,77,H)
C     Coefficient construction for loop diagram with ID 35
      CALL FFS37L1_2(PL(0,64),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,67)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,64),4,COEFS,4,4,WL(1,0,1,67))
      CALL FFS37L1_2(PL(0,67),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,68)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,67),4,COEFS,4,4,WL(1,0,1,68))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,68),4,4,3,1,1,78,H)
C     Coefficient construction for loop diagram with ID 36
      CALL FFV110L1_2(PL(0,0),W(1,1),GC_358,MDL_MT,MDL_WT,PL(0,69)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,69))
      CALL FFV87L1_2(PL(0,69),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,70)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,69),4,COEFS,4,4,WL(1,0,1,70))
      CALL FFS37L1_2(PL(0,70),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,71)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,70),4,COEFS,4,4,WL(1,0,1,71))
      CALL FFS37L1_2(PL(0,71),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,72)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,71),4,COEFS,4,4,WL(1,0,1,72))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,72),4,4,2,1,1,79,H)
C     Coefficient construction for loop diagram with ID 37
      CALL FFS37L1_2(PL(0,70),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,73)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,70),4,COEFS,4,4,WL(1,0,1,73))
      CALL FFS37L1_2(PL(0,73),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,74)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,73),4,COEFS,4,4,WL(1,0,1,74))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,74),4,4,3,1,1,80,H)
C     Coefficient construction for loop diagram with ID 38
      CALL FFS37L1_2(PL(0,2),W(1,7),GC_1842,MDL_MT,MDL_WT,PL(0,75)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,75))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,75),3,4,1,1,1,81,H)
C     Coefficient construction for loop diagram with ID 39
      CALL FFS37L1_2(PL(0,2),W(1,12),GC_546,MDL_MT,MDL_WT,PL(0,76)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,76))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,76),3,4,1,1,1,82,H)
C     Coefficient construction for loop diagram with ID 40
      CALL FFS37L1_2(PL(0,2),W(1,13),GC_546,MDL_MT,MDL_WT,PL(0,77)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,77))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,77),3,4,1,1,1,83,H)
C     Coefficient construction for loop diagram with ID 41
      CALL FFS37L1_2(PL(0,64),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,78)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,64),4,COEFS,4,4,WL(1,0,1,78))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,78),3,4,1,1,1,84,H)
C     Coefficient construction for loop diagram with ID 42
      CALL FFS37L1_2(PL(0,70),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,79)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,70),4,COEFS,4,4,WL(1,0,1,79))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,79),3,4,1,1,1,85,H)
C     Coefficient construction for loop diagram with ID 43
      CALL FFVS127L2_1(PL(0,4),W(1,2),W(1,7),GC_156,MDL_MT,MDL_WT,PL(0
     $ ,80),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,80))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,80),2,4,13,1,1,86,H)
C     Coefficient construction for loop diagram with ID 44
      CALL FFS37L2_1(PL(0,5),W(1,7),GC_1842,MDL_MT,MDL_WT,PL(0,81)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,81))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,81),3,4,1,1,1,87,H)
C     Coefficient construction for loop diagram with ID 45
      CALL FFV110L2_1(PL(0,4),W(1,2),GC_358,MDL_MT,MDL_WT,PL(0,82)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,82))
      CALL FFS37L2_1(PL(0,82),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,83)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,82),4,COEFS,4,4,WL(1,0,1,83))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,83),3,4,1,1,1,88,H)
C     Coefficient construction for loop diagram with ID 46
      CALL FFS37L2_1(PL(0,5),W(1,12),GC_546,MDL_MT,MDL_WT,PL(0,84)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,84))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,84),3,4,1,1,1,89,H)
C     Coefficient construction for loop diagram with ID 47
      CALL FFS37L2_1(PL(0,5),W(1,13),GC_546,MDL_MT,MDL_WT,PL(0,85)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,85))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,85),3,4,1,1,1,90,H)
C     Coefficient construction for loop diagram with ID 48
      CALL FFV110L2_1(PL(0,0),W(1,1),GC_358,MDL_MT,MDL_WT,PL(0,86)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,86))
      CALL FFV87L2_1(PL(0,86),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,87)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,86),4,COEFS,4,4,WL(1,0,1,87))
      CALL FFS37L2_1(PL(0,87),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,88)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,87),4,COEFS,4,4,WL(1,0,1,88))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,88),3,4,1,1,1,91,H)
C     Coefficient construction for loop diagram with ID 49
      CALL FFS37L2_1(PL(0,13),W(1,3),GC_1842,MDL_MT,MDL_WT,PL(0,89)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,89))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,89),4,4,2,1,1,92,H)
C     Coefficient construction for loop diagram with ID 50
      CALL FFS37L1_2(PL(0,1),W(1,3),GC_1842,MDL_MT,MDL_WT,PL(0,90)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,90))
      CALL FFV87L1_2(PL(0,90),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,91)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,90),4,COEFS,4,4,WL(1,0,1,91))
      CALL FFS37L1_2(PL(0,91),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,92)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,91),4,COEFS,4,4,WL(1,0,1,92))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,92),4,4,4,1,1,93,H)
C     Coefficient construction for loop diagram with ID 51
      CALL FFVS127L1_2(PL(0,15),W(1,2),W(1,4),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,93),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,93))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,93),3,4,12,1,1,94,H)
C     Coefficient construction for loop diagram with ID 52
      CALL FFS37L2_1(PL(0,5),W(1,4),GC_1842,MDL_MT,MDL_WT,PL(0,94)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,94))
      CALL FFS37L2_1(PL(0,94),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,95)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,94),4,COEFS,4,4,WL(1,0,1,95))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,95),4,4,2,1,1,95,H)
C     Coefficient construction for loop diagram with ID 53
      CALL FFS37L2_1(PL(0,82),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,96)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,82),4,COEFS,4,4,WL(1,0,1,96))
      CALL FFS37L2_1(PL(0,96),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,97)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,96),4,COEFS,4,4,WL(1,0,1,97))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,97),4,4,2,1,1,96,H)
C     Coefficient construction for loop diagram with ID 54
      CALL FFS37L1_2(PL(0,16),W(1,4),GC_1842,MDL_MT,MDL_WT,PL(0,98)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,16),4,COEFS,4,4,WL(1,0,1,98))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,98),4,4,4,1,1,97,H)
C     Coefficient construction for loop diagram with ID 55
      CALL FFV110L1_2(PL(0,15),W(1,2),GC_358,MDL_MT,MDL_WT,PL(0,99)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,99))
      CALL FFS37L1_2(PL(0,99),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,100)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,99),4,COEFS,4,4,WL(1,0,1,100))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,100),4,4,4,1,1,98,H)
C     Coefficient construction for loop diagram with ID 56
      CALL FFS37L2_1(PL(0,87),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,101)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,87),4,COEFS,4,4,WL(1,0,1,101))
      CALL FFS37L2_1(PL(0,101),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,102)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,101),4,COEFS,4,4,WL(1,0,1,102))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,102),4,4,2,1,1,99,H)
C     Coefficient construction for loop diagram with ID 57
      CALL FFS37L1_2(PL(0,69),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,103)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,69),4,COEFS,4,4,WL(1,0,1,103))
      CALL FFV87L1_2(PL(0,103),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,104)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,103),4,COEFS,4,4,WL(1,0,1,104))
      CALL FFS37L1_2(PL(0,104),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,105)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,104),4,COEFS,4,4,WL(1,0,1,105))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,105),4,4,4,1,1,100,H)
C     Coefficient construction for loop diagram with ID 58
      CALL FFS37L2_1(PL(0,18),W(1,4),GC_1842,MDL_MT,MDL_WT,PL(0,106)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,18),4,COEFS,4,4,WL(1,0,1,106))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,106),4,4,3,1,1,101,H)
C     Coefficient construction for loop diagram with ID 59
      CALL FFS37L2_1(PL(0,21),W(1,4),GC_1842,MDL_MT,MDL_WT,PL(0,107)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,107))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,107),4,4,4,1,1,102,H)
C     Coefficient construction for loop diagram with ID 60
      CALL FFVS127L1_2(PL(0,53),W(1,2),W(1,3),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,108),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,53),4,COEFS,4,4,WL(1,0,1,108))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,108),3,4,11,1,1,103,H)
C     Coefficient construction for loop diagram with ID 61
      CALL FFS37L2_1(PL(0,5),W(1,3),GC_1842,MDL_MT,MDL_WT,PL(0,109)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,109))
      CALL FFS37L2_1(PL(0,109),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,110)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,109),4,COEFS,4,4,WL(1,0,1,110))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,110),4,4,3,1,1,104,H)
C     Coefficient construction for loop diagram with ID 62
      CALL FFS37L2_1(PL(0,82),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,111)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,82),4,COEFS,4,4,WL(1,0,1,111))
      CALL FFS37L2_1(PL(0,111),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,112)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,111),4,COEFS,4,4,WL(1,0,1,112))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,112),4,4,3,1,1,105,H)
C     Coefficient construction for loop diagram with ID 63
      CALL FFS37L2_1(PL(0,4),W(1,3),GC_1842,MDL_MT,MDL_WT,PL(0,113)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,113))
      CALL FFV87L2_1(PL(0,113),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,114)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,113),4,COEFS,4,4,WL(1,0,1,114))
      CALL FFS37L2_1(PL(0,114),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,115)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,114),4,COEFS,4,4,WL(1,0,1,115))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,115),4,4,4,1,1,106,H)
C     Coefficient construction for loop diagram with ID 64
      CALL FFV110L2_1(PL(0,20),W(1,2),GC_358,MDL_MT,MDL_WT,PL(0,116)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,116))
      CALL FFS37L2_1(PL(0,116),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,117)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,116),4,COEFS,4,4,WL(1,0,1,117))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,117),4,4,4,1,1,107,H)
C     Coefficient construction for loop diagram with ID 65
      CALL FFS37L2_1(PL(0,87),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,118)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,87),4,COEFS,4,4,WL(1,0,1,118))
      CALL FFS37L2_1(PL(0,118),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,119)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,118),4,COEFS,4,4,WL(1,0,1,119))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,119),4,4,3,1,1,108,H)
C     Coefficient construction for loop diagram with ID 66
      CALL FFS37L2_1(PL(0,86),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,120)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,86),4,COEFS,4,4,WL(1,0,1,120))
      CALL FFV87L2_1(PL(0,120),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,121)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,120),4,COEFS,4,4,WL(1,0,1,121))
      CALL FFS37L2_1(PL(0,121),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,122)
     $ ,COEFS)
      CALL UPDATE_WL_3_1(WL(1,0,1,121),4,COEFS,4,4,WL(1,0,1,122))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,122),4,4,4,1,1,109,H)
C     Coefficient construction for loop diagram with ID 67
      CALL FFVS127L2_1(PL(0,51),W(1,2),W(1,3),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,123),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,51),4,COEFS,4,4,WL(1,0,1,123))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,123),3,4,11,1,1,110,H)
C     Coefficient construction for loop diagram with ID 68
      CALL FFVS127L2_1(PL(0,20),W(1,2),W(1,4),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,124),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,124))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,124),3,4,12,1,1,111,H)
C     Coefficient construction for loop diagram with ID 69
      CALL FFSS26L2_1(PL(0,5),W(1,3),W(1,4),GC_1822,MDL_MT,MDL_WT,PL(0
     $ ,125),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,125))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,125),3,4,1,1,1,112,H)
C     Coefficient construction for loop diagram with ID 70
      CALL FFVS127L2_1(PL(0,41),W(1,1),W(1,7),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,126),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,41),4,COEFS,4,4,WL(1,0,1,126))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,126),2,4,14,1,1,113,H)
C     Coefficient construction for loop diagram with ID 71
      CALL FFVS127L1_2(PL(0,49),W(1,1),W(1,4),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,127),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,49),4,COEFS,4,4,WL(1,0,1,127))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,127),3,4,10,1,1,114,H)
C     Coefficient construction for loop diagram with ID 72
      CALL FFVS127L1_2(PL(0,45),W(1,1),W(1,3),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,128),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,45),4,COEFS,4,4,WL(1,0,1,128))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,128),3,4,9,1,1,115,H)
C     Coefficient construction for loop diagram with ID 73
      CALL FFVS127L2_1(PL(0,47),W(1,1),W(1,4),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,129),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,47),4,COEFS,4,4,WL(1,0,1,129))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,129),3,4,10,1,1,116,H)
C     Coefficient construction for loop diagram with ID 74
      CALL FFVS127L2_1(PL(0,42),W(1,1),W(1,3),GC_156,MDL_MT,MDL_WT
     $ ,PL(0,130),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,42),4,COEFS,4,4,WL(1,0,1,130))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,130),3,4,9,1,1,117,H)
C     Coefficient construction for loop diagram with ID 75
      CALL FFS37L1_2(PL(0,0),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,131)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,131))
      CALL FFS37L1_2(PL(0,131),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,132)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,131),4,COEFS,4,4,WL(1,0,1,132))
      CALL FFVV186L1_2(PL(0,132),W(1,1),W(1,2),GC_360,MDL_MT,MDL_WT
     $ ,PL(0,133),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,132),4,COEFS,4,4,WL(1,0,1,133))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,133),3,4,15,1,1,118,H)
C     Coefficient construction for loop diagram with ID 76
      CALL FFS37L2_1(PL(0,0),W(1,3),GC_546,MDL_MT,MDL_WT,PL(0,134)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,134))
      CALL FFS37L2_1(PL(0,134),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,135)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,134),4,COEFS,4,4,WL(1,0,1,135))
      CALL FFVV186L2_1(PL(0,135),W(1,1),W(1,2),GC_360,MDL_MT,MDL_WT
     $ ,PL(0,136),COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,135),4,COEFS,4,4,WL(1,0,1,136))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,136),3,4,15,1,1,119,H)
C     Coefficient construction for loop diagram with ID 77
      CALL FFVV186L2_1(PL(0,0),W(1,1),W(1,2),GC_360,MDL_MT,MDL_WT,PL(0
     $ ,137),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,137))
      CALL FFS37L2_1(PL(0,137),W(1,7),GC_546,MDL_MT,MDL_WT,PL(0,138)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,137),4,COEFS,4,4,WL(1,0,1,138))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,138),2,4,16,1,1,120,H)
C     Coefficient construction for loop diagram with ID 78
      CALL FFS37L2_1(PL(0,0),W(1,4),GC_546,MDL_MT,MDL_WT,PL(0,139)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,139))
      CALL FFVVS151L2_1(PL(0,139),W(1,1),W(1,2),W(1,3),GC_161,MDL_MT
     $ ,MDL_WT,PL(0,140),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,139),4,COEFS,4,4,WL(1,0,1,140))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,140),2,4,17,1,1,121,H)
C     Coefficient construction for loop diagram with ID 79
      CALL FFVVS151L2_1(PL(0,134),W(1,1),W(1,2),W(1,4),GC_161,MDL_MT
     $ ,MDL_WT,PL(0,141),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,134),4,COEFS,4,4,WL(1,0,1,141))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,141),2,4,18,1,1,122,H)
C     Coefficient construction for loop diagram with ID 80
      CALL FFVVS151L2_1(PL(0,0),W(1,1),W(1,2),W(1,7),GC_161,MDL_MT
     $ ,MDL_WT,PL(0,142),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,142))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,142),1,4,19,1,1,123,H)
C     At this point, all loop coefficients needed for (NP=4 QCD=2
C      QED=4), i.e. of split order ID=3, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.3) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

