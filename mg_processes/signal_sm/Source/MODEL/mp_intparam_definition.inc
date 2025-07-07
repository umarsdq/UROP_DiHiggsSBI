ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        MP__G = 2 * SQRT(MP__AS*MP__PI)  ! for the first init

        MP__MDL_CW0 = MP__MDL_MW/MP__MDL_MZ

        MP__MDL_SQRT__2 = SQRT(CMPLX((2.000000E+00_16),KIND=16))

        MP__MDL_MUH0 = MP__MDL_MH/MP__MDL_SQRT__2

        MP__MDL_MW__EXP__2 = MP__MDL_MW**2

        MP__MDL_MZ__EXP__2 = MP__MDL_MZ**2

        MP__MDL_SW0 = SQRT(CMPLX((1.000000E+00_16-MP__MDL_MW__EXP__2
     $   /MP__MDL_MZ__EXP__2),KIND=16))

        MP__MDL_CF3PL1 = MP__MDL_C3PL1

        MP__MDL_CF3PL2 = MP__MDL_C3PL2

        MP__MDL_CFLL = 2.000000E+00_16*MP__MDL_CLL1221

        MP__MDL_NB__2__EXP__0_25 = 2.000000E+00_16**2.500000E-01_16

        MP__MDL_LAMBDA__EXP__2 = MP__MDL_LAMBDA**2

        MP__MDL_CTB = (-MP__MDL_CTZ+MP__MDL_CTW*MP__MDL_CW0)
     $   /MP__MDL_SW0

        MP__MDL_CW0__EXP__2 = MP__MDL_CW0**2

        MP__MDL_SW0__EXP__2 = MP__MDL_SW0**2

        MP__MDL_MH__EXP__2 = MP__MDL_MH**2

        MP__MDL_COMPLEXI = CMPLX(0.000000E+00_16,1.000000E+00_16
     $   ,KIND=16)

        MP__MDL_SW0__EXP__3 = MP__MDL_SW0**3

        MP__MDL_CW0__EXP__3 = MP__MDL_CW0**3

        MP__MDL_MT__EXP__2 = MP__MDL_MT**2

        MP__MDL_MT__EXP__3 = MP__MDL_MT**3

        MP__MDL_MUEFT__EXP__2 = MP__MDL_MUEFT**2

        MP__MDL_VEV0 = SQRT(CMPLX((1.000000E+00_16/MP__MDL_GF),KIND=16)
     $   )/MP__MDL_NB__2__EXP__0_25

        MP__MDL_VEV0__EXP__2 = MP__MDL_VEV0**2

        MP__MDL_DT = (MP__MDL_CPDC*MP__MDL_VEV0__EXP__2)/(2.000000E
     $   +00_16*MP__MDL_LAMBDA__EXP__2)

        MP__MDL_DV = ((MP__MDL_CF3PL1+MP__MDL_CF3PL2-MP__MDL_CFLL
     $   /2.000000E+00_16)*MP__MDL_VEV0__EXP__2)/MP__MDL_LAMBDA__EXP__2

        MP__MDL_DWB = (MP__MDL_CPWB*MP__MDL_VEV0__EXP__2)
     $   /MP__MDL_LAMBDA__EXP__2

        MP__MDL_DLAM = ((4.000000E+00_16*MP__MDL_CDP-MP__MDL_CPDC)
     $   *MP__MDL_VEV0__EXP__2)/(2.000000E+00_16
     $   *MP__MDL_LAMBDA__EXP__2)

        MP__MDL_EE0 = (2.000000E+00_16*MP__MDL_MW*MP__MDL_SW0)
     $   /MP__MDL_VEV0

        MP__MDL_CW = MP__MDL_CW0*(1.000000E+00_16+MP__MDL_DT/2.000000E
     $   +00_16)

        MP__MDL_MUH = (1.000000E+00_16-MP__MDL_DLAM/2.000000E+00_16)
     $   *MP__MDL_MUH0

        MP__MDL_SW = (1.000000E+00_16-(MP__MDL_CW0__EXP__2*MP__MDL_DT)
     $   /(2.000000E+00_16*MP__MDL_SW0__EXP__2))*MP__MDL_SW0

        MP__MDL_G1 = (MP__MDL_EE0*(1.000000E+00_16-MP__MDL_DV
     $   /2.000000E+00_16-MP__MDL_DT/(2.000000E+00_16
     $   *MP__MDL_SW0__EXP__2)-(MP__MDL_CW0*MP__MDL_DWB)/MP__MDL_SW0))
     $   /MP__MDL_CW0

        MP__MDL_GW = ((1.000000E+00_16-MP__MDL_DV/2.000000E+00_16)
     $   *MP__MDL_EE0)/MP__MDL_SW0

        MP__MDL_LAM = ((1.000000E+00_16-MP__MDL_DLAM-MP__MDL_DV)
     $   *MP__MDL_MH__EXP__2)/(2.000000E+00_16*MP__MDL_VEV0__EXP__2)

        MP__MDL_VEV = (1.000000E+00_16+MP__MDL_DV/2.000000E+00_16)
     $   *MP__MDL_VEV0

        MP__MDL_EE = (2.000000E+00_16*MP__MDL_MW*MP__MDL_SW*(1.000000E
     $   +00_16-MP__MDL_DV/2.000000E+00_16-(MP__MDL_CW0*MP__MDL_DWB)
     $   /MP__MDL_SW0))/MP__MDL_VEV0

        MP__MDL_EE__EXP__2 = MP__MDL_EE**2

        MP__MDL_AEW = MP__MDL_EE__EXP__2/(4.000000E+00_16*MP__PI)

        MP__AEWM1 = 1.000000E+00_16/MP__MDL_AEW

        MP__MDL_EE0__EXP__2 = MP__MDL_EE0**2

        MP__MDL_EE0__EXP__3 = MP__MDL_EE0**3

        MP__MDL_VEV0__EXP__3 = MP__MDL_VEV0**3

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      MP__AS = MP__G**2/4/MP__PI

      MP__MDL_SQRT__AS = SQRT(CMPLX((MP__AS),KIND=16))

      MP__MDL_G__EXP__2 = MP__G**2

      MP__MDL_G__EXP__3 = MP__G**3

      MP__MDL_G__EXP__4 = MP__G**4

      MP__MDL_MU_R__EXP__2 = MP__MU_R**2

C     
C     Parameters that should be updated for the loops.
C     
C     
C     Definition of the EW coupling used in the write out of aqed
C     
      MP__GAL(1) = 2 * SQRT(MP__PI/ABS(MP__AEWM1))
      MP__GAL(2) = 1D0

