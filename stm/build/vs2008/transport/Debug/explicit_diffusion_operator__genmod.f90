        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 15 15:22:21 2009
        MODULE EXPLICIT_DIFFUSION_OPERATOR__genmod
          INTERFACE 
            SUBROUTINE EXPLICIT_DIFFUSION_OPERATOR(                     &
     &EXPLICIT_DIFFUSION_TERM,CONC,CONC_PREV,MASS,MASS_PREV,AREA,AREA_LO&
     &,AREA_HI,KS_LO,KS_HI,NCELL,NVAR,TIME,STRT_FLX_BC_PREV,            &
     &END_FLX_BC_PREV,DT,DX)
              INTEGER(KIND=4), INTENT(IN) :: NVAR
              INTEGER(KIND=4), INTENT(IN) :: NCELL
              REAL(KIND=8), INTENT(OUT) :: EXPLICIT_DIFFUSION_TERM(NCELL&
     &,NVAR)
              REAL(KIND=8), INTENT(OUT) :: CONC(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: CONC_PREV(NCELL,NVAR)
              REAL(KIND=8), INTENT(OUT) :: MASS(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: MASS_PREV(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: AREA(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: AREA_LO(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: AREA_HI(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: KS_LO(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: KS_HI(NCELL,NVAR)
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(IN) :: STRT_FLX_BC_PREV
              REAL(KIND=8), INTENT(IN) :: END_FLX_BC_PREV
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: DX
            END SUBROUTINE EXPLICIT_DIFFUSION_OPERATOR
          END INTERFACE 
        END MODULE EXPLICIT_DIFFUSION_OPERATOR__genmod
