        !COMPILER-GENERATED INTERFACE MODULE: Wed Feb 03 10:50:36 2010
        MODULE SUB1__genmod
          INTERFACE 
            SUBROUTINE SUB1(INTARG,CALLBACK)
              INTEGER(KIND=4), INTENT(IN) :: INTARG
              INTERFACE 
                SUBROUTINE CALLBACK(INTARG)
                  INTEGER(KIND=4), INTENT(IN) :: INTARG
                END SUBROUTINE CALLBACK
              END INTERFACE 
            END SUBROUTINE SUB1
          END INTERFACE 
        END MODULE SUB1__genmod
