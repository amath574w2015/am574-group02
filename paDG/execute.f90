! ================================================================
! Polynomial adaptive discontinuous Galerkin (DG)
! By: Devin Light & Scott Moe
! For AMATH 574 -- March 2015
! =================================================================

PROGRAM EXECUTE
  USE testParameters
  USE gridData
  IMPLICIT NONE

  INTEGER :: nRuns,nScale,noutput,nex0
  INTEGER, PARAMETER :: inUnit=20,outUnit=21
  DOUBLE PRECISION :: cflCoeff,muMax,xLower,xUpper
  LOGICAL :: DEBUG

  INTERFACE
    SUBROUTINE DRIVER(ntest,nex0,nscale,nruns,noutput,maxCFL)
      INTEGER, INTENT(IN) :: ntest,nex0,nscale,nruns,noutput
      REAL(KIND=8), INTENT(IN) :: maxCFL
    END SUBROUTINE DRIVER
  END INTERFACE

  ! Read input parameters from $(APPSDIR)/inputs.nl
  NAMELIST /inputs/ meqn,u0,K0,rho0,xLower,xUpper,nex0,nRuns,&
                    nScale,maxPolyDegree,cflCoeff,tfinal,&
                    noutput,outdir,testID,DEBUG

  OPEN(unit=inUnit,file="_apps/inputs.nl",action="read")
  READ(inUnit,NML=inputs)

  xDomain(1) = xLower
  xDomain(2) = xUpper

  write(*,*) meqn
  muMAX  = determineCFL(maxpolyDegree,cflCoeff)

  write(*,*) '======================================================'
  write(*,*) '             BEGINNING RUN OF MODAL TESTS             '
  write(*,'(A27,F7.4)') 'muMAX=',muMAX
  write(*,*) '======================================================'

  write(*,*) '======'
  SELECT CASE(testID)
    CASE(0)
      write(*,*) 'TEST 0: Consistency test'
    CASE(1)
      write(*,*) 'TEST 1: Uniform advection of 2 Sine Waves'
    CASE(2)
      write(*,*) 'TEST 2: Uniform advection of cosinebell'
    CASE(3)
      write(*,*) 'TEST 3: Uniform advection of square wave'
  END SELECT
  write(*,*) '======'

  CALL DRIVER(testID,nex0,nScale,nRuns,noutput,muMAX)
  CLOSE(inUnit)

  WRITE(*,*) 'PROGRAM COMPLETE!'
CONTAINS
  FUNCTION determineCFL(maxPolyDegree,cflCoeff)
    ! ===============================================================
    ! Determines max CFL for SSPRK3 timestepping as a function of
    ! maximum reconstructing polynomial degree
    ! ===============================================================
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: maxPolyDegree
    DOUBLE PRECISION :: cflCoeff
    ! Outputs
    DOUBLE PRECISION determineCFL
    ! Local variables

    determineCFL=9.d0/5.d0/((maxPolyDegree+1.d0)*(maxPolyDegree+1.d0))
    SELECT CASE(maxPolyDegree)
      CASE(2)
        determineCFL = 0.209D0
      CASE(3)
        determineCFL = 0.130D0
      CASE(4)
        determineCFL = 0.089D0
      CASE(5)
        determineCFL = 0.066D0
      CASE(6)
        determineCFL = 0.051D0
      CASE(7)
        determineCFL = 0.04D0
      CASE(8)
        determineCFL = 0.033D0
      CASE(9)
        determineCFL = 0.026D0
    END SELECT
    determineCFL = determineCFL*cflCoeff

  END FUNCTION determineCFL
END PROGRAM EXECUTE
