FUNCTION fluxFunction(qvals,xpts,npts)
  ! ========================================================================
  ! User specified flux functions evaluated at xpts
  ! INPUTS: npts - number of points to evaluate at
  !         xpts(npts) - points to evaluate at
  !         qvals(npts,m) - mth equation solution values at xpts
  ! OUTPUT:
  !         fluxFunction(npts,m) - analytic flux function evaluated at xpts
  ! ========================================================================
  USE testParameters, ONLY: meqn,u0,K0,rho0
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: npts
  DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xpts
  DOUBLE PRECISION, DIMENSION(1:npts,1:meqn), INTENT(IN) :: qvals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:npts,1:meqn) :: fluxFunction
  ! Local variables
  INTEGER :: i
  DO i=1,npts
    ! Advection
    fluxFunction(i,1) = u0*qvals(i,1)

    ! Acoustics
    ! fluxFunction(i,1) = u0*qvals(i,1)+K0*qvals(i,2)
    ! fluxFunction(i,2) = qvals(i,1)/rho0+u0*qvals(i,2)
  ENDDO !i
END FUNCTION fluxFunction
