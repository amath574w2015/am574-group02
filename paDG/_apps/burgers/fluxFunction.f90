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
    ! Burgers
    fluxFunction(i,1) = 0.5d0*qvals(i,1)*qvals(i,1)
  ENDDO !i
END FUNCTION fluxFunction
