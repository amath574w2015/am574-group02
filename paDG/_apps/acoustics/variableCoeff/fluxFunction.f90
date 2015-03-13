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
  DOUBLE PRECISION rhol,cl,rhor,cr,rho,c,K

  rhol=1.0d0
  cl=1.0d0
  rhor=2.0d0
  cr=0.5d0

  DO i=1,npts
   
    if (xpts(i)<0.d0) then
       rho=rhol
       c=cl
    else
       rho=rhor
       c=cr
    endif
    K=c*c*rho
    ! Acoustics
    fluxFunction(i,1) = K*qvals(i,2)
    fluxFunction(i,2) = qvals(i,1)/rho
  ENDDO !i
END FUNCTION fluxFunction
