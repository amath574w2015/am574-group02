SUBROUTINE qinit(xvals,nx,q)
  ! ==============================================================================
  ! Computes initial conditions for q fields using L2 projection onto basis
  ! INPUTS: meqn - number of fields to evaluate
  !         nx - number of points to evaluate q at
  !
  ! OUTPUTS: q(i,neq) - initial conditions evaluated at xvals(i) for neqth field
  ! ==============================================================================
  USE testParameters, ONLY: meqn,PI,testID
  USE modalDGmod
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nx
  DOUBLE PRECISION, DIMENSION(1:nx), INTENT(IN) :: xvals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:meqn), INTENT(OUT) :: q
  ! Local variables
  INTEGER :: k
  Double Precision :: r1

  SELECT CASE(testID)
  CASE(1)
    q(:,1) = exp(-100.0d0*(xvals(:)-0.5d0)*(xvals(:)-0.5d0))
!    q(:,2) = SIN(6D0*PI*xvals(:))+SIN(8D0*PI*xvals(:))
  CASE(2)
    do k=1,nx
       r1=min(1.d0/0.15d0*abs(xvals(k)-0.5d0),1.d0)
       q(k,1) = 0.5d0*(1D0+cos(PI*r1))
    enddo
  CASE(3)
    do k=1,nx
      q(k,1) = 0.d0 
      if ((xvals(k)> 0.25d0).and.(xvals(k)<0.75d0)) then
         q(k,1)=1.d0
      endif
    enddo
  END SELECT

END SUBROUTINE qinit
