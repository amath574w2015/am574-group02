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


    do k=1,nx
      q(k,1) = 0.d0
      if ((xvals(k)> -4.d0).and.(xvals(k)<-2.d0)) then
         !q(k,1)=sqrt(1.d0-(xvals(k)+3.d0)*(xvals(k)+3.d0))
         q(k,1)=exp(-2.d0*(xvals(k)+3.d0)*(xvals(k)+3.d0))
      endif
         q(k,2) = q(k,1)
    enddo
END SUBROUTINE qinit
