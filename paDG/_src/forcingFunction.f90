FUNCTION forcingFunction(k,fluxQuad,quadWeights,dlegendreQuad,numFlux,numQuad)
  ! ========================================================================
  ! Computes the right hand side forcing term for evolution of kth coefficient
  ! INPUTS: k - which degree coefficient is being updated
  !         fluxQuad(:,m) - exact flux function for equation m eval'd at quad nodes
  !         quadWeights(:) - quadrature weights
  !         dlegendreQuad(:) - derivative of kth basis function eval'd at quad nodes
  !         numFlux(0,m) - numerical flux for mth eqn through left interface
  !         numFlux(1,m) - numerical flux for mth eqn through right interface
  ! OUTPUT:
  !         forcingFunction(m) - RHS forcing for evolution of kth coefficient
  ! ========================================================================
  USE testParameters, ONLY: meqn
  USE gridData, ONLY: nQuad
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: k
  INTEGER, INTENT(IN) :: numQuad
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights,dlegendreQuad
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn), INTENT(IN) :: fluxQuad
  DOUBLE PRECISION, DIMENSION(0:1,1:meqn), INTENT(IN) :: numFlux
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:meqn) :: forcingFunction
  ! Local variables
  INTEGER :: p,m,l
  DOUBLE PRECISION :: coeff

  forcingFunction(:)=0.d0

  coeff = 2D0*k+1D0
  DO m=1,meqn
    DO l=0,nQuad
      forcingFunction(m) = forcingFunction(m)+(quadWeights(l)*dlegendreQuad(l)*fluxQuad(l,m))
    ENDDO !l
  ENDDO !m
  forcingFunction(:) = -forcingFunction(:) + (numFlux(1,:)-(-1D0)**k*numFlux(0,:))
  forcingFunction(:) = coeff*forcingFunction(:)

END FUNCTION forcingFunction
