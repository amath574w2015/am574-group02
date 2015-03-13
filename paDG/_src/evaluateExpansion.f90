SUBROUTINE evaluateExpansion(coeffs,xivals,qvals,nxivals,maxDegree)
  ! ===========================================================================
  ! Evaluates polynomial expansion phi = \sum coeffs_k * P_k at locations given by xivals
  ! Used for outputting solution values
  ! INPUTS: xivals(1:nxivals)
  !         coeffs(0:maxPolyDegree,1:meqn)
  !
  ! OUTPUTS: qvals
  ! ===========================================================================
  USE testParameters
  USE gridData
  USE modalDGmod
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nxivals
  INTEGER, INTENT(IN) :: maxDegree
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn), INTENT(IN) :: coeffs
  DOUBLE PRECISION, DIMENSION(1:nxivals), INTENT(IN) :: xivals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nxivals,1:meqn), INTENT(OUT) :: qvals
  ! Local valriables
  INTEGER :: i,k,m
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nxivals) :: legendreVals

  DO k=0,maxDegree
    DO i=1,nxivals
      legendreVals(k,i) = legendre(xivals(i),k)
    ENDDO!i
  ENDDO !k

  DO m=1,meqn
    DO i=1,nxivals
      qvals(i,m)=0.d0
      DO k=0,maxDegree
        qvals(i,m) = qvals(i,m)+coeffs(k,m)*legendreVals(k,i)
      ENDDO !k
    ENDDO !i
  ENDDO !m

END SUBROUTINE evaluateExpansion
