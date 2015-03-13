SUBROUTINE updateCoeffs(coeffs,maxDegree)
  USE testParameters, ONLY: maxPolyDegree,meqn
  USE gridData, ONLY: nex,nQuad!,projectL2mod

  IMPLICIT NONE
  ! Inputs and Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(INOUT) :: coeffs
  INTEGER, DIMENSION(1:nex), INTENT(INOUT) :: maxDegree
  ! Local valriables
  INTEGER :: i,j,k,m,D,Dm
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn) :: coeffsold,coeffsnew
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn) :: tmpCoeffs
  DOUBLE PRECISION :: q,qp, TOL

  TOL=1.0d-8

  DO j=1,nex
      tmpCoeffs(:,:) = coeffs(:,j,:)
      D=0
      OUTER: Do m=1,meqn
      Dm=maxPolyDegree
        q=tmpCoeffs(0,m)
        qp=tmpCoeffs(0,m)
        INNER: Do k=1,maxPolyDegree
          qp=q
          q=tmpCoeffs(k,m)
          IF (((abs(q)<TOL).and.(abs(qp)<TOL)).and.(k>1)) THEN
             Dm=k
             D=max(Dm,D)
             EXIT INNER
          ENDIF
          IF (k>maxPolyDegree-1) THEN
             D=maxPolyDegree
          ENDIF
        ENDDO INNER
      ENDDO OUTER
      maxDegree(j)=D
      Do k=D+1,maxPolyDegree
        coeffs(k,j,:)=0.d0
      EndDo
  ENDDO !j

End SUBROUTINE updateCoeffs
