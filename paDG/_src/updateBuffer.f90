SUBROUTINE updateBuffer(coeffs,maxDegree)
  USE testParameters, ONLY: maxPolyDegree,meqn
  USE gridData, ONLY: nex,nQuad!,projectL2mod

  IMPLICIT NONE
  ! Inputs and Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(INOUT) :: coeffs
  INTEGER, DIMENSION(1:nex), INTENT(INOUT) :: maxDegree
  ! Local valriables
  INTEGER :: i,j,k,m,D
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn) :: coeffsold,coeffsnew
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn) :: tmpCoeffs
  DOUBLE PRECISION :: q,qp
  INTEGER, DIMENSION(1:nex) :: maxDegreetmp

D=0
D=max(D,maxDegree(nex))
D=max(D,maxDegree(1))
D=max(D,maxDegree(2))
maxDegreetmp(1)=D

  DO j=2,nex-1
     D=0
     DO k=-1,1
       D=max(D,maxDegree(j+k))
     ENDDO
     maxDegreetmp(j)=D
  EndDo

D=0
D=max(D,maxDegree(nex))
D=max(D,maxDegree(1))
D=max(D,maxDegree(nex-1))
maxDegreetmp(nex)=D

  DO j=1,nex
     maxDegree(j)=maxDegreetmp(j)
  ENDDO

END SUBROUTINE updateBuffer
   
   

