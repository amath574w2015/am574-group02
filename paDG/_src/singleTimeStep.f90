SUBROUTINE singleTimeStep(coeffs,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,dt,maxCFL,elemCenter,maxDegree,dx)
!=========================================!
!This is a subroutine that evaluates a whole
!time step of the discontinuous Galerkin 
!scheme. This version of the file is written
!to use RK3. However other time-stepping
!schemes will use similar files

!In this case we will use this RK3 
!k1 = h f(xi, yi),
!k2 = h f(xi + h / 2, yi + k1 / 2 ),
!k3 = h f(xi + h, yi - k1 + 2 k2 )
!y_{k+1}=y_k+1/6*(k1+4*k2+k3) 
!=========================================!

  USE testParameters, ONLY: maxPolyDegree,meqn
  USE gridData, ONLY: nex,nQuad!,projectL2mod

  IMPLICIT NONE
  ! Inputs and Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(INOUT) :: coeffs
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: Quadnodes
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: Quadwghts
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreQuad
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: dlegendreQuad
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemCenter
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: dx
  DOUBLE PRECISION, INTENT(INOUT) :: dt
  DOUBLE PRECISION, INTENT(IN) :: maxCFL
  INTEGER, DIMENSION(1:nex), INTENT(IN) :: maxDegree
  ! Local valriables
  INTEGER :: i,j,k,m
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn) :: k1,k2,k3
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn) :: coeffsold,coeffsnew
  DOUBLE PRECISION :: computeCFL

  INTERFACE

  SUBROUTINE  singleStage(coeffs,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,&
              istage,dt,maxCFL,elemCenter,maxDegree,dx,k1,computeCFL)
      USE testParameters, ONLY: maxPolyDegree,meqn
      USE gridData, ONLY: nex,nQuad!,projectL2mod
      IMPLICIT NONE
      ! Inputs and Outputs
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(INOUT) :: coeffs
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: Quadnodes
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: Quadwghts
      DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemCenter,dx
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(INOUT) :: k1
      DOUBLE PRECISION, INTENT(IN) :: dt,computeCFL
      INTEGER, INTENT(IN) :: istage
      DOUBLE PRECISION, INTENT(IN) :: maxCFL
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreQuad,dlegendreQuad
      INTEGER, DIMENSION(1:nex), INTENT(IN) :: maxDegree
      ! Local valriables
      INTEGER :: i,j,k,m
      DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn) :: fluxQuad
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: localCoeffs,qOut,numFlux,numFlux1
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xivals
      DOUBLE PRECISION, DIMENSION(1:meqn) :: qvalsl,qvalsr
      DOUBLE PRECISION, DIMENSION(1) :: xpt
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn) :: coeffL,coeffR
    End SUBROUTINE singleStage

  END INTERFACE
  k1(:,:,:)=0.d0
  k2(:,:,:)=0.d0
  k3(:,:,:)=0.d0
  coeffsnew(:,:,:)=0.d0

  computeCFL=0.d0
  coeffsold=coeffs
  call singleStage(coeffsold,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,1,&
       dt,maxCFL,elemCenter,maxDegree,dx,k1,computeCFL)

  DO j=1,nex
    Do m=1,meqn
      Do k=0,maxDegree(j)
        coeffsnew(k,j,m)=coeffsold(k,j,m)+0.5d0*dt*k1(k,j,m)
      ENDDO
    ENDDO
  ENDDO
  !coeffsnew=coeffsold+0.5*k1
  call singleStage(coeffsnew,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,2,&
        dt,maxCFL,elemCenter,maxDegree,dx,k2,computeCFL)
  DO j=1,nex
    Do m=1,meqn
      Do k=0,maxDegree(j)
        coeffsnew(k,j,m)=coeffsold(k,j,m)-dt*k1(k,j,m)+2.d0*dt*k2(k,j,m)
      ENDDO
    ENDDO
  ENDDO
  !coeffsnew=coeffsold-k1+2.d0*k2
  call singleStage(coeffsnew,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,3,&
        dt,maxCFL,elemCenter,maxDegree,dx,k3,computeCFL)
  if (computeCFL>maxCFL) then
     print *,'cfl3',computeCFL,maxCFL
  endif
  DO j=1,nex
    Do m=1,meqn
      Do k=0,maxDegree(j)
        coeffs(k,j,m)=coeffsold(k,j,m)+1.d0/6.d0*dt*(k1(k,j,m)+4.d0*k2(k,j,m)+k3(k,j,m))
      ENDDO
    ENDDO
  ENDDO
  !coeffs=coeffsold+1.d0/6.d0*(k1+4.d0*k2+k3)

End SUBROUTINE singleTimeStep
