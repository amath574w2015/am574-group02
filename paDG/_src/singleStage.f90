SUBROUTINE singleStage(coeffs,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,&
            istage,dt,maxCFL,elemCenter,maxDegree,dx,k1,computeCFL)
!=========================================!
!This is a subroutine that evaluates a whole
!time step of the discontinuous Galerkin 
!scheme. This version of the file is written
!to be used with RK3. However other time-stepping
!schemes can use it
!=========================================!

  USE testParameters, ONLY: maxPolyDegree,meqn
  USE gridData, ONLY: nex,nQuad!,projectL2mod

  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(IN) :: coeffs
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: Quadnodes
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: Quadwghts
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreQuad
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: dlegendreQuad
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemCenter
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: dx
  INTEGER, INTENT(IN) :: istage
  DOUBLE PRECISION, INTENT(IN) :: maxCFL
  INTEGER, DIMENSION(1:nex), INTENT(IN) :: maxDegree
  ! Output
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(INOUT) :: k1
  DOUBLE PRECISION, INTENT(INOUT) :: dt
  DOUBLE PRECISION, INTENT(INOUT) :: computeCFL
  ! Output
  ! Local valriables
  INTEGER :: i,j,k,m,degl,deg,degr
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn) :: fluxQuad
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: localCoeffs,qOut,numFlux,numFlux1
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xivals
  DOUBLE PRECISION, DIMENSION(1:meqn) :: qvalsl,qvalsr
  DOUBLE PRECISION, DIMENSION(1) :: xpt,s1
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn) :: coeffL,coeffR
  DOUBLE PRECISION :: dtmin


  INTERFACE

    FUNCTION RiemannSolve(qvalsl,qvalsr,xpt)
      USE testParameters, ONLY: meqn,u0,K0,rho0
      IMPLICIT NONE
      ! Inputs
      DOUBLE PRECISION, DIMENSION(1), INTENT(IN) :: xpt
      DOUBLE PRECISION, DIMENSION(1:meqn), INTENT(IN) :: qvalsl,qvalsr
      ! Outputs
      DOUBLE PRECISION, DIMENSION(0:1,1:meqn) :: RiemannSolve
      ! Local Variables
      Double Precision, Dimension(1,1:meqn) :: Fl,Fr
      Double Precision, Dimension(1) :: s1,s2,smax
    End FUNCTION RiemannSolve

    SUBROUTINE evaluateExpansion(coeffs,xivals,qvals,nxivals,maxDegree)
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
    End SUBROUTINE evaluateExpansion

    FUNCTION forcingFunction(k,fluxQuad,quadWeights,dlegendreQuad,numFlux)
      USE testParameters, ONLY: meqn
      USE gridData, ONLY: nQuad
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: k
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights,dlegendreQuad
      DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn), INTENT(IN) :: fluxQuad
      DOUBLE PRECISION, DIMENSION(0:1,1:meqn), INTENT(IN) :: numFlux
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:meqn) :: forcingFunction
      ! Local variables
      INTEGER :: p,m
      DOUBLE PRECISION :: coeff
    End FUNCTION forcingFunction

    FUNCTION fluxFunction(qvals,xpts,npts)
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
    End FUNCTION fluxFunction

     FUNCTION setWaveSpeed(npts,xivals,qvals)
        USE testParameters, ONLY: meqn,u0,K0,rho0
        IMPLICIT NONE
        !Inputs
        INTEGER, INTENT(IN) :: npts
        DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xivals
        DOUBLE PRECISION, DIMENSION(1:npts,1:meqn), INTENT(IN) :: qvals
        !Outputs
        Double Precision, Dimension(1:npts) :: setWaveSpeed
        Integer :: i
     End Function setWaveSpeed

  End INTERFACE

  ALLOCATE(localCoeffs(0:maxPolyDegree,1:meqn),qOut(1:nQuad+1,1:meqn),xivals(1:nQuad+1),numFlux1(0:1,1:meqn),numFlux(0:1,1:meqn))

  localCoeffs(0:maxPolyDegree,1:meqn)=0D0
  coeffL(0:maxPolyDegree,1:meqn)=0D0
  coeffR(0:maxPolyDegree,1:meqn)=0D0
 

  dtmin=100.d0
  degl=maxDegree(nex)
  deg=maxDegree(1)
  degr=maxDegree(2)
  coeffL(0:degl,:)=coeffs(0:degl,nex,:)
  DO j=1,nex-1
          degr=maxDegree(j+1)
          localCoeffs(0:deg,:) = coeffs(0:deg,j,:)
          xivals(:)=0.5D0*quadNodes(:)*dx(j)+elemCenter(j)
          CALL evaluateExpansion(localCoeffs,quadNodes(:),qOut,nQuad+1,deg)
          fluxQuad(:,:)=fluxFunction(qOut,xivals,nQuad+1)


          xpt(1)=1.d0
          CALL evaluateExpansion(coeffL,xpt,qvalsl,1,degl)
          xpt(1)=-1.d0
          CALL evaluateExpansion(localCoeffs,xpt,qvalsr,1,deg)
          xpt(1)=-0.5d0*dx(j)+elemCenter(j)
          numFlux1(:,:)=RiemannSolve(qvalsl,qvalsr,xpt)
          numFlux(0,:)=numFlux1(1,:)


          xpt(1)=-1.d0
          coeffR(0:degr,:)=coeffs(0:degr,j+1,:)
          CALL evaluateExpansion(coeffR,xpt,qvalsr,1,degr)
          xpt(1)=1.d0
          CALL evaluateExpansion(localCoeffs,xpt,qvalsl,1,deg)
          xpt(1)=0.5d0*dx(j)+elemCenter(j)
          numFlux1(:,:)=RiemannSolve(qvalsl,qvalsr,xpt)
          numFlux(1,:)=numFlux1(0,:)


          s1=setWaveSpeed(1,xpt,qvalsl)
          Do k=0,deg
              k1(k,j,:)=-1.0d0/dx(j)*forcingFunction(k,fluxQuad,Quadwghts,dlegendreQuad(k,:),numFlux)
          EndDo !k
        dtmin=min(dtmin,maxCFL/s1(1)*dx(j))
        computeCFL=max(computeCFL,s1(1)*dt/dx(j))
	degl=deg
	deg=degr
        coeffL(0:degl,:)=localCoeffs(0:degl,:)
	ENDDO !j

	j=nex
	degr=maxDegree(1)
coeffR(0:degr,:)=coeffs(0:degr,1,:)

	localCoeffs(0:deg,:) = coeffs(0:deg,j,:)
	xivals(:)=0.5D0*quadNodes(:)*dx(j)+elemCenter(j)
	CALL evaluateExpansion(localCoeffs,quadNodes(:),qOut,nQuad+1,deg)
fluxQuad(:,:)=fluxFunction(qOut,xivals,nQuad+1)


	xpt(1)=1.d0
CALL evaluateExpansion(coeffL,xpt,qvalsl,1,degl)
	xpt(1)=-1.d0
	CALL evaluateExpansion(localCoeffs,xpt,qvalsr,1,deg)
	xpt(1)=-0.5d0*dx(j)+elemCenter(j)
	numFlux1(:,:)=RiemannSolve(qvalsl,qvalsr,xpt)
numFlux(0,:)=numFlux1(1,:)


	xpt(1)=-1.d0
CALL evaluateExpansion(coeffR,xpt,qvalsr,1,degr)
	xpt(1)=1.d0
	CALL evaluateExpansion(localCoeffs,xpt,qvalsl,1,deg)
	xpt(1)=0.5d0*dx(j)+elemCenter(j)
	numFlux1(:,:)=RiemannSolve(qvalsl,qvalsr,xpt)
numFlux(1,:)=numFlux1(0,:)


s1=setWaveSpeed(1,xpt,qvalsl)
  Do k=0,deg
    k1(k,j,:)=-1.0d0/dx(j)*forcingFunction(k,fluxQuad,Quadwghts,dlegendreQuad(k,:),numFlux)
  EndDo !k
  dtmin=min(dtmin,maxCFL/s1(1)*dx(j))
  computeCFL=max(computeCFL,s1(1)*dt/dx(j))

  IF(istage.eq.1) THEN
    dt=min(dt,dtmin)
    computeCFL=maxCFL
  ENDIF
  !k1=k1*dt
End SUBROUTINE singleStage
