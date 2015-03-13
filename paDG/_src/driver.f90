! ================================================================
! Polynomial adaptive discontinuous Galerkin (DG)
! By: Devin Light & Scott Moe
! For AMATH 574 -- March 2015
! =================================================================

SUBROUTINE DRIVER(ntest,nex0,nscale,nruns,noutput,maxCFL)
! ===============================================================
! Main driver subroutine for DG simulations
! Inputs:
!   ntest  : which test is being run
!   nex0   : number of initial spatial cells
!   nscale : multiple of nex0 used for subsequent runs (nex = nex0*nscale**p)
!   nruns  : number of total runs to make
!   maxCFL : maximal CFL number to use throughout integration
! ===============================================================
  USE testParameters
  USE gridData
  USE modalDGmod

  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: ntest,nex0,nscale,nruns,noutput
  DOUBLE PRECISION, INTENT(IN) :: maxCFL
  ! Ouputs
  ! Local variables
  INTEGER :: nxi,nsteps
  DOUBLE PRECISION :: dxi,time,dt,tnext,dtout,dti,dtupdate,tupdate
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: lam
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: quadNodes,quadWeights,&
    elemCenter,dxel,xOut,dxOut,xivals
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: &
    legendreQuad,dlegendreQuad
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: &
    qOut,q0,localCoeffs
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: coeffs,coeffsPrelim
  INTEGER, DIMENSION(:), ALLOCATABLE :: localDegree

  DOUBLE PRECISION, DIMENSION(1:nRuns,1:meqn) :: e1,e2,ei,cons,tmpqMax,tmpqMin
  REAL(KIND=4), DIMENSION(1:nRuns) :: t0,tf
  REAL(KIND=4), DIMENSION(2) :: tstart,tend
  DOUBLE PRECISION, DIMENSION(1:meqn) :: qvalsl
  DOUBLE PRECISION, DIMENSION(1) :: xpt

  INTEGER :: i,j,k,l,p,iout,m
  INTEGER :: ierr

  INTERFACE
    SUBROUTINE outputData(qvals,xvals,elCent,localPolyDegree,time,iout,p)
      USE testParameters, ONLY: meqn,testID,maxPolyDegree,outdir
      USE gridData
      IMPLICIT NONE
      ! Inputs
      DOUBLE PRECISION, DIMENSION(1:nxOut,1:meqn), INTENT(IN) :: qvals
      DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xvals
      DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elCent
      INTEGER, DIMENSION(1:nex), INTENT(IN) ::localPolyDegree
      DOUBLE PRECISION, INTENT(IN) :: time
      INTEGER, INTENT(IN) :: iout,p
    END SUBROUTINE outputData
    SUBROUTINE singleTimeStep(coeffs,Quadnodes,Quadwghts,legendreQuad,dlegendreQuad,dt,maxCFL,elemCenter,maxDegree,dx)
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
      DOUBLE PRECISION, INTENT(IN) :: dt
      DOUBLE PRECISION, INTENT(IN) :: maxCFL
      INTEGER, DIMENSION(1:nex), INTENT(IN) :: maxDegree
      ! Local valriables
      INTEGER :: i,j,k,m
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn) :: k1,k2,k3
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn) :: coeffsold
    END SUBROUTINE singleTimeStep

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
      DOUBLE PRECISION :: q,qp
    End SUBROUTINE updateCoeffs

     FUNCTION setWaveSpeed(npts,xivals,qvals)
        USE testParameters, ONLY: meqn,u0,K0,rho0
        IMPLICIT NONE
        !Inputs
        INTEGER, INTENT(IN) :: npts
        DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xivals
        DOUBLE PRECISION, DIMENSION(1:npts,1:meqn), INTENT(IN) :: qvals
        Double Precision, Dimension(1:npts) :: setWaveSpeed
        Integer :: i
     End Function setWaveSpeed

  END INTERFACE

  if(nRuns.lt.1) STOP 'nRuns should be at least 1 in DRIVER()'

  PI = DACOS(-1D0)

  ! Set number of quadrature nodes per element
  nQuad = maxPolyDegree
  ! Set number of local output nodes per element
  nxi = maxPolyDegree+1
  ALLOCATE(quadNodes(0:nQuad),quadWeights(0:nQuad),&
          legendreQuad(0:maxPolyDegree,0:nQuad),&
          dlegendreQuad(0:maxPolyDegree,0:nQuad), &
          xivals(1:nxi),STAT=ierr)

  ! Precompute Gauss-Legendre quadrature nodes and weights
  WRITE(*,*) 'Precomputing quadrature nodes and weights...'
  CALL gaussQuadNodes(nQuad+1,quadNodes)
  CALL gaussQuadWeights(nQuad+1,quadNodes,quadWeights)

  ! Precompute Legendre polynomial and its derivative values at quadrature locations
  DO i=0,nQuad
    DO k=0,maxPolyDegree
      legendreQuad(k,i) = legendre(quadNodes(i),k)
      dlegendreQuad(k,i) = dlegendre(quadNodes(i),k)
    ENDDO !k
  ENDDO !i

  ! Precompute local plotting points
  dxi = 2D0/DBLE(nxi)
  xivals(1) = -1D0+0.5D0*dxi
  DO i = 2,nxi
    xivals(i) = xivals(i-1)+dxi
  ENDDO !i

  WRITE(*,*) 'Beginning runs...'
  DO p=1,nRuns
    t0(p) = ETIME(tstart)

    ! Set up physical grids
    nex = nex0*nscale**(p-1) ! Number of elements in this run

    nxOut = nxi*nex ! Total number of plotting points

    ALLOCATE(elemCenter(1:nex),dxel(1:nex),&
             coeffs(0:maxPolyDegree,1:nex,1:meqn),&
             localDegree(1:nex),&
             lam(1:nex,1),&
             coeffsPrelim(0:maxPolyDegree,1:nex,1:meqn),&
             xOut(1:nxOut),dxOut(1:nxOut),qOut(1:nxOut,1:meqn),q0(1:nxOut,1:meqn),&
             localCoeffs(0:maxPolyDegree,1:meqn),STAT=ierr)

    ! Make element grid
    CALL makeGrid(nex,dxel,elemCenter)
    ! Make output grid
    CALL makeGrid(nxOut,dxOut,xOut)

    ! Initialize q fields to maxPolyDegree accuracy
    CALL computeICs(elemCenter,localDegree,dxel,quadNodes,quadWeights,legendreQuad,coeffs)
    ! Update the local degree and a buffer region
    CALL updateCoeffs(coeffs,localDegree)
    CALL updateBuffer(coeffs,localDegree)

    ! Evalutate ICs at plotting points and also find wave speeds at cell edges
    DO j=1,nex
      localCoeffs = coeffs(0:maxPolyDegree,j,:)
      CALL evaluateExpansion(localCoeffs,xivals,qOut(1+(j-1)*nxi:j*nxi,:),nxi,maxPolyDegree)

      xpt(1)=-1.d0
      CALL evaluateExpansion(localCoeffs,xpt,qvalsl,1,localDegree(j))
      xpt(1)=-0.5d0*dxel(j)+elemCenter(j)
      lam(j,:)=setWaveSpeed(1,xpt,qvalsl)

    ENDDO !j
    q0 = qout

    ! Set qMax and qMin values
    DO m=1,meqn
      tmpqMax(p,m) = MAXVAL(q0(:,m))
      tmpqMin(p,m) = MINVAL(q0(:,m))
    ENDDO

    iout = 0
    time = 0D0
    CALL outputData(qOut,xOut,elemCenter,localDegree,0D0,iout,p)
    iout = iout+1

    IF(noutput .eq. -1) THEN
      WRITE(*,*) '=== in DRIVER: noutput should be >= 0'
      STOP
    ELSE
      nsteps = noutput*CEILING( &
              (tfinal/maxCFL)*(MAXVAL(ABS(lam(:,1)))/MINVAL(dxel))/DBLE(noutput) )
    ENDIF
    dt = tfinal/DBLE(nsteps)
    dti = tfinal/DBLE(nsteps)
    dtout = tfinal/DBLE(noutput)
    tnext = dtout

    dtupdate=(MINVAL(dxel)/MAXVAL(ABS(lam(:,1))))

    tupdate=dtupdate

    WRITE(*,*) 'Beginning time stepping...'
    DO WHILE (time - tfinal<-1.0e-8)
      dt=dti
      coeffsPrelim=coeffs
      call singleTimeStep(coeffsPrelim,Quadnodes,quadWeights,&
        legendreQuad,dlegendreQuad,dt,maxCFL,elemCenter,localDegree,dxel)
      ! Check if this is an output time
      IF((time<tnext).AND.((time+dt-tnext)>-1.d-8)) THEN
        dt=tnext-time;
        ! Redo-run to match output time
        coeffsPrelim=coeffs
        call singleTimeStep(coeffsPrelim,Quadnodes,quadWeights,&
          legendreQuad,dlegendreQuad,dt,maxCFL,elemCenter,localDegree,dxel)
        ! Evalutate approximate solution at this time
        DO j=1,nex
          localCoeffs = coeffsPrelim(0:maxPolyDegree,j,:)
          CALL evaluateExpansion(localCoeffs,xivals,qOut(1+(j-1)*nxi:j*nxi,:),nxi,maxPolyDegree)
        ENDDO !j

        ! Set qMax and qMin values
        DO m=1,meqn
          tmpqMax(p,m) = MAX(tmpqMax(p,m),MAXVAL(qOut(:,m)) )
          tmpqMin(p,m) = MIN(tmpqMin(p,m),MINVAL(qOut(:,m)) )
        ENDDO

        CALL outputData(qOut,xOut,elemCenter,localDegree,time+dt,iout,p)
        iout = iout+1
        tnext=tnext+dtout
      ENDIF
      coeffs=coeffsPrelim

      !print *, 'Super Awesome DG code dt=',dt,'time=',time
      IF((time<tupdate).AND.((time+dt-tupdate)>-1.d-8)) THEN
         ! Update local degree and buffer regions
         CALL updateCoeffs(coeffs,localDegree)
         CALL updateBuffer(coeffs,localDegree)

         ! Evalutate ICs at plotting points and also find wave speeds at cell edges
         DO j=1,nex
            localCoeffs = coeffs(0:maxPolyDegree,j,:)
            xpt(1)=-1.d0
            CALL evaluateExpansion(localCoeffs,xpt,qvalsl,1,maxPolyDegree)
            xpt(1)=-0.5d0*dxel(j)+elemCenter(j)
            lam(j,:)=setWaveSpeed(1,xpt,qvalsl)
         ENDDO !j

         dtupdate=(MINVAL(dxel)/MAXVAL(ABS(lam(:,1))))
         tupdate=time+dtupdate
         !print *,time,'maxdegree=',maxDegree
      ENDIF
      time = time + dt
    ENDDO
    tf(p) = ETIME(tend)-t0(p)
    CALL computeErrors(qOut,q0,e1,e2,ei,cons,tmpqMax,tmpqMin,tf,nxOut,nRuns,nex0,nscale,p)
    !print *,time,'maxdegree=',maxDegree

    DEALLOCATE(elemCenter,localDegree,lam,dxel,coeffs,xOut,dxOut,qOut,q0,localCoeffs,STAT=ierr)
  ENDDO !p

  CALL computeErrors(qOut,q0,e1,e2,ei,cons,tmpqMax,tmpqMin,tf,nxOut,nRuns,nex0,nscale,-1)
  DEALLOCATE(quadNodes,quadWeights,legendreQuad,dlegendreQuad,xivals,STAT=ierr)

CONTAINS
  SUBROUTINE computeErrors(qOut,q0,e1,e2,ei,cons,ovrshoot,undrshoot,tf,nxOut,nRuns,nex0,nscale,stat)
    ! =============================================================================
    ! Prints error estimates and other useful information to screen
    ! INPUTS: qOut - current estimate solution
    !         q0   - initial conditions
    !         tf(p) - cput time for pth run
    !         stat - status integer
    !         ovrshoot(p,m) - maximum overshoot in appx soln at plotting times
    !         undrshoot(p,m) - maximum undershoot in appx soln at plotting times
    ! OUTPUTS: e1(p,m) - L1 error estimate
    !          e2(p,m) - L2 error estimate
    !          ei(p,m) - Linf error estimate
    !          cons(p,m) - conservation estimate
    ! =============================================================================
    USE testParameters
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nRuns,nxOut,stat,nscale,nex0
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:meqn), INTENT(IN) :: q0,qOut
    DOUBLE PRECISION, DIMENSION(1:nRuns,1:meqn), INTENT(IN) :: ovrshoot,undrshoot
    REAL(KIND=4), DIMENSION(1:nRuns),INTENT(IN) :: tf
    ! Outputs
    DOUBLE PRECISION, DIMENSION(1:nRuns,1:meqn),INTENT(INOUT) :: e1,e2,ei,cons
    ! Local variables
    INTEGER :: p,m
    CHARACTER(len=2) :: qname
    DOUBLE PRECISION :: cnvg1,cnvg2,cnvgi

    IF(stat == -1) THEN
      ! Write error output to screen
      DO m=1,meqn
        WRITE(qname,'(a,i1)') 'q',m
        WRITE(*,*) '===================='
        WRITE(*,'(a12)') qname
        WRITE(*,*) '===================='
      WRITE(*,'(A107)') &
'nex      E1        E2       Einf      convergence rate  maximum   minimum       cons       cputime   tf'

        cnvg1 = 0D0
        cnvg2 = 0D0
        cnvgi = 0D0

        DO p=1,nRuns
          IF(p.gt.1) THEN
            cnvg1 = -log(e1(p,m)/e1(p-1,m))/log(dble(nscale))
            cnvg2 = -log(e2(p,m)/e2(p-1,m))/log(dble(nscale))
            cnvgi = -log(ei(p,m)/ei(p-1,m))/log(dble(nscale))
          ENDIF

          WRITE(*,990) nex0*nscale**(p-1), e1(p,m), e2(p,m), ei(p,m), &
                cnvg1, cnvg2, cnvgi, &
                ovrshoot(p,m), &
                undrshoot(p,m), &
                cons(p,m), tf(p),tfinal
        ENDDO!p
      ENDDO !m
      990    format(i6,3e12.4,3f5.2,3e12.4,2f8.2)
    ELSE
      ! Compute error estimates for this run
      DO m=1,meqn
        e1(stat,m) = SUM(ABS(qOut(:,m)-q0(:,m)))/DBLE(nxOut)
  			e2(stat,m) = SQRT(SUM((qOut(:,m)-q0(:,m))**2)/DBLE(nxOut))
  			ei(stat,m) = MAXVAL(ABS( qOut(:,m)-q0(:,m) ))
        cons(stat,m) = SUM(qOut(:,m)-q0(:,m))/DBLE(nxOut)
      ENDDO !m

    ENDIF

  END SUBROUTINE computeErrors

  SUBROUTINE makeGrid(nex,dxel,xCenter)
    ! =============================================================================
    ! Computes cell width and initializes cell centers and quadrature grid
    ! INPUTS:   nex - number of elements
    !           nQuad - number of quadrature nodes
    !           quadNodes(0:nQuad) - Gauss-Legendre quadrature nodes
    ! OUTPUTS:  dxel - width of elements
    !           xCenter(j) - location of jth element center
    ! =============================================================================
    USE gridData, ONLY: xDomain
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nex
    ! Outputs
    DOUBLE PRECISION, DIMENSION(1:nex), INTENT(OUT) :: dxel,xCenter

    ! Local variables
    DOUBLE PRECISION :: domainWidth
    INTEGER :: k,j

    domainWidth = xDomain(2)-xDomain(1)
    dxel(:) = domainWidth/DBLE(nex)

    xCenter(1) = xDomain(1)+0.5D0*dxel(1)
    DO j=2,nex
      xCenter(j) = xCenter(j-1)+0.5D0*(dxel(j-1)+dxel(j))
    ENDDO!j

  END SUBROUTINE makeGrid

  SUBROUTINE computeICs(elemCenter,maxDegree,dx,quadNodes,quadWeights,legendreQuad,coeffs)
    ! =============================================================================
    ! Computes initial coefficient values via L2 projection onto basis
    ! INPUTS:   maxPolyDegree - maximum degree of basis polynomial
    !           meqn - number of equations in system
    !           nex - number of elements
    !           nQuad - number of quadrature nodes
    !           quadNodes(0:nQuad) - Gauss-Legendre quadrature nodes
    !           quadWeights(0:nQuad) - Gauss-Legendre quadrature weights
    !           legendreQuad(k,i) - kth degree legendre polynomial at ith quadrature node
    !           dxel(j) - width of jth element
    !           elemCenter(j) - location of jth element center
    ! OUTPUTS:  coeffs(k,j,neq) - kth degree coefficient in jth element for neqth equation
    ! =============================================================================
    USE testParameters, ONLY: maxPolyDegree,meqn
    USE gridData, ONLY: nex,nQuad!,projectL2mod

    IMPLICIT NONE
    ! Inputs
    DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemCenter,dx
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadNodes,quadWeights
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreQuad
    ! Outputs
    INTEGER, DIMENSION(1:nex), INTENT(INOUT) :: maxDegree
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nex,1:meqn), INTENT(OUT) :: coeffs
    ! Local variables
    INTEGER :: i,j,k,m,D,Dm
    DOUBLE PRECISION, DIMENSION(0:nQuad) :: xvals
    DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn) :: qvals
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn) :: tmpCoeffs
    DOUBLE PRECISION :: q,qp,TOL

    INTERFACE
      SUBROUTINE qinit(xvals,nx,q)
        USE testParameters, ONLY: meqn,maxPolyDegree
        USE gridData, ONLY: nQuad
        INTEGER, INTENT(IN) :: nx
        DOUBLE PRECISION, DIMENSION(1:nx),INTENT(IN) :: xvals
        ! Outputs
        DOUBLE PRECISION, DIMENSION(1:nx,1:meqn),INTENT(OUT) :: q
      END SUBROUTINE qinit

      SUBROUTINE projectL2(qvals,quadWeights,legendreQuad,coeffs)
        USE testParameters, ONLY: meqn,maxPolyDegree
        USE gridData, ONLY: nQuad
        DOUBLE PRECISION, DIMENSION(0:nQuad),INTENT(IN) :: quadWeights
        DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn),INTENT(IN) :: qvals
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad),INTENT(IN) :: legendreQuad
        !Outputs
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree), INTENT(OUT) :: coeffs
      END SUBROUTINE projectL2
    END INTERFACE

    TOL=1.0d-3

    DO j=1,nex
      xvals(:) = 0.5D0*quadNodes(:)*dx(j)+elemCenter(j)

      CALL qinit(xvals,nQuad+1,qvals)
      !print *,'i=', j
      CALL projectL2(qvals,quadWeights,legendreQuad,tmpCoeffs)

      D=0
      OUTER: Do m=1,meqn
      Dm=maxPolyDegree
        q=tmpCoeffs(1,m)
        qp=tmpCoeffs(1,m)
        INNER: Do k=0,maxPolyDegree
          qp=q
          q=tmpCoeffs(k,m)
          IF (((abs(q)<TOL).and.(abs(qp)<TOL)).and.(k>0)) THEN
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
      !print *,j,D, tmpCoeffs
      coeffs(:,j,:) = tmpCoeffs(:,:)
    ENDDO !j

  END SUBROUTINE computeICs

END SUBROUTINE DRIVER
