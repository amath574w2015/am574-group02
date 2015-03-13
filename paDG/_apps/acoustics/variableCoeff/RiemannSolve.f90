FUNCTION RiemannSolve(qvalsl,qvalsr,xpt)
  !===============================!
  !This function evaluates an approximate
  !Local Lax Friedrichs flux. This is an
  !approximate Riemann solve 
  !===============================!
  USE testParameters, ONLY: meqn,u0,K0,rho0
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, DIMENSION(1), INTENT(IN) :: xpt
  DOUBLE PRECISION, DIMENSION(1,1:meqn), INTENT(IN) :: qvalsl,qvalsr
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:1,1:meqn) :: RiemannSolve
  ! Local Variables
  Double Precision, Dimension(1,1:meqn) :: Fl,Fr
  Double Precision, Dimension(1) :: s1,s2,smax
  Double Precision del0,del1,zr,zl,a0,a1,cgr,rhogr,cgl,rhogl,qm0,qm1
  Double Precision cr,rhor,cl,rhol
 
  !Interface blocks
  Interface

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
     End Function fluxFunction
     
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
  End Interface


  rhogl=1.0d0
  cgl=1.0d0
  rhogr=2.0d0
  cgr=0.5d0



    if (xpt(1)<-1.d-6) then
       rhol=rhogl
       rhor=rhogl
       cl=cgl
       cr=cgl
    else if (abs(xpt(1))<1.0d-6) then
       rhol=rhogl
       rhor=rhogr
       cl=cgl
       cr=cgr
    else 
       rhol=rhogr
       rhor=rhogr
       cl=cgr
       cr=cgr
    endif
  
  
  s1=setWaveSpeed(1,xpt,qvalsl)
  s2=setWaveSpeed(1,xpt,qvalsr)
  
  smax(1)=max(abs(s1(1)),abs(s2(1)))  

  del0=qvalsr(1,1)-qvalsl(1,1)
  del1=qvalsr(1,2)-qvalsl(1,2)

  zr = cr*rhor
  zl = cl*rhol

  a0 = (-del0 + zr*del1)/(zl+zr)
  a1 = (del0  + zl*del1)/(zl+zr)

  qm0=qvalsl(1,1)-a0*zl
  qm1=qvalsl(1,2)+a0
   
  RiemannSolve(0,1)=cl*cl*rhol*qm1
  RiemannSolve(0,2)=qm0/rhol
  RiemannSolve(1,1)=cr*cr*rhor*qm1
  RiemannSolve(1,2)=qm0/rhor

End FUNCTION RiemannSolve
