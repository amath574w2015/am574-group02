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
  
  Fl=fluxFunction(qvalsl,xpt,1)
  Fr=fluxFunction(qvalsr,xpt,1)
  
  s1=setWaveSpeed(1,xpt,qvalsl)
  s2=setWaveSpeed(1,xpt,qvalsr)
  
  smax(1)=max(abs(s1(1)),abs(s2(1)))  
    
  RiemannSolve(0,:)=0.5d0*((Fl(1,:)+Fr(1,:))+smax(1)*(qvalsl(1,:)-qvalsr(1,:)))
  RiemannSolve(1,:)=0.5d0*((Fl(1,:)+Fr(1,:))+smax(1)*(qvalsl(1,:)-qvalsr(1,:)))

End FUNCTION RiemannSolve
