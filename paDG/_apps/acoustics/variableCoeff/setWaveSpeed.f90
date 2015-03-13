FUNCTION setWaveSpeed(npts,xivals,qvals)

  !This function returns the wave speed (eigenvalues)
  !of the PDE at every point in space. This specific
  !function is written for constant coefficient advection

  USE testParameters, ONLY: meqn,u0,K0,rho0
  IMPLICIT NONE
  !Inputs
  INTEGER, INTENT(IN) :: npts
  DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xivals
  DOUBLE PRECISION, DIMENSION(1:npts,1:meqn), INTENT(IN) :: qvals

  !Outputs
  DOUBLE PRECISION, DIMENSION(1:npts) :: setWaveSpeed
  Integer :: i
  DOUBLE PRECISION rhol,cl,rhor,cr,rho,c,K

  rhol=1.0d0
  cl=1.0d0
  rhor=2.0d0
  cr=0.5d0


  Do i=1,npts
    if (xivals(i)<0.d0) then
       rho=rhol
       c=cl
    else if (abs(xivals(i))<1.d-10) then
       rho=rhor
       c=max(cr,cl)
    else 
       rho=rhor
       c=cr
    endif

     setWaveSpeed(i)=c
  EndDo

End FUNCTION setWaveSpeed
