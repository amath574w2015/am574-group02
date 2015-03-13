FUNCTION setWaveSpeed(npts) 

  !This function returns the wave speed (eigenvalues)
  !of the PDE at every point in space. This specific 
  !function is written for constant coefficient advection

  USE testParameters, ONLY: meqn,u0,K0,rho0
  IMPLICIT NONE
  !Inputs
  INTEGER, INTENT(IN) :: npts
  !Outputs
  Double Precision, Dimension(1:npts) :: setWaveSpeed
 
  Integer :: i
 
  Do i=1,npts
     setWaveSpeed(i)=u0
  EndDo

End FUNCTION setWaveSpeed

