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

  Do i=1,npts
     setWaveSpeed(i)=max(abs(u0-SQRT(K0/rho0)),abs(u0+SQRT(K0/rho0)))
  EndDo

End FUNCTION setWaveSpeed
