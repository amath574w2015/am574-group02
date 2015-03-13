MODULE gridData
  ! ========================================================================
  ! This module will contain information about the current physical domain
  ! and is used to simplify passing of this information throughout pAMR subroutines and functions
  ! ========================================================================

  IMPLICIT NONE
  INTEGER :: nex,nQuad,nxOut
  DOUBLE PRECISION, DIMENSION(1:2) :: xDomain
  SAVE

CONTAINS

END MODULE gridData
