SUBROUTINE outputData(qvals,xvals,elCent,localPolyDegree,time,iout,p)
! ==============================================================================
! outputs t,x, and q fields to output files
! INPUTS: qvals(1:nx,1:meq) - q field value at xvals locations
!         xvals(1:nx) - x locations for plotting
!         time - current output time
!         iout - current output frame
! OUTPUTS: -none-
! ==============================================================================
  USE testParameters, ONLY: meqn,testID,maxPolyDegree,outdir
  USE gridData
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:meqn), INTENT(INOUT) :: qvals
  DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xvals
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elCent
  INTEGER, DIMENSION(1:nex), INTENT(IN) :: localPolyDegree
  DOUBLE PRECISION, INTENT(IN) :: time
  INTEGER, INTENT(IN) :: iout,p
  ! Outputs
  ! Local variables
  INTEGER :: i,j,m
  CHARACTER(len=40) :: fileNameq,fileNameAux
  CHARACTER(len=4) :: ioutName,runName

  WRITE(ioutName, '(i4.4)') iout
  WRITE(runName, '(a3,i1)') 'run',p
  fileNameq = TRIM(outdir)//TRIM(runName)//'out.q'//TRIM(ioutName)
  fileNameAux = TRIM(outdir)//TRIM(runName)//'out.aux'//TRIM(ioutName)
  ! Open files
  OPEN(unit=50,file=TRIM(fileNameq),status='unknown',form='formatted')
  OPEN(unit=60,file=TRIM(fileNameAux),status='unknown',form='formatted')

  ! Print out header information for out.qxxxx
  WRITE(50,1001) maxPolyDegree,nex,xDomain(1),xDomain(2),time
1001 FORMAT(&
       i5,'                        maxPolyDegree',/,&
       i5,'                        nex',/, &
       f6.3,'                        xLeft',/,&
       f6.3,'                        xRight',/,&
       e26.16,'   time')

  WRITE(50,*) ''
  WRITE(50,'(2a26)') 'x','q'
  WRITE(50,*) ''

  DO i=1,nxOut
     DO m=1,meqn
      ! Exponents with more than 2 digits cause problems reading
      ! into matlab... reset tiny values to zero:
        IF (dabs(qvals(i,m)) .lt. 1d-99) qvals(i,m) = 0.d0
     ENDDO

     WRITE(50,1005) xvals(i),(qvals(i,m), m=1,meqn)
  ENDDO
1005    FORMAT(50e26.16)



  ! Write out auxilliary array information
  ! Headers
  WRITE(60,*) ''
  WRITE(60,'(2a26)') 'elCent','local poly degree'
  WRITE(60,*) ''

  DO i=1,nex
    WRITE(60,1006) elCent(i),localPolyDegree(i)
  ENDDO
  1006   FORMAT(1e26.16, i26.16)

  ! Close files
  CLOSE(50)
  CLOSE(60)

END SUBROUTINE outputData
