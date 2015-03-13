&inputs
    ! Equation parameters
    meqn = 2            ! Number of equations:
                        ! 1 = advection
                        ! 2 = acoustics

    u0 = 0D0            ! For advection: wave speed
                        ! For acoustics: background flow speed

    K0 = 1D0            ! Acoustics specific
    rho0 = 1D0          ! Acoustics specific

    ! Spatial parameters
    xLower = -5D0        ! Left edge of domain
    xUpper = 5D0        ! Right edge of domain
    nex0 = 100,          ! Which resolution is run first
    nRuns = 1,          ! How many runs are done
    nScale = 2,         ! Ratio between number of elements in successive runs
    maxPolyDegree = 10, ! Maximum local degree of reconstructing polynomial

    ! Time stepping paramteters
    cflCoeff = 0.6D0    ! Ratio of used CFL number to maximum stable CFL
    tfinal = 7.0D0        ! Final time of integration

    ! Outputting parameters
    noutput = 10         ! Number of times to output, including final time (must be >= 1) (automatically includes ICs)
    outdir = '_output/variableCoeff/' ! Output directory

    ! Testing parameters
    testID = 1          ! Application specific (see qinit.f90)


    ! Misc parameters
    DEBUG = .FALSE.

/
