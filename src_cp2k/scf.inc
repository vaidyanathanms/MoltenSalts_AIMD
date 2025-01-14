! This file contains the CP2K parameters for the SCF calculation

&SCF
  MAX_SCF ${SCF_NCYCLES} ! => Maximum number of SCF cycles
  EPS_SCF ${SCF_CONV} ! => Threshold for the SCF convergence
  SCF_GUESS ${SCF_GUESS} ! => Initial guess for the wave-function
                         ! Default value of initial configuration,
                         ! Initial calculation: ATOMIC
                         ! Restart calculation: RESTART

  ! Setting up the Orbital Transformation to speed up the calculation
  ! WARNING - This requires a gap HOMO-LUMO
  ! WARNING - therefore it works only for insulators
  
  &OT ON
    MINIMIZER ${SCF_MINI} ! => Which algorithm for the minimizer DIIS or CG
    PRECONDITIONER FULL_ALL ! => Preconditioner for the minimization scheme,
                            ! FULL_ALL is the most effective state selective
                            ! preconditioner and is based on diagonalization
    ENERGY_GAP 0.001 ! => Underestimate value of the Gap (HOMO-LUMO) in a.u.
                     ! to be used with the FULL_ALL preconditioner
  &END OT

  ! If after the ${SCF_NCYCLES} first SCF steps no convergence has been reached
  ! more SCF cycles can be done updating the preconditioner.
  ! detailed information is then specified in the ’OUTER_SCF’ section
  &OUTER_SCF
    MAX_SCF ${SCF_OCYCLES} ! => We update the preconditioner and start a new cycle
                           ! up to ${SCF_OCYCLES} times, thus in this case:
                           ! 500 x 100 cycles can be computed
    EPS_SCF ${SCF_CONV} ! => Convergence threshold for the extra cycles
                        ! it is safer to keep the same value
  &END OUTER_SCF


  ! To tune the printing of the restart file

  &PRINT

    &RESTART
      LOG_PRINT_KEY T ! => Printing on screen when restart file is written
      
      &EACH
        QS_SCF 0 ! => Never write restart file(s) during the SCF cycle
        MD ${OUT_STEPS} ! => Write restart file every ${OUT_STEPS} MD steps
      &END EACH
      
      ADD_LAST NUMERIC
    &END RESTART
  
  &END PRINT

&END SCF
