! Ref: https://www.ipcms.fr/wp-content/uploads/2023/09/cp2k.pdf
! This file contains the description the CP2K method to evaluate forces and energy

&FORCE_EVAL

  ! Print forces
  &PRINT
    &FORCES
    &END
  &END

  ! For NPT simulations include stress tensor
  @IF ( ${MD_ENS} == NPT_I )
    STRESS_TENSOR ANALYTICAL
  @ENDIF

  ! For first principles calculation the method is always Quickstep
  ! which is the name of the method used to expand the wave-function in CP2K

  METHOD Quickstep
  &DFT

    ! The files that will be used in the calculation
    BASIS_SET_FILE_NAME ${BASISFILE}
    POTENTIAL_FILE_NAME ${PSEUDOFILE}
    
    ! If we restart a calculation then we use the previous wave-function
    @IF ( ${RESTART} == TRUE )
      WFN_RESTART_FILE_NAME ${WAVEFILE}
    @ENDIF

    ! Setting up of the multigrids calculation
    &MGRID
      CUTOFF ${CUTOFF} ! => Cutoff of the finest grid level
      NGRIDS ${GRIDS} ! => Number of multigrids to use, default = 4
    &END MGRID

    ! Setting up of the Quickstep algorithm
    &QS
      METHOD GPW ! => This is the default value
      EPS_DEFAULT 1.0E-12 ! => Default value is 1.0E-10
      EXTRAPOLATION ASPC ! => Recommend for MD, PS otherwise
      EXTRAPOLATION_ORDER 4 ! => 3 or 4. 4 can be better but increases CPU time
    &END QS
    
    ! Setting up the Self Consistent Field calculation parameters
    ! the SCF section is inserted using an include file:
    @INCLUDE scf.inc

    ! Setting up the exchange-correlation calculation
    &XC
      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION ORIG !Standard PBE parameterization
        &END
      &END XC_FUNCTIONAL

      &VDW_POTENTIAL

        DISPERSION_FUNCTIONAL PAIR_POTENTIAL  ! Enables vdW correction
        
        &PAIR_POTENTIAL
          TYPE DFTD3  ! Uses Grimme’s D3 correction
          REFERENCE_FUNCTIONAL PBE  ! Must match the XC functional used
          PARAMETER_FILE_NAME ./dftd3.dat
          CALCULATE_C9_TERM TRUE  ! Enables three-body dispersion term
          R_CUTOFF 15.0  ! Cutoff distance for dispersion interactions (Å)
        &END PAIR_POTENTIAL

      &END VDW_POTENTIAL

    &END XC

  &END DFT


  ! The FORCE_EVAL also requires a SUBSYS section to describe the system
  ! the SUBSYS section is inserted using an include file:i
  @INCLUDE subsys.inc

&END FORCE_EVAL
