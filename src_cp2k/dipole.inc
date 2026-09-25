! Wannier centers first
&LOCALIZE T
  METHOD CRAZY
  EPS_LOCALIZATION 1.0E-6
  MAX_ITER 800
  USE_HISTORY TRUE
  &PRINT
    &WANNIER_CENTERS ON
      FORMAT XYZ
      IONS+CENTERS TRUE
      FILENAME =wannier_trajectory.xyz
      &EACH
        MD 1   # Calculate the Wannier centers every MD step
      &END EACH
     &END WANNIER_CENTERS
  &END PRINT
&END LOCALIZE

! Berry-phase total periodic dipole calculated independently
&PRINT
  &MOMENTS ON
    PERIODIC TRUE
    MAX_MOMENT 1
    FILENAME =total_dipole.dat
    &EACH
      MD 1
    &END EACH
  &END MOMENTS
&END PRINT

