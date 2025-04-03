! This file contains the description of the system in the CP2K input format

&SUBSYS

  ! Read from main file or EXT_RESTART section
  @IF ( ${RESTART} == FALSE )
  
    &CELL
      ABC [angstrom] 22.0 22.0 22.0 ! => Simulation box parameters
      PERIODIC XYZ ! => Periodicity of the system
    &END CELL
  
    &TOPOLOGY
      ! From XYZ file coordinates are always Cartesian and in angstrom
      COORDINATE XYZ
      ! Name of the file that contains the atomic coordinates
      COORD_FILE_NAME system.xyz
    &END TOPOLOGY
  
  @ENDIF
  
  &KIND Al ! => For each species a KIND section is created
    BASIS_SET DZVP-MOLOPT-SR-GTH ! => Very tricky section, thorough tests
                                ! to ensure that proper basis set is used
    POTENTIAL GTH-PBE ! => CP2K pseudo-potential: GTH, ALL or KG
                          ! the exact sequence ’GTH-BLYP-q4’ appears
                          ! in the file that contains the pseudo-pot.
                          ! in front of the name of the element
  &END KIND

  &KIND Cl
    BASIS_SET DZVP-MOLOPT-GTH
    POTENTIAL GTH-PBE
  &END KIND

  &KIND C
    BASIS_SET DZVP-MOLOPT-GTH
    POTENTIAL GTH-PBE
  &END KIND
  
  &KIND H
    BASIS_SET DZVP-MOLOPT-GTH
    POTENTIAL GTH-PBE
  &END KIND

  &KIND O
    BASIS_SET DZVP-MOLOPT-GTH
    POTENTIAL GTH-PBE
  &END KIND

  &KIND N
    BASIS_SET DZVP-MOLOPT-GTH
    POTENTIAL GTH-PBE
  &END KIND

&END SUBSYS
