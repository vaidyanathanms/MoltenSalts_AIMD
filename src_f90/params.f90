! Parameters for CP2K trajectories

MODULE PARAMS_CP2K

  USE OMP_LIB
  IMPLICIT NONE

  ! Atom, processor and time details
  INTEGER :: ntotatoms, ntotatomtypes
  INTEGER :: nframes, freqfr, nfrcntr
  REAL    :: start_time, end_time
  INTEGER :: nproc
  INTEGER :: ioncnt,  iontype
  INTEGER :: c_ioncnt, c_iontype

 
  ! Structural quantities
  INTEGER :: rdffreq,rmaxbin,npairs,rdfpaircnt
  INTEGER :: maxneighsize, neighfreq
  INTEGER :: ntotion_centers
  REAL    :: rneigh_cut,rcatan_cut,rclust_cut
  REAL    :: rvolavg,rdomcut,rbinval,rvolval

  ! All flags
  INTEGER :: box_from_file_flag, box_type_flag
  INTEGER :: rdfcalc_flag
  INTEGER :: catan_neighcalc_flag
  INTEGER :: clust_calc_flag
  INTEGER :: ion_dynflag, cion_dynflag
  INTEGER :: catan_autocfflag, catpol_autocfflag
  INTEGER :: name_to_type_map_flag
  INTEGER :: ion_diff, cion_diff
  
  ! File names and unit numbers
  CHARACTER(LEN = 256) :: ana_fname, dum_fname
  CHARACTER(LEN = 256) :: data_fname,traj_fname,log_fname,box_fname
  CHARACTER(LEN = 256) :: rdf_fname
  INTEGER, PARAMETER :: anaread = 2, logout = 3
  INTEGER, PARAMETER :: trajread = 15, boxread = 20, inpread = 100
  INTEGER, PARAMETER :: dumwrite = 50
  INTEGER, PARAMETER :: max_char = 5
  
  !Math constants
  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Trajectory file read details
  REAL :: box_xl,box_yl,box_zl, boxval
  REAL :: box_xf,box_yf,box_zf
  INTEGER*8 :: nbox_steps, timestep
  REAL*8 :: act_time

  !Required global arrays
  CHARACTER(max_char),ALLOCATABLE :: name_arr(:) 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: box_arr
  REAL,ALLOCATABLE,DIMENSION(:,:) :: masses
  REAL*8,ALLOCATABLE,DIMENSION(:,:) :: rxyz_lmp
  INTEGER,ALLOCATABLE,DIMENSION(:) :: type_arr
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: aidvals

  !Required arrays for analysis
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ionarray,counterarray
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: allionids
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_rdf
  REAL*8,ALLOCATABLE,DIMENSION(:) :: clust_avg
  REAL*8,ALLOCATABLE,DIMENSION(:,:):: rdfarray
  REAL*8,ALLOCATABLE,DIMENSION(:) :: cat_an_neighavg,an_cat_neighavg

END MODULE PARAMS_CP2K
