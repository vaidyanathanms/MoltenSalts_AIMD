! Program to analyze CP2K trajectories
! Author: Vaidyanathan M. Sethuraman
! Version: 020925

PROGRAM ANALYZE_CP2K

  USE PARAMS_CP2K
  IMPLICIT NONE

! Print headers
  PRINT *, "Static analysis of SIC system .."
  PRINT *, "Starting OMP Threads .."

!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  PRINT *, "Number of threads: ", nproc

! Call functions
  CALL READ_ANA_INP_FILE()
  CALL READ_DATAFILE()
  CALL SORTALLARRAYS()
  CALL ALLOCATE_ANALYSIS_ARRAYS()
  CALL ANALYZE_TRAJECTORYFILE()
  CALL ALLOUTPUTS()
  CALL DEALLOCATE_ARRAYS()

END PROGRAM ANALYZE_CP2K

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_INP_FILE()

  USE PARAMS_CP2K
  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j
  CHARACTER(256) :: dumchar
  CHARACTER(max_char) :: aname
  
  CALL DEFAULTVALUES()

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     ! Read file and trajectory details
     IF(dumchar == 'datafile') THEN
        
        READ(anaread,*,iostat=ierr) data_fname

     ELSEIF(dumchar == 'trajectory_file') THEN

        READ(anaread,*,iostat=ierr) traj_fname

     ELSEIF(dumchar == 'nframes') THEN

        READ(anaread,*,iostat=ierr) nframes

     ELSEIF(dumchar == 'start_time') THEN

        READ(anaread,*,iostat=ierr) start_time

     ELSEIF(dumchar == 'end_time') THEN

        READ(anaread,*,iostat=ierr) end_time

     ELSEIF(dumchar == 'freqfr') THEN

        READ(anaread,*,iostat=ierr) freqfr

     ELSEIF(dumchar == 'name_to_type_map') THEN

        READ(anaread,*,iostat=ierr) ntotatomtypes
        
        ALLOCATE(name_arr(ntotatomtypes),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate name_arr"
        ALLOCATE(type_arr(ntotatomtypes),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate type_arr"

        DO i = 1, ntotatomtypes

           READ(anaread,*,iostat=ierr) name_arr(i), type_arr(i)

        END DO
        
        name_to_type_map_flag = 1
        
     ELSEIF(dumchar == 'box_dim') THEN

        box_type_flag = 1
        READ(anaread,*,iostat=ierr) dumchar
        
        IF(trim(dumchar) == 'from_file') THEN

           READ(anaread,*,iostat=ierr) box_fname
           box_from_file_flag = 1
           
        ELSEIF(trim(dumchar) == 'fixed') THEN

           READ(anaread,*,iostat=ierr) box_xf, box_yf, box_zf

        END IF
           
     ELSEIF(dumchar == 'ion_type') THEN
        
        READ(anaread,*,iostat=ierr) iontype
        
     ELSEIF(dumchar == 'cion_type') THEN

        READ(anaread,*,iostat=ierr) c_iontype

     !Here onwards static properties
     ELSEIF(dumchar == 'compute_rdf') THEN

        rdfcalc_flag = 1
        READ(anaread,*,iostat=ierr) rdffreq, rmaxbin, rdomcut,npairs
        
        ALLOCATE(pairs_rdf(npairs,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"
      
        DO i = 1,npairs

           READ(anaread,*,iostat=ierr) pairs_rdf(i,1), pairs_rdf(i,2)

        END DO

     ELSEIF(dumchar == 'compute_clust') THEN

        clust_calc_flag = 1
        READ(anaread,*,iostat=ierr) rclust_cut
        
     ELSEIF(dumchar == 'compute_catanneigh') THEN

        catan_neighcalc_flag = 1
        READ(anaread,*,iostat=ierr) neighfreq,maxneighsize,rneigh_cut
        
     !Here onwards dynamic properties
!!$     ELSEIF(dumchar == 'compute_iondiff') THEN
!!$
!!$        READ(anaread,*,iostat=ierr) ion_diff
!!$        ion_dynflag = 1
!!$
!!$     ELSEIF(dumchar == 'compute_ciondiff') THEN
!!$
!!$        READ(anaread,*,iostat=ierr) cion_diff
!!$        cion_dynflag = 1
!!$
!!$     ELSEIF(dumchar == 'compute_catanrestime') THEN
!!$        
!!$        READ(anaread,*,iostat=ierr) rcatan_cut
!!$        catan_autocfflag = 1
!!$        ion_dynflag = 1; cion_dynflag = 1
!!$
     ! Read log filename
     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) log_fname = "log_"//trim(adjustl(traj_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  IF(box_type_flag == 0) THEN
     PRINT *, "ERROR: Box type not set: set fixed or from_file"
     STOP
  END IF

  PRINT *, "Analysis input file read finished .."

  CALL SANITY_CHECK_IONTYPES()
  
END SUBROUTINE READ_ANA_INP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE PARAMS_CP2K

  IMPLICIT NONE

  ! Frame, molecules and processor details
  nframes = 0; freqfr = 0; nfrcntr = 0
  start_time = 0.0; end_time = 0.0
  
  ! Initialize flags
  box_from_file_flag = 0
  box_type_flag = 0
  rdfcalc_flag = 0
  ion_dynflag = 0; cion_dynflag = 0
  catan_autocfflag = 0
  name_to_type_map_flag = 0
  ion_diff = 0; cion_diff = 0

  ! Initialize iontypes
  c_iontype = -1; iontype = -1

  !Initialize system quantities
  ioncnt = 0; c_ioncnt = 0

  ! Initialize distributions and frequencies
  rdffreq = 0

  ! Initialzie structural quantities
  rdomcut = 10.0;  rmaxbin = 100; rbinval = REAL(rdomcut)&
       &/REAL(rmaxbin)
  rcatan_cut = 0.0; rneigh_cut = 0.0; rclust_cut = 0.0

  ! Initialize structural averages
  rvolavg = 0

END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,j,ierr,k,u,AllocateStatus,imax
  INTEGER :: flag, cntr, nwords
  INTEGER :: molid,atype,ix,iy,iz,mtype
  REAL    :: charge,rx,ry,rz,massval
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(len=256) :: rline,dumchar
  CHARACTER(len=max_char) :: aname

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used: ", trim(adjustl(data_fname))

  ntotatoms = 0
  
  READ(inpread,*) ntotatoms
  READ(inpread,*) 

  PRINT *, "STATISTICS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_TOPO_ARRAYS()
         
  masses = -1
  
  DO j = 1,ntotatoms
     
     READ(inpread,*) aname,rx,ry,rz
     
     CALL MAP_ANAME_TO_ATYPE(aname,atype,k,0)
     
     aidvals(j,1)     = j
     aidvals(j,2)     = 1
     aidvals(j,3)     = atype
     rxyz_lmp(j,1)    = rx
     rxyz_lmp(j,2)    = ry
     rxyz_lmp(j,3)    = rz
     
     IF(masses(k,1) == -1) THEN
        
        CALL ASSIGN_MASSES(aname,mtype,j,massval)
        masses(k,1) = INT(mtype)
        masses(k,2) = massval
        
     END IF
        
  END DO

  PRINT *, "Writing mass and type data..."
  DO i = 1,ntotatomtypes
     WRITE(logout,*) name_arr(i),type_arr(i),masses(i,1),masses(i,2)
  END DO
     
  PRINT *, "Datafile read completed..."

  CLOSE(inpread)

END SUBROUTINE READ_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE MAP_ANAME_TO_ATYPE(aname,atype,kcnt,tval)

  USE PARAMS_CP2K
  IMPLICIT NONE
  
  CHARACTER(len=*), INTENT(IN) :: aname
  REAL, INTENT(IN) :: tval
  INTEGER, INTENT(OUT) :: atype,kcnt
  INTEGER :: atom_find_flag
  
  atom_find_flag = -1

  DO kcnt = 1,ntotatomtypes

     IF(trim(adjustl(aname)) == trim(adjustl(name_arr(kcnt)))) THEN
        
        atype = type_arr(kcnt)
        atom_find_flag = 1
        EXIT
        
     END IF
     
  END DO

  IF(atom_find_flag == -1) THEN

     
     PRINT *, "ERROR: Unknown atom name ", aname, " at ", tval
     STOP

  END IF
    

END SUBROUTINE MAP_ANAME_TO_ATYPE

!--------------------------------------------------------------------

SUBROUTINE ASSIGN_MASSES(aname,atype,aid,massval)

  USE PARAMS_CP2K
  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN) :: aname
  INTEGER, INTENT(IN) :: atype, aid
  REAL, INTENT(OUT) :: massval

  massval = 1.0 ! Assign unity as default
  IF(trim(adjustl(aname)) == 'Al') THEN
     
      massval = 26.981539

   ELSEIF(trim(adjustl(aname)) == 'Cl') THEN

      massval = 35.453

   ELSEIF(trim(adjustl(aname)) == 'H') THEN

      massval = 1.00784

   ELSEIF(trim(adjustl(aname)) == 'C') THEN

      massval = 12.011

   ELSEIF(trim(adjustl(aname)) == 'O') THEN

      massval = 15.999

   ELSEIF(trim(adjustl(aname)) == 'N') THEN

      massval = 14.0067

   ELSE

      PRINT *, trim(adjustl(aname)), " and ID ", aid, " not found!"
      PRINT *, "WARNING: Assigning unit mass to ", atype

   END IF
  
 END SUBROUTINE ASSIGN_MASSES

!--------------------------------------------------------------------

SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: ierr,atchk,atype,jumpfr
  INTEGER :: skip_cnt, at_cnt, mtype
  CHARACTER(len=256) :: line_read
  REAL :: eval
  CHARACTER(len=5) :: aname
  
  OPEN(unit = trajread,file =trim(traj_fname),action="read",status="ol&
       &d",iostat=ierr)
  IF(ierr /= 0) STOP "trajectory file not found"
  PRINT *, "Trajectory file used: ",trim(adjustl(traj_fname))
  WRITE(logout,*) "Trajectory file used: "&
       &,trim(adjustl(traj_fname))

  PRINT *, "Analyzing trajectory file..."

  CALL READ_AND_ASSIGN_BOX_DIM()
  CALL STRUCT_INIT()
  
  DO 

     READ(trajread,*) atchk
     READ(trajread,'(A)',iostat=ierr) line_read
     CALL PROCESS_HEADER_LINES(trim(adjustl(line_read)),timestep&
          &,act_time,eval,atchk,ierr)

     IF(mod(act_time,100.0) == 0.0) PRINT *,"Time: ", act_time
     IF(act_time .LT. start_time) THEN
        DO at_cnt = 1,atchk
           READ(trajread,*) 
        END DO
        CYCLE
     END IF
     
     IF(act_time .GT. end_time) EXIT

     IF(nfrcntr == 0) THEN
        PRINT *, "Starting time: ", act_time
        WRITE(logout,*) "Starting time: ", act_time
     END IF
     
     ! Find box-dimension for NPT systems
     CALL FIND_BOX_DIM(act_time)
     
     ! Read trajectory
     DO at_cnt = 1,atchk
        
        READ(trajread,*) aname,rxyz_lmp(at_cnt,1),rxyz_lmp(at_cnt&
             &,2),rxyz_lmp(at_cnt,3)
        
        CALL MAP_ANAME_TO_ATYPE(aname,atype,mtype,act_time)
        
        IF(atype .NE. aidvals(at_cnt,3)) THEN
           
           PRINT *, "Incorrect atom ids"
           PRINT *, timestep, act_time
           PRINT *, at_cnt,aname,atype,aidvals(at_cnt,3)
           STOP
           
        END IF
        
     END DO

     nfrcntr = nfrcntr + 1
     IF(nfrcntr == 1) PRINT *, "Beginning statics analysis..."
     CALL STRUCT_MAIN(nfrcntr)

  END DO

  PRINT *, "Trajectory read completed .."
  PRINT *, "Last frame analyzed ..", act_time
  PRINT *, "Total frames analyze ..", nfrcntr
  WRITE(logout,*) "Total frames analyze ..", nfrcntr
  WRITE(logout,*) "Last frame analyzed ..", act_time
  
!!$     DO jumpfr = 1,freqfr
!!$
!!$        READ(trajread,*)
!!$        READ(trajread,*)        
!!$        READ(trajread,*)
!!$ 
!!$        READ(trajread,*) atchk
!!$
!!$        DO at_cnt = 1,atchk+5
!!$
!!$           READ(trajread,*) 
!!$
!!$        END DO
!!$
!!$     END DO
!!$
!!$  END DO

  CLOSE(trajread)

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE PROCESS_HEADER_LINES(rline,ival,tval,enval,atomstval,ierr)

  USE PARAMS_CP2K
  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN) :: rline
  REAL, INTENT(OUT) ::  tval, enval
  INTEGER, INTENT(OUT) :: ival
  INTEGER,INTENT(IN) :: atomstval,ierr
  INTEGER :: ierr2
  CHARACTER(len=8) :: char1, char2, char3
  
  ! Check line
  IF (ierr /= 0) THEN
     PRINT *, "ERROR: reading ", rline
     STOP
  END IF

  READ(rline,'(3X,I10,8X,F20.10,5X,F20.10)',iostat=ierr2) ival, tval,&
       & enval

  ! Check for errors in header lines
  IF (ierr2 /= 0) THEN
     PRINT *, "ERROR: assigning ", rline
     STOP
  END IF
  IF(atomstval > ntotatoms) THEN
     PRINT *, "ERROR: More atoms than in the datafile: "
     PRINT *, ntotatoms, atomstval, ival, tval
     STOP
  END IF

END SUBROUTINE PROCESS_HEADER_LINES
  
!--------------------------------------------------------------------

SUBROUTINE READ_AND_ASSIGN_BOX_DIM()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: ierr, AllocateStatus, icnt, temp
  INTEGER :: dum_time
  REAL :: box_xx, box_xy, box_xz
  REAL :: box_yx, box_yy, box_yz
  REAL :: box_zx, box_zy, box_zz
  REAL :: dumvol
  CHARACTER (len=256) :: dumread

  IF (box_from_file_flag == 1) THEN

     PRINT *, "Reading box dimensions from file.."
     OPEN(unit = boxread,file =trim(box_fname),action="read",status&
          &="old",iostat=ierr)
     IF(ierr /= 0) STOP "Box file not found & `from_file` flag is ON"
     PRINT *, "Box file used: ",trim(adjustl(box_fname))
     WRITE(logout,*) "Box file used: ",trim(adjustl(box_fname))

     ! Find number of lines in the file
     READ(boxread,*)
     nbox_steps = 0
     DO 
        READ(boxread,*,iostat=ierr) dumread
        IF (ierr /= 0) EXIT
        nbox_steps = nbox_steps + 1
     END DO
     
     CALL ALLOCATE_BOX_ARRAYS()
     CLOSE(boxread)
  
     ! Re-read box file
     OPEN(unit = boxread,file =trim(box_fname),action="read",status&
          &="old",iostat=ierr)
     IF(ierr /= 0) STOP "Re-reading box-file failed!.."

     READ(boxread,*)
     DO icnt = 1, nbox_steps
        READ(boxread,*) temp, dum_time, box_xx, box_xy, box_xz,&
             & box_yx,box_yy, box_yz, box_zx, box_zy, box_zz, dumvol
        
        box_arr(icnt,1) = dum_time
        box_arr(icnt,2) = box_xx
        box_arr(icnt,3) = box_yy
        box_arr(icnt,4) = box_zz
        
     END DO

  ELSE

     CALL ALLOCATE_BOX_ARRAYS()
     PRINT *, "Fixed box dimensions..", box_xf, box_yf, box_zf
     
  END IF

  PRINT *, "Box dimensions assigned successfully.."
  
END SUBROUTINE READ_AND_ASSIGN_BOX_DIM

!--------------------------------------------------------------------

SUBROUTINE FIND_BOX_DIM(tval)

  USE PARAMS_CP2K
  IMPLICIT NONE
  
  REAL, INTENT(IN)  :: tval
  INTEGER :: icnt, boxflag

  boxflag = -1
  
  IF (box_from_file_flag == 1) THEN

     DO icnt = 1, nbox_steps

        IF(tval == box_arr(icnt,1)) THEN

           box_xl = box_arr(icnt,2)
           box_yl = box_arr(icnt,3)
           box_zl = box_arr(icnt,4)
           boxflag = 1
           
           EXIT
           
        END IF

     END DO

  ELSE

     box_xl = box_xf; box_yl = box_yf; box_zl = box_zf
     boxflag = 1

  END IF

  IF(boxflag == -1) THEN
     
     PRINT *, "ERROR: Box flag to read from file is ON!"
     PRINT *, "ERROR: Unknown time-step for analysis: ", tval
     STOP

  END IF
     
END SUBROUTINE FIND_BOX_DIM

!--------------------------------------------------------------------

SUBROUTINE STRUCT_INIT()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2,norm,acnt,fcnt,a1id,molid,flagch,flagpr,jmax
  INTEGER :: AllocateStatus

  IF(rdfcalc_flag) THEN

     rdfarray = 0.0
     rbinval = rdomcut/REAL(rmaxbin)

     DO i = 1, npairs

        t1 = 0; t2 = 0

        DO j = 1,ntotatoms

           IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
           IF(aidvals(j,3) == pairs_rdf(i,2)) t2 = t2+1

        END DO

        IF(pairs_rdf(i,1) == pairs_rdf(i,2)) THEN
           pairs_rdf(i,3) = t1*(t1-1) !g_AA(r)
        ELSE
           pairs_rdf(i,3) = t1*t2 !g_AB(r)
        END IF

     END DO

  END IF



END SUBROUTINE STRUCT_INIT

!--------------------------------------------------------------------

SUBROUTINE SANITY_CHECK_IONTYPES()

  USE PARAMS_CP2K
  IMPLICIT NONE

  IF(ion_dynflag .OR. catan_neighcalc_flag) THEN

     IF(iontype == -1) THEN

        PRINT *, "ion type undefined for neigh or diff calculation"
        STOP

     END IF

  END IF

  IF(cion_dynflag .OR. catan_neighcalc_flag) THEN

     IF(c_iontype == -1) THEN

        PRINT *, "counter-ion type undefined for neigh or diff calculation"
        STOP

     END IF

  END IF

END SUBROUTINE SANITY_CHECK_IONTYPES

!--------------------------------------------------------------------

SUBROUTINE STRUCT_MAIN(tval)

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval
  INTEGER :: t1, t2
  INTEGER :: clock_rate, clock_max

  IF(rdfcalc_flag) THEN

     IF(tval == 1) THEN

        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL COMPUTE_RDF(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for RDF analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     END IF
     IF(mod(tval-1,rdffreq)==0) CALL COMPUTE_RDF(tval)

  END IF

  IF(catan_neighcalc_flag) THEN

     IF(tval == 1) THEN

        cat_an_neighavg = 0.0; an_cat_neighavg=0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CAT_AN_NEIGHS()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for neighbor analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     END IF

     IF(mod(tval,neighfreq) == 0) CALL CAT_AN_NEIGHS()

  END IF

  IF(clust_calc_flag) THEN

     IF(tval == 1) THEN

        clust_avg = 0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CLUSTER_ANALYSIS(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for cluster analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'
     ELSE

        CALL CLUSTER_ANALYSIS(tval)

     END IF

  END IF


END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE SORTALLARRAYS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,j,a1type,cnt,AllocateStatus,ntotion_cnt,aid,molid
  INTEGER, DIMENSION(1:ntotatoms,2) :: dumsortarr,dumcionarr

  dumsortarr = -1; dumcionarr = -1
  cnt = 0
  ntotion_cnt = 0

  DO i = 1,ntotatoms

     a1type = aidvals(i,3)

     IF(a1type == iontype) THEN
        ioncnt = ioncnt + 1
        dumsortarr(ioncnt,1) = i
        dumsortarr(ioncnt,2) = a1type

     ELSEIF(a1type == c_iontype) THEN
        c_ioncnt = c_ioncnt + 1
        dumcionarr(c_ioncnt,1) = i
        dumcionarr(c_ioncnt,2) = a1type

     END IF

  END DO

  ! Always identify ion and counter-ion types
  PRINT *, "Number of atoms of ion type: ", ioncnt
  PRINT *, "Number of atoms of cntion type: ", c_ioncnt

  ALLOCATE(ionarray(ioncnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate ionarray"
  ALLOCATE(counterarray(c_ioncnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate counterarray"

  ! Load ion array

  i = 0

  DO WHILE(dumsortarr(i+1,1) .NE. -1)

     i = i + 1
     ionarray(i,1) = dumsortarr(i,1)
     ionarray(i,2) = dumsortarr(i,2)

  END DO

    IF(i .NE. ioncnt) THEN
     PRINT *, i, ioncnt
     STOP "Wrong total count in ionarray"
  END IF

  DO i = 1,ioncnt

     IF(ionarray(i,1) == -1 .OR. ionarray(i,2) == -1) THEN

        PRINT *, i,ionarray(i,1), ionarray(i,2)
        PRINT *, "Something wrong in assigning ionarray"
        STOP

     END IF

     IF(ionarray(i,2) .NE. iontype) THEN

        PRINT *, i,ionarray(i,1), ionarray(i,2)
        PRINT *, "Something wrong in ionarray type"
        STOP

     END IF

  END DO

    OPEN(unit = 93,file="iontypelist.txt",action="write",status="replace&
       &")

  WRITE(93,*) "Reference type/count: ", iontype, ioncnt

  DO i = 1,ioncnt
     WRITE(93,'(3(I0,1X))') i, ionarray(i,1), ionarray(i,2)
  END DO

  CLOSE(93)

  ! Load counterion array

  i = 0

  DO WHILE(dumcionarr(i+1,1) .NE. -1)

     i = i + 1
     counterarray(i,1) = dumcionarr(i,1)
     counterarray(i,2) = dumcionarr(i,2)

  END DO

  IF(i .NE. c_ioncnt) THEN
     PRINT *, i, c_ioncnt
     STOP "Wrong total count in counterarray"
  END IF

  DO i = 1,c_ioncnt

     IF(counterarray(i,1) == -1 .OR. counterarray(i,2) == -1) THEN

        PRINT *, i,counterarray(i,1), counterarray(i,2)
        PRINT *, "Something wrong in assigning counterarray"
        STOP

     END IF

     IF(counterarray(i,2) .NE. c_iontype) THEN

        PRINT *, i,counterarray(i,1), counterarray(i,2)
        PRINT *, "Something wrong in counterionarray type"
        STOP

     END IF

  END DO


  OPEN(unit = 93,file="cntionlist.txt",action="write",status="repl&
       &ace")

  WRITE(93,*) "Reference type/count: ", c_iontype, c_ioncnt

  DO i = 1,c_ioncnt

     WRITE(93,'(3(I0,1X))') i, counterarray(i,1), counterarray(i,2)

  END DO

  CLOSE(93)

   ! Cluster calc requires to add iontype and c_iontype in same array
  IF (clust_calc_flag) THEN

     ntotion_centers = ioncnt + c_ioncnt
     PRINT *, "Total number of ion centers: ", ntotion_centers
     cnt = 1

     ALLOCATE(allionids(ntotion_centers,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate allionids"
     ALLOCATE(clust_avg(ntotion_centers),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate clust_avg"

     allionids = 0 ! Initial allocation

     DO i = 1,ntotatoms

        a1type = aidvals(i,3)

        IF(a1type == iontype .OR. a1type == c_iontype) THEN

           allionids(cnt,1) = i
           allionids(cnt,2) = a1type
           cnt = cnt + 1

        END IF

     END DO

  ELSE

     ALLOCATE(allionids(1,2),stat = AllocateStatus)
     DEALLOCATE(allionids)
     ALLOCATE(clust_avg(1),stat = AllocateStatus)
     DEALLOCATE(clust_avg)

  END IF

END SUBROUTINE SORTALLARRAYS
  
!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF(iframe)

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER :: a1ref,a2ref
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval

  ALLOCATE(dumrdfarray(0:rmaxbin-1,npairs),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"
  dumrdfarray = 0

!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,paircnt,a1ref,a2ref) REDUCTION(+:dumrdfarray)
  DO paircnt = 1,npairs

     a1ref = pairs_rdf(paircnt,1); a2ref = pairs_rdf(paircnt,2)

     DO i = 1,ntotatoms

        a1id   = aidvals(i,1)
        a1type = aidvals(i,3)

        DO j = 1,ntotatoms

           a2id   = aidvals(j,1)
           a2type = aidvals(j,3)

           ! Remove identical IDs when computing g_AA(r)
           IF(a1id == a2id .AND. a1ref == a2ref) CYCLE


           IF(a1type == a1ref .AND. a2type == a2ref) THEN

              rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1)
              ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2)
              rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3)

              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3)

              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rzval - box_yl*ANINT(rzval/box_zl)

              rval = sqrt(rxval**2 + ryval**2 + rzval**2)
              ibin = FLOOR(rval/rbinval)

              IF(ibin .LT. rmaxbin) THEN

                 dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                      &,paircnt) + 1

              END IF

           END IF

        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j)
  DO j = 1,npairs

     DO i = 0,rmaxbin-1

        rdfarray(i,j) = rdfarray(i,j) + REAL(dumrdfarray(i,j))&
             &*rvolval/(REAL(pairs_rdf(j,3)))

     END DO

  END DO
!$OMP END DO

!$OMP END PARALLEL

  DEALLOCATE(dumrdfarray)

END SUBROUTINE COMPUTE_RDF

!--------------------------------------------------------------------

SUBROUTINE CAT_AN_NEIGHS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,neigh_cnt,tid
  INTEGER,DIMENSION(1:maxneighsize,0:nproc-1) :: cat_an_neigh_inst&
       &,an_cat_neigh_inst
  REAL :: rxval, ryval, rzval, rval

  cat_an_neigh_inst = 0; an_cat_neigh_inst = 0

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,ioncnt

     neigh_cnt = 0
     a1id = ionarray(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,c_ioncnt

        a2id = counterarray(j,1)

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1)
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2)
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rneigh_cut) THEN

           neigh_cnt = neigh_cnt + 1

        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     cat_an_neigh_inst(neigh_cnt+1,tid) = cat_an_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,c_ioncnt

     neigh_cnt = 0
     a1id = counterarray(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,ioncnt

        a2id = ionarray(j,1)

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1)
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2)
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rneigh_cut) THEN

           neigh_cnt = neigh_cnt + 1

        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     an_cat_neigh_inst(neigh_cnt+1,tid) = an_cat_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO

!$OMP DO 
  DO  i = 1,maxneighsize
     DO j = 0,nproc-1
        cat_an_neighavg(i) = cat_an_neighavg(i) + cat_an_neigh_inst(i&
             &,j)
        an_cat_neighavg(i) = an_cat_neighavg(i) + an_cat_neigh_inst(i&
             &,j)
     END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE CAT_AN_NEIGHS

!--------------------------------------------------------------------

SUBROUTINE CLUSTER_ANALYSIS(frnum)

  USE PARAMS_CP2K
  IMPLICIT NONE

!Ref Sevick et.al ., J Chem Phys 88 (2)

  INTEGER :: i,j,k,a2ptr,a1id,a2id,itype,jtype,jptr,idum,jflag,jcnt&
       &,iflag,jtot,jind,jprev
  INTEGER, DIMENSION(ntotion_centers,ntotion_centers) :: all_direct&
       &,catan_direct,all_neigh,catan_neigh
  INTEGER, DIMENSION(1:ntotion_centers) :: union_all,flag_catan,scnt&
       &,all_linked
  REAL :: rxval, ryval, rzval, rval
  INTEGER, INTENT(IN) :: frnum

!$OMP PARALLEL SHARED(catan_direct,catan_neigh)

!$OMP DO PRIVATE(i,j)
  DO i = 1,ntotion_centers

     scnt(i) = 0; all_linked(i)  = 0
     union_all(i) = -1; flag_catan(i) = -1

     DO j = 1,ntotion_centers

        IF(i .NE. j) THEN
           all_direct(i,j) = 0
           catan_direct(i,j) = 0
        END IF

        IF(i == j) THEN
           all_direct(i,j) = 1
           catan_direct(i,j) = 0
        END IF

        all_neigh(i,j) = 0
!        catan_neigh(i,j) = 0

     END DO

  END DO
!$OMP END DO

!Create Direct connectivity matrix
!allneigh - does not distinguish between Li and P neigh
!catan_neigh - neighbors with sequence cat-an-cat-an.. or an-cat-an-cat...

!$OMP DO PRIVATE(i,j,a1id,a2ptr,a2id,rxval,ryval,rzval,rval,itype&
!$OMP& ,jptr,jtype)  
  DO i = 1,ntotion_centers

     a1id = allionids(i,1)
     a2ptr = 1
     itype = aidvals(a1id,3)
     jptr  = 1
     all_neigh(i,i) = a1id

     DO j = 1,ntotion_centers

        a2id = allionids(j,1)

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1)
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2)
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rclust_cut .AND. a1id .NE. a2id) THEN

           all_direct(i,j) = 1
           all_neigh(i,j)  = a2id
!           all_neigh(i,a2ptr) = a2id
!           a2ptr = a2ptr + 1

           jtype = aidvals(a2id,3)

           IF(itype .NE. jtype) THEN

              catan_direct(i,j) = 1
!              catan_neigh(i,j)  = a2id
!              catan_neigh(i,jptr+1) = a2id
              itype = jtype
!              jptr  = jptr + 1

           END IF

        END IF

     END DO

  END DO

!$OMP END DO  

  
!Check for symmetry
  IF(frnum == 1) THEN
!$OMP DO
     DO i = 1,ntotion_centers

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) .NE. all_direct(j,i)) STOP "Unsymmetric&
                & all_direct"

          IF(all_neigh(i,j) .NE. 0) THEN

              IF(all_neigh(i,j) .NE. all_neigh(j,j) .OR. all_neigh(j&
                   &,i) .NE. all_neigh(i,i)) THEN

                 PRINT *, i,j,all_direct(i,j),all_direct(j,i)&
                      &,all_neigh(j,i),all_neigh(i,i)
                 STOP "Unsymmetric neighbor list"

              END IF

           END IF

        END DO

     END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL        

  !Intersection

  DO i = 1,ntotion_centers-1 !Ref row counter

     iflag = 0
     idum  = i

     DO WHILE(iflag == 0 .AND. union_all(i) == -1)

        jflag = 0
        k    = 1 !Column counter
        j    = idum+1 !Other row counter

        DO WHILE(jflag == 0 .AND. k .LE. ntotion_centers)

           IF((all_direct(i,k) == all_direct(j,k)).AND. all_direct(i&
                &,k)== 1) THEN

              jflag = 1
!!$              jprev = 0

              DO jcnt = 1,ntotion_centers


!!$                 IF(all_direct(j,jcnt) == 1) jprev = 1

                 !Replace highest row by union of two rows

                 all_direct(j,jcnt) = all_direct(i,jcnt) .OR.&
                      & all_direct(j,jcnt)

!!$                 IF((all_direct(j,jcnt) == 1 .AND. all_direct(i,jcnt)&
!!$                      &==1) .AND. jprev == 0) THEN
!!$                    
!!$                    all_neigh(j,jcnt) = all_neigh(i,jcnt)
!!$                    jprev = 0 !Other condition is already
!!$                    ! incorporated before
!!$                 END IF
!!$                 
              END DO

              union_all(i) = 1 !One match implies the low ranked row
              ! is present in high ranked row

           ELSE

              k = k + 1

           END IF

        END DO

        IF(union_all(i) == 1) THEN

           iflag = 1

        ELSE

           idum  = idum + 1

        END IF

        IF(idum == ntotion_centers) iflag = 1

     END DO

  END DO

!Count
  jtot = 0
!$OMP PARALLEL PRIVATE(i,j,jind) 
!$OMP DO
  DO i = 1,ntotion_centers

     IF(union_all(i) == -1) THEN

        jind = 0

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) == 1) jind = jind + 1

        END DO

        scnt(jind) = scnt(jind) + 1
        all_linked(i) = jind

     END IF

  END DO
!$OMP END DO

!$OMP DO

  DO i = 1,ntotion_centers

     clust_avg(i) = clust_avg(i) + scnt(i)

  END DO
!$OMP END DO

!$OMP END PARALLEL

  IF(frnum == 1) THEN
     OPEN(unit =90,file ="scnt.txt",action="write",status="replace")
  END IF

  jtot = 0

  DO i = 1,ntotion_centers

          IF(frnum == 1) WRITE(90,*) i,scnt(i)
     jtot = jtot + all_linked(i)

  END DO

  IF(jtot .NE. ntotion_centers) THEN

     PRINT *, "Sum of ions not equal to total ion centers"
     PRINT *, jtot, ntotion_centers
     STOP

  END IF

  IF(frnum == 1) CLOSE(90)

  IF(frnum == 1) THEN

     OPEN(unit =90,file ="all_neigh.txt",action="write",status="replace")

     DO i = 1,ntotion_centers

        IF(union_all(i) == -1) THEN

           WRITE(90,*) i,all_linked(i)

           DO j = 1,ntotion_centers

              IF(all_direct(i,j) == 1) WRITE(90,*) j,allionids(j,1),&
                   & allionids(j,2),all_direct(i,j)

           END DO

        END IF

     END DO

     CLOSE(90)

  END IF

END SUBROUTINE CLUSTER_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,ierr

  IF (nframes == 0 .AND. nfrcntr .GT. 0) nframes = nfrcntr
  PRINT *, "Number of frames from start to end: ", nframes/(freqfr+1)
  PRINT *, "Frequency of Frames: ", freqfr + 1
  PRINT *, "Total number of Frames analyzed: ", nfrcntr

  WRITE(logout,*) "Number of frames from start to end: ", nframes&
       &/(freqfr+1)
  WRITE(logout,*) "Frequency of Frames: ", freqfr+1
  WRITE(logout,*) "Total number of Frames analyzed: ", nfrcntr

  IF(rdfcalc_flag) THEN
     PRINT *, "Writing RDFs .."
     CALL OUTPUT_ALLRDF()
  END IF

  IF(catan_neighcalc_flag) THEN
     PRINT *, "Writing neighbors .."
     CALL OUTPUT_ALLNEIGHBORS()
  END IF

  IF(clust_calc_flag) THEN

     PRINT *, "Writing cluster outputs .."
     CALL OUTPUT_ALLCLUSTERS()

  END IF
  
END SUBROUTINE ALLOUTPUTS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLCLUSTERS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,ierr
  
  dum_fname = "clust_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
  
  IF(ierr /= 0) PRINT *, "Unknown clust_filename"
  
  DO i = 1,ntotion_centers
     
     WRITE(dumwrite,'(I0,1X,F14.8,1X)') i, REAL(clust_avg(i))&
          &/REAL(nframes)
     
  END DO
  CLOSE(dumwrite)

END SUBROUTINE OUTPUT_ALLCLUSTERS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLNEIGHBORS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,frnorm,ierr
  REAL :: totcat_an_neigh,totan_cat_neigh

  IF(neighfreq == 1) frnorm = nframes
  IF(neighfreq .NE. 1) frnorm = nframes/neighfreq + 1

  totcat_an_neigh = 0.0; totan_cat_neigh = 0.0

  IF(catan_neighcalc_flag) THEN
!$OMP PARALLEL DO REDUCTION(+:totcat_an_neigh,totan_cat_neigh) PRIVATE(i)

     DO i = 1,maxneighsize

        totcat_an_neigh = totcat_an_neigh + REAL(cat_an_neighavg(i))
        totan_cat_neigh = totan_cat_neigh + REAL(an_cat_neighavg(i))

     END DO

!$OMP END PARALLEL DO

     IF(catan_neighcalc_flag) THEN

        dum_fname = "catanneigh_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace")


     END IF


     DO i = 1,maxneighsize

        WRITE(dumwrite,'(I0,1X,4(F14.8,1X))') i-1,&
             & REAL(cat_an_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(cat_an_neighavg(i))/totcat_an_neigh&
             &,REAL(an_cat_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(an_cat_neighavg(i))/totan_cat_neigh

     END DO

     CLOSE(dumwrite)
     
  END IF

END SUBROUTINE OUTPUT_ALLNEIGHBORS

!--------------------------------------------------------------------


SUBROUTINE OUTPUT_ALLRDF()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm,acrnorm

  IF(rdfcalc_flag) THEN

     rdffrnorm = INT(nfrcntr/rdffreq)
     rvolavg = rvolavg/REAL(rdffrnorm)
     PRINT *, "Average volume of box", rvolavg

     IF(rdfcalc_flag) THEN
        dum_fname = "rdf_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace",iostat=ierr)

        IF(ierr /= 0) THEN
           PRINT *, "Could not open", trim(dum_fname)
        END IF

        WRITE(dumwrite,'(A,8X)',advance="no") "r"

        DO j = 1,npairs

           WRITE(dumwrite,'(I0,A1,I0,8X)',advance="no") pairs_rdf(j&
                &,1),'-',pairs_rdf(j,2)

        END DO

        WRITE(dumwrite,*)

        DO i = 0,rmaxbin-1

           rlower = real(i)*rbinval
           rupper = rlower + rbinval
           nideal = vconst*(rupper**3 - rlower**3)

           WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*rbinval&
                &*(REAL(2*i+1))

           DO j = 1,npairs

              WRITE(dumwrite,'(F16.9,1X)',advance="no")rdfarray(i,j)&
               &/(rdffrnorm*nideal)

           END DO

           WRITE(dumwrite,*)

        END DO

        CLOSE(dumwrite)

     END IF

  END IF

END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_TOPO_ARRAYS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ! Allocate global arrays

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(masses(ntotatomtypes,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate masses"

END SUBROUTINE ALLOCATE_TOPO_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_BOX_ARRAYS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ! Allocate box arrays
  IF (box_from_file_flag == 1) THEN
     ALLOCATE(box_arr(nbox_steps,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate box_arr"
  ELSE
     ALLOCATE(box_arr(1,1),stat = AllocateStatus)
     DEALLOCATE(box_arr)
  END IF

END SUBROUTINE ALLOCATE_BOX_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()

  USE PARAMS_CP2K
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ! Allocate for statics

  IF(rdfcalc_flag) THEN
     ALLOCATE(rdfarray(0:rmaxbin-1,npairs),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdfarray"
  ELSE
     ALLOCATE(rdfarray(1,1),stat = AllocateStatus)
     DEALLOCATE(rdfarray)
  END IF

  
  IF(catan_neighcalc_flag) THEN
     ALLOCATE(cat_an_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate cat_an_neighavg"
     ALLOCATE(an_cat_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate an_cat_neighavg"
  ELSE
     ALLOCATE(cat_an_neighavg(1),stat = AllocateStatus)
     ALLOCATE(an_cat_neighavg(1),stat = AllocateStatus)
     DEALLOCATE(cat_an_neighavg)
     DEALLOCATE(an_cat_neighavg)
  END IF

  
!!$! Allocate for dynamics 
!!$
!!$  IF(ion_dynflag) THEN
!!$     ALLOCATE(itrx_lmp(ioncnt,nframes),stat = AllocateStatus)
!!$     IF(AllocateStatus/=0) STOP "did not allocate itrx_lmp"
!!$     ALLOCATE(itry_lmp(ioncnt,nframes),stat = AllocateStatus)
!!$     IF(AllocateStatus/=0) STOP "did not allocate itry_lmp"
!!$     ALLOCATE(itrz_lmp(ioncnt,nframes),stat = AllocateStatus)
!!$     IF(AllocateStatus/=0) STOP "did not allocate itrz_lmp"
!!$  ELSE
!!$     ALLOCATE(itrx_lmp(1,nframes),stat = AllocateStatus)
!!$     ALLOCATE(itry_lmp(1,nframes),stat = AllocateStatus)
!!$     ALLOCATE(itrz_lmp(1,nframes),stat = AllocateStatus)
!!$     DEALLOCATE(itrx_lmp)
!!$     DEALLOCATE(itry_lmp)
!!$     DEALLOCATE(itrz_lmp)
!!$  END IF
!!$
!!$  IF(cion_dynflag) THEN
!!$     ALLOCATE(ctrx_lmp(c_ioncnt,nframes),stat = AllocateStatus)
!!$     IF(AllocateStatus/=0) STOP "did not allocate ctrx_lmp"
!!$     ALLOCATE(ctry_lmp(c_ioncnt,nframes),stat = AllocateStatus)
!!$     IF(AllocateStatus/=0) STOP "did not allocate ctry_lmp"
!!$     ALLOCATE(ctrz_lmp(c_ioncnt,nframes),stat = AllocateStatus)
!!$     IF(AllocateStatus/=0) STOP "did not allocate ctrz_lmp"
!!$  ELSE
!!$     ALLOCATE(ctrx_lmp(1,nframes),stat = AllocateStatus)
!!$     ALLOCATE(ctry_lmp(1,nframes),stat = AllocateStatus)
!!$     ALLOCATE(ctrz_lmp(1,nframes),stat = AllocateStatus)
!!$     DEALLOCATE(ctrx_lmp)
!!$     DEALLOCATE(ctry_lmp)
!!$     DEALLOCATE(ctrz_lmp)
!!$  END IF

  PRINT *, "Successfully allocated memory for analyis"

END SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMS_CP2K

  IMPLICIT NONE

  !Global arrays
  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(box_arr)

  !Statics calculations arrays
  IF(rdfcalc_flag) DEALLOCATE(rdfarray)

  !Dynamic calculations arrays
!!$  IF(ion_dynflag) THEN
!!$     DEALLOCATE(itrx_lmp)
!!$     DEALLOCATE(itry_lmp)
!!$     DEALLOCATE(itrz_lmp)
!!$  END IF
!!$
!!$  IF(cion_dynflag) THEN
!!$     DEALLOCATE(ctrx_lmp)
!!$     DEALLOCATE(ctry_lmp)
!!$     DEALLOCATE(ctrz_lmp)
!!$  END IF

END SUBROUTINE DEALLOCATE_ARRAYS

!----------------------------------------------------------------
