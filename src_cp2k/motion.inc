&MOTION

  &${RTYPE}

    @IF ( ${RTYPE} == MD ) ! Running molecular dynamics
      ENSEMBLE ${MD_ENS}
      STEPS ${MD_STEPS}
      TIMESTEP ${MD_DT}
      TEMPERATURE ${MD_TEMP}

      @IF ( ${MD_ENS} == NVT )
        &THERMOSTAT
          TYPE NOSE
          REGION GLOBAL
          
          ! Parameters of the Nosé-Hoover chain thermostat
          &NOSE
            TIMECON 50.
            LENGTH 3
            YOSHIDA 3
            MTS 2
          &END
        &END
      @ENDIF

      @IF ( ${MD_ENS} == NPT_I )
        &THERMOSTAT
          TYPE NOSE
          REGION GLOBAL
          
          ! Parameters of the Nosé-Hoover chain thermostat
          &NOSE
            TIMECON 100.
            LENGTH 3
            YOSHIDA 3
            MTS 2
          &END

        &END

        &BAROSTAT
          PRESSURE [bar] 1
          TIMECON 100
        &END


      @ENDIF

    @ENDIF

    @IF ( ${RTYPE} == GEO_OPT ) ! Running geometry optimization
      OPTIMIZER ${GEO_MINI}
      MAX_ITER ${GEO_MAXS}
    @ENDIF
  &END ${RTYPE}
  
  ! Setting up how to print results
  &PRINT
    
    &CELL
      &EACH
        ${RTYPE} ${OUT_STEPS}
      &END EACH
      ADD_LAST NUMERIC
    &END CELL

    &RESTART
      LOG_PRINT_KEY T
      &EACH
        ${RTYPE} ${OUT_STEPS} ! Restart files will be printed every ${OUT_STEPS} steps
      &END EACH
      ADD_LAST NUMERIC
    &END RESTART

    &TRAJECTORY
      LOG_PRINT_KEY T
      FORMAT ${OUT_FORM}
      UNIT ${OUT_UNIT}
      FILENAME ./${TRAJFILE}
      &EACH
        ${RTYPE} ${OUT_STEPS} ! Trajectory files will be printed every ${OUT_STEPS} steps
      &END EACH
      ADD_LAST NUMERIC
    &END TRAJECTORY

    @IF ( ${RTYPE} == MD )
      &VELOCITIES
        FORMAT ${OUT_FORM}
        UNIT ${OUT_UNIT}
        &EACH
          ${RTYPE} ${OUT_STEPS} ! Velocity files will be printed every ${OUT_STEPS} steps
        &END EACH
        ADD_LAST NUMERIC
      &END VELOCITIES

      @IF ( ${MD_ENS} == NPT_I )
        &STRESS
          &EACH
            MD 5
          &END
          FILENAME pressure
        &END STRESS
      @ENDIF

    @ENDIF

  &END PRINT

&END MOTION
