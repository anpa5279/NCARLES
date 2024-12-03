PROGRAM les_mpi

  USE pars
  USE fields
  USE con_data
  USE con_stats
  USE tracerbc, ONLY: applytracerbc

  INCLUDE 'mpif.h'

  ! INITIALIZE MPI, GET MYID, NUMPROCS
  ! TEST IF ON ROOT PROCESS

  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world,myid,ierr)
  CALL mpi_comm_size(mpi_comm_world,numprocs,ierr)

  i_root = 0
  l_root = .FALSE.
  IF(myid == i_root) l_root = .TRUE.

  l_debug = .FALSE.
  IF(idebug == 1) l_debug = .TRUE.

  ts_mpi = mpi_wtime()

  ! SET NUMBER OF X-Y SLAB CPUS
  ncpu_s   = 8
  itn      = 0
  case_inp = '30L'

  CALL get_units
  CALL gridd
  CALL setcon
  CALL set_paths
  istop = 1
  
  ! SCRATCH RUN
  IF (iti==0)  THEN
    igrdr = 2
    case = case_inp
    CALL init
    CALL setup(it)

    ! CHOOSE ROUTINE FOR GETTING INITIAL GUESS
    CALL randoc
    CALL get_max
  ELSE
    igrdr = 3
    CALL restart
    CALL get_max
    CALL setup(it)
  ENDIF

  !open input file and collect information (for now will be hardcoded instead)
  model_type = 2
  nscl = 8, nvar = (4+nscl) !number of scalars and vars


  !apply tracer intial conditions according to the input file
  CALL applytracerbc

  ! TIME LOOP
  tzero = time
  CALL get_dt(it,iti)

  DO WHILE(it<itmax)
    CALL set_sav(it,iti)

    ! UPDATE POSITION OF VORTEX
    IF(it >= new_vis .AND. ivis0 == 1) THEN
      ivis = 1
    ELSE
      ivis = 0
    ENDIF

    ! 3 STAGE RUNGE-KUTTA TIME STEPPING
    DO istage=1,3
      dtzeta = dt*zetas(istage)
      dtgama = dt*gama(istage)

      ! COMPUTE DERIVATIVES OF (U,V,W)
      CALL exchange
      CALL get_derv

      ! NEW EDDY VISCOSITY, AND BCS
      IF(iss == 0 .AND. ifree == 0) THEN
        CALL lower(it)
      ELSEIF(ifree == 1) THEN
        CALL lower_free(it)
      ENDIF

      IF(ise == numprocs-1) THEN
        CALL upper
      ENDIF

      CALL bcast_pbc
      CALL get_means(istage)

      IF(ivis == 1) THEN
        CALL iso(it)
        CALL surfvis(it)
      ENDIF

      IF(istage == 1)THEN
        CALL xy_stats
        CALL tke_budget
        CALL pbltop(itop)
      ENDIF

      ! GET RHS FOR ALL EQUATIONS
      CALL comp1(istage,it)

    ! SOLVE FOR PRESSURE
      CALL comp_p

    ! ADD PRESSURE GRADIENT AND DEALIAS
      CALL comp2

      ! update sources and sinks of chemical reactions
      CALL wrapper(model_type, tracer_matrix, required_cond_matrix, required_cond_array)


      IF(msave .AND. istage == 3) THEN
        CALL save_v(it)
      ENDIF

      IF(istage == 3) THEN
        IF(msave .AND. l_root) CALL save_c(it)
      ENDIF

      IF(micut) THEN
        CALL dealias
      ENDIF

      IF(mnout .AND. istage == 1)  THEN
        IF(l_debug) THEN
          CALL print(nprt,it,izs,ize)
        ENDIF
        IF(l_root) CALL print(6,it,1,nnz)
      ENDIF

      IF(l_root) THEN
        IF(mhis  .AND. istage == 1)  CALL write_his(itop)
        IF(mhis  .AND. istage == 1 .AND. mtape) CALL close_his
      ENDIF

    END DO

    CALL get_max
    CALL get_dt(it,iti)
  END DO !end of while loop
  close(10)

  te_mpi = mpi_wtime()

  WRITE(6,9997) (te_mpi - ts_mpi)

  CONTINUE
  CALL mpi_finalize(ierr)

! FORMAT
9997  FORMAT(' Job Execution Time = ',e15.6)

  STOP
END PROGRAM
