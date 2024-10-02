MODULE pars

  IMPLICIT NONE

! PARAMETERS ARE DESCRIBED IN PARAMS.pdf
! nslab=nproc/ncpu_s=256/24=6
! nys=ny/nslab=256/6=24

  INTEGER, PARAMETER :: flg_stokes = 1    ! stokes on or off
  INTEGER, PARAMETER :: flg_lat = 0       ! solve using lat or utau
  INTEGER, PARAMETER :: flg_reaction = 1  ! 3.1536e8 reaction model on or off
  INTEGER, PARAMETER :: chem0d = 0
  INTEGER, PARAMETER :: co2_asflux = 0    ! 2 => WB_param
  INTEGER, PARAMETER :: flg_alk = 0       ! changing alkalinity
  INTEGER, PARAMETER :: flg_npz = 1       ! changing NPZ model, mayzaud-poulet= 0 or ivlev=1
  INTEGER, PARAMETER :: iti=0, itmax=80000, imean=1, ihst=01, itape=500,        &
  itstr=1, it_his=120000, i_viz=120000

  INTEGER, PARAMETER :: nscl = 12, nvar = (4+nscl) !number of scalars and vars
  INTEGER, PARAMETER :: nxg1  = 32, nyg1  = 32, nzg1  = 32 !size of problem
  INTEGER, PARAMETER :: maxnx = 256, maxny = 256, maxnz = 256 !max size
  INTEGER, PARAMETER :: maxnz1 = maxnz + 1, maxnz2 = maxnz + 2,             &
   maxnx2 = maxnx + 2, maxny2 = maxny + 2


! INTERNAM FLAG DEFINITIONS
!       igrdr   =  3; data comes from restart file
!               =  2; data comes from initialization (random)
!               =  1; data comes from coarser grid (or otherwise)
!
!       ibcu    = 1 ; upper boundary condition set by radiation bc
!               = 0 ; fixed value = 0.
!               = -1; value defined by coarser mesh for all variables
!
!       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
!               = -1; value defined by coarser mesh for all variables
!
!       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
!               = 1 ; fixed time step set in sr. get_dt
!
!       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
!               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
!
!       ihst    = nn; frequency at which global variables are output in history file
!               < 0 ; no history files
!
!       it_his  = time step where history files start, incremented by itape
!
!       ismlt   = 1 ; use businger formulas in MO
!                 0 ; use large and everyone elses formulas in MO
!
!       iupwnd  = 0;  use skew symmetric formulas for all derivatives
!                    in scalar equations
!              = 1;  use hybrid upwind scheme for all derivatives
!                    in scalar equations
!
!      ivis0   = 0; old eddy viscosity model
!              = 1; new eddy viscosity model
!
!      new_vis = step; the iteration step for which the new model
!                      is turned on when ivis0=1
!              < 0; new model is on at all steps for ivis0=1
!
!      nscl  .ge. 1   number of scalars to be followed set in parameter statements
!                     change entries in sr. init, and sr. suft for surface bc's


  INTEGER, PARAMETER ::                                                     &
          ismlt=0, ifree=0, isfc=0, iradup=0, iupwnd=1, ibuoy=1, itcut=1,   &
          method=3, idebug=0, iz_space=0, ivis0=1, ifix_dt=0, new_vis=-1,   &
          i_dear = 0
  INTEGER, PARAMETER :: j_recl=4 !record length in "bytes" for history file
  INTEGER, PARAMETER :: k8=8 !kind parameter for integers in mpi_io routines

!
! -------- number of vars, size of problem in (x,y,z), max size, max processors
!
!
! ------------ leave following definitions as is
!
! ----------------------------------------------------------------------
  INTEGER ::                                                                &
          nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy, nnxp2, nnyp2, nnzp1,  &
          ivis, nnzm1, isize, krec, izs, ize, ixs, ixe, jxs, jxe, kxs, kxe, &
          mxs, mxe, iss, ise, iys, iye, jys, jye
  INTEGER ::                                                                &
          it_counter
!----------------------------------------------------------------------
  CHARACTER case*3
!----------------------------------------------------------------------
  INTEGER ::                                                                &
          nvel, npre, nhis1, nprt, nhisp, nvelc, nviz_z, nviz_y, nviz_x,    &
          nviz_s, kfile, jfile, ibcl, ibcu, igrdr, imach, itn, it_nxt
  LOGICAL ::                                                                &
          mnout, micut, mtape, mhis, msave, l_root, l_debug, msave_v, mviz
!----------------------------------------------------------------------
  REAL ::                                                                   &
          windm, u1xy, v1xy, t1xy(nscl), t10xy(nscl), au13m, au23m,         &
          aut3m(nscl), tsfcm(nscl), thstar(nscl), eavg(maxnz),              &
          tr_tau(0:maxnz), zi_min, dir_x, dir_y
  INTEGER ::                                                                &
          izi, iz_min
  REAL, ALLOCATABLE ::                                                      &
          wind(:,:), tau13m(:,:), tau23m(:,:), taut3m(:,:,:), t_grnd(:,:,:)
!----------------------------------------------------------------------
  REAL ::                                                                   &
         u_mn(0:maxnz1), v_mn(0:maxnz1), w_mn(0:maxnz1), t_mn(0:maxnz1,nscl)
!----------------------------------------------------------------------
  REAL ::                                                                   &
         dzw(0:maxnz2), dzu(0:maxnz2), dzw_i(0:maxnz2), dzu_i(0:maxnz2)
!----------------------------------------------------------------------
  REAL ::                                                                   &
         t_factor, t_ref
!----------------------------------------------------------------------
  REAL ::                                                                   &
         dfac(maxnz), dsl_z(0:maxnz1), xksurf, viscon, vise, almin_c,       &
         stabmin, ck, ceps, csmag, stab_c, vis_mean(0:maxnz)
  INTEGER ::                                                                &
        nmatch
!----------------------------------------------------------------------
  REAL ::                                                                   &
          zetas(3), gama(3), etas(4), dt_new, umax,vmax,wmax, wabs, emax,   &
          vismax, cfl, tzero, ucfl, vcfl, wcfl, s10_est
!----------------------------------------------------------------------
  CHARACTER*80                                                              &
          path_res, path_sav, path_his, path_prt, path_hp, path_sav_hp,     &
          path_v, path_c, path_p, path_h, path_sav_v, path_sav_c,           &
          path_sav_p, path_sav_h, bad_news
  CHARACTER case_inp*3
  CHARACTER*80 path_viz_xy, path_viz_xz, path_viz_yz, path_stuf
!----------------------------------------------------------------------
  INTEGER ::                                                                &
          myid, numprocs, i_root, ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
  INTEGER, ALLOCATABLE, DIMENSION(:) ::                                     &
          ix_s, ix_e, jx_s, jx_e, kx_s, kx_e, mx_s, mx_e, iy_s, iy_e, jy_s, &
          jy_e, is_s, is_e, iz_s, iz_e
!----------------------------------------------------------------------
  REAL ::                                                                &
          fug, khen, bet_ost, kbub, u_tau, wa, wh
!----------------------NPZ parameters that are called in the code----------------
  REAL ::                                                                &
          k_ext=0.1
END MODULE
