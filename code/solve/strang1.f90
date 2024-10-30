SUBROUTINE strang1(it)
! STRANG SPLITTING OF SCALAR REACTION TERM - 0.5*REACT, ADVECT, 0.5*REACT FOR
! FAST REACTIONS (TAU <= 1000), SEE RHS/RHS_SCL.F90 FOR SLOW REACTION SOURCES

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats
  USE reaction, ONLY: react_src

  INCLUDE 'mpif.h'

  ! TEMP SCALAR ARRAY TO HOLD RHS FROM STEP N-1 AND FIELD VARIABLES FROM STEP N
  REAL :: trhs(nnx,iys:iye,nscl,izs:ize)
  REAL, DIMENSION(nscl-1) :: tmp
  REAL :: Navg, Pavg, Zavg

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
      tmp = react_src(ix,iy,1,iz) !why even bother indexing for the first scalar when it is not even used in react_src?
        DO l=2,nscl
          trhs(ix,iy,l,iz) = tmp(l-1)
          IF(trhs(ix,iy,l,iz).le.1.0e-20)THEN
            trhs(ix,iy,l,iz) = 1.0e-20
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO iz=izs,ize
    DO l=2,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
        
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !taking averages
  !Pavg=0.0
  !Zavg=0.0
  !Navg=0.0

  !DO iz=izs,ize
  !  DO ix=1,nnx
  !    DO iy=iys,iye
  !      Pavg= Pavg+ t(ix,iy,6,iz)/(nnx**2)
  !      Zavg= Zavg+ t(ix,iy,7,iz)/(nnx**2)
  !      Navg= Navg+ t(ix,iy,8,iz)/(nnx**2)
  !      
  !    ENDDO
  !  ENDDO
  !  CALL NPZdebug(Pavg, Zavg, Navg, iz, 'end of strang1')
  !  Pavg=0.0
  !  Zavg=0.0
  !  Navg=0.0
  !ENDDO
  RETURN
END SUBROUTINE
