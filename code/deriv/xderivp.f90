SUBROUTINE xderivp(ax,trigx,xk,nnx,iys,iye)

! GET MULTIPLE X DERIVATIVES USING FFTPACK ROUTINES
! USE FFTPACK STORAGE A0, (A1,B1), (A2,B2), ...
! ASSUMES THAT FIRST WAVENUMBER XK(1) = 0.0
! ASSUMES WAVENUMBERS ARE NORMALIZED BY NUMBER OF POINTS

  REAL :: xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)

  DO iy=iys,iye
    CALL rfftf(nnx,ax(1,iy),trigx)

    ii = 1
    ax(1,iy) = 0.0
    ax(nnx,iy) = 0.0
    DO ix=2,nnx-1,2
      ii          = ii + 1 !why is this even needed
      temp        = ax(ix,iy)
      ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
      ax(ix+1,iy) = xk(ii)*temp
    ENDDO

    CALL rfftb(nnx,ax(1,iy),trigx)
  ENDDO
  RETURN
END SUBROUTINE xderivp
