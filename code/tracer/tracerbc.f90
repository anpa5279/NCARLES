module tracerbc
  use con_data, only: dz,dx,dy, zl
  use con_stats, only: z, zz
  use pars, only: flg_alk, flg_npz, iys,iye,izs,ize,izi,nnx,nnz,nny,nscl
  use fields, only: t
  use inputs

  implicit none
  integer, parameter :: flg_debug = 0
  real, dimension(nscl) :: tau,airval
  integer, dimension(nscl) :: ictype,rmodel,rdorg,rpartner,asflux
  integer, dimension(2,nscl) :: bnd
  integer, dimension(2) :: bnds
  integer, dimension(3,nscl) :: point
  integer, dimension(3) :: points
  real, dimension(nscl) :: val
  real, dimension(nscl) :: chng
  contains

! iscl      : scalar number (temperature is always iscl=1)
! tau       : reaction time scale
! ictype    : initial condition (0 = nothing, 1 = horiz. band,
!                                2 = vertical band in x, 3 = vertical band in y,
!                                4 = point source, 5 = vertical gradient,
!                                6 = horiz. gradient in x, 7 = horiz. gradient in y
!                                8= exponential decay in -z, 9= exponential decay in +z
!                                10= exponential decay in -z with a constant layer)
!             ictype does not work for iscl=1 (temperature), that is set in init/randoc.f
! val       : value of initial finite or source band/point
! np        : width of initial finite or source band
! zt        : upper/left most level or finite or source band
! bndz       :
! point     :
! rmodel    : reaction model type (0 = no reaction, 1 = single tracer decay/growth,
!                                  2 = two tracers decay/growth, 3 = carbonate chemistry)
! rdorg     : reaction decay or growth (0 = decaying tracer, 1 = growing tracer)
! rpartner  : reaction partner (iscl number of coupled tracer for reaction,
!                               0 = no coupled tracer)
! asflux    : air-sea flux boundary condition (0 = for no flux, 1 = for flux [also need
!             flag_airseaflux.eq.1 in pars.f])
! airval    : value of tracer in air (only need set for use with asflux and flag_airseaflux)

    subroutine applytracerbc(it)
      integer, intent(in) :: it
      integer :: iscl, np, zt
      real :: ta, vals, zl

            !! active tracers (temperature)
            iscl = 1; 
            ictype(iscl) = 0;   val(iscl) = 273.15 + iTsurf;      tau(iscl)    = 0;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = 0;      zt = 0;  rmodel(iscl) = 0;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            !! passive tracers
            iscl = 2;!carbon dioxide (CO2)
            ictype(iscl) = 1;   val(iscl) = c1;     tau(iscl)      = 1;
            asflux(iscl) = 1;   airval(iscl) = 8.56056;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 3;!bicarbonate (HCO3)
            ictype(iscl) = 1;   val(iscl) = c2;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 4;!carbonate (CO3)
            ictype(iscl) = 1;   val(iscl) = c3;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 5;!Boric acid (B(OH)3)
            ictype(iscl) = 1;   val(iscl) = c4;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 6; !Tetrahydroxyborate (B(OH)4)
            ictype(iscl) = 1;   val(iscl) = c5;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 7; !hydrogen ion (H+)
            ictype(iscl) = 1;   val(iscl) = c6; tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 8; !Hydroxl ion (OH-)
            ictype(iscl) = 1;   val(iscl) = c7;     tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0;

            iscl = 9; !phytoplankton (C106H175O42N16P, organic matter)
            ictype(iscl) = 1;   val(iscl) = c8;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0.02;

            iscl = 10; !zooplankton
            ictype(iscl) = 1;   val(iscl) = c9; tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0.02;

            iscl =11; !nitrate (NO3)
            ictype(iscl) = 1;   val(iscl) = c10;     tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0.02;

            iscl =12; !detritus
            ictype(iscl) = 1;   val(iscl) = c11;     tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
            chng(iscl)=0.01;

        do iscl = 2,nscl
          bnds=bnd(:,iscl); vals=val(iscl); points=point(:,iscl);
          if (it.eq.1) then !not sure why this if statement is here. it would be more efficient if this would in les_mpi line 89. 
              if (ictype(iscl).eq.1) call hbndsource(iscl,bnds,vals);
              if (ictype(iscl).eq.4) call pointsource(iscl, bnds, vals);
              if (ictype(iscl).eq.5) call vgradsource(iscl,bnds,vals);
              if (ictype(iscl).eq.8) call zdecay(iscl, chng(iscl), bnds, vals);
              if (ictype(iscl).eq.9) call nutrients(iscl, chng(iscl), bnds, vals);
              if (ictype(iscl).eq.10) call nutrients1(iscl, chng(iscl), bnds, vals);
          endif
          bnds = 0; vals = 0; points = 0;
        enddo
     

      if(flg_debug == 1) then
          open(13, file='tracerbc.txt',access='append')
          write(13,'(A)') '------------------------'
          write(13,'(A,i3)') 'RUNNING FOR IT= ',it
          write(13,'(A,f9.6)') 'Z for 5m above H', z(izi)+5.0
          write(13,'(A,f9.6)') 'Z for 5m below H', z(izi)-5.0
          close(13)
      end if

    end subroutine

    subroutine hbndsource(iscl, bnds, vals)
      integer, intent(in) :: iscl
      integer, intent(in), dimension(2) :: bnds
      real, intent(in) :: vals
      integer :: ix,iy,iz
      do iz=bnds(1),bnds(2)
         do iy=iys,iye
            do ix=1,nnx
               if ((iz >= izs) .and. (iz <= ize)) then
                     t(ix,iy,iscl,iz) = vals
               endif
            end do
         end do
      end do

    end subroutine

    subroutine vgradsource(iscl, bnds, vals)
      integer, intent(in) :: iscl
      integer, intent(in), dimension(2) :: bnds
      real, intent(in) :: vals
      integer :: ix,iy,iz,zi

      zi  = z(bnds(2))
      do iy=iys,iye
         do iz=bnds(1),bnds(2)
            do ix=1,nnx
               if ((iz >= izs) .and. (iz <= ize)) then
                  t(ix,iy,iscl,iz) = (vals/zi)*(zi-zz(iz))
               endif
            end do
         end do
      end do

    end subroutine

    subroutine zdecay(iscl, chng, bnds, vals)
      integer, intent(in) :: iscl
      real, intent(in) :: chng
      real, intent(in) :: vals
      integer, intent(in), dimension(2) :: bnds
      integer :: ix,iy,iz

      do iy=iys,iye
         do iz=bnds(1),bnds(2)
            do ix=1,nnx
              if ((iz >= izs) .and. (iz <= ize)) then
                t(ix,iy,iscl,iz) =vals*exp(chng*z(iz-1))
              endif
            end do
         end do
      end do
    end subroutine

    subroutine nutrients(iscl, chng, bnds, vals)
      integer, intent(in) :: iscl
      real, intent(in) :: chng
      real, intent(in) :: vals
      integer, intent(in), dimension(2) :: bnds
      integer :: ix,iy,iz, nnz

      do iy=iys,iye
         do iz=bnds(1),bnds(2)
            do ix=1,nnx
              if ((iz >= izs) .and. (iz <= ize)) then
                t(ix,iy,iscl,iz) =vals*exp(-chng*(z(iz)-zl))
              endif
            end do
         end do
      end do
    end subroutine

    subroutine nutrients1(iscl, chng, bnds, vals)
      integer, intent(in) :: iscl
      real, intent(in) :: chng
      real, intent(in) :: vals
      integer, intent(in), dimension(2) :: bnds
      integer :: ix,iy,iz, nnz

      do iy=iys,iye
         do iz=bnds(1),bnds(2)
            do ix=1,nnx
              if ((iz >= izs) .and. (iz <= ize)) then
                t(ix,iy,iscl,iz) =vals*exp(-chng*(z(iz+1)-zl))
              endif
            end do
         end do
      end do
    end subroutine


    subroutine pointsource(iscl, bnds, vals)
      integer, intent(in) :: iscl
      integer, intent(in), dimension(2) :: bnds
      real, intent(in) :: vals
      real :: spread=1
      integer :: ix,iy,iz,zi
      
      !point source at the surface in the middle
      do iy=iys,iye
        do iz=bnds(1),bnds(2)
            do ix=1,nnx
              if ((iz >= izs) .and. (iz <= ize)) then
                t(ix, iy, iscl,iz)=2*vals/((2*4.0*ATAN(1.0))**(3/2)*spread**3)*exp(-((ix-nnx/2)**2+(iy-nny/2)**2+(iz)**2)/(2*spread**2))
              endif
            end do
        end do
      end do
    end subroutine

    function znptobnd(zt,np)
      integer, intent(in) :: zt
      integer, intent(in) :: np
      integer, dimension(2) :: znptobnd
      integer :: iz

      ! set the first bound, and make sure it doesn't exceed dimensions
      iz = ztoiz(zt)
      znptobnd(1) = iz - int((np-1)/2)
      if (znptobnd(1) < 0) then
        znptobnd(1) = 0
      end if

      ! set the second bound based upon the first
      znptobnd(2) = znptobnd(1) + np -1

    end function

    function xnptobnd(xt,np)
      integer, intent(in) :: xt
      integer, intent(in) :: np
      integer, dimension(2) :: xnptobnd
      integer :: ix

      ! set the first bound, and make sure it doesn't exceed dimensions
      ix = xtoix(xt)
      xnptobnd(1) = ix - int((np-1)/2)
      if (xnptobnd(1) < 0) then
        xnptobnd(1) = 0
      end if

      ! set the second bound based upon the first
      xnptobnd(2) = xnptobnd(1) + np -1

    end function

    function ynptobnd(yt,np)
      integer, intent(in) :: yt
      integer, intent(in) :: np
      integer, dimension(2) :: ynptobnd
      integer :: iy

      ! set the first bound, and make sure it doesn't exceed dimensions
      iy = ytoiy(yt)
      ynptobnd(1) = iy - int((np-1)/2)
      if (ynptobnd(1) < 0) then
        ynptobnd(1) = 0
      end if

      ! set the second bound based upon the first
      ynptobnd(2) = ynptobnd(1) + np -1

    end function

    function restobnd(zt,dr)
      integer,intent(in) :: zt
      integer, intent(in) :: dr
      integer, dimension(2) :: restobnd
      integer :: iz

      iz = ztoiz(zt)
      if (dr > 0) then ! surface res
        restobnd(1) = 0
        restobnd(2) = iz
      else
        restobnd(1) = 0
        restobnd(2) = nnz
      end if
    end function

    function ztoiz(zt)
      integer, intent(in) :: zt
      integer :: ztoiz

      ! note that this will only work for equispaced z grids
      ztoiz = int(zt/dz)

    end function

    function xtoix(xt)
      integer, intent(in) :: xt
      integer :: xtoix

      ! note that this will only work for equispaced z grids
      xtoix = int(xt/dx)

    end function

    function ytoiy(yt)
      integer, intent(in) :: yt
      integer :: ytoiy

      ! note that this will only work for equispaced z grids
      ytoiy = int(yt/dy)

    end function
    
end module
