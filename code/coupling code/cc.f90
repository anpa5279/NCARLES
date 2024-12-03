module cc
  use fields, only: t
  use inputs
  use con_data, only: time,dt, dz
  use pars, only: flg_npz, nscl,chem0d,flg_alk, k_ext
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Parameters

  ! LOCATED IN MODULES/PAR.F
  ! flg_reaction : Flag to turn on or off the reaction models

  ! LOCATED IN TRACER/TRACERBC.F90
  ! rmodel       : An array of models to be used in the reactive tracers. Models are
  !                separated out in the react_src function
  ! rdorg        : An array of 0 or 1, 0 = decaying reaction, 1 = growing reaction
  ! rpartner     : An array dictating what tracers are coupled to each other
  ! tau          : The timescale to be used.


contains

  ! REACT_SRC: calculate the scalar reaction source term for a given scalar
  !            and point. This is called in rhs_scl for each scalar.
function react_src(ix, iy, iscl, iz)
   ! ix, iy, iscl, iz (in): location of interest
   real, dimension(nscl-1) :: react_src
   integer, intent(in) :: ix, iy, iscl, iz

   real, dimension(0:nscl-2) :: co2, co2tmp
   real :: temper, t_rkc, t_end
   real, dimension(0:4+nscl-1) :: workLocal
   integer :: i

   t_rkc  = time
   t_end  = time + dt * 0.5

   !!!!!!!!!!!!
   !! Zeebe et al. 2001 Carbonate Chemistry
   !!!!!!!!!!!!!
   ! c(0) = Carbon Dioxide, [CO2], t(ix,iy,2,iz)
   ! c(1) = Bicarbonate, [HCO3-], t(ix,iy,3,iz)
   ! c(2) = Carbonate, [CO32-], t(ix,iy,4,iz)
   ! c(3) = Boric Acid, [B(OH)3], t(ix,iy,5,iz)
   ! c(4) = Tetrahydroxyborate, [B(OH)4-], t(ix,iy,6,iz)
   ! c(5) = Hydrogen Ion, [H+], t(ix,iy,7,iz)
   ! c(6) = Hydroxide, [OH-], t(ix,iy,8,iz)

   workLocal(:) = 0.0
   co2(0:nscl-2) = t(ix, iy, 2:nscl, iz)

   if (chem0d == 1) then
      temper = iTsurf
   else
      temper = t(ix, iy, 1, iz)
   end if

   co2tmp = rkc_driver(t_rkc, t_end, workLocal, co2, temper, iz)
   co2 = co2tmp

   do i = 0, nscl-2
      react_src(i+1) = co2(i)
   end do

end function react_src

  function rkc_driver(t_rkc, t_end, workLocal, co2, temper, iz)
    ! Driver function for RKC integrator.
    !
    ! t_tkc       the starting time.
    ! t_end       the desired end time.
    ! task        0 to take a single integration step, 1 to integrate to tEnd.
    ! workLocal   Real work array, size 3.
    ! co2         Dependent variable array, integrated values replace initial conditions.

    nstep   = 0
    rel_tol = 1.0e-6
    abs_tol = 1.0e-10
    UROUND  = 2.22e-16
    m_max   = nint(sqrt(rel_tol / (10.0 * UROUND)))
    hmax    = abs(t_end - t_rkc)
    hmin    = 10.0 * UROUND * max(abs(t_rkc), hmax)

    if(m_max < 2)then
       m_max = 2
    endif

    do i = 0,nscl-2
       y_n(i) = co2(i)
    enddo
    

    ! calculate F_n for initial y
    F_n = dydt(t_rkc, y_n, temper, iz)
   
    ! load initial estimate for eigenvector
    if(workLocal(2) < UROUND) then
       do i = 0,nscl-2
          workLocal(4+i) = F_n(i)
       enddo
    endif

    do while (t_rkc < t_end)
       ! use time step stored in workLocal(3)

       ! estimate Jacobian spectral radius
       ! only if 25 steps passed
       ! spec_rad = workLocal(4)
       temp_arr(:)  = 0.0
       temp_arr2(:) = 0.0
       err          = 0.0

       if(mod(nstep,25) == 0)then
          workLocal(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, workLocal(4), temp_arr2, temper, iz)
       endif

       ! first step, estimate step size
       if(workLocal(2) < UROUND)then
          workLocal(2) = hmax
          if((workLocal(3) * workLocal(2)) > 1.0)then
             workLocal(2) = 1.0/workLocal(3)
          endif
          workLocal(2) = max(workLocal(2), hmin)

          do i = 0,nscl-2
             temp_arr(i) = y_n(i) + (workLocal(2) * F_n(i))
          enddo
          
          temp_arr2 = dydt(t_rkc + workLocal(2), temp_arr, temper, iz)
          err = 0.0
          do i = 0,nscl-2
             est = (temp_arr2(i) - F_n(i)) / (abs_tol + rel_tol * abs(y_n(i)))
             err = err + est*est
          enddo
          err = workLocal(2) * sqrt(err/real(nscl-2))

          if((0.1 * workLocal(2)) < (hmax * sqrt(err)))then
             workLocal(2) = max((0.1 * workLocal(2)) / sqrt(err), hmin)
          else
             workLocal(2) = hmax
          endif
       endif

       ! check if last step
       if((1.1 * workLocal(2)) .ge. abs(t_end - t_rkc))then
          workLocal(2) = abs(t_end - t_rkc)
       endif

       ! calculate number of steps
       m = 1 + nint(sqrt(1.54 * workLocal(2) * workLocal(3) + 1.0))

       if(m > m_max)then
          m = m_max
          workLocal(2) = real((m*m - 1) / (1.54*workLocal(3)))
       endif

       hmin = 10.0 * UROUND * max(abs(t_rkc), abs(t_rkc + workLocal(2)))

       ! perform tentative time step
       co2 = rkc_step(t_rkc, workLocal(2), y_n, F_n, m, temper, iz)
       
       ! calculate F_np1 with tenative y_np1
       temp_arr = dydt(t_rkc + workLocal(2), co2, temper, iz)
       
       ! estimate error
       err = 0.0
       do i = 0,nscl-2
          est = 0.0
          est = 0.8 * (y_n(i) - co2(i)) + 0.4 * workLocal(2) * (F_n(i) + temp_arr(i))
          est = est / (abs_tol + rel_tol * max(abs(co2(i)), abs(y_n(i))))
          err = err + est*est
       enddo
       err = sqrt(err / 7.0)

       if (err > 1.0) then
          ! error too large, step is rejected

          ! select smaller step size
          workLocal(2) = 0.8 * workLocal(2) / (err**(1.0/3.0))

          ! reevaluate spectral radius
          workLocal(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, workLocal(4), temp_arr2, temper, iz)
       else
          ! step accepted
          t_rkc = t_rkc + workLocal(2)
          nstep = nstep + 1

          fac   = 10.0
          temp1 = 0.0
          temp2 = 0.0
          if(workLocal(1) < UROUND)then
             temp2 = err**(1.0/3.0)
             if(0.8 < (fac * temp2))then
                fac = 0.8 /  temp2
             endif
          else
             temp1 = 0.8 * workLocal(2) * (workLocal(0)**(1.0/3.0))
             temp2 = workLocal(1) * (err**(2.0/3.0))
             if(temp1 < (fac * temp2))then
                fac = temp1 / temp2
             endif
          endif

          ! set "old" values to those for current time step
          workLocal(0) = err
          workLocal(1) = workLocal(2)

          do i = 0,nscl-2
             y_n(i) = co2(i)
             F_n(i) = temp_arr(i)
          enddo

          ! store next time step
          workLocal(2) = workLocal(2) * max(0.1, fac)
          workLocal(2) = max(hmin, min(hmax, workLocal(2)))

       endif
    enddo

    do i = 0,nscl-2
       rkc_driver(i) = co2(i)
    enddo

  end function rkc_driver

  real function rkc_spec_rad(t_rkc, hmax, co2, F, v, Fv, temper, iz)
    ! Function to estimate spectral radius.
    !
    ! t_rkc    the time.
    ! hmax     Max time step size.
    ! co2   Array of dependent variable.
    ! F        Derivative evaluated at current state
    ! v
    ! Fv

    real, intent(in) :: t_rkc
    real, intent(in) :: hmax, temper
    real, intent(inout), dimension(0:nscl-2) :: v, Fv, F
    real, intent(in), dimension(0:nscl-2) :: co2
    integer itmax, i, iter, ind
    real UROUND, small, nrm1, nrm2, dynrm, sigma
    integer, intent(in) :: iz

    UROUND  = 2.22e-16
    itmax   = 50
    small   = 1.0 / hmax
    nrm1    = 0.0
    nrm2    = 0.0
    sigma   = 0.0

    do i = 0,nscl-2
       nrm1 = nrm1 + co2(i) * co2(i)
       nrm2 = nrm2 + v(i) * v(i)
    enddo
    nrm1 = sqrt(nrm1)
    nrm2 = sqrt(nrm2)

    if((nrm1 .ne. 0.0) .and. (nrm2 .ne. 0.0))then
       dynrm = nrm1 * sqrt(UROUND)
       do i = 0,nscl-2
          v(i) = co2(i) + v(i) * (dynrm / nrm2)
       enddo
    elseif(nrm1 .ne. 0.0)then
       dynrm = nrm1 * sqrt(UROUND)
       do i = 0,nscl-2
          v(i) = co2(i) * (1.0 + sqrt(UROUND))
       enddo
    elseif(nrm2 .ne. 0.0)then
       dynrm = UROUND
       do i = 0,nscl-2
          v(i) = v(i) * (dynrm / nrm2)
       enddo
    else
       dynrm = UROUND
       do i = 0,nscl-2
          v(i) = UROUND
       enddo
    endif

    ! now iterate using nonlinear power method
    sigma = 0.0
    do iter = 1,itmax
       Fv = dydt(t_rkc, v, temper, iz)

       nrm1 = 0.0
       do i = 0,nscl-2
          nrm1 = nrm1 + ((Fv(i) - F(i)) * (Fv(i) - F(i)))
       enddo
       nrm1  = sqrt(nrm1)
       nrm2  = sigma
       sigma = nrm1 / dynrm
       if((iter .ge. 2) .and. (abs(sigma - nrm2) .le. (max(sigma, small) * 0.01)))then
          do i = 0,nscl-2
             v(i) = v(i) - co2(i)
          enddo
          rkc_spec_rad = 1.2 * sigma
       endif

       if(nrm1 .ne. 0.0)then
          do i = 0,nscl-2
             v(i) = co2(i) + ((Fv(i) - F(i)) * (dynrm / nrm1))
          enddo
       else
          ind = mod(iter, int(nscl-1))
          v(ind) = co2(ind) - (v(ind) - co2(ind))
       endif
       !CALL NPZdebug(v(4), v(5), v(6), iz, 'end of rkc_spec_rad')
       
    enddo

    rkc_spec_rad = 1.2 * sigma

  end function rkc_spec_rad

  function rkc_step(t_rkc, h, y_0, F_0, s, temper, iz)
    ! Function to take a single RKC integration step
    !
    ! t_rkc    the starting time.
    ! h        Time-step size.
    ! y_0      Initial conditions.
    ! F_0      Derivative function at initial conditions.
    ! s        number of steps.
    ! rkc_step Integrated variables

    real, intent(in) :: t_rkc
    real, intent(in) :: h, temper
    real, intent(inout), dimension(0:nscl-2) :: y_0, F_0
    integer, intent(in) :: s
    real, dimension(0:nscl-2) :: rkc_step
    real, dimension(0:nscl-2) :: y_j
    real w0, temp1, temp2, arg, w1, b_jm1, b_jm2, mu_t
    real c_jm2, c_jm1, zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2
    real zj, dzj, d2zj, b_j, gamma_t, nu, mu, c_j
    real, dimension(0:nscl-2) :: y_jm1, y_jm2
    integer i, j
    integer, intent(in) :: iz

    w0    = 1.0 + 2.0 / (13.0 * real(s * s))
    temp1 = (w0 * w0) - 1.0
    temp2 = sqrt(temp1)
    arg   = real(s) * log(w0 + temp2)
    w1    = sinh(arg) * temp1 / (cosh(arg) * real(s) * temp2 - w0 * sinh(arg))

    b_jm1 = 1.0 / (4.0 * (w0 * w0))
    b_jm2 = b_jm1

    ! calculate y_1
    mu_t = w1 * b_jm1
    do i = 0,nscl-2
       y_jm2(i) = y_0(i)
       y_jm1(i) = y_0(i) + (mu_t * h * F_0(i))
    enddo

    c_jm2 = 0.0
    c_jm1 = mu_t
    zjm1 = w0
    zjm2 = 1.0
    dzjm1 = 1.0
    dzjm2 = 0.0
    d2zjm1 = 0.0
    d2zjm2 = 0.0

    do j = 2,s

       zj = 2.0 * w0 * zjm1 - zjm2
       dzj = 2.0 * w0 * dzjm1 - dzjm2 + 2.0 * zjm1
       d2zj = 2.0 * w0 * d2zjm1 - d2zjm2 + 4.0 * dzjm1
       b_j = d2zj / (dzj * dzj)
       gamma_t = 1.0 - (zjm1 * b_jm1)

       nu = -b_j / b_jm2
       mu = 2.0 * b_j * w0 / b_jm1
       mu_t = mu * w1 / w0

       ! calculate derivative, use y array for temporary storage
       y_j = dydt(t_rkc + (h * c_jm1), y_jm1, temper, iz)

       do i = 0,nscl-2 
          y_j(i) = (1.0 - mu - nu) * y_0(i) + (mu * y_jm1(i)) + (nu * y_jm2(i)) &
               + h * mu_t * (y_j(i) - (gamma_t * F_0(i)))
       enddo
       
       c_j = (mu * c_jm1) + (nu * c_jm2) + mu_t * (1.0 - gamma_t)

       if(j < s)then
          do i = 0,nscl-2
             y_jm2(i) = y_jm1(i)
             y_jm1(i) = y_j(i)
          enddo
       endif

       c_jm2  = c_jm1
       c_jm1  = c_j
       b_jm2  = b_jm1
       b_jm1  = b_j
       zjm2   = zjm1
       zjm1   = zj
       dzjm2  = dzjm1
       dzjm1  = dzj
       d2zjm2 = d2zjm1
       d2zjm1 = d2zj
    enddo

    do i = 0,nscl-2
       rkc_step(i) = y_j(i)
    enddo

  end function rkc_step

  function dydt(t_rkc, y, temper, iz)
   ! Parameters
    real, intent(in),  dimension(0:nscl-2) :: y
    real, dimension(0:nscl-2) :: dydt, dy
    real, dimension(0:nscl-2) :: c
    real, intent(in) :: t_rkc, temper
    real K1s, K2s, Kw, Kb, Rgas, salt
    integer i
    logical reduced
    real a1, a2, a3, a4, a5, a6, a7
    real b1, b2, b3, b4, b5, b6, b7
    
    reduced = .true.
    salt   = 35.0
    do i = 0,nscl-2
       c(i) = y(i)
    enddo

    K1s = exp(-2307.1266/temper + 2.83655 - 1.5529413*log(temper) + &
         (-4.0484/temper - 0.20760841)*(salt**0.5) + 0.08468345*salt - &
         0.00654208*(salt**1.5) + log(1.0-0.001005*salt))*(1.0e6)
    K2s = exp(-3351.6106/temper - 9.226508 - 0.2005743*log(temper) + &
         (-23.9722/temper - 0.106901773)*(salt**0.5) + 0.1130822*salt - &
         0.00846934*(salt**1.5) + log(1.0-0.001005*salt))*(1.0e6)
    Kw = exp(148.96502 - 13847.26/temper - 23.65218*log(temper) + &
         (118.67/temper - 5.977 + 1.0495*log(temper))*(salt**0.5) - &
         0.01615*salt)*(1.0e6) !(DoE, 1994)
    Kb = exp((-8966.9 - 2890.53*(salt**0.5) - 77.942*salt + &
         1.728*(salt**1.5) - 0.0996*(salt**2))/temper &
         + 148.0248 + 137.1942*(salt**0.5) + 1.62142*salt - &
         (24.4344 + 25.085*(salt**0.5) + 0.2474*salt)*log(temper) + &
         0.053105*(salt**0.5)*temper)*(1.0e6) !(Dickson, 1990)
    Rgas = 0.0083143

    a1 = exp(1246.98-6.19*(10.0**4)/temper - 183.0*log(temper))
    a2 = (4.7e7)*exp(-23.3/(Rgas*temper))/(1.0e6)
    a3 = (5.0e10)/(1.0e6)
    a4 = (6.0e9)/(1.0e6)
    a5 = (1.4e-3)*(1.0e6)
    a6 = (4.58e10)*exp(-(20.8/(Rgas*temper)))/(1.0e6)
    a7 = (3.05e10)*exp(-(20.8/(Rgas*temper)))/(1.0e6)
    b1 = a1/K1s
    b2 = (Kw*a2/K1s)*(1.0e6)
    b3 = a3*K2s
    b4 = (a4*Kw/K2s)*(1.0e6)
    b5 = (a5/Kw)/(1.0e6)
    b6 = (a6*Kw/Kb)*(1.0e6)
    b7 = a7*K2s/Kb
    
   dy(0) = b1*c(1)*c(5)+b2*c(1)-a1*c(0)-a2*c(0)*c(6)
   
   dy(1) = a1*c(0)+a2*c(0)*c(6)-b1*c(1)*c(5)-b2*c(1) &
       +a3*c(2)*c(5)-b3*c(1)-a4*c(1)*c(6)+b4*c(2) &
       +a7*c(2)*c(3)-b7*c(4)*c(1)

   dy(2) = -a3*c(2)*c(5)+ b3*c(1)+a4*c(1)*c(6)-b4*c(2) &
       -a7*c(2)*c(3)+b7*c(4)*c(1)

   dy(3) = -a6*c(3)*c(6)+ b6*c(4)-a7*c(2)*c(3)+b7*c(4)*c(1)

   dy(4) = a6*c(3)*c(6)- b6*c(4)+a7*c(2)*c(3)-b7*c(4)*c(1)

   dy(5) = 0

   dy(6) = b2*c(1)-a2*c(0)*c(6)-a4*c(1)*c(6)+b4*c(2)+a5 &
      -b5*c(5)*c(6)-a6*c(3)*c(6)+b6*c(4)


    do i = 0,nscl-2
       dydt(i) = dy(i)
    enddo

  end function dydt

end module cc