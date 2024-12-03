module NPZD_Model
    implicit none
   
contains
    function dcdt(c, temper, light_intensity)
        ! NPZ from P Franks 1986 recommended by Nikki Lovenduski
        ! Parameters
        real:: vp
        real:: kn = 1.0        ! umolN/l
        real:: rm = 1.0         ! 1/d
        real:: death_rate_zoo = 0.2     ! 1/d
        real:: lambda = 0.2    ! umolN/l
        real:: death_rate_phyto = 0.1   ! 1/d
        real:: alpha = 0.3
        real:: beta = 0.6
        real:: phi = 0.4  ! 1/d
        real:: r_npzd = 0.15   ! 1/d 
        real:: intensity
        real:: a_npz = 0.6  ! 1/d
        real:: b_npz = 1.066
        real:: c_npz = 1.0 !1/Celcius
        real, intent(in):: light_intensity, temper
        real, dimension(1:4), intent(in):: c
        real:: P, Z, N, D
        real, dimension(1:4) :: dcdt

        P=c(1)
        Z=c(2)
        N=c(3)
        D=c(4)
        
        !intensity = rm * P * lambda !Mayzaud-Poulet (1/d)
        intensity = rm !Ivlev (1/d)
        vp=(a_npz*b_npz**(c_npz*temper)) !from Eppley 1972 (1/d)
            
        dcdt(1) = vp * (N / (kn + N)) * light_intensity * P - &
            intensity * (1.0 - exp(-lambda * P)) * Z - death_rate_phyto &
            * P- r_npzd * P

        dcdt(2) = beta * intensity * (1.0 - exp(-lambda * P)) * Z - death_rate_zoo * Z

        dcdt(3) = -vp * (N / (kn + N)) * light_intensity * P + alpha * &
            intensity * (1.0 - exp(-lambda * P)) * Z + death_rate_phyto * P + &
            death_rate_zoo * Z+ phi*D

        dcdt(4) = r_npzd * P + (1 - alpha - beta) * intensity * (1.0 - exp(-lambda * P)) * Z &
            - phi * D
        
        !print*, dcdt
    end function dcdt
end module NPZD_Model