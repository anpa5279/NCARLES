subroutine wrapper(model_type, tracer_matrix, required_cond_matrix, required_cond_array)
    use pars
    use fields
    use fftwk
    use con_data
    use con_stats
    integer, intent(in) :: model_type
    real, dimension(:,:,:,:), intent(inout) :: tracer_matrix
    real, dimension(:,:,:,:), intent(inout) :: required_cond_matrix
    real, dimension(:), intent(inout) :: required_cond_array
    integer :: istage
    integer, intent(in) :: ix, iy, iscl, iz

    !defining the variables in the required_array


    if model_type == 1 then
        ! Call NPZD_Model.dcdt

        use NPZD_Model, only: dcdt
        do iscl=2,nscl
            do iz=izs,ize
                do iy=iys,iye
                    do ix=1,nnx
                        rhs_scl(ix,iy,iz, iscl) = dcdt(ix,iy,iscl,iz)
                    end do  
                end do
            end do
        end do
        
    else if model_type == 2 then
        ! Call cc model

        use cc, only: react_src
        do iscl=2,nscl
            do iz=izs,ize
                do iy=iys,iye
                    do ix=1,nnx
                        rhs_scl(ix,iy,iz, iscl) = react_src(ix,iy,iscl,iz)
                    end do  
                end do
            end do
        end do

    else if model_type == 3 then
        ! call BFM17
    else if model_type == 4 then
        ! call MARBL
    end if
    tracer_matrix = tracer_matrix + dtzeta*rhs_scl


    ! Update required_cond_array if needed
    ! This part depends on how you want to use required_cond_array in your calculations

end subroutine wrapper