SUBROUTINE NPZdebug(Phyto, Zoo, Nutrients, iz, code_location)
    
    USE pars
    USE fields
    USE con_data
    USE con_stats
    
    INCLUDE 'mpif.h'
    INTEGER :: iz
    CHARACTER(len=*) :: code_location

    REAL :: Phyto, Zoo, Nutrients
    CHARACTER(len=100) :: output_filename
    LOGICAL :: file_exists

   

    ! file name
    output_filename = 'debugNPZ.txt'

    !check if file exists
    INQUIRE(file=output_filename, exist=file_exists)

    !write averages to file
    if (.not. file_exists) then
        open(10, file=output_filename, access='append')
        write(10, '(A)') 'Code: ' // trim(code_location)
        write(10, '(A, I10, /)') 'iz: ', iz
        write(10, '(A, F10.6, /)') 'P: ', Phyto
        write(10, '(A, F10.6, /)') 'Z: ', Zoo
        write(10, '(A, F10.6, /)') 'N: ', Nutrients
        write(10, '(A)') '------' // '//'
        close(10)
    else
        open(10, file=output_filename, status='old', access='append')
        write(10, '(A)') 'Code: ' // trim(code_location)
        write(10, '(A, I10, /)') 'iz: ', iz
        write(10, '(A, F10.6, /)') 'P: ', Phyto
        write(10, '(A, F10.6, /)') 'Z: ', Zoo
        write(10, '(A, F10.6, /)') 'N: ', Nutrients
        write(10, '(A)') '------' // '//'
        close(10)
    endif

    RETURN
END SUBROUTINE NPZdebug
