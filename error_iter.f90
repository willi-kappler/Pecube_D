module m_error_iter
    use m_pecube_func
    use m_pecube_config
    use m_compiler
    use m_logger

    implicit none

! Written by Willi Kappler
! willi.kappler@uni-tuebingen.de
    type observable_t
        real(8) :: age_observed
        real(8) :: age_predicted
        real(8) :: misfit
        real(8) :: px
        real(8) :: py
        real(8) :: pz
        real(8) :: z_offset
        integer(4) :: x
        integer(4) :: y
    end type

    ! Module global variables and constants
    integer(4), parameter :: file_unit = 87
    integer(4) :: io_status, num_of_observables

    type(observable_t), dimension(:), allocatable :: observable

    real(8), parameter :: FROM_KM_TO_M = 1000.0
    real(8), dimension(:, :), allocatable :: topo_data

    ! make everything private
    private

    public error_iter

    contains

    subroutine error_iter(config)
        implicit none

        type(config_t) :: config

        ! local variables
        integer(4) :: current_iteration = 1, iter_counter = 1

        ! Only one age (AHe) is allowed in the Pecube.in file.
        ! Run pecube once to get some data
        call start_pecube(config)

        do
            call log_message("Current iteration: " + current_iteration)
            call load_topography(topo_file_name(1))
            call load_observables(config)

            if (check_misfit(config%error_iter_misfit_limit)) then
                call log_message("All points are below the misfit, exit now")
                exit
            endif

            call log_message("Clean up previous run")
            call sys_command("rm -f " // trim(config%output_folder) // "/*")

            call mutate_topography(config)
            call save_topography(topo_file_name(1), config)

            call start_pecube(config)
            call logger_flush()

            current_iteration = current_iteration + 1
            iter_counter = iter_counter + 1

            if (iter_counter >= 100) then
                iter_counter = 1
                ! Change time from Pecube.in file
                config%error_iter_age_dec = config%error_iter_age_dec + 1.0
                call log_message("Change time in Pecube: " + config%error_iter_age_dec)
            endif
        end do

        call write_results("error_iter_results.txt")

        call log_message("iteration finished!")
    end subroutine error_iter

    subroutine start_pecube(config)
        implicit none

        type(config_t) :: config

        call log_message("Starting Pecube")
        call pecube_func(config)

        call sys_command("cp " // trim(config%output_folder) // "/Comparison_ages_tec_001.txt observables.txt")
        call sys_command("cp " // topo_file_name(1) // " best_topo.dat")
    end subroutine start_pecube

    subroutine load_topography(topo_filename)
        implicit none

        character(*), intent(in) :: topo_filename
        integer(4) :: x, y

        call log_message("Open topography file: " + trim(topo_filename))

        open(file_unit, file=trim(topo_filename), status="old", action="read", iostat=io_status)

        if (io_status /= 0) then
            call log_message("Error while reading file: ")
            call log_message(trim(topo_filename))
            call log_message("Error code: d1ae86a1005601273ac71ad5f377ac17")

            call clean_up()
        endif

        if (.not. allocated(topo_data)) then
            call log_message("Allocate topo_data")
            allocate(topo_data(nx0, ny0))
        endif

        do y = 1, ny0
            do x = 1, nx0
                read(file_unit, *, iostat=io_status) topo_data(x, y)

                if (io_status /= 0) then
                    call log_message("Error while reading file: ")
                    call log_message(trim(topo_filename))
                    call log_message("Error code: fbb495093af90b17158acc303a05c616")
                    call log_message("x: " +  x + ", y: " + y)
                    call log_message("Value: " + topo_data(x, y))

                    call clean_up()
                endif
            end do
        end do

        call log_message("New topography loaded, min: " + minval(topo_data) + ", max: " + maxval(topo_data))

        close(file_unit)
    end subroutine load_topography

    subroutine load_observables(config)
        implicit none

        type(config_t), intent(in) :: config

        character(255) :: line, observables_file
        integer(4) :: x, y, i, j, skip_int
        real(8) :: px, py, pz, skip_float, age_predicted, age_error, age_observed
        real(8) :: diffcoords, diffspacing, age_diff, erosion_rate

        observables_file = trim(config%output_folder) // "/Comparison_ages_tec_001.txt"

        open(file_unit, file=trim(observables_file), status="old", action="read", iostat=io_status)

        if (io_status /= 0) then
            call log_message("Error while reading file: ")
            call log_message(trim(observables_file))
            call log_message("Error code: e36edf489d307edce020967399722409")

            call clean_up()
        endif

        read(file_unit, "(a)", iostat=io_status) line

        if (io_status /= 0) then
            call log_message("Error while reading file: ")
            call log_message(trim(observables_file))
            call log_message("Error code: 75cd2517cb846e7d5382d22b2cd8e777")

            call clean_up()
        endif

        ! Input looks like this:
        ! Number of observations:          120
        read(line(25:), *, iostat=io_status) num_of_observables

        if (io_status /= 0) then
            call log_message("Error while reading file: ")
            call log_message(trim(observables_file))
            call log_message("Line: " + trim(line(25:)))
            call log_message("Error code: 45317fd210f515b4fdaf16ad276dcdd3")

            call clean_up()
        endif

        call log_message("Number of observables: " + num_of_observables)

        ! Skip header line
        read(file_unit, *, iostat=io_status)

        if (io_status /= 0) then
            call log_message("Error while reading file: ")
            call log_message(trim(observables_file))
            call log_message("Error code: 0832652050140a26f3c561053cf8c477")

            call clean_up()
        endif

        call log_message("xlonmin: " + xlonmin + ", xlatmin: " + xlatmin)
        call log_message("spacing_long: " + config%spacing_long + ", spacing_lat: " + config%spacing_lat)
        call log_message("nx0: " + nx0 + ", ny0: " + ny0)
        call log_message("nsurf: " + nsurf + ", nx0 * ny0: " + nx0 * ny0)

        ! Each line looks like this:
        ! Data Set, Longitude, Latitude, Height Obs., Height Int., AHe Type, AHe Age Obs., AHe Age Error, AHe Age Prd.

        if (.not. allocated(observable)) then
            allocate(observable(num_of_observables))
        endif

        diffspacing = hypot(config%spacing_long, config%spacing_lat)

        do i = 1, num_of_observables
            read(file_unit, *, iostat=io_status) skip_int, px, py, pz, skip_float, skip_int, age_observed, &
                age_error, age_predicted

            if (io_status /= 0) then
                call log_message("Error while reading file: ")
                call log_message(trim(config%output_folder) // "/Comparison_ages_tec_001.txt")
                call log_message("At index: " + i)
                call log_message("Error code: c1cabd2ee7c857ecba9bda3e9d490623")

                call clean_up()
            endif

            x = nint((px - xlonmin) / config%spacing_long)
            x = min(max(x, 1), nx0)
            y = nint((py - xlatmin) / config%spacing_lat)
            y = min(max(y, 1), ny0)

            do j = 1, num_of_observables
                if (i /= j) then
                    ! Multiply by 0.5 to have a buffer and ensure that the distance is really big enough
                    diffcoords = hypot(observable(j)%px - px, observable(j)%py - py) * 0.5

                    if (diffcoords < diffspacing) then
                        call log_message("Observables with the same coordinates detected!")
                        call log_message("diffcoords: " + diffcoords + ", diffspacing: " + diffspacing)
                        call log_message("1. index: " + i + ", x: " + x + ", y: " + y + ", px: " + px + ", py: " + py)
                        call log_message("2. index: " + j + ", x: " + observable(j)%x + ", y: " + observable(j)%y + ", px: " + &
                            observable(j)%px + ", py: " + observable(j)%py)
                        ! call clean_up()
                    endif
                endif
            enddo

            observable(i)%x = x
            observable(i)%y = y
            observable(i)%px = px
            observable(i)%py = py
            observable(i)%pz = pz
            observable(i)%age_predicted = age_predicted
            observable(i)%age_observed = age_observed
            age_diff = age_predicted - age_observed
            observable(i)%misfit = abs(age_diff) / age_error
            erosion_rate = last_peclet
            observable(i)%z_offset = age_diff * erosion_rate * FROM_KM_TO_M

            call log_message("Observable indices: x: " + x + ", y: " + y + ", Observable coordinates: px: " + px + &
                ", py: " + py + ", pz: " + pz)
            call log_message("age_predicted: " + age_predicted + ", age_observed: " + age_observed + ", age_error: " + age_error)
            call log_message("Index: " + i + ", z_offset: " + observable(i)%z_offset + ", misfit: " + observable(i)%misfit + &
                ", z: " + topo_data(x, y) + ", erosion rate: " + erosion_rate)
        end do
        close(file_unit)
    end subroutine load_observables

    subroutine mutate_topography(config)
        implicit none

        type(config_t), intent(in) :: config

        integer(4) :: i

        do i = 1, num_of_observables
            if (observable(i)%misfit > config%error_iter_misfit_limit) then
                call log_message("Observable " + i + ": misfit > tolerance: " + observable(i)%misfit + " > " + &
                    config%error_iter_misfit_limit + ", z_offset: " + observable(i)%z_offset)

                if (observable(i)%age_observed < observable(i)%age_predicted) then
                    call log_message("Modify topography, age_observed < age_predicted: " + observable(i)%age_observed + &
                        " < "  + observable(i)%age_predicted)

                    call mutate_point(observable(i)%x, observable(i)%y, config%error_iter_radius, &
                        observable(i)%z_offset)
                else
                    call log_message("No changes in topography, age_observed >= age_predicted: " + observable(i)%age_observed + &
                        " >= "  + observable(i)%age_predicted)
                endif
            endif
        enddo
    end subroutine mutate_topography

    subroutine mutate_point(x, y, radius, offset)
        implicit none

        integer(4), intent(in) :: x, y, radius
        real(8), intent(in) :: offset

        integer(4) :: xr, yr
        real(8) :: rr, diff

        !call log_message("x: " + x + ", y: " + y + ", radius: " + radius + ", pz: " + pz + ", offset: " + offset)

        if (radius == 0) then
            call log_message("topo before: " + topo_data(x, y))
            topo_data(x, y) = topo_data(x, y) + offset
            call log_message("topo after: " + topo_data(x, y))
        else
            do yr = y - radius, y + radius
                do xr = x - radius, x + radius
                    if (yr >= 1 .and. yr <= ny0) then
                        if (xr >= 1 .and. xr <= nx0) then
                            diff = hypot(dble(x - xr), dble(y - yr))
                            rr = exp(- diff / radius)
                            call log_message("xr: " + xr + ", yr: " + yr + ", diff: " + diff + ", rr: " + rr)
                            call log_message("topo before: " + topo_data(xr, yr))
                            topo_data(xr, yr) = topo_data(xr, yr) + (offset * rr)
                            call log_message("topo after: " + topo_data(xr, yr))
                        endif
                    endif
                enddo
            enddo
        endif

    end subroutine mutate_point

    subroutine save_topography(topo_filename, config)
        implicit none

        integer(4) :: x, y
        type(config_t), intent(in) :: config
        character(*), intent(in) :: topo_filename

        call log_message("Saving new topography data: " + trim(topo_filename))

        open(file_unit, file=trim(topo_filename), status="unknown", iostat=io_status)

        if (io_status /= 0) then
            call log_message("Error while opening file: ")
            call log_message(trim(topo_filename))
            call log_message("Error code: 858df9b09eede32f58c1044322870101")

            call clean_up()
        endif

        do y = 1, ny0
            do x = 1, nx0
                if (topo_data(x, y) < 0.0) then
                    ! If topography gets negative reset it to zero
                    call log_message("topo too low: " + topo_data(x, y) + " (x: " + x + ", y: " + y + ") will be set to 0")
                    topo_data(x, y) = 0.0
                else if (topo_data(x, y) > config%error_iter_max_topography) then
                    ! If topography gets too high reset it to the max allowed value
                    call log_message("topo too high: " + topo_data(x, y) + " (x: " + x + ", y: " + y &
                        + ") will be set to " + config%error_iter_max_topography)
                    topo_data(x, y) = config%error_iter_max_topography
                endif

                write(file_unit, "(F0.4)", iostat=io_status) topo_data(x, y)

                if (io_status /= 0) then
                    call log_message("Error while writing to file: ")
                    call log_message(trim(topo_filename))
                    call log_message("Error code: 4bcea079af481275f30c32cf7e778249")

                    call clean_up()
                endif
            end do
        end do
        close(file_unit)

        call log_message("Topo min: " + minval(topo_data))
        call log_message("Topo max: " + maxval(topo_data))
    end subroutine save_topography

    function check_misfit(tolerance)
        implicit none

        real(8), intent(in) :: tolerance
        logical :: check_misfit
        integer(4) :: i

        check_misfit = .true.

        do i = 1, num_of_observables
            if (observable(i)%age_observed < observable(i)%age_predicted) then
                if (observable(i)%misfit > tolerance) then
                    check_misfit = .false.
                    return
                endif
            endif
        enddo
    end function check_misfit

    subroutine write_results(filename)
        implicit none

        character(*), intent(in) :: filename
        integer(4) :: i

        call log_message("Saving results: " + trim(filename))

        open(file_unit, file=trim(filename), status="unknown", iostat=io_status)

        if (io_status /= 0) then
            call log_message("Error while opening file: ")
            call log_message(trim(filename))
            call log_message("Error code: 858df9b09eede32f58c1044322870101")

            call clean_up()
        endif

        write(file_unit, "(A)", iostat=io_status, advance="no") "#index, "
        write(file_unit, "(A)", iostat=io_status, advance="no") "age observed, "
        write(file_unit, "(A)", iostat=io_status, advance="no") "age predicted, "
        write(file_unit, "(A)", iostat=io_status, advance="no") "misfit, "
        write(file_unit, "(A)", iostat=io_status, advance="no") "x, "
        write(file_unit, "(A)", iostat=io_status, advance="no") "y, "
        write(file_unit, "(A)", iostat=io_status) "z"

        if (io_status /= 0) then
            call log_message("Error while writing to file: ")
            call log_message(trim(filename))
            call log_message("Error code: 4bcea079af481275f30c32cf7e778249")

            call clean_up()
        endif


        do i = 1, num_of_observables
            write(file_unit, "(I3, A)", iostat=io_status, advance="no") i, ", "
            write(file_unit, "(F0.4, A)", iostat=io_status, advance="no") observable(i)%age_observed, ", "
            write(file_unit, "(F0.4, A)", iostat=io_status, advance="no") observable(i)%age_predicted, ", "
            write(file_unit, "(F0.4, A)", iostat=io_status, advance="no") observable(i)%misfit, ", "
            write(file_unit, "(F0.4, A)", iostat=io_status, advance="no") observable(i)%px, ", "
            write(file_unit, "(F0.4, A)", iostat=io_status, advance="no") observable(i)%py, ", "
            write(file_unit, "(F0.4, A)", iostat=io_status) observable(i)%pz

            if (io_status /= 0) then
                call log_message("Error while writing to file: ")
                call log_message(trim(filename))
                call log_message("Error code: 4bcea079af481275f30c32cf7e778249")

                call clean_up()
            endif
        enddo

        close(file_unit)
    end subroutine write_results

    subroutine clean_up()
        implicit none

        logical :: file_is_opened

        if (allocated(topo_data)) then
            deallocate(topo_data)
        endif

        if (allocated(observable)) then
            deallocate(observable)
        endif

        inquire(file_unit, opened=file_is_opened)

        if (file_is_opened) then
            close(file_unit, iostat=io_status)

            if (io_status /= 0) then
                call log_message("Cloud not close file, maybe already closed")
            endif
        endif

        call logger_flush()

        stop
    end subroutine clean_up
end module m_error_iter
