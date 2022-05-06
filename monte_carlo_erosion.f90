! Written by Willi Kappler
! willi.kappler@uni-tuebingen.de

! TODO: fix bug: monte carlo only sets erosion rates for time steps 0 to n-1.
! It needs to set erosion rates for all n time steps

module m_monte_carlo_erosion
    use m_pecube_func
    use m_compiler
    use m_pecube_config
    use m_logger

    use mpi

    implicit none

    private

    ! module global constants
    integer(4), parameter :: file_unit = 101
    integer(4), parameter :: file_unit2 = 102
    integer(4), parameter :: file_unit3 = 103
    integer(4), parameter :: max_line_length = 255
    integer(4), parameter :: chronometer_name_length = 10

    character(1), parameter :: csv_delimiter = ";"


    ! type definition
    type :: tuple_age_sigma_t
        real(8) :: age
        real(8) :: sigma
    endtype

    type :: chronometer_data_t
        integer(4) :: id
        character(chronometer_name_length) :: short_name
        type(tuple_age_sigma_t), dimension(:), allocatable :: all_tuples
    endtype

    type :: good_fit_stat_t
        real(8) :: chisq_gfit
        real(8) :: ages_gfit
    endtype

    ! module global variable
    integer(4) :: io_status, num_of_data_sets_in_line, num_of_chronometers
    integer(4) :: mc_mpi_error
    ! we need 64 bit interger here:
    integer(8), dimension(:), allocatable :: number_of_good_fits
    integer(4), dimension(MPI_STATUS_SIZE) :: mc_mpi_status

    character(max_line_length) :: global_line

    real(8), dimension(:,:), allocatable :: chi_squared
    real(8), dimension(:,:), allocatable :: predicted_ages
    real(8), dimension(:), allocatable :: random_difs

    type(chronometer_data_t), dimension(:), allocatable :: all_chronometers

    public monte_carlo_erosion

    contains

    subroutine monte_carlo_erosion(config)
        implicit none

        type(config_t) :: config

        ! to compile with OpenMPI use:
        ! scons --use-mpi

        integer(4) :: seed_size, i, j, k, column1, column2, column3, num_of_entries
        integer(4) :: mc_mpi_error, return_code, cpu_id_with_min_chi_squared
        integer(4), dimension(:), allocatable :: seed

        logical :: has_good_fit, other_has_good_fit

        character(4) :: suffix
        character(64) :: folder_name

        real(8) :: chi_squared_min_sum, other_chi_squared_min_sum
        real(8) :: random_index, column4, column5, c_time, random_coin
        real(8) :: new_erosion_value
        real(8), dimension(:), allocatable :: random_erosion_list

        if (config%mc_num_of_simulations < 2) then
          call log_message("You must use at least 2 iterations")
          error stop 1
        endif

        if (config%mpi_total_num_of_cpu < 4) then
          call log_message("You must use at least 4 cpus")
          error stop 1
        endif

        ! init random seed
        call random_seed(size = seed_size)
        call log_message("seed_size: " + seed_size)
        allocate(seed(seed_size))

        open (file_unit, file="/dev/urandom", access="stream", &
             form="unformatted", action="read", status="old", iostat=io_status)

        if (io_status == 0) then
           read (file_unit) seed
           close (file_unit)
        else
            call log_message("warning, could not open /dev/urandom!")

            call cpu_time(c_time)

            call log_message("c_time: " + c_time)

            do i = 1, seed_size
                seed(i) = int(c_time * 100.0) + i
            enddo
        endif

        call log_message("seed: " + seed)

        call random_seed(put=seed)

        deallocate(seed)

        num_of_entries = int((config%mc_max_erosion_rate - config%mc_min_erosion_rate) / config%mc_erosion_step) - 1
        call log_message("num_of_entries: " + num_of_entries)
        allocate(random_erosion_list(num_of_entries + 1))

        do i = 0, num_of_entries
            random_erosion_list(i + 1) = config%mc_min_erosion_rate + (dble(i) * config%mc_erosion_step)
            call log_message("random_erosion_list(" + (i + 1) + "): " + random_erosion_list(i + 1))
        enddo

        call read_csv_input_file(config)

        write(suffix, "(I4.4)") config%mpi_current_cpu_id
        folder_name = "mpi_mc_process_" // suffix
        call log_message("create directory: '" // trim(folder_name) // "'")

        call sys_command("mkdir -p " // trim(folder_name))
        call sys_command("mkdir -p " // trim(folder_name) // "/input/Pecube-D")
        call sys_command("mkdir -p " // trim(folder_name) // "/output/Pecube-D")
        call sys_command("mkdir -p " // trim(folder_name) // "/output/mc_results")
        call sys_command("cp Pecube.in " // trim(folder_name))

        ! call sys_command("pwd >> " // config%log_filename)
        call sys_chdir(folder_name)
        ! call sys_command("pwd >> ../" // config%log_filename)

        allocate(chi_squared(num_of_chronometers, num_of_data_sets_in_line))
        allocate(predicted_ages(num_of_chronometers, num_of_data_sets_in_line))
        allocate(number_of_good_fits(num_of_data_sets_in_line))

        ! write header for the two output files
        open(file_unit, file="output/mc_results/erate_gfit.txt", status="unknown", iostat=io_status)
        if (io_status /= 0) then
            call log_message("error while opening file: output/mc_results/erate_gfit.txt")
            call log_message("error code: 29cda9b12b7e6e65395d425be3017510")
            call clean_up_stop()
        endif
        write(file_unit, *) "# number_of_good_fits(j), k, j, mc_random_erosion_rate(k)"
        close(file_unit)

        open(file_unit, file="output/mc_results/good_fit_stats.txt", status="unknown", iostat=io_status)
        if (io_status /= 0) then
            call log_message("error while opening file: output/mc_results/good_fit_stats.txt")
            call log_message("error code: 5f61996e93741163395bd03783d8113b")
            call clean_up_stop()
        endif
        write(file_unit, *) "# number_of_good_fits(j), k, j, chi_squared(k, j), predicted_ages(k, j)"
        close(file_unit)

        ! init variable with 0
        number_of_good_fits = 0

        chi_squared_min_sum = huge(chi_squared_min_sum)

        ! main loop
        do i = 1, config%mc_num_of_simulations, config%mpi_total_num_of_cpu
            call log_message("monte carlo iteration number: " + i)

            !chi_squared_min_sum = huge(chi_squared_min_sum)

            if (allocated(mc_random_erosion_rate)) then
                call log_message("monte_carlo.f90: calling pecube_func with these erosion rates: " + mc_random_erosion_rate)
            endif

            ! Clean up otuput files from previous run
            call sys_command("rm -f output/*.dat")
            call pecube_func(config)
            call log_message("monte_carlo.f90: pecube has calculated these ages: " + surface_age_info(nstep)%all_ages(1, :))

            ! calculate reduced chi-squared misfit
            do j = 1, num_of_data_sets_in_line ! u in matlab

                do k = 1, num_of_chronometers ! w in matlab
                    predicted_ages(k, j) = surface_age_info(nstep)%all_ages(1, all_chronometers(k)%id)

                    chi_squared(k, j) = ((all_chronometers(k)%all_tuples(j)%age - predicted_ages(k, j))**2) / &
                       ((all_chronometers(k)%all_tuples(j)%sigma)**2)
                enddo

                call log_message("chi_squared: (j: " + j + "), " + chi_squared(:, j))
                call log_message("predicted_ages: " + predicted_ages(:, j))

                call log_message("current number_of_good_fits(j): " + number_of_good_fits(j) + ", j: " + j)

                if (sum(chi_squared) < chi_squared_min_sum) then
                    call log_message("new minimum found: " + sum(chi_squared) + " < " + chi_squared_min_sum)
                    chi_squared_min_sum = sum(chi_squared)
                endif

                has_good_fit = .false.

                ! check if all elements in row are smaller than the threshold
                if (all(chi_squared(:, j) <= config%mc_tolerance_chi_squared)) then
                    has_good_fit = .true.
                endif

                if (config%mpi_current_cpu_id == 0) then
                    ! ask all processes if there was a good fit,
                    ! if yes increase global number of good fits

                    do k = 1, config%mpi_total_num_of_cpu - 1
                        call MPI_RECV(other_has_good_fit, 1, MPI_LOGICAL, MPI_ANY_SOURCE,&
                            MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                        if (mc_mpi_error /= MPI_SUCCESS) then
                            call log_message("mpi error: a47641b3b44cae25df50d564923cd367")
                            call clean_up_stop()
                        endif
                        if (other_has_good_fit) then
                            ! increase global number of good fits and send it to the client process
                            number_of_good_fits(j) = number_of_good_fits(j) + 1
                            call log_message("other process has good fit: " + number_of_good_fits)
                            call MPI_SEND(number_of_good_fits(j), 1, MPI_INTEGER8, mc_mpi_status(MPI_SOURCE), 0, &
                                MPI_COMM_WORLD, mc_mpi_error)
                            if (mc_mpi_error /= MPI_SUCCESS) then
                                call log_message("mpi error: 9edcf5afbc8afaa60b51c1cd588c6cf8")
                                call clean_up_stop()
                            endif
                        else
                            ! otherwise just send current value from last run
                            call MPI_SEND(number_of_good_fits(j), 1, MPI_INTEGER8, mc_mpi_status(MPI_SOURCE), 0, &
                                MPI_COMM_WORLD, mc_mpi_error)
                            if (mc_mpi_error /= MPI_SUCCESS) then
                                call log_message("mpi error: f14f8f8b1c404f7934c914c8b8d8db53")
                                call clean_up_stop()
                            endif
                        endif
                    enddo

                    if (has_good_fit) then
                        ! increase global counter for master after increasing global counter for clients
                        number_of_good_fits(j) = number_of_good_fits(j) + 1 ! cnt2 in matlab
                    endif

                else
                    ! send if we have a good fit to master process (id = 0)
                    call MPI_SEND(has_good_fit, 1, MPI_LOGICAL, 0, 0, MPI_COMM_WORLD, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: f14f8f8b1c404f7934c914c8b8d8db53")
                        call clean_up_stop()
                    endif
                    ! receive the global number of good fits from master
                    call MPI_RECV(number_of_good_fits(j), 1, MPI_INTEGER8, MPI_ANY_SOURCE, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: 3a2a15f83023f6475387fed3b6c09e9e")
                        call clean_up_stop()
                    endif
                endif


                if (has_good_fit) then
                    ! keep track of how many times a good fit is found
                    call log_message("number_of_good_fits: " + number_of_good_fits)

                    open(file_unit, file="output/mc_results/erate_gfit.txt", status="unknown", &
                        position="append", iostat=io_status)
                    if (io_status /= 0) then
                        call log_message("error while opening file: output/mc_results/erate_gfit.txt")
                        call log_message("error code: 2e04c87035fafe5d3d2e79397a290af1")
                        call clean_up_stop()
                    endif

                    do k = 1, mc_num_of_erosion_rates
                        ! the erate_gfit array stores the erosion history for each
                        ! simulation that fell below the user defined tolerance
                        ! write the erosion history for all good fit histories
                        write(file_unit, *) number_of_good_fits(j), k, j, mc_random_erosion_rate(k)
                    enddo

                    close(file_unit)

                    open(file_unit, file="output/mc_results/good_fit_stats.txt", status="unknown", &
                        position="append", iostat=io_status)
                    if (io_status /= 0) then
                        call log_message("error while opening file: output/mc_results/good_fit_stats.txt")
                        call log_message("error code: 343035d10920fff1cab149e62b109888")
                        call clean_up_stop()
                    endif

                    do k = 1, num_of_chronometers
                        ! the chisq_gfit matrix stores the chisq for each model
                        ! that fell below the tolerance
                        ! the ages_gfit matrix stores the ages for each model that
                        ! fell below the tolerance.
                        write(file_unit, *) number_of_good_fits(j), k, j, chi_squared(k, j), predicted_ages(k, j)
                    enddo

                    close(file_unit)
                else
                    open(file_unit, file="output/mc_results/erate_nofit.txt", status="unknown", position="append", iostat=io_status)
                    if (io_status /= 0) then
                        call log_message("error while opening file: output/mc_results/erate_nofit.txt")
                        call log_message("error code: 49d8238451790ff1ff69a9e0b5c35094")
                        call clean_up_stop()
                    endif

                    do k = 1, mc_num_of_erosion_rates
                        ! Write erosion history for cases with no good fit.
                        write(file_unit, *) k, j, mc_random_erosion_rate(k)
                    enddo
                endif

                ! wait for all processes to reach this point (barrier)
                call MPI_Barrier(MPI_COMM_WORLD, mc_mpi_error)
            enddo ! j = 1, num_of_data_sets_in_line

            ! begin of evolutionary / genetic algorithm
            if (config%mpi_current_cpu_id == 0) then
                cpu_id_with_min_chi_squared = 0

                do j = 1, config%mpi_total_num_of_cpu - 1
                    ! receive values from all other processes
                    call log_message("receive other_chi_squared_min_sum... " + j)
                    call MPI_RECV(other_chi_squared_min_sum, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: d6e4fba92219be2207c3a190c9bd7141")
                        call clean_up_stop()
                    endif

                    call log_message("other_chi_squared_min_sum: " + other_chi_squared_min_sum)
                    if (other_chi_squared_min_sum < chi_squared_min_sum) then
                        call log_message("new minimum found - old value: " + chi_squared_min_sum)
                        ! we found one that is better than anyone else
                        chi_squared_min_sum = other_chi_squared_min_sum
                        cpu_id_with_min_chi_squared = mc_mpi_status(MPI_SOURCE)
                        call log_message("cpu id: " + cpu_id_with_min_chi_squared)
                    endif
                enddo

                do j = 1, config%mpi_total_num_of_cpu - 1
                    ! send the cpu id to all the other processes
                    call log_message("send cpu_id_with_min_chi_squared to other pcoresses... " + j)
                    call MPI_SEND(cpu_id_with_min_chi_squared, 1, MPI_INTEGER, j, 0, &
                        MPI_COMM_WORLD, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: ffbfdbbdd31950195d307861e801ad3f")
                        call clean_up_stop()
                    endif
                enddo

                if (cpu_id_with_min_chi_squared == 0) then
                    ! no other process is better than the master process,
                    ! so we just send what the master has to the others
                    call log_message("no better configuration found. Send master configuration to others")
                    call log_message("master chi_squared_min_sum: " + chi_squared_min_sum)
                    do j = 1, config%mpi_total_num_of_cpu - 1
                        call log_message("sending master config... " + j)
                        call MPI_SEND(mc_random_erosion_rate, mc_num_of_erosion_rates, MPI_DOUBLE_PRECISION, j, 0, &
                            MPI_COMM_WORLD, mc_mpi_error)
                        if (mc_mpi_error /= MPI_SUCCESS) then
                            call log_message("mpi error: 198b165860b41b29933c250ffe68a8dc")
                            call clean_up_stop()
                        endif
                    enddo
                else
                    ! first we receive the new configuration
                    call log_message("receive new config from process: " + cpu_id_with_min_chi_squared)
                    call MPI_RECV(mc_random_erosion_rate, mc_num_of_erosion_rates, MPI_DOUBLE_PRECISION, &
                        cpu_id_with_min_chi_squared, MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: 96d04fe0183639337989a705a392cd45")
                        call clean_up_stop()
                    endif

                    call log_message("new optimal config: " + mc_random_erosion_rate)

                    ! then send this informatioin to the other processes
                    do j = 1, config%mpi_total_num_of_cpu - 1
                        call log_message("send new config to other processes... " + j)
                        call MPI_SEND(mc_random_erosion_rate, mc_num_of_erosion_rates, MPI_DOUBLE_PRECISION, j, 0, &
                            MPI_COMM_WORLD, mc_mpi_error)
                        if (mc_mpi_error /= MPI_SUCCESS) then
                            call log_message("mpi error: 9bb9768e78bbfd1d836df51130368d6d")
                            call clean_up_stop()
                        endif
                    enddo
                endif

                ! set new 100% random configuration for master node
                call log_message("setting random erosion rate values (monte carlo) for master")
                do j = 1, mc_num_of_erosion_rates

                    call random_number(random_index)
                    random_index = random_index * dble(num_of_entries)
                    mc_random_erosion_rate(j) = random_erosion_list(int(random_index) + 1)

                    if (j <= mc_max_number_of_erosion_rates) then
                        if (config%mc_fixed_erosion_rates(j) >= 0.0) then
                            mc_random_erosion_rate(j) = config%mc_fixed_erosion_rates(j)
                        endif
                    endif
                enddo
                call log_message("random erosion rate values (master): " + mc_random_erosion_rate)

            else
                ! send the value for the minimum to the master process
                call log_message("sending minimum value to the master: " + chi_squared_min_sum)
                call MPI_SEND(chi_squared_min_sum, 1, MPI_DOUBLE_PRECISION, 0, config%mpi_current_cpu_id, &
                    MPI_COMM_WORLD, mc_mpi_error)
                if (mc_mpi_error /= MPI_SUCCESS) then
                    call log_message("mpi error: ce2bd032ff192f4d205f7cbb0864c76b")
                    call clean_up_stop()
                endif

                ! receive the cpu id from the master process
                call log_message("receiving cpu id from master...")
                call MPI_RECV(cpu_id_with_min_chi_squared, 1, MPI_INTEGER, 0, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                if (mc_mpi_error /= MPI_SUCCESS) then
                    call log_message("mpi error: 29c27565d644cf3a87686e0f62af254c")
                    call clean_up_stop()
                endif
                call log_message("cpu id from master: " + cpu_id_with_min_chi_squared)

                if (cpu_id_with_min_chi_squared == 0) then
                    ! no other process is better than the master process
                    ! so we just receive what the master process has
                    call log_message("no other process was better than the master - receiving master configuration...")
                    call MPI_RECV(mc_random_erosion_rate, mc_num_of_erosion_rates, MPI_DOUBLE_PRECISION, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: 13f88a90f30b6fc4d31296216b5d9d9f")
                        call clean_up_stop()
                    endif
                else
                    if (cpu_id_with_min_chi_squared == config%mpi_current_cpu_id) then
                        ! this process is the one with the best configuration
                        ! so we send it to the master process
                        call log_message("we are the process with the best config: " + cpu_id_with_min_chi_squared + &
                            ", " + chi_squared_min_sum)
                        call log_message("sending out config to master...")
                        call MPI_SEND(mc_random_erosion_rate, mc_num_of_erosion_rates, MPI_DOUBLE_PRECISION, 0, 0, &
                            MPI_COMM_WORLD, mc_mpi_error)
                        if (mc_mpi_error /= MPI_SUCCESS) then
                            call log_message("mpi error: 7a667cf23c4c4cb24656745adbaa94fe")
                            call clean_up_stop()
                        endif
                    endif
                    ! receive new configuration from master process
                    call log_message("receiving new configuration from master")
                    call MPI_RECV(mc_random_erosion_rate, mc_num_of_erosion_rates, MPI_DOUBLE_PRECISION, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, mc_mpi_status, mc_mpi_error)
                    if (mc_mpi_error /= MPI_SUCCESS) then
                        call log_message("mpi error: ce82302d736d1b4dbbd3c4f96c40ad46")
                        call clean_up_stop()
                    endif
                    call log_message("new configuration: " + mc_random_erosion_rate)
                endif

                ! set new 100% random configuration for other nodes
                call log_message("setting random erosion rate values (monte carlo) for other nodes")
                do j = 1, mc_num_of_erosion_rates
                    !call random_number(random_index)
                    !random_index = random_index * dble(num_of_entries)
                    !mc_random_erosion_rate(j) = random_erosion_list(int(random_index) + 1)

                    new_erosion_value = mc_random_erosion_rate(j)

                    call random_number(random_coin)

                    if (random_coin < 0.33) then
                        new_erosion_value = max(new_erosion_value - config%mc_erosion_step, config%mc_min_erosion_rate)
                    else if (random_coin > 0.66) then
                        new_erosion_value = min(new_erosion_value + config%mc_erosion_step, config%mc_max_erosion_rate)
                    endif

                    if (j <= mc_max_number_of_erosion_rates) then
                        if (config%mc_fixed_erosion_rates(j) >= 0.0) then
                            new_erosion_value = config%mc_fixed_erosion_rates(j)
                        endif
                    endif

                    mc_random_erosion_rate(j) = new_erosion_value

                enddo
                call log_message("random erosion rate values (not master): " + mc_random_erosion_rate)
            endif
        enddo ! main loop, i = 1, config%mc_num_of_simulations

        call log_message("monte carlo iteration finished")

        ! wait for all processes to reach this point (barrier)
        call MPI_Barrier(MPI_COMM_WORLD, mc_mpi_error)

        if (mc_mpi_error /= MPI_SUCCESS) then
            call log_message("error: could not wait for mpi barrier")
            call log_message("error code: 53827676fbe423905d2293dc7e97a90e")
            return_code = mc_mpi_error
            call MPI_ABORT(MPI_COMM_WORLD, return_code, mc_mpi_error)
        endif

        if (config%mpi_current_cpu_id == 0) then
            ! Collect all the information
            ! an merge it into one file.
            ! This must be done by the master process (id 0)
            call sys_chdir("..")
            ! call sys_command("pwd >> " // config%log_filename)



            ! Matlab needs these files and datalayout:

            ! max_erate: mc_maximum_erosion_rate
            ! time_steps: time step years
            ! ages_gfit(w,cnt2(u),u) = predicted ages
            ! chisq_gfit(w,cnt2(u),u) = chi squared
            ! erate_gfit(ct,cnt2(u),u = monte carlo estimated erosion rates

            call log_message("merging all output files into one for matlab")

            open(file_unit2, file="erate_gfit.txt", status="unknown", iostat=io_status)
            if (io_status /= 0) then
                call log_message("error while opening file: erate_gfit.txt")
                call log_message("error code: a101227ca55278e5605bccf1a4ac785a")
                call clean_up_stop()
            endif

            write(file_unit2, *) "# number_of_good_fits(j), k, j, mc_random_erosion_rate(k)"

            open(file_unit3, file="good_fit_stats.txt", status="unknown", iostat=io_status)
            if (io_status /= 0) then
                call log_message("error while opening file: good_fit_stats.txt")
                call log_message("error code: c9a546eeb406bdddda2c92440f264b3e")
                call clean_up_stop()
            endif

            write(file_unit3, *) "# number_of_good_fits(j), k, j, chi_squared(k, j), predicted_ages(k, j)"

            do i = 0, config%mpi_total_num_of_cpu - 1
                write(suffix, "(I4.4)") i
                folder_name = "mpi_mc_process_" // suffix // "/output/mc_results/"
                call log_message("reading folder: " + folder_name)

                open(file_unit, file=trim(folder_name) // "erate_gfit.txt", iostat=io_status)
                if (io_status /= 0) then
                    call log_message("error while opening file: " // trim(folder_name) // "good_fit_stats.txt")
                    call log_message("error code: b3093947b8e26d094546d931245242ee")
                    call clean_up_stop()
                endif

                ! skip header
                read(file_unit, "(a)") global_line

                do
                    read(file_unit, *, iostat=io_status) column1, column2, column3, column4
                    if (is_iostat_end(io_status)) then
                        ! end of file reached
                        exit
                    endif
                    write(file_unit2, *) column1, column2, column3, column4
                enddo
                close(file_unit)

                open(file_unit, file=trim(folder_name) // "good_fit_stats.txt", iostat=io_status)
                if (io_status /= 0) then
                    call log_message("error while opening file: good_fit_stats.txt")
                    call log_message("error code: f32b4303cfcc34d2f43ec6ff65be2de8")
                    call clean_up_stop()
                endif

                ! skip header
                read(file_unit, "(a)") global_line

                do
                    read(file_unit, *, iostat=io_status) column1, column2, column3, column4, column5
                    if (is_iostat_end(io_status)) then
                        ! end of file reached
                        exit
                    endif
                    write(file_unit3, *) column1, column2, column3, column4, column5
                enddo
                close(file_unit)
            enddo

            close(file_unit3)
            close(file_unit2)

        endif

        ! wait for all processes to reach this point (barrier)
        ! call MPI_Barrier(MPI_COMM_WORLD, mc_mpi_error)

        call clean_up()

    end subroutine monte_carlo_erosion

    subroutine read_csv_input_file(config)
        implicit none

        type(config_t) :: config

        integer(4) :: i, j, num_of_columns

        call log_message("open csv input file: " + config%mc_csv_input_file)
        open(file_unit, file=trim(config%mc_csv_input_file), status="old", action="read", iostat=io_status)

        if (io_status /= 0) then
            call log_message("error while opening file: " + trim(config%mc_csv_input_file))
            call log_message("error code: b0e18abea50bedc5b5b513ed51fd241c")

            call clean_up_stop()
        endif

        ! skip first line (header)
        read(file_unit, "(a)", iostat=io_status) global_line

        ! read line containing actual data in order to determine number of columns
        read(file_unit, "(a)", iostat=io_status) global_line

        if (io_status /= 0) then
            call log_message("error reading csv input file: " + trim(config%mc_csv_input_file))
            call log_message("line: '" // trim(global_line) // "'")
            call log_message("error code: 1e7d8c286e6995d2e22b4e646ca43039")

            call clean_up_stop()
        endif

        num_of_columns = 1

        call log_message("global_line: " + trim(global_line))

        do i = 1, max_line_length
            if (global_line(i:i) == csv_delimiter) then
                num_of_columns = num_of_columns + 1
            endif
        enddo

        num_of_data_sets_in_line = num_of_columns / 3

        call log_message("num_of_columns: " + num_of_columns + ", num_of_data_sets_in_line: " + num_of_data_sets_in_line)

        num_of_chronometers = 1

        ! determine number of data lines in file
        do
            read(file_unit, "(a)", iostat=io_status) global_line

            if (io_status < 0) then
                ! end of file
                exit
            else if (io_status > 0) then
                call log_message("error reading csv input file: " + trim(config%mc_csv_input_file))
                call log_message("line: '" // trim(global_line) // "'")
                call log_message("error code: 354e2e99cdc06bafd51f80f8ec6cee00")

                call clean_up_stop()
            endif

            num_of_chronometers = num_of_chronometers + 1
        enddo

        call log_message("num_of_chronometers: " + num_of_chronometers)

        allocate(all_chronometers(num_of_chronometers))

        rewind(file_unit)

        ! skip first line (header)
        read(file_unit, "(a)", iostat=io_status) global_line

        do i = 1, num_of_chronometers
            call log_message("process line: " + i)

            allocate(all_chronometers(i)%all_tuples(num_of_data_sets_in_line))

            read(file_unit, "(a)", iostat=io_status) global_line

            if (io_status /= 0) then
                call log_message("error reading csv input file: " + trim(config%mc_csv_input_file))
                call log_message("line: '" // trim(global_line) // "'")
                call log_message("error code: 38d777624c59ce47474d7b3e52bd919d")

                call clean_up_stop()
            endif

            call get_line_data(config, num_of_data_sets_in_line, i)
        enddo

        open(file_unit2, file="chrono_data.txt", status="unknown", iostat=io_status)

        if (io_status /= 0) then
            call log_message("error while opening file: chrono_data.txt")
            call log_message("error code: d13cc0fb107fe952c3af9dcb5fb31c8a")

            call clean_up_stop()
        endif

        ! write header
        write(file_unit2, *) "# i, j, short_name, age, sigma"

        do i = 1, num_of_chronometers
            call log_message("chronometer id: " + all_chronometers(i)%id)
            do j = 1, num_of_data_sets_in_line
                if (config%mc_check_min_threshold) then
                    if (all_chronometers(i)%all_tuples(j)%sigma < all_chronometers(i)%all_tuples(j)%age * &
                        config%mc_min_threshold_factor) then
                        all_chronometers(i)%all_tuples(j)%sigma = all_chronometers(i)%all_tuples(j)%age * &
                            config%mc_min_threshold_factor
                        call log_message("sigma changed")
                    endif
                endif
                call log_message("tuple: " + all_chronometers(i)%all_tuples(j)%age + ", " + all_chronometers(i)%all_tuples(j)%sigma)
                write(file_unit2, *) i, j, all_chronometers(i)%short_name, all_chronometers(i)%all_tuples(j)%age, &
                    all_chronometers(i)%all_tuples(j)%sigma
            enddo
        enddo

        close(file_unit2)

        close(file_unit)

    end subroutine read_csv_input_file

    subroutine get_next_column(pos1, pos2)
        implicit none

        integer(4), intent(in) :: pos1
        integer(4), intent(out) :: pos2

        do pos2 = pos1, max_line_length
            if (global_line(pos2:pos2) == csv_delimiter) then
                exit
            endif
        enddo
    end subroutine get_next_column

    subroutine get_next_column_data(pos1, pos2, column_data)
        implicit none

        integer(4), intent(inout) :: pos1
        integer(4), intent(out) :: pos2
        character(*), intent(out) :: column_data

        call get_next_column(pos1, pos2)
        column_data = global_line(pos1:pos2-1)
        pos1 = pos2 + 1

    end subroutine get_next_column_data

    subroutine get_line_data(config, num_of_data_sets_in_line, current_line)
        implicit none

        type(config_t) :: config

        integer(4), intent(in) :: num_of_data_sets_in_line, current_line

        integer(4) :: i, pos1, pos2

        character(max_line_length) :: column_data

        type(tuple_age_sigma_t) :: tuple_age_sigma

        pos1 = 1
        pos2 = 1

        do i = 1, num_of_data_sets_in_line

            call get_next_column_data(pos1, pos2, column_data)

            all_chronometers(current_line)%id = get_chronometer_id(column_data)
            if (len(trim(column_data)) > chronometer_name_length) then
                all_chronometers(current_line)%short_name = trim(column_data(1:chronometer_name_length))
            else
                all_chronometers(current_line)%short_name = trim(column_data)
            endif

            !call log_message("chronometer_data%id: " + all_chronometers(current_line)%id)

            call get_next_column_data(pos1, pos2, column_data)

            !call log_message("column_data: " + column_data)
            read(column_data, "(D30.10)", iostat=io_status) tuple_age_sigma%age

            if (io_status /= 0) then
                call log_message("error reading csv input file: " + trim(config%mc_csv_input_file))
                call log_message("line: '" // trim(global_line) // "'")
                call log_message("column data: " + trim(column_data))
                call log_message("not a valid floating point value.")
                call log_message("perhaps replace ',' with '.'")
                call log_message("error code: 821a1de8706255bdc75c8cd8b5dfe8d9")

                call clean_up_stop()
            endif

            !call log_message("tuple_age_sigma%age: " + tuple_age_sigma%age)

            call get_next_column_data(pos1, pos2, column_data)

            !call log_message("column_data: " + column_data)
            read(column_data, "(D30.10)", iostat=io_status) tuple_age_sigma%sigma

            if (io_status /= 0) then
                call log_message("error reading csv input file: " + trim(config%mc_csv_input_file))
                call log_message("line: '" // trim(global_line) // "'")
                call log_message("column data: " + trim(column_data))
                call log_message("not a valid floating point value")
                call log_message("perhaps replace ',' with '.'")
                call log_message("error code: 1c92d9eb2c6cda1ba40c68fc88f05df3")

                call clean_up_stop()
            endif

            !call log_message("tuple_age_sigma%sigma: " + tuple_age_sigma%sigma)

            !call log_message("tuple_age_sigma: " + tuple_age_sigma)

            all_chronometers(current_line)%all_tuples(i) = tuple_age_sigma

        enddo
    end subroutine get_line_data

    function get_chronometer_id(chronometer_name)
        implicit none

        integer(4) :: get_chronometer_id, i, c

        character(max_line_length), intent(in) :: chronometer_name
        character(max_line_length) :: lower_name

        do i = 1, max_line_length
            c = iachar(chronometer_name(i:i))

            if ((c >= iachar("A")) .and. (c <= iachar("Z"))) then
                lower_name(i:i) = achar(c + 32)
            else
                lower_name(i:i) = achar(c)
            endif
        enddo

        select case(trim(lower_name))
            case ("aphe")
                get_chronometer_id = 1
            case ("zrnhe")
                get_chronometer_id = 8
            case ("apft")
                get_chronometer_id = 9
            case ("zrnft")
                get_chronometer_id = 10
            case ("msar")
                get_chronometer_id = 11
            case ("ksparar")
                get_chronometer_id = 12
            case ("btar")
                get_chronometer_id = 13
            case ("hblar")
                get_chronometer_id = 15
            case ("appb")
                get_chronometer_id = 16
            case ("btrbst")
                get_chronometer_id = 17
            case ("rtpb")
                get_chronometer_id = 18
            case ("ttnpb")
                get_chronometer_id = 19
            case ("zrnpb")
                get_chronometer_id = 20
            case ("ttnhe")
                get_chronometer_id = 21
            case default
                call log_message("unknown chronometer: " + trim(chronometer_name))
                call log_message("error code: 65c9217ba1b0d99e8896c2445f0d50db")

                call clean_up_stop()
        endselect


    end function get_chronometer_id


!    subroutine get_line_values()
!        implicit none
!
!        integer(4) :: pos1, pos2
!
!    end subroutine

    subroutine clean_up()
        implicit none

        logical :: file_is_opened

        integer(4) :: i


        if (allocated(mc_random_erosion_rate)) then
            deallocate(mc_random_erosion_rate)
        endif

        if (allocated(random_difs)) then
            deallocate(random_difs)
        endif

        if (allocated(topo_file_name)) then
            deallocate(topo_file_name)
        endif

        if (allocated(all_chronometers)) then
            do i = 1, num_of_chronometers
                if (allocated(all_chronometers(i)%all_tuples)) then
                    deallocate(all_chronometers(i)%all_tuples)
                endif
            enddo

            deallocate(all_chronometers)
        endif

        if (allocated(chi_squared)) then
            deallocate(chi_squared)
        endif

        if (allocated(predicted_ages)) then
            deallocate(predicted_ages)
        endif

        if (allocated(number_of_good_fits)) then
            deallocate(number_of_good_fits)
        endif

        inquire(file_unit, opened=file_is_opened)

        if (file_is_opened) then
            close(file_unit, iostat=io_status)

            if (io_status /= 0) then
                call log_message("cloud not close file, maybe already closed")
            endif
        endif

        inquire(file_unit2, opened=file_is_opened)

        if (file_is_opened) then
            close(file_unit2, iostat=io_status)

            if (io_status /= 0) then
                call log_message("cloud not close file, maybe already closed")
            endif
        endif

        inquire(file_unit3, opened=file_is_opened)

        if (file_is_opened) then
            close(file_unit3, iostat=io_status)

            if (io_status /= 0) then
                call log_message("cloud not close file, maybe already closed")
            endif
        endif
    end subroutine clean_up

    subroutine clean_up_stop()
        implicit none

        call clean_up()
        call MPI_Finalize(mc_mpi_error)

        error stop 1
    end subroutine clean_up_stop

end module m_monte_carlo_erosion
