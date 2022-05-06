module m_test_calculate_misfit

    use m_logger
    use m_calculate_misfit
    use m_pecube_config

    implicit none

    contains
        subroutine test_calculate_misfit()
            implicit none

            type(t_config) :: config

            integer(4), dimension(:,:), allocatable :: iconsurf
            real(8), dimension(:), allocatable :: zsurf
            real(8), dimension(:), allocatable :: track_length_mean


            print *, "Begin test_create_velocities"
            call system("date")

            config%nx = 10
            config%ny = 10
            config%spacing_long = 1
            config%spacing_lat = 1
            config%nskip = 1
            config%location_long = 0.0_8
            config%location_lat = 0.0_8
            config%length_comparison_file = "length_comparison_file.dat"
            config%num_elements_surf = 100
            config%nsurf = 100
            config%npe = 4

            allocate(iconsurf(config%npe, config%num_elements_surf))
            allocate(zsurf(config%nsurf))
            allocate(track_length_mean(10))


            call calculate_misfit(config, iconsurf, zsurf, track_length_mean)

            deallocate(track_length_mean)
            deallocate(zsurf)
            deallocate(iconsurf)

            print *, "Finished test_create_velocities"
            call system("date")
        end subroutine test_calculate_misfit

end module m_test_calculate_misfit
