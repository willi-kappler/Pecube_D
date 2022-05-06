module m_test_create_velocities

    use m_move_velocities
    use m_pecube_config

    implicit none

    contains
        subroutine test_create_velocities()
            implicit none

            character(300), dimension(:), allocatable :: velocity_files

            type(t_config) :: config

            print *, "begin test_create_velocities"
            call system("date")

            allocate(velocity_files(5))

            call create_velocities(config, "input.txt", velocity_files, 5)

            deallocate(velocity_files)

            print *, "finished test_create_velocities"
            call system("date")
        end subroutine test_create_velocities

end module m_test_create_velocities
