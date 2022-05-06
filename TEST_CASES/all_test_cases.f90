program m_test_cases

    use m_test_calculate_misfit
    use m_test_create_velocities
    use m_test_ages

    use m_pecube_config

    implicit none

    type(t_config) :: config

    print *, "Running all test cases..."
    call system("date")
    call logger_init(config)

    !call test_create_velocities()

    call test_calculate_misfit()

    call test_ages()

    print *, "Finished with all test cases"
    call system("date")
    call logger_finish()
end program
