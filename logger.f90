module m_logger
    implicit none

    private

    integer(4), parameter :: file_unit = 90
    integer(4), parameter :: max_length_of_array = 100

    integer(4) :: io_status

    integer(4), dimension(8) :: current_date_and_time

    interface operator(+)
        module procedure add_str_str
        module procedure add_str_bool
        module procedure add_str_int4
        module procedure add_str_int8
        module procedure add_str_array_int4
        module procedure add_str_array_int8
        module procedure add_str_real4
        module procedure add_str_real8
        module procedure add_str_array_real8
        module procedure add_str_array_real8_real8
    end interface operator(+)



    public log_message
    public logger_init
    public logger_finish
    public logger_flush
    public operator(+)

    contains

    subroutine log_message(message)
        implicit none

        character(*), intent(in) :: message

!        print *, "log message: ", message

        call date_and_time(values = current_date_and_time)

        write(file_unit, "('(', I4.4, '.', I2.2, '.', I2.2, ' - ', I2.2, ':', I2.2, ':', I2.2, ') ', a)") &
            current_date_and_time(1), current_date_and_time(2), current_date_and_time(3), current_date_and_time(5), &
            current_date_and_time(6), current_date_and_time(7), message


    end subroutine log_message

    subroutine logger_init(config)
        use m_pecube_config

        implicit none

        type(config_t) :: config
        character(256) :: log_filename

        call date_and_time(values = current_date_and_time)

        write(log_filename, "(I4.4, '_', I2.2, '_', I2.2, '___', I2.2, '_', I2.2, '_', I2.2, '___', I4.4, '.log')") &
            current_date_and_time(1), current_date_and_time(2), current_date_and_time(3), current_date_and_time(5), &
            current_date_and_time(6), current_date_and_time(7), config%mpi_current_cpu_id

!        print *, "log_filename: ", trim(log_filename)

        open(file_unit, file=trim(log_filename), status="unknown", iostat=io_status)

        if (io_status /= 0) then
            print *, "could not open log file: ", trim(log_filename)
            error stop 1
        endif

        config%log_filename = trim(log_filename)

    end subroutine logger_init

    subroutine logger_finish()
        implicit none

        close(file_unit)

    end subroutine logger_finish

    subroutine logger_flush()
        implicit none

        flush(file_unit)
    end subroutine

    function add_str_str(str1, str2)
        implicit none

        character(*), intent(in) :: str1
        character(*), intent(in) :: str2

        character(len(str1) + len(str2)) :: add_str_str

        add_str_str = str1 // str2
    end function add_str_str

    function add_str_int4(str1, int1)
        implicit none

        character(*), intent(in) :: str1
        integer(4), intent(in) :: int1

        character(:), allocatable :: add_str_int4
        character(12) :: value

        write(value, "(I12)") int1
        add_str_int4 = str1 // trim(adjustl(value))
    end function add_str_int4

    function add_str_int8(str1, int1)
        implicit none

        character(*), intent(in) :: str1
        integer(8), intent(in) :: int1

        character(:), allocatable :: add_str_int8
        character(20) :: value

        write(value, "(I20)") int1
        add_str_int8 = str1 // trim(adjustl(value))
    end function add_str_int8

    function add_str_bool(str1, bool1)
        implicit none

        character(*), intent(in) :: str1
        logical, intent(in) :: bool1

        character(len(str1) + 6) :: add_str_bool

        if (bool1) then
            write(add_str_bool, "(a, 'true')") str1
        else
            write(add_str_bool, "(a, 'false')") str1
        endif
    end function add_str_bool

    function add_str_array_int4(str1, array_int1)
        implicit none

        character(*), intent(in) :: str1
        integer(4), intent(in), dimension(:) :: array_int1

        character(len(str1) + (min(size(array_int1), max_length_of_array) * 14)) :: add_str_array_int4

        integer(4) :: i, end_index

        end_index = min(size(array_int1), max_length_of_array)

        write(add_str_array_int4, "(a)") str1

        do i = 1, end_index
            write(add_str_array_int4, "(a, ' ', I12)") trim(add_str_array_int4), array_int1(i)
        enddo
    end function add_str_array_int4

    function add_str_array_int8(str1, array_int1)
        implicit none

        character(*), intent(in) :: str1
        integer(8), intent(in), dimension(:) :: array_int1

        character(len(str1) + (min(size(array_int1), max_length_of_array) * 22)) :: add_str_array_int8

        integer(4) :: i, end_index

        end_index = min(size(array_int1), max_length_of_array)

        write(add_str_array_int8, "(a)") str1

        do i = 1, end_index
            write(add_str_array_int8, "(a, ' ', I20)") trim(add_str_array_int8), array_int1(i)
        enddo
    end function add_str_array_int8

    function add_str_real4(str1, real1)
        implicit none

        character(*), intent(in) :: str1
        real(4), intent(in) :: real1

        character(:), allocatable :: add_str_real4
        character(16) :: value

        write(value, "(F16.4)") real1
        add_str_real4 = str1 // trim(adjustl(value))
    end function add_str_real4

    function add_str_real8(str1, real1)
        implicit none

        character(*), intent(in) :: str1
        real(8), intent(in) :: real1

        character(:), allocatable :: add_str_real8
        character(32) :: value

        write(value, "(F32.4)") real1
        add_str_real8 = str1 // trim(adjustl(value))
    end function add_str_real8

    function add_str_array_real8(str1, array_real8)
        implicit none

        character(*), intent(in) :: str1
        real(8), intent(in), dimension(:) :: array_real8

        character(len(str1) + (min(size(array_real8), max_length_of_array) * 18)) :: add_str_array_real8

        integer(4) :: i, end_index

        end_index = min(size(array_real8), max_length_of_array)

        write(add_str_array_real8, "(a)") str1

        do i = 1, end_index
            write(add_str_array_real8, "(a, ' ', F16.4)") trim(add_str_array_real8), array_real8(i)
        enddo
    end function add_str_array_real8

    function add_str_array_real8_real8(str1, array_real8)
        implicit none

        character(*), intent(in) :: str1
        real(8), intent(in), dimension(:,:) :: array_real8

        character(len(str1) + 16) :: add_str_array_real8_real8

        ! TODO: loop over array and print values
        write(add_str_array_real8_real8, "(a, F16.4)") str1, array_real8(1,1)
    end function add_str_array_real8_real8
end module m_logger
