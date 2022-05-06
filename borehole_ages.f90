module m_borehole_ages
    use m_pecube_config
    use m_logger
    use m_data_structures

    implicit none

    integer(4), parameter :: borehole_file_unit = 89

    private

    public load_borehole_ages

  contains

    function load_borehole_ages(borehole_ages_file, world_x, world_y, world_z, borehole_ages_points)
      implicit none

      character(str_length), intent(in) :: borehole_ages_file
      real(8), intent(in) :: world_x, world_y, world_z
      type(vector3D_t), intent(out), allocatable, dimension(:) :: borehole_ages_points

      integer(4) :: file_name_length, line_counter, io_error, io_status, i
      integer(4) :: load_borehole_ages
      logical :: file_exists
      character(300) :: input_line
      real(8) :: x, y, z

      file_name_length = len_trim(borehole_ages_file)

      if (file_name_length == 0) then
        call log_message("borehole_ages_file not given in input file")
        load_borehole_ages = 0
        return
      end if

      ! Try to open and read the file:
      inquire(file=borehole_ages_file, exist=file_exists, iostat=io_error)

      if (io_error /= 0) then
          call log_message("error opening file: " // trim(borehole_ages_file))
          call log_message("io_error: " + io_error)
          call log_message("error code: b678e69822cead8421f140ef4ea7ffae")
          error stop 1
      end if

      if (.not. file_exists) then
          call log_message("error opening file: file does not exist: " // trim(borehole_ages_file))
          call log_message("error code: 310aebc856fb76315bf2c7a03a256f5b")
          error stop 1
      end if

      open(borehole_file_unit, file=borehole_ages_file, iostat=io_error, status="old", action="read")

      if (io_error /= 0) then
          call log_message("error opening file: " // trim(borehole_ages_file))
          call log_message("io_error: " + io_error)
          call log_message("error code: 42987b2ecea269a90d423b36efb19c56")
          error stop 1
      end if

      line_counter = 0

      do
          read(borehole_file_unit, *, iostat=io_status) input_line
          if (io_status < 0) then
              exit ! end of file reached
          endif

          line_counter = line_counter + 1
      enddo

      call log_message("load_borehole_ages, number of lines: " + line_counter)

      rewind(borehole_file_unit)

      allocate(borehole_ages_points(line_counter))

      ! Read in borehole coordinates in [km] and real world coordinates!
      do i = 1, line_counter
        read(borehole_file_unit, *) x, y, z
        ! Initialize last time step
        borehole_ages_points(i)%x = x - world_x
        borehole_ages_points(i)%y = y - world_y
        borehole_ages_points(i)%z = z + world_z
      enddo

      close(borehole_file_unit)

      load_borehole_ages = line_counter

    end function load_borehole_ages
end module m_borehole_ages
