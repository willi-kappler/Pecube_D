module m_dynamic_thermal_conductivity
    use m_pecube_config
    use m_logger
    use m_age_algorithms
    use m_data_structures

    implicit none

    private

    enum, bind( C )
        enumerator :: CONDUCTIVITY, TIME
    end enum

    integer(4), parameter :: thermal_file_unit = 80

    type tc_entry_t
      real(8) :: td_value
      real(8) :: z1
      real(8) :: z2
    end type tc_entry_t

    type thermal_conductivity_t
      integer(4) :: num_of_layers
      real(8) :: time
      type(tc_entry_t), dimension(:), allocatable :: layers
    end type thermal_conductivity_t

    public load_thermal_conductivity_file
    public deallocate_thermal_conductivity
    public get_conductivity
    public calculate_conductivity
    public thermal_conductivity_t

  contains

    subroutine load_thermal_conductivity_file(thermal_conductivity_file, dynamic_thermal_conductivity, &
      time_steps_pecube, model_depth, crustal_density, heat_capacity)
      implicit none

      character(str_length), intent(in) :: thermal_conductivity_file
      integer(4), intent(in) :: time_steps_pecube
      real(8), intent(in) :: model_depth, crustal_density, heat_capacity

      type(thermal_conductivity_t), intent(inout), dimension(:), allocatable :: dynamic_thermal_conductivity

      integer(4) :: file_name_length, io_error, colon_pos, num_of_layers, current_mode, current_time, i, j
      integer(4) :: time_steps
      logical :: file_exists, keyword_found
      character(str_length) :: current_line
      character(5) :: keyword
      real(8) :: time_value, conductivity_value, z, depth

      type(key_value_t), dimension(20) :: key_value_store

      time_steps = time_steps_pecube + 1

      call log_message("model depth: " + model_depth)
      call log_message("time_steps: " + time_steps)

      call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)

      file_name_length = len_trim(thermal_conductivity_file)

      ! call log_message("thermal_conductivity_file: " + thermal_conductivity_file)
      ! call log_message("len1: " + len(thermal_conductivity_file))
      ! call log_message("len2: " + len_trim(thermal_conductivity_file))
      ! call log_message("scan: " + scan(thermal_conductivity_file, "\0"))
      ! call log_message("index: " + index(thermal_conductivity_file, "\0"))

      if (file_name_length == 0) then
        call log_message("thermal_conductivity_file not given in input file")
        return
      end if

      ! Try to open and read the file:
      inquire(file=thermal_conductivity_file, exist=file_exists, iostat=io_error)

      if (io_error /= 0) then
          call log_message("error opening file: " // trim(thermal_conductivity_file))
          call log_message("io_error: " + io_error)
          call log_message("error code: 7a43711d07d370dab3f10564c5ca88a1")
          error stop 1
      end if

      if (.not. file_exists) then
          call log_message("error opening file: file does not exist: " // trim(thermal_conductivity_file))
          call log_message("error code: e3f6495d6a11a4a65fb4b3e5bb8ccc02")
          error stop 1
      end if

      open(thermal_file_unit, file=thermal_conductivity_file, iostat=io_error, status="old", action="read")

      if (io_error /= 0) then
          call log_message("error opening file: " // trim(thermal_conductivity_file))
          call log_message("io_error: " + io_error)
          call log_message("error code: 1785aef21da6363367395b9565bbf17f")
          error stop 1
      end if

      allocate(dynamic_thermal_conductivity(time_steps))

      num_of_layers = 0
      current_time = 1
      current_mode = CONDUCTIVITY

      do
        current_line(1:str_length) = " "
        read(thermal_file_unit, "(a)", iostat=io_error) current_line

        if (is_iostat_end(io_error)) then
            call log_message("end of file reached")
            exit
        else if (io_error /= 0) then
            call log_message("error reading from file: " // trim(thermal_conductivity_file))
            call log_message("line: " // trim(current_line))
            call log_message("io_error: " + io_error)
            call log_message("error code: 25ed813c791075f00e5025190a8627ee")
            call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
            close(thermal_file_unit)
            error stop 1
        end if

        if (current_line(1:1) == " ") then
          ! Ignore blank lines
          cycle
        end if

        colon_pos = scan(current_line, ":")

        if (colon_pos == 0) then
          ! No colon found, go to next line
          cycle
        else if (colon_pos < 3) then
          ! Invalid line
          call log_message("error reading from file: " // trim(thermal_conductivity_file))
          call log_message("line: " // trim(current_line))
          call log_message("keyword too short")
          call log_message("error code: 85eb373b2181b7c9cd75ac4474a428af")
          call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
          close(thermal_file_unit)
          error stop 1
        endif

        keyword = current_line(1 : colon_pos - 1)

        call log_message("keyword: " + keyword + ", colon_pos: " + colon_pos)

        select case (current_mode)
          case (CONDUCTIVITY)
            if (keyword == "time") then
              current_mode = TIME
              read (current_line(colon_pos + 1:), *) time_value
              call log_message("time_value: " + time_value)
              dynamic_thermal_conductivity(current_time)%time = time_value
              dynamic_thermal_conductivity(current_time)%num_of_layers = num_of_layers
              allocate(dynamic_thermal_conductivity(current_time)%layers(num_of_layers))
            else if (keyword(1:1) == "k") then
              num_of_layers = num_of_layers + 1
              call log_message("num_of_layers: " + num_of_layers)
              read (current_line(colon_pos + 1:), *) conductivity_value
              call log_message("time_value: " + time_value)
              key_value_store(num_of_layers)%key = keyword
              ! Calculate diffusivity from conductivity
              ! thermal conductivity [W/m K]
              ! specific heat capacity [J/kg K]
              ! rock density [kg/m^3]
              key_value_store(num_of_layers)%value = conductivity_value * sec_per_mil_year / &
                ((crustal_density * heat_capacity) * 1000000.0)
            else
              call log_message("error reading from file: " // trim(thermal_conductivity_file))
              call log_message("line: " // trim(current_line))
              call log_message("invalid keyword: " + keyword)
              call log_message("error code: 9bfd3f8b711e5f1662fd7436ba11e986")
              call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
              close(thermal_file_unit)
              error stop 1
            end if
          case (TIME)
            if (keyword == "time") then
              current_time = current_time + 1
              read (current_line(colon_pos + 1:), *) time_value
              dynamic_thermal_conductivity(current_time)%time = time_value
              dynamic_thermal_conductivity(current_time)%num_of_layers = num_of_layers
              allocate(dynamic_thermal_conductivity(current_time)%layers(num_of_layers))
              call log_message("time_value: " + time_value)
            else if (keyword(1:1) == "k") then
              keyword_found = .false.

              do i = 1, num_of_layers
                if (key_value_store(i)%key == keyword) then
                  keyword_found = .true.
                  exit
                end if
              end do

              if (.not. keyword_found) then
                call log_message("error reading from file: " // trim(thermal_conductivity_file))
                call log_message("line: " // trim(current_line))
                call log_message("keyword not found: " + keyword)
                call log_message("error code: 9248c4d30c2f83eaeeaca0789928fc97")
                call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
                close(thermal_file_unit)
                error stop 1
              end if

              z = 0.0
              depth = model_depth

              read (current_line(colon_pos + 1:), *, iostat=io_error) z, depth

              if (io_error /= 0) then
                ! Try one value:
                read (current_line(colon_pos + 1:), *, iostat=io_error) z
                if (io_error /= 0) then
                  ! OK, looks like an input error...
                  call log_message("error reading from file: " // trim(thermal_conductivity_file))
                  call log_message("line: " // trim(current_line))
                  call log_message("numeric value invalid!")
                  call log_message("error code: fbdb94682412df1a870e514eae85438b")
                  call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
                  close(thermal_file_unit)
                  error stop 1
                end if
              end if

              dynamic_thermal_conductivity(current_time)%layers(i)%td_value = key_value_store(i)%value
              dynamic_thermal_conductivity(current_time)%layers(i)%z1 = z
              dynamic_thermal_conductivity(current_time)%layers(i)%z2 = min(z + depth, depth)

            else
              call log_message("error reading from file: " // trim(thermal_conductivity_file))
              call log_message("line: " // trim(current_line))
              call log_message("invalid keyword: " + keyword)
              call log_message("error code: d0d98434f1319384dda636b1a0ac08b1")
              call deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
              close(thermal_file_unit)
              error stop 1
            end if
        end select

      end do

      call log_message("dynamic thermal conductivity values:")
      do i = 1, time_steps
        call log_message("time: " + dynamic_thermal_conductivity(i)%time + ", num_of_layers: " + dynamic_thermal_conductivity(i)%num_of_layers)
        do j = 1, dynamic_thermal_conductivity(i)%num_of_layers
          call log_message("td_value: " + dynamic_thermal_conductivity(i)%layers(j)%td_value)
          call log_message("z1: " + dynamic_thermal_conductivity(i)%layers(j)%z1)
          call log_message("z2: " + dynamic_thermal_conductivity(i)%layers(j)%z2)
        end do
      end do

      close(thermal_file_unit)

    end subroutine load_thermal_conductivity_file

    subroutine calculate_conductivity(current_step, current_dynamic_thermal_conductivity, dynamic_thermal_conductivity, ftime)
      implicit none

      integer(4), intent(in) :: current_step
      real(8), intent(in) :: ftime

      type(thermal_conductivity_t), intent(inout) :: current_dynamic_thermal_conductivity
      type(thermal_conductivity_t), intent(in), dimension(:), allocatable :: dynamic_thermal_conductivity

      integer(4) :: i

      if (current_step < 1) then
        do i = 1, current_dynamic_thermal_conductivity%num_of_layers
          current_dynamic_thermal_conductivity%layers(i)%td_value = dynamic_thermal_conductivity(1)%layers(i)%td_value
          current_dynamic_thermal_conductivity%layers(i)%z1 = dynamic_thermal_conductivity(1)%layers(i)%z1
          current_dynamic_thermal_conductivity%layers(i)%z2 = dynamic_thermal_conductivity(1)%layers(i)%z2
        enddo
      else
        do i = 1, current_dynamic_thermal_conductivity%num_of_layers
          current_dynamic_thermal_conductivity%layers(i)%td_value = dynamic_thermal_conductivity(current_step + 1)%layers(i)%td_value
          current_dynamic_thermal_conductivity%layers(i)%z1 = dynamic_thermal_conductivity(current_step)%layers(i)%z1 + &
            (ftime * (dynamic_thermal_conductivity(current_step + 1)%layers(i)%z1 - dynamic_thermal_conductivity(current_step)%layers(i)%z1))
          if (i == current_dynamic_thermal_conductivity%num_of_layers) then
            current_dynamic_thermal_conductivity%layers(i)%z2 = dynamic_thermal_conductivity(current_step + 1)%layers(i)%z2
          else
            current_dynamic_thermal_conductivity%layers(i)%z2 = dynamic_thermal_conductivity(current_step)%layers(i)%z2 + &
              (ftime * (dynamic_thermal_conductivity(current_step + 1)%layers(i)%z2 - dynamic_thermal_conductivity(current_step)%layers(i)%z2))
          endif
        enddo
      endif
    end subroutine calculate_conductivity

    subroutine get_conductivity(dynamic_thermal_conductivity, z, current_diffusivity)
      implicit none

      real(8), intent(in) :: z
      real(8), intent(out) :: current_diffusivity

      type(thermal_conductivity_t), intent(in) :: dynamic_thermal_conductivity

      integer(4) :: i
      real(8) :: z1, z2

      do i = 1, dynamic_thermal_conductivity%num_of_layers
        z1 = dynamic_thermal_conductivity%layers(i)%z1
        z2 = dynamic_thermal_conductivity%layers(i)%z2
        if (z2 > z1) then
          if ((z1 <= z) .and. (z <= z2)) then
            current_diffusivity = dynamic_thermal_conductivity%layers(i)%td_value
          end if
        end if
      end do
    end subroutine get_conductivity

    subroutine deallocate_thermal_conductivity(dynamic_thermal_conductivity, time_steps)
      implicit none

      integer(4), intent(in) :: time_steps

      type(thermal_conductivity_t), dimension(:), allocatable :: dynamic_thermal_conductivity

      integer(4) :: i

      if (allocated(dynamic_thermal_conductivity)) then
        do i = 1, time_steps
          if (allocated(dynamic_thermal_conductivity(i)%layers)) then
            deallocate(dynamic_thermal_conductivity(i)%layers)
          end if
        end do

        deallocate(dynamic_thermal_conductivity)
      end if
    end subroutine deallocate_thermal_conductivity

end module m_dynamic_thermal_conductivity
