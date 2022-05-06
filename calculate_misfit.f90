module m_calculate_misfit
    use m_logger
    use m_pecube_config

    implicit none
    save

    ! make everything private
    private

    ! internal constants
    integer(4), parameter :: global_file_unit = 80
    integer(4) :: num_of_observables

    enum, bind( C )
      enumerator :: TYPE_AFT = 1, TYPE_ZFT, END_TYPE
    end enum

    type :: sample_t
      integer(4) :: s_type
      real(8) :: longitude
      real(8) :: latitude
      real(8) :: elevation
      real(8) :: length
      real(8) :: length_error
      real(8) :: w1
      real(8) :: w2
      real(8) :: w3
      real(8) :: w4
      real(8) :: new_height
      real(8) :: length_predicted
      real(8) :: misfit
      integer(4) :: index
    end type sample_t

    type(sample_t), dimension(:), allocatable :: all_samples

    ! export identifyer
    public calculate_misfit

  contains

    subroutine load_length_data(config)
      implicit none

      type(config_t), intent(in) :: config

      integer(4) :: iostatus, i, pos, i1, j1, ndx, ndy, name_length

      real(8) :: r, s

      character(6) :: s_type
      character(255) :: line, token

      ndx = int(config%spacing_long * config%nskip)
      ndy = int(config%spacing_lat * config%nskip)

      name_length = 1
      do name_length = 1, str_length
        if (config%length_comparison_file(name_length:name_length) == char(0)) then
          exit
        endif
      enddo

      call log_message("calculate_misfit.f90: name_length: " + name_length)

      if (name_length == 1) then
        call log_message("calculate_misfit.f90: no filename given (length_comparison_file)")
        return
      endif

      call log_message("calculate_misfit.f90: trying to open file: 'input/" + trim(config%length_comparison_file) + "'")

      open (global_file_unit, file='input/' // trim(config%length_comparison_file), status='old', iostat=iostatus)
      if (iostatus /= 0) then
        call log_message("could not open input file!")
        return
      endif

      num_of_observables = 0

      do
        read(global_file_unit, *, iostat=iostatus) line

        if (is_iostat_end(iostatus)) then
          exit ! end of file reached
        endif

        num_of_observables = num_of_observables + 1
      enddo

      call log_message("number of lines to read: " + num_of_observables)

      allocate(all_samples(num_of_observables))

      rewind(global_file_unit)

      do i = 1, num_of_observables
        read(global_file_unit, "(A)", iostat=iostatus) line

        pos = 1
        call extract_token(pos, line, token)
        s_type = trim(token)

        select case(s_type)
          case("AFT")
            all_samples(i)%s_type = TYPE_AFT
          case("ZFT")
            all_samples(i)%s_type = TYPE_ZFT
          case default
            call log_message("unknown type: " + s_type)
            all_samples(i)%s_type = END_TYPE
            continue
        end select

        call extract_token(pos, line, token)
        read(token, "(F15.5)") all_samples(i)%longitude

        call extract_token(pos, line, token)
        read(token, "(F15.5)") all_samples(i)%latitude

        call extract_token(pos, line, token)
        read(token, "(F15.5)") all_samples(i)%elevation

        call extract_token(pos, line, token)
        read(token, "(F15.5)") all_samples(i)%length

        call extract_token(pos, line, token)
        read(token, "(F15.5)") all_samples(i)%length_error

        ! Calculate position of sample inside the internal model
        i1 = int((all_samples(i)%longitude - config%location_long) / ndx) + 1
        if (i1 == config%nx) i1 = config%nx - 1

        j1 = int((all_samples(i)%latitude - config%location_lat) / ndy) + 1
        if (j1 == config%ny) j1 = config%ny - 1

        all_samples(i)%index = i1 + (j1 - 1) * (config%nx - 1)

        if ((all_samples(i)%index < 1) .or. (all_samples(i)%index > config%num_elements_surf)) then
            call log_message("i1: " + i1 + ", j1: " + j1 + ", num_elements_surf: " + config%num_elements_surf)
            call log_message("index: " + all_samples(i)%index + ", nx: " + config%nx + ", ny: " + config%ny)
            call log_message("longitude: " + all_samples(i)%longitude + ", lon1: " + config%location_long)
            call log_message("latitude: " + all_samples(i)%latitude + ", la1t: " + config%location_lat)

            call log_message("error in file: " // config%length_comparison_file)
            call log_message("ensure that longitude > lon1 and latitude > lat1")

            deallocate(all_samples)
            return
        end if

        r = 2.0 * (all_samples(i)%longitude - (i1 - 1) * ndx - config%location_long) / ndx
        r = r - 1.0
        s = 2.0 * (all_samples(i)%latitude - (j1 - 1) * ndy - config%location_lat) / ndy
        s = s - 1.0

        all_samples(i)%w1 = (1.0 - r) * (1.0 - s) / 4.0
        all_samples(i)%w2 = (1.0 + r) * (1.0 - s) / 4.0
        all_samples(i)%w3 = (1.0 + r) * (1.0 + s) / 4.0
        all_samples(i)%w4 = (1.0 - r) * (1.0 + s) / 4.0

      enddo

      close(global_file_unit)

      do i = 1, num_of_observables
        print *, "type:",  all_samples(i)%s_type, &
                 ", long:", all_samples(i)%longitude, &
                 ", lat:", all_samples(i)%latitude, &
                 ", elevation:", all_samples(i)%elevation, &
                 ", length", all_samples(i)%length, &
                 ", length_error:", all_samples(i)%length_error
        print *, "index:", all_samples(i)%index, &
                 "w1:", all_samples(i)%w1, &
                 "w2:", all_samples(i)%w2, &
                 "w3:", all_samples(i)%w3, &
                 "w4:", all_samples(i)%w4
      enddo
    end subroutine load_length_data

    subroutine calculate_misfit(config, iconsurf, zsurf, track_length_mean)
      implicit none

      type(config_t), intent(in) :: config
      integer(4), dimension(:,:), intent(in) :: iconsurf
      real(8), dimension(:), intent(in) :: zsurf
      real(8), dimension(:), intent(in) :: track_length_mean


      integer(4) :: i, index, ico1, ico2, ico3, ico4
      real(8), dimension(:), allocatable :: misfit_per_age
      real(8) :: total_misfit

      call log_message("calculating misfit (calculate_misfit.f90)")

      call load_length_data(config)

      if (.not. allocated(all_samples)) then
        call log_message("error loading input file!")
        return
      endif

      allocate(misfit_per_age(END_TYPE))

      misfit_per_age = 0.0

      do i = 1, num_of_observables
        index = all_samples(i)%index
        ico1 = iconsurf(1, index)
        ico2 = iconsurf(2, index)
        ico3 = iconsurf(3, index)
        ico4 = iconsurf(4, index)

        if ((ico1 > config%nsurf) .or. &
            (ico2 > config%nsurf) .or. &
            (ico3 > config%nsurf) .or. &
            (ico4 > config%nsurf)) then
          call log_message('Error: Surface node connectivity value out of range')
          call log_message('Check that thermochronological data file coordinates are same system (degrees/utm) as model')
          stop
        endif

        all_samples(i)%new_height = ((all_samples(i)%w1 * zsurf(ico1)) + &
          (all_samples(i)%w2 * zsurf(ico2)) + &
          (all_samples(i)%w3 * zsurf(ico3)) + &
          (all_samples(i)%w4 * zsurf(ico4))) * 1000.0

        all_samples(i)%length_predicted = (all_samples(i)%w1 * track_length_mean(ico1)) + &
          (all_samples(i)%w2 * track_length_mean(ico2)) + &
          (all_samples(i)%w3 * track_length_mean(ico3)) + &
          (all_samples(i)%w4 * track_length_mean(ico4))

        all_samples(i)%misfit = (all_samples(i)%length - all_samples(i)%length_predicted)**2 / &
          all_samples(i)%length_error**2

        misfit_per_age(all_samples(i)%s_type) = misfit_per_age(all_samples(i)%s_type) + all_samples(i)%misfit

      enddo

      total_misfit = 0.0_8

      do i = 1, END_TYPE - 1
        misfit_per_age(i) = sqrt(misfit_per_age(i))
        total_misfit = total_misfit + misfit_per_age(i)
        print *, "i, ", i, ", misfit per age: ", misfit_per_age(i)
      enddo

      print *, "total_misfit: ", total_misfit

      deallocate(misfit_per_age)
      deallocate(all_samples)
    end subroutine calculate_misfit

    subroutine extract_token(pos, input, token)
      implicit none

      integer(4), intent(inout) :: pos

      character(*), intent(in) :: input
      character(*), intent(out) :: token

      integer(4) :: i, mode

      i = 1
      token = " "

      !print *, "start pos: ", pos

      ! ignore all white spaces at the beginning
      mode = 1

      do
        if (mode == 1) then ! skip white space
          if (input(pos:pos) == " ") then
            pos = pos + 1
            continue
          else
            ! first non-white space character found: we need to copy it
            mode = 2
          endif
        else if (mode == 2) then ! copy character
          if (input(pos:pos) == " ") then
            ! first white space character found after token: return
            ! print *, "pos:", pos, ", token: '" // trim(token) // "'"
            return
          else
            ! copy all non-white space character to token string
            token(i:i) = input(pos:pos)
            pos = pos + 1
            i = i + 1
          end if
        else
          print *, "unknown mode: ", mode
          return
        endif
      enddo
    end subroutine extract_token

end module m_calculate_misfit
