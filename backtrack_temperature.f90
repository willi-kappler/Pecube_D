module m_backtrack_temperature
contains

  subroutine backtrack_temperature(start_step, moving_points, model_points, time_temperature_history, &
        velocity_info, output_dir, config)
    use m_logger
    use m_data_structures
    use m_find_temperature
    use m_pecube_config

    implicit none

    integer(4), parameter :: file_unit_temperature = 84
    integer(4), parameter :: file_unit_temperature_history = 85
    integer(4), parameter :: file_unit_velocity_info = 87
    integer(4), intent(in) :: start_step, moving_points, model_points
    integer(4) :: i, j, k, l, current_step, ntime, sub_step, nnode

    character(*), intent(in) :: time_temperature_history, velocity_info, output_dir
    character(4) :: file_id2

    real(8), dimension(:, :), allocatable :: temperature_history_sub
    real(8), dimension(:), allocatable :: tx, ty, tz, t_in
    real(8) :: dt_sub, time_history_value, tnow

    logical :: file_exists

    type(config_t), intent(in) :: config

    type(vector3D_t), dimension(:, :), allocatable :: point_pos
    ! velocity doesn't change between sub time steps
    type(vector3D_t), dimension(:), allocatable :: point_velo

    ! In order to work this backtracking needs the velocity (uplift rate)
    ! at the beginning and at the end:

    ! 130  1  0  1  1  0 0  0  0  0  0  0  0 0 0 0 0
    ! 129  1  0  1  1  0.2  0  0  0  0  0  0 0 0 0 0
    ! 12   1  0  1  1  0.2  0  0  0  0  0  0 0 0 0 0
    ! 11   1  0  1  1  0.3  0  0  0  0  0  0 0 0 0 0
    ! 0.0  1  0  1  1  0.3  0  0  0  0  0  0 0 0 0 0


    inquire(file=trim(time_temperature_history), exist=file_exists)

    if (file_exists .and. config%use_cached_files) then
      call log_message("Using already pre-calculated time_temperature_history files")
      return
    endif

    call log_message("backtrack_temperature.f90, start_step: " + start_step + ", moving_points: " + moving_points + &
      ", model_points: " + model_points)
    !call log_message("backtrack_temperature.f90, time_temperature_history: " + time_temperature_history)
    !call log_message("backtrack_temperature.f90, velocity_info: " + velocity_info)
    !call log_message("backtrack_temperature.f90, output_dir: " + output_dir)

    ! This file contains the whole history for time step "start_step"
    ! # point id, x, y, z, vx, vy, vz
    open(file_unit_velocity_info, file=trim(velocity_info), status="old", form="unformatted", access="stream")

    read(file_unit_velocity_info) j, k

    !call log_message("backtrack_temperature.f90, from velocity file: j: " + j + ", k: " + k)

    if (j /= start_step) then
      call log_message("Error reading in velocity info: start_step != j")
      call log_message("start_step: " + start_step + ", j: " + j)
      call log_message("error code: 4ac7975fe6392c51d4c47d2de976b1fd")
      error stop 1
    endif

    if (k /= moving_points) then
      call log_message("Error reading in velocity info: moving_points != k")
      call log_message("moving_points: " + moving_points + ", k: " + k)
      call log_message("error code: 0676aad49739803f4f05ff67ee6559d5")
      error stop 1
    endif

    open(file_unit_temperature_history, file=trim(time_temperature_history), status="unknown", &
        form="unformatted", access="stream")

    write(file_unit_temperature_history) moving_points, start_step

    allocate(tx(model_points), ty(model_points), tz(model_points), t_in(model_points))

    !call log_message("backtrack_temperature.f90, start back track loop")

    ! Go back in time
    do current_step = start_step, 1, -1
      !call log_message("backtrack_temperature.f90, current_step: " + current_step)

      read(file_unit_velocity_info) j

      if (j /= current_step) then
        call log_message("Error reading in velocity info: current_step != j")
        call log_message("current_step: " + current_step + ", j: " + j)
        call log_message("error code: fdb5a493fcd1756e39117ee90502582a")
        error stop 1
      endif

      write(file_id2, "(i4.4)") current_step

      ! #node id, x, y, z, t
      open(file_unit_temperature, file=trim(output_dir)//"/temperature_field_sub_"//file_id2//".bin", &
          status="old", form="unformatted", access="stream")
      ! Skip header
      !read(file_unit_temperature, *)

      read(file_unit_temperature) ntime, j, nnode

      !call log_message("backtrack_temperature.f90, ntime: " + ntime)

      if (j /= current_step) then
        call log_message("Error reading in temperature info: current_step != j")
        call log_message("current_step: " + current_step + ", j: " + j)
        call log_message("error code: 93629db3e700c2e2a3d459598ddf3444")
        error stop 1
      endif

      if (nnode /= model_points) then
        call log_message("Error reading in temperature info: nnode != model_pos")
        call log_message("nnode: " + nnode + ", model_points: " + model_points)
        call log_message("error code: 0ec5bfade8f85b6128be346514da6777")
        error stop 1
      endif

      write(file_unit_temperature_history) current_step, ntime
      allocate(temperature_history_sub(ntime, moving_points))
      allocate(point_pos(ntime + 1, moving_points), point_velo(moving_points))

      !call log_message("backtrack_temperature.f90, start point loop")

      do i = 1, moving_points
        ! Read in the already back tracked positions and velocities
        read(file_unit_velocity_info) k, point_pos(1, i)%x, point_pos(1, i)%y, point_pos(1, i)%z, &
             point_velo(i)%x, point_velo(i)%y, point_velo(i)%z

        if (k /= i) then
          call log_message("Error reading in velocity info: k != i")
          call log_message("k: " + k + ", i: " + i)
          call log_message("error code: e23b8613cb8241fad4908ea9ac2d94e2")
          error stop 1
        endif
      enddo ! i = 1, moving_points

      !call log_message("backtrack_temperature.f90, start sub step loop")

      do sub_step = 1, ntime
        read(file_unit_temperature) dt_sub, k, time_history_value

        !call log_message("dt_sub: " + dt_sub + ", k: " + k + ", time_history_value: " + time_history_value)

        if (k /= sub_step) then
          call log_message("Error reading in temperature info: k != sub_step")
          call log_message("k: " + k + ", sub_step: " + sub_step)
          call log_message("error code: f0f6b80d99d6317c32521913ea2f4f13")
          error stop 1
        endif

        write(file_unit_temperature_history) sub_step, time_history_value

        do i = 1, model_points
          read(file_unit_temperature) l, tx(i), ty(i), tz(i), t_in(i)

          ! call log_message("i: " + i + ", tx: " + tx(i) + ", ty: " + ty(i) + ", tz: " + tz(i) + ", t_in: " + t_in(i))

          if (l /= i) then
            call log_message("Error reading in temperature info: l != i")
            call log_message("l: " + l + ", i: " + i)
            call log_message("error code: 695365b60e9e50d901c6aa728331cde2")
            error stop 1
          endif
        enddo ! i = 1, model_points

        ! Calculate temperature for selected points:
        !call log_message("sub_step: " + sub_step)
        !call find_temperature(0.0_8, 0.0_8, 0.0_8, tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 0): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 32.0_8, tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 32): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 33.0_8, tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 33): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 34.0_8, tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 34): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 34.5_8, tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 34.5): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 35.0_8, tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 35.0): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 110.0_8 - 0.2_8, tx, ty, tz, t_in, model_points, radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 110 - 0.2): " + tnow)
        !call find_temperature(0.0_8, 0.0_8, 110.0_8 - 0.0_8, tx, ty, tz, t_in, model_points, radius, tnow, .true.)
        !call log_message("temperature at (0, 0, 110 - 0): " + tnow)

        ! Now calculate the temperature for each surface node position
        !$omp parallel default( shared ) &
        !$omp private( i, tnow )
        !$omp do schedule( static )
        do i = 1, moving_points
          call find_temperature(point_pos(sub_step, i)%x, point_pos(sub_step, i)%y, point_pos(sub_step, i)%z, &
              tx, ty, tz, t_in, model_points, config%temperature_radius, tnow, .false.)
          temperature_history_sub(sub_step, i) = tnow
          !if (i == 1) then
          !call log_message("backtrack_temperature.f90: i: " + i + ", temperature: " + tnow + ", pz: " + point_pos(sub_step, i)%z)
          !end if
          ! Move surface point in sub time step
          point_pos(sub_step + 1, i)%x = point_pos(sub_step, i)%x + (point_velo(i)%x * dt_sub)
          point_pos(sub_step + 1, i)%y = point_pos(sub_step, i)%y + (point_velo(i)%y * dt_sub)
          point_pos(sub_step + 1, i)%z = point_pos(sub_step, i)%z + (point_velo(i)%z * dt_sub)
          !if (i == 1) then
          !call log_message("backtrack_temperature.f90: pz after move: " + point_pos(sub_step + 1, i)%z + ", vz: " + point_velo(i)%z)
          !end if
        enddo ! i = 1, moving_points
        !$omp end do
        !$omp end parallel

      enddo ! sub_step = 1, ntime

      !call log_message("backtrack_temperature.f90, start temperature sub step loop")

      ! Write pre-calculated temperature value for each node in each sub-time step to the output file
      do sub_step = 1, ntime
        do i = 1, moving_points
          write(file_unit_temperature_history) i, temperature_history_sub(sub_step, i), &
            point_pos(sub_step, i)%x, point_pos(sub_step, i)%y, point_pos(sub_step, i)%z, &
            point_velo(i)%x, point_velo(i)%y, point_velo(i)%z
        enddo
      enddo

      !call log_message("backtrack_temperature.f90, finish temperature sub step loop")

      deallocate(temperature_history_sub)

      !call log_message("backtrack_temperature.f90, deallocated temperature_history_sub")

      deallocate(point_pos, point_velo)

      !call log_message("backtrack_temperature.f90, deallocated point_pos, point_velo")

      close(file_unit_temperature)
    enddo ! current_step = start_step, 1, -1

    !call log_message("backtrack_temperature.f90, Finish step back track loop")

    deallocate(tx, ty, tz, t_in)
    close(file_unit_temperature_history)
    close(file_unit_velocity_info)
  end subroutine backtrack_temperature

end module m_backtrack_temperature
