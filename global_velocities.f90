!> Module for reading in velocities from input files
module m_global_velocities
    use m_compiler
    use m_logger
    use m_data_structures

    implicit none
    save

    ! make everything private
    private

    !> @struct vec
    !! @brief contains position [m] and velocities [mm/year]
    type :: vec_t
        real(8) :: px, py, pz, vx, vy, vz
    endtype vec_t

    !> Holds all velocities from current file
    type(vec_t), dimension(:), allocatable :: global_velocities

    !> Holds all velocities values for a given index
    type(vec_t), dimension(:,:,:), allocatable :: global_velocities_grid

    !> Holds all surface z values in [km] for a given index
    real(8), dimension(:,:), allocatable :: global_surface_grid

    !> Wether velocity file was successfully loaded or not
    logical :: file_loaded = .FALSE.

    !> Number of entries in current velocity file
    integer(4) :: number_of_elements = 0

    !> Number of surface points
    integer(4) :: number_of_surface_points

    !> Velocity grid step size in [km]
    real(8) :: global_gridStepX, global_gridStepY, global_gridStepZ

    !> Surface topography in [km]
    type(vector3D_t), dimension(:), allocatable :: global_surface

    !> Velocity grid dimensions [number of points]
    integer(4) :: global_sizeX, global_sizeY, global_sizeZ

    !> Minimum allowed velocity value
    real(8) :: vx_min_global, vy_min_global, vz_min_global

    !> Maximum allowed velocity value
    real(8) :: vx_max_global, vy_max_global, vz_max_global

    !> constants
    real(8), parameter :: FROM_MYRS_TO_YRS = 1e6
    real(8), parameter :: FROM_KM_TO_MM = 1e6
    real(8), parameter :: FROM_MM_RO_KM = 1e-6
    real(8), parameter :: eps = 0.0001

    integer(4), parameter :: file_unit_input = 60
    integer(4), parameter :: file_unit_cached = 61

    interface operator(+)
        module procedure add_str_vec
    end interface operator (+)





    !> Exported subroutines:
    public get_global_velocity
    public load_velocities
    public check_get_global_velocity

    contains
        !> Load all velocity entries from given filename
        !! @param filename name of the velocity file to load
        !! @param xlonmin
        !! @param xlatmin
        !! @param zl
        !! @param sizeX
        !! @param sizeY
        !! @param sizeZ
        !! @param gridStepX
        !! @param gridStepY
        subroutine load_velocities(filename, xlonmin, xlatmin, zl, sizeX, sizeY, sizeZ, &
                                   gridStepX, gridStepY, xsurf, ysurf, zsurf, nsurf, &
                                   vx_min, vy_min, vz_min, vx_max, vy_max, vz_max, use_new_velocity)
            implicit none

            character(300), intent(in) :: filename
            character(300) :: input_line

            integer(4) :: i, j, k, io_status, line_counter, id, index_x, index_y, index_z
            integer(4) :: fill_index, file_x, file_y, file_z
            integer(4), intent(in) :: sizeX, sizeY, sizeZ, nsurf
            integer(4), dimension(:,:,:), allocatable :: sample_counter

            real(8) :: px, py, pz, vx, vy, vz, sz, factor
            real(8) :: gridStepZ, fill_value, lowest_x_velocity
            real(8), intent(in) :: xlonmin, xlatmin, zl, gridStepX, gridStepY, xsurf(nsurf), ysurf(nsurf), zsurf(nsurf)
            real(8), intent(in) :: vx_min, vy_min, vz_min
            real(8), intent(in) :: vx_max, vy_max, vz_max

            type(vec_t) :: one_line, vec1
            logical :: file_exists
            logical, intent(in) :: use_new_velocity

            gridStepZ = zl / sizeZ

            call log_message("global_velocities.f90:")
            call log_message("number of nodes, x:" + sizeX + ", y:" + sizeY + ", z:" + sizeZ)
            call log_message("step in km, x:" + gridStepX + ", y:" + gridStepY + ", z:" + gridStepZ)
            call log_message("number of surface (z) values: " + nsurf)
            call log_message("xlonmin: " + xlonmin + ", xlatmin: " + xlatmin + ", zl: " + zl)
            call log_message("filename: " + filename)

            ! store current values for module global usage
            global_gridStepX = gridStepX
            global_gridStepY = gridStepY
            global_gridStepZ = gridStepZ

            global_sizeX = sizeX + 8
            global_sizeY = sizeY + 8
            global_sizeZ = sizeZ + 8

            ! store minimum and maximum velocity in global variables

            vx_min_global = vx_min
            vy_min_global = vy_min
            vz_min_global = vz_min

            vx_max_global = vx_max
            vy_max_global = vy_max
            vz_max_global = vz_max

            call log_message("minimum velocity values (vx, vy, vz): " + vx_min_global + ", " + vy_min_global + ", " + vz_min_global)
            call log_message("maximum velocity values (vx, vy, vz): " + vx_max_global + ", " + vy_max_global + ", " + vz_max_global)

            ! Deallocate 'global_velocities_grid' if it already has been allocated in last call
            if (allocated(global_velocities_grid)) then
                deallocate(global_velocities_grid)
            endif

            allocate(global_velocities_grid(global_sizeX, global_sizeY, global_sizeZ))

            global_velocities_grid = vec_t(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

            ! Deallocate 'global_surface' if it already has been allocated in last call
            if (allocated(global_surface)) then
                deallocate(global_surface)
            endif

            allocate(global_surface(nsurf))

            global_surface = vector3D_t(0.0, 0.0, 0.0)

            ! Deallocate 'global_surface_grid' if it already has been allocated in last call
            if (allocated(global_surface_grid)) then
                deallocate(global_surface_grid)
            endif

            allocate(global_surface_grid(global_sizeX, global_sizeY))

            global_surface_grid = 0.0

            ! copy surface x, y, z values to module global variable
            do i = 1, nsurf
                global_surface(i) = vector3D_t(xsurf(i), ysurf(i), zsurf(i) - zl)
            enddo

            number_of_surface_points = nsurf

            file_loaded = .FALSE.

            if (filename == "NIL") then
                call log_message("global_velocities.f90, load_velocities, 'NIL' found, no file will be loaded")
            else
                ! Determine file name lenght

                call log_message("global_velocities.f90, load_velocities, filename: " // trim(filename))

                ! Check if pre-calculated version of the file exists:
                inquire(file="input/cached_"//trim(filename), exist=file_exists)

                if (file_exists) then ! File exists, so we can use it
                    call log_message("read cached velocity file")
                    open(file_unit_cached, file="input/cached_"//trim(filename), status="old", &
                        form="unformatted", access="stream")
                    ! Read in binray file content

                    read(file_unit_cached) file_x, file_y, file_z

                    if (global_sizeX /= file_x) then
                      call log_message("error in cached velocity file:")
                      call log_message("sizeX: " + global_sizeX + " != file_x: " + file_x)
                      call log_message("error code: 19d170c8ae29d9a160a3f5f47c22ebdb")
                      error stop 1
                    endif

                    if (global_sizeY /= file_y) then
                      call log_message("error in cached velocity file:")
                      call log_message("sizeY: " + global_sizeY + " != file_y: " + file_y)
                      call log_message("error code: db59dbffa8fd5a741c475b7d7a9d0a93")
                      error stop 1
                    endif

                    if (global_sizeZ /= file_z) then
                      call log_message("error in cached velocity file:")
                      call log_message("sizeZ: " + global_sizeZ + " != file_z: " + file_z)
                      call log_message("error code: 891822b9a6ce34ad6edb334ce7e64c70")
                      error stop 1
                    endif

                    do i=1, global_sizeX
                        do j=1, global_sizeY
                            do k=1, global_sizeZ
                                read(file_unit_cached) global_velocities_grid(i, j, k)
                            enddo
                            read(file_unit_cached) global_surface_grid(i, j)
                        enddo
                    enddo

                    close(file_unit_cached)
                    file_loaded = .TRUE.
                    return
                endif

                ! Cached version does not exist, so we have to calculate the grid first
                ! Open file and check if it is readable
                call log_message("cached velocity files does not exist, calculating it...")
                open(file_unit_input, file="input/"//trim(filename),status='old',iostat=io_status)
                if (io_status /= 0) then
                    call log_message("could not open velocity file: input/"//trim(filename))
                    error stop 1
                endif

                ! count the number of lines in file.
                ! this gives us the number of entries in that file and we
                ! can allocate our array
                line_counter = 0

                ! the first line just contains the header, so we can read it in and forget it
                read(file_unit_input, *) input_line

                do
                    read(file_unit_input, *, iostat=io_status) input_line
                    if (io_status < 0) then
                        exit ! end of file reached
                    endif

                    line_counter = line_counter + 1
                enddo

                call log_message("global_velocities.f90, load_velocities, number of lines: " + line_counter)

                number_of_elements = line_counter

                ! Deallocate 'global_velocities' if it already has been allocated in last call
                if (allocated(global_velocities)) then
                    deallocate(global_velocities)
                endif

                allocate(global_velocities(number_of_elements))

                ! start rading input file form the beginning
                rewind(file_unit_input)

                ! the first line just contains the header, so we can read it in and ignore it
                read(file_unit_input, *) input_line

                do i = 1, number_of_elements
                    read(file_unit_input, *) id, px, py, pz, vx, vy, vz

                    if (vx < vx_min_global) then
                        !call log_message("vx too small, old value: " + vx + ", new value: " + vx_min_global)
                        vx = vx_min_global + eps
                    else if (vx > vx_max_global) then
                        !call log_message("vx too big, old value: " + vx + ", new value: " + vx_max_global)
                        vx = vx_max_global - eps
                    endif

                    if (vy < vy_min_global) then
                        !call log_message("vy too small, old value: " + vy + ", new value: " + vy_min_global)
                        vy = vy_min_global + eps
                    else if (vy > vy_max_global) then
                        !call log_message("vy too big, old value: " + vy + ", new value: " + vy_max_global)
                        vy = vy_max_global - eps
                    endif

                    if (vz < vz_min_global) then
                        !call log_message("vz too small, old value: " + vz + ", new value: " + vz_min_global)
                        vz = vz_min_global + eps
                    else if (vz > vz_max_global) then
                        !call log_message("vz too big, old value: " + vz + ", new value: " + vz_max_global)
                        vz = vz_max_global - eps
                    endif

                    ! position value is in km
                    ! velocity is in mm/year
                    one_line = vec_t(px - xlonmin, py - xlatmin, pz + zl, vx, vy, vz)

                    global_velocities(i) = one_line

!                        call log_message("global_velocities.f90, load_velocities, global_velocities: " + i + ", " + &
!                        global_velocities(i))
                enddo

                close(file_unit_input)

                call logger_flush()

                if (use_new_velocity) then
                  allocate(sample_counter(global_sizeX, global_sizeY, global_sizeZ))

                  sample_counter = 0

                  call log_message("New velocity calculation start")
                  do i = 1, number_of_elements
                    index_x = floor(global_velocities(i)%px / gridStepX)
                    index_y = floor(global_velocities(i)%py / gridStepY)
                    index_z = floor(global_velocities(i)%pz / gridStepZ)

                    call calculate_velo(i, index_x, index_y, index_z, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x + 1, index_y, index_z, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x, index_y + 1, index_z, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x + 1, index_y + 1, index_z, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x, index_y, index_z + 1, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x + 1, index_y, index_z + 1, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x, index_y + 1, index_z + 1, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                    call calculate_velo(i, index_x + 1, index_y + 1, index_z + 1, sample_counter, &
                      gridStepX, gridStepY, gridStepZ)

                  enddo ! i = 1, number_of_elements

                  do i = 1, global_sizeX
                    do j = 1, global_sizeY
                      do k = 1, global_sizeZ
                        if (sample_counter(i, j, k) == 0) then
                          ! TODO: Go through all global_velocities_grid and check if the sample_counter is zero
                          ! If yes, try to interpolate from neighbours

                          global_velocities_grid(i, j, k)%px = dble(i - 1) * gridStepX
                          global_velocities_grid(i, j, k)%py = dble(j - 1) * gridStepY
                          global_velocities_grid(i, j, k)%pz = dble(k - 1) * gridStepZ
                        else if (sample_counter(i, j, k) > 1) then
                          factor = dble(sample_counter(i, j, k))
                          vec1 = global_velocities_grid(i, j, k)
                          global_velocities_grid(i, j, k)%vx = vec1%vx / factor
                          global_velocities_grid(i, j, k)%vy = vec1%vy / factor
                          global_velocities_grid(i, j, k)%vz = vec1%vz / factor
                        endif

                        if (global_velocities_grid(i, j, k)%vx < vx_min_global) then
                          global_velocities_grid(i, j, k)%vx = vx_min_global + eps
                        else if (global_velocities_grid(i, j, k)%vx > vx_max_global) then
                          global_velocities_grid(i, j, k)%vx = vx_max_global - eps
                        endif

                        if (global_velocities_grid(i, j, k)%vy < vy_min_global) then
                          global_velocities_grid(i, j, k)%vy = vy_min_global + eps
                        else if (global_velocities_grid(i, j, k)%vy > vy_max_global) then
                          global_velocities_grid(i, j, k)%vy = vy_max_global - eps
                        endif

                        if (global_velocities_grid(i, j, k)%vz < vz_min_global) then
                          global_velocities_grid(i, j, k)%vz = vz_min_global + eps
                        else if (global_velocities_grid(i, j, k)%vz > vz_max_global) then
                          global_velocities_grid(i, j, k)%vz = vz_max_global - eps
                        endif
                      enddo
                    enddo
                  enddo

                  deallocate(sample_counter)
                  call log_message("New velocity calculation end")

                  call log_message("Surface calculation start")
                  ! fill grid with values:
                  !$omp parallel default( shared ) &
                  !$omp private( i, j, k, px, py, pz, vx, vy, vz, sz ) &
                  !$omp private( fill_index, fill_value, lowest_x_velocity )
                  !$omp do schedule( static )
                  do k = 0, global_sizeZ - 1
                      pz = k * gridStepZ ! z values start from 0, since the pecube internal model also start from z = 0
                      do j = 0, global_sizeY - 1
                          py = j * gridStepY
                          lowest_x_velocity = 0.0

                          do i = 0, global_sizeX - 1
                              px = i * gridStepX
                              call get_surface_from_neighbours(px, py, sz)
                              global_surface_grid(i + 1, j + 1) = sz
                          enddo ! i

                      enddo ! j
                  enddo ! k
                  !$omp end do
                  !$omp end parallel
                else
                  ! fill grid with values:
                  !$omp parallel default( shared ) &
                  !$omp private( i, j, k, px, py, pz, vx, vy, vz, sz ) &
                  !$omp private( fill_index, fill_value, lowest_x_velocity )
                  !$omp do schedule( static )
                  do k = 0, global_sizeZ - 1
                      pz = k * gridStepZ ! z values start from 0, since the pecube internal model also start from z = 0
                      do j = 0, global_sizeY - 1
                          py = j * gridStepY
                          lowest_x_velocity = 0.0

                          do i = 0, global_sizeX - 1
                              px = i * gridStepX
                              call get_velo_from_neighbours(px, py, pz, vx, vy, vz)
                              global_velocities_grid(i + 1, j + 1, k + 1) = vec_t(px, py, pz, vx, vy, vz)

                              call get_surface_from_neighbours(px, py, sz)
                              global_surface_grid(i + 1, j + 1) = sz
                          enddo ! i

                      enddo ! j
                  enddo ! k
                  !$omp end do
                  !$omp end parallel
                endif ! use_new_velocity
                call log_message("Surface calculation end")

                ! Write grid to cache file
                call log_message("write result to cache file")
                open(file_unit_cached, file="input/cached_"//trim(filename), &
                    status="new", form="unformatted", access="stream")

                write(file_unit_cached) global_sizeX, global_sizeY, global_sizeZ

                do i = 1, global_sizeX
                    do j = 1, global_sizeY
                        do k = 1, global_sizeZ
                            write(file_unit_cached) global_velocities_grid(i, j, k)
                        enddo
                        write(file_unit_cached) global_surface_grid(i, j)
                    enddo
                enddo
                close(file_unit_cached)

                file_loaded = .TRUE.
            endif
        end subroutine load_velocities

        subroutine calculate_velo(i, index_in_x, index_in_y, index_in_z, sample_counter, &
            gridStepX, gridStepY, gridStepZ)
          implicit none

          integer(4), intent(in) :: i, index_in_x, index_in_y, index_in_z
          integer(4), intent(inout), dimension(:,:,:) :: &
            sample_counter(global_sizeX, global_sizeY, global_sizeZ)
          integer(4) :: index_x, index_y, index_z

          real(8), intent(in) :: gridStepX, gridStepY, gridStepZ
          real(8) :: px, py, pz

          type(vec_t) :: vec1, vec2

          index_x = index_in_x
          index_y = index_in_y
          index_z = index_in_z

          if (index_x <= 0) then
            index_x = 1
          else if (index_x > global_sizeX) then
            index_x = global_sizeX
          endif

          if (index_y <= 0) then
            index_y = 1
          else if (index_y > global_sizeY) then
            index_y = global_sizeY
          endif

          if (index_z <= 0) then
            index_z = 1
          else if (index_z > global_sizeZ) then
            index_z = global_sizeZ
          endif

          px = dble(index_x - 1) * gridStepX
          py = dble(index_y - 1) * gridStepY
          pz = dble(index_z - 1) * gridStepZ

          vec1 = vec_t(px, py, pz, 0.0, 0.0, 0.0)
          call interpolate_velo(global_velocities(i), vec1)

          if (sample_counter(index_x, index_y, index_z) == 0) then
            global_velocities_grid(index_x, index_y, index_z) = &
              vec_t(px, py, pz, vec1%vx, vec1%vy, vec1%vz)

            sample_counter(index_x, index_y, index_z) = 1
          else
            vec2 = global_velocities_grid(index_x, index_y, index_z)
            global_velocities_grid(index_x, index_y, index_z)%vx = vec1%vx + vec2%vx
            global_velocities_grid(index_x, index_y, index_z)%vy = vec1%vy + vec2%vy
            global_velocities_grid(index_x, index_y, index_z)%vz = vec1%vz + vec2%vz

            sample_counter(index_x, index_y, index_z) = sample_counter(index_x, index_y, index_z) + 1
          endif
        end subroutine calculate_velo

        !> Find four nearest neighbours of given point
        !! and set velocity
        !! @param px x coordinate of point
        !! @param py y coordinate of point
        !! @param pz z coordinate of point
        !! @param vx output x velocity
        !! @param vy output y velocity
        !! @param vz output z velocity
        subroutine get_velo_from_neighbours(px, py, pz, vx, vy, vz)
            implicit none

            real(8), intent(in) :: px, py, pz
            real(8), intent(out) :: vx, vy, vz

            real(8) :: d1, d2, d3, d4, d, factor1, factor2, factor3, factor4, factor_sum
            real(8) :: radius
            integer(4) :: i, index1, index2, index3, index4
            type(vec_t) :: p1, p2, p3, p4

            ! radius in km
            radius = 10

            ! set distances to the maximum floating point value
            d1 = huge(d1)
            d2 = huge(d2)
            d3 = huge(d3)
            d4 = huge(d4)

            do i=1,number_of_elements
                d = dist(px, py, pz, global_velocities(i)%px, global_velocities(i)%py, global_velocities(i)%pz)

                if (d < d1) then
                    d4 = d3
                    d3 = d2
                    d2 = d1
                    d1 = d

                    p4 = p3
                    p3 = p2
                    p2 = p1
                    p1 = global_velocities(i)

                    index4 = index3
                    index3 = index2
                    index2 = index1
                    index1 = i
                else if (d < d2) then
                    ! d1 is the distance to the closest point
                    d4 = d3
                    d3 = d2
                    d2 = d

                    p4 = p3
                    p3 = p2
                    p2 = global_velocities(i)

                    index4 = index3
                    index3 = index2
                    index2 = i
                else if (d < d3) then
                    d4 = d3
                    d3 = d

                    p4 = p3
                    p3 = global_velocities(i)

                    index4 = index3
                    index3 = i
                else if (d < d4) then
                    d4 = d

                    p4 = global_velocities(i)

                    index4 = i
                endif
            enddo

            factor1 = exp(-d1/radius)
            factor2 = exp(-d2/radius)
            factor3 = exp(-d3/radius)
            factor4 = exp(-d4/radius)

            factor_sum = factor1 + factor2 + factor3 + factor4

            if (factor_sum < 0.001) then
              vx = 0.0
              vy = 0.0
              vz = 0.0
            else
              vx = (p1%vx*factor1 + p2%vx*factor2 + p3%vx*factor3 + p4%vx*factor4) / factor_sum
              vy = (p1%vy*factor1 + p2%vy*factor2 + p3%vy*factor3 + p4%vy*factor4) / factor_sum
              vz = (p1%vz*factor1 + p2%vz*factor2 + p3%vz*factor3 + p4%vz*factor4) / factor_sum
            endif

            ! Fix undefined points below the surface
            if ((p4%pz > pz) .and. (abs(p4%px - px) < radius) .and. (abs(p4%py - py) < radius) .and. (factor4 < 0.1)) then
              vx = p4%vx
              vy = p4%vy
              vz = p4%vz
            end if

            if (vx /= vx) then
              call log_message("global_velocities.f90: x velocity is NaN")
              call velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
                  d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
            endif

            if (vy /= vy) then
              call log_message("global_velocities.f90: y velocity is NaN")
              call velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
                  d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
            endif

            if (vz /= vz) then
              call log_message("global_velocities.f90: z velocity is NaN")
              call velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
                  d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
            endif

            if (is_velocity_out_of_bounds(vx, vx_min_global, vx_max_global)) then
              call log_message("x velocity out of range, program will stop now")
              call log_message("valid range: " + vx_min_global + ", " + vx_max_global)
              call log_message("global_velocities.f90, pos 1: " + vx)
              call velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
                  d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
            endif

            if (is_velocity_out_of_bounds(vy, vy_min_global, vy_max_global)) then
              call log_message("y velocity out of range, program will stop now")
              call log_message("valid range: " + vy_min_global + ", " + vy_max_global)
              call log_message("global_velocities.f90, pos 2: " + vy)
              call velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
                  d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
            endif

            if (is_velocity_out_of_bounds(vz, vz_min_global, vz_max_global)) then
              call log_message("z velocity out of range, program will stop now")
              call log_message("valid range: " + vz_min_global + ", " + vz_max_global)
              call log_message("global_velocities.f90, pos 3: " + vz)
              call velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
                  d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
            endif
        end subroutine get_velo_from_neighbours

        subroutine interpolate_velo(vec1, vec2)
          implicit none

          type(vec_t), intent(in) :: vec1
          type(vec_t), intent(inout) :: vec2

          real(8) :: distance, factor
          real(8), parameter :: radius = 10 ! radius in [km]

          distance = dist(vec1%px, vec1%py, vec1%pz, vec2%px, vec2%py, vec2%pz)

          factor = exp(-distance / radius)

          vec2%vx = vec1%vx * factor
          vec2%vy = vec1%vy * factor
          vec2%vz = vec1%vz * factor
        end subroutine interpolate_velo

        !> Calculate the norm of the given velocity vector
        !! @param velo velocity vector
        real(8) function normVelo(velo)
            implicit none

            type(vec_t), intent(in) :: velo

            normVelo = sqrt((velo%vx**2) + (velo%vy**2) + (velo%vz**2))
        end function normVelo

        !> Calculate distance from given point
        !! @param x1 x coordinate of point (from geometry)
        !! @param y1 y coordinate of point (from geometry)
        !! @param z1 z coordinate of point (from geometry)
        !! @param x2, y2, z2 current point from global_velocities array
        real(8) function dist(x1, y1, z1, x2, y2, z2)
            implicit none

            real(8), intent(in) :: x1, y1, z1, x2, y2, z2

            dist = sqrt(((x1 - x2)**2) + ((y1 - y2)**2) + ((z1 - z2)**2))
        end function dist

        subroutine get_surface_from_neighbours(px, py, sz)
            real(8), intent(in) :: px, py
            real(8), intent(out) :: sz
            real(8) :: radius, d, factor, factor_sum
            integer(4) :: i

            ! radius in km
            radius = 10
            sz = 0.0
            factor_sum = 0.0

            do i = 1, number_of_surface_points
                d = dist(px, py, 0.0_8, global_surface(i)%x, global_surface(i)%y, 0.0_8)
                factor = exp(-d/radius)
                factor_sum = factor_sum + factor
                sz = sz + (factor * global_surface(i)%z)
            enddo

            sz = sz / factor_sum

        end subroutine get_surface_from_neighbours

        subroutine linear_interpolation(p1, p2, val, p3, axis)
            implicit none

            type(vec_t), intent(in) :: p1, p2
            type(vec_t), intent(out) :: p3

            integer(4), intent(in) :: axis

            real(8), intent(in) :: val
            real(8) :: factor

            if (axis == 1) then ! along x axis
                factor = (val - p1%px) / (p2%px - p1%px)
                p3%px = val
                p3%py = p1%py
                p3%pz = p1%pz
            else if (axis == 2) then ! along y axis
                factor = (val - p1%py) / (p2%py - p1%py)
                p3%px = p1%px
                p3%py = val
                p3%pz = p1%pz
            else if (axis == 3) then ! along z axis
                factor = (val - p1%pz) / (p2%pz - p1%pz)
                p3%px = p1%px
                p3%py = p1%py
                p3%pz = val
            else
                call log_message("global_velocities.f90, linear_interpolation, unknown axis: " + axis)
                error stop 1
            endif

            p3%vx = p1%vx + ((p2%vx - p1%vx) * factor)
            p3%vy = p1%vy + ((p2%vy - p1%vy) * factor)
            p3%vz = p1%vz + ((p2%vz - p1%vz) * factor)
        end subroutine linear_interpolation

        subroutine get_velo_from_grid(px, py, pz, vx, vy, vz)
            implicit none

            real(8), intent(in) :: px, py, pz
            real(8), intent(out) :: vx, vy, vz

            integer(4) :: indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1

            type(vec_t) :: c000, c001, c010, c011, c100, c101, c110, c111
            type(vec_t) :: c00z, c01z, c10z, c11z
            type(vec_t) :: c0yz, c1yz
            type(vec_t) :: cxyz

            real(8) :: sz, pz_corrected

            ! trilinear interpolation.
            ! see http://en.wikipedia.org/wiki/Trilinear_interpolation

            indexX0 = floor(px / global_gridStepX) + 1

            if (indexX0 < 1) then
                vx = 0
                vy = 0
                vz = 0
                return
            else if (indexX0 > global_sizeX - 1) then
                vx = 0
                vy = 0
                vz = 0
                return
            endif

            indexX1 = indexX0 + 1

            indexY0 = floor(py / global_gridStepY) + 1

            if (indexY0 < 1) then
                vx = 0
                vy = 0
                vz = 0
                return
            else if (indexY0 > global_sizeY - 1) then
                vx = 0
                vy = 0
                vz = 0
                return
            endif

            indexY1 = indexY0 + 1

            sz = global_surface_grid(indexX0, indexY0)

            pz_corrected = pz - sz

            indexZ0 = floor((pz_corrected) / global_gridStepZ) + 1

            if (indexZ0 < 1) then
                vx = 0
                vy = 0
                vz = 0
                return
            else if (indexZ0 > global_sizeZ - 1) then
                vx = 0
                vy = 0
                vz = 0
                return
            endif

            indexZ1 = indexZ0 + 1

            c000 = global_velocities_grid(indexX0, indexY0, indexZ0)
            c001 = global_velocities_grid(indexX0, indexY0, indexZ1)
            c010 = global_velocities_grid(indexX0, indexY1, indexZ0)
            c011 = global_velocities_grid(indexX0, indexY1, indexZ1)
            c100 = global_velocities_grid(indexX1, indexY0, indexZ0)
            c101 = global_velocities_grid(indexX1, indexY0, indexZ1)
            c110 = global_velocities_grid(indexX1, indexY1, indexZ0)
            c111 = global_velocities_grid(indexX1, indexY1, indexZ1)

            ! interpolate along z axis
            call linear_interpolation(c000, c001, pz_corrected, c00z, 3)
            call linear_interpolation(c010, c011, pz_corrected, c01z, 3)
            call linear_interpolation(c100, c101, pz_corrected, c10z, 3)
            call linear_interpolation(c110, c111, pz_corrected, c11z, 3)

            ! interpolate along y axis
            call linear_interpolation(c00z, c01z, py, c0yz, 2)
            call linear_interpolation(c10z, c11z, py, c1yz, 2)

            ! interpolate along x axis
            call linear_interpolation(c0yz, c1yz, px, cxyz, 1)

            vx = cxyz%vx
            vy = cxyz%vy
            vz = cxyz%vz

            if (vx /= vx) then
                call velocity_is_nan(vx, vy, vz, px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
                  c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
            endif

            if (vy /= vy) then
              call velocity_is_nan(vx, vy, vz, px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
                c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
            endif

            if (vz /= vz) then
              call velocity_is_nan(vx, vy, vz, px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
                c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
            endif



            if (is_velocity_out_of_bounds(vx, vx_min_global, vx_max_global)) then
              call log_message("x velocity out of range, program will stop now")
              call log_message("valid range: " + vx_min_global + ", " + vx_max_global)
              call log_message("global_velocities.f90, pos 5: " + vx)
              call velocity_out_of_bounds(px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
                  c000, c001, c010, c001, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
            endif

            if (is_velocity_out_of_bounds(vy, vy_min_global, vy_max_global)) then
              call log_message("y velocity out of range, program will stop now")
              call log_message("valid range: " + vy_min_global + ", " + vy_max_global)
              call log_message("global_velocities.f90, pos 6: " + vy)
              call velocity_out_of_bounds(px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
                  c000, c001, c010, c001, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
            endif

            if (is_velocity_out_of_bounds(vz, vz_min_global, vz_max_global)) then
              call log_message("z velocity out of range, program will stop now")
              call log_message("valid range: " + vz_min_global + ", " + vz_max_global)
              call log_message("global_velocities.f90, pos 5: " + vz)
              call velocity_out_of_bounds(px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
                  c000, c001, c010, c001, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
            endif

        end subroutine get_velo_from_grid

        !> Return the velocitiy for given point (back to 'geometry')
        !! @param px x coordinate of point
        !! @param py y coordinate of point
        !! @param pz z coordinate of point
        !! @param vx velocity in x
        !! @param vy velocity in y
        !! @param vz velocity in z
        subroutine get_global_velocity(px, py, pz, vx, vy, vz)
            implicit none

            real(8), intent(in) :: px, py, pz
            real(8), intent(out) :: vx, vy, vz

            if (file_loaded) then
                if (.not. allocated(global_velocities_grid)) then
                    call log_message("global_velocities.f90, get_global_velocity, array 'global_velocities' not allocated!")
                    error stop 1
                endif

                call get_velo_from_grid(px, py, pz, vx, vy, vz)
            else
                ! velocity file not loaded yet, so we just set all values to zero
                call log_message("global_velocities.f90: velocity file not loaded yet, set all to zero.")
                vx = 0.0
                vy = 0.0
                vz = 0.0
            endif
        end subroutine get_global_velocity

        subroutine check_get_global_velocity()
            implicit none

            integer(4), parameter :: file_unit = 83
            integer(4) :: i
            real(8) :: px, py, pz, vx, vy, vz

            open (file_unit, file="debug1.dat", status="new")

            do i=1,number_of_elements
                px = global_velocities(i)%px
                py = global_velocities(i)%py
                pz = global_velocities(i)%pz
                call get_global_velocity(px, py, pz, vx, vy, vz)
                write (file_unit, *) px, py, pz, vx, vy, vz
            enddo

            close(file_unit)

            error stop 1
        end subroutine check_get_global_velocity

        function add_str_vec(str1, vec1)
            implicit none

            character(*), intent(in) :: str1
            type(vec_t), intent(in) :: vec1

            character(len(str1) + (16 * 6)) :: add_str_vec

            write(add_str_vec, "(a, F16.4, F16.4, F16.4, F16.4, F16.4, F16.4)") str1, vec1%px, vec1%py, vec1%pz, vec1%vx, vec1%vy, vec1%vz
        end function add_str_vec

        function is_velocity_out_of_bounds(v, vmin, vmax)
          implicit none

          real(8), intent(in) :: v, vmin, vmax

          logical :: is_velocity_out_of_bounds

          is_velocity_out_of_bounds = (v < vmin) .or. (v > vmax)

        end function is_velocity_out_of_bounds

        subroutine velocity_is_nan(vx, vy, vz, px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
            c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
          implicit none

          real(8), intent(in) :: vx, vy, vz, px, py, pz, sz, pz_corrected
          integer(4), intent(in) :: indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1
          type(vec_t), intent(in) :: c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz

          call log_message("global_velocities.f90: velocity is NaN")
          call log_message("vx: " + vx + ", vy: " + vy + ", vz: " + vz)
          call log_message("px, py, pz, sz, pz_corrected: " + px + ", " + py + ", " + pz + ", " + sz + ", " + pz_corrected)
          call log_message("indexX0, indexX1: " + indexX0 + ", " + indexX1)
          call log_message("indexY0, indexY1: " + indexY0 + ", " + indexY1)
          call log_message("indexZ0, indexZ1: " + indexZ0 + ", " + indexZ1)
          call log_message("c000: " + c000)
          call log_message("c001: " + c001)
          call log_message("c010: " + c010)
          call log_message("c011: " + c011)
          call log_message("c100: " + c100)
          call log_message("c101: " + c101)
          call log_message("c110: " + c110)
          call log_message("c111: " + c111)
          call log_message("c00z: " + c00z)
          call log_message("c01z: " + c01z)
          call log_message("c10z: " + c10z)
          call log_message("c11z: " + c11z)
          call log_message("c0yz: " + c0yz)
          call log_message("c1yz: " + c1yz)
          call log_message("cxyz: " + cxyz)

          error stop 1
        end subroutine velocity_is_nan

        subroutine velocity_out_of_bounds(px, py, pz, sz, pz_corrected, indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1, &
            c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz)
          implicit none

          real(8), intent(in) :: px, py, pz, sz, pz_corrected
          integer(4), intent(in) :: indexX0, indexX1, indexY0, indexY1, indexZ0, indexZ1
          type(vec_t), intent(in) :: c000, c001, c010, c011, c100, c101, c110, c111, c00z, c01z, c10z, c11z, c0yz, c1yz, cxyz

          call log_message("px, py, pz, sz, pz_corrected: " + px + ", " + py + ", " + pz + ", " + sz + ", " + pz_corrected)
          call log_message("indexX0, indexX1: " + indexX0 + ", " + indexX1)
          call log_message("indexY0, indexY1: " + indexY0 + ", " + indexY1)
          call log_message("indexZ0, indexZ1: " + indexZ0 + ", " + indexZ1)
          call log_message("c000: " + c000)
          call log_message("c001: " + c001)
          call log_message("c010: " + c010)
          call log_message("c011: " + c011)
          call log_message("c100: " + c100)
          call log_message("c101: " + c101)
          call log_message("c110: " + c110)
          call log_message("c111: " + c111)
          call log_message("c00z: " + c00z)
          call log_message("c01z: " + c01z)
          call log_message("c10z: " + c10z)
          call log_message("c11z: " + c11z)
          call log_message("c0yz: " + c0yz)
          call log_message("c1yz: " + c1yz)
          call log_message("cxyz: " + cxyz)

          error stop 1
        end subroutine velocity_out_of_bounds

        subroutine velocity_out_of_bounds2(px, py, pz, factor1, factor2, factor3, factor4, &
            d1, d2, d3, d4, p1, p2, p3, p4, index1, index2, index3, index4)
          implicit none

          real(8), intent(in) :: px, py, pz, factor1, factor2, factor3, factor4
          real(8), intent(in) :: d1, d2, d3, d4
          type(vec_t), intent(in) :: p1, p2, p3, p4
          integer(4), intent(in) :: index1, index2, index3, index4

          call log_message("px, py, pz: " + px + ", " + py + ", " + pz)
          call log_message("factor 1, 2, 3, 4: " + factor1 + ", " + factor2 + ", " + factor3 + ", " + factor4)
          call log_message("d 1, 2, 3, 4: " + d1 + ", " + d2 + ", " + d3 + ", " + d4)
          call log_message("index 1, 2, 3, 4: " + index1 + ", " + index2 + ", " + index3 + ", " + index4)
          call log_message("p1: " + p1)
          call log_message("p2: " + p2)
          call log_message("p3: " + p3)
          call log_message("p4: " + p4)
          error stop 1
        end subroutine
end module m_global_velocities
