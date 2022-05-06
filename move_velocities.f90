module m_move_velocities
    use m_compiler
    use m_logger
    use m_pecube_config

    implicit none
    save

    ! make everything private
    private

    integer(4), parameter :: str_len = 300

    !> @struct point_info_t
    !! @brief contains x, z, and point id values from move input file
    type :: point_info_t
        integer(4) :: pt_id
        real(8) :: x
        real(8) :: z
    end type point_info_t

    !> @struct bin_tree_t
    !! @brief contains sorted point ids from velocity file
    type :: bin_tree_t
        integer(4) :: pt_id
        real(8) :: x
        real(8) :: z
        type(bin_tree_t), pointer :: left ! smaller
        type(bin_tree_t), pointer :: right ! greater
    end type bin_tree_t

    !> @struct node_list_t
    !! @brief contains output node information
    type :: node_list_t
        real(8) :: x_out
        real(8) :: z_out
        real(8) :: vx
        real(8) :: vz
        logical :: valid
    end type node_list_t

    !> @struct line_item_t
    !! @brief contains year and move velocity file
    type :: line_item_t
        real(8) :: year
        character(str_len) file_name
    end type line_item_t

    public create_velocities

    contains

        subroutine create_velocities(config, input_file_name, velocity_files, num_of_steps)
            implicit none

            ! openssl rand -hex 8

            type(config_t), intent(in) :: config

            integer(4), parameter :: input_file = 60
            integer(4), parameter :: last_move_file = 61
            integer(4), parameter :: current_move_file = 62
            integer(4), parameter :: output_file = 63
            integer(4), parameter :: top_nodes_file = 64

            character(*), intent(in) :: input_file_name
            character(str_len), dimension(*), intent(out) :: velocity_files
            character(str_len) :: output_file_name, top_nodes_file_name
            character(str_len) :: current_move_file_name, last_move_file_name

            integer(4), intent(in) :: num_of_steps
            integer(4) :: io_status, node_index
            integer(4) :: line_counter
            integer(4) :: number_of_nodes, num_of_samples
            integer(4) :: extension_pos, pt_id2, i, j, k
            integer(4) :: output_nodes_capacity, output_nodes_length

            real(8) :: depth_in_KM, step_y_in_KM
            real(8) :: interp_step_x, xmin, xmax, zmin
            real(8) :: current_time_in_Myrs, last_time_in_Myrs
            real(8) :: time_diff, x1, z1, x2, z2
            real(8) :: x_out, y_out, z_out, value_range
            real(8) :: vx, vz, vx1, vz1, vx2, vz2
            real(8) :: vx_interp, vz_interp, x_factor
            real(8) :: x_interp, z_interp

            real(8), parameter :: FROM_MYRS_TO_YRS = 1e6
            real(8), parameter :: FROM_KM_TO_MM = 1e6
            real(8), parameter :: FROM_MM_TO_KM = 1e-6

            logical :: output_exist
            logical :: top_node_exist

            type(line_item_t), dimension(:), allocatable :: input_items

            type(point_info_t), dimension(:), allocatable :: last_move_file_points
            type(point_info_t), dimension(:), allocatable :: current_move_file_points

            type(bin_tree_t), pointer :: last_move_file_tree, current_move_file_tree

            type(node_list_t), dimension(:), allocatable :: node_list, valid_nodes
            type(node_list_t), dimension(:), allocatable :: output_nodes, tmp_nodes

            type(node_list_t) :: last_node

            call log_message("create_velocities: trying to open file: '" // trim(input_File_name) // "'")

            call open_and_check_file_read(input_file, input_file_name)

            ! The input file must have the following format:

            ! first line: depth in y dimension in km
            ! second line: number of nodes in y dimension
            ! third line: number of samples to take for interpolation
            ! fourth line: step size in x direction for interpolation
            ! fifth line to end of file: time in Myrs and filename
            !
            ! example:
            !
            ! 10.0
            ! 5
            ! 10
            ! 1.0
            ! 50.0 file1.dat
            ! 25.0 file2.dat
            ! 12.56 file3.dat
            ! 7.1 file4.dat
            ! 1.8 file5.dat
            ! 0.0 file6.dat

            read(input_file, *) depth_in_KM
            read(input_file, *) number_of_nodes
            step_y_in_KM = depth_in_KM / (dble(number_of_nodes) - 1)
            read(input_file, *) num_of_samples
            read(input_file, *) interp_step_x

            call log_message("depth_in_KM: " + depth_in_KM + ", number_of_nodes: " + number_of_nodes + ", step_y_in_KM: " + step_y_in_KM)
            call log_message("num_of_samples: " + num_of_samples + ", interp_step_x: " + interp_step_x)

            close(input_file)

            allocate(input_items(num_of_steps))

            call open_and_check_file_read(input_file, input_file_name)

            ! skip header lines
            call skip_lines(input_file, input_file_name, 4)

            do i = 1, num_of_steps
                read(input_file, *, iostat=io_status) input_items(i)%year, input_items(i)%file_name

                call check_read_write(input_file_name, i, io_status)

                if (io_status < 0) then
                    exit ! EOF reached
                end if

                velocity_files(i) = input_items(i)%file_name

                input_items(i)%file_name = "input/"//trim(input_items(i)%file_name)

            end do

            close(input_file)



            if (config%just_velocity) then
              ! Everything is done, we just needed the file names (velocity_files)
              return
            end if



            last_time_in_Myrs = input_items(1)%year
            last_move_file_name = input_items(1)%file_name

            xmin = huge(xmin)
            xmax = -huge(xmax)
            zmin = huge(zmin)

            do i = 1, num_of_steps
                if (input_items(i)%file_name(1:3) == 'Nil') then
                    if (i < num_of_steps) then
                        last_move_file_name = input_items(i + 1)%file_name
                    end if
                    continue
                end if

                current_time_in_Myrs = input_items(i)%year
                current_move_file_name = input_items(i)%file_name

                extension_pos = 0

                do j = 1, str_len
                    if (current_move_file_name(j:j) == '.') then
                        extension_pos = j - 1
                    end if
                end do

                call log_message("processing time: " + current_time_in_Myrs + ", file name: '" // &
                trim(current_move_file_name) // "'")

                ! find position of file name extension if any

                ! call log_message("extension_pos: " + extension_pos)
                ! call log_message("end_of_string: " + end_of_string)

                if (extension_pos > 1) then
                    ! file name with extension
                    output_file_name = current_move_file_name(1:extension_pos)
                    top_nodes_file_name = trim(output_file_name)//"_vel_top.dat"
                    output_file_name = trim(output_file_name)//"_vel_new.dat"
                else
                    ! file name without extension
                    output_file_name = trim(current_move_file_name)
                    top_nodes_file_name = trim(output_file_name)//"_vel_top.dat"
                    output_file_name = trim(output_file_name)//"_vel_new.dat"
                end if

                ! export to pecube main, remove leading "input/"
                ! TODO: this is currently hard coded
                velocity_files(i) = trim(output_file_name(7:300))

                ! call log_message("last_move_file_name: " + last_move_file_name)
                ! call log_message("current_move_file_name: " + current_move_file_name)
                ! call log_message("output_file_name: " + output_file_name)
                ! call log_message("top_nodes_file_name: " + top_nodes_file_name)

                ! TODO: check if files already exist...

                inquire(file=output_file_name, exist=output_exist)
                inquire(file=top_nodes_file_name, exist=top_node_exist)

                if (output_exist .and. top_node_exist .and. config%use_cached_files) then
                    call log_message("these two files already exist, skipping:")
                    call log_message(trim(output_file_name))
                    call log_message(trim(top_nodes_file_name))
                    cycle
                end if


                call open_and_check_file(output_file, output_file_name)

                call open_and_check_file(top_nodes_file, top_nodes_file_name)

                write (output_file, "(A)") "# particle id,x,y,z,vx,vy,vz,color; units: km, mm/year"
                write (top_nodes_file, "(A)") "# x,y,z,vx,vy,vz; units: km, mm/year"

                time_diff = (last_time_in_Myrs - current_time_in_Myrs) * FROM_MYRS_TO_YRS

                call log_message("time_diff: " + time_diff)

                ! read in data from last move file

                call get_num_of_lines_in_file(last_move_file_name, line_counter)

                allocate(last_move_file_points(line_counter))

                ! call log_message("number of lines in previous move file: " + line_counter)

                call read_in_point_data(last_move_file, last_move_file_name, last_move_file_points, &
                                        line_counter, last_move_file_tree)

                ! read in data from current move file

                call get_num_of_lines_in_file(current_move_file_name, line_counter)

                allocate(current_move_file_points(line_counter))

                ! call log_message("number of lines in current move file: " + line_counter)

                call read_in_point_data(current_move_file, current_move_file_name, current_move_file_points, &
                                        line_counter, current_move_file_tree)

                ! process data
                allocate(node_list(line_counter))

                do j = 1, line_counter
                    pt_id2 = current_move_file_points(j)%pt_id
                    x2 = current_move_file_points(j)%x
                    z2 = current_move_file_points(j)%z

                    if (.not. get_x_z_from_tree(last_move_file_tree, pt_id2, x1, z1)) then
                        call log_message("could not find point ID in previous file!")
                        call log_message("file1: '" // trim(last_move_file_name) // "'")
                        call log_message("file2: '" // trim(current_move_file_name) // "'")
                        call log_message("id: " + pt_id2 + ", x: " + x2 + ", z: " + z2)
                        call log_message("")
                        continue
                    end if

                    ! unit of velocity is mm / year

                    vx = 0.0
                    vz = 0.0

                    if (time_diff /= 0.0) then
                        vx = FROM_KM_TO_MM * (x2 - x1) / time_diff
                        vz = FROM_KM_TO_MM * (z2 - z1) / time_diff
                    end if

                    x_out = x2
                    z_out = z2

                    if (x_out < xmin) then
                        xmin = x_out
                    elseif (x_out > xmax) then
                        xmax = x_out
                    end if

                    if (z_out < zmin) then
                        zmin = z_out
                    end if

                    node_list(j)%x_out = x_out
                    node_list(j)%z_out = z_out
                    node_list(j)%vx = vx
                    node_list(j)%vz = vz
                    node_list(j)%valid = .true.

                    do k = 1, number_of_nodes
                        y_out = dble(k * step_y_in_KM)
                        write (output_file, "(I10, F20.8, F20.8, F20.8, F20.8, F20.8, F20.8)") pt_id2, x_out, y_out, z_out, vx, 0.0, vz
                    end do
                end do

                close (output_file)

                ! find top nodes:
                value_range = xmax - xmin

                call log_message("value_range: " + value_range)

                allocate(valid_nodes(num_of_samples + 1))

                do j = 1, num_of_samples + 1
                    valid_nodes(j)%valid = .false.
                end do

                do j = 1, line_counter
                    node_index = int(((node_list(j)%x_out - xmin) / value_range) * num_of_samples) + 1
                    if (.not. valid_nodes(node_index)%valid) then
                        valid_nodes(node_index) = node_list(j)
                    else if (node_list(j)%z_out > valid_nodes(node_index)%z_out) then
                        valid_nodes(node_index) = node_list(j)
                    end if
                end do

                last_node%valid = .false.
                last_node%x_out = 0.0
                last_node%z_out = 0.0
                last_node%vx = 0.0
                last_node%vz = 0.0

                output_nodes_capacity = 2
                output_nodes_length = 0
                allocate(output_nodes(output_nodes_capacity))

                do j = 1, num_of_samples + 1
                    if (valid_nodes(j)%valid) then
                        x2 = valid_nodes(j)%x_out
                        z2 = valid_nodes(j)%z_out
                        vx2 = valid_nodes(j)%vx
                        vz2 = valid_nodes(j)%vz

                        if (last_node%valid) then
                            x1 = last_node%x_out
                            z1 = last_node%z_out
                            vx1 = last_node%vx
                            vz1 = last_node%vz

                            x_interp = x1
                            z_interp = z1
                            vx_interp = vx1
                            vz_interp = vz1

                            do while(x_interp < x2)
                                output_nodes_length = output_nodes_length + 1

                                if (output_nodes_length > output_nodes_capacity) then
                                    !call log_message("output_nodes_length: " + output_nodes_length)
                                    !call log_message("output_nodes_capacity: " + output_nodes_capacity)

                                    allocate(tmp_nodes(output_nodes_capacity))

                                    tmp_nodes = output_nodes

                                    !call log_message("size output_nodes: " + size(output_nodes))

                                    deallocate(output_nodes)

                                    allocate(output_nodes(output_nodes_capacity * 2))

                                    !call log_message("size output_nodes: " + size(output_nodes))

                                    do k = 1, output_nodes_capacity
                                        output_nodes(k) = tmp_nodes(k)
                                    end do

                                    output_nodes_capacity = output_nodes_capacity * 2

                                    !call log_message("new output_nodes_capacity: " + output_nodes_capacity)
                                    deallocate(tmp_nodes)
                                end if

                                output_nodes(output_nodes_length)%x_out = x_interp
                                output_nodes(output_nodes_length)%z_out = z_interp
                                output_nodes(output_nodes_length)%vx = vx_interp
                                output_nodes(output_nodes_length)%vz = vz_interp

                                x_interp = x_interp + interp_step_x
                                x_factor = (x_interp - x1) / (x2 - x1)
                                z_interp = z1 + (x_factor * (z2 - z1))
                                vx_interp = vx1 + (x_factor * (vx2 - vx1))
                                vz_interp = vz1 + (x_factor * (vz2 - vz1))
                            end do

                        end if

                        last_node = valid_nodes(j)
                    end if
                end do

                call log_message("number of output nodes: " + output_nodes_length)
                call log_message("capacity: " + output_nodes_capacity)

                if (output_nodes_length == 0) then
                  call log_message("no output nodes!")
                  error stop 1
                end if

                do j = 1, output_nodes_length
                    x1 = output_nodes(j)%x_out
                    z1 = output_nodes(j)%z_out
                    vx1 = output_nodes(j)%vx
                    vz1 = output_nodes(j)%vz

                    write(top_nodes_file, "(F20.8, F20.8, F20.8, F20.8, F20.8, F20.8)") x1, 0.0, z1, vx1, 0.0, vz1
                end do

                close (top_nodes_file)

                last_move_file_name = current_move_file_name
                last_time_in_Myrs = current_time_in_Myrs

                call clean_up_id_tree(last_move_file_tree)
                call clean_up_id_tree(current_move_file_tree)

                deallocate(last_move_file_points)
                deallocate(current_move_file_points)
                deallocate(node_list)
                deallocate(valid_nodes)
                deallocate(output_nodes)

                call logger_flush()
            end do ! i = 1, num_of_steps

            call log_message("you need these values for pecube:")
            call log_message("xmin: " + xmin + ", zmin: " + zmin)

        end subroutine create_velocities

        subroutine open_and_check_file_read(file_unit, file_name)
            implicit none

            character(*), intent(in) :: file_name

            integer(4), intent(in) :: file_unit
            integer(4) :: io_stat

            open (file_unit, file=trim(file_name), status="old", iostat=io_stat)
            if (io_stat /= 0) then
                call log_message("could not open file: " // file_name)
                call log_message("code: c4a24658a41e8e01, io_stat: " + io_stat)
                error stop 1
            end if
        end subroutine open_and_check_file_read

        subroutine open_and_check_file(file_unit, file_name)
            implicit none

            character(*), intent(in) :: file_name

            integer(4), intent(in) :: file_unit
            integer(4) :: io_stat

            open (file_unit, file=trim(file_name), status="unknown", iostat=io_stat)
            if (io_stat /= 0) then
                call log_message("could not open file: " // file_name)
                call log_message("code: c4a24658a41e8e01, io_stat: " + io_stat)
                error stop 1
            end if
        end subroutine open_and_check_file

        subroutine check_read_write(file_name, line_number, io_stat)
            implicit none

            character(*), intent(in) :: file_name

            integer(4), intent(in) :: line_number
            integer(4), intent(in) :: io_stat

            if (io_stat > 0) then
                call log_message("io error in input file '" // trim(file_name) // "'")
                call log_message("in line: " + line_number)
                call log_message("code: f2d404b50b5ee138, io_stat: " + io_stat)
                error stop 1
            end if
        end subroutine check_read_write

        subroutine get_num_of_lines_in_file(file_name, line_counter)
            implicit none

            character(*), intent(in) :: file_name
            integer(4), intent(out) :: line_counter

            integer(4), parameter :: file_unit = 91
            integer(4) :: io_stat

            call open_and_check_file_read(file_unit, file_name)

            line_counter = 0

            do
                read (file_unit, *, iostat=io_stat)

                call check_read_write(file_name, line_counter, io_stat)

                if (io_stat < 0) then
                    exit
                end if

                line_counter = line_counter + 1
            end do

            close (file_unit)

            ! call log_message("get_num_of_lines_in_file, line_counter: " + line_counter)
        end subroutine get_num_of_lines_in_file

        subroutine skip_lines(file_unit, file_name, num_of_lines_to_skip)
            implicit none

            character(*), intent(in) :: file_name

            integer(4), intent(in) :: file_unit
            integer(4), intent(in) :: num_of_lines_to_skip
            integer(4) :: io_stat, i

            do i = 1, num_of_lines_to_skip
                read (file_unit, *, iostat=io_stat)

                call check_read_write(file_name, i, io_stat)

                if (io_stat < 0) then
                    exit
                end if
            end do
        end subroutine skip_lines

        recursive function check_add_id_tree(tree, pt_id, x, z) result(ret)
            implicit none

            type(bin_tree_t), pointer :: tree

            integer(4), intent(in) :: pt_id

            logical :: ret

            real(8), intent(in) :: x, z

            if (associated(tree)) then
                if (pt_id < tree%pt_id) then
                    ret = check_add_id_tree(tree%left, pt_id, x, z)
                else if (pt_id > tree%pt_id) then
                    ret = check_add_id_tree(tree%right, pt_id, x, z)
                else
                    ret = .true.
                end if
            else
                allocate(tree)
                tree%pt_id = pt_id
                tree%x = x
                tree%z = z
                tree%left => null()
                tree%right => null()
                ret = .false.
            end if
        end function check_add_id_tree

        recursive subroutine clean_up_id_tree(tree)
            implicit none

            type(bin_tree_t), pointer :: tree

            if (associated(tree)) then
                call clean_up_id_tree(tree%left)
                call clean_up_id_tree(tree%right)
                deallocate(tree)
            end if
        end subroutine clean_up_id_tree

        recursive function get_x_z_from_tree(tree, pt_id, x, z) result(ret)
            implicit none

            type(bin_tree_t), pointer :: tree

            integer(4), intent(in) :: pt_id

            logical :: ret

            real(8), intent(out) :: x, z

            ret = .false.

            if (associated(tree)) then
                if (pt_id == tree%pt_id) then
                    x = tree%x
                    z = tree%z
                    ret = .true.
                elseif (pt_id < tree%pt_id) then
                    ret = get_x_z_from_tree(tree%left, pt_id, x, z)
                elseif (pt_id > tree%pt_id) then
                    ret = get_x_z_from_tree(tree%right, pt_id, x, z)
                end if
            else
                ret = .false.
            end if
        end function get_x_z_from_tree

        subroutine read_in_point_data(move_file, move_file_name, move_file_points, line_counter, pt_id_tree)
            implicit none

            integer(4), intent(in) :: move_file, line_counter

            character(*), intent(in) :: move_file_name

            type(point_info_t), dimension(:), intent(out) :: move_file_points

            type(bin_tree_t), pointer :: pt_id_tree

            integer(4) :: i, io_stat, pt_color, pt_id

            real(8) :: x, y, z

            pt_id_tree => null()

            call open_and_check_file_read(move_file, move_file_name)

            do i = 1, line_counter
                read(move_file, *, iostat=io_stat) x, y, z, pt_color, pt_id

                if (check_add_id_tree(pt_id_tree, pt_id, x, z)) then
                    call log_message("point id not unique: " + pt_id)
                    call log_message("values: x: " + x + ", y: " + y + ", z: " + z)
                    call log_message("code: cde4da37a2604ce6")
                    call clean_up_id_tree(pt_id_tree)
                    error stop 1
                end if

                move_file_points(i)%x = x
                move_file_points(i)%z = z
                move_file_points(i)%pt_id = pt_id

                call check_read_write(move_file_name, i, io_stat)
            end do

            close (move_file)
        end subroutine read_in_point_data

end module m_move_velocities
