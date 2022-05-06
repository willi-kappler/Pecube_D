!> Module for reading in temperature from input files
    module m_global_temperature
        use m_compiler
        use m_logger

        implicit none
        save


        ! make everything private
        private

        real(8), dimension(:), allocatable :: x, y, z, t

        public load_temperatures

        contains
            subroutine load_temperatures(nnode, xp, yp, zp, tp, temperature_file)
                implicit none

                integer(4), intent(in) :: nnode

                real(8), intent(in) :: xp(nnode), yp(nnode), zp(nnode)
                real(8), intent(out) :: tp(nnode)

                character(300), intent(in) :: temperature_file
                character(300) :: input_line

                integer(4) :: i, num_of_char, io_status, line_counter

                tp = 0.0

                call log_message("global_temperature.f90: temperature_file: " // trim(temperature_file))

                num_of_char = 1

                do i=1,300
                    if (temperature_file(i:i) == ' ') then
                        num_of_char = i
                        exit
                    endif
                enddo

                open(60,file="input/"//temperature_file(1:num_of_char),status='old',iostat=io_status)

                if (io_status /= 0) then
                    call log_message("could not open temperature file: input/"//temperature_file(1:num_of_char))
                    error stop 1
                endif

                ! count the number of lines in file.
                ! this gives us the number of entries in that file and we
                ! can allocate our array
                line_counter = 0

                ! the first line just contains the header, so we can read it in and forget it
                read(60, *) input_line

                do
                    read(60, *, iostat=io_status) input_line
                    if (io_status < 0) then
                        exit ! end of file reached
                    endif

                    line_counter = line_counter + 1
                enddo

                call log_message("global_temperature.f90, load_temperatures, number of lines: " + line_counter)
                call log_message("global_temperature.f90, load_temperatures, number of nodes: " + nnode)

                allocate(x(line_counter), y(line_counter), z(line_counter), t(line_counter))

                ! start rading input file form the beginning
                rewind(60)

                ! the first line just contains the header, so we can read it in and forget it
                read(60, *) input_line

                do i = 1, line_counter
                    read(60, *) x(i), y(i), z(i), t(i)
                enddo

                close(60)

                ! interpolate onto the pecube grid

                call log_message("start interpolation...")

                do i = 1, nnode
                    call interpolate(xp(i), yp(i), zp(i), tp(i), line_counter)
                    ! call log_message(i, xp(i), yp(i), zp(i), tp(i))
                enddo

                call log_message("interpolation finished")

                deallocate(x, y, z, t)

                ! call log_message("global_temperature.f90:")
                ! call log_message(tp)

            end subroutine load_temperatures

            subroutine interpolate(xp, yp, zp, tp, number_of_elements)
                implicit none

                integer(4), intent(in) :: number_of_elements
                real(8), intent(in) :: xp, yp, zp
                real(8), intent(out) :: tp

                real(8) :: d1, d2, d3, d4, d, factor1, factor2, factor3, factor4, factor_sum
                real(8) :: radius
                integer(4) :: i, index1, index2, index3, index4

                ! radius in km
                radius = 10

                ! set distances to the maximum floating point value
                d1 = huge(d1)
                d2 = huge(d2)
                d3 = huge(d3)
                d4 = huge(d4)

                index1 = 1
                index2 = 1
                index3 = 1
                index4 = 1

                ! find the four nearest neighbour
                do i = 1, number_of_elements
                    d = sqrt(((x(i) - xp)**2) + ((y(i) - yp)**2) + ((z(i) - zp)**2))
                    if (d < d1) then
                        d4 = d3
                        d3 = d2
                        d2 = d1
                        d1 = d

                        index4 = index3
                        index3 = index2
                        index2 = index1
                        index1 = i
                    else if (d < d2) then
                        d4 = d3
                        d3 = d2
                        d2 = d

                        index4 = index3
                        index3 = index2
                        index2 = i
                    else if (d < d3) then
                        d4 = d3
                        d3 = d

                        index4 = index3
                        index3 = i
                    else if (d < d4) then
                        d4 = d

                        index4 = i
                    endif
                enddo

                ! multiplication factor with exponential fall-off
                factor1 = exp(-d1/radius)
                factor2 = exp(-d2/radius)
                factor3 = exp(-d3/radius)
                factor4 = exp(-d4/radius)

                factor_sum = factor1 + factor2 + factor3 + factor4

                tp = (t(index1)*factor1 + t(index2)*factor2 + t(index3)*factor3 + t(index4)*factor4) / factor_sum

            end subroutine interpolate

    end module m_global_temperature
