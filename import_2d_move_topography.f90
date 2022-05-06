!> Module for reading in 2d move topography files
        module m_import_2d_move_topography
            use m_logger
            use m_data_structures

            implicit none
            save

            ! make everything private
            private

            !> exported subroutines:
            public fileImport

            contains
                !> Load all coordinate entries from the 2d move file
                !! @param
                subroutine fileImport(currentTimeStep, zValues, numOfCols, numOfRows, stepX, xlonmin)
                    implicit none

                    integer(4), parameter :: file_unit = 82
                    integer(4), intent(in) :: currentTimeStep, numOfCols, numOfRows
                    real(4), dimension(:,:,:), intent(out) :: zValues
                    real(4), intent(in) :: stepX, xlonmin

                    integer(4) :: ioStatus, lineCounter, id, i, j, currentIndex, lastIndex
                    integer(4) :: numOfColumns, lastCharType
                    integer(4), parameter :: UNSET = 0, WHITESPACE = 1, NON_WHITESPACE = 2
                    character(300) :: inputLine
                    character :: currentChar
                    character(4) :: indexStr
                    type(vector3D_t), dimension(:), allocatable :: vector3D_ts
                    type(vector3D_t) :: tmp
                    real(8) :: zValue, modelX
                    logical(1) :: orderChanged

                    lineCounter = 0
                    numOfColumns = 0
                    lastCharType = UNSET

                    !call log_message("xlonmin  [km] " + xlonmin)

                    read(45, "(A)", iostat=ioStatus) inputLine

                    do i = 1,60
                        currentChar = inputLine(i:i)

                        if ((currentChar == ' ').or.(currentChar == '\t').or.(iachar(currentChar) == 9)) then
                            lastCharType = WHITESPACE
                        else
                            if (lastCharType == WHITESPACE) then
                                numOfColumns = numOfColumns + 1
                            endif
                            lastCharType = NON_WHITESPACE
                        endif
                    enddo

                    numOfColumns = numOfColumns + 1
                    !call log_message("numOfColumns: " + numOfColumns)

                    if ((numOfColumns < 3).or.(numOfColumns > 4)) then
                        call log_message("unsupported number of columns")
                        error stop 1
                    endif

                    do
                        if (ioStatus < 0) then
                            exit ! end of file reached
                        endif

                        lineCounter = lineCounter + 1

                        read(45, *, iostat=ioStatus) inputLine
                    enddo

                    !call log_message("lineCounter: " + lineCounter)

                    allocate(vector3D_ts(lineCounter))

                    rewind(45)

                    do i = 1,lineCounter
                        if (numOfColumns == 4) then
                            ! x, y, z, id
                            read(45, *) vector3D_ts(i)%x, vector3D_ts(i)%y, vector3D_ts(i)%z, id
                        else if (numOfColumns == 3) then
                            ! x, y, z
                            read(45, *) vector3D_ts(i)%x, vector3D_ts(i)%y, vector3D_ts(i)%z
                        endif
                    enddo

                    ! bubble sort, for now :-(
                    ! needs to be changed !
                    do
                        orderChanged = .FALSE.
                        do i = 1, lineCounter - 1
                            if (vector3D_ts(i)%x > vector3D_ts(i+1)%x) then
                                orderChanged = .TRUE.
                                tmp = vector3D_ts(i+1)
                                vector3D_ts(i+1) = vector3D_ts(i)
                                vector3D_ts(i) = tmp
                            endif
                        enddo

                        if (.not.orderChanged) then
                            exit
                        endif
                    enddo

                    currentIndex = 2
                    lastIndex = 1

                    write (indexStr,'(I4.4)') currentTimeStep

                    open(file_unit, file="interpol_surface_"//indexStr//".dat", status="unknown")

                    !modelX = vector3D_ts(1)%x - stepX - xlonmin
                    ! Internal model must start at origin specified in input file
                    modelX = xlonmin

                    !call log_message("modelX [km]: " + modelX)
                    !call log_message("vector3D_ts(1)%x [km]: " + vector3D_ts(1)%x)
                    !call log_message("stepX [km]: " + stepX)
                    !call log_message("xlonmin [km]: " + xlonmin)
                    !call log_message("")

                    do i = 1, numOfCols
                        modelX = modelX + stepX

                        do
                            if (modelx > vector3D_ts(currentIndex)%x) then
                                if (currentIndex < lineCounter) then
                                    lastIndex = currentIndex
                                    currentIndex = currentIndex + 1
                                else
                                    ! no more vector3D_ts in array, so we have to use the last one
                                    lastIndex = lineCounter - 1
                                    currentIndex = lineCounter
                                    vector3D_ts(lastIndex)%z = vector3D_ts(currentIndex)%z
                                    exit
                                endif
                            else
                                exit
                            endif
                        enddo

                        ! linear interpolation
                        zValue = vector3D_ts(lastIndex)%z + ((vector3D_ts(currentIndex)%z - vector3D_ts(lastIndex)%z) * (modelX - vector3D_ts(lastIndex)%x) / (vector3D_ts(currentIndex)%x - vector3D_ts(lastIndex)%x))

                        do j = 1, numOfRows ! just copy the values along the y axis
                            zValues(currentTimeStep, i, j) = real(zValue * 1000.0) ! convert from km to m
                        enddo

                        write(file_unit, *) modelX, zValue, vector3D_ts(currentIndex)%z
                    enddo

                    close(file_unit)
                end subroutine fileImport
        end module m_import_2d_move_topography
