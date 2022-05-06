module m_find_upstream_points
use m_pdfmaker_for_pecube

contains

      subroutine find_upstream_points (fnme,num_topo_files,xoutlet,youtlet,nstep,nx,ny,dx,dy,&
                 run,age_info,age_type,nsurf,xdepth,ydepth,edot,&
                 xlonmin,xlatmin,&
                 num_basins,pages,perates,nrun)

      ! This subroutine is passed in x, y, and z coordinates and creates a flow path network
      ! The closest grid point to the specified basin outlet point is found and used for
      ! finding all points upstream of the basin outlet
      ! Once all points upstream are found, pdfmaker_for_pecube.f90 is called to calculate
      ! and output the PDF for the specified basin
      !
      ! The UPSTREAMJ and UPSTREAMI arrays hold the x and y indices, respectively, of the points which
      ! flow into the associated point
      ! The surrounding points to a center point appear as follows:
      ! *******************************
      ! ******** D7   D8   D1 *********
      ! ******** D6   C    D2 *********
      ! ******** D5   D4   D3 *********
      ! *******************************
      ! So, the point D7 with center point (C) of (j,i) would be held in the UPSTREAM arrays as
      ! UPSTREAMJ(j,i,7) and UPSTREAMI(j,i,7). If it flows into the center point then the
      ! values of j-1 and i-1 are stored in the arrays, respectively
      ! The point j,i can have at most 8 values in the UPSTREAM arrays because at most 8 points
      ! can surround the center point (ie. (j-1,i-1); (j-1,i); (j+1,i+1); etc)
      ! This creates a flow network that can easily be tracked when given a point on the grid
      ! The subroutine getCatchment recursively finds the upstream points and appends a list
      ! of the x values, y values, ages, and erosion rates

        use m_bivar
        use m_deplss
        use m_getCatchment
        use m_logger
        use m_data_structures

      IMPLICIT NONE

      integer(4) num_topo_files,i,nx,ny,nstep,j,k,m,nsurf,num_basins,iwk
      character fnme(num_topo_files)*300,run*100
      real*8 xoutlet(num_basins),youtlet(num_basins),dx,dy,distance
      real*8 DI,DIS,CAREA
      real(8), allocatable, dimension(:,:) :: Z
      integer age_type(num_basins)
      integer(4), allocatable, dimension(:,:,:) :: UPSTREAMI, UPSTREAMJ
      real*8 x(nx,ny),y(nx,ny),min_distance
      integer xmin,ymin
      real(8), intent(in) :: edot(nstep,nsurf)
      real(8) :: ages(nstep,nsurf)
      real*8 xdepth(nstep,nsurf),ydepth(nstep,nsurf)
      integer md,ncp,num_xchars,num_ychars
      real*8 x2(nx*ny),y2(nx*ny),new_ages(nx*ny)
      real*8 ages2(nsurf),new_edot(nx*ny)
      real*8 pdf_ages(nx*ny),edot2(nsurf),pages(num_basins,nx*ny)
      real*8 new_ages2(nx,ny),new_edot2(nx,ny),pdf_edot(nx*ny),perates(num_basins,nx*ny)
      real*8 ZD,ZDD,CHAN,wk
      integer counter,c,MFP,nrun,ios,l
      character x_basin_char*100,y_basin_char*100,t4*4
      logical contained(nx,ny)
      real*8 xlonmin,xlatmin
      real*8 xstore(nx*ny),ystore(nx*ny),zstore(nx*ny)

      real(8) :: minval_ages, maxval_ages

      type(age_info_t) age_info(nstep)

      md=1
      ncp=4

      ! 2013.07.09, WK: put stuff on heap instead of stack
      allocate(Z(nx, ny))
      allocate(UPSTREAMI(nx, ny, 8))
      allocate(UPSTREAMJ(nx, ny, 8))


      ! Iterates through all time steps
      do i=1,nstep
        if (i.lt.num_topo_files) then
          call log_message("read in file: '"//fnme(i+1)//"'")
          open (203,file='input/'//fnme(i+1),status="old",iostat=ios)              ! Opens the topography file for the time step
          if (ios.ne.0) then
            open (203,file=fnme(i+1),status='old')
          endif
        else
          call log_message("read in file: '"//fnme(num_topo_files)//"'")
          open (203,file='input/'//fnme(num_topo_files),status="old",iostat=ios)
          if (ios.ne.0) then
            open (203,file=fnme(num_topo_files),status='old')
          endif
        endif
        call log_message("read in " + (nx * ny) + " data values for z level")
        call log_message("(" + nx + ", " + ny + ")")
        read (203,*) (Z(:,k),k=ny,1,-1)                                     ! Reads in all elevation values
        close (203)

        CAREA = 9e20                                            ! Maximum cross-grading area (m^2). Hard coded not sure if value is valid for all models
        DI = dx*1000.                                           ! Grid cell size (m)
        DIS = dx*1000.                                          ! Grid cell size (m)
        ZD = -999.                                              ! Value for missing data point in input file
        ZDD = -999.                                             ! Value to be output for missing data point
        MFP = 1                                                 ! Method of catchment area computations: 1=D8 or Rho8, 2=FD8/FRho8 (Multiple drainage
                                                                ! direction method using slope weighting algorithm), 3=Stream tube method (currently not implemented)
        CHAN = 1.0                                              ! Scaling factor to reduce elevations to absolute elevations in meters (typically 1.0 or 0.01)
        UPSTREAMI = 0                                      ! Initialize upstream y value indices to 0
        UPSTREAMJ = 0                                      ! Initialize upstream x value indices to 0

        call log_message("Z(1,1): " + Z(1,1))
        call log_message("Z(nx, ny): " + Z(nx,ny))

        ! Subroutine creates flow path network
        ! Stores x and y indices of surrounding points (max of 8) that flow into
        ! center point

        call log_message("nx: " + nx + ", ny: " + ny + ", CAREA: " + CAREA + ", DI: " + DI)

        call log_message("DIS: " + DIS + ", ZDD: " + ZDD + ", ZD: " + ZD + ", CHAN: " + CHAN)

        call log_message("MFP: " + MFP + ", UPSTREAMI(1, 1, 1): " + UPSTREAMI(1, 1, 1) + &
          ", UPSTREAMI(nx, ny, 8): " + UPSTREAMI(nx, ny, 8))

        call log_message("UPSTREAMJ(1, 1, 1): " + UPSTREAMJ(1, 1, 1) + ", UPSTREAMJ(nx, ny, 8): " + &
          UPSTREAMJ(nx, ny, 8))


        call DEPLSS(ny, nx, nx*ny, CAREA, DI, DIS, ZDD, ZD, CHAN, MFP, &
             Z, UPSTREAMI, UPSTREAMJ)
        !call DEPLSS(ny, & ! i4
        !            nx, & ! i4
        !            nx*ny, & ! i4
        !            CAREA, & ! r8
        !            DI, & ! r8
        !            DIS, & ! r8
        !            ZDD, & ! r8
        !            ZD,& ! r8
        !            CHAN, & ! r8
        !            MFP, & ! i4
        !            Z, & ! r8
        !            UPSTREAMI, & ! i4
        !            UPSTREAMJ) ! i4

        ! Interpolates the current calculated erosion rate distrbution (with skipping factor)
        ! onto the full resolution mesh (with no skipping factor)
        wk = 0.
        iwk = 0
        new_edot = 0.
        edot2 = edot(i,:)

        call idbvip (nsurf,xdepth(i,:),ydepth(i,:),edot2,nx*ny,&
                     x2,y2,new_edot)

        x = 0.
        y = 0.
        x2 = 0.
        y2 = 0.
        counter=1
        do k=1,ny
          do j=1,nx
            x(j,ny-k+1) = xlonmin+dx*(nx-1)*float(j-1)/float(nx-1)
            y(j,ny-k+1) = xlatmin+dy*(ny-1)*float(k-1)/float(ny-1)
            !x(j,k)=xlonmin+(xlonmax-xlonmin)*(x(j,k)-xmin2)/(xmax-xmin2)
            !y(j,k)=xlatmin+(xlatmax-xlatmin)*(y(j,k)-ymin2)/(ymax-ymin2)
            x2(counter) = x(j,ny-k+1)
            y2(counter) = y(j,ny-k+1)
            counter = counter + 1
          enddo
        enddo

        do m=1,num_basins

          do l = 1, nstep
            if (age_type(m).ge.1 .and. age_type(m).le.7) then
              ages(l,:) = age_info(l)%all_ages(:,age_type(m))
            else if (age_type(m).eq.8) then
              ages(l,:) = age_info(l)%all_ages(:,9)
            else if (age_type(m).eq.9) then
              ages(l,:) = age_info(l)%all_ages(:,8)
            else if (age_type(m).eq.10) then
              ages(l,:) = age_info(l)%all_ages(:,10)
            else if (age_type(m).eq.11) then
              ages(l,:) = age_info(l)%all_ages(:,11)
            endif
          enddo

          minval_ages = minval(ages)
          maxval_ages = maxval(ages)

          if (minval_ages == 0.0 .and. maxval_ages == 0.0) then
             call log_message("find_upstream_points.f90: minval(ages): " + minval_ages)
             call log_message("find_upstream_points.f90: maxval(ages): " + maxval_ages)

             call log_message("Invalid values for ages, will exit now")
             stop
          endif

          ! Interpolates the current calculated age distribution set (with skipping factor)
          ! onto the full resolution mesh (with no skipping factor)
          wk = 0.
          iwk = 0
          new_ages = 0.
          ages2 = ages(i,:)

          do j=1,nsurf
             if (ages2(j) /= ages2(j)) then
                call log_message("ages2 is Nan, j: " + j)
                stop
             endif
          enddo

          call log_message("find_upstream_points.f90: minval(ages2): " + minval(ages2))
          call log_message("find_upstream_points.f90: maxval(ages2): " + maxval(ages2))

          call idbvip (nsurf,xdepth(i,:),ydepth(i,:),ages2,nx*ny,&
                       x2,y2,new_ages)

          ! Converts the interpolated 1-D arrays of ages and erosion rates to 2-D
          ! to be used when finding the upstream points since the upstream indices
          ! are stored as (j,i) pairs in UPSTREAMJ and UPSTREAMI
          counter=1
          do j=1,ny
            do k=1,nx
              new_ages2(k,j) = new_ages(counter)
              new_edot2(k,j) = new_edot(counter)
              counter = counter + 1
            enddo
          enddo

          ! Calculates the x and y values for the full resolution Pecube model (skipping factor of 1)
          ! The 2-D arrays are used for ease in obtaining the x and y coordinates when finding
          ! upstream points of basin outlet
          ! The 1-D arrays are used for the interpolations (see below) of the ages and erosion rates
          ! onto the full resolution model
          ! Distance is always >= 0 so initialize to -1 to check if it is the first distance calculation
          ! Finds the closest point on the Pecube grid to the basin outlet specified
          ! This closest point is used for the upstream calculations
          min_distance = 100000.
          xmin = 0
          ymin = 0
          distance = 0.
          do k=1,ny
            do j=1,nx
              distance = sqrt ( ( xoutlet(m) - x(j,ny-k+1) )**2. + ( youtlet(m) - y(j,ny-k+1) )**2. )
              if ( ( distance .lt. min_distance ) ) then
                min_distance = distance
                ! The xmin and ymin variables store the indices of the x and y values of the closest point
                xmin = j
                ymin = ny-k+1
              endif
            enddo
          enddo

          ! If the distance between the specified basin outlet and the closest point is more than
          ! the distance between two nodes in the full resolution mesh, than skip that outlet
          if (min_distance .gt. sqrt(dx**2.+dy**2.)) then
            call log_message('Error: The basin outlet (' + xoutlet(m) + ',' + youtlet(m) + &
                ') is too far from its closest point in the model')
            call log_message('Skipping outlet for timestep ' + i)
          else
            ! xoutlet, youtlet, and t4 are used in the getCatchment.f90 subroutine
            ! for opening a new file for output of the x, y, and z positions
            ! of the basin points
            do j=1,100
              x_basin_char(j:j)=' '
              y_basin_char(j:j)=' '
            enddo
            write (x_basin_char,'(f100.2)') xoutlet(m)
            write (y_basin_char,'(f100.2)') youtlet(m)
            do j=1,100
              if (x_basin_char(j:j).ne.'0' .and. x_basin_char(j:j).ne.'1' .and. x_basin_char(j:j).ne.'2' &
                  .and. x_basin_char(j:j).ne.'3' .and. x_basin_char(j:j).ne.'4' .and. x_basin_char(j:j).ne.'5' &
                  .and. x_basin_char(j:j).ne.'6' .and. x_basin_char(j:j).ne.'7' .and. x_basin_char(j:j).ne.'8' &
                  .and. x_basin_char(j:j).ne.'9' .and. x_basin_char(j:j).ne.'.') num_xchars=j+1
              if (y_basin_char(j:j).ne.'0' .and. y_basin_char(j:j).ne.'1' .and. y_basin_char(j:j).ne.'2' &
                  .and. y_basin_char(j:j).ne.'3' .and. y_basin_char(j:j).ne.'4' .and. y_basin_char(j:j).ne.'5' &
                  .and. y_basin_char(j:j).ne.'6' .and. y_basin_char(j:j).ne.'7' .and. y_basin_char(j:j).ne.'8' &
                  .and. y_basin_char(j:j).ne.'9' .and. y_basin_char(j:j).ne.'.') num_ychars=j+1
            enddo

            write (t4,'(i4)') i

            if (i.lt.1000) t4(1:1)='0'
            if (i.lt.100) t4(1:2)='00'
            if (i.lt.10) t4(1:3)='000'

            ! The variable 'c' keeps track of the number of upstream points
            ! This subroutine compiles a list of all points upstream of the
            ! basin outlet and their cooresponding x values, y values, ages
            ! and erosion rates
            c=1
            call getCatchment(UPSTREAMI,UPSTREAMJ,nx,ny,new_ages2,x,y,Z,&
                              xmin,ymin,pdf_ages,new_edot2,pdf_edot,c,contained,&
                              xstore,ystore,zstore)

            open (105,file=run(1:nrun)//'/Timestep_'//t4//'_Basin_X_'//x_basin_char(num_xchars:100)//'_Y_'//y_basin_char(num_ychars:100)//'.dat',status='unknown')
            write (105,*) 'TITLE = "Pecube User Defined Basin"'
            write (105,*) 'VARIABLES = "x (km)" "y (km)" "z (km)"'
            write (105,*) 'ZONE T = "Basin"'
            write (105,*) 'I=',c,', J=1, K=1, ZONETYPE=Ordered'
            write (105,*) 'DATAPACKING=POINT'
            write (105,*) 'DT=(DOUBLE DOUBLE DOUBLE)'

            do j=1,c
              write (105,'(3f15.5)') xstore(j),ystore(j),zstore(j)
            enddo
            close (105)

            ! Subroutine to calculate and output PDF for the basin
            call pdfmaker_for_pecube(pdf_ages,c,age_type(m),i,-1,run,pdf_edot,&
                                     x_basin_char,y_basin_char,num_xchars,num_ychars,nrun)

            pages(m,:)=pdf_ages
            perates(m,:)=pdf_edot

          endif
        enddo
      enddo

      deallocate(Z)
      deallocate(UPSTREAMI)
      deallocate(UPSTREAMJ)

      return
      end subroutine
end module m_find_upstream_points
