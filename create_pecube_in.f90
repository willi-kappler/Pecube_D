module m_create_pecube_in
    contains

! Subroutine that reads in input parameters from Pecube.in, modifies variables,
! and writes variables to scratch file for Pecube.f90 to read in

      subroutine create_pecube_in(config)

      use m_import_2d_move_topography
      use m_logger
      use m_pecube_config

      type(config_t), intent(inout) :: config


! The user can modify this file at will but the output it produces (Pecube.in)
! must obey rules that are described in the user guide and in Pecube.f90

      real,dimension(:,:,:),allocatable::zNZ
      real,dimension(:,:),allocatable::z
      real,dimension(:),allocatable::timek,Peclet,topomag,topooffset,x1f
      real,dimension(:),allocatable::y1f,x2f,y2f,def,dif,zoff
      real,dimension(:),allocatable::Pe2
      real,dimension(:),allocatable::theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,thermflag
      integer,dimension(:),allocatable::nfnme,age_flags,nobsfile
      real,dimension(:),allocatable::zmin,zmax
      real*8,dimension(:),allocatable :: age_obs,dage_obs
      character,dimension(:),allocatable::fnme*300,obsfile*300
      real*4 dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7
      real*8 isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp
      real*8 istc,shtc,lhtc,ghtc,thtc,isrho,shrho,lhrho,ghrho,thrho
      real*8 ishc,shhc,lhhc,ghhc,thhc
      character run*255,line*1024,topo_prefix*50,file_iter*4,det_calc*300
      integer num_topo_files,topo_file_flag,prefix_length,He_flag
      real xlonobs,xlatobs,heightobs
      integer shearint,nrun,ios,node_thresh,num_of_char
      real*8 fric
      character cascadedir*100

      integer(4), dimension(:), allocatable :: geoflag
      real(8) :: yrsec
      real(4) :: dx, dy, xlon, xlat
      real(8), dimension(:), allocatable :: intrusion_start, intrusion_end, intrusion_temperature
      ! temperature hold period, suggested by Byron Adams
      ! Maintain intrusion temperature for some time
      real(8), dimension(:), allocatable :: temperature_hold_period

      character(300) :: move_file, skip_line
      character(300), dimension(:), allocatable :: temperature_files

      ! may be uninitialized
      num_of_char = 0
      xl = 0.0
      yl = 0.0
      prefix_length = 0

      dx = 0.0
      dy = 0.0
      xlon = 0.0
      xlat = 0.0

      open (54,file='Pecube.in',status='old')

      ! remove comment lines and empty lines

      open (55,status='scratch')
    1 read (54,'(a1024)',end=2) line
      if (line(1:1).ne.'$'.and. line(1:1).ne.' ') write (55,'(a)') line
      goto 1
    2 close (54)
      rewind (55)

! skip first lines since they define the pecube run mode

    do
        read(55, "(a)") skip_line
        if (skip_line(1:19) == "pecube_end_run_mode") then
            exit
        endif
    enddo



! run is the name of the run (assumes a directory of that name exists)
! should be 5 character long
      do i=1,100
        run(i:i)=' '
      enddo
      read (55,'(a100)') run
      do i=1,100
        if (run(i:i).ne.' ') nrun=i
      enddo
!         run='output/Pecube-D'

    config%output_folder = trim(run)

! Number of topography files to be loaded (added 2/7/2007 by cspath)
    call log_message("create_pecube_in.f90: read number of topography files")
    read(55,*) num_topo_files

! Reads the topography flag value (1=User lists all files individually,
! 2=User sets prefix and Pecube iterates through all files with that prefix)
    call log_message("create_pecube_in.f90: read topo file flag")
    read(55,*) topo_file_flag

! Checks that the program should be reading in topography file names
    if (num_topo_files .ne. 0) then
          allocate(fnme(num_topo_files),nfnme(num_topo_files))
          if (topo_file_flag == 1 .or. topo_file_flag == 3) then
! Reads in the names of the topography files (added 2/7/2007 by cspath)
            do k=1,num_topo_files
              do i=1,300
                fnme(k)(i:i)=' '
              enddo
              read (55,'(a)') fnme(k)
              do i=1,300
                if (fnme(k)(i:i).ne.' ') then
                  nfnme(k)=i
                else
                  exit
                end if
              enddo
            enddo

          else
! Reads in topography prefix name
            read(55,*) topo_prefix
            do i=1,50
              if (topo_prefix(i:i).ne.' ') prefix_length=i
            enddo
          endif

    else ! num_topo_files == 0, so we just read in NIL
      allocate(fnme(1),nfnme(1))
      do i=1,300
        fnme(1)(i:i)=' '
      enddo
      read (55,'(a)') fnme(1)
      do i=1,300
        if (fnme(1)(i:i).ne.' ') nfnme(1)=i
      enddo
    endif

! nx0 and ny0 are the number of points in the input topographic file
! dx and dy are the data spacings in the x- and y-direction in the input topographic
! file (in degrees or meters)
! nskip is the number of points that are skipped when sampling the initial
! topo file
! xlat and xlon are the latitude/longitude of the bottom left corner of the
! data set (in deg. or meters)
! Reads in coordinate flag (coordflag)
! If coordflag=1, the input type is degrees
! If coordflag=2, the input type is utm
      call log_message("create_pecube_in.f90: read coord flags (Input 5, 6, 7, 8, 9)")
      read (55,*) coordflag
      read (55,*) nx0,ny0
      read (55,*) dx,dy
      read (55,*) nskip
      read (55,*) xlon,xlat

! Convert necessary inputs from meters to kilometers if coordinate system is utm
      if(coordflag.eq.2) then
        dx=dx/1000.
        dy=dy/1000.
        xlon=xlon/1000.
        xlat=xlat/1000.
        call log_message("dx: " + dx + ", dy: " + dy + ", xlon: " + xlon + ", xlat: " + xlat)
       endif


! nstep is the number of time step
! tau is the erosion time scale (assuming an exponential decvrease in topography) (in Myr)
      call log_message("create_pecube_in.f90: read number of time steps (Input 10, 11)")
      read (55,*) nstep,tau

! Compares the number of topography files to the number of time steps
! If there are less topography files than nstep+1, then it assumes the same topography
! for the last time steps
! If there are more topography files, then it quits
! Added 2/6/2007 by cspath
    if ((nstep+1) .gt. num_topo_files) then
      call log_message('Warning: Number of topography files less than number of time steps + 1.')
      call log_message('Assuming same topography for last ' + ((nstep+1)-num_topo_files) + ' files.')
    else if ((nstep+1) .lt. num_topo_files) then
      call log_message('Error: Number of topography files exceeds number of time steps + 1')
      call log_message('Reduce files to or below amount of time steps')
      stop
    endif

! for each time step + 1
! timek is the time (in Myr)
! Peclet is velocity (erosion rate) in mm / year
! topomag is the topographic amplification factor at this step
! topooffset is the topographic offset

      allocate (timek(nstep+1),Peclet(nstep+1),topomag(nstep+1),topooffset(nstep+1))
      allocate (x1f(nstep+1),y1f(nstep+1),x2f(nstep+1),y2f(nstep+1),def(nstep+1),dif(nstep+1))
      allocate (geoflag(nstep+1),theta(nstep+1),phi(nstep+1),thermflag(nstep+1))
      allocate (mft_ratein(nstep+1),mbt_ratein(nstep+1),mct_ratein(nstep+1),stf_ratein(nstep+1))
      allocate (Pe2(nstep+1))
      allocate (intrusion_start(nstep + 1), intrusion_end(nstep + 1), intrusion_temperature(nstep + 1))
      allocate (temperature_hold_period(nstep + 1))

! 2012.07.14, WK: initialize all values with zero

      x1f=0.0
      y1f=0.0
      x2f=0.0
      y2f=0.0
      def=0.0
      dif=0.0
      theta=0.0
      phi=0.0
      mft_ratein=0.0
      mbt_ratein=0.0
      mct_ratein=0.0
      stf_ratein=0.0
      Pe2=0.0

! Reads in time step information from Pecube.in
! Reads in 12 values no matter what geometry is set
! Checks the geometry and then determines which values are needed for later calculations
! Sets all unnecessary variables to 0
        call log_message("create_pecube_in.f90: read time step information (Input 12)")
        do istep = 1, nstep + 1
          !call log_message("reading step: " + istep)
          !if (istep > 1) then
          !    call log_message("timek(istep - 1): " + timek(istep - 1))
          !    call log_message("topomag(istep - 1): " + topomag(istep - 1))
          !    call log_message("topooffset(istep - 1): " + topooffset(istep - 1))
          !    call log_message("thermflag(istep - 1): " + thermflag(istep - 1))
          !    call log_message("geoflag(istep - 1): " + geoflag(istep - 1))
          !endif
          read (55,*) timek(istep), &
            topomag(istep), &
            topooffset(istep), &
            thermflag(istep), &
            geoflag(istep), &
            Peclet(istep), &
            dummy1, &
            dummy2, &
            dummy3, &
            dummy4, &
            dummy5, &
            dummy6, &
            dummy7, &
            intrusion_start(istep), & ! Added new parameter for intrusion
            intrusion_end(istep), &
            intrusion_temperature(istep), &
            temperature_hold_period(istep)

            if (geoflag(istep) .eq. 1) then
              x1f(istep)=0.
              y1f(istep)=0.
              x2f(istep)=0.
              y2f(istep)=0.
              def(istep)=0.
              dif(istep)=0.
              theta(istep)=0.
              phi(istep)=0.
              mft_ratein(istep)=0.
              mbt_ratein(istep)=0.
              mct_ratein(istep)=0.
              stf_ratein(istep)=0.
              Pe2(istep)=0.
            else if (geoflag(istep) .eq. 2) then
              theta(istep)=dummy1
              phi(istep)=dummy2
              x1f(istep)=0.
              y1f(istep)=0.
              x2f(istep)=0.
              y2f(istep)=0.
              def(istep)=0.
              dif(istep)=0.
              mft_ratein(istep)=0.
              mbt_ratein(istep)=0.
              mct_ratein(istep)=0.
              stf_ratein(istep)=0.
              Pe2(istep)=0.
            else if (geoflag(istep) .eq. 3) then
              x1f(istep)=dummy1
              y1f(istep)=dummy2
              x2f(istep)=dummy3
              y2f(istep)=dummy4
              def(istep)=dummy5
              dif(istep)=dummy6
              Pe2(istep)=dummy7
              theta(istep)=0.
              phi(istep)=0.
              mft_ratein(istep)=0.
              mbt_ratein(istep)=0.
              mct_ratein(istep)=0.
              stf_ratein(istep)=0.
            else if (geoflag(istep) .eq. 4) then
              mft_ratein(istep)=dummy1
              mbt_ratein(istep)=dummy2
              mct_ratein(istep)=dummy3
              stf_ratein(istep)=dummy4
              x1f(istep)=0.
              y1f(istep)=0.
              x2f(istep)=0.
              y2f(istep)=0.
              def(istep)=0.
              dif(istep)=0.
              theta(istep)=0.
              phi(istep)=0.
              Pe2(istep)=0.
            else if (geoflag(istep) .eq. 5) then
              x1f(istep)=dummy1
              y1f(istep)=dummy2
              x2f(istep)=dummy3
              y2f(istep)=dummy4
              def(istep)=dummy5
              dif(istep)=dummy6
              theta(istep)=0.
              phi(istep)=0.
              mft_ratein(istep)=0.
              mbt_ratein(istep)=0.
              mct_ratein(istep)=0.
              stf_ratein(istep)=0.
              Pe2(istep)=0.
            else if ((geoflag(istep) == 6) .or. (geoflag(istep) == 7) .or. (geoflag(istep) == 8)) then
              x1f(istep)=dummy1
              y1f(istep)=dummy2
              x2f(istep)=dummy3
              y2f(istep)=dummy4
              def(istep)=dummy5
              dif(istep)=dummy6
              Pe2(istep)=dummy7 ! Pe2 is uplift rate outer ellipse
              theta(istep)=0.
              phi(istep)=0.
              mft_ratein(istep)=0.
              mbt_ratein(istep)=0.
              mct_ratein(istep)=0.
              stf_ratein(istep)=0.
            endif

    enddo ! istep = 1, nstep + 1


! 2012.12.14, WK
! read in the ranges for allowed velocity values inside the velocity field

    call log_message("read in min and max values for velocity files (Input 12.1)")
    read (55, *) vx_min, vx_max
    read (55, *) vy_min, vy_max
    read (55, *) vz_min, vz_max

    call log_message("vx min: " + vx_min + ", vx max: " + vx_max)
    call log_message("vy min: " + vy_min + ", vy max: " + vy_max)
    call log_message("vz min: " + vz_min + ", vz max: " + vz_max)

! converts geological time into model time
        do istep=nstep+1,1,-1
          timek(istep)=timek(1)-timek(istep)
        enddo

! isostasy flag (0 no isostasy, 1 isostasy on)
! ** NOTE: rhoc and rhom are now defined below in the thermal properties section **
! young is the elastic plate young modulus (in Pa)
! poisson is poisson's ratio (dimensionless)
! thickness is the elastic thickness of the plate (in km)
! nxiso and nyiso are the resolutions in the x- and y-directions of the grid on
! which the isostatic (flexural) calculations are performed (including the FFT)
! note that these numbers must be powers of two.
      call log_message("create_pecube_in.f90: read isostasy flag, young modules, etc. (Input 13)")
      read (55,*) isoflag,young,poisson,thickness,nxiso,nyiso

! crustal thickness is the averaged crustal thickness (i.e. the depth at which the
! temperature is assumed to be constant) (in km)
! nzin is the number of points in the z-direction
! tcond is the thermal conductivity (in W/m K)
! heatcap is the heat capacity (in J/kg K)
! rhoc and rhom are the densities for the crust and mantle, respectively (in kg/m3)
! ** NOTE: these values are also used in the isostatic calculations! **
! -REMOVED- diffusivity is the heat diffusivity (in km2/Myr) -REMOVED-
! ** NOTE: Diffusivity is now calculated below, rather than input **
! tmax is the basal temperature (in C)
! tmsl is the temperature at the top of the model (at z=0)
! tlapse is the lapse rate (or change of temperature with height in the atmosphere)
! (in C/km)
! heatproduction is the rate of heat production (in uW/m^3)

      call log_message("create_pecube_in.f90: read crustal thickness, temperature, etc. (Input 14)")
      read (55,*) crustal_thickness,nzin,tcond,heatcap,rhoc,rhom
      read (55,*) tmax,tmsl,tlapse,hpc,efold,hpm,shearint,fric

! Thermal model parameters for the Nepal model geometry
! Each line lists the values below for the Indian shield, Sub-Himalaya
! Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
! Line one is the volumetric heat production [uW/m^3]
! Line two is the thermal conductivity [W/m K]
! Line three is the rock density [kg/m^3]
! Line four is the specific heat capacity [J/kg K]
      call log_message("create_pecube_in.f90: read thermal model parameters (Input 15)")
      read(55,*) ishp,shhp,lhhp,ghhp,thhp
      read(55,*) istc,shtc,lhtc,ghtc,thtc
      read(55,*) isrho,shrho,lhrho,ghrho,thrho
      read(55,*) ishc,shhc,lhhc,ghhc,thhc

! Convert input values to units used in Pecube
      ! crustal thickness (km --> m)
      crustal_thickness = real(crustal_thickness * 1.d3)
      ! number of seconds in 1 My
      yrsec=3.15569259747d7
      ! heat production (uW/m^3-->C/My)
      heatproduction=real(((hpc*1.d-6)/(rhoc*heatcap))*(yrsec*1.e6))
      ishp=((ishp*1.d-6)/(isrho*ishc))*(yrsec*1.e6)                             ! Indian Shield HP
      shhp=((shhp*1.d-6)/(shrho*shhc))*(yrsec*1.e6)                             ! Sub-Himalaya HP
      lhhp=((lhhp*1.d-6)/(lhrho*lhhc))*(yrsec*1.e6)                             ! Lesser Himalaya HP
      ghhp=((ghhp*1.d-6)/(ghrho*ghhc))*(yrsec*1.e6)                             ! Greater Himalaya HP
      thhp=((thhp*1.d-6)/(thrho*thhc))*(yrsec*1.e6)                             ! Tethyan Himalaya HP
      ! thermal diffusivity (km^2/My)
      diffusivity = real((tcond/(rhoc*heatcap))*(1.d-3)**2*yrsec*1.e6)
      isdiff=(istc/(isrho*ishc))*(1.d-3)**2*yrsec*1.e6                          ! Indian Shield diffusivity
      shdiff=(shtc/(shrho*shhc))*(1.d-3)**2*yrsec*1.e6                          ! Sub-Himalaya diffusivity
      lhdiff=(lhtc/(lhrho*lhhc))*(1.d-3)**2*yrsec*1.e6                          ! Lesser Himalaya diffusivity
      ghdiff=(ghtc/(ghrho*ghhc))*(1.d-3)**2*yrsec*1.e6                          ! Greater Himalaya diffusivity
      thdiff=(thtc/(thrho*thhc))*(1.d-3)**2*yrsec*1.e6                          ! Tethyan Himalaya diffusivity

! obsfile is the name of the observation file
      call log_message("create_pecube_in.f90: read observation file (Input 16)")
        read (55,*) num_obsfiles
        allocate (obsfile(num_obsfiles),nobsfile(num_obsfiles))
        do j=1,num_obsfiles
          do i=1,300
            obsfile(j)(i:i)=' '
          enddo
          read (55,'(a)') obsfile(j)
          do i=1,300
            if (obsfile(j)(i:i).ne.' ') nobsfile(j)=i
          enddo
        enddo

! Read in flags for which ages to calculate and output
      allocate (age_flags(config%num_of_age_flags))
      call log_message("create_pecube_in.f90: read age flags (Input 17)")
      read (55, *) age_flags(1:config%num_of_age_flags)

      ! Type of detrital calculation, if any
      ! Options are for no PDF calculation, Cascade catchments, and user specified catchments for non-cascade models
      call log_message("create_pecube_in.f90: read det calc (Input 18)")
      read (55,*) det_calc

      ! Minimum number of nodes for a Cascade model catchment to be output
      call log_message("create_pecube_in.f90: read node_thresh (Input 19)")
      read (55,*) node_thresh

      ! Directory where the Cascade output files are located for use in Pecube
      ! Only necessary if option '1' is selected for det_calc
      do i=1,100
         cascadedir(i:i)=' '
      enddo

      call log_message("create_pecube_in.f90: read cascade dir (Input 20)")
      read (55,'(a100)') cascadedir

      ! WK, 2014.03.05, read in temperature file name
      call log_message("create_pecube_in.f90: read in temperature files (Input 21)")

      allocate(temperature_files(nstep + 1))

      do i = 1, nstep + 1
        read (55, '(A)') temperature_files(i)
        if (temperature_files(1)(1:3) == "Nil") then
          exit ! No temperature file given
        endif
        call log_message("create_pecube_in.f90: temperature file from input file: " + temperature_files(i))
      enddo

      call log_message("create_pecube_in.f90: read in 2D move file (Input 22)")
      read(55, '(A)') move_file

      close (55)

! ----------
! WK: end of config file 'Pecube.in'
! ----------


! Opens all topography files and stores the values in 3D array (zNZ)
! Modified for multiple topography file inputs on 2/7/2007 by cspath

      call log_message("create_pecube_in.f90: read all topography files")

    allocate (zNZ(nstep+1,nx0,ny0))

    zNZ = 0.0

    if (topo_file_flag.eq.1) then ! all file names are listed
      if (num_topo_files.eq.0) then
            zNZ=0.0
      else
        do k=1, nstep+1
          if(k.le.num_topo_files) then
              open (45,file='input/'//fnme(k)(1:nfnme(k)),status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=fnme(k)(1:nfnme(k)),status='old')
              endif
              call log_message("create_pecube_in.f90: read topo file: "//fnme(k)(1:nfnme(k)))
              do m=1,ny0
                do n=1,nx0
                  read (45,*) zNZ(k,n,m)
                enddo
              enddo
          else
          open (45,file='input/'//fnme(num_topo_files)(1:nfnme(num_topo_files)),status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=fnme(num_topo_files)(1:nfnme(num_topo_files)),status='old')
              endif
              call log_message("create_pecube_in.f90: read topo file: "//fnme(num_topo_files)(1:nfnme(num_topo_files)))
              do m=1,ny0
                do n=1,nx0
                  read (45,*) zNZ(k,n,m)
                enddo
              enddo
          endif
              close (45)
        enddo
      endif
    elseif (topo_file_flag.eq.2) then ! only prefix is given
      if (num_topo_files.eq.0) then
        zNZ=0.0
      else
        do k=0, nstep
          if (k.le.(num_topo_files-1)) then
          write (file_iter,'(i4)') k
          if (k.lt.1000) file_iter(1:1)='0'
              if (k.lt.100) file_iter(1:2)='00'
              if (k.lt.10) file_iter(1:3)='000'
              open (45,file='input/'//topo_prefix(1:prefix_length)//file_iter//'.dat',status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=topo_prefix(1:prefix_length)//file_iter//'.dat',status='old')
              endif
              fnme(k+1)=topo_prefix(1:prefix_length)//file_iter//'.dat'
              read (45,*) zNZ(k+1,:,:)
          else
              open (45,file='input/'//topo_prefix(1:prefix_length)//file_iter//'.dat',status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=topo_prefix(1:prefix_length)//file_iter//'.dat',status='old')
              endif
              read (45,*) zNZ(k+1,:,:)
          endif
          close (45)
        enddo
      endif
    ! 2011.09.02, WK read in data file exported from 2d move
    elseif (topo_file_flag == 3) then ! use topo file exported from 2d move
      if (num_topo_files.eq.0) then
        zNZ=0.0
      else
        do k=1, nstep + 1
            if (k <= num_topo_files) then
                open (45,file='input/'//fnme(k)(1:nfnme(k)),status='old',iostat=ios)
                if (ios.ne.0) then ! look in current directory, if file not found in "input/"
                    open (45,file=fnme(k)(1:nfnme(k)),status='old')
                    call log_message("2D move topography file: " // fnme(k)(1:nfnme(k)))
                else
                    call log_message("2D move topography file: input/" // fnme(k)(1:nfnme(k)))
                endif

                ! import_2d_move_topography.f90
                call fileImport(k, zNZ, nx0, ny0, dx, xlon)

                close (45)
            else ! just copy the vales from the last time step
                do m = 1, ny0
                    do n = 1, nx0
                        zNz(k,n,m) = zNz(num_topo_files,n,m)
                    enddo
                enddo
            endif
        enddo

        ! write all values to file for debugging
        open (45, file="znz.dat", status="unknown")

        do m = 1, ny0
            do n = 1, nx0
                write (45, *) n, m, zNz(1, n, m)
            enddo
        enddo

        close(45)
!        stop

      endif
    endif


! call log_message("finish")
! stop


! Stores certain values of the topography files input based on nx, ny, nskip
! Modified for multiple topography file input on 2/7/2007 by cspath
      nx=(nx0-1)/nskip+1
      ny=(ny0-1)/nskip+1
      allocate (z(nstep+1,nx*ny),zoff(nx*ny))

      do k=1, nstep+1
        ij=0
        do j=1,ny0,nskip
          do i=1,nx0,nskip
            ij=ij+1
            z(k, ij) = zNZ(k,i,j)
            ! call log_message("Indices for skipping factor: x:" + i + ", y:" + j)
            ! if (z(k, ij) /= z(k, ij)) then
            !   call log_message("create_pecube_in.f90, z(k, ij) is Nan, k: " + k + ", ij: " + ij)
            !   call log_message("zNZ(k,i,j): " + zNZ(k, i, j) + ", i: " + i + ", j: " + j)
            !   stop
            ! endif
          enddo
        enddo
      enddo

! Modified for multiple topography inputs on 2/7/2007 by cspath
      allocate(zmin(nstep+1),zmax(nstep+1))
      do k=1,nstep+1
        zmin(k)=minval(z(k,:))
        zmax(k)=maxval(z(k,:))
      enddo


! Calculates xl,yl based on coordinate flag (1=Degrees, 2=UTM) (added by cspath)
      if (coordflag .eq. 1) then
          xl=dx*(nx-1)*nskip*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
          yl=dy*(ny-1)*nskip*111.11
      else if (coordflag .eq. 2) then
          xl=dx*(nx-1)*nskip
          yl=dy*(ny-1)*nskip
      else
          call log_message('Coordinate flag is set to an invalid value')
      endif

      zl=crustal_thickness/1.e3

! Modified for multiple topography inputs on 2/7/2007 by cspath
    !call log_message('z1:' + z(1,1))
      do k=1,nstep+1
        z(k,:)=(z(k,:)-zmin(k))/crustal_thickness*zl

        ! do j = 1, nx*ny
        !   if (z(k,j) /= z(k,j)) then
        !     call log_message("create_pecube_in.f90,  z(k,j) is Nan, k: " + k + ", j: " + j)
        !     call log_message("zmin(k): " + zmin(k) + ", crustal_thickness: " + crustal_thickness + ", zl: " + zl)
        !     stop
        !   endif
        ! enddo

      enddo
    !call log_message('z2:' + z(1,1))

      ! Commented out because this was leading to incorrect heat
      !   production values.
      ! dwhipp - 10/07
      !heatproduction=heatproduction*diffusivity/zl**2*tmax

! Converts values to kilometers based on coordinate flag (added by cspath)
! Note: If coordinate system is utm, the values inputted here should be in km as specified in Pecube.in
      if (coordflag .eq. 1) then
        do istep=1,nstep+1
            x1f(istep)=(x1f(istep)-xlon)*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
            y1f(istep)=(y1f(istep)-xlat)*111.11
            x2f(istep)=(x2f(istep)-xlon)*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
            y2f(istep)=(y2f(istep)-xlat)*111.11
        enddo
      else if (coordflag .eq. 2) then
        do istep=1,nstep+1
            x1f(istep)=x1f(istep)-xlon
            y1f(istep)=y1f(istep)-xlat
            x2f(istep)=x2f(istep)-xlon
            y2f(istep)=y2f(istep)-xlat
        enddo
      endif

! Modified for multiple topography inputs on 2/7/2007 by cspath
      do k=1,nstep+1
          zmax(k)=maxval(z(k,:))
      enddo

      open (7, file = trim(config%output_folder) // "/Pecube.dat", status = "unknown")
      !open (7,status='scratch')

      ilog=0
      iterative=1
      interpol=1

      nxs=nx/2
! 2011.07.25, WK: non standard format specifier
      write (7,'(a,i10)') run,nrun
      write (7,*) num_topo_files,det_calc

      do i=1,config%num_of_age_flags
          write (7,*) age_flags(i)
      end do

      write (7,'(A)') (fnme(k),k=1,num_topo_files)
      config%num_elements_surf = (nx-1)*(ny-1)
      config%nsurf = nx * ny
      config%npe = 4
      write (7,*) config%npe,config%nsurf,nzin,config%num_elements_surf,zl,diffusivity,heatproduction,efold,shearint,fric
      ! Write out Nepal model geometry info
      ! dwhipp 11/07
      write (7,*) isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,lhhp
      write (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol
      write (7,*) isoflag,tau,rhoc,rhom,heatcap
      write (7,*) nx,ny,nxiso,nyiso,nx0,ny0,dx,dy
      write (7,*) xl/(nx-1)*1.d3,yl/(ny-1)*1.d3,young,poisson,thickness*1.d3
      ! xlonmin, xlonmax, xlatmin, xlatmax
      write (7,*) xlon,xlon+dx*(nx-1)*nskip,xlat,xlat+dy*(ny-1)*nskip
      ! 2012.12.14, WK: ranges for velocity files
      write (7,*) vx_min, vy_min, vz_min
      write (7,*) vx_max, vy_max, vz_max

    write (7,'(A)') trim(move_file)

        do j=1,ny
          do i=1,nx
            x=xl*float(i-1)/float(nx-1)
            y=yl*float(j-1)/float(ny-1)
            write (7,*) x,y
          enddo
        enddo

        do j=1,ny-1
          do i=1,nx-1
            icon1=(j-1)*nx+i
            icon2=icon1+1
            icon3=icon1+nx+1
            icon4=icon1+nx
            write (7,*) icon1,icon2,icon3,icon4
          enddo
        enddo

! Modified for multiple topography inputs on 2/7/2007 by cspath
      do k=1,num_topo_files
        zmax(k)=maxval(z(k,:))
      enddo

! Modified for multiple topography inputs on 2/7/2007 by cspath
! Writes all time step information to scratch file to be read by Pecube.f90
    do istep = 0, nstep
        write (7, *) timek(istep + 1), Peclet(istep + 1), 0, &
                    x1f(istep + 1), y1f(istep + 1), x2f(istep + 1), y2f(istep + 1), &
                    def(istep + 1), dif(istep + 1), thermflag(istep + 1), geoflag(istep + 1), theta(istep + 1), phi(istep + 1),&
                    mft_ratein(istep + 1), mbt_ratein(istep + 1) ,mct_ratein(istep + 1), stf_ratein(istep + 1), Pe2(istep + 1)
        write (7, *) intrusion_start(istep + 1), intrusion_end(istep + 1), intrusion_temperature(istep + 1)
        write (7, *) temperature_hold_period(istep + 1)
        write (7, *) (z(istep+1,k)*topomag(istep+1)+topooffset(istep+1),k=1,nx*ny)
        if (temperature_files(1)(1:3) == "Nil") then
          write (7, *) "Nil"
        else
          write (7, *) trim(temperature_files(istep + 1))
        endif

        ! do k = 1, nx*ny
        !   if (z(istep+1,k) /= z(istep+1,k)) then
        !     call log_message("create_pecube_in.f90, z(istep+1,k) is Nan, istep: " + istep + ", k: " + k)
        !     stop
        !   endif
        ! enddo

        ! if (topomag(istep+1) /= topomag(istep+1)) then
        !   call log_message("create_pecube_in.f90,  topomag(istep+1) is Nan, istep: " + istep)
        !   stop
        ! endif

        ! if (topooffset(istep+1) /= topooffset(istep+1)) then
        !   call log_message("create_pecube_in.f90, topooffset(istep+1) is Nan, istep: " + istep)
        !   stop
        ! endif

    enddo

    if (det_calc .eq. '1') write (7,*) node_thresh
    if (det_calc .eq. '1') write (7,*) cascadedir


    config%nx = nx
    config%ny = ny
    config%spacing_long = dx
    config%spacing_lat = dy
    config%nskip = nskip
    config%location_long = xlon
    config%location_lat = xlat

! observations
      write (7,*) num_obsfiles
      allocate (age_obs(config%num_of_age_flags),dage_obs(config%num_of_age_flags))
      age_obs = -1.0
      dage_obs = -1.0

      if (num_obsfiles.eq.0) then
        nobs=0
        write (7,*) nobs
      else if (obsfile(1)(1:nobsfile(1)).eq.'Nil') then
        nobs=0
        write (7,*) nobs
      else
        do j = 1, num_obsfiles
          open (8,file='input/'//obsfile(j)(1:nobsfile(j)),status='old',iostat=ios)
          if (ios.ne.0) then
            open (8,file=obsfile(j)(1:nobsfile(j)),status='old')
          endif

          call log_message("create_pecube_in.f90: read obs file: "//obsfile(j)(1:nobsfile(j)))

          call log_message("reading in " + (config%num_of_age_flags - 10 + 4) + &
              " entrys per line (number of age flags: " + config%num_of_age_flags + ")")

          read (8,*) nobs
          write (7,*) nobs
          do i = 1, nobs

            !call log_message("obs: read line " + i + " of " + nobs)
            read (8,*)  xlonobs,xlatobs,heightobs,He_flag,&
                        age_obs(1),dage_obs(1),&
                        age_obs(9),dage_obs(9),&
                        age_obs(8),dage_obs(8),&
                        (age_obs(k),dage_obs(k), k = 10, config%num_of_age_flags)

            do k = 1, config%num_of_age_flags
              if ((dage_obs(k) == 0.0) .and. (age_obs(k) > 0.0)) then
                call log_message("error: dage_obs is zero, index: " + k + ", age_obs: " + age_obs(k))
                call log_message("xlonobs: " + xlonobs + ", xlatobs: " + xlatobs + ", heightobs: " + heightobs)
                call log_message("line: " + i + " of " + nobs)
                call log_message("will exit now")
                stop
              endif
            enddo

            age_obs(He_flag) = age_obs(1)
            dage_obs(He_flag) = dage_obs(1)

            i1=int((xlonobs-xlon)/(dx*nskip))+1
            if (i1.eq.nx) i1=nx-1

            j1=int((xlatobs-xlat)/(dy*nskip))+1
            if (j1.eq.ny) j1=ny-1

            ieobs=i1+(j1-1)*(nx-1)

            if ((ieobs < 1) .or. (ieobs > config%num_elements_surf)) then
                call log_message("i1: " + i1 + ", j1: " + j1 + ", num_elements_surf: " + config%num_elements_surf)
                call log_message("ieobs: " + ieobs + ", nx: " + nx + ", ny: " + ny)
                call log_message("xlonobs: " + xlonobs + ", xlon: " + xlon)
                call log_message("xlatobs: " + xlatobs + ", xlat: " + xlat)

                call log_message("error in file: "//obsfile(j)(1:nobsfile(j)))
                call log_message("ensure that xlonobs > xlon and xlatobs > xlat")

                stop
            end if

            r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip)
            r=-1.+2.*r
            s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip)
            s=-1.+2.*s
            wobs1=(1.-r)*(1.-s)/4.
            wobs2=(1.+r)*(1.-s)/4.
            wobs3=(1.+r)*(1.+s)/4.
            wobs4=(1.-r)*(1.+s)/4.

            write (7,*) He_flag
            write (7,*) xlonobs,xlatobs,heightobs,wobs1,wobs2,wobs3,wobs4,ieobs,&
                        age_obs(He_flag),dage_obs(He_flag)

            do k = 8,config%num_of_age_flags
                write(7, *) age_obs(k), dage_obs(k)
            end do

          enddo
          close(8)
        enddo
      endif

      close(7)

      deallocate (zNZ,z,zoff,zmin,zmax)
      deallocate (timek,Peclet,topomag,topooffset)
      deallocate (x1f,x2f,y1f,y2f,def,dif)
      deallocate (geoflag,theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,thermflag)
      deallocate (Pe2, intrusion_start, intrusion_end, intrusion_temperature)
      deallocate (temperature_hold_period)
      deallocate (fnme,nfnme)
      deallocate (age_flags)
      deallocate (age_obs,dage_obs)
      deallocate (obsfile,nobsfile)
      deallocate (temperature_files)

      return
      end subroutine create_pecube_in
end module m_create_pecube_in
