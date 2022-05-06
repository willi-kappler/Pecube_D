module m_pecube_func
    use m_data_structures

    implicit none

    private

    integer(4) :: nx0, ny0, nstep, nsurf
    integer(4) :: num_threads, num_topo_files
    integer(4) :: num_obsfiles, mc_num_of_erosion_rates
    integer(4), dimension(:), allocatable :: age_flags

    real(8) :: xlonmin, xlatmin, xstep, ystep, error_iter_misfit, last_peclet
    real(8), dimension(:), allocatable :: xsurf, ysurf, zsurf
    real(8), dimension(:), allocatable :: mc_random_erosion_rate
    real(8), dimension(:,:), allocatable :: edot_store

    character(300), dimension(:), allocatable :: topo_file_name

    type(age_info_t), dimension(:), allocatable :: surface_age_info

    ! visible outside this function:
    public pecube_func

    public nx0, ny0, nstep, nsurf
    public num_threads, num_topo_files
    public num_obsfiles

    public age_flags

    public xsurf, ysurf, zsurf
    public xlonmin, xlatmin, xstep, ystep

    public topo_file_name, error_iter_misfit

    public mc_random_erosion_rate, mc_num_of_erosion_rates

    public surface_age_info

    ! public edot_store

    public last_peclet

    contains

    subroutine pecube_func(config)
      use m_compiler
      use omp_lib
      use m_global_velocities
      use m_make_matrix
      use m_find_upstream_points
      use m_interpolate
      use m_ages_header
      use m_bivar
      use m_catchments_output
      use m_create_pecube_in
      use m_find_dt
      use m_find_neighbours
      use m_solve_iterative
      use m_tec_mat_output
      use m_erates
      use m_isostatic_rebound
      use m_pdfmaker_for_data
      use m_detrital_mc
      use m_global_temperature
      use m_move_velocities
      use m_pecube_config
      use m_calculate_misfit
      use m_dynamic_thermal_conductivity
      use m_borehole_ages
      use m_backtrack_temperature
      use m_move_points
      use m_calculate_ages
      use m_export_surface_line

      implicit none

      type(config_t), intent(inout) :: config

! WK: data structure for position and temperature history
!        type :: type_pos_t_hist
!            real(8) :: x, y, z, t
!        end type type_pos_t_hist

! WK: array holding some sample points for position and temperature history
!        type(type_pos_t_hist), dimension(:), allocatable :: pos_t_hist

! WK: constant holding the number of samples to take
!        integer(4), parameter :: num_of_samples = 1000

! WK: variable for random number routine:
!        real(8) :: random_pos

! WK: variable for length (x), width (y) and depth (z)
!        real(8) :: length, width, depth

      type(system_clock_info_t), dimension(9) :: clock_info

      real(8) :: totalTime

! WK: statistics to find out how often tracked particles are outside the model
        integer(4) :: stat1, stat2, stat3


      real*8,dimension(:),allocatable :: x,y,z,xp,yp,zp,dummy,dummyp,header_info
      real*8,dimension(:),allocatable :: t,tp,f
      integer,dimension(:,:),allocatable :: icon
      real*8,dimension(:,:,:),allocatable :: ael
      real*8,dimension(:,:),allocatable :: bel
      integer, dimension(:),allocatable :: kfix,ielsurf
      real*8,dimension(:),allocatable :: zsurfp
      real*8,dimension(:),allocatable :: topoa,topob,rsurf,rebound
      integer,dimension(:,:),allocatable :: iconsurf,neighbour
      real*8,dimension(:,:),allocatable :: xdepth,ydepth,zdepth
      real(8), dimension(:), allocatable :: edot,bg_edot
      real(8), dimension(:), allocatable :: topo_edot
      real*8 mft_ratein,mbt_ratein,mct_ratein,stf_ratein,tb_velo
      real*8 iout,PecletN,PecletO
      real*8 isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp
      real*8 cur_depth,xy_mean,fact,zmin,zmax
      real*8 scale_fact,zmin_orig
      integer ie,je,num_therm_his,niter,plane_store1,plane_store2
      integer,dimension(:),allocatable :: jrec
      integer::counter,nz,nzin,num_left
      character run*255,det_calc*300,comp*3
      logical depth_flag
      real*8,dimension(:),allocatable :: age_obs,dage_obs
      real*8,dimension(:),allocatable :: xmisfit_tot,age_prd,sum_xmisfit
      real*8,dimension(:,:),allocatable :: xmisfit
      integer He_flag,node_thresh
      real*8 dx,dy,x_basin,y_basin
      real*8,dimension(:),allocatable :: pdf_ages,error
      real*8,dimension(:),allocatable :: x_basin1,y_basin1
      real*8,dimension(:),allocatable :: x_basin2,y_basin2
      character,dimension(:),allocatable :: data_file1*300
      character data_file*300,data_comp*3,cascadedir*100
      character(300) move_file, temperature_file
      character(300), dimension(:), allocatable :: velocity_files
      integer,dimension(:),allocatable :: age_type1,age_type2
      integer dataCount,nilCount,age_type,foundtwo,shearint,nrun,ios
      real*8,dimension(:,:),allocatable :: pages,perates
      character,dimension(:),allocatable :: data_comp1*3,data_comp2*3
      real*8 friction
      real*8 Peclet2





      logical :: has_dynamic_thermal_conductivity, file_exists

      real(8) :: swap_z, ftime
      real(8) :: vx_min, vy_min, vz_min
      real(8) :: vx_max, vy_max, vz_max
      real(8) :: diffusivity, diffusivityO
      real(8) :: area_width, area_height, pecube_width, pecube_height
      real(8) :: ftldmean1, ftldmean2, ftldmean3, latitude_step, longitude_step

      real(8) :: alpha, def, dif, dt, efold, eps, ftimep
      real(8) :: heatproduction, hei, heightobs
      real(8) :: peclet, phi, poisson, rhoc, rhom, heatcap, tau, thickness, time, timesurf
      real(8) :: tlapse, tmax, tmsl, tsurf, wobs1, wobs2, wobs3, wobs4
      real(8) :: x1, x10, x1f, x2, x2f, x3, x4, x5, x6, x7, x8, x9, tfinal, theta
      real(8) :: timesurfp, xlatmax, xlatobs, xlonmax, xlonobs
      real(8) :: xmax, xmin, xxx, y1f, y2f, ymax, ymin, young, yyy
      real(8) :: zh, zl, zsurf1, zsurf2, thermflag
      real(8) :: intrusion_start, intrusion_end, intrusion_temperature
      real(8) :: prev_intrusion_start, prev_intrusion_end, prev_intrusion_temperature
      real(8) :: temperature_hold_period, manual_dt, actual_time

      real(8), dimension(:), allocatable :: surf_longitude, surf_latitude

      integer(4) :: geoflag, index1, index2, in

      integer(4) :: i, i1, i2, i3, i4, i5, i6, ic, ieobs, iesurf, ij, ilog, inp
      integer(4) :: interpol, iobs, irec, isoflag, istatic, istep, iterative
      integer(4) :: itime, j, k, kstep, last_val, m, mpe, n, nelem
      integer(4) :: nelemsurf, kk, nnode, nobs, npe, ntime
      integer(4) :: num_basins, number, nx, ny, nxiso, nyiso
      integer(4) :: observ_x, observ_y, p1, p2, p3, p4, number_of_borehole_ages
      integer(4) :: start_step, current_step

      integer(4), dimension(:), allocatable :: dummy_values, therm_his_val
      integer(4), dimension(:), allocatable :: sub_step_sum, num_of_sub_steps

      integer(4), parameter :: file_unit_closure = 81
      integer(4), parameter :: file_unit_length_dist = 82
      integer(4), parameter :: file_unit_surface = 83
      integer(4), parameter :: file_unit_temperature = 84
      integer(4), parameter :: file_unit_temperature_history = 85
      integer(4), parameter :: file_unit_velocity_info = 87
      integer(4), parameter :: file_unit_borehole_velocity_info = 88
      integer(4), parameter :: file_unit_scratch1 = 89

      type(thermal_conductivity_t), dimension(:), allocatable :: dynamic_thermal_conductivity
      type(thermal_conductivity_t) :: current_dynamic_thermal_conductivity

      type(vector3D_t), dimension(:), allocatable :: borehole_ages_points
      type(vector3D_t), dimension(:), allocatable :: surface_pos
      type(vector3D_t), dimension(:), allocatable :: model_pos, model_velo

      type(velocity_info_t), dimension(:), allocatable :: velo_info
      type(age_info_t), dimension(:), allocatable :: borehole_age_info

      character(4) :: file_id1
      character(512) :: file1, file2



! may be not initialized:
    counter       = 0
    ftime         = 0.0_8
    zmin_orig     = 0.0_8
    xy_mean       = 0.0_8
    PecletO       = 0.0_8
    PecletN       = 0.0_8
    fact          = 0.0_8
    cur_depth     = 0.0_8
    diffusivity   = 0.0_8
    diffusivityO  = 0.0_8

    in            = 0
    plane_store1  = 0
    plane_store2  = 0
    tfinal        = 0

    eps           = tiny(eps)

    stat1         = 0
    stat2         = 0
    stat3         = 0

    prev_intrusion_start       = 0.0
    prev_intrusion_end         = 0.0
    prev_intrusion_temperature = 0.0

    ! call cpu_time (times(1))
    call system_clock(clock_info(1)%sys_count, clock_info(1)%sys_count_rate, clock_info(1)%sys_count_max)


! Pecube is a Finite Element solver of the 3D, transient heat transfer
! equation that allows for conduction, vertical advection and production of
! heat. Pecube also allows for the surface geometry to vary with time.

! This version is an improvement (we hope) of the original version
! available by download from the author's website in that it reads an
! input (topography) file as well as a set of cooling ages (He and FT in
! apatite of known location and calculates synthetic ages at the same
! locations for comparison.

! A flexural isostatic model has also been included to
! calculate the effect of isostasy in amplifying the exhumation caused
! by erosional unloading

! The tectonic uplift has been replaced by movement along a listric
! thrust fault

! To understand how this version works, read the information in the input
! file named Pecube.in

! Note that this input file is read in the subroutine create_pecube_in which
! generates another (scratch) file defined as unit 7 from which Pecube reads
! in the controlling data/information. Note that if you know what you are doing
! you can by pass this operation and create your own input file to be read directly
! by Pecube

! WK, 2013.11.13
! set number of age parameters
      config%num_of_age_flags = 26

      allocate(dummy_values(config%num_of_age_flags))

      !stop

      !call log_message("Test cases")
      !call all_test_cases

      call log_message('Reading input')
      call logger_flush()

      call create_pecube_in(config)

      call log_message("Use new velocity interpolation: " + config%use_new_velocity)

!      call cpu_time (times(2))
      call system_clock(clock_info(2)%sys_count, clock_info(2)%sys_count_rate, clock_info(2)%sys_count_max)
      call logger_flush()

! opens input files

      open (7, file = trim(config%output_folder) // "/Pecube.dat", status = "unknown")

! read in general information

! first line:
! run: string that will determine the name of the folder where the input file
!      will be copied to and where the output files (Pecube.out and Pecube.ptt) will be stored.
! num_topo_files: number of topography files to read

      if (allocated(age_flags)) then
          deallocate(age_flags)
      endif
      allocate (age_flags(config%num_of_age_flags))

! 2011.07.25, WK: format specifier not correct
      read (7,'(a,i10)') run, nrun
      read (7,*) num_topo_files,det_calc       ! modified 01/08

      call log_message("run: " + run(1:nrun))
      call log_message("nrun: " + nrun)
      call log_message("output_folder: " + trim(config%output_folder))

      if (nrun == 0) then
        call log_message("Output folder is empty, exit now")
        error stop 1
      endif

      ! Create output folder if it does not exist
      call sys_command("mkdir -p " // trim(config%output_folder))

      do i = 1, config%num_of_age_flags
          read (7,*) age_flags(i)
      enddo

      if (allocated(topo_file_name)) then
          deallocate(topo_file_name)
      endif

      allocate (topo_file_name(num_topo_files))
      read (7,'(A)') (topo_file_name(k),k=1,num_topo_files)

! second line
! npe: number of nodes per surface (2D) elements (3 = triangular elements - 4 = rectangular elements)
! nsurf: number of nodes defining the surface geometry
! nzin: number of nodes in the vertical direction
! nelemsurf: number of surface 2D elements
! zl: thickness of crustal layer (= depth at which t is fixed) (in km), total model thickness
! diffusivity: thermal diffusivity (in km^2/Myr)
! heatproduction: heat production (in degC/Myr)
! efold: e-folding depth for decrease in heat production
! Note: HP is constant above msl and decreases exponentially below
! 1/e decrease in HP at efold

      read (7,*) npe,nsurf,nzin,nelemsurf,zl,diffusivity,heatproduction,efold,&
                 shearint,friction

! Read in Nepal model geometry info
! This line lists the thermal diffusivity and heat production values for the
! Indian shield, Sub-Himalaya Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
      read (7,*) isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp
      mpe=2*npe

      allocate (iconsurf(npe,nelemsurf),ielsurf(nsurf))
      allocate (neighbour(npe,nelemsurf))

      if (allocated(xsurf)) then
        ! Deallocate from previous run
        deallocate(xsurf)
      endif

      if (allocated(ysurf)) then
        ! Deallocate from previous run
        deallocate(ysurf)
      endif

      allocate (xsurf(nsurf), ysurf(nsurf), zsurf(nsurf),zsurfp(nsurf))
      allocate (topoa(nsurf),topob(nsurf),rsurf(nsurf))

      allocate (surf_longitude(nsurf))
      allocate (surf_latitude(nsurf))

! Initialize all memory with zero:

xsurf = 0.0
ysurf = 0.0
zsurf = 0.0
zsurfp = 0.0
topoa = 0.0
topob = 0.0
rsurf = 0.0
surf_longitude = 0.0
surf_latitude = 0.0


! second line:
! tmax: basal temperature (at z=-zl) (in degC)
! tmsl: temperature at mean sea level (z=0) (in degC)
! tlapse: lapse rate (in degC/km)
! nstep: number of stages (the surface topo and exhumatoin rate can be specified at each time stage)
! ilog: redundant (dont use)
! iterative: iterative solver method (1 = Gauss Siedel with overrelaxation - 2 = Conjugate gradient)
! interpol: interpolation used (1 = linear - 2 = quadratic, not recommended)

! isoflag: flag for isostasy (0 no isostasy; 1 isostasy)
! tau: erosion time scale which determines the exponential rate qt which topography
! is transformed from one step to the other
! rhoc crustal density
! rhom mantle density
! nx, ny, nxiso,nyiso are the x and y discretization of the FE grid and the isostatic grid
! note that we have lost a bit of generality in Pecube when implementing the isostatic response
! in that we need (for computational efficiency) the FE grid to be rectangular
! note also that nxiso and nyiso need to be powers of 2
! xstep, ystep are the grid spacing in meter
! young and poisson are young modulus and poisson's ratio
! thickness is effective elastic thickness

      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol
      read (7,*) isoflag,tau,rhoc,rhom,heatcap
      read (7,*) nx,ny,nxiso,nyiso,nx0,ny0,dx,dy
      read (7,*) xstep,ystep,young,poisson,thickness
      if (ilog.eq.1) open (9,file='Pecube.log',status='unknown')

      call log_message("xstep: " + xstep + ", ystep: " + ystep + ", dx: " + dx + ", dy: " + dy)

! Note: the *_store arrays currently have a hard coded 1st column length of 500
! This is because this length is determined by how many sub-time intervals Pecube
! calculates for every time step
      allocate(xdepth(nstep,nsurf),ydepth(nstep,nsurf),zdepth(nstep,nsurf))
      allocate(therm_his_val(nstep+1))

      if (allocated(edot_store)) then
        ! Deallocate from previous error iteration run
        deallocate(edot_store)
      endif

      allocate(edot_store(nstep, nsurf))

      allocate(velocity_files(nstep))
      allocate(sub_step_sum(nstep), num_of_sub_steps(nstep))

! Initialize:
      xdepth           = 0.0
      ydepth           = 0.0
      zdepth           = 0.0

      sub_step_sum     = 0
      num_of_sub_steps = 0

! read in nodal geometry
! xlonmin, xlonmax, xlatmin, xlatmax are all in [km]
      read (7, *) xlonmin, xlonmax, xlatmin, xlatmax

! 2012.12.14, WK:
! read in minimum and maximum allowed velocity values for the velocity field:

      read (7,*) vx_min, vy_min, vz_min
      read (7,*) vx_max, vy_max, vz_max

      call log_message("v ranges: " + vx_min + ", " + vy_min + ", " + vz_min + ", " + vx_max + ", " + vy_max + ", " + vz_max)

      read (7, '(A)') move_file
      if (move_file(1:3) == 'Nil') then
          do i = 1, nstep
              velocity_files(i) = 'Nil'
          enddo
      else
          ! TODO: Fix hard coded input folder
          call create_velocities(config, "input/"//trim(move_file), velocity_files, nstep)
      endif

! xsurf, ysurf: x- and y-locations of the surface nodes
! iconsurf: connectivity matrix describing the 2D surface elements

      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

! WK, 2015.02.19, Monte Carlo

    if (config%mc_use_monte_carlo) then
        if (.not. allocated(mc_random_erosion_rate)) then
            if (allocated(config%mc_time_slices)) then
              mc_num_of_erosion_rates = size(config%mc_time_slices)
            else
              mc_num_of_erosion_rates = nstep
            endif


            allocate(mc_random_erosion_rate(mc_num_of_erosion_rates))
            call random_number(mc_random_erosion_rate)

            call log_message("setting random erosion rate values for monte carlo")

            ! initialize all erosion rates with random numbers
            do i = 1, mc_num_of_erosion_rates
                if (i <= mc_max_number_of_erosion_rates) then
                  if (config%mc_fixed_erosion_rates(i) >= 0.0) then
                    mc_random_erosion_rate(i) = config%mc_fixed_erosion_rates(i)
                  else
                    mc_random_erosion_rate(i) = mc_random_erosion_rate(i) * config%mc_max_erosion_rate
                  endif
                else
                    mc_random_erosion_rate(i) = mc_random_erosion_rate(i) * config%mc_max_erosion_rate
                endif
            enddo

            call log_message("pecube_func.f90: random erosion rate values: " + mc_random_erosion_rate)
        endif

        ! Delete old *.dat files from previous monte carlo run
        call execute_command_line("rm " // run(1:nrun) // "/*.dat")
    endif


! Check for dynamic thermal conductivity
    call load_thermal_conductivity_file(config%thermal_conductivity_file, dynamic_thermal_conductivity, nstep, zl, rhoc, heatcap)
    has_dynamic_thermal_conductivity = allocated(dynamic_thermal_conductivity)

    if (has_dynamic_thermal_conductivity) then
      current_dynamic_thermal_conductivity%num_of_layers = dynamic_thermal_conductivity(1)%num_of_layers
      allocate(current_dynamic_thermal_conductivity%layers(dynamic_thermal_conductivity(1)%num_of_layers))
    end if

! Check for borehole age calculation
    number_of_borehole_ages = load_borehole_ages(config%borehole_ages_file, xlonmin, xlatmin, zl, borehole_ages_points)

    if (number_of_borehole_ages > 0) then
      call log_message("number of borehole ages points: " + number_of_borehole_ages)
      allocate(model_pos(number_of_borehole_ages), model_velo(number_of_borehole_ages))

      model_pos = vector3D_t(0.0, 0.0, 0.0)
      model_velo = vector3D_t(0.0, 0.0, 0.0)
    endif


! WK, 2012.02.15, support for OpenMP
      call omp_set_num_threads(num_threads)

      call log_message('Number of threads: ' + num_threads)
      call log_message('Advecting rocks')

      xmin = minval(xsurf)
      xmax = maxval(xsurf)
      ymin = minval(ysurf)
      ymax = maxval(ysurf)

      allocate(velo_info(nstep))

      do i = 1, nstep
        allocate(velo_info(i)%vz(nsurf))
      enddo










    do istep = nstep, 1, -1
      timesurfp = 0.
      rewind (7)
      read (7,*)

      read (7,*) i1,det_calc,dummy_values

      read (7,'(A)') (topo_file_name(k),k=1,num_topo_files)
      read (7,*) i1,i2,i3,i4,x1,x2,x3,x4,x5,x6
      read (7,*) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
      read (7,*) x1,x2,x3,i1,i2,i3,i4
      read (7,*) x1,x2,x3,x4
      read (7,*) i1,i2,i3,i4,i5,i6,x1,x2
      read (7,*) x1,x2,x3,x4,x5
      read (7,*) x1,x2,x3,x4
      read (7,*) vx_min, vy_min, vz_min
      read (7,*) vx_max, vy_max, vz_max

! WK, 2014.03.05: read in temperature file name
      read (7, '(A)') move_file

      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)

      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)
      ! WK: read in first istep-1 velocity lines and discard them
      ! Peclet is erosion rate (mm/year or km/mil year)
      ! Peclet2 is uplift rate of outer ellipse
      do kstep=0,istep-1
          read (7, *) timesurfp,Peclet,iout,x1f,y1f,x2f,y2f,def,dif,thermflag,geoflag, &
                   theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2
          read (7, *) intrusion_start, intrusion_end, intrusion_temperature, temperature_hold_period
          read (7, *) (zsurfp(i),i=1,nsurf)
          read (7, *) temperature_file
      enddo
      ! WK: we need velocity information for step istep
      read (7,*) timesurf,Peclet,iout,x1f,y1f,x2f,y2f,def,dif,thermflag,geoflag,&
                 theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2
      read (7, *) intrusion_start, intrusion_end, intrusion_temperature, temperature_hold_period
      read (7, *) (zsurf(i),i=1,nsurf)
      read (7, *) temperature_file


      call log_message("pecube_func.f90: current temperature file: " + trim(temperature_file))

      if ((istep == nstep) .and. (config%error_iter_age_dec > 0.0)) then
        call log_message("Reducing time, timesurf: " + timesurf)
        timesurf = timesurf - config%error_iter_age_dec
        call log_message("New reduced time, timesurf: " + timesurf)
      endif

      if (timesurf < timesurfp) then
        call log_message("Time values must be decreasing, timesurf < timesurfp: " + timesurf + " < " + timesurfp)
        stop 1
      endif

      !call log_message("zsurfp: " + zsurfp(1) + ", " + zsurfp(2) + ", " + zsurfp(3))
      !call log_message("zsurf: " + zsurf(1) + ", " + zsurf(2) + ", " + zsurf(3))

      ! do i = 1, nsurf
      !   if (zsurf(i) /= zsurf(i)) then
      !     call log_message("pecube_func.f90, zsurf(i) is Nan, i: " + i + ", istep: " + istep)
      !   endif
      ! enddo

! Calculate proper z-node spacing if nz is zero
! cspath and dwhipp 11/07
      if (nzin == 0) then                                                     ! If number of input z levels is zero, then code will calculate nz and spacing
        zmin       = minval(zsurf)                                            !   of z levels.  z levels will be spaced 1: 1 with the input x and y spacing
        zmin_orig  = zmin                                                     !   down to ~5 km below the model surface, 3:     1 down to ~15 km below the surface
        zmax       = maxval(zsurf)                                            !   and ~9:                                  1 for the rest of the model.
        xy_mean    = (xstep+ystep)/2                                          ! This first portion of the variable z spacing code determines the number of
        cur_depth  = zmin+zl                                                  !   z levels for the new geometry.
        nz         = 1
        depth_flag = .true.
        do while (cur_depth.gt.0.)                                            ! Work down from min elevation to base determining number of z levels needed
          nz=nz+1                                                             ! While still in this loop, increment nz
          if (cur_depth.gt.(zmin+zl)-5.) then                                 ! If within 5 km of model top surface, space node planes at xy_mean (1:1)
            cur_depth=cur_depth-(xy_mean/1000.)                               ! Subtract off new node plane spacing from remaining depth range
            plane_store1=nz-1                                                 ! Store number of planes used at 1:1 spacing
          else if (cur_depth.gt.zmin+zl-15) then                              ! If within 15 km of model surface, space node planes at 3*xy_mean (3:1)
            cur_depth=cur_depth-3*(xy_mean/1000.)                             ! Subtract off new node plane spacing from remaining depth range
            plane_store2=nz-1                                                 ! Store number of planes used at 3:1 spacing
          else                                                                ! If greater than 15 km from model surface, use ~9:1 node plane spacing
            if (depth_flag) then                                              ! Calculate node plane spacing on first click through this condition
              depth_flag=.false.                                              ! Set depth_flag to false to avoid repeating this calculation
              num_left=int((cur_depth)/(9*xy_mean/1000.))+1                   ! Number of remaining node planes is equal to the remaining model depth over
              fact=cur_depth/(num_left*(xy_mean/1000.))                       !   the 9:1 spacing increment, plus one.  This yields spacing of <=9:1.
            endif
            if (int(cur_depth-fact*(xy_mean/1000.)) == 0.) then               ! If near the base of the model, set the cur_depth to zero
              cur_depth=0.
            else                                                              ! Subtract off new node spacing from remaining depth range
              cur_depth=cur_depth-fact*(xy_mean/1000.)
            endif
          endif
        enddo
      else if (nzin.gt.0) then                                                ! If number of input z levels is positive, then use that number for nz
        nz=nzin
      else                                                                    ! Stop if input nz value is negative
        call log_message('Error in Pecube.in: nz must be zero or a positive integer')
        stop
      endif ! nzin

    !call log_message("zmin: " + zmin + ", zmax: " + zmax)

! Stores the thermal flag values into array to be used to determine
! what thermal histories to write out at the end
! Moved so thermal history flag values are stored before time stepping (cspath 10/07)
      therm_his_val(istep) = int(thermflag)
      call log_message("pecube_func.f90, thermflag: " + thermflag)

! isostatic rebound
      topoa = zsurfp
      topob = zsurf

      if (isoflag.eq.1) then
        call isostatic_rebound (topoa,topob,nsurf,rsurf, &
                              rhoc,rhom,nx,ny,nxiso,nyiso, &
                              xstep,ystep,young,poisson,thickness)
      else
        rsurf=0.
      endif

      if (istep >= nstep) tfinal = timesurf


    ! monte carlo simulation
    ! change Peclet here

      if (config%mc_use_monte_carlo) then
          if (allocated(config%mc_time_slices)) then
            actual_time = tfinal - timesurf

            do i = 1, size(config%mc_time_slices)
              if (actual_time >= config%mc_time_slices(i)%s_begin .and. actual_time <= config%mc_time_slices(i)%s_end) then
                Peclet = mc_random_erosion_rate(i)
              endif
            enddo
          else
            if (istep >= nstep) then
              Peclet = mc_random_erosion_rate(nstep)
            else
              Peclet = mc_random_erosion_rate(istep + 1)
            endif
          endif
      endif
    !call log_message("pecube_func.f90: Peclet after mc: " + Peclet)

    ! Use correct Peclet number and diffusivity for Nepal model geometry
    ! Added by dwhipp 11/07
      if (geoflag == 4) then
        ! DAVE: Should include uplift in velo!
        tb_velo      = (mft_ratein+mbt_ratein+mct_ratein-stf_ratein)
        diffusivityO = diffusivity
        diffusivity  = max(isdiff,shdiff,lhdiff,ghdiff,thdiff)
        if (20.-tb_velo.ge.tb_velo) then                                      ! Max velocity is underthrusting in Nepal model geometry
          PecletN = 20.-tb_velo
        else                                                                  ! Max velocity is overthrusting in Nepal model geometry
          PecletN = tb_velo
        endif
        PecletO = Peclet
        Peclet  = PecletN
      endif

      !call log_message("before find_dt, dt: " + dt + ", ntime: " + ntime + ", istep: " + istep)
      call find_dt(zl,diffusivity,nsurf,zsurf,zsurfp, &
                   Peclet,timesurf,timesurfp,istep,eps, &
                   dt,ntime,istatic)
      !call log_message("after find_dt, dt: " + dt + ", ntime: " + ntime)

      !call log_message("Manual dt: " + manual_dt)

      if (manual_dt > 0.0) then
        call log_message("Use manual dt for this time step: " + dt + " -> " + manual_dt)
        dt = manual_dt
        call log_message("ntime old: " + ntime)
        ntime = int((timesurf - timesurfp) / dt) + 1
        call log_message("ntime new: " + ntime)
      endif

      ! Reset Peclet for Nepal model geometry
      ! dwhipp 11/07
      if (geoflag == 4) then
        Peclet      = PecletO
        diffusivity = diffusivityO
      endif


      ! WK: store position of surface nodes for later use
      do i = 1, nsurf
        xdepth(istep, i) = xsurf(i)
        ydepth(istep, i) = ysurf(i)
        zdepth(istep, i) = zsurf(i) + zl ! Adds the model thickness to the extra elevation

        ! if (xdepth(istep, i) /= xdepth(istep, i)) then
        !   call log_message("pecube_func.f90, xdepth(istep, i) is NaN, istep: " + istep + ", i: " + i)
        ! endif
        !
        ! if (ydepth(istep, i) /= ydepth(istep, i)) then
        !   call log_message("pecube_func.f90, ydepth(istep, i) is NaN, istep: " + istep + ", i: " + i)
        ! endif
        !
        ! if (zdepth(istep, i) /= zdepth(istep, i)) then
        !   call log_message("pecube_func.f90, zdepth(istep, i) is NaN, istep: " + istep + ", i: " + i)
        !   call log_message("zsurf(i): " + zsurf(i) + ", zl: " + zl)
        ! endif

      enddo

      ! WK: store velocity information for each time step to use it later
      velo_info(istep)%zl = zl
      velo_info(istep)%x1f = x1f
      velo_info(istep)%y1f = y1f
      velo_info(istep)%x2f = x2f
      velo_info(istep)%y2f = y2f
      velo_info(istep)%def = def
      velo_info(istep)%dif = dif
      velo_info(istep)%theta = theta
      velo_info(istep)%phi = phi
      velo_info(istep)%mft_ratein = mft_ratein
      velo_info(istep)%mbt_ratein = mbt_ratein
      velo_info(istep)%mct_ratein = mct_ratein
      velo_info(istep)%stf_ratein = stf_ratein
      velo_info(istep)%Peclet = Peclet
      velo_info(istep)%Peclet2 = Peclet2
      velo_info(istep)%dt = dt
      velo_info(istep)%geoflag = geoflag
      velo_info(istep)%ntime = ntime
      velo_info(istep)%vz = 0.0

      zsurfp    = zsurf
      timesurfp = timesurf

      !call log_message("timesurf: " + timesurf + ", istep: " + istep + ", dt: " + dt)

    enddo ! istep = nstep, 1, -1















    ! WK: Pre-calculate positions and velocities for all time steps
    ! Backtrack for each time step from the surface:
    ! 1. time step: just the surface
    ! 2. time step: surface and 1 step back
    ! 3. time step: surface and 2 steps back
    ! 4. time step: surface and 3 steps back
    ! etc.
    !

    call log_message("Pre-calculating position and velocities of surface nodes")

    allocate(surface_pos(nsurf))

    do start_step = 1, nstep
      call log_message("Time step: " + start_step)
      write(file_id1, "(i4.4)") start_step

      call log_message("velocity_info_"//file_id1//".bin")
      inquire(file=run(1:nrun)//"/velocity_info_"//file_id1//".bin", exist=file_exists)

      if (file_exists .and. config%use_cached_files) then
        call log_message("Using already pre-calculated velocity_info files")
        cycle
      endif

      open(file_unit_velocity_info, file=run(1:nrun)//"/velocity_info_"//file_id1//".bin", &
            status="unknown", form="unformatted", access="stream")
      write(file_unit_velocity_info) start_step, nsurf

      do i = 1, nsurf
        surface_pos(i)%x = xdepth(start_step, i)
        surface_pos(i)%y = ydepth(start_step, i)
        surface_pos(i)%z = zdepth(start_step, i)

        ! if (surface_pos(i)%x /= surface_pos(i)%x) then
        !   call log_message("pecube_func.f90, surface_pos(i)%x is Nan, i:" + i + ", start_step: " + start_step)
        ! endif
        !
        ! if (surface_pos(i)%y /= surface_pos(i)%y) then
        !   call log_message("pecube_func.f90, surface_pos(i)%y is Nan, i: " + i + ", start_step: " + start_step)
        ! endif
        !
        ! if (surface_pos(i)%z /= surface_pos(i)%z) then
        !   call log_message("pecube_func.f90, surface_pos(i)%z is Nan, i: " + i + ", start_step: " + start_step)
        ! endif
      enddo

      ! Go back one step after another and calculate the positions and velocity for each step
      do current_step = start_step, 1, -1
        !call log_message("current_step: " + current_step)

        write(file_unit_velocity_info) current_step

        if (velo_info(current_step)%geoflag == 8) then
          call load_velocities(velocity_files(current_step), xlonmin, xlatmin, &
                velo_info(current_step)%zl, nx0, ny0, nz, dx, dy, &
                xdepth(current_step, :), ydepth(current_step, :), zdepth(current_step, :), &
                nsurf, vx_min, vy_min, vz_min, vx_max, vy_max, vz_max, config%use_new_velocity)
        endif

        do i = 1, nsurf
          call move_points(surface_pos(i)%x, surface_pos(i)%y, surface_pos(i)%z, &
                velo_info(current_step), file_unit_velocity_info, i, xmax, ymax)
        enddo
      enddo ! current_step = start_step, 1, -1
      close(file_unit_velocity_info)
    enddo ! start_step = 1, nstep

    ! If bore hole points are defined:
    if (number_of_borehole_ages > 0) then
      do start_step = 1, nstep
        !call log_message("Time step: " + start_step)
        write(file_id1, "(i4.4)") start_step

        inquire(file=run(1:nrun)//"/borehole_velocity_info_"//file_id1//".bin", exist=file_exists)

        if (file_exists .and. config%use_cached_files) then
          call log_message("Using already pre-calculated borehole_velocity_info files")
          cycle
        endif

        open(file_unit_borehole_velocity_info, file=run(1:nrun)//"/borehole_velocity_info_"//file_id1//".bin", &
              status="unknown", form="unformatted", access="stream")
        write(file_unit_borehole_velocity_info) start_step, number_of_borehole_ages

        do i = 1, number_of_borehole_ages
          model_pos(i)%x = borehole_ages_points(i)%x
          model_pos(i)%y = borehole_ages_points(i)%y
          model_pos(i)%z = borehole_ages_points(i)%z

          ! if (model_pos(i)%x /= model_pos(i)%x) then
          !   call log_message("model_pos(i)%x is Nan, i: " + i)
          ! endif
          !
          ! if (model_pos(i)%y /= model_pos(i)%y) then
          !   call log_message("model_pos(i)%y is Nan, i: " + i)
          ! endif
          !
          ! if (model_pos(i)%z /= model_pos(i)%z) then
          !   call log_message("model_pos(i)%z is Nan, i: " + i)
          ! endif

        enddo

        ! Go back one step after another and calculate the positions and velocity for each step
        do current_step = start_step, 1, -1
          write(file_unit_borehole_velocity_info) current_step

          if (velo_info(current_step)%geoflag == 8) then
            call load_velocities(velocity_files(current_step), xlonmin, xlatmin, &
                  velo_info(current_step)%zl, nx0, ny0, nz, dx, dy, &
                  xdepth(current_step, :), ydepth(current_step, :), zdepth(current_step, :), &
                  nsurf, vx_min, vy_min, vz_min, vx_max, vy_max, vz_max, config%use_new_velocity)
          endif

          do i = 1, number_of_borehole_ages
            call move_points(model_pos(i)%x, model_pos(i)%y, model_pos(i)%z, &
                  velo_info(current_step), file_unit_borehole_velocity_info, i, xmax, ymax)
          enddo
        enddo ! current_step = start_step, 1, -1
        close(file_unit_borehole_velocity_info)
      enddo ! start_step = 1, nstep
    endif ! number_of_borehole_ages > 0





    call log_message("Finished velocity pre-calculation")
















! reset the input file to its proper position

      call log_message('initializing')
      call logger_flush()

      rewind (7)
! 2011.07.25, WK: format specifier was not correct, now fixed!
      read (7,'(a,i10)') run, nrun
      read (7,*) num_topo_files,det_calc

      do i=1,config%num_of_age_flags
        read (7,*) age_flags(i)
      end do

      read (7,'(A)') (topo_file_name(k),k=1,num_topo_files)
      read (7,*) npe,nsurf,nzin,nelemsurf,zl,diffusivity,heatproduction,efold,&
                 shearint,friction

      ! Read in Nepal model geometry info
      ! This line lists the thermal diffusivity and heat production values for the
      ! Indian shield, Sub-Himalaya Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
      read (7,*) isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp

      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol
      nstep=abs(nstep)
      read (7,*) isoflag,tau,rhoc,rhom,heatcap
      read (7,*) nx,ny,nxiso,nyiso,nx0,ny0,dx,dy
      read (7,*) xstep,ystep,young,poisson,thickness
      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax
      read (7,*) vx_min, vy_min, vz_min
      read (7,*) vx_max, vy_max, vz_max

! WK, 2014.03.05: read in temperature file name
      read (7, '(A)') move_file

      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

      nnode = nsurf*nz
      nelem = nelemsurf*(nz-1)
      if (ilog.eq.1) call log_message('nnode/nelem= ' + nnode + "" + nelem)

! opens output files

      allocate (x(nnode),y(nnode),z(nnode),t(nnode))
      allocate (xp(nnode),yp(nnode),zp(nnode),tp(nnode))
      allocate (dummy(nnode),dummyp(nnode))
      allocate (icon(mpe,nelem))
      allocate (kfix(nnode))
      allocate (f(nnode),rebound(nnode))
      allocate (jrec(nstep+1))

! build 3D element connectivity

      ie=0
        do iesurf=1,nelemsurf
          do k=1,(nz-1)
          ie=ie+1
            do kk=1,npe
              icon(kk,ie)=(iconsurf(kk,iesurf)-1)*nz+k
              icon(kk+npe,ie)=icon(kk,ie)+1
            enddo
          enddo
        enddo

      if (ie.ne.nelem) then
        stop 'nelem mismatch'
      endif

! allocate reduced number of elemental matrices

      allocate (ael(mpe,mpe,nelem),bel(mpe,nelem))

      if (ilog.eq.1) then
        call log_message('icon')
        do ie=1,nelem
          call log_message("" + icon(:,ie))
        enddo
      endif

! finds neighbour connectivity matrix

      call find_neighbours (iconsurf,neighbour,npe,nelemsurf,nsurf)
      ielsurf = 1

! initialize global parameters
! alpha is the time integration parameter

      alpha    = 0.5
      time     = 0.
      timesurf = 0.
      irec     = 0
      jrec     = 0
      niter    = 0

! allocate erosion rate arrays - dwhipp (08/07)
      allocate(edot(nsurf),bg_edot(nsurf),topo_edot(nsurf))

! Initialize erosion rates:

      edot = 0.0
      bg_edot = 0.0
      topo_edot = 0.0

! allocate heat production and thermal conductivity arrays
! dwhipp (10/07)

























!*********************************************
!*********************************************
!** begining of surface stepping (stageing) **
!*********************************************
!*********************************************
      call log_message("Start of time stepping")

      do istep = 0, nstep
        write(file_id1, "(i4.4)") istep

        call log_message("Step " + istep + ", niter: " + niter)
        call logger_flush()

        if (istep > 0) then
          open(file_unit_temperature, file=run(1:nrun)//"/temperature_field_sub_"//file_id1//".bin", &
              status="unknown", form="unformatted", access="stream")
        endif
        ! write(file_unit_temperature, *) "#id, x, y, z, t"

        if (istep.ne.0) zsurfp=zsurf
        timesurfp = timesurf

! read the step information

! timesurf is the time since the begining of the experiment
! Peclet is the exhumation velocity since the end of the last stage (mm/year or km/mil year) from Pecube.in
! iout indicates if the T information is to be saved at the end of the stage

      read (7, *) timesurf,Peclet,iout,x1f,y1f,x2f,y2f,def,dif,thermflag,geoflag,theta,phi, &
                 mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2

      read (7, *) intrusion_start, intrusion_end, intrusion_temperature, temperature_hold_period

      if (intrusion_temperature > 0.0) then
        call log_message("Temperature hold period:" + temperature_hold_period + " [mil year]")
        prev_intrusion_temperature = 0.0
      endif

      if ((istep == nstep) .and. (config%error_iter_age_dec > 0.0)) then
        call log_message("Reducing time, timesurf: " + timesurf)
        timesurf = timesurf - config%error_iter_age_dec
        call log_message("New reduced time, timesurf: " + timesurf)
      endif

      call log_message('timesurf: ' + timesurf)

      if (istep == 0) then
        if (timesurf.gt.eps) then
        call log_message('first topography record must be at time zero ...')
        stop
        endif
      endif
      read (7, *) (zsurf(i),i=1,nsurf)

      read (7, *) temperature_file
      call log_message("pecube_func.f90: temperature file: " + trim(temperature_file))

      !call log_message("zsurfp: " + zsurfp(1) + ", " + zsurfp(2) + ", " + zsurfp(3))
      !call log_message("zsurf: " + zsurf(1) + ", " + zsurf(2) + ", " + zsurf(3))

! Stores the thermal flag values into array to be used to determine
! what thermal histories to write out at the end
        if(istep /= 0) then
          therm_his_val(istep) = int(thermflag)
        endif

! initial (or zeroth) step

        if (istep == 0) then
          zsurfp = zsurf
          in     = 0
          if (nzin.gt.0) then                                                   ! If the input nz value is positive, build geometry with constant z node plane
            do i=1,nsurf                                                        ! spacing equal to the total model thickness divided by the input nz
              do k=1,nz
                in   = in+1
                fact = float(k-1)/float(nz-1)
                if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
                xp(in)   = xsurf(i)
                x(in)    = xp(in)
                yp(in)   = ysurf(i)
                y(in)    = yp(in)
                zh       = zsurf(i)+zl
                zp(in)   = zh*fact
                kfix(in) = 0
                ! Modified to identify top and bottom fixed T B/Cs (original
                ! version below) - dwhipp (09/07)
                !if (k.eq.1 .or. k.eq.nz) kfix(in)=1
                if (k == 1) kfix(in)=1
                if (k == nz) kfix(in)=2
              enddo
            enddo
          else                                                                  ! If the input nz value is zero, build geometry with variable z node plane
            zmin = minval(zsurf)                                                ! spacing.  The spacing is such that from the top surface of the model down
            zmax = maxval(zsurf)                                                ! to 5 km below it, the spacing is 1: 1 with the x and y node spacing (xy_mean)
            do i=1,nsurf                                                        ! Below that, down to 15 km below the surface, the spacing is ~3:1 (3*xy_mean)
              do k=1,nz                                                         ! Below that, the spacing is ~9:1 down to the base of the model
                in=in+1                                                         ! Set the topography scaling factor to zero if the top surface is flat
                if ((zmax-zmin).eq.0.) then
                  scale_fact=0.
                else                                                            ! Define scaling factor to shift nodes beneath the topography
                  scale_fact=(zsurf(i)-zmin)*(1/(zsurf(i)+zl))
                endif
                if (k.gt.nz-(plane_store1)) then
                  if (k.eq.nz) then
                    scale_fact=0.
                    cur_depth=zsurf(i)+zl
                  else
                    cur_depth=cur_depth+(xy_mean/1000.)
                  endif
                else if (k.gt.nz-(plane_store2)) then
                  cur_depth=cur_depth+3*(xy_mean/1000.)
                else
                  if (k.eq.1) then
                    cur_depth=0.
                  else
                    cur_depth=cur_depth+fact*(xy_mean/1000.)
                  endif
                endif
                xp(in)   = xsurf(i)
                x(in)    = xp(in)
                yp(in)   = ysurf(i)
                y(in)    = yp(in)
                zp(in)   = cur_depth+cur_depth*scale_fact
                kfix(in) = 0
                ! Modified to identify top and bottom fixed T B/Cs (original
                !   version below) - dwhipp (09/07)
                !if (k.eq.1 .or. k.eq.nz) kfix(in)=1
                if (k.eq.1) kfix(in)=1
                if (k.eq.nz) kfix(in)=2
              enddo
            enddo
          endif ! nzin.gt.0

    !call log_message("zmin: " + zmin + ", zmax: " + zmax)

! calculates initial temperature

          call log_message("initial temperature calculation...")
          call log_message("nsurf: " + nsurf + ", nz: " + nz + ", nnode: " + nnode)


          do i = 1, nsurf
            tsurf = -zsurf(i) * tlapse + tmsl
            do k = 1, nz
              in     = (i - 1) * nz + k
              zh     = zsurf(i) + zl
              tp(in) = tsurf + (tmax - tsurf) * (zh - zp(in)) / zh
            enddo
          enddo

          t = tp

        endif ! istep == 0

        ! WK, 2019.04.29: load temperature filed from file and use it
        if (temperature_file(1:3) /= "Nil") then
          call load_temperatures(nnode, xp, yp, zp, tp, temperature_file)
          t = tp
        endif

        if(istep > 0) then
          jrec(istep+1)=jrec(istep)
        endif

! calculates time step length to satisfy stability and accuracy criteria
! note that an additional criterion must be fullfilled for stability
! but it is independent of the time step: (dx u / kappa) < 1 (the local
! Peclet number must be smaller than 1) where dx is the size of an element
! in a direction where a velocity u is applied, with kappa being the thermal
! diffusivity


    ! monte carlo simulation
    ! change Peclet here
    if (config%mc_use_monte_carlo) then
        if (allocated(config%mc_time_slices)) then
          actual_time = tfinal - timesurf

          do i = 1, size(config%mc_time_slices)
            if (actual_time >= config%mc_time_slices(i)%s_begin .and. actual_time <= config%mc_time_slices(i)%s_end) then
                Peclet = mc_random_erosion_rate(i)
            endif
          enddo
        else
          if (istep >= nstep) then
              Peclet = mc_random_erosion_rate(nstep)
          else
              Peclet = mc_random_erosion_rate(istep + 1)
          endif
        endif
    endif
    !call log_message("pecube_func.f90: Peclet after mc: " + Peclet)



      ! Use correct Peclet number for Nepal model geometry
      ! Added by dwhipp 11/07
      if (geoflag == 4) then
        tb_velo=(mft_ratein+mbt_ratein+mct_ratein-stf_ratein)
        if (20.-tb_velo.ge.tb_velo) then                                      ! Max velocity is underthrusting in Nepal model geometry
          PecletN = 20.-tb_velo
        else                                                                  ! Max velocity is overthrusting in Nepal model geometry
          PecletN = tb_velo
        endif
        PecletO = Peclet
        Peclet  = PecletN
      endif

      !call log_message("before find_dt, dt: " + dt + ", ntime: " + ntime)
      call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                    Peclet,timesurf,timesurfp,istep,eps, &
                    dt,ntime,istatic)
      !call log_message("after find_dt, dt: " + dt + ", ntime: " + ntime)

      !call log_message("Manual dt: " + manual_dt)

      if (manual_dt > 0.0) then
        call log_message("Use manual dt for this time step: " + dt + " -> " + manual_dt)
        dt = manual_dt
        call log_message("ntime old: " + ntime)
        ntime = int((timesurf - timesurfp) / dt) + 1
        call log_message("ntime new: " + ntime)
      endif

      call log_message("dt: " + dt + ", Peclet: (mm/year or km/mil year) " + Peclet + ", PecletO: " + PecletO + ", PecletN: " + PecletN)
      call log_message ("ntime: " + ntime + ", niter: " + niter)

      if (istep > 0) then
        write(file_unit_temperature) ntime, istep, nnode
      endif

      last_val = max(istep, 1)

      if (geoflag == 8) then
          call load_velocities(velocity_files(last_val), xlonmin, xlatmin, zl, nx0, ny0, nz, dx, dy, &
                 xdepth(last_val, :), ydepth(last_val, :), zdepth(last_val, :), nsurf, &
                 vx_min, vy_min, vz_min, vx_max, vy_max, vz_max, config%use_new_velocity)
      endif

      num_of_sub_steps(last_val) = ntime


























!******************************
!* beginning of time substeps *
!******************************

        do itime = 1, ntime
          call system_clock(clock_info(3)%sys_count, clock_info(3)%sys_count_rate, clock_info(3)%sys_count_max)
          time   = time+dt
          ftime  = float(itime)/float(ntime)
          ftimep = float(itime-1)/float(ntime)
          ftime  = (1.-exp(-ftime*tfinal/tau))/(1.-exp(-tfinal/tau))
          ftimep = (1.-exp(-ftimep*tfinal/tau))/(1.-exp(-tfinal/tau))
          call log_message('Doing time step ' + itime + ' of ' + ntime + &
                   ' in stage ' + istep + ' of ' + nstep + ' (' + niter + ')')
          call log_message("time: " + time + ", ftime: " + ftime)
          ! Reset Peclet for Nepal model geometry
          ! dwhipp 11/07
          if (geoflag == 4) Peclet=PecletO

! build new node geometry

          z = zp
          do i = 1, nsurf
            in     = i*nz
            inp    = i*nz-1
            zsurf1 = zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime+zl
            zsurf2 = zsurfp(i)+(zsurf(i)-zsurfp(i))*ftimep+zl
            ! calculate rate of surface change
            ! dwhipp 08/07
            if (istep == 0) then
              topo_edot(i) = 0.0
            else
              topo_edot(i) = real((zsurf1-zsurf2)/dt)
            endif
            topoa(i) = zsurf2
            topob(i) = zsurf1
            z(in)    = z(in)+zsurf1-zsurf2
            z(inp)   = z(inp)+(zsurf1-zsurf2)/2.
          enddo

        !call log_message('Min-Max Topo: ' + (minval(topoa)-zl) + "" + (maxval(topoa)-zl))

! isostatic rebound

        !call log_message("pecube_func.f90: isostatic rebound")

          if (isoflag.eq.1) then
            call isostatic_rebound (topoa,topob,nsurf,rsurf, &
                                       rhoc,rhom,nx,ny,nxiso,nyiso, &
                                       xstep,ystep,young,poisson,thickness)
          else
            !call log_message('Min-Max Rebound: ' + minval(rsurf) + "" + maxval(rsurf))
            rsurf = 0.
          endif

          do i = 1, nsurf
            do k = 1, nz
              rebound(k+(i-1)*nz) = rsurf(i)
            enddo
          enddo

          if (in.ne.nnode) then
            stop 'nnode mismatch'
          endif


! calculate current dynamic thermal conductivity
          if (has_dynamic_thermal_conductivity) then
            call calculate_conductivity(istep, current_dynamic_thermal_conductivity, dynamic_thermal_conductivity, ftime)

            do i = 1, current_dynamic_thermal_conductivity%num_of_layers
              call log_message("current dynamic thermal conductivity: " + current_dynamic_thermal_conductivity%layers(i)%td_value)
              call log_message("current dynamic thermal conductivity z1: " + current_dynamic_thermal_conductivity%layers(i)%z1)
              call log_message("current dynamic thermal conductivity z2: " + current_dynamic_thermal_conductivity%layers(i)%z2)
            enddo
          endif

! build local FE matrices

        !call log_message('Building matrix ')

        !call log_message("pecube_func.f90: make matrix")

        ! call cpu_time (times(7))
          call system_clock(clock_info(7)%sys_count, clock_info(7)%sys_count_rate, clock_info(7)%sys_count_max)
          do je=1,nelem

              do i=1,mpe/2
                  index1 = icon(i,je)
                  index2 = icon(i+(mpe/2),je)
                  if (z(index1) > z(index2)) then
                      call log_message("cell z order incorrect, je: " + je + ", i: " + i)
                      call log_message("index1, index2: " + index1 + "" + index2)
                      call log_message("z(index1): " + z(index1))
                      call log_message("x(index1), y(index1): " + x(index1) + ", " + y(index1))
                      call log_message("z(index2): " + z(index2))
                      call log_message("x(index2), y(index2): " + x(index2) + ", " + y(index2))
                      swap_z    = z(index1)
                      z(index1) = z(index2)
                      z(index2) = swap_z
                  endif

                  if (zp(index1) > zp(index2)) then
                      call log_message("cell zp order incorrect, je: " + je)
                      call log_message("index1, index2: " + index1 + "" + index2)
                      call log_message("zp(index1): " + zp(index1))
                      call log_message("xp(index1), yp(index1): " + xp(index1) + ", " + yp(index1))
                      call log_message("zp(index2): " + zp(index2) )
                      call log_message("xp(index2), yp(index2): " + xp(index2) + ", " + yp(index2))
                      swap_z     = zp(index1)
                      zp(index1) = zp(index2)
                      zp(index2) = swap_z
                  endif
              enddo

              ! tp is input here
              call make_matrix (mpe,ael(1,1,je),bel(1,je),icon(1,je),x,y,z,xp,&
                                 yp,zp,kfix,diffusivity,heatproduction,alpha,tp,nnode,&
                                 istatic,tlapse,tmsl,efold,je,&
                                 isdiff,shdiff,lhdiff,ghdiff,thdiff,&
                                 ishp,shhp,lhhp,ghhp,thhp,shearint,friction,&
                                 current_dynamic_thermal_conductivity, &
                                 has_dynamic_thermal_conductivity,velo_info(last_val))

          enddo ! je = 1, nelem

          call system_clock(clock_info(8)%sys_count, clock_info(8)%sys_count_rate, clock_info(8)%sys_count_max)

! build global RHS vector

        f = 0.d0
        do je = 1, nelem
          do k = 1, mpe
            ic    = icon(k, je)
            f(ic) = f(ic) + bel(k, je)
          enddo
        enddo

! solve global FE equations

        !call log_message('Solving matrix ')

        !call log_message("pecube_func.f90: solve iterative")

        ! call cpu_time (times(5))
        call system_clock(clock_info(5)%sys_count, clock_info(5)%sys_count_rate, clock_info(5)%sys_count_max)
        if (iterative == 1) then
          ! t is output (modified here):
          !call log_message("niter before solve_iterative: " + niter)
          call solve_iterative (mpe,ael,f,t,kfix,icon,nnode,nelem,niter,z,tlapse,tmsl,zl)
          !call log_message("niter after solve_iterative: " + niter)
        else
          stop 'solution strategy not implemented'
        endif







        !call cpu_time (times(6))
        call system_clock(clock_info(6)%sys_count, clock_info(6)%sys_count_rate, clock_info(6)%sys_count_max)

! sends short result to the screen

      !call log_message('Min-Max Temperature: ' + minval(t) + "" + maxval(t))

! stretch old grid

      if (nzin.gt.0) then                                                       ! Linearly shift nodes if the input nz value is positive
        do i=1,nsurf
          do k=1,nz
            in   = (i-1)*nz+k
            fact = float(k-1)/float(nz-1)
            if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
            zsurf1 = zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime
            zh     = zsurf1+zl
            zp(in) = zh*fact
          enddo
        enddo
      else                                                                      ! Shift variably spaced z node planes for new topography
        zmin = minval(zsurfp+(zsurf-zsurfp)*ftime)                              ! Calculate new zmin and zmax for new topo
        zmax = maxval(zsurfp+(zsurf-zsurfp)*ftime)
        do i=1,nsurf
          do k=1,nz
            in = (i-1)*nz+k
            if ((zmax-zmin).eq.0.) then                                         ! Set scale_fact to zero if surface topo is flat
              if (zmin.eq.0.) then
                scale_fact = 0.
              else
                scale_fact = (zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)*(1/zl)        ! Set scale factor to allow for topography shifts
              endif
            else                                                                ! Set scale factor to allow for topography shifts and relief scaling
                scale_fact = ((zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)-zmin_orig)*(1/zl)
            endif
            if (k.gt.nz-(plane_store1)) then
              if (k.eq.nz) then
                cur_depth  = (zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)+zl
                scale_fact = 0.
              else
                cur_depth=cur_depth+(xy_mean/1000.)
              endif
            else if (k.gt.nz-(plane_store2)) then
              cur_depth = cur_depth+3.*(xy_mean/1000.)
            else
              if (k.eq.1) then
                cur_depth = 0.
              else
                cur_depth = cur_depth+fact*(xy_mean/1000.)
              endif
            endif
            zp(in) = cur_depth+cur_depth*scale_fact
          enddo
        enddo
      endif ! nzin > 0

    !call log_message("zmin: " + zmin + ", zmax: " + zmax)

      if (itime == ntime) then
        ! Only on the last sub time step
        do i = 1, nsurf
          bg_edot(i) = Peclet * velo_info(last_val)%vz(i)
          edot(i) = bg_edot(i) - topo_edot(i)
        enddo
      endif

      !call log_message("intrusion_start: " + intrusion_start)
      !call log_message("intrusion_end: " + intrusion_end)
      !call log_message("intrusion_temperature: " + intrusion_temperature)
      !call log_message("prev_intrusion_temperature: " + prev_intrusion_temperature)
      !call log_message("temperature_hold_period: " + temperature_hold_period)

      if (temperature_hold_period > 0.0) then
        if (prev_intrusion_temperature > 0.0) then
          ! Move intrusion position according to the uplift
          intrusion_start = prev_intrusion_start + (dt * Peclet)
          intrusion_end = prev_intrusion_end + (dt * Peclet)
          intrusion_temperature = prev_intrusion_temperature
          temperature_hold_period = temperature_hold_period - dt
        endif
      else
          ! Turn off intrusion
          intrusion_temperature = 0.0
          prev_intrusion_temperature = 0.0
          temperature_hold_period = 0.0
          ! call log_message("Turn off intrusion")
      endif

      !call log_message("intrusion_start: " + intrusion_start)
      !call log_message("intrusion_end: " + intrusion_end)
      !call log_message("intrusion_temperature: " + intrusion_temperature)
      !call log_message("prev_intrusion_temperature: " + prev_intrusion_temperature)
      !call log_message("temperature_hold_period: " + temperature_hold_period)

      ! WK: temperature intrusion
      if (intrusion_temperature > 0.0) then
        ! Intrusion is active
        !call log_message("Intrusion is active, temperature: " + intrusion_temperature)

         do i = 1, nnode
           if ((z(i) >= intrusion_start) .and. (z(i) <= intrusion_end)) then
             t(i) = intrusion_temperature
           endif
         enddo

         if (temperature_hold_period > 0.0) then
            prev_intrusion_start = intrusion_start
            prev_intrusion_end = intrusion_end
            prev_intrusion_temperature = intrusion_temperature
         endif
      endif





! interpolate result onto undeformed mesh
      do i = 1, nsurf
        ij = (i-1)*nz+1
        ! t is input, tp is output here:
        call interpolate (t(ij),tp(ij),z(ij),zp(ij),nz)
      enddo

      ! WK: write temperature field to file
      ! We need to write the whole temperature field for all time steps
      ! since we need all the information for the back tracking in each time stlep

      if (istep > 0) then
        !write(file_unit_temperature) dt / dble(ntime), itime, tfinal - time
        write(file_unit_temperature) dt, itime, tfinal - time

        write(file_id1, "(i4.4)") istep

        !open(file_unit_scratch1, file="temperature_field_" // file_id1 // ".txt", status="unknown")
        do i = 1, nnode
          write(file_unit_temperature) i, x(i), y(i), z(i), t(i)
          !write(file_unit_scratch1, *) x(i), y(i), z(i), t(i)
        enddo
        !close(file_unit_scratch1)
      endif

      ! call cpu_time (times(4))
      call system_clock(clock_info(4)%sys_count, clock_info(4)%sys_count_rate, clock_info(4)%sys_count_max)

      call log_message('This time step (sys clock): ' + (dble(clock_info(4)%sys_count - clock_info(3)%sys_count) / &
          dble(clock_info(4)%sys_count_rate)))

!************************
!* end of sub time stepping *
!************************

    enddo ! itime = 1, ntime, sub time step loop



















      ! Write out Tecplot formatted output files
      if (has_dynamic_thermal_conductivity) then
        call calculate_conductivity(istep, current_dynamic_thermal_conductivity, dynamic_thermal_conductivity, 1.0_8)
      endif

      call tec_mat_output (x, y, z, t, nnode, icon, nelem, mpe, run, &
                     file_id1, xlonmin, xlatmin, nrun, current_dynamic_thermal_conductivity, &
                     has_dynamic_thermal_conductivity, velo_info(last_val) )

      ! Call subroutine to write out surface erosion rates
      ! dwhipp - (08/07)
      call erates (nsurf,xxx,xlonmin,xlonmax,xsurf,xmin,xmax,yyy,&
                   xlatmin,xlatmax,ysurf,ymin,ymax,zsurf,run,edot,file_id1,&
                   bg_edot,topo_edot,nrun,nx,nelemsurf)

! Stores current elevation values for later output
      if (istep /= 0) edot_store(istep,:) = edot

      if (istep > 0) then
        close(file_unit_temperature)
      endif

!*****************************
!** end of (outer) big time tepping **
!*****************************

      enddo ! istep = 0, nstep, big time step loop



















      ! WK: calculate total number of sub time steps
      do start_step = 1, nstep
        do current_step = start_step, 1, -1
          sub_step_sum(start_step) = sub_step_sum(start_step) + num_of_sub_steps(current_step)
        enddo
        call log_message("Number of total sub time steps: " + sub_step_sum(start_step) + &
                " in step: " + start_step + ", sub steps: " + num_of_sub_steps(start_step))
      enddo

      ! WK: Calculate the temperature for each back tracked surface node
      call log_message("Temperature calculation for each back tracked surface node (this will take some time...)")
      do start_step = 1, nstep
        call log_message("Time step: " + start_step)
        call logger_flush()

        write(file_id1, "(i4.4)") start_step

        file1 = run(1:nrun)//"/time_temperature_history_"//file_id1//".bin"
        file2 = run(1:nrun)//"/velocity_info_"//file_id1//".bin"
        call backtrack_temperature(start_step, nsurf, nnode, file1, file2, run, &
          config)

        if (number_of_borehole_ages > 0) then
          file1 = run(1:nrun)//"/borehole_time_temperature_history_"//file_id1//".bin"
          file2 = run(1:nrun)//"/borehole_velocity_info_"//file_id1//".bin"
          call backtrack_temperature(start_step, number_of_borehole_ages, nnode, &
            file1, file2, run, config)
        endif
      enddo ! start_step = 1, nstep

      call log_message("Finished temperature calculation")

      call log_message("xmin: " + minval(x) + ", xmax: " + maxval(x))
      call log_message("ymin: " + minval(y) + ", ymax: " + maxval(y))
      call log_message("zmin: " + minval(z) + ", zmax: " + maxval(z))

! Deallocations
! Moved 09/07 by cspath

      deallocate (ielsurf,neighbour)
      deallocate (topoa,topob,rsurf,rebound)
      deallocate (x,y,z,t)
      deallocate (xp,yp,zp,tp,dummy,dummyp)
      deallocate (icon,ael,bel)
      deallocate (kfix,f)
      deallocate (edot)


! Finds the number of requested thermal history outputs (specified in Pecube.in)
! Added 10/07 by cspath
      num_therm_his = 0
      do k = 1, nstep
        if (therm_his_val(k) == 1) num_therm_his = num_therm_his + 1
      enddo

    call log_message("pecube_func.f90: num_therm_his: " + num_therm_his)
      if (ilog == 1) close (9)





















! calculate FT ages
! Made age arrays 2,3-D to hold the age values of the surface nodes for
! different topography files as they are tracked
! fte is the Helium ages array, fta is fission track ages
! ftm is the muscovite ages
! The fte array is 3-D because of the extra Helium ages calculated within Mad_He
! in addition to finding ages for every time step

      call log_message('Calculating ages and thermal history...')
      call log_message('Please be patient. This may take awhile.')
      call log_message("age flags: " + age_flags)


      ! Initialize RDAAM calculation if needed. This has to be done only once since the
      ! values of the configuration parameters don't change.
      ! You can either have RDAAM or ZRDAAM but not both in the same simulation.
      if (has_RDAAM_ApUThHe(age_flags)) then
        if (has_RDAAM_ZirUThHe(age_flags)) then
          call log_message("Only one is allowed: RDAAM or ZRDAAM, but not both at the same time!")
          stop
        endif
        call rdaam_init(config%RDAAM_grain_radius, config%RDAAM_ppm_U, config%RDAAM_ppm_Th, config%RDAAM_ppm_Sm)
      else if (has_RDAAM_ZirUThHe(age_flags)) then
        call zrdaam_init(config%RDAAM_grain_radius, config%RDAAM_ppm_U, config%RDAAM_ppm_Th, config%RDAAM_ppm_Sm)
      endif


      allocate (header_info(6))

      if (allocated(surface_age_info)) then
        call deallocate_age_info(surface_age_info, nstep)
      endif

      call allocate_age_info(surface_age_info, nstep, nsurf, config%num_of_age_flags)

      if (number_of_borehole_ages > 0) then
        call allocate_age_info(borehole_age_info, nstep, number_of_borehole_ages, config%num_of_age_flags)
      endif


      call log_message("pecube_func.f90: therm_his_val: " + therm_his_val)

      ! Request from Byron Adams: write closure temperature to file
      ! Added 2014.06.04, Willi Kappler

      open (file_unit_closure, file=run(1:nrun)//"/closure_temp.dat", status="unknown")

      write (file_unit_closure, *) "Step, Surface, Method, Temperature, Closure Temperature"







      do m = 1, nstep
        if (therm_his_val(m) == 1) then
          call log_message("Timestep: " + m)
          call logger_flush()

          write(file_id1, "(I4.4)") m

          file1 = run(1:nrun)//"/time_temperature_history_"//file_id1
          call calculate_ages(config, sub_step_sum(m), file1, m, nsurf, &
              age_flags, surface_age_info(m), header_info, zl, config%export_surface_history)

          call log_message("pecube_func.f90, calculated ages: " + surface_age_info(m)%all_ages(1, :))

          if (number_of_borehole_ages > 0) then
            file1 = run(1:nrun)//"/borehole_time_temperature_history_"//file_id1
            call calculate_ages(config, sub_step_sum(m), file1, m, &
                number_of_borehole_ages, age_flags, borehole_age_info(m), &
                header_info, zl, config%export_borehole_history)
          endif

        else ! therm_his_val
           call log_message("pecube_func.f90, input 12d is zero in Pecube.in! (m = " + m + ")")
        endif
      enddo ! m = 1, nstep

      close (file_unit_closure) ! Closure temperature

      ! WK: write length distribution to file
      ! TODO: put this into the ages file below (ages_header)
      open (file_unit_length_dist, file=run(1:nrun)//"/length_distribution.dat", status="unknown")

      write(file_unit_length_dist, *) "time step, node, pecube x, pecube y, long, lat, mean track length, track length standard deviation"

      area_width     = abs(xlonmax - xlonmin)
      area_height    = abs(xlatmax - xlatmin)
      pecube_width   = xmax - xmin
      pecube_height  = ymax - ymin

      longitude_step = area_width / nx
      latitude_step  = area_height / ny

      call log_message("area_width: " + area_width + ", area_height: " + area_height)
      call log_message("pecube_width: " + pecube_width + ", pecube_height: " + pecube_height)
      call log_message("longitude_step: " + longitude_step + ", latitude_step: " + latitude_step)

      do i = 1, nsurf
        surf_longitude(i) = xlonmin + area_width * (xsurf(i) - xmin) / pecube_width
        surf_latitude(i)  = xlatmin + area_height * (ysurf(i) - ymin) / pecube_height
        write(file_unit_length_dist, *) nstep, i, xsurf(i), ysurf(i), surf_longitude(i), &
            surf_latitude(i), surface_age_info(nstep)%ftldmean(i), surface_age_info(nstep)%ftldsd(i)
      enddo

      close(file_unit_length_dist)

      ! Write out ages for all timesteps
      file1 = trim(run)//"/Ages_tec"
      call ages_header_surface (nsurf, xlonmin, xsurf, xlatmin, ysurf, zl,&
          xdepth, ydepth, zdepth, surface_age_info, file1, nstep, &
          header_info, therm_his_val, age_flags, nx, nelemsurf, &
          config%num_of_age_flags)

      if (number_of_borehole_ages > 0) then
        file1 = trim(run)//"/Ages_borehole_tec"
        call ages_header_borehole (number_of_borehole_ages, xlonmin, xlatmin, &
            zl, borehole_ages_points, borehole_age_info, file1, nstep, header_info, &
            therm_his_val, age_flags, config%num_of_age_flags)
      endif

      ! Checks if user specified to have detrital age calculations
      ! If det_calc is 1, then PDFs for all catchments that drain to edge of model are calculated
      ! Note: The catchments_output subroutine will only work for models from cascade because it uses
      ! the catchments found from a cascade output file (topo_tec*.dat)

      call log_message("--- detridal age ---")

      if (det_calc .eq. '1') then
        call log_message('Calculating Cascade catchment ages and PDFs...')
        read (7,*) node_thresh
        read (7,'(a100)') cascadedir
        call catchments_output (surface_age_info,nstep,nsurf,run,xdepth,ydepth,det_calc,&
                                edot_store,nrun,cascadedir,node_thresh,age_flags,header_info,&
                                topo_file_name,num_topo_files,&
                                config%num_of_age_flags)
      elseif (det_calc .ne. '0') then  ! If user specified a file with specific basin outlets
        call log_message('Calculating PDFs for user defined basins...')
        open (104,file='input/'//det_calc,status='old',iostat=ios)
        if (ios.ne.0) then
          open (104,file=det_calc,status='old')
        endif
        num_basins=0
        do
          read (104,*,end=20)
          num_basins=num_basins+1
        enddo
20      rewind(104)
        allocate (x_basin1(num_basins),y_basin1(num_basins),data_file1(num_basins),age_type1(num_basins))
        allocate (x_basin2(num_basins),y_basin2(num_basins),age_type2(num_basins))
        allocate (data_comp1(num_basins),data_comp2(num_basins))
        nilCount=0
        dataCount=0
        do j=1,num_basins      ! Iterates through basin list file
          read (104,*) x_basin,y_basin,age_type,data_file,data_comp   ! x position of outlet, y position of outlet, age type (1-8), associated data file or Nil
          if (data_file.eq.'Nil') then
            nilCount             = nilCount+1
            x_basin2(nilCount)   = x_basin
            y_basin2(nilCount)   = y_basin
            age_type2(nilCount)  = age_type
            data_comp2(nilCount) = data_comp
          else
            dataCount             = dataCount+1
            x_basin1(dataCount)   = x_basin
            y_basin1(dataCount)   = y_basin
            age_type1(dataCount)  = age_type
            data_file1(dataCount) = data_file
            data_comp1(dataCount) = data_comp
          endif
        enddo
        close (104)

        allocate (pages(nilCount,nx0*ny0),perates(nilCount,nx0*ny0))
        perates = 0.0_8
        pages = 0.0_8

        if (nilCount .ne. 0) then
          call find_upstream_points (topo_file_name,num_topo_files,x_basin2,y_basin2,nstep,nx0,ny0,dx,dy,&
               run,surface_age_info,age_type2,nsurf,xdepth,ydepth,edot_store,&
               xlonmin,xlatmin,nilCount,pages,perates,nrun)
        endif
        if (dataCount .ne. 0) then    ! If data file specified, read in data and error and output PDF
          do i=1,dataCount
            open (105,file='input/'//data_file1(i),status='old',iostat=ios)
            if (ios.ne.0) then
              open (105,file=data_file1(i),status='old')
            endif
            number=0
            do
              read (105,*,end=21)
              number=number+1
            enddo
21          rewind(105)
            allocate (pdf_ages(number),error(number))
            do j=1,number
              read (105,*) pdf_ages(j),error(j)
            enddo
            call pdfmaker_for_data (pdf_ages,error,number,age_type1(i),run,x_basin1(i),y_basin1(i),nrun)
            if (data_comp1(i) .eq. 'yes' .or. data_comp1(i) .eq. 'Yes') then
              foundtwo=0
              do j=1,nilCount
                if (x_basin2(j).eq.x_basin1(i) .and. y_basin2(j).eq.y_basin1(i) .and. age_type2(j).eq.age_type1(i)) then
                  if (data_comp2(j) .eq. 'yes' .or. data_comp2(j) .eq. 'Yes') then
                    foundtwo=1
                    call detrital_mc (data_file1(i),number,pdf_ages,error,&
                         pages(j,:),perates(j,:),nx0*ny0,run,x_basin1(i),y_basin1(i),nrun)
                  endif
                endif
              enddo
              if (foundtwo .eq. 0) then
! 2011.07.25, WK: format specifier not correct!
                call log_message('Warning: Could not find second basin data for basin:' + x_basin1(i) + ',' + y_basin1(i) + &
                    ' to run Monte Carlo test')
                call log_message('Check that the two basins have the same x position, y position, age type,')
                call log_message('and that both have "Yes" for the Monte Carlo test to run')
                call log_message('Skipping Monte Carlo test for this basin')
              endif
            endif
            deallocate (pdf_ages,error)
            close (105)
          enddo ! I = 1, dataCount
        endif
!         enddo
        deallocate (pages,perates)

        deallocate (x_basin1)
        deallocate (y_basin1)
        deallocate (data_file1)
        deallocate (age_type1)
        deallocate (x_basin2)
        deallocate (y_basin2)
        deallocate (age_type2)
        deallocate (data_comp1)
        deallocate (data_comp2)


      endif ! det_calc

      call logger_flush()


! calculate misfit
! Modified to calculate misfit values for all different types of ages
! Also allows for multiple data files to be entered in Pecube.in
! Will give a Comparison file for each data file and a total
! of all data files in Comparison.txt

! TODO: WK: Add misfit for fission track length


        open(file_unit_length_dist, file=run(1:nrun)//"/length_distribution2.dat", status="unknown")
        write(file_unit_length_dist, *) "long, lat, mean track length"

        call calculate_misfit(config, iconsurf, zsurf, surface_age_info(nstep)%ftldmean)

        call log_message("--- calculate misfit ---")

        allocate (age_obs(config%num_of_age_flags),&
          dage_obs(config%num_of_age_flags),&
          age_prd(config%num_of_age_flags))

        age_obs = 0.0
        dage_obs = 0.0
        age_prd = 0.0


        read (7,*) num_obsfiles
        allocate (xmisfit(num_obsfiles,config%num_of_age_flags),xmisfit_tot(num_obsfiles))

        xmisfit = 0.0
        xmisfit_tot = 0.0

        allocate (sum_xmisfit(8))
        sum_xmisfit = 0.0

        if (num_obsfiles.ne.0) then
        open (14,file=run(1:nrun)//'/Comparison_stats.txt',status='unknown')
        do j = 1, num_obsfiles

          write (comp,'(i3)') j
          if (j.lt.100) comp(1:1)='0'
          if (j.lt.10) comp(1:2)='00'

          call log_message("pecube_func.f90: write misfit to file: "//"Comparison_ages_tec_"//comp//".txt")
          open (13,file=run(1:nrun)//'/Comparison_ages_tec_'//comp//'.txt',status='unknown')
          xmisfit(j, :) = 0.0

          ! nobs is from obervation file specfied in Pecube.in
          read (7,*) nobs
          write (13,*) 'Number of observations: ',nobs

          call log_message('Number of observations: ' + nobs)


          write (13,'(A55)',ADVANCE="no") 'Data Set, Longitude, Latitude, Height Obs., Height Int., '
          write (13,'(A49)',ADVANCE="no") 'AHe Type, AHe Age Obs., AHe Age Error, AHe Age Prd., '
          write (13,'(A40)',ADVANCE="no") 'AFT Age Obs., AFT Age Error, AFT Age Prd., '
          write (13,'(A40)',ADVANCE="no") 'ZHe Age Obs., ZHe Age Error, ZHe Age Prd., '
          write (13,'(A40)',ADVANCE="no") 'ZFT Age Obs., ZFT Age Error, ZFT Age Prd., '
          write (13,'(A40)',ADVANCE="no") 'MAr Age Obs., MAr Age Error, MAr Age Prd.'
          write (13,*)

          do iobs=1,nobs
            read (7,*) He_flag
            read (7,*) xlonobs,xlatobs,heightobs,wobs1,wobs2,wobs3,wobs4,ieobs,&
                      age_obs(He_flag),dage_obs(He_flag)

              do k = 8, config%num_of_age_flags
                  read(7,*) age_obs(k), dage_obs(k)
              end do

              !call log_message("He_flag: " + He_flag + ", ieobs: " + ieobs + ", num_of_age_flags: " + config%num_of_age_flags)
              !call log_message("npe: " + npe + ", nelemsurf: " + nelemsurf)


              ! interpolate track length to observed coordinates:
              ! xlonobs, xlatobs

              observ_x = floor((xlonobs - xlonmin) / longitude_step) + 1
              observ_y = floor((xlatobs - xlatmin) / latitude_step) + 1

              p1 = (observ_y * nx) + observ_x
              p2 = (observ_y * nx) + observ_x + 1
              p3 = ((observ_y + 1) * nx) + observ_x
              p4 = ((observ_y + 1) * nx) + observ_x + 1

              !call log_message("nx: " + nx + ", ny: " + ny + ", nx0: " + nx0 + ", ny0: " + ny0)
              !call log_message("xlonobs: " + xlonobs + ", xlatobs: " + xlatobs)
              !call log_message("observ_x: " + observ_x +  ", observ_y: " + observ_y + ", nsurf: " + nsurf)
              !call log_message("p1: " + p1 + ", p2: " + p2 + ", p3: " + p3 + ", p4: " + p4)

              if ((p1 >= 1 .and. p1 <= nsurf) .and. &
                  (p2 >= 1 .and. p2 <= nsurf) .and. &
                  (p3 >= 1 .and. p3 <= nsurf) .and. &
                  (p4 >= 1 .and. p4 <= nsurf)) then

                ftldmean1 = linear_interpolation(surf_longitude(p1), surf_longitude(p2), &
                    surface_age_info(nstep)%ftldmean(p1), surface_age_info(nstep)%ftldmean(p2), xlonobs)

                ftldmean2 = linear_interpolation(surf_longitude(p3), surf_longitude(p4), &
                    surface_age_info(nstep)%ftldmean(p3), surface_age_info(nstep)%ftldmean(p4), xlonobs)

                ftldmean3 = linear_interpolation(surf_latitude(p1), surf_latitude(p3), &
                    ftldmean1, ftldmean2, xlatobs)

                write(file_unit_length_dist, *) xlonobs, xlatobs, ftldmean3

              endif

            if (age_flags(He_flag) == 0 .and. age_obs(He_flag) > 0) then
              call log_message('Warning: Apatite Helium age, type:' + He_flag + 'not being predicted; cannot compare to observed age')
              call log_message('For comparison enter 1 for prediction calculation in Pecube.in')
            endif

            if (iconsurf(1,ieobs) > nsurf .or. iconsurf(2,ieobs) > nsurf .or. iconsurf(3,ieobs) > nsurf .or. iconsurf(4,ieobs) > nsurf) then
              call log_message('Error: Surface node connectivity value out of range')
              call log_message('Check that thermochronological data file coordinates are same system (degrees/utm) as model')
              stop
            endif

            hei=wobs1*zsurf(iconsurf(1,ieobs)) &
              +wobs2*zsurf(iconsurf(2,ieobs)) &
              +wobs3*zsurf(iconsurf(3,ieobs)) &
              +wobs4*zsurf(iconsurf(4,ieobs))
            hei=hei*1000.

            ! j = 1 .. num_obsfiles
            write (13,'(i5,4f12.4)',ADVANCE="no") j,xlonobs,xlatobs,heightobs,hei

                age_prd(He_flag)=wobs1*surface_age_info(nstep)%all_ages(iconsurf(1,ieobs),He_flag) &
                      +wobs2*surface_age_info(nstep)%all_ages(iconsurf(2,ieobs),He_flag) &
                      +wobs3*surface_age_info(nstep)%all_ages(iconsurf(3,ieobs),He_flag) &
                      +wobs4*surface_age_info(nstep)%all_ages(iconsurf(4,ieobs),He_flag)

                age_prd(8)=wobs1*surface_age_info(nstep)%all_ages(iconsurf(1,ieobs),8) &
                      +wobs2*surface_age_info(nstep)%all_ages(iconsurf(2,ieobs),8) &
                      +wobs3*surface_age_info(nstep)%all_ages(iconsurf(3,ieobs),8) &
                      +wobs4*surface_age_info(nstep)%all_ages(iconsurf(4,ieobs),8)

                age_prd(9)=wobs1*surface_age_info(nstep)%all_ages(iconsurf(1,ieobs), 9) &
                      +wobs2*surface_age_info(nstep)%all_ages(iconsurf(2,ieobs), 9) &
                      +wobs3*surface_age_info(nstep)%all_ages(iconsurf(3,ieobs), 9) &
                      +wobs4*surface_age_info(nstep)%all_ages(iconsurf(4,ieobs), 9)

                age_prd(10)=wobs1*surface_age_info(nstep)%all_ages(iconsurf(1,ieobs), 10) &
                      +wobs2*surface_age_info(nstep)%all_ages(iconsurf(2,ieobs), 10) &
                      +wobs3*surface_age_info(nstep)%all_ages(iconsurf(3,ieobs), 10) &
                      +wobs4*surface_age_info(nstep)%all_ages(iconsurf(4,ieobs), 10)

                age_prd(11)=wobs1*surface_age_info(nstep)%all_ages(iconsurf(1,ieobs), 11) &
                      +wobs2*surface_age_info(nstep)%all_ages(iconsurf(2,ieobs), 11) &
                      +wobs3*surface_age_info(nstep)%all_ages(iconsurf(3,ieobs), 11) &
                      +wobs4*surface_age_info(nstep)%all_ages(iconsurf(4,ieobs), 11)

                if (age_prd(He_flag).eq.0) age_prd(He_flag)=-999
                if (age_prd(8).eq.0) age_prd(8)=-999
                if (age_prd(9).eq.0) age_prd(9)=-999
                if (age_prd(10).eq.0) age_prd(10)=-999
                if (age_prd(11).eq.0) age_prd(11)=-999

                write (13,'(i7,3f14.4)',ADVANCE="no") He_flag,age_obs(He_flag),dage_obs(He_flag),age_prd(He_flag)
                write (13,'(3f13.4)',ADVANCE="no") age_obs(9),dage_obs(9),age_prd(9)
                write (13,'(3f13.4)',ADVANCE="no") age_obs(8),dage_obs(8),age_prd(8)
                write (13,'(3f13.4)',ADVANCE="no") age_obs(10),dage_obs(10),age_prd(10)
                write (13,'(3f14.4)',ADVANCE="no") age_obs(11),dage_obs(11),age_prd(11)

              if (age_obs(He_flag) > 0.0 .and. age_flags(He_flag) == 1) then
                !call log_message("xmisfit(j, He_flag), before: " + xmisfit(j, He_flag))
                !call log_message("age_obs(He_flag)" + age_obs(He_flag))
                !call log_message("age_prd(He_flag)" + age_prd(He_flag))
                !call log_message("dage_obs(He_flag)" + dage_obs(He_flag))
                xmisfit(j, He_flag) = xmisfit(j, He_flag) + abs((age_obs(He_flag) - age_prd(He_flag)) / dage_obs(He_flag))
                !call log_message("xmisfit(j, He_flag), after: " + xmisfit(j, He_flag))
              endif

              if (age_obs(9) > 0.0 .and. age_flags(9) == 1) then
                xmisfit(j,9) = xmisfit(j,9) + abs((age_obs(9)-age_prd(9)) / dage_obs(9))
              endif

              if (age_obs(8) > 0.0 .and. age_flags(8) == 1) then
                xmisfit(j,8) = xmisfit(j,8) + abs((age_obs(8)-age_prd(8)) / dage_obs(8))
              endif

              if (age_obs(10) > 0.0 .and. age_flags(10) == 1) then
                xmisfit(j,10) = xmisfit(j,10) + abs((age_obs(10)-age_prd(10)) / dage_obs(10))
              endif

              if (age_obs(11) > 0.0 .and. age_flags(11) == 1) then
                xmisfit(j,11) = xmisfit(j,11) + abs((age_obs(11)-age_prd(11)) / dage_obs(11))
              endif

              ! TODO: add age observation for new ages

              write (13,*)

            enddo ! iobs=1,nobs

          if (age_flags(8) > 0 .and. age_obs(8) > 0) then
            call log_message('Warning: Zircon Helium age not being predicted; cannot compare to observed age')
            call log_message('For comparison enter 1 for prediction calculation in Pecube.in')
          endif
          if (age_flags(9) > 0 .and. age_obs(9) > 0) then
            call log_message('Warning: Apatite Fission Track age not being predicted; cannot compare to observed age')
            call log_message('For comparison enter 1 for prediction calculation in Pecube.in')
          endif
          if (age_flags(10) > 0 .and. age_obs(10) > 0) then
            call log_message('Warning: Zircon Fission Track age not being predicted; cannot compare to observed age')
            call log_message('For comparison enter 1 for prediction calculation in Pecube.in')
          endif
          if (age_flags(11) > 0 .and. age_obs(11) > 0) then
            call log_message('Warning: Muscovite age not being predicted; cannot compare to observed age')
            call log_message('For comparison enter 1 for prediction calculation in Pecube.in')
          endif

          !call log_message("xmisfit: " + xmisfit(j, :))

          xmisfit(j, 1) = xmisfit(j, 1) + xmisfit(j, 2) +xmisfit(j, 3) + xmisfit(j, 4) + xmisfit(j, 5) + xmisfit(j, 6) + xmisfit(j, 7)
          xmisfit_tot(j) = xmisfit(j, 1) + xmisfit(j, 8) + xmisfit(j, 9) + xmisfit(j, 10) + xmisfit(j, 11)

          if (nobs > 0) xmisfit_tot(j)=sqrt(xmisfit_tot(j))

          if (xmisfit(j,1) > 0.0) then
              write (14,*) "Misfit AHe:",xmisfit(j,1)
          endif

          if (xmisfit(j,9) > 0.0) then
              write (14,*) "Misfit AFt:",xmisfit(j,9)
          endif

          if (xmisfit(j,8) > 0.0) then
              write (14,*) "Misfit ZHe:",xmisfit(j,8)
          endif

          if (xmisfit(j,10) > 0.0) then
              write (14,*) "Misfit ZFt:",xmisfit(j,10)
          endif

          if (xmisfit(j,11) > 0.0) then
              write (14,*) "Misfit MAr:",xmisfit(j,11)
          endif

          write (14,*) 'Total Misfit: ',xmisfit_tot(j)
          write (14,*)

          close (13)
        enddo ! j = 1, num_obsfiles

        write (14,*) 'Total Misfits for all data sets:'

        do n = 1, num_obsfiles
          sum_xmisfit(1) = sum_xmisfit(1) + xmisfit(n,1)
          sum_xmisfit(2) = sum_xmisfit(2) + xmisfit(n,9)
          sum_xmisfit(3) = sum_xmisfit(3) + xmisfit(n,8)
          sum_xmisfit(4) = sum_xmisfit(4) + xmisfit(n,10)
          sum_xmisfit(5) = sum_xmisfit(5) + xmisfit(n,11)
          sum_xmisfit(6) = sum_xmisfit(6) + xmisfit_tot(n)
          sum_xmisfit(8) = sum_xmisfit(8) + (xmisfit_tot(n) * xmisfit_tot(n))
        enddo

        sum_xmisfit(7) = sum_xmisfit(6) / nobs
        sum_xmisfit(8) = sqrt(sum_xmisfit(8) / nobs)

        write (14,*) 'Misfit AHe:', sum_xmisfit(1), ', norm: ', (sum_xmisfit(1) / nobs)
        write (14,*) 'Misfit AFt:', sum_xmisfit(2), ', norm: ', (sum_xmisfit(2) / nobs)
        write (14,*) 'Misfit ZHe:', sum_xmisfit(3), ', norm: ', (sum_xmisfit(3) / nobs)
        write (14,*) 'Misfit ZFt:', sum_xmisfit(4), ', norm: ', (sum_xmisfit(4) / nobs)
        write (14,*) 'Misfit MAr:', sum_xmisfit(5), ', norm: ', (sum_xmisfit(5) / nobs)
        write (14,*) 'Total Misfit:', sum_xmisfit(6)
        write (14,*) 'Average Misfit1:', sum_xmisfit(7)
        write (14,*) 'Average Misfit2:', sum_xmisfit(8)

        call log_message("Misfit overview: ")
        call log_message("AHe: " + sum_xmisfit(1) + ", norm: " + (sum_xmisfit(1) / nobs) )
        call log_message("Aft: " + sum_xmisfit(2) + ", norm: " + (sum_xmisfit(2) / nobs) )
        call log_message("ZHe: " + sum_xmisfit(3) + ", norm: " + (sum_xmisfit(3) / nobs) )
        call log_message("ZFt: " + sum_xmisfit(4) + ", norm: " + (sum_xmisfit(4) / nobs) )
        call log_message("MAr: " + sum_xmisfit(5) + ", norm: " + (sum_xmisfit(5) / nobs) )
        call log_message("Total Misfit: " + sum_xmisfit(6))
        call log_message("Misfit avg1: " + sum_xmisfit(7) + ", avg2: " + sum_xmisfit(8))

        error_iter_misfit = sum_xmisfit(1) / nobs

        ! Needed for error_iteration
        last_peclet = velo_info(nstep)%Peclet

        deallocate (sum_xmisfit)
        endif
        close(7)
        close(14)

        close(file_unit_length_dist)

        deallocate (age_obs, dage_obs, age_prd)
        deallocate (xmisfit, xmisfit_tot)









! terminates the job
      deallocate (iconsurf)
      deallocate (header_info,therm_his_val)
      deallocate (xdepth,ydepth,zdepth,zsurfp)
      deallocate (zsurf)
      deallocate (dummy_values)

      if (config%pecube_run_mode /= pecube_run_mode_error_iter) then
        deallocate (edot_store)
        deallocate (xsurf, ysurf)
      endif

      deallocate (bg_edot)
      deallocate (velocity_files)
      deallocate (jrec)
      deallocate (topo_edot)
      deallocate (surf_longitude)
      deallocate (surf_latitude)

      if (number_of_borehole_ages > 0) then
        deallocate (borehole_ages_points)
        deallocate (model_pos, model_velo)
      endif

      deallocate(surface_pos)

      do i = 1, nstep
        deallocate(velo_info(i)%vz)
      enddo

      deallocate(velo_info)
      deallocate(sub_step_sum, num_of_sub_steps)

      call deallocate_thermal_conductivity(dynamic_thermal_conductivity, nstep)

      if (has_dynamic_thermal_conductivity) then
        deallocate(current_dynamic_thermal_conductivity%layers)
      end if

      call log_message("--- finished ---")

      ! call cpu_time (times(9))
      call system_clock(clock_info(9)%sys_count, clock_info(9)%sys_count_rate, clock_info(9)%sys_count_max)

      totalTime = dble(clock_info(9)%sys_count - clock_info(1)%sys_count) / dble(clock_info(9)%sys_count_rate)
      call log_message('Total times : ' + totalTime + " seconds")
      call log_message('Total times : ' + (totalTime / 60.0) + " minutes")
      call log_message('Total times : ' + (totalTime / 3600.0) + " hours")
      call logger_flush()
    end subroutine pecube_func

    function linear_interpolation(x1, x2, y1, y2, x)
      implicit none

      real(8), intent(in) :: x1
      real(8), intent(in) :: x2
      real(8), intent(in) :: y1
      real(8), intent(in) :: y2
      real(8), intent(in) :: x

      real(8) :: linear_interpolation

      real(8) :: m, b

      m = (y2 - y1) / (x2 - x1)
      b = ((y1 * x2) - (y2 * x1)) / (x2 - x1)

      linear_interpolation = (m * x) + b
    end function linear_interpolation
end module m_pecube_func
