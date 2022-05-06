module m_pecube_config
    implicit none

    private

    ! global pecube constants
    integer(4), parameter :: num_of_time_step_data = 8
    integer(4), parameter :: str_length = 255
    integer(4), parameter :: pecube_run_mode_normal = 0
    integer(4), parameter :: pecube_run_mode_error_iter = 1
    integer(4), parameter :: pecube_run_mode_monte_carlo = 2
    integer(4), parameter :: mc_max_number_of_erosion_rates = 30

    !> @struct time_step_t
    !! @brief contains configuration data for pecube each time step
    type :: time_step_t
        real(8) :: time
        real(8) :: amplification_factor
        real(8) :: vertical_offset
        logical :: temperature_history
        integer(4) :: kinematic_type
        real(8), dimension(num_of_time_step_data) :: data
    end type time_step_t

    !> @struct time_slice_t
    !! @brief contains begin and end of time slice for monte carlo intrusion simulation
    type :: time_slice_t
        real(8) :: s_begin
        real(8) :: s_end
    end type time_slice_t

    !> @struct config_t
    !! @brief contains configuration data for pecube from the input configuration file
    type :: config_t
        integer(4) :: mpi_total_num_of_cpu
        integer(4) :: mpi_current_cpu_id
        character(str_length) :: log_filename
        integer(4) :: pecube_run_mode
        logical :: just_velocity
        character(str_length) :: thermal_conductivity_file
        character(str_length) :: borehole_ages_file
        logical:: use_aft_ketcham

        ! number of time steps is stored in the global variable 'nstep'
        character(str_length) :: output_folder
        integer(4) :: topography_file_mode
        integer(4) :: topography_type
        character(str_length) :: topography_file_name
        integer(4) :: coordinate_system
        integer(4) :: nx
        integer(4) :: ny
        real(8) :: spacing_long
        real(8) :: spacing_lat
        integer(4) :: nskip
        real(8) :: location_long
        real(8) :: location_lat
        real(8) :: erosial_time_step
        type(time_step_t), dimension(:), allocatable :: time_step
        real(8) :: vx_min
        real(8) :: vx_max
        real(8) :: vy_min
        real(8) :: vy_max
        real(8) :: vz_min
        real(8) :: vz_max
        logical :: isostacy
        real(8) :: young_modules
        real(8) :: poisson_ratio
        real(8) :: elastic_plate
        integer(4) :: fft_grid_x
        integer(4) :: fft_grid_y
        real(8) :: model_thickness
        integer(4) :: number_of_z_planes
        real(8) :: thermal_conductivity
        real(8) :: specific_heat_capacity
        real(8) :: crustal_density
        real(8) :: mantle_density
        real(8) :: base_temperature
        real(8) :: z0_temperature
        real(8) :: atmospheric_lapse_rate
        real(8) :: crustal_heat_production
        real(8) :: e_fold_depth
        real(8) :: mantle_heat_production
        logical :: brittle_shear_heating
        real(8), dimension(5) :: nepal_model1
        real(8), dimension(5) :: nepal_model2
        real(8), dimension(5) :: nepal_model3
        real(8), dimension(5) :: nepal_model4
        ! thermocron file names, TODO
        ! age calculation parameters, TODO
        character(str_length) :: detridal_age
        integer(4) :: min_nodes
        character(str_length) :: cascade_out
        character(str_length) :: temperature_file
        character(str_length) :: move_file
        character(str_length) :: length_comparison_file
        integer(4) :: num_elements_surf
        integer(4) :: nsurf
        integer(4) :: npe
        integer(4) :: num_of_age_flags


        ! error iteration options
        integer(4) :: error_iter_radius
        real(8) :: error_iter_misfit_limit
        ! Age decrement for each Pecube run
        real(8) :: error_iter_age_dec
        real(8) :: error_iter_max_topography ! [m]

        ! monte carlo erosion rate options
        logical :: mc_use_monte_carlo
        real(8) :: mc_min_erosion_rate
        real(8) :: mc_max_erosion_rate
        integer(4) :: mc_num_of_simulations
        real(8) :: mc_tolerance_chi_squared
        logical :: mc_check_min_threshold
        real(8) :: mc_min_threshold_factor
        real(8) :: mc_erosion_step
        real(8), dimension(mc_max_number_of_erosion_rates) :: mc_fixed_erosion_rates
        character(str_length) :: mc_csv_input_file
        type(time_slice_t), dimension(:), allocatable :: mc_time_slices

        ! export time-temperature history to readable text file:
        logical :: export_surface_history
        logical :: export_borehole_history

        ! radius for exponential fall-off in find_temperature
        real(8) :: temperature_radius

        ! use new faster algorithm for velocity interpolation
        ! on regular grid
        logical :: use_new_velocity

        ! use cached files for velocity_info and borehole
        logical :: use_cached_files

        ! write VTK / Paraview output
        logical :: output_vtk

        ! RDAAM settings
        real(8) :: RDAAM_grain_radius
        real(8) :: RDAAM_ppm_U
        real(8) :: RDAAM_ppm_Th
        real(8) :: RDAAM_ppm_Sm
    end type config_t

    ! type(config_t) :: pecube_config

    ! public pecube_config
    public str_length
    public num_of_time_step_data
    public pecube_run_mode_normal
    public pecube_run_mode_error_iter
    public pecube_run_mode_monte_carlo
    public config_t
    public mc_max_number_of_erosion_rates
    public time_slice_t

end module m_pecube_config
