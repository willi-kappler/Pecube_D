! Written by Willi Kappler

module m_read_config_file
    use m_compiler
    use m_logger
    use m_pecube_config


    ! TODO: read variable thermal conductivity
    ! diffusivity = real((tcond/(rhoc*heatcap))*(1.d-3)**2*yrsec*1.e6)
    ! find_dt - not needed ?
    ! make_matrix



    ! Generate error code: openssl rand -hex 16

    implicit none
    save

    ! make everything private
    private

    ! parser internal constants
    integer(4), parameter :: global_file_unit = 80
    integer(4), parameter :: parse_mode_normal = 1
    integer(4), parameter :: parse_mode_time_step = 2
    integer(4), parameter :: parse_mode_thermocron = 3
    integer(4), parameter :: parse_mode_stop = 4

    ! keyword for configuration options
    character(*), parameter :: keyword_output_folder = "output_folder:" ! Input 1
    character(*), parameter :: keyword_topo_file_mode = "topo_file_mode:" ! Input 2
    character(*), parameter :: keyword_topo_type = "topo_type:" ! Input 3
    character(*), parameter :: keyword_topo_file_name = "topo_file_name:" ! Input 4
    character(*), parameter :: keyword_coordinate_system = "coordinate_system:" ! Input 5
    character(*), parameter :: keyword_nx = "nx:" ! Input 6
    character(*), parameter :: keyword_ny = "ny:" ! Input 6
    character(*), parameter :: keyword_spacing_long = "spacing_long:" ! Input 7
    character(*), parameter :: keyword_spacing_lat = "spacing_lat:" ! Input 7
    character(*), parameter :: keyword_nskip = "nskip:" ! Input 8
    character(*), parameter :: keyword_location_long = "location_long:" ! Input 9
    character(*), parameter :: keyword_location_lat = "location_lat:" ! Input 9
    character(*), parameter :: keyword_erosial_time_step = "erosial_time_step:" ! Input 10
    character(*), parameter :: keyword_erosial_time_scale = "erosial_time_scale:" ! Input 11

    ! Input 12
    character(*), parameter :: keyword_begin_time_step = "begin_time_step"
    ! time, amplification factor (af) = 1.0, vertical offset (vo) = 0.0, output (out) = on,
    ! kinematic field (kf) = 1, erosion rate (er) = 1.0
    !
    !
    character(*), parameter :: keyword_end_time_step = "end_time_step"

    ! Input 12.1
    character(*), parameter :: keyword_vx_min= "vx_min:"
    character(*), parameter :: keyword_vx_max = "vx_max:"
    character(*), parameter :: keyword_vy_min = "vy_min:"
    character(*), parameter :: keyword_vy_max = "vy_max:"
    character(*), parameter :: keyword_vz_min = "vz_min:"
    character(*), parameter :: keyword_vz_max = "vz_max:"

    ! Input 13
    character(*), parameter :: keyword_isostacy = "isostacy:"
    character(*), parameter :: keyword_young_modules = "young_modules:"
    character(*), parameter :: keyword_poisson_ratio = "poisson_ratio:"
    character(*), parameter :: keyword_elastic_plate = "elastic_plate:"
    character(*), parameter :: keyword_fft_grid_x = "fft_grid_x:"
    character(*), parameter :: keyword_fft_grid_y = "fft_grid_y:"

    ! Input 14
    character(*), parameter :: keyword_model_thickness = "model_thickness:"
    character(*), parameter :: keyword_number_of_z_planes = "number_of_z_planes:"
    character(*), parameter :: keyword_thermal_conductivity = "thermal_conductivity:"
    character(*), parameter :: keyword_specific_heat_capacity = "specific_heat_capacity:"
    character(*), parameter :: keyword_crustal_density = "crustal_density:"
    character(*), parameter :: keyword_mantle_density = "mantle_density:"
    character(*), parameter :: keyword_base_temperature = "base_temperature:"
    character(*), parameter :: keyword_z0_temperature = "z0_temperature:"
    character(*), parameter :: keyword_atmospheric_lapse_rate = "atmospheric_lapse_rate:"
    character(*), parameter :: keyword_crustal_heat_production = "crustal_heat_production:"
    character(*), parameter :: keyword_e_fold_depth = "e_fold_depth:"
    character(*), parameter :: keyword_mantle_heat_production = "mantle_heat_production:"
    character(*), parameter :: keyword_brittle_shear_heating = "brittle_shear_heating:"

    ! Input 15
    character(*), parameter :: keyword_nepal_model1 = "nepal_model1:"
    character(*), parameter :: keyword_nepal_model2 = "nepal_model2:"
    character(*), parameter :: keyword_nepal_model3 = "nepal_model3:"
    character(*), parameter :: keyword_nepal_model4 = "nepal_model4:"

    ! Input 16
    character(*), parameter :: keyword_begin_thermocron = "begin_thermocron"
    character(*), parameter :: keyword_end_thermocron = "end_thermocron"

    character(*), parameter :: keyword_age_calculation = "age_calculation:" ! Input 17
    character(*), parameter :: keyword_detridal_age = "detridal_age:" ! Input 18
    character(*), parameter :: keyword_min_nodes = "min_nodes:" ! Input 19
    character(*), parameter :: keyword_cascade_out = "cascade_out:" ! Input 20
    character(*), parameter :: keyword_temperature_file = "temperature_file:" ! Input 21
    character(*), parameter :: keyword_2d_move_file = "2d_move_file:" ! Input 22

    ! pecube run mode configuration
    character(*), parameter :: keyword_pecube_run_mode = "pecube_run_mode:"
    character(*), parameter :: keyword_pecube_end_run_mode = "pecube_end_run_mode"

    !! age calculation
    character(*), parameter :: keyword_use_aft_ketcham = "use_aft_ketcham:"

    ! error iteration options
    character(*), parameter :: keyword_error_iter_radius = "error_iter_radius:"
    character(*), parameter :: keyword_error_iter_misfit_limit = "error_iter_misfit_limit:"
    character(*), parameter :: keyword_error_iter_max_topography = "error_iter_max_topography:"

    ! monte carlo erosion rate options
    character(*), parameter :: keyword_mc_max_erosion_rate = "mc_max_erosion_rate:"
    character(*), parameter :: keyword_mc_min_erosion_rate = "mc_min_erosion_rate:"
    character(*), parameter :: keyword_mc_num_of_simulations = "mc_num_of_simulations:"
    character(*), parameter :: keyword_mc_tolerance_chi_squared = "mc_tolerance_chi_squared:"
    character(*), parameter :: keyword_mc_check_min_threshold = "mc_check_min_threshold:"
    character(*), parameter :: keyword_mc_min_threshold_factor = "mc_min_threshold_factor:"
    character(*), parameter :: keyword_mc_erosion_step = "mc_erosion_step:"
    character(*), parameter :: keyword_mc_csv_input_file = "mc_csv_input_file:"
    character(*), parameter :: keyword_mc_fix_erosion_rate = "mc_fix_erosion_rate:"
    character(*), parameter :: keyword_mc_time_slices = "mc_time_slices:"
    character(*), parameter :: keyword_length_comparison_file = "length_comparison_file:"

    ! RDAAM settings:
    character(*), parameter :: keyword_RDAAM_grain_radius = "RDAAM_grain_radius:"
    character(*), parameter :: keyword_RDAAM_ppm_U = "RDAAM_ppm_U:"
    character(*), parameter :: keyword_RDAAM_ppm_Th = "RDAAM_ppm_Th:"
    character(*), parameter :: keyword_RDAAM_ppm_Sm = "RDAAM_ppm_Sm:"

    ! character(*), parameter :: keyword_ = ""

    ! monte carlo intrusion options

    ! TODO: Add settings

    ! use 2D Move input to calculate velocity field or just load the velocity files:
    character(*), parameter :: keyword_just_velocity = "just_velocity:"

    ! dynamic thermal conductivity
    character(*), parameter :: keyword_thermal_conductivity_file = "thermal_conductivity_file:"

    ! age calculation for points inside the model
    character(*), parameter :: keyword_borehole_ages_file = "borehole_ages_file:"

    ! export time-temperature history to readable text file:
    character(*), parameter :: keyword_export_surface_history = "export_surface_history:"
    character(*), parameter :: keyword_export_borehole_history = "export_borehole_history:"

    ! radius of exponential fall-off for find_temperature
    character(*), parameter :: keyword_temperature_radius = "temperature_radius:"

    character(*), parameter :: keyword_use_new_velocity = "use_new_velocity:"
    character(*), parameter :: keyword_use_cached_files = "use_cached_files:"

    ! Write output for VTK / ParaView
    character(*), parameter :: keyword_output_vtk = "output_vtk:"

    ! keywords for configuration values

    ! Input 2
    character(*), parameter :: keyword_no_topo = "no_topo"
    character(*), parameter :: keyword_same_topo = "same_topo"
    character(*), parameter :: keyword_new_topo = "new_topo"

    ! Input 3
    character(*), parameter :: keyword_all_files = "all_files"
    character(*), parameter :: keyword_prefix = "prefix"
    character(*), parameter :: keyword_2d_move = "2d_move"

    ! Input 4
    character(*), parameter :: keyword_nil = "nil"

    ! Input 5
    character(*), parameter :: keyword_degrees = "degrees"
    character(*), parameter :: keyword_utm = "utm"

    ! Input 17
    character(*), parameter :: keyword_ahe0 = "ahe0"
    character(*), parameter :: keyword_ahe1 = "ahe1"
    character(*), parameter :: keyword_ahe2 = "ahe2"
    character(*), parameter :: keyword_ahe3 = "ahe3"
    character(*), parameter :: keyword_ahe4 = "ahe4"
    character(*), parameter :: keyword_ahe5 = "ahe5"
    character(*), parameter :: keyword_ahe6 = "ahe6"
    character(*), parameter :: keyword_aft = "aft"
    character(*), parameter :: keyword_zhe = "zhe"
    character(*), parameter :: keyword_zft = "zft"
    character(*), parameter :: keyword_musco = "musco"
    character(*), parameter :: keyword_ar_k_feldspar = "ar_k_feldspar"
    character(*), parameter :: keyword_ar_biotite = "ar_biotite"
    character(*), parameter :: keyword_ar_musco = "ar_musco"
    character(*), parameter :: keyword_ar_horn = "ar_horn"
    character(*), parameter :: keyword_a_uth = "a_uth"
    character(*), parameter :: keyword_biotite = "biotite"
    character(*), parameter :: keyword_ruite = "ruite"
    character(*), parameter :: keyword_titanite_upb = "titanite_upb"
    character(*), parameter :: keyword_zir_upb = "zir_upb"
    character(*), parameter :: keyword_titanite_uth = "titanite_uth"

    character(*), parameter :: keyword_normal_mode = "normal_mode"
    character(*), parameter :: keyword_error_iteration = "error_iteration"
    character(*), parameter :: keyword_monte_carlo_erosion = "monte_carlo_erosion"
    character(*), parameter :: keyword_monte_carlo_intrusion = "monte_carlo_intrusion"

    character(str_length) :: global_current_line
    character(str_length) :: global_last_comment
    character(str_length) :: global_line_value
    character(str_length) :: global_current_keyword

    integer(4) :: global_current_line_number
    integer(4) :: global_io_error
    integer(4) :: global_parse_mode


    ! export identifyer
    public read_config_file
    public num_of_time_step_data

    contains

    function check_for_keyword(keyword)
        implicit none

        character(*), intent(in) :: keyword

        logical :: check_for_keyword

        check_for_keyword = (global_current_line(1:index(global_current_line, " ")) == keyword)

        if (check_for_keyword) then
            global_current_keyword = keyword
        endif
    end function check_for_keyword

    function extract_value()
        implicit none

        character(:), allocatable :: extract_value

        extract_value = trim(adjustl(global_current_line((len(trim(global_current_keyword)) + 1):)))
    end function extract_value

    function check_for_value(keyword)
        implicit none

        character(*), intent(in) :: keyword

        logical :: check_for_value

        check_for_value = (trim(global_line_value) == keyword)

    end function check_for_value

    function extract_int_value()
        implicit none

        integer(4) :: extract_int_value

        global_line_value = extract_value()

        read(global_line_value, *, iostat=global_io_error) extract_int_value
        if (global_io_error /= 0) then
            call log_message("error reading value: " // trim(global_current_keyword))
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 6f51f76828571a132bf587921fa4af28")
            close(global_file_unit)
            error stop 1
        end if
    end function extract_int_value

    function extract_float_value()
        implicit none

        real(8) :: extract_float_value

        global_line_value = extract_value()

        read(global_line_value, *, iostat=global_io_error) extract_float_value
        if (global_io_error /= 0) then
            call log_message("error reading value: " // trim(global_current_keyword))
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 1676174598b9226c086d9c0961944cea")
            close(global_file_unit)
            error stop 1
        end if
    end function extract_float_value

    subroutine extract_5_float_values(result_array)
        implicit none

        real(8), dimension(5) :: result_array

        global_line_value = extract_value()

        read(global_line_value, *, iostat=global_io_error) result_array
        if (global_io_error /= 0) then
            call log_message("error reading value: " // global_current_keyword)
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 64546850e861da2cb5e80c6c765be9ba")
            close(global_file_unit)
            error stop 1
        end if
    end subroutine extract_5_float_values

    ! subroutine extract_4_float_values(result_array)
    !     implicit none

    !     real(8), dimension(4) :: result_array

    !     global_line_value = extract_value()

    !     read(global_line_value, *, iostat=global_io_error) result_array
    !     if (global_io_error /= 0) then
    !         call log_message("error reading value: " // global_current_keyword)
    !         call log_message("'" // trim(global_line_value) // "'")
    !         call log_message("line: " + global_current_line_number)
    !         call log_message("last comment: " // trim(global_last_comment))
    !         call log_message("parse mode: " + global_parse_mode)
    !         call log_message("error code: f76d1206ecc9e79aa66664cd84f636cd")
    !         close(global_file_unit)
    !         error stop 1
    !     end if
    ! end subroutine extract_4_float_values

    subroutine extract_int_float_values(int_value, float_value)
        implicit none

        integer(4) :: int_value
        real(8) :: float_value

        global_line_value = extract_value()

        read(global_line_value, *, iostat=global_io_error) int_value, float_value
        if (global_io_error /= 0) then
            call log_message("error reading value: " // global_current_keyword)
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 66b7e2faf2592cc6de56f371698efa72")
            close(global_file_unit)
            error stop 1
        end if
    end subroutine extract_int_float_values

    subroutine extract_option_bool(variable)
        implicit none

        logical, intent(out) :: variable

        integer(4) :: int_value

        call extract_option_2("off", "on", int_value)

        if (int_value == 0) then
            variable = .false.
        else if (int_value == 1) then
            variable = .true.
        else
            call log_message("impossible bool value: " + int_value)
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("allowed values: off, on")
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 695a7d74b5900d2d53c78d73e407c303")
            close(global_file_unit)
            error stop 1
        end if
    end subroutine extract_option_bool

    subroutine extract_option_2(val1, val2, variable)
        implicit none

        character(*), intent(in) :: val1
        character(*), intent(in) :: val2

        integer(4), intent(out) :: variable

        global_line_value = extract_value()
        if (check_for_value(val1)) then
            variable = 0
        else if (check_for_value(val2)) then
            variable = 1
        else
            call log_message("unknown value for: " + trim(global_current_keyword))
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " + trim(global_last_comment))
            call log_message("allowed values: " + val1 + ", " + val2)
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 21a8211a824a81f525ff2399fbf270a8")
            close(global_file_unit)
            error stop 1
        end if

    end subroutine extract_option_2

    subroutine extract_option_3(val1, val2, val3, variable)
        implicit none

        character(*), intent(in) :: val1
        character(*), intent(in) :: val2
        character(*), intent(in) :: val3

        integer(4), intent(out) :: variable

        global_line_value = extract_value()
        if (check_for_value(val1)) then
            variable = 0
        else if (check_for_value(val2)) then
            variable = 1
        else if (check_for_value(val3)) then
            variable = 2
        else
            call log_message("unknown value for: " + trim(global_current_keyword))
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("allowed values: " + val1 + ", " + val2 + ", " + val3)
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 832831c982118619bb9a8eb555775db6")
            close(global_file_unit)
            error stop 1
        end if

    end subroutine extract_option_3

    subroutine extract_option_4(val1, val2, val3, val4, variable)
        implicit none

        character(*), intent(in) :: val1
        character(*), intent(in) :: val2
        character(*), intent(in) :: val3
        character(*), intent(in) :: val4

        integer(4), intent(out) :: variable

        global_line_value = extract_value()
        if (check_for_value(val1)) then
            variable = 0
        else if (check_for_value(val2)) then
            variable = 1
        else if (check_for_value(val3)) then
            variable = 2
        else if (check_for_value(val4)) then
            variable = 3
        else
            call log_message("unknown value for: " + trim(global_current_keyword))
            call log_message("'" // trim(global_line_value) // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("allowed values: " + val1 + ", " + val2 + ", " + val3 + ", " + val4)
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: c69553e917586a305993eaed0e1210a5")
            close(global_file_unit)
            error stop 1
        end if

    end subroutine extract_option_4

    function extract_list()
        implicit none

        integer(4) :: num_of_commas, prev_pos, i, j
        character(20), dimension(:), allocatable :: extract_list

        global_line_value = extract_value()

        num_of_commas = 0

        do i = 1, str_length
            if (global_line_value(i:i) == ',') then
                num_of_commas = num_of_commas + 1
            endif
        enddo

        if (num_of_commas == 0) then
            call log_message("comma is missing in input file: " + trim(global_current_keyword))
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 693ab48e0e831c5238ad7210baaf1a48")
            close(global_file_unit)
            error stop 1
        endif

        allocate(extract_list(num_of_commas))

        prev_pos = 1
        j = 1

        do i = 1, str_length
            if (global_line_value(i:i) == ',') then
                extract_list(j) = global_line_value(prev_pos:i-1)
                prev_pos = i + 1
                j = j + 1
            endif
        enddo
    end function extract_list

    function extract_time_slice(string)
        implicit none

        integer(4) :: pos
        character(20), intent(in) :: string
        type(time_slice_t) :: extract_time_slice

        pos = index(string, "-")

        if (pos == 0) then
            call log_message("error reading value: " // trim(global_current_keyword))
            call log_message("'" // trim(string) // "'")
            call log_message("'-' missing")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: d68db23eaf4f7dfd9fe6f3236e823b9e")
            close(global_file_unit)
            error stop 1
        endif

        read(string(1:pos-1), *, iostat=global_io_error) extract_time_slice%s_begin
        if (global_io_error /= 0) then
            call log_message("error reading value: " // trim(global_current_keyword))
            call log_message("'" // trim(string) // "'")
            call log_message("begin is not a real number")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: f4f513c19cb5714a69746cdc0b83179e")
            close(global_file_unit)
            error stop 1
        end if

        read(string(pos+1:), *, iostat=global_io_error) extract_time_slice%s_end
        if (global_io_error /= 0) then
            call log_message("error reading value: " // trim(global_current_keyword))
            call log_message("'" // trim(string) // "'")
            call log_message("end is not a real number")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 6188f7326a30175f74126f008d5fd29b")
            close(global_file_unit)
            error stop 1
        end if

        if (extract_time_slice%s_begin >= extract_time_slice%s_end) then
            call log_message("error reading value: " // trim(global_current_keyword))
            call log_message("'" // trim(string) // "'")
            call log_message("begin >= end")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: b0a874ef102cf47c94b59117d0894959")
            close(global_file_unit)
            error stop 1
        end if

    end function extract_time_slice

    function get_time_slices(string_list)
        implicit none

        character(20), dimension(:), allocatable, intent(in) :: string_list
        type(time_slice_t), dimension(:), allocatable :: get_time_slices

        integer(4) :: num_of_slices, i

        num_of_slices = size(string_list)

        allocate(get_time_slices(num_of_slices))

        do i = 1, num_of_slices
            get_time_slices(i) = extract_time_slice(string_list(i))
        enddo

    end function get_time_slices

    subroutine handle_parse_mode_normal(config)
        implicit none

        type(config_t), intent(inout) :: config

        integer(4) :: int_value
        real(8) :: float_value

        character(20), dimension(:), allocatable :: string_list

        if (check_for_keyword(keyword_output_folder)) then
            config%output_folder = extract_value()
        else if (check_for_keyword(keyword_topo_file_mode)) then
            call extract_option_3(keyword_no_topo, &
                keyword_same_topo, keyword_new_topo, &
                config%topography_file_mode)
        else if (check_for_keyword(keyword_topo_type)) then
            call extract_option_3(keyword_all_files, &
                keyword_prefix, keyword_2d_move, &
                config%topography_type)
        else if (check_for_keyword(keyword_topo_file_name)) then
            config%topography_file_name = extract_value()
        else if (check_for_keyword(keyword_coordinate_system)) then
            call extract_option_2(keyword_degrees, &
                keyword_utm, config%coordinate_system)
        else if (check_for_keyword(keyword_nx)) then
            config%nx = extract_int_value()
        else if (check_for_keyword(keyword_ny)) then
            config%ny = extract_int_value()
        else if (check_for_keyword(keyword_spacing_long)) then
            config%spacing_long = extract_float_value()
        else if (check_for_keyword(keyword_spacing_lat)) then
            config%spacing_lat = extract_float_value()
        else if (check_for_keyword(keyword_nskip)) then
            config%nskip = extract_int_value()
        else if (check_for_keyword(keyword_location_long)) then
            config%location_long = extract_float_value()
        else if (check_for_keyword(keyword_location_lat)) then
            config%location_lat = extract_float_value()
        else if (check_for_keyword(keyword_erosial_time_step)) then
            config%erosial_time_step = extract_float_value()
        else if (check_for_keyword(keyword_begin_time_step)) then
            call log_message("begin parse mode time step")
            global_parse_mode = parse_mode_time_step
        else if (check_for_keyword(keyword_vx_min)) then
            config%vx_min = extract_float_value()
        else if (check_for_keyword(keyword_vx_max)) then
            config%vx_max = extract_float_value()
        else if (check_for_keyword(keyword_vy_min)) then
            config%vy_min = extract_float_value()
        else if (check_for_keyword(keyword_vy_max)) then
            config%vy_max = extract_float_value()
        else if (check_for_keyword(keyword_vz_min)) then
            config%vy_min = extract_float_value()
        else if (check_for_keyword(keyword_vz_max)) then
            config%vz_max = extract_float_value()
        else if (check_for_keyword(keyword_isostacy)) then
            call extract_option_bool(config%isostacy)
        else if (check_for_keyword(keyword_young_modules)) then
            config%young_modules = extract_float_value()
        else if (check_for_keyword(keyword_poisson_ratio)) then
            config%poisson_ratio = extract_float_value()
        else if (check_for_keyword(keyword_elastic_plate)) then
            config%elastic_plate = extract_float_value()
        else if (check_for_keyword(keyword_fft_grid_x)) then
            config%fft_grid_x = extract_int_value()
        else if (check_for_keyword(keyword_fft_grid_y)) then
            config%fft_grid_y = extract_int_value()
        else if (check_for_keyword(keyword_model_thickness)) then
            config%model_thickness = extract_float_value()
        else if (check_for_keyword(keyword_number_of_z_planes)) then
            config%number_of_z_planes = extract_int_value()
        else if (check_for_keyword(keyword_thermal_conductivity)) then
            config%thermal_conductivity = extract_float_value()
        else if (check_for_keyword(keyword_specific_heat_capacity)) then
            config%specific_heat_capacity = extract_float_value()
        else if (check_for_keyword(keyword_crustal_density)) then
            config%crustal_density = extract_float_value()
        else if (check_for_keyword(keyword_mantle_density)) then
            config%mantle_density = extract_float_value()
        else if (check_for_keyword(keyword_base_temperature)) then
            config%base_temperature = extract_float_value()
        else if (check_for_keyword(keyword_z0_temperature)) then
            config%z0_temperature = extract_float_value()
        else if (check_for_keyword(keyword_atmospheric_lapse_rate)) then
            config%atmospheric_lapse_rate = extract_float_value()
        else if (check_for_keyword(keyword_crustal_heat_production)) then
            config%crustal_heat_production = extract_float_value()
        else if (check_for_keyword(keyword_e_fold_depth)) then
            config%e_fold_depth = extract_float_value()
        else if (check_for_keyword(keyword_mantle_heat_production)) then
            config%mantle_heat_production = extract_float_value()
        else if (check_for_keyword(keyword_brittle_shear_heating)) then
            call extract_option_bool(config%brittle_shear_heating)
        else if (check_for_keyword(keyword_nepal_model1)) then
            call extract_5_float_values(config%nepal_model1)
        else if (check_for_keyword(keyword_nepal_model2)) then
            call extract_5_float_values(config%nepal_model2)
        else if (check_for_keyword(keyword_nepal_model3)) then
            call extract_5_float_values(config%nepal_model3)
        else if (check_for_keyword(keyword_nepal_model4)) then
            call extract_5_float_values(config%nepal_model4)
        else if (check_for_keyword(keyword_begin_thermocron)) then
            call log_message("begin parse mode thermocron")
            global_parse_mode = parse_mode_thermocron
        else if (check_for_keyword(keyword_age_calculation)) then
            call log_message("parse age calculation options, TODO...")
        else if (check_for_keyword(keyword_detridal_age)) then
            config%detridal_age = extract_value()
        else if (check_for_keyword(keyword_min_nodes)) then
            config%min_nodes = extract_int_value()
        else if (check_for_keyword(keyword_cascade_out)) then
            config%cascade_out = extract_value()
        else if (check_for_keyword(keyword_temperature_file)) then
            config%temperature_file = extract_value()
        else if (check_for_keyword(keyword_2d_move_file)) then
            config%move_file = extract_value()
        else if (check_for_keyword(keyword_pecube_run_mode)) then
            call extract_option_4(keyword_normal_mode, &
                keyword_error_iteration, keyword_monte_carlo_erosion, &
                keyword_monte_carlo_intrusion, config%pecube_run_mode)
        else if (check_for_keyword(keyword_pecube_end_run_mode)) then
            call log_message("parse mode stop")
            global_parse_mode = parse_mode_stop
        else if (check_for_keyword(keyword_error_iter_radius)) then
            config%error_iter_radius = extract_int_value()
        else if (check_for_keyword(keyword_error_iter_misfit_limit)) then
            config%error_iter_misfit_limit = extract_float_value()
        else if (check_for_keyword(keyword_error_iter_max_topography)) then
            config%error_iter_max_topography = extract_float_value()
        else if (check_for_keyword(keyword_mc_min_erosion_rate)) then
            config%mc_min_erosion_rate = extract_float_value()
        else if (check_for_keyword(keyword_mc_max_erosion_rate)) then
            config%mc_max_erosion_rate = extract_float_value()
        else if (check_for_keyword(keyword_mc_num_of_simulations)) then
            config%mc_num_of_simulations = extract_int_value()
        else if (check_for_keyword(keyword_mc_tolerance_chi_squared)) then
            config%mc_tolerance_chi_squared = extract_float_value()
        else if (check_for_keyword(keyword_mc_check_min_threshold)) then
            call extract_option_bool(config%mc_check_min_threshold)
        else if (check_for_keyword(keyword_mc_min_threshold_factor)) then
            config%mc_min_threshold_factor = extract_float_value()
        else if (check_for_keyword(keyword_mc_erosion_step)) then
            config%mc_erosion_step = extract_float_value()
        else if (check_for_keyword(keyword_mc_csv_input_file)) then
            config%mc_csv_input_file = trim(extract_value())
        else if (check_for_keyword(keyword_mc_fix_erosion_rate)) then
          call extract_int_float_values(int_value, float_value)
          if (int_value >= 1 .and. int_value <= mc_max_number_of_erosion_rates) then
            config%mc_fixed_erosion_rates(int_value) = float_value
          else
            call log_message("erosion rate index out of range:")
            call log_message("index: " + int_value)
            call log_message("erosion rate: " + float_value)
            call log_message("index must be between 1 and " + mc_max_number_of_erosion_rates)
            call log_message("'" // global_current_line // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 1b7be42025025ae23da85a6ca9529be7")
            close(global_file_unit)
            error stop 1
          end if
        else if (check_for_keyword(keyword_mc_time_slices)) then
            string_list = extract_list()
            config%mc_time_slices = get_time_slices(string_list)
            deallocate(string_list)
        else if (check_for_keyword(keyword_length_comparison_file)) then
          config%length_comparison_file = trim(extract_value())
        else if (check_for_keyword(keyword_just_velocity)) then
          call extract_option_bool(config%just_velocity)
        else if (check_for_keyword(keyword_thermal_conductivity_file)) then
          config%thermal_conductivity_file = trim(extract_value())
        else if (check_for_keyword(keyword_borehole_ages_file)) then
          config%borehole_ages_file = trim(extract_value())
        else if (check_for_keyword(keyword_export_surface_history)) then
          call extract_option_bool(config%export_surface_history)
        else if (check_for_keyword(keyword_export_borehole_history)) then
          call extract_option_bool(config%export_borehole_history)
        else if (check_for_keyword(keyword_temperature_radius)) then
          config%temperature_radius = extract_float_value()
        else if (check_for_keyword(keyword_use_new_velocity)) then
          call extract_option_bool(config%use_new_velocity)
        else if (check_for_keyword(keyword_use_cached_files)) then
          call extract_option_bool(config%use_cached_files)
        else if (check_for_keyword(keyword_output_vtk)) then
          call extract_option_bool(config%output_vtk)
        else if (check_for_keyword(keyword_use_aft_ketcham)) then
          call extract_option_bool(config%use_aft_ketcham)
        else if (check_for_keyword(keyword_RDAAM_grain_radius)) then
            config%RDAAM_grain_radius = extract_float_value()
        else if (check_for_keyword(keyword_RDAAM_ppm_U)) then
            config%RDAAM_ppm_U = extract_float_value()
        else if (check_for_keyword(keyword_RDAAM_ppm_Th)) then
            config%RDAAM_ppm_Th = extract_float_value()
        else if (check_for_keyword(keyword_RDAAM_ppm_Sm)) then
            config%RDAAM_ppm_Sm = extract_float_value()
!        else if (check_for_keyword(keyword_)) then
!        else if (check_for_keyword(keyword_)) then
        else
            call log_message("unknown keyword:")
            call log_message("'" // global_current_line // "'")
            call log_message("line: " + global_current_line_number)
            call log_message("last comment: " // trim(global_last_comment))
            call log_message("parse mode: " + global_parse_mode)
            call log_message("error code: 96dad153ed7d44ebd6c1d7e1f0fb1f86")
            close(global_file_unit)
            error stop 1
        end if
    end subroutine handle_parse_mode_normal

    subroutine handle_parse_mode_time_step
        implicit none

        call log_message("parse mode: time step")

        if (check_for_keyword(keyword_end_time_step)) then
            call log_message("end parse mode time step")
            global_parse_mode = parse_mode_normal
        else
            call log_message("global_current_line: " + global_current_line)
        end if
    end subroutine handle_parse_mode_time_step

    subroutine handle_parse_mode_thermocron
        implicit none

        call log_message("parse mode: thermocron")

        if (check_for_keyword(keyword_end_thermocron)) then
            call log_message("end parse mode thermocron")
            global_parse_mode = parse_mode_normal
        else
            call log_message("global_current_line: " + global_current_line)
        end if
    end subroutine handle_parse_mode_thermocron


    subroutine read_config_file(config_file_name, config)
        implicit none

        character(*), intent(in) :: config_file_name
        type(config_t), intent(inout) :: config

        logical :: file_exists

        ! Initialize Pecube configuration with default values
        config%mc_fixed_erosion_rates = -100.0
        config%thermal_conductivity_file = ""
        config%borehole_ages_file = ""
        config%export_surface_history = .false.
        config%export_borehole_history = .false.
        config%temperature_radius = 1.0
        config%length_comparison_file = ""
        config%use_new_velocity = .false.
        config%use_cached_files = .false.
        config%error_iter_radius = 10
        config%error_iter_misfit_limit = 5.0
        ! max z in [m]
        config%error_iter_max_topography = 10000.0
        config%output_vtk = .false.
        config%error_iter_age_dec = 0.0
        config%mc_min_erosion_rate = 0.0
        config%mc_erosion_step = 0.1
        config%use_aft_ketcham = .true.
        config%RDAAM_grain_radius = 100.0
        config%RDAAM_ppm_U = 10.0
        config%RDAAM_ppm_Th = 40.0
        config%RDAAM_ppm_Sm = 0.0

        if (allocated(config%mc_time_slices)) then
            deallocate(config%mc_time_slices)
        endif

        call log_message("Trying to open input file: " + config_file_name)

        inquire(file=config_file_name, exist=file_exists, iostat=global_io_error)

        if (global_io_error /= 0) then
            call log_message("error opening file: " // trim(config_file_name))
            call log_message("io_error: " + global_io_error)
            call log_message("error code: c4e7744497ddb8551da996708cdff4ec")
            error stop 1
        end if

        if (.not. file_exists) then
            call log_message("error opening file: file does not exist: " // trim(config_file_name))
            call log_message("error code: d9c0b9c8883003afccf7626ccac27b0f")
            error stop 1
        end if

        open(global_file_unit, file=config_file_name, iostat=global_io_error, status="old", action="read")

        if (global_io_error /= 0) then
            call log_message("error opening file: " // trim(config_file_name))
            call log_message("io_error: " + global_io_error)
            call log_message("error code: 174903e2cd234dbdd6f2e5eb4ea3ac1a")
            error stop 1
        end if

        global_current_line_number = 0
        global_parse_mode = parse_mode_normal

        do
            global_current_line(1:str_length) = " "
            global_line_value(1:str_length) = " "
            read(global_file_unit, "(a)", iostat=global_io_error) global_current_line

            if (is_iostat_end(global_io_error)) then
                call log_message("end of file reached")
                call log_message("parse mode: " + global_parse_mode)
                call log_message("total number of lines: " + global_current_line_number)
                exit
            else if (global_io_error /= 0) then
                call log_message("error reading from file: " // trim(config_file_name))
                call log_message("line number: " + global_current_line_number)
                call log_message("line: " // trim(global_current_line))
                call log_message("last comment: " // trim(global_last_comment))
                call log_message("io_error: " + global_io_error)
                call log_message("parse mode: " + global_parse_mode)
                call log_message("error code: 80a988eb4b8585f17d8ee1102ea4e7eb")
                close(global_file_unit)
                error stop 1
            end if

            global_current_line_number = global_current_line_number + 1

            if ((global_current_line(1:1) == "#") .or. (global_current_line(1:1) == "$")) then
                ! ignore comments
                global_last_comment = global_current_line
                cycle
            else if (global_current_line(1:1) == " ") then
                ! ignore blank lines
                cycle
            else
                if (global_parse_mode == parse_mode_normal) then
                    call handle_parse_mode_normal(config)
                    cycle
                else if (global_parse_mode == parse_mode_time_step) then
                    call handle_parse_mode_time_step()
                    cycle
                else if (global_parse_mode == parse_mode_thermocron) then
                    call handle_parse_mode_thermocron()
                    cycle
                else if (global_parse_mode == parse_mode_stop) then
                    exit
                else
                    call log_message("unknown parse mode: " + global_parse_mode)
                    call log_message("'" // global_current_line // "'")
                    call log_message("line: " + global_current_line_number)
                    call log_message("last comment: " // trim(global_last_comment))
                    call log_message("error code: a5930ff58f456fd90c35c1ebc4fdc329")
                    close(global_file_unit)
                    error stop 1
                end if
            end if
        end do

        close(global_file_unit)

        if ((global_parse_mode /= parse_mode_normal) .and. (global_parse_mode /= parse_mode_stop)) then
            call log_message("not normal parse mode: " + global_parse_mode)
            call log_message("missing keyword: end_")
            call log_message("error code: c25cb76c9d795e511a61390289086ac4")
            error stop 1
        end if

        ! Sanity checks for values:

        if (config%error_iter_radius < 0) then
            config%error_iter_radius = 0
        endif

        if (config%mc_min_erosion_rate >= config%mc_max_erosion_rate) then
            call log_message("Min erosion rate >= max erosion rate!")
            call log_message("mc_min_erosion_rate: " + config%mc_min_erosion_rate)
            call log_message("mc_max_erosion_rate: " + config%mc_max_erosion_rate)
            error stop 1
        endif

        call log_message("configuration:")
        call printconfig_t(config)

    end subroutine read_config_file

    subroutine printconfig_t(config)
        implicit none

        integer(4) :: i
        type(config_t) :: config

        call log_message("pecube run mode: " + config%pecube_run_mode)
        if (config%pecube_run_mode == pecube_run_mode_normal) then
            call log_message("normal mode")
        else if (config%pecube_run_mode == pecube_run_mode_error_iter) then
            call log_message("error iteration mode")
            call log_message("error_iter_radius: " + config%error_iter_radius)
            call log_message("error_iter_misfit_limit: " + config%error_iter_misfit_limit)
            call log_message("error_iter_max_topography [m]: " + config%error_iter_max_topography)
        else if (config%pecube_run_mode == pecube_run_mode_monte_carlo) then
            call log_message("monte carlo mode")
            call log_message("mc_min_erosion_rate: " + config%mc_min_erosion_rate)
            call log_message("mc_max_erosion_rate: " + config%mc_max_erosion_rate)
            call log_message("mc_num_of_simulations: " + config%mc_num_of_simulations)
            call log_message("mc_tolerance_chi_squared: " + config%mc_tolerance_chi_squared)
            call log_message("mc_check_min_threshold: " + config%mc_check_min_threshold)
            call log_message("mc_min_threshold_factor: " + config%mc_min_threshold_factor)
            call log_message("mc_erosion_step: " + config%mc_erosion_step)
            call log_message("mc_csv_input_file: " + trim(config%mc_csv_input_file))
        else
            call log_message("unknown run mode:" + config%pecube_run_mode)
        endif

        if (allocated(config%mc_time_slices)) then
            do i = 1, size(config%mc_time_slices)
                call log_message("mc time slice: " + config%mc_time_slices(i)%s_begin + " - " + &
                    config%mc_time_slices(i)%s_end)
            enddo
        endif


        call log_message("RDAAM_grain_radius: " + config%RDAAM_grain_radius)
        call log_message("RDAAM_ppm_Sm: " + config%RDAAM_ppm_Sm)
        call log_message("RDAAM_ppm_Th: " + config%RDAAM_ppm_Th)
        call log_message("RDAAM_ppm_U: " + config%RDAAM_ppm_U)

        return

        ! TODO

        call log_message("output_folder: " + trim(config%output_folder))
        call log_message("topography_file_mode: " + config%topography_file_mode)
        call log_message("topography_type: " + config%topography_type)
        call log_message("topography_file_name: " + trim(config%topography_file_name))
        call log_message("coordinate_system: " + config%coordinate_system)
        call log_message("nx: " + config%nx)
        call log_message("ny: " + config%ny)
        call log_message("spacing_long: " + config%spacing_long)
        call log_message("spacing_lat: " + config%spacing_lat)
        call log_message("nskip: " + config%nskip)
        call log_message("location_long: " + config%location_long)
        call log_message("location_lat: " + config%location_lat)
        call log_message("erosial_time_step: " + config%erosial_time_step)
        call log_message("time step: TODO")
        call log_message("vx_min: " + config%vx_min)
        call log_message("vx_max: " + config%vx_max)
        call log_message("vy_min: " + config%vy_min)
        call log_message("vy_max: " + config%vy_max)
        call log_message("vz_min: " + config%vz_min)
        call log_message("vz_max: " + config%vz_max)
        call log_message("isostacy: " + config%isostacy)
        call log_message("young_modules: " + config%young_modules)
        call log_message("poisson_ratio: " + config%poisson_ratio)
        call log_message("elastic_plate: " + config%elastic_plate)
        call log_message("fft_grid_x: " + config%fft_grid_x)
        call log_message("fft_grid_y: " + config%fft_grid_y)
        call log_message("model_thickness: " + config%model_thickness)
        call log_message("number_of_z_planes: " + config%number_of_z_planes)
        call log_message("specific_heat_capacity: " + config%specific_heat_capacity)
        call log_message("thermal_conductivity: " + config%thermal_conductivity)
        call log_message("specific_heat_capacity: " + config%specific_heat_capacity)
        call log_message("crustal_density: " + config%crustal_density)
        call log_message("mantle_density: " + config%mantle_density)
        call log_message("base_temperature: " + config%base_temperature)
        call log_message("z0_temperature: " + config%z0_temperature)
        call log_message("atmospheric_lapse_rate: " + config%atmospheric_lapse_rate)
        call log_message("crustal_heat_production: " + config%crustal_heat_production)
        call log_message("e_fold_depth: " + config%e_fold_depth)
        call log_message("mantle_heat_production: " + config%mantle_heat_production)
        call log_message("brittle_shear_heating: " + config%brittle_shear_heating)
        call log_message("nepal_model1: " + config%nepal_model1)
        call log_message("nepal_model2: " + config%nepal_model2)
        call log_message("nepal_model3: " + config%nepal_model3)
        call log_message("nepal_model4: " + config%nepal_model4)
        call log_message("thermocron data, TODO")
        call log_message("age_calculation, TODO")
        call log_message("detridal_age: " + config%detridal_age)
        call log_message("min_nodes: " + config%min_nodes)
        call log_message("cascade_out: " + config%cascade_out)
        call log_message("temperature_file: " + config%temperature_file)
        call log_message("2d_move: " + config%move_file)
    end subroutine printconfig_t

end module m_read_config_file
