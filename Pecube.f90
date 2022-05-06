! Pecube-D   - University of Tuebingen, Germany, Version
!
!
! Todd Ehlers, May 2015.  todd.ehlers@uni-tuebingen.de
!
! This version of pecube is based on the original version by Jean Braun.  However, it
! it includes many changes that have been implemented by Todd Ehlers work group over the
! years.  These changes include:
! a.  Calculation of predicted ages for different thermo- geochronometer systems.
! b.  Many different options for user defined velocity input fields (e.g. McQuarrie and
!     Ehlers, 2015).
! c.  Different output format options.
! d.  calculation of detrital cooling ages for user defined sample points on the topog.
!      (e.g. Whipp et al., 2009)
! e.  Iterative Inversion of cooling ages for topographic change scenarios (e.g. Olen et al. 2012)
! f.  Monte Carlo Inversion of cooling ages to identify the denudation histories that
!     that can produce observed ages.  (e.g. Thiede & Ehlers, 2013)


! Please quote the following reference(s) when publishing results obtained
! with this versionPecube:

! Braun , J., 2003. Pecube: A new finite element code to solve the 3D heat
!  transport equation including the effects of a time-varying, finite
!  amplitude surface topography.  Computers and Geosciences, v.29, pp.787-794.

! Braun , J., 2002. Quantifying the effect of recent relief changes on age-elevation
!  relationships.  Earth and Planetary Science Letters, v.200, pp.331-343.

! Braun, J., 2002. Estimating exhumation rate and relief evolution by spectral
!  analysis of age-elevation datasets. Terra Nova, v.14, pp.210-214.

! Reference to use concerning program changes made by T. Ehlers group:

! Olen, S., Ehlers, T.A., Densmore, M.S., 2012, Limits to reconstructing paleotopography
! from ther- mochronometer data, J. Geophysical Res – Earth Surface, v. 117,
! doi:10.1029/2011/ JF001985

! Whipp, D.M. Jr., Ehlers, T.A., Braun, J., Spath, C.D., 2009, Effects of exhumation
! kinematics and topo- graphic evolution on detrital thermochronometer data,
! J. Geophysical Res. – Earth Surface, V. 114, F04021, doi:10.1029/2008JF001195.

! Thiede, R.C., Ehlers, T.A., 2013, Large spatial and temporal variations in Himalayan
! denudation, Earth and Planetary Science Letters, 371-372, pp. 278-293.

! McQuarrie, N., and Ehlers, T.A., 2015 Influence of thrust belt geometry
! and shortening rate on thermochronometer cooling ages: Insights from the
! Bhutan Himalaya, Tectonics. 34, doi:10.1002/2014TC003783.





! References for the age calculation (provided by Paul Eizenhöfer):
!
! AHe: Farley, K. A. (2000). Helium diffusion from apatite: General behavior as illustrated by Durango fluorapatite. Journal of Geophysical Research: Solid Earth, 105(B2), 2903-2914.
!
! AFT: Crowley, K. D., Cameron, M., & Schaefer, R. L. (1991). Experimental studies of annealing of etched fission tracks in fluorapatite. Geochimica et Cosmochimica Acta, 55(5), 1449-1465.
!
! ZHe: Reiners, P. W., Spell, T. L., Nicolescu, S., & Zanetti, K. A. (2004). Zircon (U-Th)/He thermochronometry: He diffusion and comparisons with 40Ar/39Ar dating. Geochimica et cosmochimica acta, 68(8), 1857-1887.
!
! ZFT: Brandon, M. T., Roden-Tice, M. K., & Garver, J. I. (1998). Late Cenozoic exhumation of the Cascadia accretionary wedge in the Olympic Mountains, northwest Washington State. Geological Society of America Bulletin, 110(8), 985-1009.
!
! MAr: Hames, W. E., & Bowring, S. A. (1994). An empirical evaluation of the argon diffusion geometry in muscovite. Earth and Planetary Science Letters, 124(1-4), 161-169.




program Pecube
    use m_read_config_file
    use m_pecube_func
    use m_error_iter
    use m_monte_carlo_erosion
    use m_logger
    use m_pecube_config
    use m_compiler
    use m_version
    use m_data_structures

    use mpi

    implicit none

    type(config_t) :: config

    integer(4) :: pecube_mpi_error, return_code
    logical :: file_exists

    num_threads = 8

    call MPI_INIT(pecube_mpi_error)

    if (pecube_mpi_error /= MPI_SUCCESS) then
        print *, "error: could not set up MPI"
        print *, "error code: 5f1cd18df3dcf57ebd5455988d8aa140"
        return_code = pecube_mpi_error
        call MPI_ABORT(MPI_COMM_WORLD, return_code, pecube_mpi_error)
        stop
    endif

    call MPI_COMM_SIZE(MPI_COMM_WORLD, config%mpi_total_num_of_cpu, pecube_mpi_error)

    if (pecube_mpi_error /= MPI_SUCCESS) then
        print *, "error: could not get the total number of cpus"
        print *, "error code: ea2b6e161636da729750b50c4b27d956"
        return_code = pecube_mpi_error
        call MPI_ABORT(MPI_COMM_WORLD, return_code, pecube_mpi_error)
    endif

    call MPI_COMM_RANK(MPI_COMM_WORLD, config%mpi_current_cpu_id, pecube_mpi_error)

    if (pecube_mpi_error /= MPI_SUCCESS) then
        print *, "error: could not get the current cpu id"
        print *, "error code: cd74ad0d7f502754adb36397bf8becdb"
        return_code = pecube_mpi_error
        call MPI_ABORT(MPI_COMM_WORLD, return_code, pecube_mpi_error)
    endif

    call logger_init(config)

    call log_message('Begin Pecube Execution')
    call print_program_version()

    call log_message("number of mpi cpus: " + config%mpi_total_num_of_cpu)
    call log_message("current mpi id: " + config%mpi_current_cpu_id)
    call log_message("log file name: " + trim(config%log_filename))
    call log_message("System information:")
    call sys_command("date >> " // trim(config%log_filename))
    call sys_command("pwd >> " // trim(config%log_filename))
    call sys_command("whoami >> " // trim(config%log_filename))
    call sys_command("uname -a >> " // trim(config%log_filename))
    call log_message("End system information")

    inquire (file="Pecube.in", exist = file_exists)

    if (.not. file_exists) then
        call log_message("The file 'Pecube.in' does not exist! Exit now")
        stop
    endif

    call read_config_file("Pecube.in", config)

    config%mc_use_monte_carlo = .false.

    select case (config%pecube_run_mode)
        case (pecube_run_mode_normal)
            call log_message("normal pecube simulation")
            call pecube_func(config)
        case (pecube_run_mode_error_iter)
            call log_message("error iteration simulation")
            call error_iter(config)
        case (pecube_run_mode_monte_carlo)
            call log_message("monte carlo simulation")
            config%mc_use_monte_carlo = .true.
            call monte_carlo_erosion(config)
        case default
            call log_message("unknown run mode")
    end select

    ! clean up allocated memory
    if (allocated(topo_file_name)) then
        deallocate(topo_file_name)
    endif

    if (allocated(surface_age_info)) then
      call deallocate_age_info(surface_age_info, nstep)
    endif

    if (allocated(age_flags)) then
        deallocate (age_flags)
    endif

    call logger_finish()

    call MPI_Finalize(pecube_mpi_error)

    ! end of pecube
end program pecube
