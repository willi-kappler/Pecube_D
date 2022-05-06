# SConstruct build file for scons
# Version 0.1 (2014.08.20), written by Willi Kappler

import subprocess
import platform
import os

# Add our own command line option for debug builds
AddOption("--debug-build", action="store_true", dest="debug_build")
AddOption("--use-intel", action="store_true", dest="use_intel")
AddOption("--use-gnu", action="store_true", dest="use_gnu")
AddOption("--use-pgi", action="store_true", dest="use_pgi")
AddOption("--use-mpi", action="store_true", dest="use_mpi")
AddOption("--profile-build", action="store_true", dest="profile_build")

# We append compiler specific files if needed
sourceFiles = []

# Define system environmen with all shell variables.
env = Environment(ENV = os.environ, F90PATH = ".")


# Helper function to set all fortran compiler setting
def set_fortran_compiler(
        env,
        debug_build,
        profile_build,
        compiler,
        compile_flags,
        debug_flags,
        optimize_flags,
        link_flags,
        link_flags_debug,
        profile_flags):

    env.Replace(LINK = compiler)
    env.Replace(F77 = compiler)
    env.Replace(F90 = compiler)
    env.Replace(F95 = compiler)

    env.Append(FORTRANFLAGS = compile_flags)
    env.Append(F77FLAGS = compile_flags)
    env.Append(F90FLAGS = compile_flags)
    env.Append(F95FLAGS = compile_flags)

    if debug_build:
        env.Append(FORTRANFLAGS = debug_flags)
        env.Append(F77FLAGS = debug_flags)
        env.Append(F90FLAGS = debug_flags)
        env.Append(F95FLAGS = debug_flags)
        env.Append(LINKFLAGS = link_flags_debug)
    else:
        env.Append(FORTRANFLAGS = optimize_flags)
        env.Append(F77FLAGS = optimize_flags)
        env.Append(F90FLAGS = optimize_flags)
        env.Append(F95FLAGS = optimize_flags)
        env.Append(LINKFLAGS = link_flags)

    if profile_build:
        env.Append(FORTRANFLAGS = profile_flags)
        env.Append(F77FLAGS = profile_flags)
        env.Append(F90FLAGS = profile_flags)
        env.Append(F95FLAGS = profile_flags)
        env.Append(LINKFLAGS = profile_flags)

    env.Append(LIBS = "stdc++")

    print("FORTRAN: {}".format(env["FORTRAN"]))
    print("LINK: {}, LINKFLAGS: {}, LIBS: {}".format(env["LINK"], env["LINKFLAGS"], env["LIBS"]))
    print("F90FLAGS: {}".format(env["F90FLAGS"]))

# Only change environmemt if we do a real build
if not (GetOption("help") or GetOption("clean")):
    debug_build = GetOption("debug_build")
    use_intel = GetOption("use_intel")
    use_gnu = GetOption("use_gnu")
    use_pgi = GetOption("use_pgi")
    use_mpi = GetOption("use_mpi")
    profile_build = GetOption("profile_build")

    if use_intel:
        env.Replace(FORTRAN = "ifort")
    elif use_gnu:
        env.Replace(FORTRAN = "gfortran")
    elif use_pgi:
        env.Replace(FORTRAN = "pgfortran")
    elif use_mpi:
        env.Replace(FORTRAN = "mpifort")


    if env["CXX"] == "g++":
        print("g++ found, setting options...")

        env.Append(CCFLAGS = "-Wall -Wextra -Werror -Wno-error=maybe-uninitialized")

        if debug_build:
            env.Append(CCFLAGS = "-g3")
        else:
            env.Append(CCFLAGS = "-O3 -g3")

        if profile_build:
            env.Append(CCFLAGS = "-pg")

        print("CCFLAGS: {}".format(env["CCFLAGS"]))

    if env["FORTRAN"] == "ifort":
        print("using the intel fortran compiler, setting options...")

        set_fortran_compiler(
            env,
            debug_build,
            profile_build,
            "ifort",
            "-openmp -warn all",
            "-g -check all -debug all",
            "-fast",
            "-openmp -fast",
            "-openmp -g -check all -debug all",
            "-p"
            )

        sourceFiles.append("intel_compiler.f90")

    elif env["FORTRAN"] == "gfortran":
        print("using the gnu fortran compiler, setting options...")

        set_fortran_compiler(
            env,
            debug_build,
            profile_build,
            "gfortran",
            "-Werror -Wall -Wunderflow -Wno-error=maybe-uninitialized -ffree-line-length-none -fopenmp",
            "-g3 -fbacktrace -fdump-core -fcheck=all",
            "-O3 -g3",
            "-fopenmp -O3 -g3",
            "-fopenmp -g3",
            "-pg"
            )

        sourceFiles.append("gnu_compiler.f90")
        sourceFiles.append("fake_mpi.f90")

    elif env["FORTRAN"] == "pgfortran":
        print("using the portland group fortran compiler, setting options...")

        set_fortran_compiler(
            env,
            debug_build,
            profile_build,
            "pgfortran",
            "",
            "-g",
            "-fastsse -O4",
            "-fastsse -O4",
            "-g",
            "-pg"
            )

        sourceFiles.append("portland_compiler.f90")

    elif env["FORTRAN"] == "mpifort":
        print("using the mpi fortran wrapper, setting options...")

        set_fortran_compiler(
            env,
            debug_build,
            profile_build,
            "mpifort",
            "-Werror -Wall -Wunderflow -Wno-error=maybe-uninitialized -ffree-line-length-none -fopenmp",
            "-g3 -fbacktrace -fdump-core -fcheck=all",
            "-O3 -g3",
            "-fopenmp -O3 -g3",
            "-fopenmp -g3",
            "-pg"
            )

        sourceFiles.append("gnu_compiler.f90")

    currentSystem = platform.system()

    if currentSystem == "Linux":
        print("System is Linux, running updateVersion...")
        subprocess.call("bash updateVersion.sh", shell=True)
    elif currentSystem == "Windows":
        print("System is Windows")
    elif currentSystem == "MacOS":
        print("System is MacOS")

    SetOption('num_jobs', 8)

sourceFiles = sourceFiles + [
    "aft.cc",
    "age_algorithms.f90",
    "ages_header.f90",
    "backtrack_temperature.f90",
    "bivar.f90",
    "borehole_ages.f90",
    "calculate_ages.f90",
    "calculate_misfit.f90",
    "catchments_output.f90",
    "create_pecube_in.f90",
    "define_proc.f90",
    "detrital_mc.f90",
    "data_structures.f90",
    "dynamic_thermal_conductivity.f90",
    "erates.f90",
    "error_iter.f90",
    "export_surface_line.f90",
    "find_dt.f90",
    "find_neighbours.f90",
    "find_temperature.f90",
    "find_upstream_points.f90",
    "find_velo.f90",
    "four1.f90",
    "geometry.f90",
    "getCatchment.f90",
    "global_temperature.f90",
    "global_velocities.f90",
    "import_2d_move_topography.f90",
    "indexx.f90",
    "init_random_seed.f90",
    "interpolate.f90",
    "isostatic_rebound.f90",
    "ketch.c",
    "kptwo.f90",
    "logger.f90",
    "Mad_He.f90",
    "Mad_He2.f90",
    "Mad_Trax_Zircon.f90",
    "Mad_Trax.f90",
    "make_matrix.f90",
    "monte_carlo_erosion.f90",
    "move_points.f90",
    "move_velocities.f90",
    "Muscovite.f90",
    "nr.f90",
    "nrtype.f90",
    "nrutil.f90",
    "pdfmaker_for_cascade.f90",
    "pdfmaker_for_data.f90",
    "pdfmaker_for_pecube.f90",
    "pecube_config.f90",
    "pecube_func.f90",
    "Pecube.f90",
    "RDAAM.cpp",
    "read_config_file.f90",
    "read_time_temperature_history.f90",
    "realft.f90",
    "sinft.f90",
    "solve_iterative.f90",
    "sort2.f90",
    "tapesg_subroutines.f90",
    "tec_mat_output.f90",
    "tridag.f90",
    "version.f90",
    "ZFT.f90"
]

t = env.Program(target="pecube", source=sourceFiles)
Clean(t, ["m_compiler.mod", "gnu_compiler.o", "intel_compiler.o", "mpi.mod", "fake_mpi.o", "m_error_iter2.mod", "m_find_element.mod", "m_interpol3d.mod"])
Default(t)
