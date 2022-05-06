module m_compiler
    use ifport

    implicit none

    !> make everything private
    private

    !> Exported subroutines:
    public sys_command
    public sys_chdir

    contains

        subroutine sys_command(cmd)
            implicit none

            character(*) :: cmd
            integer(4) :: res

            res = system(cmd)

        end subroutine sys_command

        subroutine sys_chdir(dir)
            implicit none

            character(*) :: dir
            integer(4) :: res

            res = chdir(dir)

        end subroutine sys_chdir
end module m_compiler
