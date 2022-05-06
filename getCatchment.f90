module m_getCatchment
contains

recursive subroutine getCatchment(UPSTREAMI,UPSTREAMJ,nx,ny,new_ages,x,y,z,&
                                  xmin,ymin,pdf_ages,new_edot,pdf_edot,c,&
                                  contained,xstore,ystore,zstore)

  ! This subroutine uses recursion to compile the upstream points 
  ! It is passed in a center point and finds all immediately surrounding
  ! points that flow into it. The ages and erosion rates of these points are added to
  ! a list. Then, all the surrouding points that flow into the center are checked for
  ! points that flow into them and so on until the entire flow path is found.

  IMPLICIT NONE

  integer i,nx,ny,c
  integer*4 xindex,yindex,xmin,ymin
  integer*4 UPSTREAMI(nx,ny,8),UPSTREAMJ(nx,ny,8)
  real*8 x(nx,ny),y(nx,ny),z(nx,ny),new_ages(nx,ny)
  real*8 pdf_ages(nx*ny),pdf_edot(nx*ny),new_edot(nx,ny)
  logical contained(nx,ny)
  real*8 xstore(nx*ny),ystore(nx*ny),zstore(nx*ny)

  ! If this is the first call to getCatchment then open a new file for the basin,
  ! store the basin outlet point, and write out the basin outlet point
  ! The 'contained' array is used to prevent storing a point more than once; this could happen
  ! when a point is upstream of multiple points, thus being counted two or more times
  if (c.eq.1) then
    xstore(c)=x(xmin,ymin)
    ystore(c)=y(xmin,ymin)
    zstore(c)=z(xmin,ymin)/1000.
    contained(:,:) = .false.
    contained(xmin,ymin) = .true.
    pdf_ages(c) = new_ages(xmin,ymin)
    pdf_edot(c) = new_edot(xmin,ymin)
  endif

  ! Iterates through all 8 points surrounding the center point
  do i=1,8
    if ((UPSTREAMI(xmin,ymin,i).ne.0) .or. (UPSTREAMJ(xmin,ymin,i).ne.0)) then                    ! Checks if the current surrounding point flows into the center point

      xindex = UPSTREAMJ(xmin,ymin,i)
      yindex = UPSTREAMI(xmin,ymin,i)

      if (.not.contained(xindex,yindex)) then                  ! If the point is not already in the list of upstream points
        c = c + 1
        xstore(c)=x(xindex,yindex)
        ystore(c)=y(xindex,yindex)
        zstore(c)=z(xindex,yindex)/1000.

        pdf_ages(c) = new_ages(xindex,yindex)
        pdf_edot(c) = new_edot(xindex,yindex)
        contained(xindex,yindex) = .true.

        ! Recursive call to check upstream point for points upstream of it
        call getCatchment(UPSTREAMI,UPSTREAMJ,nx,ny,new_ages,x,y,z,xindex,yindex,pdf_ages,&
                          new_edot,pdf_edot,c,contained,xstore,ystore,zstore)
      endif
    endif
  enddo

end subroutine getCatchment
end module m_getCatchment

