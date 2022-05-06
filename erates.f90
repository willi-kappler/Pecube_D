module m_erates
contains

!       This subroutine generates maps of the surface erosion rates for each
!       model time step (nstep)
!       dwhipp - 02/07

        subroutine erates (nsurf,xxx,xlonmin,xlonmax,xsurf,xmin,xmax,yyy,&
                           xlatmin,xlatmax,ysurf,ymin,ymax,zsurf,run,edot,file_id1,&
                           bg_edot,topo_edot,nrun,nx,nelemsurf)

    real*8 xsurf(nsurf),ysurf(nsurf),zsurf(nsurf)
    real(8) edot(nsurf),bg_edot(nsurf),topo_edot(nsurf)
    real*8 xlonmin,xlonmax,xmin,xmax,xlatmin,xlatmax,xxx,yyy,ymin,ymax
    character run*100,file_id1*4
    integer nelemsurf,counter

        open (68,file=run(1:nrun)//'/erates_tec'//file_id1//'.dat',status='unknown')

!   Writes Tecplot header in Ages_tec.dat file
        write (68,*) 'TITLE = "Pecube Erosion Rates"'
        write (68,'(a126)') 'VARIABLES = "x (km)" "y (km)" "z (km)" "Total Erosion rate &
                      & (mm/yr)" "Background Erosion Rate (mm/yr)" "Relief Change (mm/yr)"'
        write (68,*) 'ZONE T="Erosion Rates"'
        write(68,'(A2,i10)',advance="no") 'n=',nsurf
        write(68,'(A4,i10)',advance="no") ', e=',nelemsurf
        write(68,*) ', et=quadrilateral, f=fepoint'
!         write (68,*) 'I=',nsurf,', J=1, K=1, ZONETYPE=Ordered'
!         write (68,*) 'DATAPACKING=POINT'
!         write (68,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'

!   Loop to write latitude,longitude,elevation,AHe age, and AFT age
!   This loop is taken directly from the Pecube.f90 file
        do i=1,nsurf
          xxx=xlonmin+(xlonmax-xlonmin)*(xsurf(i)-xmin)/(xmax-xmin)
          yyy=xlatmin+(xlatmax-xlatmin)*(ysurf(i)-ymin)/(ymax-ymin)
          write (68,'(6f12.4)') xxx,yyy,zsurf(i),edot(i),bg_edot(i),topo_edot(i)
        enddo

      counter=0
      do i=1,nelemsurf
        counter=counter+1
        if (mod(counter,nx).eq.0) counter=counter+1
        write (68,*) counter+nx,counter+nx+1,counter+1,counter
      enddo


! close file
        close(68)

        return
        end subroutine erates
end module m_erates

