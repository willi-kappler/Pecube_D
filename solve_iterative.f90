module m_solve_iterative
use m_compiler
use m_logger

contains

subroutine solve_iterative (n,ael,f,x,kfix,icon, &
                                    nnode,nelem,niter, &
                                    znow,tlapse,tmsl,zl)

      implicit none

! Gauss-Seidel iterative solver

! WK: x is temperature here!

      integer(4), intent(in) :: n, nelem, nnode
      real(8), intent(in) :: tlapse, tmsl, zl

      real(8) :: ael(n,n,nelem),f(nnode),x(nnode)
      real(8) :: znow(nnode)
      integer(4) :: icon(n,nelem),kfix(nnode)

      integer,dimension(:),allocatable :: jcmax
      integer,dimension(:,:),allocatable :: jc,jcloc,kc
      real(8),dimension(:),allocatable :: diag, x0
      real(8) :: sum0, sum1, sumref, sup, beta, sumrefall
      real(8) :: tol

      ! real(8), parameter :: min_temp = -100.0_8, max_temp = 10000.0_8
      integer(4) :: k, knode, jelem
      integer(4) :: inode, it, itmax, je, ielemloc, jnode
      integer(4) :: melem, ic, ielem, niter

      !call log_message("solve_iterative.f90: date and time")

! 2012.09.13, WK: check if numbers for temperature are in range

!    do i=1,nnode
!        if ((x(i) < min_temp) .or. (x(i) > max_temp)) then
!            call log_message( "x(i) out of range: ", x(i), i
!            stop
!        endif
!    enddo

    !call log_message("solve_iterative.f90: min temp, max temp: " + minval(x) + "" + maxval(x))
    !call log_message("solve_iterative.f90: min ael, max ael: " + minval(ael) + "" + maxval(ael))
    !call log_message("solve_iterative.f90: min f, max f: " + minval(f) + "" + maxval(f))
    !call log_message("solve_iterative.f90: min znow, max znow: " + minval(znow) + "" + maxval(znow))
    !call log_message("solve_iterative.f90: min kfix, max kfix: " + minval(kfix) + "" + maxval(kfix))
    !call log_message("solve_iterative.f90: min icon, max icon: " + minval(icon) + "" + maxval(icon))

! prepare inverse connectivity arrays

      allocate (jcmax(nnode),diag(nnode),x0(nnode))

      jcmax=0

! 2012.07.15, WK: initialize memory
        diag = 0.0
        x0 = 0.0

        ! WK: can be parallized with OMP
        do jelem=1,nelem
          do inode=1,n
              ic=icon(inode,jelem)
              jcmax(ic)=jcmax(ic)+1
          enddo
        enddo

      melem=maxval(jcmax)

      allocate (jc(melem,nnode),jcloc(melem,nnode),kc(melem,nnode))

! 2012.07.15, WK: initialize memory
      jc = 0
      jcloc = 0
      kc = 0

      jcmax=0
        ! WK: can be parallized with OMP
        do jelem=1,nelem
          do inode=1,n
              ic=icon(inode,jelem)
              jcmax(ic)=jcmax(ic)+1
              jc(jcmax(ic),ic)=jelem
              jcloc(jcmax(ic),ic)=jelem
              kc(jcmax(ic),ic)=inode
          enddo
        enddo

! prepare diagonal

      diag=0.
        ! WK: can be parallized with OMP
        do je=1,nelem
          do k=1,n
              ic=icon(k,je)
              diag(ic)=diag(ic)+ael(k,k,je)
          enddo
        enddo


! start of iterations

      itmax=100000
      tol=1.e-6
      beta=1.3
      x0=x


      it=0
1     it=it+1
      if (it.gt.itmax) then
          print*,'too many iterations'
          stop
      endif


        ! WK: can be parallized with OMP
        do k=1,nnode
          sup=f(k)
          do jelem=1,jcmax(k)
            ielem=jc(jelem,k)
            ielemloc=jcloc(jelem,k)
            knode=kc(jelem,k)
            do jnode=1,n
              ic=icon(jnode,ielem)
              sup=sup-ael(knode,jnode,ielemloc)*x(ic)
            enddo
          enddo
          ! Edited to reflect new code where top and bottom fixed T B/Cs are
          ! identified separately (old version below)
          ! dwhipp (09/07)
          if (kfix(k).eq.0) then
            x(k)=x0(k)+beta*sup/diag(k)
          endif
          if (kfix(k).eq.2) then
            x(k)=(-znow(k)+zl)*tlapse+tmsl
          endif
        enddo


! check for convergence

      sum1=0.
      sum0=0.
        do inode=1,nnode
            sum1=max(sum1,abs(x(inode)))
            sum0=max(sum0,abs((x(inode)-x0(inode))))
        enddo

      sumref=sum0/sum1
      x0=x

      sumrefall=sumref

      if(sumrefall.gt.tol) then
        goto 1
      endif

      deallocate (jcmax,jc,jcloc,kc,diag,x0)

      niter=it

      return
      end subroutine

end module m_solve_iterative
