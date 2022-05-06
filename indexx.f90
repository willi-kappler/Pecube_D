module m_indexx
contains

! Subroutine indexx.f90 from Numerical Recipes in Fortran 90
! Slightly modified to pass in the arr array size
!
! dwhipp - 04/08

      SUBROUTINE indexx(arr,n1,index)
      USE nrutil, ONLY : arth,assert_eq,nrerror,swap
      IMPLICIT NONE
      REAL*8,DIMENSION(n1),INTENT(IN) :: arr
      INTEGER,DIMENSION(n1),INTENT(OUT) :: index
      INTEGER,PARAMETER :: NN=15, NSTACK=50
      REAL*8 :: a
      INTEGER :: n,k,i,j,indext,jstack,l,r,n1
      INTEGER,DIMENSION(NSTACK) :: istack
      n=assert_eq(size(index),size(arr),'indexx')
      index=arth(1,1,n)
      jstack=0
      l=1
      r=n
      do
        if (r-l < NN) then
          do j=l+1,r
            indext=index(j)
            a=arr(indext)
            do i=j-1,l,-1
              if (arr(index(i)) <= a) exit
              index(i+1)=index(i)
            end do
            index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
            do
              i=i+1
              if (arr(index(i)) >= a) exit
            end do
            do
              j=j-1
              if (arr(index(j)) <= a) exit
            end do
            if (j < i) exit
            call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          end if
        end if
      end do
      CONTAINS
!BL
      SUBROUTINE icomp_xchg(i,j)
      INTEGER,INTENT(INOUT) :: i,j
      INTEGER :: swp
      if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
      end if
      END SUBROUTINE icomp_xchg
      END SUBROUTINE indexx


      SUBROUTINE indexx_i4b(iarr,index)
      USE nrutil, ONLY : arth,assert_eq,nrerror,swap
      IMPLICIT NONE
      INTEGER,DIMENSION(:),INTENT(IN) :: iarr
      INTEGER,DIMENSION(:),INTENT(OUT) :: index
      INTEGER,PARAMETER :: NN=15, NSTACK=50
      INTEGER :: a
      INTEGER :: n,k,i,j,indext,jstack,l,r
      INTEGER,DIMENSION(NSTACK) :: istack
      n=assert_eq(size(index),size(iarr),'indexx_sp')
      index=arth(1,1,n)
      jstack=0
      l=1
      r=n
      do
        if (r-l < NN) then
          do j=l+1,r
            indext=index(j)
            a=iarr(indext)
            do i=j-1,1,-1
              if (iarr(index(i)) <= a) exit
              index(i+1)=index(i)
            end do
            index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=iarr(indext)
          do
            do
              i=i+1
              if (iarr(index(i)) >= a) exit
            end do
            do
              j=j-1
              if (iarr(index(j)) <= a) exit
            end do
            if (j < i) exit
            call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          end if
        end if
      end do
      CONTAINS
!BL
      SUBROUTINE icomp_xchg(i,j)
      INTEGER,INTENT(INOUT) :: i,j
      INTEGER :: swp
      if (iarr(j) < iarr(i)) then
        swp=i
        i=j
        j=swp
      end if
      END SUBROUTINE icomp_xchg
      END SUBROUTINE indexx_i4b

end module m_indexx

