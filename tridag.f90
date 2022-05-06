module m_tridag
contains

      SUBROUTINE tridag (a,b,c,r,u,n)

      ! (C) Copr. 1986-92 Numerical Recipes Software

      implicit none

      INTEGER n
      REAL(8) :: a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      REAL(8) :: bet
      REAL(8), dimension(:), allocatable :: gam

      allocate (gam(n))

      if(b(1).eq.0.) stop 'in tridag'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
          print*,'tridag failed'
          stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      deallocate (gam)

      return

      END subroutine

! SUBROUTINE tridag (a,b,c,r,u,n,k,i)
! 
!   ! (C) Copr. 1986-92 Numerical Recipes Software
! 
!   implicit none
! 
!   INTEGER n,k,i
!   REAL a(i,n),b(i,n),c(i,n),r(i,n),u(i,n)
!   INTEGER j
!   REAL bet
!   REAL,dimension(:),allocatable :: gam
!   allocate (gam(n))
!   if(b(k,1).eq.0.) stop 'in tridag'
!   bet=b(k,1)
!   u(k,1)=r(k,1)/bet
!   do 11 j=2,n
!     gam(j)=c(k,j-1)/bet
!     bet=b(k,j)-a(k,j)*gam(j)
!     if(bet.eq.0.) then
!       print*,'tridag failed'
!       stop
!     endif
!     u(k,j)=(r(k,j)-a(k,j)*u(k,j-1))/bet
! 11 continue
!   do 12 j=n-1,1,-1
!     u(k,j)=u(k,j)-gam(j+1)*u(k,j+1)
! 12 continue
!   deallocate (gam)
! 
!   return
! 
! END


end module m_tridag

