module m_define_proc
contains

subroutine define_proc (nproc,nsurf,nz,proc)

      integer proc(nsurf*nz)

        do iprc=0,nproc-1
        node1=(iprc*nsurf)/nproc+1
        node2=((iprc+1)*nsurf)/nproc
        if (iprc.eq.nproc) node2=nsurf
          do i=node1,node2
            do k=1,nz
            proc((i-1)*nz+k)=iprc
            enddo
          enddo
        enddo

      return
      end subroutine define_proc
end module m_define_proc

