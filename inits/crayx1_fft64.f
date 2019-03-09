













      subroutine cray_1d(x, rank, n1, n2, n3, idir)

      implicit none

      integer :: rank, n1, n2, n3, idir
      complex :: x(n1)

      write(0,'("Dummy Cray X1 1D FFT - error")')
      call stop_all_cpus

      return
      end

      subroutine cray_2d(x, rank, n1, n2, n3, idir)

      implicit none

      integer :: rank, n1, n2, n3, idir
      complex :: x(n1,n2)

      write(0,'("Dummy Cray X1 2D FFT - error")')
      call stop_all_cpus

      return
      end

      subroutine cray_3d( x, rank, n1, n2, n3, idir )

      implicit none

      integer :: rank, n1, n2, n3, idir
      complex :: x(n1,n2,n3)

      write(0,'("Dummy Cray X1 3D FFT - error")')
      call stop_all_cpus

      return
      end

