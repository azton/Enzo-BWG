














      subroutine ibm_1d(x, rank, n1, n2, n3, idir)

      implicit none

      integer :: rank, n1, n2, n3, idir
      complex :: x(n1)

      write(0,'("Dummy IBM 1D FFT - error")')
      CALL f_error("essl_3d_fft.src",441)

      return
      end

      subroutine ibm_2d(x, rank, n1, n2, n3, idir)

      implicit none

      integer :: rank, n1, n2, n3, idir
      complex :: x(n1,n2)

      write(0,'("Dummy IBM 2D FFT - error")')
      CALL f_error("essl_3d_fft.src",454)

      return
      end

      subroutine ibm_3d( x, rank, n1, n2, n3, idir )

      implicit none

      integer :: rank, n1, n2, n3, idir
      complex :: x(n1,n2,n3)

      write(0,'("Dummy IBM 3D FFT - error")')
      CALL f_error("essl_3d_fft.src",467)

      return
      end

