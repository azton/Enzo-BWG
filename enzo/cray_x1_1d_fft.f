














      subroutine cray_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("Dummy Cray X1 1D FFT - error")')
      CALL f_error("cray_x1_1d_fft.src",91)

      return
      end

