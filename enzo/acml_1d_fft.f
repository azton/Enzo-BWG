














      subroutine acml_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("ACML stride 1 FFT error")')
      CALL f_error("acml_1d_fft.src",102)

      return
      end

