

      subroutine ibm_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("IBM stride 1 FFT error")')
      CALL f_error("essl_1d_fft.src",148)

      return
      end

