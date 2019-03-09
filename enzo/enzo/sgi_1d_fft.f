

      subroutine sgi_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("SGI stride 1 FFT error")')
      CALL f_error("sgi_1d_fft.src",105)
      
      return
      end
      
