













      subroutine ibm_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("IBM stride 1 FFT error")')
      call stop_all_cpus

      return
      end

