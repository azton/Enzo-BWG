












! Complex to complex

! isign = 0   initialize coeffts
! isign = -1  forward normal
! isign = +1  inverse normal
! isign = -2  forward normal in, bit-rev out
! isign = +2  inverse input bit-rev, out normal

! cfft1d( r, n, isign, wsave )
! zfft1d( r, n, isign, wsave )

! r(n)  complex / double complex
! n     integer must be power of 2
! wsave complex / double complex  array((3*n)/2)



      subroutine mkl_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      write(0,'("MKL stride 1 FFT error")')
      call stop_all_cpus

      return
      end

