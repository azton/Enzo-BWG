












      module enzo_precision


      integer*8     :: ii
      integer*8     :: ll


      real*8        :: sp
      real*8        :: dp
      complex*16    :: cp
      complex*16    :: zp

      integer, parameter :: enzo_int = kind(ii)
      integer, parameter :: enzo_ill = kind(ll)
      integer, parameter :: enzo_fpr = kind(sp)
      integer, parameter :: enzo_fpc = kind(cp)
      integer, parameter :: enzo_xpr = kind(dp)
      integer, parameter :: enzo_xpc = kind(zp)

      end module enzo_precision
