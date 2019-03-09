












      subroutine nr_1d(x, rank, n1, n2, n3, dir)

      implicit none

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1)

      integer :: n(3)
      real :: factor

!     write(*,*) 'NR_1D ',rank,n1,n2,n3,dir

      if( rank /= 1 ) then
        write(0,*) 'NR_1D rank != 1'
        call stop_all_cpus
      end if

      if( n2 /= 1 ) then
        write(0,*) 'NR_1D dim2 > 1'
        call stop_all_cpus
      end if

      if( n3 /= 1 ) then
        write(0,*) 'NR_1D dim3 > 1'
        call stop_all_cpus
      end if

      n(1) = n1
      n(2) = n2
      n(3) = n3

      factor = 1.0/real(n1)

      if( dir == -1 ) then
        call fourn(x, n, rank, dir)
      else
        call fourn(x, n, rank, dir)
        x = x * factor
      end if

      return
      end


      subroutine nr_2d(x, rank, n1, n2, n3, dir)

      implicit none

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1,n2)

      integer :: n(3)
      real :: factor

!     write(*,*) 'NR_2D ',rank,n1,n2,n3,dir

      if( rank /= 2 ) then
        write(0,*) 'NR_2D rank != 2'
        call stop_all_cpus
      end if

      if( n3 /= 1 ) then
        write(0,*) 'NR_2D dim3 > 1'
        call stop_all_cpus
      end if

      n(1) = n1
      n(2) = n2
      n(3) = n3

      factor = 1.0/real(n1*n2)

      if( dir == -1 ) then
        call fourn(x, n, rank, dir)
      else
        call fourn(x, n, rank, dir)
        x = x * factor
      end if

      return
      end


      subroutine nr_3d(x, rank, n1, n2, n3, dir)

      implicit none

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1,n2,n3)

      integer :: n(3)
      real :: factor

!     write(*,*) 'NR_3D ',rank,n1,n2,n3,dir

      if( rank /= 3 ) then
        write(0,*) 'NR_3D rank != 3'
        call stop_all_cpus
      end if

      n(1) = n1
      n(2) = n2
      n(3) = n3

      factor = 1.0/real(n1*n2*n3)

      if( dir == -1 ) then
        call fourn(x, n, rank, dir)
      else
        call fourn(x, n, rank, dir)
        x = x * factor
      end if

      return
      end
