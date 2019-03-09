







      subroutine open_mpi_error_file( fn, iunit, type )

      implicit none


      character*(*) fn, type
      integer       iunit

      character*4   post

      integer*4     id, ierr
      integer       i
      character*32  fnx

      id = 0

      write(post,'(i4)') id
      do i=1,4
      if(post(i:i).eq.' ') post(i:i)='0'
      end do
!     write(*,'(i4,4x,a4)') id,post

      fnx=fn // '_' // post

      open(unit=iunit,file=fnx,status=type,position='append')

      return
      end

      subroutine close_mpi_error_file( iunit )

      implicit none

      integer iunit

      close(unit=iunit)

      return
      end

