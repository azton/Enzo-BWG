












c=======================================================================
c//////////////////////////  SUBROUTINE MG_RELAX  \\\\\\\\\\\\\\\\\\\\\c
      subroutine mg_relax(solution, rhs, ndim, dim1, dim2, dim3)
c
c  MULTIGRID: RELAX SOLUTION WITH DIFFERENCED POISSON OPERATOR
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:  Oliver Hahn, 
c  date:       February, 2010
c
c  PURPOSE:
c
c  INPUTS:
c     solution     - solution field
c     rhs          - right hand side
c     dim1-3       - dimensions
c     ndim         - rank of fields
c
c  OUTPUT ARGUMENTS: 
c     solution     - solution field
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer*8 ndim, dim1, dim2, dim3
      real*8    solution(dim1, dim2, dim3), rhs(dim1, dim2, dim3)
c
c  locals
c
      integer*8 i, j, k, ipass, istart, jstart, kstart
      real*8    h1, h2, h3, coef1, coef2, coef3
      
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 2nd order 7-point Laplacian --
c

c
c     Precompute some things
c
      h1 = 1._RKIND/REAL(dim1-1,RKIND)
      if (ndim .ge. 2) h2 = h1/REAL(dim2-1,RKIND)
      if (ndim .ge. 3) h3 = h2/REAL(dim3-1,RKIND)
      coef1 = 1._RKIND/2._RKIND
      coef2 = 1._RKIND/4._RKIND
      coef3 = 1._RKIND/6._RKIND
      
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do ipass=1, 2
            do i=ipass+1, dim1-1, 2
               solution(i,1,1) = coef1*(
     &                           solution(i-1,1,1)+solution(i+1,1,1) -
     &                           h1*rhs(i,1,1))
            enddo
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         jstart = 1
         do ipass=1, 2
            istart = jstart
            do j=2, dim2-1
               do i=istart+1, dim1-1, 2
                  solution(i,j,1) = coef2*(
     &                    solution(i-1,j  ,1) + solution(i+1,j  ,1) +
     &                    solution(i  ,j-1,1) + solution(i  ,j+1,1) -
     &                    h2*rhs(i,j,1))
               enddo
               istart = 3-istart
            enddo
            jstart = 3-jstart
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         kstart = 1
! Cache aware Gauss-Seidel method: Data is passed through cache only once.         
         ! Red cells (even-even, odd-odd) on the first slab
         jstart = kstart
         istart = jstart
         do j = 2, dim2-1
            do i = istart+1, dim1-1, 2
               solution(i,j,2) = coef3*(
     &              solution(i-1,j  ,2) + solution(i+1,j  ,2) +
     &              solution(i  ,j-1,2) + solution(i  ,j+1,2) +
     &              solution(i  ,j  ,1) + solution(i  ,j  ,3) -
     &              h3*rhs(i,j,2))
            enddo
            istart = 3-istart
         enddo
         jstart = 3-jstart

         ! Fused loop with red/black cells
         do k = 3, dim3-1
            istart = jstart
            do j = 2, dim2-1
               do i = istart+1, dim1-1, 2
                  ! Red on this slab
                  solution(i,j,k) = coef3*(
     &                 solution(i-1,j  ,k  ) + solution(i+1,j  ,k  ) +
     &                 solution(i  ,j-1,k  ) + solution(i  ,j+1,k  ) +
     &                 solution(i  ,j  ,k-1) + solution(i  ,j  ,k+1) -
     &                 h3*rhs(i,j,k))
                  ! Black on previous slab
                  solution(i,j,k-1) = coef3*(
     &                 solution(i-1,j,k-1) + solution(i+1,j,k-1) +
     &                 solution(i,j-1,k-1) + solution(i,j+1,k-1) +
     &                 solution(i,j  ,k-2) + solution(i,j  ,k) -
     &                 h3*rhs(i,j,k-1))
               enddo
               istart = 3-istart
            enddo
            jstart = 3-jstart
         enddo

         ! Black cells (even-odd, odd-even) on the last slab
         istart = jstart
         do j = 2, dim2-1
            do i = istart+1, dim1-1, 2
               solution(i,j,dim3-1) = coef3*(
     &            solution(i-1,j,dim3-1) + solution(i+1,j,dim3-1)+
     &            solution(i,j-1,dim3-1) + solution(i,j+1,dim3-1)+
     &            solution(i,j  ,dim3-2) + solution(i,j  ,dim3) -
     &            h3*rhs(i,j,dim3-1))
            enddo
            istart = 3-istart
         enddo
      endif
c
      return
      end