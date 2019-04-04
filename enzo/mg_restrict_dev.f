












c=======================================================================
c/////////////////////////  SUBROUTINE MG_RESTRICT  \\\\\\\\\\\\\\\\\\\c
      subroutine mg_restrict(source, dest, ndim, sdim1, sdim2, sdim3,
     &                       ddim1, ddim2, ddim3)
c
c  MULTIGRID: RESTRICT FROM SOURCE TO DEST
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c
c  OUTPUT ARGUMENTS: 
c     dest         - restricted field
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
      integer*8 ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim
      real*8    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      integer*8 i, j, k, i1, j1, k1
      real*8  fact1, fact2, fact3, x, y, z, dxm, dym, dzm, dx0, dy0,
     &        dz0, dxp, dyp, dzp, coef1, coef3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     These coefficients modify the restricted values for 1 and 2D problems
c       to speed up convergence.
c
      coef1 = 2._RKIND
      coef3 = 0.52_RKIND
c
c     Precompute some things
c
      fact1 = REAL(sdim1-1,RKIND)/REAL(ddim1-1,RKIND)
      if (ndim .ge. 2) fact2 = REAL(sdim2-1,RKIND)/REAL(ddim2-1,RKIND)
      if (ndim .ge. 3) fact3 = REAL(sdim3-1,RKIND)/REAL(ddim3-1,RKIND)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=2, ddim1-1
            x = REAL(i-1,RKIND)*fact1 + 0.5_RKIND
            i1 = int(x,IKIND) + 1
            dxm = 0.5_RKIND*(   -x+REAL(i1,RKIND))
            dxp = 0.5_RKIND*(1._RKIND+x-REAL(i1,RKIND))
            dx0 = 1._RKIND - dxp - dxm
            dest(i,1,1) = source(i1-1,1,1)*dxm + 
     &                    source(i1  ,1,1)*dx0 +
     &                    source(i1+1,1,1)*dxp
            dest(i,1,1) = coef1*dest(i,1,1)
         enddo
         dest(    1,1,1) = source(    1,1,1)
         dest(ddim1,1,1) = source(sdim1,1,1)
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=2, ddim2-1
            y = REAL(j-1,RKIND)*fact2 + 0.5_RKIND
            j1 = int(y,IKIND) + 1
            dym = 0.5_RKIND*(   -y+REAL(j1,RKIND))
            dyp = 0.5_RKIND*(1._RKIND+y-REAL(j1,RKIND))
            dy0 = 1._RKIND - dyp - dym
            do i=2, ddim1-1
               x = REAL(i-1,RKIND)*fact1 + 0.5_RKIND
               i1 = int(x,IKIND) + 1
               dxm = 0.5_RKIND*(   -x+REAL(i1,RKIND))
               dxp = 0.5_RKIND*(1._RKIND+x-REAL(i1,RKIND))
               dx0 = 1._RKIND - dxp - dxm
               dest(i,j,1) = source(i1-1,j1-1,1)*dxm*dym + 
     &                       source(i1  ,j1-1,1)*dx0*dym +
     &                       source(i1+1,j1-1,1)*dxp*dym +
     &                       source(i1-1,j1  ,1)*dxm*dy0 + 
     &                       source(i1  ,j1  ,1)*dx0*dy0 +
     &                       source(i1+1,j1  ,1)*dxp*dy0 +
     &                       source(i1-1,j1+1,1)*dxm*dyp + 
     &                       source(i1  ,j1+1,1)*dx0*dyp +
     &                       source(i1+1,j1+1,1)*dxp*dyp
            enddo
            dest(1    ,j,1) = source(1    ,j1,1)
            dest(ddim1,j,1) = source(sdim1,j1,1)
         enddo
         do i=1, ddim1
            i1 = min(max(int(REAL(i-1,RKIND)*fact1 + 0.5_RKIND,IKIND)+1, 
     &           1), sdim1)
            dest(i,    1,1) = source(i1,    1,1)
            dest(i,ddim2,1) = source(i1,sdim2,1)
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=2, ddim3-1
            z = REAL(k-1,RKIND)*fact3 + 0.5_RKIND
            k1 = int(z,IKIND) + 1
            dzm = 0.5_RKIND*(   -z+REAL(k1,RKIND))**2
            dzp = 0.5_RKIND*(1._RKIND+z-REAL(k1,RKIND))**2
            dz0 = 1._RKIND - dzp - dzm
            do j=2, ddim2-1
               y = REAL(j-1,RKIND)*fact2 + 0.5_RKIND
               j1 = int(y,IKIND) + 1
               dym = 0.5_RKIND*(   -y+REAL(j1,RKIND))**2
               dyp = 0.5_RKIND*(1._RKIND+y-REAL(j1,RKIND))**2
               dy0 = 1._RKIND - dyp - dym
               do i=2, ddim1-1
                  x = REAL(i-1,RKIND)*fact1 + 0.5_RKIND
                  i1 = int(x,IKIND) + 1
                  dxm = 0.5_RKIND*(   -x+REAL(i1,RKIND))**2
                  dxp = 0.5_RKIND*(1._RKIND+x-REAL(i1,RKIND))**2
                  dx0 = 1._RKIND - dxp - dxm
                  dest(i,j,k) = source(i1-1,j1-1,k1-1)*dxm*dym*dzm + 
     &                          source(i1  ,j1-1,k1-1)*dx0*dym*dzm +
     &                          source(i1+1,j1-1,k1-1)*dxp*dym*dzm +
     &                          source(i1-1,j1  ,k1-1)*dxm*dy0*dzm + 
     &                          source(i1  ,j1  ,k1-1)*dx0*dy0*dzm +
     &                          source(i1+1,j1  ,k1-1)*dxp*dy0*dzm +
     &                          source(i1-1,j1+1,k1-1)*dxm*dyp*dzm + 
     &                          source(i1  ,j1+1,k1-1)*dx0*dyp*dzm +
     &                          source(i1+1,j1+1,k1-1)*dxp*dyp*dzm +

     &                          source(i1-1,j1-1,k1  )*dxm*dym*dz0 + 
     &                          source(i1  ,j1-1,k1  )*dx0*dym*dz0 +
     &                          source(i1+1,j1-1,k1  )*dxp*dym*dz0 +
     &                          source(i1-1,j1  ,k1  )*dxm*dy0*dz0 + 
     &                          source(i1  ,j1  ,k1  )*dx0*dy0*dz0 +
     &                          source(i1+1,j1  ,k1  )*dxp*dy0*dz0 +
     &                          source(i1-1,j1+1,k1  )*dxm*dyp*dz0 + 
     &                          source(i1  ,j1+1,k1  )*dx0*dyp*dz0 +
     &                          source(i1+1,j1+1,k1  )*dxp*dyp*dz0 +

     &                          source(i1-1,j1-1,k1+1)*dxm*dym*dzp + 
     &                          source(i1  ,j1-1,k1+1)*dx0*dym*dzp +
     &                          source(i1+1,j1-1,k1+1)*dxp*dym*dzp +
     &                          source(i1-1,j1  ,k1+1)*dxm*dy0*dzp + 
     &                          source(i1  ,j1  ,k1+1)*dx0*dy0*dzp +
     &                          source(i1+1,j1  ,k1+1)*dxp*dy0*dzp +
     &                          source(i1-1,j1+1,k1+1)*dxm*dyp*dzp + 
     &                          source(i1  ,j1+1,k1+1)*dx0*dyp*dzp +
     &                          source(i1+1,j1+1,k1+1)*dxp*dyp*dzp
                  dest(i,j,k) = coef3*dest(i,j,k)
               enddo
               dest(1    ,j,k) = source(1    ,j1,k1)
               dest(ddim1,j,k) = source(sdim1,j1,k1)
            enddo
            do i=1, ddim1
               i1 = min(max(int(REAL(i-1,RKIND)*fact1 + 0.5_RKIND,IKIND) 
     &              + 1, 1), sdim1)
               dest(i,    1,k) = source(i1,    1,k1)
               dest(i,ddim2,k) = source(i1,sdim2,k1)
            enddo
         enddo
         do j=1, ddim2
            j1 = min(max(int(REAL(j-1,RKIND)*fact2 + 0.5_RKIND,IKIND)+1, 
     &           1), sdim2)
            do i=1, ddim1
               i1 = min(max(int(REAL(i-1,RKIND)*fact1 + 0.5_RKIND,IKIND) 
     &              + 1, 1), sdim1)
               dest(i,j,    1) = source(i1,j1,    1)
               dest(i,j,ddim3) = source(i1,j1,sdim3)
            enddo
         enddo
      endif
c
      return
      end