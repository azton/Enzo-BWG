



















c=======================================================================
c//////////////////////  SUBROUTINE COOL_TIME \ \\\\\\\\\\\\\\\\\\\\\\\c
      subroutine cool_time(
     &                d, e, ge, u, v, w, cooltime,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, idim,
     &                is, js, ks, ie, je, ke, 
     &                dt, aye, temstart, temend, fh, utem,
     &                eta1, eta2, gamma, coola)
c
c  CALCULATE COOLING TIME IN CODE UNITS (NON-MULTI SPECIES VERSION)
c
c  written by: Greg Bryan
c  date:       March, 1997
c  modified1:
c
c  PURPOSE:
c    returns the cooling time in code units on the grid
c
c  INPUTS:
c    is,ie   - start and end indicies of active region (zero-based!)
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  Arguments
c
      integer in, jn, kn, is, js, ks, ie, je, ke,
     &        idual, iexpand, nratec, idim, imethod
      real    dt, aye, temstart, temend, fh, utem,
     &        eta1, eta2, gamma, coola(nratec)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        cooltime(in,jn,kn)
c
c  Parameters
c
      integer ijk
      parameter (ijk = 4103)
c
c  Locals
c
      integer i, j, k
      real    energy
c
c  Slice locals
c 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), edot(ijk), tgas(ijk), tgasold(ijk),
     &     cool(ijk)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Loop over slices (in the k-direction)
c
      do k = ks+1, ke+1
       do j = js+1, je+1
c
c        Set time elapsed to zero for each cell (and convert dens to abs)
c
         do i = is+1, ie+1
            d(i,j,k) = d(i,j,k)/aye**3
         enddo
c
c        Compute the cooling rate on a slice
c
         call cool1d(d, e, ge, u, v, w,
     &               in, jn, kn, nratec, idual, idim, imethod, 1,
     &               is, ie, j, k, 
     &               temstart, temend, fh, utem,
     &               eta1, eta2, gamma, coola,
     &               indixe, t1, t2, logtem, tdef, edot,
     &               tgas, tgasold, p2d, cool)
c
c        Compute the cooling time on the slice
c
         do i = is+1, ie+1
            energy = max(p2d(i)/d(i,j,k)/(gamma-1.0), 1.0e-20)
            cooltime(i,j,k) = abs(energy/edot(i))
         enddo
c
c       Set the density back to comoving
c
         do i = is+1, ie+1
            d(i,j,k) = d(i,j,k)*aye**3
         enddo
c
c     Next slice
c
        enddo
      enddo
c
      return
      end

