






c=======================================================================
c///////////////////////  SUBROUTINE ZEUS_MAIN  \\\\\\\\\\\\\\\\\\\\\\\c
      subroutine zeus_main(d, e, u, v, w, p, C1, C2,
     &                  gravity, gr_xacc, gr_yacc, gr_zacc,
     &                  gamma, dt, nhy, dx, dy, dz,
     &                  rank, in, jn, kn, start, end,
     &                  gridvel, iflatten, ipresfree,
     &                  idiff, isteepen, idual, igamfield, eta1, eta2,
     &                  nsubgrids, lface, rface,
     &                  fistart, fiend, fjstart, fjend,
     &                  array, dindex, eindex,
     &                  uindex, vindex, windex, geindex, tmp,
     &                  ncolor, colorpt, coloff, colindex,
     &                  bottom, minsupecoef)
c
c  written by:
c  date:       
c  modified1: Alexei Kritsuk, July 2003; Corrected 2D permutation scheme.
c
c  PURPOSE:
c
c  EXTERNALS:
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     end     - array (of dimension 3) specifying the end of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     gridvel - bulk grid velocity (vector of length 3)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (not currently used)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     bottom  - true (1) if this is the lowest level
c     minsupecoef - coefficient for minimum pressure support (0 - not used)
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer gravity, idiff, idual, iflatten, isteepen, nhy, rank,
     &        ipresfree, end(3), in, jn, kn, nsubgrids, start(3),
     &        ncolor, coloff(ncolor), bottom, igamfield
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real d(in,jn,kn), e(in,jn,kn), u(in,jn,kn),
     &     v(in,jn,kn), w(in,jn,kn), p(in,jn,kn), gamma(in,jn,kn),
     &     gr_xacc(in,jn,kn), gr_yacc(in,jn,kn), gr_zacc(in,jn,kn),
     &     dx(in), dy(jn), dz(kn)
      real dt, eta1, eta2, gridvel(3), pmin, C1, C2
      real array(1), colorpt(1), minsupecoef
c
c  Locals
c
      integer i, ie, is, j, je, js, k, ks, ke, n, ixyz
      real tmp(*)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Error check
c
      if (max(in,jn,kn) .gt. 4103) then
         write(6,*) 'ZEUS_MAIN: A grid dimension is too long.'
         write(6,*) '   (increase max_any_single_direction.)'
         CALL f_error("zeus_main.src",92)
      endif
c
c  Convert arguments to usable form
c
      is = start(1) + 1
      js = start(2) + 1
      ks = start(3) + 1
      ie = end(1) + 1
      je = end(2) + 1
      ke = end(3) + 1
c
c  If DEFAULT_GHOST_ZONES is set to 4, then use the extra space
c
      if (is .eq. 5) then
         is = is - 1
         ie = ie + 1
      endif
      if (js .eq. 5) then
         js = js - 1
         je = je + 1
      endif
      if (ks .eq. 5) then
         ks = ks - 1
         ke = ke + 1
      endif
c
c     Set minimum pressure (better if it were a parameter)
c
      pmin = 1.0e-20
c
      if (rank .eq. 3) then
         do k=ks-3,ke+3
            do j=js-3,je+3
               do i=is-3,ie+3
                  if (abs(u(i,j,k)) .gt. dx(i)/dt .or.
     &                abs(v(i,j,k)) .gt. dy(j)/dt .or.
     &                abs(w(i,j,k)) .gt. dz(k)/dt    ) then
                     write(6,*) 'zpre',i,j,k,ie,je,ke
                     write(6,*) u(i,j,k),v(i,j,k),w(i,j,k)
                     write(6,*) d(i,j,k),e(i,j,k)
                     write(6,*) gr_xacc(i,j,k),gr_yacc(i,j,k),
     &                          gr_zacc(i,j,k)
                     write(6,*) dx(i),dt
                     CALL f_error("zeus_main.src",136)
                  endif
               enddo
            enddo
         enddo
      endif
c     
c     1) Add source terms
c
c      j=(js+je)/2
c      k=(ks+ke)/2
c      write(30,*) 'zeus_a:'
c      do i=is,ie
c         write(30,1040) i,d(i,j,k),u(i,j,k),v(i,j,k),e(i,j,k)
c 1040    format(i4,1p,5(g12.4,1x))
c      enddo
      call zeus_source(d, e, u, v, w, p, in, jn, kn, rank, igamfield,
     &                       is, ie, js, je, ks, ke, C1, C2, ipresfree,
     &                       gamma, dt, pmin, dx, dy, dz,
     &                       gravity, gr_xacc, gr_yacc, gr_zacc, 
     &                       bottom, minsupecoef)
c
c     2) Transport step
c
      ixyz = mod(nhy,rank)
!      ixyz = mod(nhy,3)
      do n=ixyz,ixyz+rank-1
!      do n=ixyz,ixyz+2
c
      if (rank .eq. 3) then
         do k=ks,ke
            do j=js,je
               do i=is,ie
                  if (abs(u(i,j,k)) .gt. dx(i)/dt .or.
     &                abs(v(i,j,k)) .gt. dy(j)/dt .or.
     &                abs(w(i,j,k)) .gt. dz(k)/dt    ) then
                     write(6,*) 'zpost',i,j,k,ie,je,ke,n,ixyz
                     write(6,*) u(i,j,k),v(i,j,k),w(i,j,k)
                     write(6,*) d(i,j,k),e(i,j,k)
                     write(6,*) gr_xacc(i,j,k),gr_yacc(i,j,k),
     &                          gr_zacc(i,j,k)
                     write(6,*) dx(i),dt
                     CALL f_error("zeus_main.src",178)
                  endif
               enddo
            enddo
         enddo
      endif
c
c        Transport step - x direction
c
         if (mod(n,rank) .eq. 0) then
!         if (mod(n,3) .eq. 0)
            call zeus_xtransport(d, e, u, v, w, in, jn, kn, rank,
     &                           is, ie, js, je, ks, ke,
     &                           dt, dx, p, bottom,
     &                           nsubgrids, lface, rface,
     &                           fistart, fiend, fjstart, fjend,
     &                           dindex, eindex, geindex,
     &                           uindex, vindex, windex, array,
     &                           ncolor, colorpt, coloff, colindex)
         endif
c
c        Transport step - y direction
c
         if (mod(n,rank) .eq. 1 .and. rank .gt. 1) then
!         if (mod(n,3) .eq. 2 .and. rank .gt. 1)
            call zeus_ytransport(d, e, u, v, w, in, jn, kn, rank,
     &                           is, ie, js, je, ks, ke,
     &                           dt, dy, p, bottom,
     &                           nsubgrids, lface, rface,
     &                           fistart, fiend, fjstart, fjend,
     &                           dindex, eindex, geindex,
     &                           uindex, vindex, windex, array,
     &                           ncolor, colorpt, coloff, colindex)
         endif
c
c        Transport step - z direction
c
         if (mod(n,rank) .eq. 2 .and. rank .gt. 2) then
!         if (mod(n,3) .eq. 1 .and. rank .gt. 2)
            call zeus_ztransport(d, e, u, v, w, in, jn, kn, rank,
     &                           is, ie, js, je, ks, ke,
     &                           dt, dz, p, bottom,
     &                           nsubgrids, lface, rface,
     &                           fistart, fiend, fjstart, fjend,
     &                           dindex, eindex, geindex,
     &                           uindex, vindex, windex, array,
     &                           ncolor, colorpt, coloff, colindex)
         endif
c
      enddo
c
      return
      end
