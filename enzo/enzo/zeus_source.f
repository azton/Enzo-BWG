






c=======================================================================
c//////////////////////  SUBROUTINE ZEUS_SOURCE  \\\\\\\\\\\\\\\\\\\\\\c
      subroutine zeus_source(d, e, u, v, w, p, in, jn, kn, rank, 
     &                          igamfield,
     &                       is, ie, js, je, ks, ke, C1, C2, ipresfree,
     &                       gamma, dt, pmin, dx, dy, dz,
     &                       gravity, gr_xacc, gr_yacc, gr_zacc, 
     &                       bottom, minsupecoef)
c
c  SOURCE TERMS FOR ZEUS HYDRO (CARTESIAN ONLY)
c
c  written by: Greg Bryan (implemented from Stone & Norman, ApJS 80, 753)
c  date:       February, 1997
c  modified1: 
c
c  PURPOSE:
c     Adds the source substeps
c
c  EXTERNALS:
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     i,j,kn  - dimensions of field arrays
c     igamfield - indicates if gamma should be a field
c     ipresfree - pressure-free flag (0 = off, 1 = on, i.e. p=0)
c     rank    - dimension of problem (not currently used)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     C1,C2   - Linear and quadratic artifificla viscosity parameters
c     bottom  - true (1) if this is the lowest level
c     minsupecoef - coefficient for minimum pressure support
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
c-----------------------------------------------------------------------
c
c     Arguments
c
      integer in, jn, kn, rank, is, ie, js, je, ks, ke, gravity, 
     &        bottom, igamfield, ipresfree, nstep
      real    d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn),
     &        w(in,jn,kn), p(in,jn,kn), dx(in), dy(jn), dz(kn),
     &        gr_xacc(in,jn,kn), gr_yacc(in,jn,kn), gr_zacc(in,jn,kn)
      real    gamma(in,jn,kn), dt, pmin, C1, C2, C3, minsupecoef
c
c     Parameters
c
      integer ijk
      parameter (ijk = 4103)
c
c     Locals
c
      integer i, j, k, jsm1, ksm1, jep1, kep1, jsm2, ksm2, jep2, kep2, n
      real    alpha, q(ijk), div(ijk), deltav, deltavmax, e1, gamma1
      real    dt1, mom(ijk), momf(ijk), dt2, pfparam, u1, C1a
      integer nmax
      parameter (pfparam = 0.001, nmax=200)
c
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c     Compute varients on start indexes
c
      jsm1 = max(js-1, 1)
      ksm1 = max(ks-1, 1)
      jep1 = min(je+1, jn)
      kep1 = min(ke+1, kn)
      jsm2 = max(js-2, 1)
      ksm2 = max(ks-2, 1)
      jep2 = min(je+2, jn)
      kep2 = min(ke+2, kn)
c
      gamma1 = gamma(1,1,1)  ! if gamma is a scalar
c
c     Compute the pressure
c
      do k = 1, kn
         do j = 1, jn
            if (ipresfree .eq. 1) then
c
c              Pressurefree - set pressure field to zero
c
               do i = 1, in
                  p(i,j,k) = 0
               enddo
            else
c
               if (igamfield .eq. 1) then
c
c                 Compute pressure with variable gamma
c
                  do i = 1, in
                     e1 = e(i,j,k)
                     e1 = max(e(i,j,k), minsupecoef*d(i,j,k))
                     p(i,j,k) = max((gamma(i,j,k)-1.0)*
     &                              d(i,j,k)*e1, pmin)
                  enddo
               else
c
c                 Compute pressure with constant gamma
c
                  do i = 1, in
                     e1 = e(i,j,k)
                     e1 = max(e(i,j,k), minsupecoef*d(i,j,k))
                     p(i,j,k) = max((gamma1-1.0)*d(i,j,k)*e1, pmin)
                     if (e(i,j,k) .le. 0.0 .or. d(i,j,k) .le. 0.0) then
                        write(6,*) 'zeus_source1:',
     &                             e(i,j,k),d(i,j,k),i,j,k
                        CALL f_error("zeus_source.src",122)
                     endif
                  enddo
               endif
c
            endif
c
c        Next j,k
c
         enddo
      enddo
c
c
c     1) Substep 1 -- pressure and gravity terms
c
      do k = ksm2, kn
         do j = jsm2, jn
c
c           Update velocities with compression term
c              (limit increase to preventy pressure driven instability)
c
            do i = is-2, ie+3
               deltav =         dt*(p(i-1,j,k)-p(i,j,k))/
     &                  (dx(i)*0.5*(d(i-1,j,k)+d(i,j,k)))
c               if (abs(deltav) .gt. vlim*dx(i)/dt) write(6,*)
c     &            'z_lim_i',deltav,dx(i),dt,p(i-1,j,k),p(i,j,k),
c     &            i,j,k,ie,je,ke,d(i-1,j,k),d(i,j,k),e(i-1,j,k),
c     &            e(i,j,k)
c               deltav = sign(min(abs(deltav),vlim*dx(i)/dt), deltav)
               u(i,j,k) = u(i,j,k) + deltav
            enddo
            if (rank .gt. 1) then
              do i = is-2, ie+3
                 deltav =         dt*(p(i,j-1,k)-p(i,j,k))/
     &                    (dy(j)*0.5*(d(i,j-1,k)+d(i,j,k)))
c               if (abs(deltav) .gt. vlim*dy(j)/dt) write(6,*)
c     &            'z_lim_j',deltav,dy(j),dt,p(i,j-1,k),p(i,j,k),
c     &            i,j,k,ie,je,ke,d(i-1,j,k),d(i,j,k),e(i,j-1,k),
c     &            e(i,j,k)
c                 deltav = sign(min(abs(deltav),vlim*dy(j)/dt), deltav)
                 v(i,j,k) = v(i,j,k) + deltav
              enddo
            endif
            if (rank .gt. 2) then
              do i = is-2, ie+3
                 deltav =         dt*(p(i,j,k-1)-p(i,j,k))/
     &                    (dz(k)*0.5*(d(i,j,k-1)+d(i,j,k)))
c               if (abs(deltav) .gt. vlim*dz(k)/dt) write(6,*)
c     &            'z_lim_k',deltav,dz(k),dt,p(i,j,k-1),p(i,j,k),
c     &            i,j,k,ie,je,ke,d(i,j,k-1),d(i,j,k),e(i,j,k-1),
c     &            e(i,j,k)
c                 deltav = sign(min(abs(deltav),vlim*dz(k)/dt), deltav)
                 w(i,j,k) = w(i,j,k) + deltav
              enddo
            endif
c
c           Update velocities with acceleration
c
            if (gravity .eq. 1) then
               do i = is-2, ie+3
                  u(i,j,k) = u(i,j,k) + dt*gr_xacc(i,j,k)
               enddo
               if (rank .gt. 1) then
                  do i = is-2, ie+3
                     v(i,j,k) = v(i,j,k) + dt*gr_yacc(i,j,k)
                  enddo
               endif
               if (rank .gt. 2) then
                  do i = is-2, ie+3
                     w(i,j,k) = w(i,j,k) + dt*gr_zacc(i,j,k)
                  enddo
               endif
            endif
c
         enddo
      enddo
c
c     2) Substep 2 -- artificial viscosity
c
      nstep = 5
      dt1 = dt/real(nstep)
      C3 = -1.0
      C1a = C1
      if (bottom .eq. 1) C1a = C1*5.0

c
      do n = 1, nstep
      do k = ksm2, kn
         do j = jsm2, jn
c
c           a) Quadratic viscosity
c
            do i = is-2, ie+2
               if ((u(i+1,j,k)-u(i,j,k)) .lt. 0.0) then
                  q(i) = C2*d(i,j,k)*(u(i+1,j,k)-u(i,j,k))**2
               else
                  q(i) = 0.0
               endif
            enddo
c
c           b) linear viscosity
c
            if (C1 .ne. 0.0) then
               if (igamfield .eq. 1) then
                  do i = is-2, ie+2
                     q(i) = q(i) + C1*d(i,j,k)*(u(i+1,j,k)-u(i,j,k))*
     &                          sqrt(gamma(i,j,k)*p(i,j,k)/d(i,j,k))
                  enddo
               else
                  do i = is-2, ie+2
                     q(i) = q(i) + C1*d(i,j,k)*(u(i+1,j,k)-u(i,j,k))*
     &                          sqrt(gamma1*p(i,j,k)/d(i,j,k))
c     &                          0.5*(u(i+1,j,k)+u(i,j,k))
                  enddo
               endif
            endif
c
            q(is-3) = q(is-2)
            q(ie+3) = q(ie+2)

c
c           update velocity1 and energy
c
            do i = is-2, ie+2
               e(i,j,k) = e(i,j,k) + dt1*q(i)/d(i,j,k)*
     &                               (u(i,j,k)-u(i+1,j,k))/dx(i)
            enddo
c
            do i = is-2, ie+3
               u(i,j,k) = u(i,j,k) + dt1*(q(i-1)-q(i))/
     &                               (dx(i)*0.5*(d(i,j,k)+d(i-1,j,k)))
            enddo
c
         enddo
      enddo
c
c           update velocity2 and energy
c
      if (rank .gt. 1) then
         do k = ksm2, kn
            do i = is-2, in
c
               do j = js-2, je+2
                  if ((v(i,j+1,k)-v(i,j,k)) .lt. 0.0) then
                     q(j) = C2*d(i,j,k)*(v(i,j+1,k)-v(i,j,k))**2
                  else
                     q(j) = 0.0
                  endif
               enddo
c
               if (C1 .ne. 0.0) then
                  if (igamfield .eq. 1) then
                     do j = js-2, je+2
                        q(j) = q(j) + C1*d(i,j,k)*(v(i,j+1,k)-v(i,j,k))*
     &                             sqrt(gamma(i,j,k)*p(i,j,k)/d(i,j,k))
                     enddo
                  else
                     do j = js-2, je+2
                        q(j) = q(j) + C1*d(i,j,k)*(v(i,j+1,k)-v(i,j,k))*
     &                             sqrt(gamma1*p(i,j,k)/d(i,j,k))
c     &                          0.5*(v(i,j+1,k)+v(i,j,k))
                     enddo
                  endif
               endif
c
               q(js-3) = q(js-2)
               q(je+3) = q(je+2)
c
c
               do j = js-2, je+2
                  e(i,j,k) = e(i,j,k) + dt1*q(j)/d(i,j,k)*
     &                               (v(i,j,k)-v(i,j+1,k))/dy(j)
               enddo
c
               do j = js-2, je+3
                  v(i,j,k) = v(i,j,k) + dt1*(q(j-1)-q(j))/
     &                               (dy(j)*0.5*(d(i,j,k)+d(i,j-1,k)))
               enddo
            enddo
         enddo
      endif
c
c           update velocity3 and energy
c
      if (rank .gt. 2) then
         do j = jsm2, jn
            do i = is-2, in
c
               do k = ks-2, ke+2
                  if ((w(i,j,k+1)-w(i,j,k)) .lt. 0.0) then
                     q(k) = C2*d(i,j,k)*(w(i,j,k+1)-w(i,j,k))**2
                  else
                     q(k) = 0.0
                  endif
               enddo
c
               if (C1 .ne. 0.0) then
                  if (igamfield .eq. 1) then
                     do k = ks-2, ke+2
                        q(k) = q(k) + C1*d(i,j,k)*(w(i,j,k+1)-w(i,j,k))*
     &                             sqrt(gamma(i,j,k)*p(i,j,k)/d(i,j,k))
                     enddo
                  else
                     do k = ks-2, ke+2
                        q(k) = q(k) + C1*d(i,j,k)*(w(i,j,k+1)-w(i,j,k))*
     &                             sqrt(gamma1*p(i,j,k)/d(i,j,k))
c     &                          0.5*(w(i,j,k+1)+w(i,j,k))
                     enddo
                  endif
               endif
c
               q(ks-3) = q(ks-2)
               q(ke+3) = q(ke+2)
c
c
               do k = ks-2, ke+2
                  e(i,j,k) = e(i,j,k) + dt1*q(k)/d(i,j,k)*
     &                               (w(i,j,k)-w(i,j,k+1))/dz(k)
               enddo
c
               do k = ks-2, ke+3
                  w(i,j,k) = w(i,j,k) + dt1*(q(k-1)-q(k))/
     &                               (dz(k)*0.5*(d(i,j,k)+d(i,j,k-1)))
               enddo
c
            enddo
         enddo
      endif
      enddo
c
c
c
c     3) Substep 3 -- compression term
c
      do k = ksm2, kep1
         do j = jsm2, jep1
c
c           Compute the divergence (should use old u,v,w?)
c
            do i = is-2, ie+1
               div(i) = (u(i+1,j,k) - u(i,j,k))/dx(i)
            enddo
c
            if (rank .gt. 1) then
               do i = is-2, ie+1
                  div(i) = div(i) + (v(i,j+1,k) - v(i,j,k))/dy(j)
               enddo
            endif
c
            if (rank .gt. 2) then
               do i = is-2, ie+1
                  div(i) = div(i) + (w(i,j,k+1) - w(i,j,k))/dz(k)
               enddo
            endif
c
c           Update energy 
c
            do i = is-2, ie+1
               if (igamfield .eq. 1) then
                  alpha = 0.5*dt*(gamma(i,j,k) - 1.0)*div(i)
               else
                  alpha = 0.5*dt*(gamma1 - 1.0)*div(i)
               endif
               e(i,j,k) = e(i,j,k) * (1.0 - alpha)/(1.0 + alpha)
c               if (d(i,j,k) .gt. 1000) write(20,111) i,j,k,d(i,j,k),
c     &                 u(i,j,k),v(i,j,k),w(i,j,k),dt
c 111           format(3i4, 1p, 5(1x,e11.3))
               if (e(i,j,k) .le. 0.0 .or. d(i,j,k) .le. 0.0) then
                  write(6,*) 'zeus_div',e(i,j,k),alpha
                  write(6,*) div(i),i,j,k
                  write(6,*) d(i,j,k),p(i,j,k),dt,dx(i)
                  write(6,*) p(i+1,j,k),p(i,j+1,k),p(i,j,k+1)
                  write(6,*) p(i-1,j,k),p(i,j-1,k),p(i,j,k-1)
                  write(6,*) d(i-1,j,k),d(i,j-1,k),d(i,j,k-1)
                  write(6,*) u(i-1,j,k),u(i,j,k),u(i+1,j,k)
                  write(6,*) v(i,j-1,k),v(i,j,k),v(i,j+1,k)
                  write(6,*) w(i,j,k-1),w(i,j,k),w(i,j,k+1)
                  write(6,*) gr_xacc(i-1,j,k),gr_xacc(i  ,j,k),
     &                                        gr_xacc(i+1,j,k)
                  write(6,*) gr_yacc(i,j-1,k),gr_yacc(i,j  ,k),
     &                                        gr_yacc(i,j+1,k)
                  write(6,*) gr_zacc(i,j,k-1),gr_zacc(i,j,k  ),
     &                                        gr_zacc(i,j,k+1)
                  do n=is,ie
                     write(6,*) d(n,j,k),e(n,j,k),u(n,j,k),
     &                          v(n,j,k),w(n,j,k)
                  enddo
                  CALL f_error("zeus_source.src",625)
               endif
            enddo
c
         enddo
      enddo
c
c     Compute the pressure
c
c
      return
      end
