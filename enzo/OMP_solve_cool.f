




















c=======================================================================
c//////////////////////  SUBROUTINE SOLVE_COOL  \\\\\\\\\\\\\\\\\\\\\\\c
      subroutine solve_cool(
     &                d, e, ge, u, v, w,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, idim,
     &                is, js, ks, ie, je, ke, 
     &                dt, aye, temstart, temend, fh,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, coola)
c
c  SOLVE RADIATIVE COOLING/HEATING EQUATIONS (NON-MULTI SPECIES VERSION)
c
c  written by: Greg Bryan
c  date:       March, 1997
c  modified1:  Robert Harkness
c  date:       November 2003
c              Tighten up convergence criteria
c
c  PURPOSE:
c    Solve the equilbrium energy cooling.
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
      real    dt, aye, temstart, temend, fh,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma, coola(nratec)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn)
c
c  Parameters
c
      integer itmax, ijk
      parameter (itmax = 5000, ijk = 4103)


      real tolerance
      parameter (tolerance = 1.0e-10)

c  Locals
c
      integer i, j, k, iter
      real energy, dt2, ttmin
c
c  Slice locals
c 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), tgas(ijk), tgasold(ijk),
     &     cool(ijk), edot(ijk)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Loop over slices (in the k-direction)
c

!$omp parallel
!$omp-  shared(is, ie, js, je, ks, ke, in, jn, kn)
!$omp-  shared(d, e, ge, u, v, w)
!$omp-  shared(idual, idim, nratec, iexpand, imethod)
!$omp-  shared(aye, dt, temstart, temend, fh)
!$omp-  shared(utem, uxyz, uaye, urho, utim, eta1, eta2, gamma, coola)
!$omp-  private(t1, t2, logtem, tdef, p2d, dtit, ttot, tgas, tgasold)
!$omp-  private(cool, edot, indixe)
!$omp-  private(dt2, ttmin, energy)
!$omp-  private(i, j, k, iter)
!$omp-  default(none)

!$omp do

      do k = ks+1, ke+1
       do j = js+1, je+1
c
c        Set time elapsed to zero for each cell (and convert dens to abs)
c
         do i = is+1, ie+1
            ttot(i) = 0.0
            d(i,j,k) = d(i,j,k)/aye**3
         enddo
c
c        Loop over cooling subcycles
c
         do iter = 1, itmax
c
c        Compute the cooling rate on a slice
c
            call cool1d(d, e, ge, u, v, w,
     &                  in, jn, kn, nratec, idual, idim, imethod, iter,
     &                  is, ie, j, k, 
     &                  temstart, temend, fh, utem,
     &                  eta1, eta2, gamma, coola,
     &                  indixe, t1, t2, logtem, tdef, edot,
     &                  tgas, tgasold, p2d, cool
     &                     )
c
c         Compute maximim timstep that keeps any fractional change to 10%
c              (maximum timestep is 1/2 hydro step).
c         Then, use each cell~s individual timestep tp update it~s energy
c         Find minimum elapsed timestep (ttot)
c
          dt2 = dt/2.0
          ttmin = 1.0e+20
c
          do i = is+1, ie+1
             if (tgas(i) .le. temstart .and. edot(i) .lt. 0.0) 
     &              edot(i) = 1.0e-20
             if (abs(edot(i)) .lt. 1.0e-20) edot(i) = 1.0e-20
             energy = max(p2d(i)/d(i,j,k)/(gamma-1.0), 1.0e-20)
             dtit(i) = min(abs(0.1*energy/edot(i)), 
     &                      dt-ttot(i), dt2)
c


c
c            Update total and gas energy
c
             e(i,j,k)  = e(i,j,k) + edot(i)*dtit(i)

             if (idual .eq. 1) then
                ge(i,j,k) = ge(i,j,k) + edot(i)*dtit(i)
c               ge(i,j,k) = max(ge(i,j,k) + edot(i)*dtit(i),
c     &                      0.5*ge(i,j,k))
c            if (ge(i,j,k) .le. 1.0e-20) ge(i,j,k) = (energy + 
c     &           edot(i)*dtit(i))

                if (ge(i,j,k) .le. 0.0)
     &          write(3,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
             endif
c
c            Update time-stepping
c
             ttot(i) = ttot(i) + dtit(i)
             ttmin = min(ttot(i), ttmin)
          enddo
c
c         If the all cells are done then skip out of loop
c
          if (abs(dt-ttmin) .lt. tolerance*dt) go to 8888
c
c         Next iteration
c
         enddo
c
8888     continue
c      
c       Print warning if iteration count exceeds maximum
c
        if (iter .gt. itmax) then
           write(6,*) 'SOLVE_COOL iter exceeds ',itmax,' at j,k =',j,k
           CALL f_error("OMP_solve_cool.src",191)
        endif
        if (iter .gt. itmax/2) write(6,*) 'COOL iter,j,k =',iter,j,k
c
c       Set the density back to comoving
c
        do i = is+1, ie+1
           d(i,j,k) = d(i,j,k)*aye**3
        enddo
c
c     Next j,k row
c
       enddo
      enddo

!$omp end do
!$omp end parallel

c
      return
      end
