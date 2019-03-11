



















c=======================================================================
c///////////////////////  SUBROUTINE EXPAND_TERMS  \\\\\\\\\\\\\\\\\\\\c
      subroutine expand_terms(rank, isize, idual, coef, imethod, gamma,
     &                        p, pdual, d, e, ge, u, v, w,
     &                        dold, eold, geold, uold, vold, wold)
c
c  ADDS THE COMOVING EXPANSION TERMS TO THE PHYSICAL VARIABLES
c
c     written by: Greg Bryan
c     date:       February, 1996
c     modified1:
c
c  PURPOSE:
c         Note: p is modified on exit
c
c  INPUTS:
c    isize   - size of fields (treat as 1d)
c    idual   - dual energy flag (1 = on, 0 = off)
c    coef    - coefficent (dt * adot / a)
c    d       - density field
c    p       - pressure field (from total energy - 0.5v^2)
c    pdual   - pressure field (from gas energy, if available)
c    e,ge    - total energy and gas energy (specific)
c    u,v,w   - velocities
c
c  OUTPUTS:
c    d,e,ge,u,v,w - output fields
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer isize, idual, imethod, rank
      real    gamma, coef
      real    d(isize), p(isize), e(isize), ge(isize),
     &        u(isize), v(isize), w(isize), pdual(isize)
      real    dold(isize), eold(isize), geold(isize),
     &        uold(isize), vold(isize), wold(isize)
c
c     Locals
c
      integer i
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c   METHOD1 is sortof time-centered
c   METHOD2 is time-backward
c   METHOD3 is semi-implicit
c
c   (If this is changed, you must also change the pressure computation
c    in Grid_ComovingExpandsionTerms.C)
c
c
c     Do gas energy first (if necessary); the term is 3*p/d
c
      if (idual .eq. 1) then
         do i = 1, isize
            ge(i) = ge(i)*(1.0-coef)/(1.0+coef)
c            if (ge(i)-coef*(3.0-2.0/(gamma-1.0))*pdual(i)/d(i).le.0)
c     &         write(6,*) 'get:',i,ge(i),p(i),d(i)
c
c  this line should be there if gamma != 5/3:
c            ge(i) = ge(i) - coef*(3.0 - 2.0/(gamma-1.0))*pdual(i)/d(i)
c
         enddo
      endif
c
c     Now do total energy; the term is 3*p/d + v^2
c       (for zeus method (imethod=2), only use 3*p/d term)
c
      if (imethod .eq. 2) then
         do i = 1, isize
            e(i) = max(e(i) - coef*6.0*p(i)/(d(i)+dold(i)), 0.5*e(i))
         enddo
      else
         do i = 1, isize
            e(i) = e(i)*(1.0-coef)/(1.0+coef)
            e(i) = max(e(i) - coef*(3.0 - 2.0/(gamma-1.0))*p(i)/d(i),
     &                 0.5*e(i))
         enddo
      endif
c
c     Velocity terms
c
c
c        i) sortof time-centered: */
c
c
c        iii) time-forward */
c
c
c        iii) semi-implicit way: */
c
      do i = 1, isize
                          u(i) = u(i)*(1.0 - 0.5*coef)/(1.0 + 0.5*coef)
         if (rank .gt. 1) v(i) = v(i)*(1.0 - 0.5*coef)/(1.0 + 0.5*coef)
         if (rank .gt. 2) w(i) = w(i)*(1.0 - 0.5*coef)/(1.0 + 0.5*coef)
      enddo
c
c
      return
      end
