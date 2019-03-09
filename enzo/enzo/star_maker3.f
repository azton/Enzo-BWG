







c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\c
      subroutine star_maker3(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, zfield1, zfield2,
     &                      dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, masseff, smthresh, level, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf)
c
c  CREATES GALAXY PARTICLES
c
c  written by: Brian O'Shea
c  date:       13 November 2002
c	This file was originally a copy of star_maker2.src,
c	which was originally written by Chris Loken.  As of
c	today, this is intended to be the unigrid version of
c	star_maker2, so the jeans mass and stochastic star
c	formation have been completely removed.  See
c	star_maker2.src for changes previous to 13 Nov. 2002.
c
c  modified1:  20 Dec 2002 by BWO
c       Stochastic star formation is added again.
c       The particle masses are averaged over several cells to avoid
c       the "runaway star particle" phenomenon
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold (some number * avg. density)
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    level - current level of refinement
c    procnum - processor number (for output)
c
c  OUTPUTS:
c
c    np   - number of particles created
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    nmax     - particle array size specified by calling routine
c
c
c-----------------------------------------------------------------------
       implicit none
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      integer procnum
      real    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real    r(nx,ny,nz), cooltime(nx,ny,nz)
      real    metal(nx,ny,nz), zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real    dt, dx, z
      real    d1, x1, v1, t1
      real*8 xstart, ystart, zstart, t
      real*8 xp(nmax), yp(nmax), zp(nmax)
      real    up(nmax), vp(nmax), wp(nmax)
      real    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      real    odthresh, masseff, smthresh, mintdyn
c
      real   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      integer  i, j, k, ii
      real   div, tdyn, dtot
      real   pi, G, sndspdC
      real   isosndsp2, starmass, starfraction, bmass, jeanmass
      double precision msolar
      parameter (pi=3.14159265, G=6.67e-8, sndspdC=1.3095e8,
     &           msolar=1.989e33)
c
      ii = 0

!     print*,'star_maker3: imetal is:',imetal

c
c  for each zone, : "star" particle is created if answers to all the
c  following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density greater than a critical density ?
c    is the flow convergent ?
c    is the cooling time less than a dynamical time ? 
c    is the gas mass greater than the Jeans mass?
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              1) is this finest level of refinement?
c
               if (r(i,j,k) .ne. 0.0) goto 10
c
c              2) is density greater than threshold?
c
               if (d(i,j,k) .lt. odthresh) goto 10
c
c              3) is divergence negative?
c                 (the first calculation is face centered for ZEUS, 
c                  the second is cell-centered for PPM)
c
               if (imethod .eq. 2) then
                  div = u(i+1,j  ,k  ) - u(i,j,k)
     &                + v(i  ,j+1,k  ) - v(i,j,k)
     &                + w(i  ,j  ,k+1) - w(i,j,k)
               else
                  div = u(i+1,j  ,k  ) - u(i-1,j  ,k  )
     &                + v(i  ,j+1,k  ) - v(i  ,j-1,k  )
     &                + w(i  ,j  ,k+1) - w(i  ,j  ,k-1)
               endif
               if (div .ge. 0.0) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3.0*pi/32.0/G/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6.0*sqrt(d(i,j,k)*dble(d1))) *
     &                    dble(pi * isosndsp2 / G)**1.5 / msolar

c
c  THIS IS COMMENTED OUT - NO JEANS MASS CRITERION IN THIS ALGORITHM!!!
c  BWO, 13 NOV 02 (fix 3 dec 02)
c               if (bmass .lt. jeanmass) goto 10
c
c              6) Check to see if star is above threshold (given
c                 in units of M_solar)
c
               starfraction = min(masseff*dt/tdyn, 0.9)
               tdyn = max(tdyn, mintdyn*3.15e7/t1)

c
c  STOCHASTIC STAR FORMATION HAS BEEN ADDED AGAIN - BWO 20 Dec 2002
c
c
c
c                 Keep global count of "unfullfilled" star formation
c                 and when total is larger than threshold, then create
c                 a star particle with the threshold mass or 1/2 the
c                 gas in the cell, whichever is smaller.
c
               if (starfraction*bmass .lt. smthresh) then
                  sformsum = sformsum + starfraction*bmass
                  if (sformsum .lt. smthresh) goto 10
                  starfraction = min(smthresh/bmass, 0.5)
                  sformsum = sformsum - starfraction*bmass
               endif

c
c              Create a star particle
c
               ii = ii + 1
               mp(ii)  = starfraction * d(i,j,k)
               tcp(ii) = t
               tdp(ii) = tdyn
               xp(ii) = xstart + (float(i)-0.5)*dx
               yp(ii) = ystart + (float(j)-0.5)*dx
               zp(ii) = zstart + (float(k)-0.5)*dx
c
c              Star velocities averaged over multiple cells to
c              avoid "runaway star particle" phenomenon
c              imethod = 2 is zeus, otherwise PPM

               if (imethod .eq. 2) then
                  up(ii) = 1.0e-20
                  vp(ii) = 1.0e-20
                  wp(ii) = 1.0e-20

               else

                  up(ii) = 1.0e-20
                  vp(ii) = 1.0e-20
                  wp(ii) = 1.0e-20

               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
!                 write(*,'("Setting metal fraction")')
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
!                 write(*,'("Zero metal fraction")')
                  metalf(ii) = 0.0
               endif
c
c              Remove mass from grid
c
               d(i,j,k) = (1.0 - starfraction)*d(i,j,k)
c
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii)
c
 1000          format(i5,1x,6(1pe10.3,1x))
c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20

10          continue

            enddo
         enddo
      enddo
 20   continue
c	
      if (ii .ge. nmax) then
         write(6,*) 'star_maker2: reached max new particle count'
         CALL f_error("star_maker3.src",262)
      endif
      np = ii
c
c      if (np .ne. 0) then
c         write(6,*) 'Stars created: number,time,level: ', np, t, level
c      endif
cc
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\c
      subroutine star_feedback3(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w,
     &                      metal, zfield1, zfield2,
     &                      idual, imetal, imulti_metals, imethod, 
     &                      dt, r, dx, t, z,
     &                      d1, x1, v1, t1, sn_param, m_eject, yield,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, justburn)
c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Brian O,Shea
c  date:       13 November 2002
c  	this is a copy of star_maker2.src, and is designated
c	as of this writing to be for unigrid runs with star 
c	formation, and for experimentation.  See star_maker2.src
c	for the changelog prior to 13 November 2002.
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metallicity density field
c    r     - refinement field (0 if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imulti_metals - flag to use multi metals zfield 1 and 2
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
       implicit none
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod
      real    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real    r(nx,ny,nz), ge(nx,ny,nz)
      real    metal(nx,ny,nz)
      real    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real    dt, dx, z
      real    d1, x1, v1, t1, justburn
      real*8 xstart, ystart, zstart, t
      real*8 xp(npart), yp(npart), zp(npart)
      real    up(npart), vp(npart), wp(npart)
      real    mp(npart), tdp(npart), tcp(npart), metalf(npart)
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      integer i, j, k, n
      real mform, tfactor, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio
      parameter (clight = 3e10, msolar_e51 = 1800.0)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback2: start'
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0) then
c


c        The star particle creation algorithm partnered with this 
c          feedback algorithm creates a star particle instantaneously.
c          However, we do feedback as if the star particles are created 
c          over a long period of time (in code units), so the particle
c          actually loses mass over time in an exponentially decaying
c          way.

c        Determine how much of a given star particle would have been 
c          turned into stars during this timestep.  Then calculate the mass
c          which should have formed during this timestel dt using the integral
c          form of the Cen & Ostriker formula.

            xv1 = (t      - tcp(n))/tdp(n)
            if (xv1 .gt. 12.0) goto 10     ! t-tcp >> tdp so ignore
            xv2 = (t + dt - tcp(n))/tdp(n)

c        First calculate the initial mass of the star particle 
c          in question.
            minitial = mp(n) / 
     &              (1.0 - m_eject*(1.0 - (1.0 + xv1)*exp(-xv1)))
c
c       Then, calculate the amount of mass that would have formed in
c         this timestep.
c
            mform = minitial * ((1.0 + xv1)*exp(-xv1) - 
     &                          (1.0 + xv2)*exp(-xv2))
            mform = max(min(mform, mp(n)), 0.0)
c
c         Compute index of the cell that the star particle
c           resides in.
c 
            i = int((xp(n) - xstart)/dx) + 1
            j = int((yp(n) - ystart)/dx) + 1
            k = int((zp(n) - zstart)/dx) + 1
c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'warning: star particle out of grid',i,j,k
               goto 100
            endif
c
c          skip if very little mass is formed.
c
            if (mform/d(i,j,k) .lt. 1.0e-10) goto 10
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c
c           Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3
c
c           Calculate how much of the star formation in this
c           timestep would have gone into supernova energy.
c
            energy = sn_param * mform * (clight/v1)**2 / 
     &                (d(i,j,k)+mform*m_eject)
c
c             (aug 7/00 addition: share ejected energy between gas
c              being ejected and that remaining in star particle)
c
c
c           Add energy to energy field
c
            dratio = d(i,j,k)/(d(i,j,k) + mform * m_eject)
            te(i,j,k) = te(i,j,k)*dratio + energy
            if (idual .eq. 1) 
     &          ge(i,j,k) = ge(i,j,k)*dratio + energy

c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
            if (imetal .eq. 1) then
c
c           "Cen method".  This takes into account gas recycling.
c
               metal(i,j,k) = (metal(i,j,k)*d(i,j,k) + mform * 
     &              (yield * (1.0-metalf(n)) + m_eject * metalf(n)))
     &              / (d(i,j,k)+mform*m_eject)  ! metal is a fraction

               if (imulti_metals .eq. 1) then
                  zfield1(i,j,k) = metal(i,j,k)
                  zfield2(i,j,k) = metal(i,j,k)
               endif
c
            endif
c
c           Mass and momentum feedback
c
            u(i,j,k) = u(i,j,k)*d(i,j,k) + mform * m_eject * up(n)
            v(i,j,k) = v(i,j,k)*d(i,j,k) + mform * m_eject * vp(n)
            w(i,j,k) = w(i,j,k)*d(i,j,k) + mform * m_eject * wp(n)
            d(i,j,k) = d(i,j,k) + mform * m_eject
            u(i,j,k) = u(i,j,k)/d(i,j,k)
            v(i,j,k) = v(i,j,k)/d(i,j,k)
            w(i,j,k) = w(i,j,k)/d(i,j,k)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
            if (imethod .ne. 2 .and. idual .eq. 1) then
               te(i,j,k) = 0.5*(u(i,j,k)**2 + v(i,j,k)**2 + 
     &                          w(i,j,k)**2) + ge(i,j,k)
            endif
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback2: end'
      return
      end
