
























c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\c
      subroutine star_maker3mom(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, zfield1, zfield2,
     &                      dx, t, z, procnum, 
     &                      dunits, x1, vunits, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, masseff, smthresh, level, np, 
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      imetalSNIa, metalSNIa, metalfSNIa,
     &                      exptime)

c
c  CREATES STAR PARTICLES FOR KINETIC FEEDBACK
c
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
c    dunits,x1,vunits,t1 - factors to convert d,dx,v,t to physical units
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
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
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
c    metalfSNIa - metallicity fraction of particle (from SN Ia) ! MKRJ
c
c
c-----------------------------------------------------------------------
       implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c-----------------------------------------------------------------------
c
c  Arguments
c
      integer*8 nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      integer*8 procnum, imetalSNIa
      real*8    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      real*8    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real*8    r(nx,ny,nz), cooltime(nx,ny,nz)
      real*8    metal(nx,ny,nz), zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real*8    dt, dx, z, exptime
      real*8    dunits, x1, vunits, t1
      real*8 xstart, ystart, zstart, t
      real*8 xp(nmax), yp(nmax), zp(nmax)
      real*8    up(nmax), vp(nmax), wp(nmax)
      real*8    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      real*8    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      real*8    odthresh, masseff, smthresh, mintdyn
c
      real*8   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      integer*8  i, j, k, ii
      real*8   div, tdyn, dtot
      real*8   pi, G, sndspdC
      real*8   isosndsp2, starmass, starfraction, bmass, jeanmass
      real*8 msolar
      parameter (pi=3.141592653589793d0, G=6.67428d-8, 
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=1.9891d33)
c
      ii = np

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
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) is density greater than threshold?

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
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*dunits
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4_RKIND) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(dunits)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(dunits))) *
     &                    dble(pi * isosndsp2 / G)**1.5_RKIND / msolar

c
c  THIS IS COMMENTED OUT - NO JEANS MASS CRITERION IN THIS ALGORITHM!!!
c  BWO, 13 NOV 02 (fix 3 dec 02)
c               if (bmass .lt. jeanmass) goto 10
c
c              6) Check to see if star is above threshold (given
c                 in units of M_solar)
c
               starfraction = min(masseff*dt/tdyn, 0.9_RKIND)
               tdyn = max(tdyn, mintdyn*3.15e7_RKIND/t1)

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
                  starfraction = min(smthresh/bmass, 0.5_RKIND)
                  sformsum = sformsum - starfraction*bmass
               endif
c
c              Create a star particle
c
               ii = ii + 1
               mp(ii)  = starfraction * d(i,j,k)
               tcp(ii) = t
               tdp(ii) = tdyn
c              If discrete explosions are used, then use tdp as
c              a flag indicating whether the particle has done
c              feedback rather than dynamical time field
               if (exptime .ge. 0._RKIND) then
                  tdp(ii) = 1._RKIND
               endif
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
c
c              Star velocities averaged over multiple cells to
c              avoid "runaway star particle" phenomenon
c              imethod = 2 is zeus, otherwise PPM

               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
!                 write(*,'("Setting metal fraction")')
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
!                 write(*,'("Zero metal fraction")')
                  metalf(ii) = 0._RKIND
               endif
c
c              MKRJ 2/20/08 Do the same for particle metal fraction from SN Ia
c
               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)    ! in here metal is a fraction
               endif
c
c              Remove mass from grid
c
               d(i,j,k) = (1._RKIND - starfraction)*d(i,j,k)
c
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*dunits,z,metalf(ii)
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
         write(6,*) 'star_maker3: reached max new particle count'
         CALL f_error("star_maker3mom.src",267)
      endif
      np = ii
c
c      if (np .ne. 0) then
c         write(6,*) 'Stars created: number,time,level: ', np, t, level
c      endif
c
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\c
      subroutine star_feedback3mom(nx, ny, nz,
     &               d, mu, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, r, dx, t, z,
     &               dunits, x1, vunits, t1, sn_param, m_eject, yield,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type, justburn,
     &               kinf_in,exptime_in)

c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Greg Bryan, Christine Simpson
c  date:      July 2015
c  	this is a copy of star_maker3.F that has been modifed to 
c       include the deposition of kinetic as well as thermal energy
c       and metals into a patch surrounding the star particle via a 
c       CIC approach.
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
c    dunits,x1,vunits,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imulti_metals - flag to use multi metals zfield 1 and 2
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created or
c               if exptime >= 0, a flag for whether the discrete 
c               explosion has occurred
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c    kinf_fixed - fixed kinetic energy fraction for use when var is
c                 false
c    var      - variable kinetic energy flag
c    exptime  - delay time for discrete explosions (if set to -1,
c               continuous energy injection used)
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
       implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c-----------------------------------------------------------------------
c
c  Arguments
c
      integer*8 nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod,
     &      distrad, diststep, distcells
      real*8    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      real*8    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real*8    r(nx,ny,nz), ge(nx,ny,nz)
      real*8    metal(nx,ny,nz)
      real*8    mu(nx,ny,nz)
      real*8    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real*8    dt, dx, z
      real*8    dunits, x1, vunits, t1, justburn
      real*8 xstart, ystart, zstart, t
      real*8 xp(npart), yp(npart), zp(npart)
      real*8    up(npart), vp(npart), wp(npart)
      real*8    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      integer*8 type(npart)
      real*8 	  kinf,kinf_in,exptime_in,exptime
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      integer*8 n, ic, jc, kc, ip, jp, kp,
     &     iface, jface, kface, i, j, k
      real*8 mform, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, 
     &     dist_mass_cells, dist_mom_cells, mass_per_cell, te_per_cell,
     &     thermal_energy_per_cell, kinetic_energy_per_cell,
     &	   energy_per_cell, mom_per_cell, face_shift,
     &     xface, yface, zface, xpos, ypos, zpos,
     &     energy_before, energy_after, mass_before, mass_after,
     &     kin_energy_before, kin_energy_after,
     &     asum, bsum, csum, delta_ke, ke_injected, ke_after
      real*8 dxf, dyf, dzf, dxc, dyc, dzc, xfc, yfc, zfc,
     &     xfcshift,yfcshift,zfcshift
      real*8 u1(4,4,4), v1(4,4,4), w1(4,4,4), d1(4,4,4),
     &     te1(4,4,4), ge1(4,4,4), metal1(4,4,4),
     &     ke_before(4,4,4)
      real*8 Zsol, d_ave, num_d, energy_51, eunits, R_PDS, t_PDS, 
     &	   R_resolve, KE_rem, e_const, ergs_51, realmass, mu_cell, fbuff
      double precision msolar, mH
      parameter (clight = 2.99792458d10, msolar_e51 = 1800._RKIND, 
c     &           msolar = 1.9891d33,
     &	         mH = 1.67262171d-24, e_const = 2.71828_RKIND,
     &	         ergs_51 = 52.57_RKIND)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
      write(6,*) 'star_feedback3: start',npart
      do n=1_IKIND, npart
         if (tcp(n) .gt. 0._RKIND .and. mp(n) .gt. 0._RKIND .and. 
     &        type(n) .eq. 2_IKIND) then
c
c        The star particle creation algorithm partnered with this 
c        feedback algorithm creates a star particle instantaneously.
c        This feedback routine can either be used to also inject energy
c        instantaneously or to do feedback over a longer period of time
c        as done in star_maker3 modeled on Cen & Ostriker's method that
c        accounts for the true (unsimulated) longer formation time of
c        a stellar population.  The instananeous injection mode is 
c        intended for use with low mass star particles that produce a
c        handful of SN.
c
c     Do instantaneous injeciton
c
            if (exptime_in .ge. 0._RKIND) then
c              convert explosion time to code units
               exptime = exptime_in * 1.e6_RKIND * 3.15e7_RKIND / t1  
               write(6,*) 'mform, minitial = ',mform,minitial
               write(6,*) 'exptime,t,tcp,tdp = ',exptime,t,tcp(n),tdp(n)
               if (tdp(n) .eq. 0._RKIND
     &              .or. abs(t-tcp(n)) .lt. exptime) goto 10
               tdp(n) = 0._RKIND
               minitial = mp(n)
               mform = mp(n)
               goto 5
            endif

c     Determine how much of a given star particle would have been 
c     turned into stars during this timestep.  Then calculate the mass
c     which should have formed during this timestel dt using the integral
c     form of the Cen & Ostriker formula.

            xv1 = (t      - tcp(n))/tdp(n)
            write(6,*) t,tcp(n),tdp(n),xv1
            if (xv1 .lt. 0._RKIND) goto 10
            if (xv1 .gt. 12._RKIND) goto 10 ! t-tcp >> tdp so ignore
            xv2 = (t + dt - tcp(n))/tdp(n)

c     First calculate the initial mass of the star particle 
c     in question.
            minitial = mp(n) / 
     &           (1._RKIND - 
     &            m_eject*(1._RKIND - (1._RKIND + xv1)*exp(-xv1)))
c     
c     Then, calculate the amount of mass that would have formed in
c     this timestep.
c     
            mform = minitial * ((1._RKIND + xv1)*exp(-xv1) - 
     &           (1._RKIND + xv2)*exp(-xv2))
            mform = max(min(mform, mp(n)), 0._RKIND)
            
c     Compute index of the cell that the star particle
c     resides in.
c 
 5          ip = int((xp(n) - xstart)/dx) + 1_IKIND
            jp = int((yp(n) - ystart)/dx) + 1_IKIND
            kp = int((zp(n) - zstart)/dx) + 1_IKIND
c
c     Assuming 3x3x3 cube, calculate cell distribution
c
            distrad = 3_IKIND
            dist_mass_cells = distrad**3_IKIND
c
c          skip if very little mass is formed.
c
            if (mform/dist_mass_cells/d(ip,jp,kp) .lt. 1e-10_RKIND) 
     &           goto 10
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c     
c     Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3_IKIND
c     
c           Calculate mass per cell ejected
c     
            mass_per_cell = mform * m_eject / dist_mass_cells
c
c           Calculate how much of the star formation in this
c           timestep would have gone into supernova energy.
c            
            energy = sn_param * mform * (clight/vunits)**2_IKIND
            energy_per_cell = energy / dist_mass_cells

            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*nz) then
               write(6,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
               goto 100
            endif
c
c	Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + 2._RKIND
c
c         check bounds - if star particle is near grid edge
c         then shift center of feedback region
c
            if (xfc .lt. xstart+fbuff*dx .or. 
     &          xfc .gt. xstart+dx*nx-fbuff*dx .or.
     &          yfc .lt. ystart+fbuff*dx .or. 
     &          yfc .gt. ystart+dx*ny-fbuff*dx .or.
     &          zfc .lt. zstart+fbuff*dx .or. 
     &          zfc .gt. zstart+dx*nz-fbuff*dx) then
               write(6,*) 'warning1: star feedback zone shifted',
     &              xfc,yfc,zfc, xstart, ystart, zstart,fbuff
c               
	       xfcshift = xfc
	       yfcshift = yfc 
	       zfcshift	= zfc
c
	       xfc = max(xfc,xstart+fbuff*dx) 
	       yfc = max(yfc,ystart+fbuff*dx) 
	       zfc = max(zfc,zstart+fbuff*dx) 
c              
	       xfc = min(xfc,xstart+dx*nx-fbuff*dx) 
	       yfc = min(yfc,ystart+dx*ny-fbuff*dx) 
	       zfc = min(zfc,zstart+dx*nz-fbuff*dx)
c
	       xfcshift = xfcshift - xfc
	       yfcshift = yfcshift - yfc 
	       zfcshift	= zfcshift - zfc
c
	       write(6,*) 'warning2: star feedback zone shifted',
     &              xfc,yfc,zfc,dx,xfcshift,yfcshift,zfcshift, 
     &		    nx, ny,nz
            endif
c
c           If using zeus, then velocities are face-centered so shift
c
            face_shift = 0._RKIND
            if (imethod .eq. 2) face_shift = 0.5_RKIND
c
c         Compute index of the first cell to add momentum,
c             accounting for possible face-centering
c
            xface = (xfc - xstart)/dx - 0.5_RKIND - face_shift
            yface = (yfc - ystart)/dx - 0.5_RKIND - face_shift
            zface = (zfc - zstart)/dx - 0.5_RKIND - face_shift
c
            iface = int(xface + 0.5_RKIND)
            jface = int(yface + 0.5_RKIND)
            kface = int(zface + 0.5_RKIND)
c
            dxf = real(iface) + 0.5_RKIND - xface
            dyf = real(jface) + 0.5_RKIND - yface
            dzf = real(kface) + 0.5_RKIND - zface
c
c         Compute index of the first cell to add mass, assuming cell-centering
c 
            xpos = (xfc - xstart)/dx - 0.5_RKIND
            ypos = (yfc - ystart)/dx - 0.5_RKIND
            zpos = (zfc - zstart)/dx - 0.5_RKIND
c
            ic = int(xpos + 0.5_RKIND)
            jc = int(ypos + 0.5_RKIND)
            kc = int(zpos + 0.5_RKIND)
c
            dxc = real(ic) + 0.5_RKIND - xpos
            dyc = real(jc) + 0.5_RKIND - ypos
            dzc = real(kc) + 0.5_RKIND - zpos
c
c     Use fixed kinf value, unless kinf < 0 - then compute
c     variable kinf
            kinf = kinf_in
            if (kinf .lt. 0._RKIND) then
c     Calculate variable kinf based on R_PDS from
c     Cioffi et al. 1988
c
               realmass = mform * dunits * dx * x1 *dx *x1 *dx *x1
               realmass = realmass /msolar
               energy_51 = sn_param * msolar_e51 * realmass
c
               Zsol = 0._RKIND
               num_d = 0._RKIND
               d_ave = 0._RKIND
c     
c     Compute average density and metallcity around cell containing
c     particle
c     
c               write(6,*) 'averaging over patch'
               do k = -1_IKIND,1_IKIND
                  do j = -1_IKIND,1_IKIND
                     do i = -1_IKIND,1_IKIND
                        mu_cell = mu(ic+i,jc+j,kc+k)
                        Zsol = Zsol + 
     &                       metal(ic+i,jc+j,kc+k)/0.02_RKIND
                        num_d = num_d + 
     &                       d(ic+i,jc+j,kc+k)*dunits/mu_cell/mH
                        d_ave = d_ave + d(ic+i,jc+j,kc+k)*dunits
                     enddo
                  enddo
               enddo
c               
               Zsol = Zsol/27._RKIND
               num_d = num_d/27._RKIND
               d_ave = d_ave/27._RKIND
c   
c     compute time and radius of transition to PDS phase for gas
c     with the computed properties
c     Note: t_PDS has units of 1e3yrs and R_PDS has units of 3.0857d18
c          
c     For metal poor gas
c     
               if (Zsol .lt. 0.01_RKIND) then
                  t_PDS = 3.06e2_RKIND * 
     &                 (energy_51**(1._RKIND/8._RKIND)) *
     &                 (num_d**(-3._RKIND/4._RKIND))
                  R_PDS = 49.3_RKIND * (energy_51**0.25_RKIND) * 
     &                 (num_d**(-0.5_RKIND))
               endif
c
c     For metal rich gas
c
               if (Zsol .ge. 0.01_RKIND) then
                  t_PDS = 26.5_RKIND * 
     &                 (energy_51**(3._RKIND/14._RKIND)) * 
     &                 (Zsol**(-5._RKIND/14._RKIND)) * 
     &                 (num_d**(-4._RKIND/7._RKIND))
                  R_PDS = 18.5_RKIND * energy_51**(2._RKIND/7._RKIND) 
     &                    * num_d**(-3._RKIND/7._RKIND)
     &                    * Zsol**(-1._RKIND/7._RKIND)
               endif
               R_resolve = dx*x1/3.0857d18
               if (R_PDS .gt. 4.5*R_resolve) then
                  kinf = 0._RKIND
                  write(6,*) 'Sedov phase resolved:',Zsol,num_d,
     &                 energy_51,R_PDS,R_resolve, realmass
               else                  
                  kinf = 3.97133e-6_RKIND * (d_ave/mH) 
     &                 * (R_resolve**(-2._RKIND)) 
     &                 * (R_PDS**7._RKIND) * (t_PDS**(-2._RKIND)) 
     &                 * (energy_51**(-1._RKIND))

                  write(6,*) 'Sedov phase not resolved:',kinf,Zsol,
     &                 num_d, d_ave, 
     &                 energy_51,t_PDS,R_PDS,R_resolve
                  
               endif
            endif
c
c     Check kinf.  If it is too small, set it to 0
c
            if (kinf .le. 1e-10_RKIND) then 
               write(6,*) 'kinf = ',kinf
               write(6,*) 'kinf too small!  Setting to 0'
               kinf = 0._RKIND
            endif

c
c     Finished computing kinf; now compute the amount of thermal energy
c     that needs to be added to each cell
            thermal_energy_per_cell = (1._RKIND - kinf) * 
     &                                energy_per_cell

c            write(6,*) 'thermal_energy_per_cell:',
c     &           thermal_energy_per_cell,kinf,energy_per_cell
c
c     Inject energy, mass & metals
c
c     Zero local dummy field and kinetic energy field
c
            do k = 1_IKIND, 4_IKIND
               do j = 1_IKIND, 4_IKIND
                  do i = 1_IKIND, 4_IKIND
                     u1(i,j,k) = 0._RKIND
                     v1(i,j,k) = 0._RKIND
                     w1(i,j,k) = 0._RKIND
                     d1(i,j,k) = 0._RKIND
                     ge1(i,j,k) = 0._RKIND
                     te1(i,j,k) = 0._RKIND
                     metal1(i,j,k) = 0._RKIND
                  enddo
               enddo
            enddo
c
c     Compute the kinetic energy in the affected region
c     before momentum is added (except for ZEUS).  This
c     is needed at the end of the calculation to update the
c     total energy (te) field.
c
            if (imethod .ne. 2_IKIND) then
               do k = -1_IKIND, +2_IKIND
                  do j = -1_IKIND, +2_IKIND
                     do i = -1_IKIND, +2_IKIND
                        ke_before(i+2, j+2, k+2) = 
     &                       0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                                (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                 v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                 w(ic+i ,jc+j ,kc+k)**2_IKIND)                     
                     enddo
                  enddo
               enddo
            endif
c
c     First convert velocities to momenta and transform
c     into frame comoving with particle
c
            call momentum(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, 
     &                    +1_IKIND)
c
c           Sum mass and energy before
c
            call sum_mass_kinetic_energy(u, v, w, d, ge, te, nx, ny, nz,
     &           iface, jface, kface, ic, jc, kc,
     &           mass_before, kin_energy_before, idual, imethod)
c            write(6,*) "m,e before:", mass_before, kinetic_energy_before
c
c           Now add mass and momentum terms (normalization 1.0) to
c              local dummy fields
c
            call add_feedback(u1, v1, w1, d1, ge1, te1, metal1, 4_IKIND, 
     &                        4_IKIND, 4_IKIND, 
     &                        2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND, 
     &                        2_IKIND, 2_IKIND,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n),
     &                        mass_per_cell, 1._RKIND, 0._RKIND)
c
c           The kinetic energy after the feedback event is a quadratic equation
c           with delta_p, where delta_p is the total amount of momentum added
c           to the feedback region:
c           E_k,a = asum + bsum*delta_p + csum*delta_p*delta_p
c
c           Sum a, b and c terms to get momentum normalization
c
            if (kinf .gt. 0._RKIND) then
               call sum_abc(u, v, w, d, ge, u1, v1, w1, d1, 
     &              nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &              asum, bsum, csum)
               write(6,*) 'abc1',asum, bsum, csum, kinf
               asum = asum - (kin_energy_before + kinf*energy)
c               write(6,*) 'kin_energy_before,kinf*energy,energy:',
c     &              kin_energy_before,kinf*energy,energy
               write(6,*) 'abc2',asum, bsum, csum
c
c           Calculate momentum contribution
c
               mom_per_cell = (-bsum + 
     &              sqrt(bsum**2_IKIND - 4._RKIND*asum*csum))/
     &              (2._RKIND*csum)
            else
c
c           If the kinetic fraction is set to 0, then 0 momentum is added
c           to each cell.  Note that this is different from adding 0 kinetic
c           energy to 
c
               mom_per_cell = 0._RKIND
            endif
c           
            write(6,*) 'mass_per_cell:',mass_per_cell
            write(6,*) 'mom_per_cell:',mom_per_cell
            write(6,*) 'thermal_energy_per_cell:',
     &           thermal_energy_per_cell
c
c           Now add mass and momentum, using three-point CIC
c           Note that if the te field is being used, it is 
c           only updated with the thermal energy
c
      write(6,*) 'momentum per ptcl: ', mom_per_cell*dunits*x1**3._RKIND*x1/t1

            call add_feedback(u, v, w, d, ge, te, metal, nx, ny, nz, 
     &                        ic, jc, kc, iface, jface, kface, 
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf(n),
     &                        mass_per_cell, mom_per_cell, 
     &			      thermal_energy_per_cell)
c
c           Sum mass and energy before
c
            call sum_mass_kinetic_energy(u, v, w, d, ge, te, nx, ny, nz,
     &           iface, jface, kface, ic, jc, kc,
     &           mass_after, kin_energy_after, idual, imethod)
c
            write(6,*) 'energy:',energy
            write(6,*) 'mass (b,a,preda):', mass_before, mass_after, 
     &           mass_before + mform*m_eject,mform*m_eject
            write(6,*) 'kinetic energy (b,a,preda):',kin_energy_before,
     &           kin_energy_after, kin_energy_before + kinf*energy, 
     &           kinf
c
c           Error checks to see if the kinetic energy added in the particle frame
c           matches the values set by the input parameters
c
            if (abs(kin_energy_after - kin_energy_before - kinf*energy)
     &           /(kinf*energy) > 0.01_RKIND .and. 
     &           kinf .ne. 0._RKIND) then
c
               write(6,*) 'Kinetic energy added to mesh does not match', 
     &              ' kinetic energy feedback parameters'
               write(6,*) 'kinf:',kinf
               write(6,*) 'kin_energy_after - kin_energy_before',
     &              kin_energy_after - kin_energy_before
               write(6,*) 'energy*fkin',kinf*energy
               write(6,*) '%diff',abs(kin_energy_after - 
     &              kin_energy_before)
     &              /(kin_energy_before)
               CALL f_error("star_maker3mom.src",840)
            endif

            if (abs(mass_after - mass_before 
     &              - mass_per_cell*dist_mass_cells)
     &           /(mass_per_cell*dist_mass_cells) > 0.01_RKIND) then
c
               write(6,*) 'Mass added to mesh does not match', 
     &              ' kinetic energy feedback parameters'
               write(6,*) 'kinf:',kinf
               write(6,*) 'mass_after, mass_before:',
     &                     mass_after, mass_before
               write(6,*) 'kin_energy_after - kin_energy_before',
     &                     kin_energy_after - kin_energy_before
               write(6,*) 'energy*fkin',kinf*energy
               write(6,*) '%diff',abs(mass_after - 
     &              mass_before)
     &              /(mass_before)
               CALL f_error("star_maker3mom.src",858)
            endif



c
c           Convert momenta back to velocities and transform back to lab frame
c
            call momentum(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal,
     &                    -1_IKIND)
c
c           Add the increase in the kinetic energy to the total energy field
c             (unless we're using Zeus).  If using dual energy formalism, we
c             might want to enforce consistency.
c     
            if (imethod .ne. 2_IKIND) then
               ke_injected = 0._RKIND
               do k = -1_IKIND, +2_IKIND
                  do j = -1_IKIND, +2_IKIND
                     do i = -1_IKIND, +2_IKIND
                        ke_after = 0.5_RKIND*
     &                        d(ic+i, jc+j, kc+k) *
     &                       (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                        v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                        w(ic+i ,jc+j ,kc+k)**2_IKIND) 
                        delta_ke = ke_after - ke_before(i+2,j+2,k+2)
                        te(ic+i ,jc+j ,kc+k) = 
     &                       te(ic+i ,jc+j ,kc+k) + 
     &                       delta_ke/d(ic+i, jc+j, kc+k)
                        ke_injected = ke_injected + delta_ke
                     enddo
                  enddo
               enddo
            endif
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback3: end'
      return
      end
c
c ==========================================================
c
c     Convert velocities to momentum and back
c
      subroutine momentum(u, v, w, d, metal, up, vp, wp, 
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, idir)
c
      implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c     Arguments
c
      integer*8 nx, ny, nz, ic, jc, kc
      integer*8 iface, jface, kface, imethod, imetal, idir
      real*8    d(nx, ny, nz), metal(nx,ny,nz)
      real*8    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real*8    up, vp, wp
c
c     Locals
c
      integer*8 i, j, k
c
c     Error check
c
      if (idir .ne. -1_IKIND .and. idir .ne. 1_IKIND) then
         write(6,*) 'incorrect idir value in momentum call'
         CALL f_error("star_maker3mom.src",933)
      endif
c
c     Loop over velocities, multiplying by densities (or dividing if
c       converting back)
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
c
c              idir = +1: convert vel -> mom
c               
               if (idir. eq. +1_IKIND) then
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i ,jc+j ,kc+k) = 
     &                    (u(iface+i ,jc+j ,kc+k) - up)
     &                    * 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                    d(iface+i+1_IKIND, jc+j, kc+k))
                     v(ic+i ,jface+j ,kc+k) = 
     &                    (v(ic+i ,jface+j ,kc+k) - vp)
     &                    * 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                    d(ic+i, jface+j+1_IKIND, kc+k))
                     w(ic+i ,jc+j ,kface+k) = 
     &                    (w(ic+i ,jc+j ,kface+k) - wp)
     &                    * 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
     &                    d(ic+i, jc+j, kface+k+1_IKIND))
                  else
                     u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                                      d(ic+i, jc+j, kc+k)
                     v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                                      d(ic+i, jc+j, kc+k)
                     w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                                      d(ic+i, jc+j, kc+k)
                  endif
                  if (imetal .eq. 1_IKIND) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) *
     &                    d(ic+i,jc+j,kc+k)
                  endif
c
c              if idir = -1: convert mom -> vel
c               
               else
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i ,jc+j ,kc+k) = u(iface+i ,jc+j ,kc+k)
     &                    /( 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                             d(iface+i+1_IKIND, jc+j, kc+k))) + up
                     v(ic+i ,jface+j ,kc+k) = v(ic+i ,jface+j ,kc+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                             d(ic+i, jface+j+1_IKIND, kc+k))) + vp
                     w(ic+i ,jc+j ,kface+k) = w(ic+i ,jc+j ,kface+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
     &                             d(ic+i, jc+j, kface+k+1_IKIND))) + wp
                  else
                     u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + up
                     v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + vp
                     w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + wp
                  endif
                  if (imetal .eq. 1_IKIND) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) /
     &                    d(ic+i,jc+j,kc+k)
                  endif
               endif
c
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine sum_mass_kinetic_energy(pu, pv, pw, d, ge, te, 
     &                           nx, ny, nz,
     &                           iface, jface, kface, ic, jc, kc,
     &                           mass_sum, kin_energy_sum,
     &                           idual, imethod)
c
      implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c     Arguments
c
      integer*8 nx, ny, nz, iface, jface, kface, ic, jc, kc, idual
      integer*8 imethod
      real*8    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      real*8    d(nx,ny,nz), ge(nx,ny,nz), te(nx,ny,nz)
      real*8    mass_sum, kin_energy_sum
c
c     Locals
c
      integer*8 i, j, k
      real*8    mass_term, mom_term, kin_energy
c
      mass_sum = 0._RKIND
      kin_energy_sum = 0._RKIND
c
c     Sum mass and energy      
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
               mass_term = d(ic   +i ,jc   +j ,kc   +k)
               mom_term = pu(iface+i ,jc   +j, kc   +k)**2_IKIND + 
     &                    pv(ic   +i ,jface+j, kc   +k)**2_IKIND + 
     &                    pw(ic   +i ,jc   +j, kface+k)**2_IKIND
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2._RKIND * mass_term)
               mass_sum = mass_sum + mass_term
               kin_energy_sum = kin_energy_sum + kin_energy
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine sum_abc(pu, pv, pw, d, ge, pu1, pv1, pw1, d1, 
     &                   nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &                   asum, bsum, csum)
c
      implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c     Arguments
c
      integer*8 nx, ny, nz, iface, jface, kface, ic, jc, kc
      real*8    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      real*8    d(nx,ny,nz), ge(nx,ny,nz)
      real*8    pu1(4,4,4), pv1(4,4,4), pw1(4,4,4), d1(4,4,4)
      real*8    asum, bsum, csum
c
c     Locals
c
      integer*8 i, j, k, istart, jstart, kstart
      real*8    mass_term, mom_term, aterm, bterm, cterm
      real*8    mass_sum, energy_sum
c
      asum = 0._RKIND
      bsum = 0._RKIND
      csum = 0._RKIND
      mass_sum = 0._RKIND
      energy_sum = 0._RKIND
c
      istart = 2_IKIND
      jstart = 2_IKIND
      kstart = 2_IKIND
c
c     Loop over all affected cells
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
               mass_term = d(ic   +i ,jc   +j ,kc+k)
               mom_term = pu(iface+i ,jc   +j, kc+k   )**2_IKIND + 
     &                    pv(ic   +i ,jface+j, kc+k   )**2_IKIND +
     &                    pw(ic   +i ,jc   +j, kface+k)**2_IKIND

               mass_term = mass_term + d1(istart+i,jstart+j,kstart+k)
               asum = asum + mom_term/(2._RKIND * mass_term)
               bterm = 
     &           pu(iface+i,jc+j,kc+k)*pu1(istart+i,jstart+j,kstart+k) +
     &           pv(ic+i,jface+j,kc+k)*pv1(istart+i,jstart+j,kstart+k) +
     &           pw(ic+i,jc+j,kface+k)*pw1(istart+i,jstart+j,kstart+k)
               bsum = bsum + bterm/mass_term ! factors of 2 cancel
               cterm = 
     &           pu1(istart+i,jstart+j,kstart+k)**2_IKIND +
     &           pv1(istart+i,jstart+j,kstart+k)**2_IKIND +
     &           pw1(istart+i,jstart+j,kstart+k)**2_IKIND
c               
               csum = csum + cterm/(2._RKIND*mass_term)
            enddo
         enddo
      enddo
c
c      write(6,*) "m,e in sum:", mass_sum, energy_sum
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine add_feedback(pu, pv, pw, d, ge, te, metal, nx, ny, nz, 
     &                        ic, jc, kc, iface, jface, kface,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf,
     &                        mass_per_cell, mom_per_cell, 
     &                        therm_per_cell)
c
      implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c     Arguments
c
      integer*8 nx, ny, nz, ic, jc, kc, idual
      integer*8 iface, jface, kface, imethod, imetal, imulti_metals
      real*8    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      real*8    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      real*8    mass_per_cell, mom_per_cell, therm_per_cell
      real*8    m_eject, yield, metalf
      real*8    dxf, dyf, dzf, dxc, dyc, dzc,
     &        dxf1, dyf1, dzf1, dxc1, dyc1, dzc1
c
c     Locals
c
      integer*8 i, j, k, i1, j1, k1
      real*8    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &     dratio, tot_mass
c
c     Error check
c
      if (imulti_metals .eq. 1_IKIND) then
         write(6,*) "momentum: not supported"
         CALL f_error("star_maker3mom.src",1163)
      endif
c
c     Loop "cells" in particle-frame
c
      tot_mass = 0._RKIND
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
c     For each particle "cell", do CIC-like deposit
c     compute zone and face centered weight factors
c
               do i1 = i, i+1_IKIND
                  dxf1 = dxf
                  dxc1 = dxc
                  if (i1 .eq. i+1_IKIND) dxf1 = 1._RKIND - dxf
                  if (i1 .eq. i+1_IKIND) dxc1 = 1._RKIND - dxc
                  do j1 = j, j+1_IKIND
                     dyf1 = dyf
                     dyc1 = dyc
                     if (j1 .eq. j+1_IKIND) dyf1 = 1._RKIND - dyf
                     if (j1 .eq. j+1_IKIND) dyc1 = 1._RKIND - dyc
                     do k1 = k, k+1_IKIND
                        dzf1 = dzf
                        dzc1 = dzc
                        if (k1 .eq. k+1_IKIND) dzf1 = 1._RKIND - dzf
                        if (k1 .eq. k+1_IKIND) dzc1 = 1._RKIND - dzc

                        delta_mass =   mass_per_cell*dxc1*dyc1*dzc1
                        delta_pu   =  i*mom_per_cell*dxf1*dyc1*dzc1
                        delta_pv   =  j*mom_per_cell*dxc1*dyf1*dzc1
                        delta_pw   =  k*mom_per_cell*dxc1*dyc1*dzf1
                        delta_therm = therm_per_cell*dxc1*dyc1*dzc1
c
c     Add mass, momentum
c     (add thermal energy here)
c
                        dratio = d(ic+i1, jc+j1, kc+k1)/
     &                       (d(ic+i1, jc+j1, kc+k1) + delta_mass)
                        d(ic+i1 ,jc+j1 ,kc+k1) = d(ic+i1 ,jc+j1 ,kc+k1)
     &                       + delta_mass
                        pu(iface+i1 ,jc+j1 ,kc+k1) = 
     &                       pu(iface+i1 ,jc+j1 ,kc+k1)+ delta_pu
                        pv(ic+i1 ,jface+j1 ,kc+k1) = 
     &                       pv(ic+i1 ,jface+j1 ,kc+k1) + delta_pv
                        pw(ic+i1 ,jc+j1 ,kface+k1) = 
     &                       pw(ic+i1 ,jc+j1 ,kface+k1)+ delta_pw
c                        
                        tot_mass = tot_mass + delta_mass
c
c     Add thermal energy
c     
                        te(ic+i1 ,jc+j1 ,kc+k1) = 
     &                       te(ic+i1 ,jc+j1 ,kc+k1)*dratio + 
     &                       delta_therm / d(ic+i1, jc+j1, kc+k1)
c                        
                        if (idual .eq. 1_IKIND)
     &                       ge(ic+i1 ,jc+j1 ,kc+k1) = 
     &                       ge(ic+i1, jc+j1, kc+k1)*dratio +
     &                       delta_therm / d(ic+i1, jc+j1, kc+k1)
                        
c     
c     Metal feedback (note that in this function gas metal is
c     a fraction (rho_metal/rho_gas) rather than a density.
c     The conversion has been done in the handling routine)
c     
                        if (imetal .eq. 1_IKIND) then
c     
c     "Cen method".  This takes into account gas recycling.
c     (metal is not a fraction here)
c     
                           metal(ic+i1, jc+j1, kc+k1) = 
     &                          metal(ic+i1, jc+j1, kc+k1) + 
     &                          (delta_mass/m_eject)
     &                          * (yield * (1._RKIND-metalf) + 
     &                          m_eject * metalf)
                           
                        endif
c     
c     End loop over CIC-deposit
c     
                     enddo
                  enddo
               enddo
c     
c     End loop over "cells" in particle-frame
c     
            enddo
         enddo
      enddo
c
c     write(6,*) 'tot_mass:',tot_mass
c     
      return
      end