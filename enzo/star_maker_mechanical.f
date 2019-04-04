























c
c==================================================================================
c////////////////////// SUBROUTINE STAR_MAKER_MECHANICAL \\\\\\\\\\\\\\\\\\\\\\\\\C
      subroutine star_maker_mechanical(nx, ny, nz,
     &            d, dm, temp, u, v, w,  
     &            dt, r, metal, zfield1, zfield2,
     &            dx, t, z, procnum, 
     &            dunits, x1, vunits, t1, 
     &            nmax, xstart, ystart, zstart,
     &            ibuff, imetal,imethod, mintdyn,
     &            odthresh, level, np, xp,yp,zp,
     &            up, vp, wp, mp, tdp, tcp, metalf)
c
c  CREATES STAR PARTICLES FOR MECHANICAL FEEDBACK
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
c
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
c-----------------------------------------------------------------------------------
        implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c-----------------------------------------------------------------------------------
c
c Arguments
        integer*8 nx, ny, nz, ibuff, nmax, np, level
        integer*8 procnum, imetalSNIa, imetal, imethod
        real*8    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
        real*8    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
        real*8    r(nx,ny,nz)
        real*8    metal(nx,ny,nz), zfield1(nx,ny,nz), zfield2(nx,ny,nz)
        real*8    dt, dx, z, odthresh, mintdyn
        real*8    dunits, x1, vunits, t1
        real*8 xstart, ystart, zstart, t
        real*8 xp(nmax), yp(nmax), zp(nmax)
        real*8    up(nmax), vp(nmax), wp(nmax)
        real*8    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
c
c LOCALS
c
        integer*8  i, j, k, ii
        real*8   div, tdyn, grad_rho, tau_cell, phi_cell, hcell
        real*8   psi_cell, f_shielded,dtot, zcell, dcell
        real*8 msolar, pi, G, sndspdC
c
        ii = np
        msolar = 1.9891d33
        sndspdC = 1.3095e8_RKIND
        pi = 3.141592653589793
        G = 6.67428d-8
c
c This is a nearly parameter-free star formation routine.  The tunable parameter
c is odthresh, which sets an overdensity requirement before formation can occur.
c
c for each zone, a star particle is created if the answers to the following are 
c affirmative:
c   
c     is this the finest level of refinement in this grid?
c     is the density > odthresh?
c     is the gas self shielded to some fraction > 0.0?
c     is the flow converging?
c
        do k=1+ibuff, nz-ibuff
          do j = 1+ibuff, ny-ibuff
            do i = 1+ibuff, nx-ibuff
c         1: Is this the finest level?
              if (r(i,j,k) .ne. 0._RKIND) goto 20
c         2: is density > odthresh?
              if (d(i,j,k) .lt. odthresh) goto 20
c         3: is gas self shielded according to Krumholz & Gnedin 2011
              hcell = dx*x1
              dcell = d(i,j,k) * dunits
          zcell = metal(i,j,k)/d(i,j,k) ! metallicity fraction
          grad_rho = abs(d(i+1, j, k) - d(i-1, j, k)
     &                  + d(i, j+1, k) -d(i, j-1, k)
     &                  + d(i, j, k+1) - d(i,j,k-1))
        tau_cell = 434.8_RKIND * dcell * (dx *x1+ dcell/(grad_rho))
        phi_cell = 0.756_RKIND * 
     &              (1+ 3.1_RKIND*zcell/0.02_RKIND)**0.365_RKIND
        psi_cell = (0.6_RKIND * tau_cell 
     &               * (0.01+zcell/0.02_RKIND))
     &               / (log(1_RKIND+0.06_RKIND 
     &               * phi_cell + 0.01_RKIND
     &               * (phi_cell)**2_RKIND))
        f_shielded = 1_RKIND - 3_RKIND/(1_RKIND+ 4_RKIND*psi_cell)
        if (f_shielded .lt. 0._RKIND) goto 20
c         4: Is the flow converging (negative divergence)?
              if (imethod .eq. 2) then
                div = u(i+1,j  ,k  ) - u(i,j,k)
     &                + v(i  ,j+1,k  ) - v(i,j,k)
     &                + w(i  ,j  ,k+1) - w(i,j,k)
              else
                div = u(i+1,j  ,k  ) - u(i-1,j  ,k  )
     &                + v(i  ,j+1,k  ) - v(i  ,j-1,k  )
     &                + w(i  ,j  ,k+1) - w(i  ,j  ,k-1)
              endif
              if (div .ge. 0._RKIND) goto 20
c         If we got to this point, its time to make a star
              ii = ii + 1
c              write(6,*) 'n_created: ', ii-np
c              write(6,*) 'n_max: ',nmax
              mp(ii) = 0.1_RKIND * d(i,j,k)
              tcp(ii) = t
              dtot = ( d(i,j,k) + dm(i,j,k) )*dunits
              tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1
              tdp = max(tdyn, mintdyn*3.15e7_RKIND/t1)
              xp(ii) = xstart+(REAL(i,RKIND)-0.5_RKIND)*dx
              yp(ii) = ystart+(REAL(i,RKIND)-0.5_RKIND)*dx
              zp(ii) = zstart+(REAL(i,RKIND)-0.5_RKIND)*dx
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
c       Remove mass from grid 
              d(i,j,k) = d(i,j,k) - mp(ii)
c              write(6,*) 'star created mass: ', 
c     &        mp(ii)*dunits*(x1*dx)**3*msolar
              if (ii .eq. nmax) goto 30                       
 20     continue              
            enddo !loop i
          enddo ! loop j
        enddo ! loop k
 30     continue
        if (ii .ge. nmax) then
          write(6,*) 'star_maker_mech: reached max new particle count'
          CALL f_error("star_maker_mechanical.src",181)
        endif
        np = ii
      return
      end
c
c
c
c===================================================================================
c////////////////////// SUBROUTINE STAR FEEDBACK MECHANICAL \\\\\\\\\\\\\\\\\\\\\\\c
c
c    Routine to handle mechanical feedback from stars.  Creation is handled as
c    star_maker3mom, i.e., Cen & Ostriker stochastic formation but with stars given
c    the bulk velocity of progenitor gas.  For brevity, use
c
c    StarParticleCreation = 14 (or 16384 in enzo-dev)
c
c    to use the maker from star_maker3mom
c
c    Supernova are handled discretely, drawing the probability of a SN II, Ia.
c    Stars stay active indefinitely, feeding back winds based on the age of the
c    particle:
c    young particles have strong & fast winds from OB stars, old particles
c    have slow winds from AGB.
c
c    This routine currently uses the fact that PPM uses cell-centered qtys
c    to treat neighbor cells as points (ala reverse SPH), so may have very
c    unintended consequences with Zeus hydro.
c
c    Particle creation is handled via method 3 or 14: cen & ostriker +
c    stochastic star formation,

c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\c
      subroutine star_feedback_mechanical(nx, ny, nz,
     &               d, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod,
     &               dt, r, dx, t, z, h, omegaM, omegaL,
     &               dunits, x1, vunits, t1,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type,
     &               star_winds, single_sn)

c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Azton Wells
c  date:      Jan 2019
c
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
c    h     - hubble constant at z=0 km/s/mpc
c    O_m   - omega_matter at z=0
c    O_l   - omega_lambda at z=0
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
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c    star_winds - flag to use continuous winds (stellar mass loss)
c    discrete_sn - flag to use discrete sn or continuous injection
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
      integer nx, ny, nz, ibuff, npart, idual, imetal,
     &      imulti_metals, imethod, single_sn,
     &      distrad, diststep, distcells, star_winds
      real*8    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      real*8    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      real*8    r(nx,ny,nz), ge(nx,ny,nz)
      real*8    metal(nx,ny,nz)
      real*8    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      real*8    dt, dx, z, h, omegaM, omegaL
      real*8    dunits, x1, vunits, t1, justburn
      real*8 xstart, ystart, zstart, t
      real*8 xp(npart), yp(npart), zp(npart)
      real*8    up(npart), vp(npart), wp(npart)
      real*8    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      integer type(npart)
c
c
c  Locals
c
c
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      integer*8 n, ic, jc, kc, ip, jp, kp,
     &            i, j, k, a, ax, p, one_event
      integer*8 n_sn_ii, n_sn_ia, iface, jface, kface
      real*8 mform, zcell, f_factor,
     &         odthresh, div, grad_rho, dcell, acell,
     &         tau_cell, phi_cell, psi_cell, f_shielded,
     &         t_ff, m_form, modxba, num, denom, rmp, rpm,
     &         sum_weights, sum_factor, wind_factor, psi_k, temp,
     &         hcell, metal_ejecta
      real*8 R_ii, R_ia, m_eject, m_winds, wind_metal, snii_metals,
     &          snia_metals, e_sn, p_ej, v_wind,
     &          e_wind, P_ii, P_ia, adotx, m_units,
     &          m_deposited, e_deposited, p_units, e_units
      real*8 mag_xba(3,3,3)
      real*8 f(3,2)
      real*8 xba(3,3,3,3), p_deposited(3), ahats(3,3,3,3)
      real*8 nbor(3,3,3,3)
      real*8 weightsFinal(3,3,3,3)
      real*8 weights(3,3,3)
      real*8 pm(3,3,3,3,2)
      real*8 dxf, dyf, dzf, dxc, dyc, dzc, xfc, yfc, zfc,
     &     xfcshift,yfcshift,zfcshift, xpos, ypos, zpos,
     &     xface, yface, zface, face_shift
      real*8 Zsol, energy_51, eunits, e_const, ergs_51, fbuff, pi
      real*8 msolar, msolar_e51, mH, g, mstar, age

      real*8 tempx, tempy, tempz
c    external random
      real random
      call random_seed
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
      if (imethod .ne. 0) then
c        write (6,*) 'mechanical feedback only works with PPM hydro solver! method= ',imethod
        goto 10
      endif
      one_event = 0._RKIND ! flag for testing: =1 implies only one sn at first step
      msolar = 1.9891d33 ! grams per m_sol
      pi = 3.14159265358979323846_RKIND
      g = 6.67428d-8
      acell = 1._RKIND/(1_RKIND + z)
      m_units = dunits * (x1)**3*dx**3 / msolar     ! Msun
      e_units = dunits * x1**5*dx**3  /t1**2 ! ergs
      p_units = sqrt(m_units*msolar*e_units)/1e5/msolar ! Msun*km/s
c      write(6,*) 'x units: ', x1
c      write(6,*) 'time units: ',t1
c      write(6,*) 'rho units: ', dunits
c      write(6,*) 'mass_units: ', m_units
c      write(6,*) 'energy_units: ', e_units
c      write(6,*) 'momenta_units: ', p_units
c      write(6,*) 'alt p_units: ', dunits * x1**4/t1
      do n=1_IKIND, npart
         if (tcp(n) .gt. 0._RKIND .and. mp(n) .gt. 0._RKIND .and.
     &        type(n) .eq. 2_IKIND) then
c   center of feedback zone
            ip = int((xp(n) - xstart)/dx) + 1_IKIND
            jp = int((yp(n) - ystart)/dx) + 1_IKIND
            kp = int((zp(n) - zstart)/dx) + 1_IKIND

            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*nz) then
               write(6,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart, ip,jp,kp
               goto 100
            endif
c
c	Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + 2._RKIND
            mstar = mp(n) * m_units
            write(6,*) 'mp pre',mstar
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
            xpos = (xfc - xstart)/dx - 0.5_RKIND
            ypos = (yfc - ystart)/dx - 0.5_RKIND
            zpos = (zfc - zstart)/dx - 0.5_RKIND
            ic = int(xpos + 0.5_RKIND)
            jc = int(ypos + 0.5_RKIND)
            kc = int(zpos + 0.5_RKIND)
c#ifdef 
c            write(6,*) 'new index: ',ic,jc,kc
c#endif
c
c    7/?  Transform to comoving coords
c
c
            call mech_momentum(u,v,w,d,metal,up(n), vp(n), wp(n),
     &                        nx, ny, nz, ic, jc, kc,
     &                        iface, jface, kface,
     &                        imethod, imetal, +1_IKIND)

c            temp = 0.0_RKIND
c            do i = 1_IKIND, 4_IKIND
c              do j = 1_IKIND, 4_IKIND
c                do k = 1_IKIND, 4_IKIND
c                  temp = temp + (u(ip-i+1, jp-j+1, kp-k+1)
c     &                            *p_units)**2
c                  temp = temp + (v(ip-i+1, jp-j+1, kp-k+1)
c     &                            *p_units)**2
c                  temp = temp + (w(ip-i+1, jp-j+1, kp-k+1)
c     &                            *p_units)**2
c                enddo
c              enddo
c            enddo
c            write(6,*) 'area momenta pre=', sqrt(temp)
            face_shift = 0._RKIND
            if (imethod .eq. 2) face_shift = 0.5_RKIND
            dxc = real(ic) + 0.5_RKIND - xpos
            dyc = real(jc) + 0.5_RKIND - ypos
            dzc = real(kc) + 0.5_RKIND - zpos
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

c    1/?  Check for continual formation:
c              a.)  require gas still be dense enough for
c                   self shielding using local Sobolev approximation
c                   and that overdensity is > n.
c
c              b.)  make sure divergence of gas flow is still negative
c              c.)  need M_jeans < m_crit
c
            mform = 0._RKIND
            zcell = metal(ip, jp, kp)
            dcell = d(ip,jp,kp) * dunits

c            write(6,*) 'cell_mass: ', d(ic,jc,kc)*m_units
            m_deposited = 0._RKIND
            e_deposited = 0._RKIND
            do i = 1_IKIND, 3_IKIND
              p_deposited(i) = 0._RKIND
            enddo
c
c           a.) divergence of the gas
c
            if (imethod .eq. 2) then
                  div = u(ip+1,jp  ,kp  ) - u(ip,jp,kp)
     &                + v(ip  ,jp+1,kp  ) - v(ip,jp,kp)
     &                + w(ip  ,jp  ,kp+1) - w(ip,jp,kp)
            else
                  div = u(ip+1,jp  ,kp  ) - u(ip-1,jp  ,kp  )
     &                + v(ip  ,jp+1,kp  ) - v(ip  ,jp-1,kp  )
     &                + w(ip  ,jp  ,kp+1) - w(ip  ,jp  ,kp-1)
            endif

c
c           c.) calculate shielded fraction via Krumholz & Gnedin 2011
c
            hcell = dx*x1
            grad_rho = abs(d(ip+1, jp, kp) - d(ip-1, jp, kp)
     &                  + d(ip, jp+1, kp) -d(ip, jp-1, kp)
     &                  + d(ip, jp, kp+1) - d(ip,jp,kp-1))
            tau_cell = 434.8_RKIND * dcell * (dx *x1+ dcell/(grad_rho))
            phi_cell = 0.756_RKIND * (1+ 3.1_RKIND*zcell/0.02_RKIND)**0.365_RKIND
            psi_cell = (0.6_RKIND * tau_cell * (0.01+zcell/0.02_RKIND))
     &               / (log(1_RKIND+0.06_RKIND * phi_cell + 0.01_RKIND
     &               * (phi_cell)**2_RKIND))
            f_shielded = 1_RKIND - 3_RKIND/(1_RKIND+ 4_RKIND*psi_cell)
            if (div .le. 0._RKIND .and. d(ip,jp,kp) .ge. odthresh
     &        .and. f_shielded .ge. 0._RKIND) then
              t_ff = sqrt(3._RKIND*pi / (32._RKIND * g * dcell))/t1
c              write(6,*) 'f, tff, d, dt: ',f_shielded,
c     &              t_ff*t1, d(ic,jc,kc)*m_units, dt
              m_form = dt * f_shielded
     &                  *d(ic,jc,kc)*m_units / (t_ff)
            else
              m_form = 0._RKIND
            endif
c           add mass to star particle
           write(6,*) 'mass formed: ', m_form
           m_form = m_form /m_units
            if (m_form > 0._RKIND) mp(n) = mp(n) + m_form
            mstar = mp(n) * m_units
c           remove mass from the grid cell.
            d(ip, jp, kp) = d(ip, jp, kp) - m_form
c
c
c    2/?  Construct the 3x3 mask to distribute feedback.
c                 make arrays of neighbor positions (x,y,z),
c                 weights.  The weights estimate the amount of
c                 feedback that will be deposited to cell 'b'
c                 and is calculated as in Hopkins, 2017
c
c
c           get a list of positions of neighbor "particles"
c    __ Construct Neighbor particle mask that acts as coupled particles for deposition
c

            do k = 1_IKIND, 3_IKIND
              do j = 1_IKIND, 3_IKIND
                do i = 1_IKIND, 3_IKIND
                  nbor(i,j,k, 1) = xp(n) + dx * (i - 2)
                  nbor(i,j,k, 2) = yp(n) + dx * (j - 2)
                  nbor(i,j,k, 3) = zp(n) + dx * (k - 2)
                  xba(i,j,k, 1) = dx * (i - 2)
                  xba(i,j,k, 2) = dx * (j - 2)
                  xba(i,j,k, 3) = dx * (k - 2)
                  mag_xba(i,j,k) = 0._RKIND
                enddo
              enddo
            enddo
c
c
c    __ Construct Wieghts vector.  Each entry is the fraction of energy that goes into that
c               coupling particle
c
c           __ Area vector for A.dot(xba)
c                 generate vector that has scalar xba^2
            do k = 1_IKIND, 3_IKIND
              do j = 1_IKIND, 3_IKIND
                do i = 1_IKIND, 3_IKIND
                  do ax = 1_IKIND, 3_IKIND
                    mag_xba(i,j,k) = mag_xba(i,j,k)
     &                        + xba(i,j,k,ax)*xba(i,j,k,ax)
                  enddo
                enddo
              enddo
            enddo
c
c         Make ahat vector; entries point from xp to fb particle
c
            do k = 1_IKIND, 3_IKIND
              do j = 1_IKIND, 3_IKIND
                do i = 1_IKIND, 3_IKIND
                  do ax = 1_IKIND, 3_IKIND
                    ahats(i,j,k,ax) = 4.8_RKIND*pi/26._RKIND
     &                             * xba(i,j,k,ax)/sqrt(mag_xba(i,j,k))
                  enddo
                enddo
              enddo
            enddo
c
c
c    __ Make array of scalar weights that determine fraction of
c         feedback that goes into the neighbor particles.
c
            sum_weights= 0.0_RKIND
            do k = 1_IKIND, 3_IKIND
              do j = 1_IKIND, 3_IKIND
                do i = 1_IKIND, 3_IKIND
                  if (i .eq. 2_IKIND
     &                  .and. j .eq. 2_IKIND
     &                  .and. k .eq. 2_IKIND) then
                    weights(i,j,k) = 0._RKIND
                  else
                    temp = 0._RKIND
                    do ax = 1_IKIND, 3_IKIND
                      temp = temp + ahats(i,j,k,ax)
     &                      * xba(i,j,k,ax)/sqrt(mag_xba(i,j,k))
                    enddo
                    weights(i,j,k) = sqrt(1._RKIND
     &                            + temp/(pi*mag_xba(i,j,k)))
                    sum_weights = sum_weights + weights(i,j,k)
                  endif
                enddo
              enddo
            enddo
            do k= 1_IKIND, 3_IKIND
              do j= 1_IKIND, 3_IKIND
                do i = 1_IKIND,3_IKIND
                  weights(i,j,k) = weights(i,j,k) / sum_weights
                enddo
              enddo
            enddo
c            tempx = 0._RKIND
            tempy = 0._RKIND
            tempz = 0._RKIND
c            write(6,*) 'Final Weights:'
            do k= 1_IKIND, 3_IKIND
              do j= 1_IKIND, 3_IKIND
                do i = 1_IKIND,3_IKIND
                  do ax = 1_IKIND, 3_IKIND
                    if (i == 2_IKIND .and. j == 2_IKIND
     &                .and. k == 2_IKIND) then
                      weightsFinal(i,j,k,ax) = 0.0_RKIND
                      continue
                    else
                      weightsFinal(i,j,k,ax) =
     &                      xba(i,j,k, ax)/sqrt(mag_xba(i,j,k))
     &                            *weights(i,j,k)
                    endif
                  enddo
c                  write(6,*) "WF: ",weightsFinal(i,j,k,1), i,j,k
c                  tempx = tempx + weightsFinal(i,j,k,1)
                  tempy = tempy + weightsFinal(i,j,k,2)
                  tempz = tempz + weightsFinal(i,j,k,3)
                enddo
              enddo
            enddo
c            write(6,*) 'Sum weights: ', tempx*m_units,
c     &                      tempy*m_units , tempz*m_units
c
c               Although Hopkins 2017 has tensor corrections using plus/minus
c                 vectors, this formulation does not require it, since
c                 the coupled particles are symmetric about the
c                 star particle.
c
c
c    3/?  if age is low enough, check for supernova
c              rates are per solar mass per Myr
            R_ii = 0._RKIND
            R_ia = 0._RKIND
            age = (t-tcp(n))*t1/3.14e13_RKIND
c            write(6,*) 'age: ',age
            if (age .lt. 0._RKIND) goto 100
            if (single_sn .eq. 1_IKIND) then
              if (age .lt. 37.53_RKIND) then
                R_ia = 0._RKIND
                if (age .lt. 10.37_RKIND
     &               .and. age .gt. 3.401_RKIND) then
                  R_ii = 0.0005408_RKIND
                endif
                if (age .gt. 10.37_RKIND) then
                         R_ii = 0.0002516_RKIND
                endif
              endif
              if (age .ge. 37.53_RKIND) then
                R_ia = 5.3e-8_RKIND * 1.6e-5_RKIND
     &                   * exp(-1.0*(age-50._RKIND)/2._RKIND)
              endif
              if (age .lt. 3.401_RKIND) then
                    R_ia = 0._RKIND
              endif
c
c    4/?  Draw from random numbers to see if SN happen or not this timestep
c
            random = rand()
            random = random
            n_sn_ii = 0._RKIND
            n_sn_ia = 0._RKIND
            if (one_event .gt. 0._RKIND) then
              if (tdp(n) .gt. one_event) goto 110
              if (tdp(n) .lt. one_event 
     &                .and. t .ge. tdp(n)*0.05) then
                n_sn_ii = 1._RKIND
                tdp(n) = tdp(n) + 1._RKIND
                goto 130
              endif
            endif
            P_ii = mstar * R_ii*dt*t1/3.14e13_RKIND
            if (P_ii .gt. 1._RKIND) then
              write(6,*) 'P_ii: ', P_ii
              write(6,*) 'P_ii > 1.0!!'
              write(6,*) 
     &          'Need to reduce timesteps or reduce particle mass!!'
c
c       Allow a SNe or two to go off with no restrictions
c         but after that, require that p_ii < 1 s.t. only 1 sne per timestep!
c
              if (t-tcp(n)*t1/3.14e13_RKIND .gt. 3.7) then 
                CALL f_error("star_maker_mechanical.src",716) 
              endif
            endif
            if (random .le. P_ii) then
               n_sn_ii = 1._RKIND !anint(mp(n) * R_ii * dt * t1/3.14e13_RKIND)
            else
               n_sn_ii = 0._RKIND
            endif
            if (age .ge. 37.53_RKIND) then
              random = rand()
              random = random
              if (random .le. mstar * R_ia * dt * t1/3.14e13_RKIND) then
                n_sn_ia = 1._RKIND !anint(mp(n) * R_ia * dt * t1/3.14e13_RKIND)
              else
                n_sn_ia = 0._RKIND
              endif
            else
              n_sn_ia = 0._RKIND
            endif
 130        if (n_sn_ia .eq. 0._RKIND
     &              .and. n_sn_ii .eq. 0._RKIND) then
c              write(6,*) 'no supernova this step'
              goto 110
            else
              write(6,*) 'SUPERNOVA!!',n, n_sn_ii, n_sn_ia, 
     &          age
            endif
c
c    5/?  Calculate mass ejected, metal ejected and energy from sn
c
            m_eject = 0._RKIND

c
c         i.) mass in Msun of metal from type II
            snii_metals = n_sn_ii * (1.91_RKIND
     &                   + 0.0479*MAX(zcell/0.02_RKIND, 1.65_RKIND))
            m_eject = m_eject + n_sn_ii * 10.5_RKIND
c         ii.) mass in Msun of metals from type Ia
            snia_metals = n_sn_ia * 1.4_RKIND
            m_eject = m_eject + n_sn_ia * 1.4_RKIND
c         iii.) energy in supernova (simple assuming 1e51 ergs/sn)
            e_sn = 1.0e51_RKIND * (n_sn_ia + n_sn_ii) / e_units
c
c    6/? send to add_feedback subroutine with supernova feedback values
c
            m_eject = m_eject/m_units
            metal_ejecta = (snia_metals + snii_metals)
     &                    /m_units
            p_ej = sqrt(2._RKIND * e_sn
     &              * m_eject)
            call mech_add_feedback (nbor, weights,weightsFinal,xba,
     &                        u, v, w,d,ge,te,metal,
     &                        nx,ny,nz,ic,jc,kc,iface,
     &                        jface,kface,dxf,dyf,
     &                        dzf,dxc,dyc,dzc,imethod,
     &                        imetal, idual, m_eject,
     &                        metal_ejecta, p_ej, e_sn,
     &                        m_deposited, e_deposited, p_deposited,
     &                        p_units,zcell, e_units,
     &                        m_units,dunits,dx, 0_IKIND,x1)
            endif
            temp = 1.998e33 * 1e5
c
c    6/?  Calculate mass loss from winds, if applicable
c
c         a.) Mass loss; msun / Gyr
c
c        write(6,*) 'calculating winds'
 110        if (star_winds .eq. 1_IKIND) then
              if (age < 1._RKIND) then
                    wind_factor = 4.763_RKIND
     &                    * min((0.01_RKIND + zcell/0.02_RKIND), 1.0)
     &                           * age
              endif
              if (age > 1._RKIND .and. age < 3.5_RKIND) then
                    wind_factor = 4.763_RKIND *
     &                    min((0.01_RKIND + zcell/0.02_RKIND), 1.0)
     &                            *(age)
     &                           **(1.45_RKIND + 0.08_RKIND
     &                           * min(log(zcell/0.02_RKIND), 1.0))
              endif
              if (age > 3.5_RKIND .and. age < 100_RKIND) then
                    wind_factor =
     &                  29.4 * (age / 3.5_RKIND)**(-3.25_RKIND)
     &                           + 0.0042_RKIND
              endif
              if (age > 100._RKIND) then
                    wind_factor = 0.42 * (age / 1000)**(-1.1)
     &                           / (19.81/log(age))
              endif
c
c              write(6,*) 'wind factor: ',wind_factor
c              write(6,*) 'zcell: ',zcell/0.02
                m_winds = mp(n)
     &                    * wind_factor * dt * t1 / 3.14e16_RKIND
                m_eject = m_winds + m_eject
c         b.) winds metals
                wind_metal = (0.0278_RKIND + 0.0041_RKIND
     &          * min(max(zcell/0.02_RKIND, 1.65_RKIND), 5.0)) * m_winds
c              write(6,*) 'wind metal: ', wind_metal*m_units
c
c         c.) energy in winds
c              i.) l_kinetic
              if (age .lt. 100_RKIND) then
                  psi_k = 5.94e4_RKIND
     &                  / (1._RKIND+ age/2.5_RKIND)**(1.4_RKIND)
     &                  + (t/50._RKIND)**5._RKIND + 4.83_RKIND
              endif
              if (age .ge. 100_RKIND) then
                  psi_k = 4.83_RKIND
              endif
c              ii.) v_wind
c
c              v_wind = sqrt(2_RKIND * psi_k * 10**12_RKIND)
              e_wind = psi_k * 1e12_RKIND * m_winds * m_units
     &                /e_units
              p_ej = sqrt(2._RKIND * e_wind / m_winds)
              write(6,*) 'e_winds: ', e_wind*e_units
c              write(6,*) 'winds p: ', p_ej*p_units
              if (e_wind .lt. 1e-50_RKIND) goto 170
              call mech_add_feedback(nbor, weights,weightsFinal,xba,
     &                        u, v, w,d,ge,te,metal,
     &                        nx,ny,nz,ic,jc,kc,iface,
     &                        jface,kface,dxf,dyf,
     &                        dzf,dxc,dyc,dzc,imethod,
     &                        imetal, idual, m_winds,
     &                        wind_metal, p_ej, e_wind,
     &                        m_deposited, e_deposited, p_deposited,
     &                        p_units,zcell, e_units,
     &                        m_units,dunits,dx, 1_IKIND, x1)
            endif ! star_winds == 1
c    remove mass from particle
 170           mp(n) = mp(n) - m_eject
            temp = 0.0_RKIND
            do i = -8_IKIND, 8_IKIND
              do j = -8_IKIND, 8_IKIND
                do k = -8_IKIND, 8_IKIND
                    temp = temp + (u(ip-i+1, jp-j+1, kp-k+1)
     &                                    *p_units)**2
                    temp = temp + (v(ip-i+1, jp-j+1, kp-k+1)
     &                                    *p_units)**2
                    temp = temp + (w(ip-i+1, jp-j+1, kp-k+1)
     &                                    *p_units)**2
                enddo
              enddo
            enddo
c            write(6,*) 'post metal: ', metal(ic,jc,kc)
c            write(6,*) 'area momenta post=', sqrt(temp)
c
c
            call mech_momentum(u,v,w,d,metal, up(n), vp(n), wp(n),
     &                  nx, ny, nz, ic, jc, kc,
     &                  iface, jface, kface,
     &                  imethod, imetal, -1_IKIND)
c
c       Error checking:
c
c          write(6,*) 'mp post: ', mp(n)*m_units
          temp = 0.0_RKIND
          do i = -8_IKIND, 8_IKIND
            do j = -8_IKIND, 8_IKIND
              do k = -8_IKIND, 8_IKIND
                temp = temp + (u(ip-i+1, jp-j+1, kp-k+1)
     &                    * d(ip-i+1, jp-j+1, kp-k+1)
     &                    * p_units)**2
                temp = temp + (v(ip-i+1, jp-j+1, kp-k+1)
     &                    * d(ip-i+1, jp-j+1, kp-k+1)
     &                                    *p_units)**2
                temp = temp + (w(ip-i+1, jp-j+1, kp-k+1)
     &                    * d(ip-i+1, jp-j+1, kp-k+1)
     &                                    *p_units)**2
              enddo
            enddo
          enddo
          write(6,*) 'm_deposited/m_ejected', m_deposited/m_eject
          write(6,*) 'e_deposited/e_ejected: ',
     &                  e_deposited/(e_sn+e_wind)
          write(6,*) '|P|: ', (p_deposited(1) + p_deposited(2)
     &                + p_deposited(3))*p_units
          endif ! star particle calculation
 100    continue
        enddo    ! end loop over particles
 10   return
      end
c
c
c=================================================================================
c                 Adding Feedback
c=================================================================================
c
cc     Compute and deposit momentum.  This follows the
c              mechanical feedback methods present in Hopkins 2017.
c              Each neighbor cell is taken as a point of gas at cell
c              center to calculate the coupling from stellar feedback.
c              Takes in a list of "particle" positions creating a 3x3 cloud
c              along with weight vectors for each "particle".
c              Each cloud particle is then distributed to the mesh using
c              cloud-in-cell
c
c       nbor - position of neighbor cloud particles that have coupled to feedback
c       weightsFinal - array of weights determining coupling to neighbor particles
c       m_eject - mass of ejecta
c       e_sn    - energy of ejecta
c       p_ej    - total ejected momenta e_sn= p_ej^2/2 m_eject
c       d, metal, te, ge - density, metal, total energy and thermal energy fields
c       nx, ny, nz -number of entries in d, metal, te, ge
c       x1, t1, dunits - conversion of distance, time, density to physical units
c       imetal - metallicity field flag
c       mcell - mass of the feedback particle is ~mass of cell
c       mdep, edep, pdep - running sum of deposited momenta, energy and mass
      subroutine mech_add_feedback(nbor, weights, weightsFinal, xba,
     &                        pu, pv, pw,
     &                        d, ge, te, metal, nx, ny, nz,
     &                        ic, jc, kc, iface, jface, kface,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, idual,
     &                        m_eject, metal_ejecta, p_ej, e_sn,
     &                        m_deposited, e_deposited,p_deposited,
     &                        p_units,zcell,e_units,
     &                        m_units,d_units,dx, winds, xunits)
c
      implicit none


      integer, parameter :: PKIND=8



      integer, parameter :: RKIND=8



      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8

c
c     Arguments
c
      integer*8 nx, ny, nz, ic, jc, kc, idual, a, winds
      integer*8 iface, jface, kface, imethod, imetal, imulti_metals
      real*8    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      real*8    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      real*8    metal_ejecta, m_eject, p_ej, e_sn, delta_metal, delta_m
      real*8    delta_e, p_units, zcell, e_units, t1, nb, dx
      real*8    dxf, dyf, dzf, dxc, dyc, dzc, tempx,xunits,
     &        dxf1, dyf1, dzf1, dxc1, dyc1, dzc1
c
c     Locals
c
      integer*8 i, j, k, i1, j1, k1, n_sn
      real*8    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &          dratio, tot_mass, temp, scalar_weight, p_t,p_rat,
     &          m_deposited, e_deposited, fz,m_units, d_units, m_cell
      real*8  r_cool(3,3,3)
      real*8  delta_p(3), p_deposited(3)
      real*8  weightsFinal(3,3,3,3), nbor(3,3,3,3), xba(3,3,3,3)
      real*8  weights(3,3,3)
      real*8 r_factor(3,3,3), mod_xba(3,3,3)
c
c     Error check
c
      if (imulti_metals .eq. 1_IKIND) then
         write(6,*) "momentum: not supported"
      endif
c
c     Loop "cells" in particle-frame
c
      tot_mass = 0._RKIND
      temp = zcell / 0.02_RKIND
      if (temp .lt. 0.01_RKIND) then
          fz = 2._RKIND
      else
          fz = (zcell/0.02_RKIND)**(-0.14_RKIND)
      endif
      nb = d(ic,jc,kc)*d_units/1.67262171e-24_RKIND
      p_t = (4.8 *1.0/nb**(1._RKIND/7._RKIND)
     &            *  (e_sn*e_units/1e51_RKIND)**(13._RKIND/14._RKIND)
     &            * fz)
     &            *1e5_RKIND / p_units 
      p_rat = p_t/p_ej      
      m_cell = d(ic, jc, kc)
      do k = -1_IKIND, 1_IKIND
        do j = -1_IKIND, 1_IKIND
          do i = -1_IKIND, 1_IKIND
                ! P_t in hopkins is in units Msun*km/s.

            r_cool (i+2, j+2, k+2) = 28.4*3.018e18*nb**(-3.0/7.0)
     &                * (e_sn/1e51*e_units)**(2.0/7.0)*fz
     &                /xunits
            mod_xba(i+2,j+2,k+2) = 
     &              sqrt(xba(i+2,j+2,k+2, 1)*xba(i+2,j+2,k+2,1)
     &                  + xba(i+2,j+2,k+2,2)*xba(i+2,j+2,k+2,2)
     &                  + xba(i+2,j+2,k+2,3)*xba(i+2,j+2,k+2,3))
            r_factor (i+2, j+2, k+2) = 
     &              mod_xba(i+2,j+2, k+2)
     $                    / r_cool(i+2,j+2,k+2)

          enddo
        enddo
      enddo

c
      tempx = 0._RKIND
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
c     calculate mass, energy, momentum, etc for this cell particle
c
            scalar_weight= weights(i+2, j+2, k+2)
            delta_m = m_eject*scalar_weight
            delta_e = e_sn * scalar_weight
c            if (delta_e .lt. 0_RKIND)
c            write(6,*) 'delta_e: ',delta_e*e_units
            delta_metal = metal_ejecta * scalar_weight
              do a = 1_IKIND, 3_IKIND
                temp = sqrt(1._RKIND+m_cell/delta_m)
c                write(6,*) 'temp, p_t/p_ej ',temp, p_t/p_ej
                delta_p(a) = p_ej*min(temp, p_rat)
     &                        *weightsFinal(i+2,j+2,k+2,a)
 120              continue
                enddo

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
                        delta_mass =   delta_m*dxc1*dyc1*dzc1
                        delta_pu   =  delta_p(1)*dxf1*dyc1*dzc1
                        delta_pv   =  delta_p(2)*dxc1*dyf1*dzc1
                        delta_pw   =  delta_p(3)*dxc1*dyc1*dzf1
                        delta_therm = delta_e*dxc1*dyc1*dzc1
c                        write(6,*) 'total energy to deposit:',
c     &                  delta_therm
c
c     Check for energy conservation: if pu^2+pv^2+pw^2/2m > dTherm, rescale?
c                        
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
                        m_deposited = m_deposited+delta_mass
                        p_deposited(1) = delta_pu + p_deposited(1)
                        p_deposited(2) = delta_pv + p_deposited(2)
                        p_deposited(3) = delta_pw + p_deposited(3)
c
c     Add increase in kinetic energy to the total energy field
c
                      temp = (delta_pu**2 + delta_pv**2 + delta_pw**2)
     &                        / (2._RKIND*d(ic+i1, jc+j1, kc+k1))
                      if (delta_therm- temp .lt. 0._RKIND) then
c                        write(6,*) 'KE, dtherm', temp*e_units,
c     &                         delta_therm*e_units, i,j,k
c                        CALL f_error("star_maker_mechanical.src",1086)
                        delta_therm = 0._RKIND
                      else
                        delta_therm = delta_therm - temp
                      endif
          ! Reduce added thermal energy by u > u*r/r_cool if r > r_cool
c                      if (mod_xba(i+2,j+2,k+2) .gt.
c     &                        r_cool(i+2,j+2,k+2)) then
c                       delta_therm = delta_therm 
c     &                    * (mod_xba(i+2,j+2, k+2)
c     &                      / r_cool(i+2,j+2,k+2))**(-6.5_RKIND)
c                      endif
          
          ! Add kinetic and thermal energy to total energy
c            write(6,*) 'e_cell, KE, dtherm',te(ic+i1 ,jc+j1 ,kc+k1),
c     &                       temp,
c     &                       delta_therm
                      te(ic+i1 ,jc+j1 ,kc+k1) =
     &                       te(ic+i1 ,jc+j1 ,kc+k1)*dratio
     &                       + temp/d(ic+i1, jc+j1, kc+k1)
     &                       + delta_therm/d(ic+i1, jc+j1, kc+k1)
c
          ! Add thermal energy to gas energy
                        if (idual .eq. 1_IKIND)
     &                       ge(ic+i1 ,jc+j1 ,kc+k1) =
     &                       ge(ic+i1, jc+j1, kc+k1)*dratio +
     &                       delta_therm / d(ic+i1, jc+j1, kc+k1)
             e_deposited = e_deposited + delta_therm + temp

c     Adding increase in kinetic energy to the total energy field
c     Metal feedback (note that in this function gas metal is
c     a fraction (rho_metal/rho_gas) rather than a density.
c     The conversion has been done in the handling routine)
c
                        if (imetal .eq. 1_IKIND) then
                           metal(ic+i1, jc+j1, kc+k1) =
     &                          metal(ic+i1, jc+j1, kc+k1) +
     &                          delta_metal/d(ic+i1,jc+j1,kc+k1)

                        endif
c
c     End loop over CIC-deposit
c
                     enddo
                  enddo
               enddo
c               write(6,*) 'D_pu = w dp xba(1): ', weightsFinal(i,j,k,1),
c     &               delta_p(1)*p_units
c               write(6,*) 'delta_p(1)=sum(d_pu): ',
c     &              delta_p(1)*p_units, i,j,k
c               write(6,*) 'added: ', tempx*p_units
c               write(6,*) 'total added, x,y,z: ', 
c     &                p_deposited(1)*p_units, p_deposited(2)*p_units,
c     &                p_deposited(3)*p_units
c
c     End loop over "cells" in particle-frame
c
            enddo
         enddo
      enddo
c      write(6,*) 'mass deposited: ', m_deposited*m_units
c
c

        return
        end
c ==============================================================================
c
c ==============================================================================
c
c ==========================================================
c
c     Convert velocities to momentum and back
c
      subroutine mech_momentum(u, v, w, d, metal, up, vp, wp,
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