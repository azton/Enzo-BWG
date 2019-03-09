







!=======================================================================
!//////////////////////  SUBROUTINE MULTI_COOL  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine multi_cool(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, imcool, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht,
     &                dt, aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma, 
     &                icmbTfloor, iClHeat, iClMMW,
     &                clMetNorm, clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating, clMMW)


!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
!  modified4: November, 2003 by Robert Harkness; tighten convergence
!
!  PURPOSE:
!    Solve the energy cooling equations.
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE

!  Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, imethod,
     &        idual, iexpand, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, imcool,
     &        iradtrans
      real    dt, aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      real    photogamma(in,jn,kn)
      real    inutot(nfreq), avgsighp, avgsighep, avgsighe2p

!  Cloudy cooling data

      integer icmbTfloor, iClHeat, iClMMW, clGridRank, clDataSize
      integer clGridDim(clGridRank)
      real clMetNorm, clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)),
     &     clPar5(clGridDim(5))
      real clCooling(clDataSize), clHeating(clDataSize),
     &     clMMW(clDataSize)

!  STATIC global data

      integer nratec, ih2co, ipiht

      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)
      real    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp

!  Parameters

      integer itmax, ijk
      parameter (itmax = 10000, ijk = 4103)


      real tolerance
      parameter (tolerance = 1.0e-10)

      double precision mh
      parameter (mh = 1.67d-24)

!  Local scalars

      integer i, j, k, n, iter
      real dom, energy
      real dt2, ttmin, comp1, comp2

!  Row locals

      logical itmask(ijk) 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), tgas(ijk), tgasold(ijk)
      double precision edot(ijk)

!  Cooling/heating row locals

      double precision ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk)
      real hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set units

      dom      = urho*(aye**3)/mh

!     Convert densities from comoving to proper

!$omp parallel
!$omp-  shared(is, ie, js, je, ks, ke, in, jn, kn)
!$omp-  shared(d, e, ge, u, v, w, de)
!$omp-  shared(HI, HII, HeI, HeII, HeIII)
!$omp-  shared(nratec, iexpand, imethod,idual, ispecies, imetal, imcool)
!$omp-  shared(ih2co, ipiht, dt, aye, temstart, temend, idim)
!$omp-  shared(utem, uxyz, uaye, urho, utim, eta1, eta2, gamma)
!$omp-  shared(ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa)
!$omp-  shared(ciHeISa, ciHeIIa, reHIIa, reHeII1a)
!$omp-  shared(reHeII2a, reHeIIIa, brema, compa)
!$omp-  shared(comp_xraya, comp_temp, piHI, piHeI, piHeII)
!$omp-  shared(HM, H2I, H2II, DI, DII, HDI, metal)
!$omp-  shared(hyd01ka, h2k01a, vibha, rotha, rotla)
!$omp-  shared(gpldla, gphdla, hdltea, hdlowa)
!$omp-  shared(inutot, iradtype, nfreq, imetalregen)
!$omp-  shared(iradshield, avgsighp, avgsighep, avgsighe2p)
!$omp-  shared(gaHIa, gaH2a, gaHea, gaHpa, gaela)
!$omp-  shared(icmbTfloor, iClHeat, iClMMW)
!$omp-  shared(clMetNorm, clEleFra, clGridRank, clGridDim)
!$omp-  shared(clPar1, clPar2, clPar3, clPar4, clPar5)
!$omp-  shared(clDataSize, clCooling, clHeating, clMMW)
!$omp-  shared(iradtrans, photogamma)
!$omp-  private(i, j, k, n, iter)
!$omp-  private(dom, energy, dt2, ttmin, comp1, comp2)
!$omp-  private(indixe, t1, t2, logtem, tdef, p2d, dtit)
!$omp-  private(ttot, tgas, tgasold, edot, itmask)
!$omp-  private(ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII)
!$omp-  private(reHII, reHeII1, reHeII2, reHeIII, brem)
!$omp-  private(hyd01k, h2k01, vibh, roth, rotl)
!$omp-  private(gpldl, gphdl, hdlte, hdlow)
!$omp-  default(none)

!$omp do
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)/aye**3
               de(i,j,k)    = de(i,j,k)/aye**3
               HI(i,j,k)    = HI(i,j,k)/aye**3
               HII(i,j,k)   = HII(i,j,k)/aye**3
               HeI(i,j,k)   = HeI(i,j,k)/aye**3
               HeII(i,j,k)  = HeII(i,j,k)/aye**3
               HeIII(i,j,k) = HeIII(i,j,k)/aye**3
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)/aye**3
                  H2I(i,j,k)  = H2I(i,j,k)/aye**3
                  H2II(i,j,k) = H2II(i,j,k)/aye**3
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)/aye**3
                  DII(i,j,k) = DII(i,j,k)/aye**3
                  HDI(i,j,k) = HDI(i,j,k)/aye**3
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)/aye**3
               enddo
            endif
         enddo
      enddo
!$omp end do
 
!     Solve energy cooling (subcycle)

!     Loop over rows of cells

!$omp do
      do k = ks+1, ke+1
       do j = js+1, je+1

!       tolerance = 1.0e-06 * dt

        do i = is+1, ie+1
           itmask(i) = .true.
        end do

!       Set time elapsed to zero for each cell

        do i = is+1, ie+1
           ttot(i) = 0.0
        enddo

!       Loop over cooling subcycles
     
        do iter = 1, itmax

!       Compute the cooling rate on this row

          call cool1d_multi(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idim,
     &                is, ie, j, k, ih2co, ipiht, iter,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, 
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma, 
     &                icmbTfloor, iClHeat, iClMMW,
     &                clMetNorm, clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating, clMMW,
     &                itmask
     &                     )

!         Compute maximum timstep that keeps any fractional change to 10%
!         (maximum timestep is 1/2 hydro step).
!         Then, use each cell~s individual timestep tp update it~s energy
!         Find minimum elapsed timestep (ttot)

          dt2 = dt/2.0
          ttmin = 1.0e+20

          do i = is+1, ie+1

!            Set energy of this cell (the gamma used here is the right
!            one even for H2 since p2d is calculated with this gamma).

             if (tgas(i) .le. temstart .and. edot(i) .lt. 0.0) 
     &              edot(i) = 1.0e-20*1.0e-3

             if (abs(edot(i)) .lt. 1.0e-20) edot(i) = 1.0e-20

             energy = max(p2d(i)/(gamma-1.0), 1.0e-20)

c            energy = max(ge(i,j,k)*d(i,j,k), p2d(i)/(gamma-1.0), 
c    &                    1.0e-20)
c            if (energy .lt. 1.0e-20) energy = d(i,j,k)*(e(i,j,k) - 
c    &              0.5*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
c            energy = p(i,j,k)/(gamma-1.0)


!            Compute timestep for 10% change

c            if (iter .gt. 100) then
c               dtit(i) = min(real(abs(0.01*energy/edot(i))), 
c    &                        dt-ttot(i), dt2)
c            else
c               dtit(i) = min(real(abs(0.1*energy/edot(i))), 
c    &                        dt-ttot(i), dt2)
c            endif

             dtit(i) = min(real(abs(0.1*energy/edot(i))),
     &                     dt-ttot(i), dt2)

             if ( dt-ttot(i) <= tolerance*dt ) then
                itmask(i) = .false.
             end if

             if ( itmask(i) ) then


 1000        format(4(i4,1x),1p,10(e14.3))


!            Update total and gas energy

             e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dtit(i)


             if (idual .eq. 1) then
                ge(i,j,k) = ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i)

c               ge(i,j,k) = max(ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i),
c     &                      0.5*ge(i,j,k))
c            if (ge(i,j,k) .le. 1.0e-20) ge(i,j,k) = (energy + 
c     &           edot(i)*dtit(i))/d(i,j,k)


             endif

!            Update time-stepping

             ttot(i) = ttot(i) + dtit(i)
             ttmin = min(ttot(i), ttmin)

             end if  ! test of itmask(i)

          enddo  ! end of loop over cells(i)


!         If the all cells are done then skip out of loop

          if (abs(dt-ttmin) .lt. tolerance*dt) go to 8888

         enddo  ! end of iteration loop

 8888    continue
      
!       Abort if iteration count exceeds maximum

         if (iter .gt. itmax) then
            write(6,*) 'MULTI_COOL iter > ',itmax,' at j,k =',j,k
            write(0,*) 'FATAL error (2) in OMP_MULTI_COOL'
            write(0,'(" dt = ",1pe10.3," ttmin = ",1pe10.3)') dt, ttmin
            write(0,'((16(1pe8.1)))') (dtit(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (ttot(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (edot(i),i=is+1,ie+1)
            CALL f_error("OMP_NEW_multi_cool.src",447)
         endif

         if (iter .gt. itmax/2) then
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
         end if

!        call open_mpi_error_file( 'F30', 30, 'unknown' )
!        write(30,'("J=",i4,4x,"K=",i4)') j,k
!        write(30,'((12(1pd10.3)))') (ge(i,j,k),i=is+1,ie+1)
!        call close_mpi_error_file( 30 )


!      Next j,k row

       enddo

      enddo
!$omp end do

!     Convert densities back to comoving from proper

!$omp do
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*aye**3
               de(i,j,k)    = de(i,j,k)*aye**3
               HI(i,j,k)    = HI(i,j,k)*aye**3
               HII(i,j,k)   = HII(i,j,k)*aye**3
               HeI(i,j,k)   = HeI(i,j,k)*aye**3
               HeII(i,j,k)  = HeII(i,j,k)*aye**3
               HeIII(i,j,k) = HeIII(i,j,k)*aye**3
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*aye**3
                  H2I(i,j,k)  = H2I(i,j,k)*aye**3
                  H2II(i,j,k) = H2II(i,j,k)*aye**3
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*aye**3
                  DII(i,j,k) = DII(i,j,k)*aye**3
                  HDI(i,j,k) = HDI(i,j,k)*aye**3
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*aye**3
               enddo
            endif
         enddo
      enddo
!$omp end do
!$omp end parallel

      return
      end
