




















!=======================================================================
!//////////////////////  SUBROUTINE COOL1D_MULTI  \\\\\\\\\\\\\\\\\\\\\
      subroutine cool1d_multi(
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

!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
!  modified4: January, 2011 by DRR, added photo-heating
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

      integer in, jn, kn, is, ie, j, k, imethod, idim,
     &        idual, iexpand, ispecies, imetal, 
     &        nfreq, iradshield, iradtype, imetalregen, iter, imcool,
     &        iradtrans

      real    aye,
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

!  Cloudy cooling data arguments

      integer icmbTfloor, iClHeat, iClMMW
      integer clGridRank, clDataSize
      integer clGridDim(clGridRank)

      real clMetNorm, clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)),
     &     clPar5(clGridDim(5))
      real clCooling(clDataSize), clHeating(clDataSize),
     &     clMMW(clDataSize)

!     WORK arrays defined in multi_cool

      double precision ceHI(in), ceHeI(in), ceHeII(in),
     &     ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &     reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &     brem(in)
      real hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &     gpldl(in), gphdl(in), hdlte(in), hdlow(in)
      real gaHI(in), gaH2(in), gaHe(in), gaHp(in), gael(in),
     &     galdl(in)

      logical itmask(in)
      integer indixe(in)
      real    t1(in), t2(in), logtem(in), tdef(in),
     &        tgas(in), tgasold(in), p2d(in)
      double precision edot(in)

      real    comp1, comp2
      real    inutot(nfreq), avgsighp, avgsighep, avgsighe2p

!     STATIC global data

      integer nratec, ih2co, ipiht
      real    temstart, temend

      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec),
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec),
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec)

      real    compa, comp_xraya, comp_temp, piHI, piHeI, piHeII

      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec),
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)

!  Parameters

      double precision mh
      real    ZSOLAR
      parameter (mh = 1.67d-24)      !DPC
      parameter (ZSOLAR = 0.02041)

!  Locals

      integer i, j1, iradfield
      real dom, qq, vibl, logtem0, logtem9, dlogtem, energy, zr
      real dt2, ttmin, scoef, acoef, HIdot,
     &     hdlte1, hdlow1, gamma2, x, nH2, nother, fudge, fH2,
     &     gphdl1, factor
      double precision coolunit, dbase1, tbase1, xbase1, rtunits



!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1)

!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1.0/(aye*uaye) - 1.0
      fudge    = 1.0
      iradfield = -1

!     Set compton cooling coefficients (and temperature)

      if (iexpand .eq. 1) then
         comp1 = compa * (1.0 + zr)**4
         comp2 = 2.73  * (1.0 + zr)
      else
         comp1 = 1.0e-20
         comp2 = 1.0e-20
      endif

!     Compute Pressure

      if (imethod .eq. 2) then

!        Zeus - e() is really gas energy

         do i = is+1, ie+1
            if ( itmask(i) ) then
            p2d(i) = (gamma - 1.0)*d(i,j,k)*e(i,j,k)
            end if
         enddo
      else
         if (idual .eq. 1) then

!           PPM with dual energy -- use gas energy

            do i = is+1, ie+1
               if ( itmask(i) ) then
               p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i,j,k)
               end if
            enddo
         else

!           PPM without dual energy -- use total energy

            do i = is+1, ie+1
               if ( itmask(i) ) then
               p2d(i) = e(i,j,k) - 0.5*u(i,j,k)**2
               if (idim .gt. 1) p2d(i) = p2d(i) - 0.5*v(i,j,k)**2
               if (idim .gt. 2) p2d(i) = p2d(i) - 0.5*w(i,j,k)**2
               p2d(i) = max((gamma - 1.0)*d(i,j,k)*p2d(i), 1.0e-20)
               end if
            enddo
         endif
      endif

!     Compute temperature

      do i = is+1, ie+1
         if ( itmask(i) ) then
         tgas(i) = 
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &            HI(i,j,k) + HII(i,j,k) + de(i,j,k)
         end if
      enddo

!          (include molecular hydrogen, but ignore deuterium)

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            tgas(i) = tgas(i) +
     &            HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2.0
            end if
         enddo
      endif

      do i = is+1, ie+1
         if ( itmask(i) ) then
         tgas(i) = max(p2d(i)*utem/tgas(i), temstart)
         end if
      enddo

!     Correct temperature for gamma from H2

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            nH2 = 0.5*(H2I(i,j,k) + H2II(i,j,k))
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            if (nH2/nother .gt. 1.0e-3) then
               x = 6100/tgas(i) ! not quite self-consistent
               if (x .gt. 10.0) then
                  gamma2 = 0.5*5.0
               else
                  gamma2 = 0.5*(5.0 + 2.0*x**2 * exp(x)/(exp(x)-1)**2)
               endif
            else
               gamma2 = 2.5
            endif
            gamma2 = 1.0 + (nH2 + nother)/
     &                     (nH2*gamma2 + nother/(gamma-1.0))
            tgas(i) = tgas(i) * (gamma2 - 1.0)/(gamma - 1.0)
            end if
         enddo
      endif

!     If this is the first time through, just set tgasold to tgas

      if (iter .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            tgasold(i) = tgas(i)
            end if
         enddo
      endif

!     --- 6 species cooling ---

      do i = is+1, ie+1
         if ( itmask(i) ) then

!        Compute log temperature and truncate if above/below table max/min

         logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Compute index into the table and precompute parts of linear interp

         indixe(i) = min(nratec-1,
     &                   max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)

!        Lookup cooling values and do a linear temperature in log(T)

         ceHI(i) = ceHIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))/tdef(i)
         ceHeI(i) = ceHeIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))/tdef(i)
         ceHeII(i) = ceHeIIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))/tdef(i)
         ciHI(i) = ciHIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))/tdef(i)
         ciHeI(i) = ciHeIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))/tdef(i)
         ciHeIS(i) = ciHeISa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))/tdef(i)
         ciHeII(i) = ciHeIIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))/tdef(i)
         reHII(i) = reHIIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))/tdef(i)
         reHeII1(i)=reHeII1a(indixe(i)) + (logtem(i) - t1(i))
     &        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))/tdef(i)
         reHeII2(i)=reHeII2a(indixe(i)) + (logtem(i) - t1(i))
     &        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))/tdef(i)
         reHeIII(i)=reHeIIIa(indixe(i)) + (logtem(i) - t1(i))
     &        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))/tdef(i)
         brem(i) = brema(indixe(i)) + (logtem(i) - t1(i))
     &         *(brema(indixe(i)+1) -brema(indixe(i)))/tdef(i)

         end if
      enddo

!     Compute the cooling function

      do i = is+1, ie+1
         if ( itmask(i) ) then
         edot(i) = (

!                    Collisional excitations

     &             - ceHI  (i)*HI  (i,j,k)*de(i,j,k)             ! ce of HI
     &             - ceHeI (i)*HeII(i,j,k)*de(i,j,k)**2*dom/4.0  ! ce of HeI
     &             - ceHeII(i)*HeII(i,j,k)*de(i,j,k)/4.0         ! ce of HeII

!                    Collisional ionizations

     &             - ciHI  (i)*HI  (i,j,k)*de(i,j,k)             ! ci of HI
     &             - ciHeI (i)*HeI (i,j,k)*de(i,j,k)/4.0         ! ci of HeI
     &             - ciHeII(i)*HeII(i,j,k)*de(i,j,k)/4.0         ! ci of HeII
     &             - ciHeIS(i)*HeII(i,j,k)*de(i,j,k)**2*dom/4.0  ! ci of HeIS

!                    Recombinations

     &             - reHII  (i)*HII  (i,j,k)*de(i,j,k)          ! re of HII
     &             - reHeII1(i)*HeII (i,j,k)*de(i,j,k)/4.0      ! re of HeII
     &             - reHeII2(i)*HeII (i,j,k)*de(i,j,k)/4.0      ! re of HeII
     &             - reHeIII(i)*HeIII(i,j,k)*de(i,j,k)/4.0      ! re of HeIII

!                    Compton cooling or heating

     &             - comp1*(tgas(i)-comp2)*de(i,j,k)/dom

!                    X-ray compton heating

     &             - comp_xraya * (tgas(i)-comp_temp)*de(i,j,k)/dom

!                    Bremsstrahlung

     &             - brem(i)*(HII(i,j,k)+HeII(i,j,k)/4.0+HeIII(i,j,k))
     &                        *de(i,j,k)
     &                 )

         end if
      enddo
     
!     --- H2 cooling ---

      if (ispecies .gt. 1) then

         do i = is+1, ie+1
            if ( itmask(i) ) then
            gaHI(i) = gaHIa(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaHIa(indixe(i)+1) - gaHIa(indixe(i)))/tdef(i)
            gaH2(i) = gaH2a(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaH2a(indixe(i)+1) - gaH2a(indixe(i)))/tdef(i)
            gaHe(i) = gaHea(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaHea(indixe(i)+1) - gaHea(indixe(i)))/tdef(i)
            gaHp(i) = gaHpa(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaHpa(indixe(i)+1) - gaHpa(indixe(i)))/tdef(i)
            gael(i) = gaela(indixe(i)) + (logtem(i) - t1(i))
     &         *(gaela(indixe(i)+1) - gaHIa(indixe(i)))/tdef(i)
            gphdl(i) = gphdla(indixe(i)) + (logtem(i) - t1(i))
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))/tdef(i)
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then

            galdl(i) = gaHI(i) * HI(i,j,k)  + gaH2(i) * H2I(i,j,k)
     &               + gaHe(i) * HeI(i,j,k) + gaHp(i) * HII(i,j,k)
     &               + gael(i) * de(i,j,k)
            gphdl1 = gphdl(i)/dom
            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1.0 + gphdl1/galdl(i)) / (2.0*dom)

            end if
         enddo

      endif

!     --- Cooling from HD ---

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
c CMB cooling floor
               if (tgas(i) .gt. comp2) then
                  hdlte(i) = hdltea(indixe(i)) + (logtem(i) - t1(i))
     &             *(hdltea(indixe(i)+1) - hdltea(indixe(i)))/tdef(i)
                  hdlow(i) = hdlowa(indixe(i)) + (logtem(i) - t1(i))
     &             *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))/tdef(i)
               else
                  hdlte(i) = 1.0e-20
                  hdlow(i) = 1.0e-20
               end if
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
c  old (incorrect) way:
c               hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2.0)
c               hdlow1 = max(hdlow(i), 1.0e-20)
c               edot(i) = edot(i) - HDI(i,j,k)*
c     .                     (hdlte1/(1.0 + hdlte1/hdlow1)/(2.0*dom))
c  new (correct) way: (april 4, 2007)
               hdlte1 = hdlte(i)/(HI(i,j,k)*dom)
               hdlow1 = max(hdlow(i), 1.0e-20)
               edot(i) = edot(i) - HDI(i,j,k)*
     .                     (hdlte(i)/(1.0 + hdlte1/hdlow1)) / (3.0*dom)
            end if
         enddo
      endif

!     --- Compute (external) radiative heating terms ---

!                       Photoionization heating

      if (iradshield .eq. 0) then

!        regular version

         if (iradtype .eq. 8) then

!           1) heating assuming high energy photons produces secondary
!              electrons which do the heating (Shull & Steenberg, 1985).

            do i = is+1, ie+1
               if ( itmask(i) ) then
               x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4)
               factor = 0.9971*(1.0-(1.0-x**0.2663)**1.3163)
               edot(i) = edot(i) + float(ipiht)*factor*(
     &                + piHI  *HI  (i,j,k)         ! pi of HI
     &                + piHeI *HeI (i,j,k)*0.25     ! pi of HeI
     &                + piHeII*HeII(i,j,k)*0.25     ! pi of HeII
     &              )/dom
               end if
            enddo

         else

!           2) standard heating

            do i = is+1, ie+1
               if ( itmask(i) ) then
               edot(i) = edot(i) + float(ipiht)*(
     &                + piHI  *HI  (i,j,k)         ! pi of HI
     &                + piHeI *HeI (i,j,k)*0.25     ! pi of HeI
     &                + piHeII*HeII(i,j,k)*0.25     ! pi of HeII
     &              )/dom
               end if
            enddo

         endif

      else

!        version with approximate self-shielding

         do i = is+1, ie+1
            if ( itmask(i) ) then
            edot(i) = edot(i) + float(ipiht)*(
     &                + piHI  *HI  (i,j,k)*
     &                   exp(-avgsighp*HI(i,j,k)*dom)
     &                + piHeI *HeI (i,j,k)*0.25*
     &                   exp(-avgsighep*0.25*HeI(i,j,k)*dom)
     &                + piHeII*HeII(i,j,k)*0.25*
     &                   exp(-avgsighe2p*0.25*HeII(i,j,k)*dom)
     &           )/dom
            end if
         enddo

      endif

!     Photoheating from radiative transfer

      if (iradtrans .eq. 1) then
c     for radiative transfer, convert from eV/s*TimeUnits to the coolunits use
         rtunits = 1.60217646d-12/utim/coolunit/dom
         do i = is+1, ie+1
            if (itmask(i)) then
               edot(i) = edot(i) + float(ipiht) * photogamma(i,j,k) * 
     &              rtunits * HI(i,j,k)
            endif
         enddo
      endif


!     --- Cooling/heating due to metals ---


!     --- Cloudy metal cooling and heating ---

      if (imcool .eq. 3) then

         call cool1D_cloudy(d, de, HI, HII, HeI, HeII, HeIII,
     &        HM, H2I, H2II, DI, DII, HDI, metal,
     &        in, jn, kn, is, ie, j, k,
     &        logtem, edot, comp2, ispecies, dom, zr,
     &        icmbTfloor, iClHeat, iClMMW,
     &        clMetNorm, clEleFra, clGridRank, clGridDim,
     &        clPar1, clPar2, clPar3, clPar4, clPar5,
     &        clDataSize, clCooling, clHeating, clMMW,
     &        itmask)

      endif

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) ) then
         tgasold(i) = tgas(i)
         end if
      enddo

      return
      end