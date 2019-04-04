




















!=======================================================================
!/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\
      subroutine solve_rate_cool(d, e, ge, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, imcool, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht,
     &                dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, fh, dtoh,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k22a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
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
     &                iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                kdissH2I, photogamma, 
     &                icmbTfloor, iClHeat, iClMMW,
     &                clMetNorm, clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating, clMMW)
!
!  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
!  modified2:  October, 1996 by GB; adapted to AMR
!  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
!  modified4:  June,    2005 by GB to solve rate & cool at same time
!  modified5:  September, 2009 by BDS to include cloudy cooling
!  modified6:  January, 2011 by DRR to include JHW's radiative transfer
!
!  PURPOSE:
!    Solve the multi-species rate and cool equations.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!    HI,HII   - H density fields (neutral & ionized)
!    HeI/II/III - He density fields
!    DI/II    - D density fields (neutral & ionized)
!    HDI      - neutral HD molecule density field
!    HM       - H- density field
!    H2I      - H_2 (molecular H) density field
!    H2II     - H_2+ density field
!    kph*     - photoionization fields
!    gamma*   - photoheating fields
!
!    is,ie    - start and end indices of active region (zero based)
!    idual    - dual energy formalism flag (0 = off, 1 = on)
!    iexpand  - comoving coordinates flag (0 = off, 1 = on)
!    idim     - dimensionality (rank) of problem
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    iradshield - flag for crude radiative shielding correction (rad type 12)
!    iradtype - type of radiative field (only used if = 8)
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!    imcool   - flag if there is metal cooling
!    imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
!    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
!    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
!    iradtrans - flag to include radiative transfer (1 = on, 0 = off)
!
!    fh       - Hydrogen mass fraction (typically 0.76)
!    dtoh     - Deuterium to H mass ratio
!    dt       - timestep to integrate over
!    aye      - expansion factor (in code units)
!
!    utim     - time units (i.e. code units to CGS conversion factor)
!    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
!    urho     - density units
!    uxyz     - length units
!    utem     - temperature(-like) units
!
!    temstart, temend - start and end of temperature range for rate table
!    nratec   - dimensions of chemical rate arrays (functions of temperature)
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    iClMMW     - flag to include addition to mean molecular weight from metals
!    clMetNorm  - parameter to convert metal density to metallicty (see CloudyCoolingData.h)
!    clEleFra   - parameter to account for additional electrons from metals 
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3, clPar4, clPar5 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clCooling  - cloudy cooling data
!    clHeating  - cloudy heating data
!    clMMW      - cloudy mean molecular weight data
!
!  OUTPUTS:
!    update chemical rate densities (HI, HII, etc)
!
!  PARAMETERS:
!    itmax   - maximum allowed sub-cycle iterations
!    mh      - H mass in cgs units
!
!-----------------------------------------------------------------------

      implicit NONE

!  General Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, imcool,
     &        iradtrans, irt_honly
      real    dt, aye, temstart, temend, eta1, eta2, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh

!  Density, energy and velocity fields fields

      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn)

!  Radiation fields

      real kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &     kdissH2I(in,jn,kn)
      real photogamma(in,jn,kn)

!  Cooling tables (coolings rates as a function of temperature)

      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)
      real    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p

!  Chemistry tables (rates as a function of temperature)

      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!  Cloudy cooling data

      integer icmbTfloor, iClHeat, iClMMW, clGridRank, clDataSize
      integer clGridDim(clGridRank)
      real clMetNorm, clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)), 
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)), 
     &     clPar5(clGridDim(5))
      real clCooling(clDataSize), clHeating(clDataSize), 
     &     clMMW(clDataSize)

!  Parameters

      integer itmax, ijk
      parameter (itmax = 10000, ijk = 4103)


      real tolerance
      parameter (tolerance = 1.0e-10)

      double precision mh
      parameter (mh = 1.67d-24)

!  Locals

      integer i, j, k, iter
      real ttmin, dom, energy, comp1, comp2
      double precision coolunit, dbase1, tbase1, xbase1

!  row temporaries

      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), 
     &     dtit(ijk), ttot(ijk), p2d(ijk), tgas(ijk), tgasold(ijk)

!  Rate equation row temporaries

      real HIp(ijk), HIIp(ijk), HeIp(ijk), HeIIp(ijk), HeIIIp(ijk),
     &     HMp(ijk), H2Ip(ijk), H2IIp(ijk),
     &     dep(ijk), dedot(ijk),HIdot(ijk), dedot_prev(ijk),
     &     DIp(ijk), DIIp(ijk), HDIp(ijk), HIdot_prev(ijk),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)
      real k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k50(ijk), k51(ijk), k52(ijk), k53(ijk), k54(ijk),
     &     k55(ijk), k56(ijk), k13dd(ijk, 7)

!  Cooling/heating row locals

      double precision ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk), edot(ijk)
      real hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk)

!  Iteration mask

      logical itmask(ijk)
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)

!  Convert densities from comoving to proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**(-3))

      call ceiling_species(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal)

!  Loop over zones, and do an entire i-column in one go

!$omp parallel
!$omp-  shared(d, e, ge, u, v, w, de)
!$omp-  shared(HI, HII, HeI, HeII, HeIII)
!$omp-  shared(in, jn, kn, nratec, iexpand, imethod)
!$omp-  shared(idual, ispecies, imetal, imcool, idim)
!$omp-  shared(is, js, ks, ie, je, ke, ih2co, ipiht)
!$omp-  shared(dt, aye, temstart, temend)
!$omp-  shared(utem, uxyz, uaye, urho, utim)
!$omp-  shared(eta1, eta2, gamma, fh, dtoh)
!$omp-  shared(k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a)
!$omp-  shared(k11a, k12a, k13a, k13dda, k14a, k15a)
!$omp-  shared(k16a, k17a, k18a, k19a, k22a)
!$omp-  shared(k24, k25, k26, k27, k28, k29, k30, k31)
!$omp-  shared(k50a, k51a, k52a, k53a, k54a, k55a, k56a)
!$omp-  shared(ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa)
!$omp-  shared(ciHeISa, ciHeIIa, reHIIa, reHeII1a)
!$omp-  shared(reHeII2a, reHeIIIa, brema, compa)
!$omp-  shared(comp_xraya, comp_temp, piHI, piHeI, piHeII)
!$omp-  shared(HM, H2I, H2II, DI, DII, HDI, metal)
!$omp-  shared(hyd01ka, h2k01a, vibha, rotha, rotla)
!$omp-  shared(gpldla, gphdla, hdltea, hdlowa)
!$omp-  shared(inutot, iradtype, nfreq, imetalregen)
!$omp-  shared(iradshield, avgsighp, avgsighep, avgsighe2p)
!$omp-  shared(dom, tbase1, xbase1, dbase1, coolunit)
!$omp-  shared(gaHIa, gaH2a, gaHea, gaHpa, gaela)
!$omp-  shared(icmbTfloor, iClHeat, iClMMW)
!$omp-  shared(clMetNorm, clEleFra, clGridRank, clGridDim)
!$omp-  shared(clPar1, clPar2, clPar3, clPar4, clPar5)
!$omp-  shared(clDataSize, clCooling, clHeating, clMMW)
!$omp-  shared(iradtrans, irt_honly)
!$omp-  shared(kphHI, kphHeI, kphHeII, kdissH2I, photogamma)
!$omp-  private(i, j, k, iter, itmask)
!$omp-  private(ttmin, energy, comp1, comp2)
!$omp-  private(indixe, t1, t2, logtem, tdef)
!$omp-  private(dtit, ttot, p2d, tgas, tgasold)
!$omp-  private(HIp, HIIp, HeIp, HeIIp, HeIIIp)
!$omp-  private(HMp, H2Ip, H2IIp)
!$omp-  private(dep, dedot,HIdot, dedot_prev)
!$omp-  private(DIp, DIIp, HDIp, HIdot_prev)
!$omp-  private(k24shield, k25shield, k26shield)
!$omp-  private(k1 , k2 , k3 , k4 , k5 )
!$omp-  private(k6 , k7 , k8 , k9 , k10)
!$omp-  private(k11, k12, k13, k14, k15)
!$omp-  private(k16, k17, k18, k19, k22)
!$omp-  private(k50, k51, k52, k53, k54)
!$omp-  private(k55, k56, k13dd)
!$omp-  private(ceHI, ceHeI, ceHeII)
!$omp-  private(ciHI, ciHeI, ciHeIS, ciHeII)
!$omp-  private(reHII, reHeII1, reHeII2, reHeIII)
!$omp-  private(brem, edot)
!$omp-  private(hyd01k, h2k01, vibh, roth, rotl)
!$omp-  private(gpldl, gphdl, hdlte, hdlow)
!$omp-  default(none)

!$omp do
      do k = ks+1, ke+1
      do j = js+1, je+1

!       tolerance = 1.0e-06 * dt

         do i = is+1, ie+1
            itmask(i) = .true.
         end do

!        Set time elapsed to zero for each cell in 1D section

         do i = is+1, ie+1
            ttot(i) = 0.0
         enddo

!        ------------------ Loop over subcycles ----------------

         do iter = 1, itmax

!           Compute the cooling rate (and tgas) for this row

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

!        Look-up rates as a function of temperature for 1D set of zones
!         (maybe should add itmask to this call)

            call lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &               is, ie, ijk, iradtype, iradshield, in, jn, kn,
     &               ispecies, tgas, HI, HII, HeI, HeII,
     &               k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &               k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &               k17a, k18a, k19a, k22a,
     &               k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &               avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &               k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &               k11, k12, k13, k14, k15, k16, k17, k18,
     &               k19, k22, k24, k25, k26, 
     &               k50, k51, k52, k53, k54, k55,
     &               k56, k13dd, k24shield, k25shield, k26shield,
     &               t1, t2, tdef, logtem, indixe, 
     &               dom, coolunit, tbase1, itmask)

!           Compute dedot and HIdot, the rates of change of de and HI
!             (should add itmask to this call)

            call rate_timestep(dedot, HIdot, ispecies,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask)

!           Find timestep that keeps relative chemical changes below 10%

            do i = is+1, ie+1
               if (itmask(i)) then
!              Bound from below to prevent numerical errors
               
	       if (abs(dedot(i)) .lt. 1.0e-20) 
     &             dedot(i) = min(1.0e-20,de(i,j,k))
	       if (abs(HIdot(i)) .lt. 1.0e-20)
     &             HIdot(i) = min(1.0e-20,HI(i,j,k))

!              If the net rate is almost perfectly balanced then set
!                  it to zero (since it is zero to available precision)

               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)),
     &                 abs(k2(i)*HII(i,j,k)*de(i,j,k)))/
     &             max(abs(dedot(i)),abs(HIdot(i))) .gt. 1.0e6) then
                  dedot(i) = 1.0e-20
                  HIdot(i) = 1.0e-20
               endif

!              If the iteration count is high then take the smaller of
!                the calculated dedot and last time step's actual dedot.
!                This is intended to get around the problem of a low
!                electron or HI fraction which is in equilibrium with high
!                individual terms (which all nearly cancel).

               if (iter .gt. 50) then
                  dedot(i) = min(abs(dedot(i)), abs(dedot_prev(i)))
                  HIdot(i) = min(abs(HIdot(i)), abs(HIdot_prev(i)))
               endif

!              compute minimum rate timestep

               dtit(i) = min(abs(0.1*de(i,j,k)/dedot(i)), 
     &                       abs(0.1*HI(i,j,k)/HIdot(i)),
     &                       dt-ttot(i), 0.5*dt)

!              Output some debugging information if required

            else               ! itmask
               dtit(i) = dt
            endif
            enddo               ! end loop over i

!           Compute maximum timestep for cooling/heating

            do i = is+1, ie+1
               if (itmask(i)) then
!              Set energy per unit volume of this cell based in the pressure
!              (the gamma used here is the right one even for H2 since p2d 
!               is calculated with this gamma).

               energy = max(p2d(i)/(gamma-1.0), 1.0e-20)

!              This is an alternate energy calculation, based directly on
!              the code's specific energy field, which differs from the above
!              only if using the dual energy formalism.

!              energy = max(ge(i,j,k)*d(i,j,k), p2d(i)/(gamma-1.0), 
!     &                    1.0e-20)
!              if (energy .lt. 1.0e-20) energy = d(i,j,k)*(e(i,j,k) - 
!     &              0.5*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))

!              If the temperature is at the bottom of the temperature look-up 
!              table and edot < 0, then shut off the cooling.

               if (tgas(i) .le. temstart .and. edot(i) .lt. 0.0) 
     &              edot(i) = 1.0e-20
	       if (abs(edot(i)) .lt. 1.0e-20) edot(i) = 1.0e-20

!              Compute timestep for 10% change

!              if (iter .gt. 100) then
!                 dtit(i) = min(real(abs(0.01*energy/edot(i))), 
!     &                        dt-ttot(i), dtit(i))
!              else
                  dtit(i) = min(real(abs(0.1*energy/edot(i))), 
     &                        dt-ttot(i), dtit(i))
!              endif


!              If the timestep is too small, then output some debugging info


 2000          format(4(i4,1x),1p,10(e14.3))
            endif   ! itmask
            enddo   ! end loop over i

!           Update total and gas energy

            do i = is+1, ie+1
               if (itmask(i)) then
               e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dtit(i)
               if (e(i,j,k) .ne. e(i,j,k)) 
     &              write(3,*) edot(i),d(i,j,k),dtit(i)

!              If using the dual energy formalism, there are 2 energy fields

               if (idual .eq. 1) then
                  ge(i,j,k) = ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i)

!                 Alternate energy update schemes (not currently used)

!                 ge(i,j,k) = max(ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i),
!     &                      0.5*ge(i,j,k))
!                 if (ge(i,j,k) .le. 1.0e-20) ge(i,j,k) = (energy + 
!     &              edot(i)*dtit(i))/d(i,j,k)
                  if (ge(i,j,k) .le. 0.0) write(3,*)
     &                 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
               endif
            endif               ! itmask
            enddo

!           Solve rate equations with backward differencing ---

            call step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask)

!           Add the timestep to the elapsed time for each cell and find
!            minimum elapsed time step in this row

            ttmin = 1.0e+20
            do i = is+1, ie+1
               ttot(i) = min(ttot(i) + dtit(i), dt)
               if (abs(dt-ttot(i)) .lt. 0.001*dt) itmask(i) = .false.
               if (ttot(i).lt.ttmin) ttmin = ttot(i)
            enddo

!           If all cells are done (on this slice), then exit

            if (abs(dt-ttmin) .lt. tolerance*dt) go to 9999

!           Next subcycle iteration

         enddo

 9999    continue

!       Abort if iteration count exceeds maximum

         if (iter .gt. itmax) then
            write(0,*) 'inside if statement solve rate cool:',is,ie
            write(6,*) 'MULTI_COOL iter > ',itmax,' at j,k =',j,k
            write(0,*) 'FATAL error (2) in MULTI_COOL'
            write(0,'(" dt = ",1pe10.3," ttmin = ",1pe10.3)') dt, ttmin
            write(0,'((16(1pe8.1)))') (dtit(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (ttot(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (edot(i),i=is+1,ie+1)
            write(0,'((16(l3)))') (itmask(i),i=is+1,ie+1)
            CALL f_error("OMP_solve_rate_cool.src",613)
         endif

         if (iter .gt. itmax/2) then
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
         end if
!     
!     Next j,k
!     
       enddo
      enddo
!$omp end do
!$omp end parallel

!     Convert densities back to comoving from proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**3)

!     Correct the species to ensure consistency (i.e. type conservation)

      call make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, d,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, fh, dtoh)

      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real     d(in,jn,kn),metal(in,jn,kn)
      real    factor

!     locals

      integer i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

!$omp parallel
!$omp-  shared(d, de, HI, HII, HeI, HeII, HeIII)
!$omp-  shared(HM, H2I, H2II, DI, DII, HDI, metal)
!$omp-  shared(is, ie, js, je, ks, ke)
!$omp-  shared(in, jn, kn, ispecies, imetal, factor)
!$omp-  private(i, j, k)
!$omp-  default(none)

!$omp do
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*factor
               de(i,j,k)    = de(i,j,k)*factor
               HI(i,j,k)    = HI(i,j,k)*factor
               HII(i,j,k)   = HII(i,j,k)*factor
               HeI(i,j,k)   = HeI(i,j,k)*factor
               HeII(i,j,k)  = HeII(i,j,k)*factor
               HeIII(i,j,k) = HeIII(i,j,k)*factor
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*factor
               enddo
            endif
         enddo
      enddo
!$omp end do
!$omp end parallel

      return
      end


c -----------------------------------------------------------
!   This routine ensures that the species aren't below tiny.

      subroutine ceiling_species(de, HI, HII, HeI, HeII, HeIII,
     &                           HM, H2I, H2II, DI, DII, HDI, metal,
     &                           is, ie, js, je, ks, ke,
     &                           in, jn, kn, ispecies, imetal)
c -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real metal(in,jn,kn)

!     locals

      integer i, j, k
      real small

      small = 1.0e-20

!$omp parallel
!$omp-  shared(de, HI, HII, HeI, HeII, HeIII)
!$omp-  shared(HM, H2I, H2II, DI, DII, HDI, metal)
!$omp-  shared(is, ie, js, je, ks, ke)
!$omp-  shared(in, jn, kn, ispecies, imetal, small)
!$omp-  private(i, j, k)
!$omp-  default(none)

!$omp do
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               de(i,j,k)    = max(de(i,j,k), small)
               HI(i,j,k)    = max(HI(i,j,k), small)
               HII(i,j,k)   = max(HII(i,j,k), small)
               HeI(i,j,k)   = max(HeI(i,j,k), small)
               HeII(i,j,k)  = max(HeII(i,j,k), small)
               HeIII(i,j,k) = max(HeIII(i,j,k), small)
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = max(HM(i,j,k), small)
                  H2I(i,j,k)  = max(H2I(i,j,k), small)
                  H2II(i,j,k) = max(H2II(i,j,k), small)
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = max(DI(i,j,k), small)
                  DII(i,j,k) = max(DII(i,j,k), small)
                  HDI(i,j,k) = max(HDI(i,j,k), small)
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = max(metal(i,j,k), small)
               enddo
            endif
         enddo
      enddo
!$omp end do
!$omp end parallel
      
      return
      end

! -----------------------------------------------------------
! This routine uses the temperature to look up the chemical
!   rates which are tabulated in a log table as a function
!   of temperature.

      subroutine lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &                is, ie, ijk, iradtype, iradshield, in, jn, kn,
     &                ispecies, tgas1d, HI, HII, HeI, HeII,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask)
! -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer is, ie, ijk, iradtype, iradshield, nratec, 
     &        in, jn, kn, ispecies, j, k
      real temstart, temend, tgas1d(in), dom
      double precision coolunit, tbase1
      logical itmask(ijk)

!     Chemistry rates as a function of temperature

      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7),
     &     k24, k25, k26,
     &     avgsighp, avgsighep, avgsighe2p, piHI, piHeI

!     Density fields

      real    HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn)

!     Returned rate values

      real k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k50(ijk), k51(ijk), k52(ijk), k53(ijk), k54(ijk),
     &     k55(ijk), k56(ijk), k13dd(ijk, 7),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)

!     1D temporaries (passed in)

      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk)

!     Parameters

      double precision everg, e24, e26
      parameter(everg = 1.60184d-12, e24 = 13.6d0, e26 = 24.6d0)

!     locals

      integer i, n1
      real factor, x, logtem0, logtem9, dlogtem, nh

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/real(nratec-1)

      do i = is+1, ie+1
         if (itmask(i)) then
!        Compute temp-centered temperature (and log)

!        logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
         logtem(i) = log(tgas1d(i))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Find index into tble and precompute interpolation values

         indixe(i) = min(nratec-1,
     &                max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)

!        Do linear table lookup (in log temperature)

         k1(i) = k1a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k1a(indixe(i)+1) -k1a(indixe(i)))/tdef(i)
         k2(i) = k2a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k2a(indixe(i)+1) -k2a(indixe(i)))/tdef(i)
         k3(i) = k3a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k3a(indixe(i)+1) -k3a(indixe(i)))/tdef(i)
         k4(i) = k4a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k4a(indixe(i)+1) -k4a(indixe(i)))/tdef(i)
         k5(i) = k5a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k5a(indixe(i)+1) -k5a(indixe(i)))/tdef(i)
         k6(i) = k6a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k6a(indixe(i)+1) -k6a(indixe(i)))/tdef(i)
      endif
      enddo

!     Look-up for 9-species model

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k7(i) = k7a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k7a(indixe(i)+1) -k7a(indixe(i)))/tdef(i)
            k8(i) = k8a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k8a(indixe(i)+1) -k8a(indixe(i)))/tdef(i)
            k9(i) = k9a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k9a(indixe(i)+1) -k9a(indixe(i)))/tdef(i)
            k10(i) = k10a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k10a(indixe(i)+1) -k10a(indixe(i)))/tdef(i)
            k11(i) = k11a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k11a(indixe(i)+1) -k11a(indixe(i)))/tdef(i)
            k12(i) = k12a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k12a(indixe(i)+1) -k12a(indixe(i)))/tdef(i)
            k13(i) = k13a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k13a(indixe(i)+1) -k13a(indixe(i)))/tdef(i)
            k14(i) = k14a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k14a(indixe(i)+1) -k14a(indixe(i)))/tdef(i)
            k15(i) = k15a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k15a(indixe(i)+1) -k15a(indixe(i)))/tdef(i)
            k16(i) = k16a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k16a(indixe(i)+1) -k16a(indixe(i)))/tdef(i)
            k17(i) = k17a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k17a(indixe(i)+1) -k17a(indixe(i)))/tdef(i)
            k18(i) = k18a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k18a(indixe(i)+1) -k18a(indixe(i)))/tdef(i)
            k19(i) = k19a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k19a(indixe(i)+1) -k19a(indixe(i)))/tdef(i)
            k22(i) = k22a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k22a(indixe(i)+1) -k22a(indixe(i)))/tdef(i)
         endif
         enddo

         do n1 = 1, 7
            do i = is+1, ie+1
            if (itmask(i)) then
               k13dd(i,n1) = k13dda(indixe(i),n1) + (logtem(i) - t1(i))
     &             *(k13dda(indixe(i)+1,n1) - 
     &               k13dda(indixe(i)  ,n1) )/tdef(i)
            endif
            enddo
         enddo
      endif

!     Look-up for 12-species model

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k50(i) = k50a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k50a(indixe(i)+1) -k50a(indixe(i)))/tdef(i)
            k51(i) = k51a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k51a(indixe(i)+1) -k51a(indixe(i)))/tdef(i)
            k52(i) = k52a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k52a(indixe(i)+1) -k52a(indixe(i)))/tdef(i)
            k53(i) = k53a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k53a(indixe(i)+1) -k53a(indixe(i)))/tdef(i)
            k54(i) = k54a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k54a(indixe(i)+1) -k54a(indixe(i)))/tdef(i)
            k55(i) = k55a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k55a(indixe(i)+1) -k55a(indixe(i)))/tdef(i)
            k56(i) = k56a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k56a(indixe(i)+1) -k56a(indixe(i)))/tdef(i)
            endif
         enddo
      endif

!        Include approximate self-shielding factors if requested

      do i = is+1, ie+1
         if (itmask(i)) then
         k24shield(i) = k24
         k25shield(i) = k25
         k26shield(i) = k26
         endif
      enddo
      if (iradshield .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k24shield(i) = k24shield(i)*exp(-HI(i,j,k)*avgsighp*dom)
            k25shield(i) = k25shield(i)*exp(-HeII(i,j,k)*avgsighe2p*dom)
            k26shield(i) = k26shield(i)*exp(-HeI(i,j,k)*avgsighep*dom)
            endif
         enddo
      endif

!        If using a high-energy radiation field, then account for
!          effects of secondary elections (Shull * Steenberg 1985)
!          (see calc_rate.src)

      if (iradtype .eq. 8) then
         do i = is+1, ie+1
            if (itmask(i)) then
            x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4)
            factor = 0.3908*(1.0 - x**0.4092)**1.7592
            k24shield(i) = k24shield(i) + 
     &         factor*(piHI + 0.08*piHeI)/(e24*everg) * coolunit*tbase1
            factor = 0.0554*(1.0 - x**0.4614)**1.6660
            k26shield(i) = k26shield(i) + 
     &         factor*(piHI/0.08 + piHeI)/(e26*everg) * coolunit*tbase1
            endif
         enddo
      endif


!           If using H2, and using the density-dependent collisional
!             H2 dissociation rate, then replace the the density-independant
!                k13 rate with the new one.
!         May/00: there appears to be a problem with the density-dependent
!             collisional rates.  Currently turned off until further notice.

            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  if (itmask(i)) then
                  nh = min(HI(i,j,k)*dom, 1.0e9)
                  k13(i) = 1.0e-20
                  if (tgas1d(i) .ge. 500.0 .and.
     &                tgas1d(i) .lt. 1.0e6) then
                     k13(i) = k13dd(i,1)-k13dd(i,2)/
     &                          (1.0+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1.0+(nh/k13dd(i,6))**k13dd(i,7))
                     k13(i) = max(10.0**k13(i), 1.0e-20)
                  endif
                  endif
               enddo
            endif

      return
      end

! -------------------------------------------------------------------
!  This routine calculates the electron and HI rates of change in
!    order to determine the maximum permitted timestep

      subroutine rate_timestep(dedot, HIdot, ispecies,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, 
     &                     kdissH2I, itmask)
! -------------------------------------------------------------------

      implicit NONE

!     arguments

      integer ispecies, is, ie, j, k, in, jn, kn, iradtrans, irt_honly
      real dedot(in), HIdot(in)
      logical itmask(in)

!     Density fields

      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)

!     Radiation fields

      real    kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

!     Rate values

      real k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k50(in), k51(in), k52(in), k53(in), k54(in),
     &     k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!     locals

      integer i

      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then

!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k3(i)*HeI(i,j,k)*de(i,j,k)/4.0
     &               + k5(i)*HeII(i,j,k)*de(i,j,k)/4.0
     &               - k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k4(i)*HeII(i,j,k)*de(i,j,k)/4.0
     &               - k6(i)*HeIII(i,j,k)*de(i,j,k)/4.0
     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4.0
     &               + k26shield(i)*HeI(i,j,k)/4.0)

!     Compute the HI density rate-of-change

            HIdot(i) =
     &               - k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k2(i)*HII(i,j,k)*de(i,j,k)
     &               -      k24shield(i)*HI(i,j,k)

         endif                  ! itmask
         enddo
      else

!         Include molecular hydrogen rates for HIdot

         do i = is+1, ie+1
            if (itmask(i)) then
            HIdot(i) = 
     &               -    k1(i) *de(i,j,k)    *HI(i,j,k)  
     &               -    k7(i) *de(i,j,k)    *HI(i,j,k)
     &               -    k8(i) *HM(i,j,k)    *HI(i,j,k)
     &               -    k9(i) *HII(i,j,k)   *HI(i,j,k)
     &               -    k10(i)*H2II(i,j,k)  *HI(i,j,k)/2.
     &               - 2.*k22(i)*HI(i,j,k)**2 *HI(i,j,k)
     &               +    k2(i) *HII(i,j,k)   *de(i,j,k) 
     &               + 2.*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2.
     &               +    k11(i)*HII(i,j,k)   *H2I(i,j,k)/2.
     &               + 2.*k12(i)*de(i,j,k)    *H2I(i,j,k)/2.
     &               +    k14(i)*HM(i,j,k)    *de(i,j,k)
     &               +    k15(i)*HM(i,j,k)    *HI(i,j,k)
     &               + 2.*k16(i)*HM(i,j,k)    *HII(i,j,k)
     &               + 2.*k18(i)*H2II(i,j,k)  *de(i,j,k)/2.
     &               +    k19(i)*H2II(i,j,k)  *HM(i,j,k)/2.
     &               -      k24shield(i)*HI(i,j,k)

!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i) * HI(i,j,k)   * de(i,j,k)
     &               + k3(i) * HeI(i,j,k)  * de(i,j,k)/4.
     &               + k5(i) * HeII(i,j,k) * de(i,j,k)/4.
     &               + k8(i) * HM(i,j,k)   * HI(i,j,k)
     &               + k15(i)* HM(i,j,k)   * HI(i,j,k)
     &               + k17(i)* HM(i,j,k)   * HII(i,j,k)
     &               + k14(i)* HM(i,j,k)   * de(i,j,k)
     &               - k2(i) * HII(i,j,k)  * de(i,j,k)
     &               - k4(i) * HeII(i,j,k) * de(i,j,k)/4.
     &               - k6(i) * HeIII(i,j,k)* de(i,j,k)/4.
     &               - k7(i) * HI(i,j,k)   * de(i,j,k)
     &               - k18(i)* H2II(i,j,k) * de(i,j,k)/2.

     &               + (k24shield(i)*HI(i,j,k)
     &               +  k25shield(i)*HeII(i,j,k)/4.0
     &               +  k26shield(i)*HeI(i,j,k)/4.0)
         endif                  ! itmask
         enddo
      endif

!     Add photo-ionization rates if needed

      if (iradtrans .eq. 1) then
         if (irt_honly .eq. 0) then
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     &                 + kphHeI(i,j,k) * HeI(i,j,k) / 4.0
     &                 + kphHeII(i,j,k) *HeII(i,j,k) / 4.0
               endif
            enddo
         else
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
               endif
            enddo
         endif
      endif

      return
      end


! -----------------------------------------------------------
!  This routine uses a backward-finite difference time integrator
!   to advance the rate equations by one (sub-)cycle (dtit).

      subroutine step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII,
     &                     kdissH2I, itmask)
c -------------------------------------------------------------------

      implicit NONE

!     arguments

      integer ispecies, in, jn, kn, is, ie, j, k, iradtrans, irt_honly
      real    dtit(in), dedot_prev(in), HIdot_prev(in)
      logical itmask(in)

!     Density fields

      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)

!     Radiation fields

      real    kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

!     Rate values

      real k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k50(in), k51(in), k52(in), k53(in), k54(in),
     &     k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!     temporaries (passed in)

      real HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &     HMp(in), H2Ip(in), H2IIp(in), dep(in),
     &     DIp(in), DIIp(in), HDIp(in)

!     locals

      integer i
      real scoef, acoef

!   A) the 6-species integrator
!      
      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI

            scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
            acoef  = k1(i)*de(i,j,k)
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
            HIp(i)  = (scoef*dtit(i) + HI(i,j,k))/(1.0 + acoef*dtit(i))

!        2) HII
c 
            scoef  = k1(i)*HIp(i)*de(i,j,k)
     &                   + k24shield(i)*HIp(i)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*HIp(i)
            acoef  = k2(i)*de (i,j,k)
            HIIp(i) = (scoef*dtit(i) + HII(i,j,k))/(1.0 + acoef*dtit(i))

!                 3) Electron density

            scoef = 0.0
     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4.
     &                 + k26shield(i)*HeI(i,j,k)/4.
            if (iradtrans .eq. 1 .and. irt_honly .eq. 0) 
     &           scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
     &           + kphHeI(i,j,k)  * HeI(i,j,k)/4.
     &           + kphHeII(i,j,k) * HeII(i,j,k)/4.
            if (iradtrans .eq. 1 .and. irt_honly .eq. 1) 
     &           scoef = scoef + kphHI(i,j,k)   * HI(i,j,k)
            acoef = -(k1(i)*HI(i,j,k)      - k2(i)*HII(i,j,k)
     &              + k3(i)*HeI(i,j,k)/4.  - k6(i)*HeIII(i,j,k)/4.0
     &              + k5(i)*HeII(i,j,k)/4. - k4(i)*HeII(i,j,k)/4.0)
            dep(i)   = (scoef*dtit(i) + de(i,j,k))
     &                     / (1.0 + acoef*dtit(i))

         endif                  ! itmask
         enddo

      endif                     ! (ispecies .eq. 1)

!  --- (B) Do helium chemistry in any case: (for all ispecies values) ---

      do i = is+1, ie+1
         if (itmask(i)) then

!        4) HeI

         scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
         acoef  = k3(i)*de(i,j,k)
     &                + k26shield(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        acoef = acoef + kphHeI(i,j,k)
         HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &              / ( 1. + acoef*dtit(i) )

!        5) HeII

         scoef  = k3(i)*HeIp(i)*de(i,j,k)
     &          + k6(i)*HeIII(i,j,k)*de(i,j,k)
     &          + k26shield(i)*HeIp(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        scoef = scoef + kphHeI(i,j,k)*HeIp(i)
         acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)
     &          + k25shield(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        acoef = acoef + kphHeII(i,j,k)
         HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &              / ( 1. + acoef*dtit(i) )

!       6) HeIII

         scoef   = k5(i)*HeIIp(i)*de(i,j,k)
     &           + k25shield(i)*HeIIp(i)
         if (iradtrans.eq.1 .and. irt_honly.eq.0) 
     $        scoef = scoef + kphHeII(i,j,k) * HeIIp(i)
         acoef   = k6(i)*de(i,j,k)
         HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                / ( 1. + acoef*dtit(i) )

      endif                     ! itmask
      enddo

c --- (C) Now do extra 3-species for molecular hydrogen ---

      if (ispecies .gt. 1) then

!        First, do HI/HII with molecular hydrogen terms

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI
!     
            scoef  =    k2(i) * HII(i,j,k) * de(i,j,k) 
     &             + 2.*k13(i)* HI(i,j,k)  * H2I(i,j,k)/2.
     &             +    k11(i)* HII(i,j,k) * H2I(i,j,k)/2.
     &             + 2.*k12(i)* de(i,j,k)  * H2I(i,j,k)/2.
     &             +    k14(i)* HM(i,j,k)  * de(i,j,k)
     &             +    k15(i)* HM(i,j,k)  * HI(i,j,k)
     &             + 2.*k16(i)* HM(i,j,k)  * HII(i,j,k)
     &             + 2.*k18(i)* H2II(i,j,k)* de(i,j,k)/2.
     &             +    k19(i)* H2II(i,j,k)* HM(i,j,k)/2.
     &             + 2.*k31   * H2I(i,j,k)/2.
            if (iradtrans.eq.1)
     &           scoef = scoef + 2.*kdissH2I(i,j,k)* H2I(i,j,k)/2.
            acoef  =    k1(i) * de(i,j,k)
     &             +    k7(i) * de(i,j,k)  
     &             +    k8(i) * HM(i,j,k)
     &             +    k9(i) * HII(i,j,k)
     &             +    k10(i)* H2II(i,j,k)/2.
     &             + 2.*k22(i)* HI(i,j,k)**2
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1. + acoef*dtit(i) )

!          2) HII

            scoef  =    k1(i)  * HI(i,j,k) * de(i,j,k)
     &             +    k10(i) * H2II(i,j,k)*HI(i,j,k)/2.
     &             + k24shield(i)*HI(i,j,k)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*HI(i,j,k)
            acoef  =    k2(i)  * de(i,j,k)
     &             +    k9(i)  * HI(i,j,k)
     &             +    k11(i) * H2I(i,j,k)/2.
     &             +    k16(i) * HM(i,j,k)
     &             +    k17(i) * HM(i,j,k)
            HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1. + acoef*dtit(i) )
!     
!          3) electrons:

            scoef =   k8(i) * HM(i,j,k) * HI(i,j,k)
     &             +  k15(i)* HM(i,j,k) * HI(i,j,k)
     &             +  k17(i)* HM(i,j,k) * HII(i,j,k)
!                  
     &             + k24shield(i)*HIp(i)
     &             + k25shield(i)*HeIIp(i)/4.
     &             + k26shield(i)*HeIp(i)/4.
            if (iradtrans .eq. 1 .and. irt_honly.eq.0) 
     &           scoef = scoef + kphHI(i,j,k)   * HIp(i)
     &           + kphHeI(i,j,k)  * HeIp(i)/4.
     &           + kphHeII(i,j,k) * HeIIp(i)/4.
            if (iradtrans .eq. 1 .and. irt_honly.eq.1) 
     &           scoef = scoef + kphHI(i,j,k)   * HIp(i)
            acoef = - (k1(i) *HI(i,j,k)    - k2(i)*HII(i,j,k)
     &              +  k3(i) *HeI(i,j,k)/4. - k6(i)*HeIII(i,j,k)/4.
     &              +  k5(i) *HeII(i,j,k)/4.- k4(i)*HeII(i,j,k)/4.
     &              +  k14(i)*HM(i,j,k)
     &              -  k7(i) *HI(i,j,k)
     &              -  k18(i)*H2II(i,j,k)/2.0)
            dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                / ( 1. + acoef*dtit(i) )

!           7) H2

            scoef = 2.*(k8(i)  * HM(i,j,k)   * HI(i,j,k)
     &            +     k10(i) * H2II(i,j,k) * HI(i,j,k)/2.
     &            +     k19(i) * H2II(i,j,k) * HM(i,j,k)/2.
     &            +     k22(i) * HI(i,j,k)**3)
            acoef = ( k13(i)*HI(i,j,k) + k11(i)*HII(i,j,k)
     &              + k12(i)*de(i,j,k) )
     &              + k29 + k31
            if (iradtrans.eq.1) acoef = acoef + kdissH2I(i,j,k)

            H2Ip(i) = ( scoef*dtit(i) + H2I(i,j,k) )
     &                / ( 1. + acoef*dtit(i) )

!           8) H-

            HMp(i) = ( k7(i)*HIp(i)*dep(i) )
     &             / ( (k8(i)+k15(i))*HIp(i)
     &             + ( k16(i)+k17(i))*HIIp(i)+k14(i)*dep(i)
     &             + k27
     &                 )

!           9) H2+

            H2IIp(i) = 2.*( k9 (i)*HIp(i)*HIIp(i)
     &                    + k11(i)*H2Ip(i)/2.0*HIIp(i)
     &                    + k17(i)*HMp(i)*HIIp(i)
     &                    + k29*H2Ip(i)
     &                    )
     &                 /  ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                    + k19(i)*HMp(i)
     &                    + (k28+k30)
     &                    )

         endif                  ! itmask
         enddo
!     
      endif                     ! H2

!  --- (D) Now do extra 3-species for molecular HD ---
!     
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
!     
!         1) DI
!     
            scoef =   (     k2(i) * DII(i,j,k) * de(i,j,k)
     &                 +    k51(i)* DII(i,j,k) * HI(i,j,k)
     &                 + 2.*k55(i)* HDI(i,j,k) * HI(i,j,k)/3.
     &                 )
            acoef  =    k1(i) * de(i,j,k)
     &             +    k50(i) * HII(i,j,k)
     &             +    k54(i) * H2I(i,j,k)/2.
     &             +    k56(i) * HM(i,j,k)
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
            DIp(i)    = ( scoef*dtit(i) + DI(i,j,k) ) / 
     &                  ( 1. + acoef*dtit(i) )

!         2) DII
c 
            scoef =   (   k1(i)  * DI(i,j,k) * de(i,j,k)
     &            +       k50(i) * HII(i,j,k)* DI(i,j,k)
     &            +    2.*k53(i) * HII(i,j,k)* HDI(i,j,k)/3.
     &            )
     &            + k24shield(i)*DI(i,j,k)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*DI(i,j,k)
            acoef =    k2(i)  * de(i,j,k)
     &            +    k51(i) * HI(i,j,k)
     &            +    k52(i) * H2I(i,j,k)/2.

            DIIp(i)   = ( scoef*dtit(i) + DII(i,j,k) )
     &                 / ( 1. + acoef*dtit(i) )

!          3) HDI
c 
            scoef = 3.*(k52(i) * DII(i,j,k)* H2I(i,j,k)/2./2.
     &             +    k54(i) * DI(i,j,k) * H2I(i,j,k)/2./2.
     &             + 2.*k56(i) * DI(i,j,k) * HM(i,j,k)/2.
     &                 )
            acoef  =    k53(i) * HII(i,j,k)
     &             +    k55(i) * HI(i,j,k)

            HDIp(i)   = ( scoef*dtit(i) + HDI(i,j,k) )
     &                 / ( 1. + acoef*dtit(i) )

         endif                  ! itmask
         enddo
      endif

!   --- (E) Set densities from 1D temps to 3D fields ---

      do i = is+1, ie+1
         if (itmask(i)) then
         HIdot_prev(i) = abs(HI(i,j,k)-HIp(i))/max(dtit(i),1.0e-20)
         HI(i,j,k)    = max(HIp(i), 1.0e-20)
         HII(i,j,k)   = max(HIIp(i), 1.0e-20)
         HeI(i,j,k)   = max(HeIp(i), 1.0e-20)
         HeII(i,j,k)  = max(HeIIp(i), 1.0e-20)
         HeIII(i,j,k) = max(HeIIIp(i), 1.0e-20)

!        de(i,j,k)    = dep(i)

!        Use charge conservation to determine electron fraction

         dedot_prev(i) = de(i,j,k)
         de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4. + HeIII(i,j,k)/2.
         if (ispecies .gt. 1) 
     &      de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2.
         dedot_prev(i) = abs(de(i,j,k)-dedot_prev(i))/
     &                         max(dtit(i),1.0e-20)

         if (ispecies .gt. 1) then
            HM(i,j,k)    = max(HMp(i), 1.0e-20)
            H2I(i,j,k)   = max(H2Ip(i), 1.0e-20)
            H2II(i,j,k)  = max(H2IIp(i), 1.0e-20)
         endif

         if (ispecies .gt. 2) then
            DI(i,j,k)    = max(DIp(i), 1.0e-20)
            DII(i,j,k)   = max(DIIp(i), 1.0e-20)
            HDI(i,j,k)   = max(HDIp(i), 1.0e-20)
         endif
      endif                     ! itmask
!     
      enddo                     ! end loop over i

      return
      end

! ------------------------------------------------------------------
!   This routine correct the highest abundence species to
!     insure conservation of particle number and charge.

      subroutine make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, fh, dtoh)
! -------------------------------------------------------------------

      implicit NONE

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)  
      real    fh, dtoh

!     Parameters

      integer ijk
      parameter (ijk = 4103)

!     locals

      integer i, j, k
      real totalH(ijk), totalHe(ijk), 
     &     totalD, correctH, correctHe, correctD

!     Loop over all zones

!$omp parallel
!$omp-  shared(de, HI, HII, HeI, HeII, HeIII)
!$omp-  shared(HM, H2I, H2II, DI, DII, HDI, d)
!$omp-  shared(is, js, ks, ie, je, ke)
!$omp-  shared(in, jn, kn, ispecies, fh, dtoh)
!$omp-  private(i, j, k)
!$omp-  private(totalH, totalHe, totalD)
!$omp-  private(correctH, correctHe, correctD)
!$omp-  default(none)

!$omp do
      do k = ks+1, ke+1
      do j = js+1, je+1

!     Compute total densities of H and He
!         (ensure non-negativity)

      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
      enddo

!     include molecular hydrogen

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
         enddo
      endif

!     Correct densities by keeping fractions the same

      do i = is+1, ie+1
         correctH = fh*d(i,j,k)/totalH(i)
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH

         correctHe = (1.0 - fh)*d(i,j,k)/totalHe(i)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe

!     Correct molecular hydrogen-related fractions

         if (ispecies .gt. 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) + 2.0/3.0*HDI(i,j,k)
            correctD = fh*dtoh*d(i,j,k)/totalD
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4. + HeIII(i,j,k)/2.
         if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &             - HM(i,j,k) + H2II(i,j,k)/2.
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k
!$omp end do
!$omp end parallel

      return
      end