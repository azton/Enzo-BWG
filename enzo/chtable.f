




















c=======================================================================
c//////////////////////  SUBROUTINE CHTABLE  \\\\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine chtable(NFREQ, FREQDEL, SIGH, SIGHE, SIGHE2,
     &                   SPSTAR1, SPSTAR2)

c
c  GENERATES THE LOOK-UP TABLE USED FOR MANY ATOMIC RATES
c       (CURRENTLY ONLY RADIATIVE)
c
c  written by: Renyue Cen
c  date:       
c  modified1:  September, 1999 by Greg Bryan; converted to AMR
c  modified2:  February, 2004 by Robert Harkness
c              Remove obsolete syntax
c
c  PURPOSE:
C     THIS ROUTINE CREATES A LOOK UP TABLE FOR COLLISIONAL IONIZATION
C     AND RECOMBINATION COEFFICIENTS, AND COOLING AND HEATING TERMS,
C     EXCLUDING PHOTOIONIZATION AND PHOTOHEATING, COMPTON HEATING/COOLING
C   
C     THE LOOKUP TABLE IS A FUNCTION OF TEMPERATURE ONLY
c
c  INPUTS:
c    NFREQ    - Number of frequency bins
c    FREQDEL  - space between frequency bins, in log10(eV)
c
c  OUTPUTS:
c    SIGH     - HI photo-ionization heating cross-section
c    SIGHE    - HeI photo-ionization heating cross-section
c    SIGHE2   - HeII photo-ionization heating cross-section
c    SPSTAR1  - normalized shape of stellar radiation field
c    SPSTAR2  - normalized shape of quasar radiation field
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  Arguments
c
      integer NFREQ
      real    FREQDEL, SIGH(NFREQ), SIGHE(NFREQ), SIGHE2(NFREQ),
     &        SPSTAR1(NFREQ), SPSTAR2(NFREQ)
c
c  Parameters
c
      real    EV2HZ, EV2ERG, PI
c
c  Locals
c
      integer N, NTI
      real    E0, sigma0, ya, P, yw, y0, y1, x, y, Fy, FNU, FNU2,
     &        SPTMO, SUM, TBB, RAT, FLAMRAT, FNURAT, ALA, A,
     &        FNUPOWBK1, FNUPOWBK2, SPITP(9), FNUITP(9), WORK(9),
     &        FNU1
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set some cosntants
c
      PI        = 4.*ATAN(1.)
      EV2HZ     = 2.415E14 
      EV2ERG    = 1.6E-12
C
c-----------------------------------------------------------------------
C     CREATE THE CROSS SECTIONS
C
      E0 = 13.61 
      sigma0 = 9.492e-16
      ya = 1.469
      P  = 3.188
      yw = 2.039
      y0 = 0.4434
      y1 = 2.136
      DO N=1,NFREQ
         FNU  = 10.0**((N-1.0)*FREQDEL)
C  
C        i) HYDROGEN I PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU .GT. 13.6) THEN
           FNU1    = SQRT(FNU/13.6 -1.0)
           SIGH(N) = 1.18E-11*FNU**(-4)*EXP(-4.*ATAN(FNU1)/FNU1)
     .                       /(1.0 -EXP(-2.*PI/FNU1))
         ELSE 
           IF(FNU .EQ. 13.6) THEN
             SIGH(N) = 2.16E-13*FNU**(-4)
           ELSE
             SIGH(N) = 0.
           ENDIF
         ENDIF
C
C
C        ii) HELIUM I PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU.GE.24.59) THEN
c          SIGHE(N) = 1.13E-14*(FNU**(-2.05) -9.775*FNU**(-3.05))
           x       = FNU/E0-y0
           y       = sqrt(x*x+y1*y1)
           Fy      = ((x-1.0)**2+yw**2)*y**(0.5*P-5.5)
     .             *(1.0+sqrt(y/ya))**(-P)
           SIGHE(N) = sigma0*Fy
         ELSE
           SIGHE(N) = 0.
         ENDIF
C
C
C        iii) HELIUM II PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU.GT.54.4) THEN
           FNU2      = SQRT(FNU/54.4 -1.0)
           SIGHE2(N) = 7.55E-10*FNU**(-4)*EXP(-4.*ATAN(FNU2)/FNU2)
     .                         /(1.0 -EXP(-2.*PI/FNU2))
         ELSE
           IF(FNU.EQ.54.4) THEN
             SIGHE2(N) = 1.38E-11*FNU**(-4)
           ELSE
             SIGHE2(N) = 0.
           ENDIF
         ENDIF
      ENDDO
C
c-----------------------------------------------------------------------
C     GENERATE PHOTOIONIZATION SPECTRUM FOR MASSIVE STARS
C
      FNUITP(1)= 10.00
      SPITP(1) = 0.291

      FNUITP(2)= 12.60
      SPITP(2) = 0.287

      FNUITP(3)= 15.85
      SPITP(3) = 0.283

      FNUITP(4)= 19.95
      SPITP(4) = 0.280

      FNUITP(5)= 25.12
      SPITP(5) = 0.210

      FNUITP(6)= 30.63
      SPITP(6) = 0.130

      FNUITP(7)= 39.81
      SPITP(7) = 0.081

      FNUITP(8)= 50.12
      SPITP(8) = 0.018

      FNUITP(9)= 60.00
      SPITP(9) = 0.010
C
      DO N=1,NFREQ
         FNU = 10.0**((N-1.0)*FREQDEL)
         IF(FNU.LT.50.12) THEN
           CALL IN_EXTP(FNUITP,SPITP,9,FNU,SPTMO,WORK)
           SPSTAR1(N) = SPTMO
         ELSE
           SPSTAR1(N) = 10.0**(-(FNU-60.00)/(60.00-50.12)
     .                        *(LOG10(0.018)-LOG10(0.010))
     .                        + LOG10(0.010))
         ENDIF
      ENDDO
C
      SUM = 0.0
      DO N=1,NFREQ
         FNU = 10.0**((N-1.0)*FREQDEL)
         IF(FNU.GE.13.6) THEN
           SUM = SUM + SPSTAR1(N)*EV2HZ*FNU
     .                *(10.0**(0.5*FREQDEL)-10**(-0.5*FREQDEL))
         ENDIF
      ENDDO
C
      OPEN(15,FILE='STAR.sp')
      DO N=1,NFREQ
         FNU = 10.0**((N-1.0)*FREQDEL)
C
C        SPSTAR1(N) IS NORMALIZED SUCH THAT 
C        \int_(13.6 eV)^{\infty} SPSTAR1(N) dNU = 1 erg
C        WHERE dNU is in hz
         SPSTAR1(N) = SPSTAR1(N)/SUM
         WRITE(15,505)N,FNU,SPSTAR1(N),SPSTAR1(N)*FNU
     .                  ,SIGH(N),SIGHE(N),SIGHE2(N)
      ENDDO
      CLOSE(15)
C
  505 FORMAT(I4,6(1X,E12.4))
C
c-----------------------------------------------------------------------
C     GENERATE PHOTOIONIZATION SPECTRUM FOR QUASAR
C
C     BLACK BODY TEMPERATURE IN KELVIN
c     TBB     = 2.6e4
      TBB     = 3.0e4
C     BLACK BODY TEMPERATURE IN eV
      TBB     = Tbb/1.16e4
      RAT     = 0.3
      FLAMRAT = 5450.0
c     h\nu IN eV
      FNURAT  = 6.625e-27*3.e10/(FLAMRAT*1.e-8)/(1.6e-12)
c
      ALA     = -1.3
c
      A       = RAT/FNURAT**(3.0-ALA)*(EXP(FNURAT/TBB)-1.0) 
c
c
      FNUPOWBK1= 1.E3
      FNUPOWBK2= 1.E5
c
      DO N=1,NFREQ
         FNU = 10.0**((N-1.0)*FREQDEL)
         IF(FNU.LT.FNUPOWBK1) THEN
           SPSTAR2(N) = A*FNU**3/(exp(min(50.0,FNU/Tbb))-1.0) 
     .                 + FNU**ala
         ELSE
           IF(FNU.LT.FNUPOWBK2) THEN
             SPSTAR2(N) = FNU**ala*(1.0+(FNU/FNUPOWBK1)**0.7)/2.0
           ELSE
             SPSTAR2(N) = FNU**ala*(1.0+(FNU/FNUPOWBK1)**0.7)/2.0
     .                /(1.0+(FNU/FNUPOWBK2)**1.0)*2.0
             ENDIF
         ENDIF
      ENDDO
C
      SUM = 0.0
      DO N=1,NFREQ
         FNU = 10.0**((N-1.0)*FREQDEL)
         IF(FNU.GE.13.6) THEN
           SUM = SUM + SPSTAR2(N)*EV2HZ*FNU
     .                *(10.0**(0.5*FREQDEL)-10**(-0.5*FREQDEL))
         ENDIF
      ENDDO
C
      OPEN(15,FILE='QUASAR.sp')
      DO N=1,NFREQ
         FNU        = 10.0**((N-1.0)*FREQDEL)
         SPSTAR2(N) = SPSTAR2(N)/SUM
C
C        SPSTAR2(N) IS NORMALIZED SUCH THAT 
C        \int_(13.6 eV)^{\infty} SPSTAR2(N) dNU = 1 erg
C        WHERE dNU is in hz
         WRITE(15,505)N,FNU,SPSTAR2(N),SPSTAR2(N)*FNU
     .                 ,SIGH(N),SIGHE(N),SIGHE2(N)
      ENDDO
      CLOSE(15)
C
C
c-----------------------------------------------------------------------
C     NOW GENERATE TABLE AS A FUNCTION OF TEMPERATURE 
C     FOR IONIZATION COEFFICIENTS AND COOLING/HEATING 
C     FOR A PLASMA OF HYDROGEN AND HELIUM WITH PRIMORDIAL COMPOSITION 
C
C
C
      RETURN
      END
C
C
C
C
C     NAME IN_EXTP
      SUBROUTINE IN_EXTP(XN,YN,N,X,Y,XX)
C
C     INPUT  : XN(N),YN(N),N,XX(N),X
C              WHERE XX(N) IS A WORKING ARRAY
C     OUTPUT : Y
C
C     PURPOSE: GIVEN N VALUES AT N POINTS, THIS ROUTINE INTERPOLATES OT
C              EXTRAPOLATES TO GIVE THE VALUE Y AT X
C
C     WRITTEN: BY RENYUE CEN, 11/7/92
C
C
      DIMENSION XN(N),YN(N),XX(N)
      DIMENSION XNTMP(10),YNTMP(10)
C
C
      IF(X.GE.XN(2) .AND. X.LE.XN(N-1)) THEN
        DIRV1 = (YN(2)-YN(1))/(XN(2)-XN(1))
        DIRVN = (YN(N)-YN(N-1))/(XN(N)-XN(N-1))
        CALL SPLINE(XN,YN,N,DIRV1,DIRVN,XX)
        CALL SPLINT(XN,YN,XX,N,X,Y)
      ELSE
        IF(X.LT.XN(2)) THEN
          IF(N.GE.10) THEN
            NPOLY = 10
          ELSE
            NPOLY = N
          ENDIF
          DO I=1,NPOLY
             XNTMP(I) = XN(I)
             YNTMP(I) = YN(I)
          ENDDO
          CALL RATINT(XNTMP,YNTMP,NPOLY,X,Y,DY)
        ELSE
          IF(N.GE.10) THEN
            NPOLY  = 10
            NSHIFT = N - 10
            DO I=1,NPOLY
               XNTMP(I) = XN(NSHIFT+I)
               YNTMP(I) = YN(NSHIFT+I)
            ENDDO
          ELSE
            NPOLY  = N
            DO I=1,NPOLY
               XNTMP(I) = XN(I)
               YNTMP(I) = YN(I)
            ENDDO
          ENDIF
          CALL RATINT(XNTMP,YNTMP,NPOLY,X,Y,DY)
        ENDIF
      ENDIF
c
c
      RETURN
      END
c
c
c
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=1000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
c
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)

      if (h.eq.0.) then
        write(0,'("Bad XA input in SPLINT")')
        CALL f_error("chtable.src",509)
      end if

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
C
C
      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)

          if (dd.eq.0.) then
            write(0,'("DD = 0 in RATINT")')
            CALL f_error("chtable.src",549)
          end if

          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
