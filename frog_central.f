C To compile -
C     f77 -o frog -O -u frog.f
C
C Program to integrate test particle orbits 
C in potential of a central star.
C
      CALL get_data
C          --------
      CALL integrate
C          ---------
      CALL check
C          -----
      END

CC***********************************************************************
C      SUBROUTINE accel(j)
CC***********************************************************************
C      INCLUDE 'frog.h'
C      INTEGER i,j
C      REAL*8 r,ar, ar_km_s2
C      
C      DO 10 i=1,nbods
CC Acceleration due to Sun
C         r=SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
C         ar=-G*msun/(r*r)
CC
CC        acceleration in km/s^2: 1kpc=3.0857e16km
C         ar_km_s2= ar/kpc2km
C         ax(i)=ar_km_s2*x(i)/r
C         ay(i)=ar_km_s2*y(i)/r
C         az(i)=ar_km_s2*z(i)/r
CC
C         pot(i)=-G*msun/r
CC        print *, 'r=',r,'  ar=',ar,'  pot(i)=',pot(i)
C        write(*,991) 'ar=',ar,'  ax=',ax(i),'  ay=',ay(i),'  az=',az(i)
C 991    format(4(A,e12.4))
CC
C 10   CONTINUE
C      RETURN
C      END

C***********************************************************************
      SUBROUTINE accel(j)
C***********************************************************************
C     Pal5: 
      INCLUDE 'frog.h'
      INTEGER i,j
      REAL*8 r,aror,c,Mhalo,phi0,p,mr
      SAVE c, Mhalo, phi0
      c = rvir/rs
      Mhalo = Mvir/(dlog(c+1.)-c/(c+1.))
      phi0 = G*Mhalo/rs
      
      DO 10 i=1,nbods
C NFW potential
         r = SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         p = r/rs
         mr = Mhalo*(dlog(p+1.)-p/(p+1.))
C truncation of NFW halo
         if (r.gt.rvir) mr = Mhalo*(dlog(c+1.)-c/(c+1.))
         aror = -(G*mr/r**2)
         pot(i) = -(phi0/p)*dlog(p+1.)
C
         aror_km_s2 = aror/kpc2km
         ax(i)=aror_km_s2*x(i)/r
         ay(i)=aror_km_s2*y(i)/r
         az(i)=aror_km_s2*z(i)/r
C        write(*,991) 'ar=',ar,'  ax=',ax(i),'  ay=',ay(i),'  az=',az(i)
C 991    format(4(A,e12.4))
C
 10   CONTINUE
      RETURN
      END


C***********************************************************************
      SUBROUTINE check
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 eps
      PARAMETER(eps=0.01d0)
      
      DO 10 i=1,nbods
         IF(maxe(i).GT.eps)WRITE(6,*)i,maxe(i)
 10   CONTINUE

      RETURN
      END

C***********************************************************************
      SUBROUTINE get_data
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i

      OPEN(UNIT=8,FILE='frogin_central',STATUS='OLD')
      READ(8,*)nbods
      READ(8,*)nsteps,nout
      READ(8,*)dt
      CLOSE(8)

C Read the cartesian coordinage of Pal5 postion.
      OPEN(UNIT=8,FILE='Pal5_initial',STATUS='OLD')
      DO 10 i=1,nbods
         READ(8,*)t_age,x(i),y(i),z(i),vx(i),vy(i),vz(i)
 10   CONTINUE
      CLOSE(8)

      RETURN
      END

C***********************************************************************
      SUBROUTINE integrate
C***********************************************************************
C
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1,dt2

      t=zero
      dt1=dt * Myr2sec
      dt2=dt/2. * Myr2sec

      i=0
      CALL accel(i)
C          -----
C
      CALL out
C          ---
      CALL stepvel(dt2)
C          -------
      DO 10 i=1,nsteps
C
         IF(MOD(i,1000).EQ.0)WRITE(6,*)i
C
         CALL steppos(dt1)
C             -------
         CALL accel(i)
C             -----
         IF(MOD(i,nout).EQ.0)THEN
            CALL stepvel(dt2)
C                -------
            CALL out
C               ----
c            IF(i+nout.GT.nsteps) CALL out_final
cC          --- Added by Joo H.
            CALL stepvel(dt2)
C                -------
         ELSE
            CALL stepvel(dt1)
C                -------
         ENDIF

 10   CONTINUE
      
      RETURN
      END

C***********************************************************************
      SUBROUTINE out
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e
      CHARACTER*10 nstring,sstring
      CHARACTER*30 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc

      IF(firstc)THEN
         firstc=.FALSE.
C
         DO 5 i=1,nbods
            sstring(1:4)='part'
            sstring(5:5)=nstring(1+i/1000:1+i/1000)
            istring=1+MOD(i,1000)/100
            sstring(6:6)=nstring(istring:istring)
            istring=1+MOD(i,100)/10
            sstring(7:7)=nstring(istring:istring)
            istring=1+MOD(i,10)
            sstring(8:8)=nstring(istring:istring)
C
            fname=sstring(1:8)
            OPEN(UNIT=10+i,FILE=fname,STATUS='UNKNOWN')

C The particles should conserve their total energy
            maxe(i)=zero
C Initial energy (JH)
            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
 5       CONTINUE
      ENDIF

      DO 10 i=1,nbods
C The particles should conserve their total energy         
C The energy at a certain time step which should be the same as the
C initial one(JH)
         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
C Fractional change
         e=(e-e0(i))/e0(i)
         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(e)
         WRITE(10+i,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),e
 10   CONTINUE

 99   FORMAT(8(1pe12.4))
      CLOSE(10+i)
      RETURN
      END

cC***********************************************************************
c      SUBROUTINE out_final
cC     output final position of particles   --- added by Joo H.
cC***********************************************************************
c      INCLUDE 'frog.h'
c      REAL*8 e
c
cC print out the initial position (at -2.95Gyr)
c      OPEN(UNIT=991,FILE='Pal5_initial',STATUS='UNKNOWN')
c      DO i=1,nbods
c         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
c         e=(e-e0(i))/e0(i)
c         WRITE(991,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),e
c      ENDDO
c 99   FORMAT(8(1pe12.4))
cC      CLOSE(999)
c      RETURN
c      END
    


C***********************************************************************
      SUBROUTINE steppos(dt1)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1

      DO 10 i=1,nbods
         x(i)=x(i)+dt1*vx(i)/kpc2km
         y(i)=y(i)+dt1*vy(i)/kpc2km
         z(i)=z(i)+dt1*vz(i)/kpc2km
 10   CONTINUE

C      t=t+dt1
      t=t+dt

      RETURN
      END

C***********************************************************************
      SUBROUTINE stepvel(dt1)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1

      DO 10 i=1,nbods
         vx(i)=vx(i)+dt1*ax(i)
         vy(i)=vy(i)+dt1*ay(i)
         vz(i)=vz(i)+dt1*az(i)
 10   CONTINUE

      RETURN
      END
