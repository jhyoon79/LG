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

C***********************************************************************
      SUBROUTINE accel(j)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,j
      REAL*8 r,aror,c,Mhalo,phi0,p,mr,aror_km_s2
      SAVE c, Mhalo, phi0
      c = rvir/rs
      Mhalo = Mvir/(dlog(c+1.)-c/(c+1.))
      phi0 = G*Mhalo/rs
      
C test particles orbit
      DO 10 i=1,nbods
C NFW potential
         r = SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         p = r/rs
         mr = Mhalo*(dlog(p+1.)-p/(p+1.))
C truncation of NFW halo
         if (r.gt.rvir) mr = Mhalo*(dlog(c+1.)-c/(c+1.))
         aror = -(G*mr/r**2.)
         pot(i) = -(phi0/p)*dlog(p+1.)
         aror_km_s2 = aror/kpc2km
         ax(i)=aror_km_s2*x(i)/r
         ay(i)=aror_km_s2*y(i)/r
         az(i)=aror_km_s2*z(i)/r
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

      OPEN(UNIT=8,FILE='frogin_backward_cen',STATUS='OLD')
      READ(8,*)nbods
      READ(8,*)nsteps,nout,nout_snap
      READ(8,*)dt
      CLOSE(8)

C Read in test particle positions and velocities
      OPEN(UNIT=8,FILE='nbody_SCFCEN12796',STATUS='OLD')
      DO 10 i=1,nbods
         READ(8,*)tmp,x(i),y(i),z(i),vx(i),vy(i),vz(i)
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
      CALL out_snap
C          ---
      CALL stepvel(dt2)
C          -------
      DO 10 i=1,nsteps
C
         IF(MOD(i,1000).EQ.0)WRITE(6,*)i,' steps'
C   Integrate the orbit
         CALL steppos(dt1)
C             -------
         CALL accel(i)
C             -----
         IF(MOD(i,nout).EQ.0)THEN
            CALL stepvel(dt2)
C                -------
C               ----
c            out_snap_dt = dt*i
c            IF(MOD(out_snap_dt,nout_snap).EQ.0) THEN
            out_snap_dt = 12796+dt*i
            IF(MOD(i,nout_snap).EQ.0) THEN
               CALL out_snap
            ENDIF
C   output snapshot in every 10Myr
c            IF(i+nout.GT.nsteps) CALL out_final
C          --- Added by Joo H.
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
      SUBROUTINE out_snap
C output snapshot of orbits
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e,J,de
      CHARACTER*20 nstring,sstring
      CHARACTER*40 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc

      IF(firstc)THEN
         firstc=.FALSE.
         DO 5 i=1,nbods
C The particles should conserve their total energy
            maxe(i)=zero
C Initial energy --- added by J.H. Yoon
            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
 5       CONTINUE
      ENDIF

      sstring(1:4)='snap'
      istring=1+MOD(out_snap_dt,100000)/10000
      sstring(5:5)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,10000)/1000
      sstring(6:6)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,1000)/100
      sstring(7:7)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,100)/10
      sstring(8:8)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,10)/1
      sstring(9:9)=nstring(istring:istring)

      fname='backward/SCFCEN_'//sstring(1:9)
      OPEN(UNIT=10,FILE=Nbody_dir//fname,STATUS='UNKNOWN')

      DO 10 i=1,nbods
C The particles should conserve their total energy         
C The energy at a certain time step which should be the same as the
C initial one --- added by J.H. Yoon
         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2 +
     &           (z(i)*vx(i)-x(i)*vz(i))**2 +
     &           (x(i)*vy(i)-y(i)*vx(i))**2 )
C Fractional change
         de=(e-e0(i))/e0(i)
         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(de)
         WRITE(10,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),de,e,J
 10   CONTINUE

 99   FORMAT(10(1pe14.6))
      CLOSE(10)
      RETURN
      END


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
