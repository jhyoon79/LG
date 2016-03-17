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
         aror = -(G*mr/r**2.)
         pot(i) = -(phi0/p)*dlog(p+1.)
         aror_km_s2 = aror/kpc2km
         ax(i)=aror_km_s2*x(i)/r
         ay(i)=aror_km_s2*y(i)/r
         az(i)=aror_km_s2*z(i)/r
c        print *, 'aft',i,ax(i)
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

      OPEN(UNIT=8,FILE='frogin_satellite',STATUS='OLD')
      READ(8,*)nbods
      READ(8,*)nsteps,nout,nout_snap
      READ(8,*)dt
      CLOSE(8)

C Read in test particle positions and velocities
      OPEN(UNIT=8,FILE='frog_satellite.dat',STATUS='OLD')
      DO 10 i=1,nbods
         READ(8,*)x(i),y(i),z(i),vx(i),vy(i),vz(i)
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
      CALL out
c      CALL out_snap
C          ---
      CALL stepvel(dt2)
C          -------
      DO 10 i=1,nsteps
C
         IF(MOD(i,1000).EQ.0)WRITE(6,*)i,' steps'
C   Find the pericenter
C   Integrate the orbit
         CALL steppos(dt1)
C             -------
         CALL accel(i)
C             -----
         IF(MOD(i,nout).EQ.0)THEN
            CALL stepvel(dt2)
C                -------
            CALL out
            out_snap_dt = dt*i
c            IF(MOD(out_snap_dt,nout_snap).EQ.0) THEN
c               CALL out_snap
c            ENDIF
C   output snapshot in every 10Myr
c            IF(i+nout.GT.nsteps) CALL out_final
C          --- Added by Joo H.
            CALL stepvel(dt2)
C                -------
         ELSE
            CALL stepvel(dt1)
C                -------
         ENDIF
         CALL find_peri(i)
 10   CONTINUE

      CALL out_peri
      
      RETURN
      END


C***********************************************************************
      SUBROUTINE find_peri(i)
C   Find the pericenter
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER j
      REAL*8 r

      DO j=1,nbods 
         IF(i.eq.1) r_peri(j)=99999.
         r = sqrt(x(j)**2.+y(j)**2.+z(j)**2.)
         IF(r.LE.r_peri(j))THEN 
            r_peri(j) = r
            t_peri(j) = t
            x_peri(j) = x(j)
            y_peri(j) = y(j)
            z_peri(j) = z(j)
            Vx_peri(j) = Vx(j)
            Vy_peri(j) = Vy(j)
            Vz_peri(j) = Vz(j)
            pot_peri(j) = pot(j)
            e_peri(j)=(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))/2.+pot(j)
            J_peri(j)=SQRT( (y(j)*vz(j)-z(j)*vy(j))**2 +
     &                      (z(j)*vx(j)-x(j)*vz(j))**2 +
     &                      (x(j)*vy(j)-y(j)*vx(j))**2 )
            Jz_peri(j)=x(j)*vy(j)-y(j)*vx(j)
         ENDIF
      ENDDO

      RETURN
      END


C***********************************************************************
      SUBROUTINE out_peri
C     output the position of pericenter  --- added by JH Yoon
C***********************************************************************
      INCLUDE 'frog.h'
      REAL*8 e
      INTEGER i

      OPEN(UNIT=997,FILE='satellite_peri',STATUS='UNKNOWN')
      DO i=1,nbods
         WRITE(997,99)t_peri(i),x_peri(i),y_peri(i),z_peri(i),
     &       Vx_peri(i),Vy_peri(i),Vz_peri(i),e_peri(i),
     &       J_peri(i),Jz_peri(i)
      ENDDO
 99   FORMAT(10(1pe12.4))
      CLOSE(997)

      RETURN
      END


C***********************************************************************
      SUBROUTINE out_snap
C output snapshot of orbits
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e,J,Jz,de
C   , E_total
      CHARACTER*10 nstring,sstring
      CHARACTER*30 fname
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

      sstring(1:4)='sate'
      sstring(5:5)=nstring(1+out_snap_dt/1000:1+out_snap_dt/1000)
      istring=1+MOD(out_snap_dt,1000)/100
      sstring(6:6)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,100)/10
      sstring(7:7)=nstring(istring:istring)
      sstring(8:8)='0'

      fname='satellite/'//sstring(1:8)
      OPEN(UNIT=10,FILE=fname,STATUS='UNKNOWN')

      DO 10 i=1,nbods
C The particles should conserve their total energy         
C The energy at a certain time step which should be the same as the
C initial one --- added by J.H. Yoon
         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2 +
     &           (z(i)*vx(i)-x(i)*vz(i))**2 +
     &           (x(i)*vy(i)-y(i)*vx(i))**2 )
         Jz=x(i)*vy(i)-y(i)*vx(i)
C Fractional change
         E_total(i) = e
         de=(e-e0(i))/e0(i)
         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(de)
         WRITE(10,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),
     &               de,E_total(i),J,Jz
 10   CONTINUE

 99   FORMAT(11(1pe12.4))
      CLOSE(10)
      RETURN
      END


C***********************************************************************
      SUBROUTINE out
C   print out the position & velocity of the first particle(Pal5 center)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e,J,dJ,Jz,dJz,J0,Jz0
c,e0(0:nmax)
      CHARACTER*10 nstring,sstring
      CHARACTER*30 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc,J0,Jz0
c,e0

      IF(firstc)THEN
         firstc=.FALSE.
C
         DO 5 i=1,nbods
c         DO 5 i=1,1
            sstring(1:3)='sat'
            sstring(4:4)=nstring(1+i/100:1+i/100)
            istring=1+MOD(i,100)/10
            sstring(5:5)=nstring(istring:istring)
            istring=1+MOD(i,10)
            sstring(6:6)=nstring(istring:istring)
C
            fname=sstring(1:6)
            fname='satellite/'//sstring(1:6)
            OPEN(UNIT=10+i,FILE=fname,STATUS='UNKNOWN')

C The particles should conserve their total energy
            maxe(i)=zero
            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
            J0=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2. +
     &               (z(i)*vx(i)-x(i)*vz(i))**2. +
     &               (x(i)*vy(i)-y(i)*vx(i))**2. )
            Jz0=x(i)*vy(i)-y(i)*vx(i)
 5       CONTINUE
      ENDIF

      DO 10 i=1,nbods
c      DO 10 i=1,1
C The particles should conserve their total energy         
         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2. +
     &           (z(i)*vx(i)-x(i)*vz(i))**2. +
     &           (x(i)*vy(i)-y(i)*vx(i))**2. )
         Jz=x(i)*vy(i)-y(i)*vx(i)
         dJ = (J-J0)/J0
         dJz = (Jz-Jz0)/Jz0
C Fractional change
         de=(e-e0(i))/e0(i)
         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(de)
         WRITE(10+i,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),de,e,dJ,dJz
 10   CONTINUE

 99   FORMAT(11(1pe12.4))

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
