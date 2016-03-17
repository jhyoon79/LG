C================================================
C This is to find the orbits of all Via Lactea subhalos
C================================================

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
      REAL*8 r_subhalo,aror_subhalo,p_subhalo,mr_subhalo
      SAVE c, Mhalo, phi0
      c = rvir/rs
      Mhalo = Mvir/(dlog(c+1.)-c/(c+1.))
      phi0 = G*Mhalo/rs
      
C subhalos orbit
      DO 20 i=1,nsubhalos
C NFW potential
         r_subhalo=SQRT(x_subhalo(i)**2.+y_subhalo(i)**2.
     &             +z_subhalo(i)**2.)
         p_subhalo=r_subhalo/rs
         mr_subhalo=Mhalo*(dlog(p_subhalo+1.)-p_subhalo/(p_subhalo+1.))
         aror_subhalo=-(G*mr_subhalo/r_subhalo**2.)
         pot_subhalo(i)=-(phi0/p_subhalo)*dlog(p_subhalo+1.)
C
         aror_subhalo_km_s2 = aror_subhalo/kpc2km
         ax_subhalo(i)=aror_subhalo_km_s2*x_subhalo(i)/r_subhalo
         ay_subhalo(i)=aror_subhalo_km_s2*y_subhalo(i)/r_subhalo
         az_subhalo(i)=aror_subhalo_km_s2*z_subhalo(i)/r_subhalo
 20   CONTINUE

cC test particles orbit
c      DO 10 i=1,nbods
cC NFW potential
c         r = SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
c         p = r/rs
c         mr = Mhalo*(dlog(p+1.)-p/(p+1.))
c         aror = -(G*mr/r**2.)
c         pot(i) = -(phi0/p)*dlog(p+1.)
cC Acceleration by subhalos
c         CALL accel_by_subhalo(i)
cC Acceleration by hosthalo+subhalo
c         aror_km_s2 = aror/kpc2km
cc        print *,'host',Mvir,Mhalo,mr,p
cc        print *, 'pre',i,aror_km_s2*x(i)/r,ax_by_subhalo(i)
cc         ax(i)=aror_km_s2*x(i)/r
cc         ay(i)=aror_km_s2*y(i)/r
cc         az(i)=aror_km_s2*z(i)/r
c         ax(i)=aror_km_s2*x(i)/r + ax_by_subhalo(i)
c         ay(i)=aror_km_s2*y(i)/r + ay_by_subhalo(i)
c         az(i)=aror_km_s2*z(i)/r + az_by_subhalo(i)
cc        print *, 'aft',i,ax(i)
c 10   CONTINUE
      RETURN
      END


cC***********************************************************************
c      SUBROUTINE accel_by_subhalo(j)
cC***********************************************************************
c      INCLUDE 'frog.h'
c      INTEGER i,j
c      REAL*8 c_subhalo, Mhalo_subhalo, phi0_subhalo
c      REAL*8 r_subhalo_sat,p_subhalo_sat,mr_subhalo_sat,
c     &      aror_subhalo_sat,aror_subhalo_sat_km_s2
c
c      DO i=1,nsubhalos
cC acceleration by subhalos
c         c_subhalo = rvir_subhalo(i)/rs_subhalo(i)
c         Mhalo_subhalo = Mvir_subhalo(i) / 
c     &        (dlog(c_subhalo+1.)-c_subhalo/(c_subhalo+1.))
c         phi0_subhalo = G*Mhalo_subhalo/rs_subhalo(i)
c         r_subhalo_sat = SQRT((x(j)-x_subhalo(i))**2.
c     &                       +(y(j)-y_subhalo(i))**2.
c     &                       +(z(j)-z_subhalo(i))**2.)
c         p_subhalo_sat=r_subhalo_sat/rs_subhalo(i)
c         mr_subhalo_sat=Mhalo_subhalo*(dlog(p_subhalo_sat+1.)
c     &                 -p_subhalo_sat/(p_subhalo_sat+1.))
c         IF(r_subhalo_sat.ge.rvir_subhalo(i))
c     &               mr_subhalo_sat=Mvir_subhalo(i)
c         aror_subhalo_sat = -(G*mr_subhalo_sat/r_subhalo_sat**2.)
c         aror_subhalo_sat_km_s2 = aror_subhalo_sat/kpc2km
cc      print *,'mass',Mvir_subhalo(i),Mhalo_subhalo,mr_subhalo_sat
cc      print *,'r_sub_sat',p_subhalo_sat,r_subhalo_sat
cc      print *,'a',aror_subhalo_sat_km_s2
c         ax_each=aror_subhalo_sat_km_s2*(x(j)-x_subhalo(i))
c     &           /r_subhalo_sat
c         ay_each=aror_subhalo_sat_km_s2*(y(j)-y_subhalo(i))
c     &           /r_subhalo_sat
c         az_each=aror_subhalo_sat_km_s2*(z(j)-z_subhalo(i))
c     &           /r_subhalo_sat
c         ax_by_subhalo(j)=ax_by_subhalo(j)+ax_each
c         ay_by_subhalo(j)=ay_by_subhalo(j)+ay_each
c         az_by_subhalo(j)=az_by_subhalo(j)+az_each
cc        print *,'accel_sub',ax_each, ax_by_subhalo(j)
c      ENDDO
c      RETURN
c      END

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

      OPEN(UNIT=8,FILE='frogin',STATUS='OLD')
      READ(8,*)nbods
      READ(8,*)nsteps,nout,nout_snap
      READ(8,*)dt
      READ(8,*)nsubhalos
      CLOSE(8)

C Read in test particle positions and velocities
c      OPEN(UNIT=8,FILE='frog_Pal5.dat',STATUS='OLD')
c      DO 10 i=1,nbods
c         READ(8,*)x(i),y(i),z(i),vx(i),vy(i),vz(i)
c 10   CONTINUE
c      CLOSE(8)

C Read in subhalo positions, velocities, and properties
      OPEN(UNIT=9,FILE='frog_VLsubhalo.dat',STATUS='OLD')
      DO 20 i=1,nsubhalos
         READ(9,*) x_subhalo(i),y_subhalo(i),z_subhalo(i),
     &         vx_subhalo(i),vy_subhalo(i),vz_subhalo(i)
 20   CONTINUE
      CLOSE(9)

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
c      CALL out_snap
      CALL out_snap_subhalo
c      CALL out
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
c            CALL out
C               ----
            out_snap_dt = dt*i
            IF(MOD(out_snap_dt,nout_snap).EQ.0) THEN
c              CALL out_snap
               CALL out_snap_subhalo
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
c         CALL find_peri(i)
 10   CONTINUE

c      CALL out_peri
      
      RETURN
      END


cC***********************************************************************
c      SUBROUTINE find_peri(i)
cC   Find the pericenter
cC***********************************************************************
c      INCLUDE 'frog.h'
c      INTEGER j
c      REAL*8 r
c
c      DO j=1,nbods 
c         IF(i.eq.1) r_peri(j)=99999.
c         r = sqrt(x(j)**2.+y(j)**2.+z(j)**2.)
c         IF(r.LE.r_peri(j))THEN 
c            r_peri(j) = r
c            t_peri(j) = t
c            x_peri(j) = x(j)
c            y_peri(j) = y(j)
c            z_peri(j) = z(j)
c            Vx_peri(j) = Vx(j)
c            Vy_peri(j) = Vy(j)
c            Vz_peri(j) = Vz(j)
c            pot_peri(j) = pot(j)
c            e_peri(j)=(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))/2.+pot(j)
c            J_peri(j)=SQRT( (y(j)*vz(j)-z(j)*vy(j))**2 +
c     &                      (z(j)*vx(j)-x(j)*vz(j))**2 +
c     &                      (x(j)*vy(j)-y(j)*vx(j))**2 )
c            Jz_peri(j)=x(j)*vy(j)-y(j)*vx(j)
c         ENDIF
c      ENDDO
c
c      RETURN
c      END
c
c
cC***********************************************************************
c      SUBROUTINE out_peri
cC     output the position of pericenter  --- added by JH Yoon
cC***********************************************************************
c      INCLUDE 'frog.h'
c      REAL*8 e
c      INTEGER i
c
c      OPEN(UNIT=997,FILE='part_peri',STATUS='UNKNOWN')
c      DO i=1,nbods
c         WRITE(997,99)t_peri(i),x_peri(i),y_peri(i),z_peri(i),
c     &       Vx_peri(i),Vy_peri(i),Vz_peri(i),e_peri(i),
c     &       J_peri(i),Jz_peri(i)
c      ENDDO
c 99   FORMAT(10(1pe12.4))
c      CLOSE(997)
c
c      RETURN
c      END


cC***********************************************************************
c      SUBROUTINE out_snap
cC output snapshot of orbits
cC***********************************************************************
c      INCLUDE 'frog.h'
c      INTEGER i,istring
c      REAL*8 e,J,Jz,de
cC   , E_total
c      CHARACTER*10 nstring,sstring
c      CHARACTER*30 fname
c      LOGICAL firstc
c      DATA nstring/'0123456789'/
c      DATA firstc/.TRUE./
c      SAVE nstring,firstc
c
c      IF(firstc)THEN
c         firstc=.FALSE.
c         DO 5 i=1,nbods
cC The particles should conserve their total energy
c            maxe(i)=zero
cC Initial energy --- added by J.H. Yoon
c            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
c 5       CONTINUE
c      ENDIF
c
c      sstring(1:4)='snap'
c      sstring(5:5)=nstring(1+out_snap_dt/1000:1+out_snap_dt/1000)
c      istring=1+MOD(out_snap_dt,1000)/100
c      sstring(6:6)=nstring(istring:istring)
c      istring=1+MOD(out_snap_dt,100)/10
c      sstring(7:7)=nstring(istring:istring)
c      sstring(8:8)='0'
c
c      fname='snapshot/'//sstring(1:8)
c      OPEN(UNIT=10,FILE=fname,STATUS='UNKNOWN')
c
c      DO 10 i=1,nbods
cC The particles should conserve their total energy         
cC The energy at a certain time step which should be the same as the
cC initial one --- added by J.H. Yoon
c         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
c         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2 +
c     &           (z(i)*vx(i)-x(i)*vz(i))**2 +
c     &           (x(i)*vy(i)-y(i)*vx(i))**2 )
c         Jz=x(i)*vy(i)-y(i)*vx(i)
cC Fractional change
c         E_total(i) = e
c         de=(e-e0(i))/e0(i)
c         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(de)
c         WRITE(10,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),
c     &               de,E_total(i),J,Jz
c 10   CONTINUE
c
c 99   FORMAT(11(1pe12.4))
c      CLOSE(10)
c      RETURN
c      END


C***********************************************************************
      SUBROUTINE out_snap_subhalo
C output snapshot of orbits
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e_subhalo,de_subhalo
C   , E_total
      CHARACTER*20 nstring,sstring
      CHARACTER*40 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc

      IF(firstc)THEN
         firstc=.FALSE.
         DO 5 i=1,nsubhalos
C The particles should conserve their total energy
            maxe_subhalo(i)=zero
C Initial energy --- added by J.H. Yoon
            e0_subhalo(i)=(vx_subhalo(i)**2.+vy_subhalo(i)**2.
     &                    +vz_subhalo(i)**2.)/2.+pot_subhalo(i)
 5       CONTINUE
      ENDIF

      sstring(1:7)='subhalo'
      sstring(8:8)=nstring(1+out_snap_dt/1000:1+out_snap_dt/1000)
      istring=1+MOD(out_snap_dt,1000)/100
      sstring(9:9)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,100)/10
      sstring(10:10)=nstring(istring:istring)
      sstring(11:11)='0'

      fname='snapshot_VLsubhalo/'//sstring(1:11)
      OPEN(UNIT=10,FILE=fname,STATUS='UNKNOWN')

      DO 10 i=1,nsubhalos
         e_subhalo=(vx_subhalo(i)**2.+vy_subhalo(i)**2.
     &             +vz_subhalo(i)**2.)/2.+pot_subhalo(i)
c         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2 +
c     &           (z(i)*vx(i)-x(i)*vz(i))**2 +
c     &           (x(i)*vy(i)-y(i)*vx(i))**2 )
c         Jz=x(i)*vy(i)-y(i)*vx(i)
C Fractional change
c         E_total(i) = e
         de_subhalo=(e_subhalo-e0_subhalo(i))/e0_subhalo(i)
         IF(ABS(e_subhalo).GT.maxe_subhalo(i))maxe_subhalo(i)
     &                                       =ABS(de_subhalo)
         WRITE(10,99)t,x_subhalo(i),y_subhalo(i),z_subhalo(i),
     &               vx_subhalo(i),vy_subhalo(i),vz_subhalo(i),
     &               de_subhalo
 10   CONTINUE

 99   FORMAT(8(1pe12.4))
      CLOSE(10)
      RETURN
      END


cC***********************************************************************
c      SUBROUTINE out
cC   print out the position & velocity of the first particle(Pal5 center)
cC***********************************************************************
c      INCLUDE 'frog.h'
c      INTEGER i,istring
c      REAL*8 e,J,dJ,Jz,dJz,J0,Jz0
cc,e0(0:nmax)
c      CHARACTER*10 nstring,sstring
c      CHARACTER*7 fname
c      LOGICAL firstc
c      DATA nstring/'0123456789'/
c      DATA firstc/.TRUE./
c      SAVE nstring,firstc,J0,Jz0
cc,e0
c
c      IF(firstc)THEN
c         firstc=.FALSE.
cC
cc         DO 5 i=1,nbods
c         DO 5 i=1,1
c            sstring(1:4)='part'
c            sstring(5:5)=nstring(1+i/100:1+i/100)
c            istring=1+MOD(i,100)/10
c            sstring(6:6)=nstring(istring:istring)
c            istring=1+MOD(i,10)
c            sstring(7:7)=nstring(istring:istring)
cC
c            fname=sstring(1:7)
c            OPEN(UNIT=10+i,FILE=fname,STATUS='UNKNOWN')
c
cC The particles should conserve their total energy
c            maxe(i)=zero
c            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
c            J0=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2. +
c     &               (z(i)*vx(i)-x(i)*vz(i))**2. +
c     &               (x(i)*vy(i)-y(i)*vx(i))**2. )
c            Jz0=x(i)*vy(i)-y(i)*vx(i)
c 5       CONTINUE
c      ENDIF
c
cc      DO 10 i=1,nbods
c      DO 10 i=1,1
cC The particles should conserve their total energy         
c         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
c         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2. +
c     &           (z(i)*vx(i)-x(i)*vz(i))**2. +
c     &           (x(i)*vy(i)-y(i)*vx(i))**2. )
c         Jz=x(i)*vy(i)-y(i)*vx(i)
c         dJ = (J-J0)/J0
c         dJz = (Jz-Jz0)/Jz0
cC Fractional change
c         de=(e-e0(i))/e0(i)
c         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(de)
c         WRITE(10+i,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),de,e,dJ,dJz
c 10   CONTINUE
c
c 99   FORMAT(11(1pe12.4))
c
c      RETURN
c      END


C***********************************************************************
      SUBROUTINE steppos(dt1)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1

c      DO 10 i=1,nbods
c         x(i)=x(i)+dt1*vx(i)/kpc2km
c         y(i)=y(i)+dt1*vy(i)/kpc2km
c         z(i)=z(i)+dt1*vz(i)/kpc2km
c 10   CONTINUE

      DO 20 i=1,nsubhalos
         x_subhalo(i)=x_subhalo(i)+dt1*vx_subhalo(i)/kpc2km
         y_subhalo(i)=y_subhalo(i)+dt1*vy_subhalo(i)/kpc2km
         z_subhalo(i)=z_subhalo(i)+dt1*vz_subhalo(i)/kpc2km
 20   CONTINUE

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

c      DO 10 i=1,nbods
c         vx(i)=vx(i)+dt1*ax(i)
c         vy(i)=vy(i)+dt1*ay(i)
c         vz(i)=vz(i)+dt1*az(i)
c 10   CONTINUE

      DO 20 i=1,nsubhalos
         vx_subhalo(i)=vx_subhalo(i)+dt1*ax_subhalo(i)
         vy_subhalo(i)=vy_subhalo(i)+dt1*ay_subhalo(i)
         vz_subhalo(i)=vz_subhalo(i)+dt1*az_subhalo(i)
 20   CONTINUE

      RETURN
      END
