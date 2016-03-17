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
      
C test particles orbit
      DO 10 i=1,nbods_new
C NFW potential
         r = SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         p = r/rs
         mr = Mhalo*(dlog(p+1.)-p/(p+1.))
C truncation of NFW halo
         if (r.gt.rvir) mr = Mhalo*(dlog(c+1.)-c/(c+1.))
         aror = -(G*mr/r**2.)
         pot(i) = -(phi0/p)*dlog(p+1.)
C Acceleration by subhalos
         CALL accel_by_subhalo(i)
C Acceleration by hosthalo+subhalo
         aror_km_s2 = aror/kpc2km
c         ax(i)=aror_km_s2*x(i)/r
c         ay(i)=aror_km_s2*y(i)/r
c         az(i)=aror_km_s2*z(i)/r
         ax(i)=aror_km_s2*x(i)/r + ax_by_subhalo(i)
         ay(i)=aror_km_s2*y(i)/r + ay_by_subhalo(i)
         az(i)=aror_km_s2*z(i)/r + az_by_subhalo(i)
c        print *,'ax MW',i,r,ax_by_subhalo(i),ax(i)
 10   CONTINUE

C subhalos orbit
      DO 20 i=1,nsubhalos
C NFW potential
         r_subhalo=SQRT(x_subhalo(i)**2.+y_subhalo(i)**2.
     &             +z_subhalo(i)**2.)
         p_subhalo=r_subhalo/rs
         mr_subhalo=Mhalo*(dlog(p_subhalo+1.)-p_subhalo/(p_subhalo+1.))
C truncation of NFW halo
         if (r.gt.rvir) mr_subhalo = Mhalo*(dlog(c+1.)-c/(c+1.))
         aror_subhalo=-(G*mr_subhalo/r_subhalo**2.)
         pot_subhalo(i)=-(phi0/p_subhalo)*dlog(p_subhalo+1.)
C
         aror_subhalo_km_s2 = aror_subhalo/kpc2km
         ax_subhalo(i)=aror_subhalo_km_s2*x_subhalo(i)/r_subhalo
         ay_subhalo(i)=aror_subhalo_km_s2*y_subhalo(i)/r_subhalo
         az_subhalo(i)=aror_subhalo_km_s2*z_subhalo(i)/r_subhalo
 20   CONTINUE

      RETURN
      END


C***********************************************************************
      SUBROUTINE accel_by_subhalo(j)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,j
      REAL*8 c_subhalo, Mhalo_subhalo, phi0_subhalo
      REAL*8 r_subhalo_sat,p_subhalo_sat,mr_subhalo_sat,
     &      aror_subhalo_sat,aror_subhalo_sat_km_s2

      ax_by_subhalo(j)=0.
      ay_by_subhalo(j)=0.
      az_by_subhalo(j)=0.

      DO i=1,nsubhalos
C acceleration by subhalos
         c_subhalo = rvir_subhalo(i)/rs_subhalo(i)
         Mhalo_subhalo = Mvir_subhalo(i) / 
     &        (dlog(c_subhalo+1.)-c_subhalo/(c_subhalo+1.))
         phi0_subhalo = G*Mhalo_subhalo/rs_subhalo(i)
         r_subhalo_sat = SQRT((x(j)-x_subhalo(i))**2.
     &                       +(y(j)-y_subhalo(i))**2.
     &                       +(z(j)-z_subhalo(i))**2.)

         p_subhalo_sat=r_subhalo_sat/rs_subhalo(i)
         mr_subhalo_sat=Mhalo_subhalo*(dlog(p_subhalo_sat+1.)
     &                 -p_subhalo_sat/(p_subhalo_sat+1.))
C truncation of NFW subhalo
c         IF(r_subhalo_sat.gt.rvir_subhalo(i))
c     &               mr_subhalo_sat=Mvir_subhalo(i)
         IF(r_subhalo_sat.gt.rvir_subhalo(i)) mr_subhalo_sat = 
     &       Mhalo_subhalo*(dlog(c_subhalo+1.)-c_subhalo/(c_subhalo+1.))
         aror_subhalo_sat = -(G*mr_subhalo_sat/r_subhalo_sat**2.)
         aror_subhalo_sat_km_s2 = aror_subhalo_sat/kpc2km
         ax_each=aror_subhalo_sat_km_s2*(x(j)-x_subhalo(i))
     &           /r_subhalo_sat
         ay_each=aror_subhalo_sat_km_s2*(y(j)-y_subhalo(i))
     &           /r_subhalo_sat
         az_each=aror_subhalo_sat_km_s2*(z(j)-z_subhalo(i))
     &           /r_subhalo_sat
         ax_by_subhalo(j)=ax_by_subhalo(j)+ax_each
         ay_by_subhalo(j)=ay_by_subhalo(j)+ay_each
         az_by_subhalo(j)=az_by_subhalo(j)+az_each
c        print *,'ax by sub',i,r_subhalo_sat,ax_each,ax_by_subhalo(j)
      ENDDO
      RETURN
      END

C***********************************************************************
      SUBROUTINE check
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 eps
      PARAMETER(eps=0.01d0)
      
      DO 10 i=1,nbods_new
         IF(maxe(i).GT.eps)WRITE(6,*)i,maxe(i)
 10   CONTINUE

      RETURN
      END


C***********************************************************************
      SUBROUTINE get_data
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring,fstep
      REAL*8 m,tmp_x,tmp_y,tmp_z,tmp_vx,tmp_vy,tmp_vz,ru,mu,vu,tu
      REAL*8 tmp_pot_int,tmp_pot_ext,tunbound
      CHARACTER*10 tmp1
      CHARACTER*20 nstring,sstring,fname,fname_next
      DATA nstring/'0123456789'/
      LOGICAL firstc
      DATA firstc/.TRUE./
      DATA fstep/0/
      SAVE fstep,firstc,ru,mu,vu,tu

      IF(firstc)THEN
         OPEN(UNIT=8,FILE='frogin',STATUS='OLD')
         READ(8,*)nbods
         READ(8,*)nsteps,nout,nout_snap
         READ(8,*)dt
         READ(8,*)nsubhalos
         CLOSE(8)

C Read Nbody SCFPAR
         OPEN(UNIT=7,FILE=Nbody_dir//'SCFPAR',STATUS='OLD')
         DO 10 i=1,16
            READ(7,*)tmp1
 10      CONTINUE
         READ(7,*)ru
         READ(7,*)mu
         CLOSE(7)
         tu = sqrt((ru*kpc2km)**3./(mu*G*kpc2km))/(Myr2sec*1e3)
         vu = (ru*kpc2km)/tu/(Myr2sec*1e3)
         DO i=1,nbods
            pre_tunbound(i)=0
         ENDDO
         nbods_new = 0
         fstep = 1
      ENDIF

C Read in test particle positions and velocities from Nbody
      sstring(1:4)='SNAP'
      istring=1+MOD(fstep,1000)/100
      sstring(5:5)=nstring(istring:istring)
      istring=1+MOD(fstep,100)/10
      sstring(6:6)=nstring(istring:istring)
      istring=1+MOD(fstep,10)/1
      sstring(7:7)=nstring(istring:istring)
      fname=sstring(1:7)
      OPEN(UNIT=8,FILE=Nbody_dir//fname,STATUS='OLD')
      READ(8,*)tmp1,t_nbody
      t_nbody = t_nbody*tu*1000.
        print *,fname,t_nbody

      DO 15 i=1,nbods
         READ(8,*)m,tmp_x,tmp_y,tmp_z,tmp_vx,tmp_vy,tmp_vz,
     +       tmp_pot_int,tmp_pot_ext,tunbound
         if(tunbound.ne.0.and.pre_tunbound(i).eq.0)then
           nbods_new = nbods_new + 1
           x(nbods_new) = tmp_x*ru
           y(nbods_new) = tmp_y*ru
           z(nbods_new) = tmp_z*ru
           vx(nbods_new) = tmp_vx*vu
           vy(nbods_new) = tmp_vy*vu
           vz(nbods_new) = tmp_vz*vu
           pot_int(nbods_new) = tmp_pot_int*vu**2.
           pot_ext(nbods_new) = tmp_pot_ext*vu**2.
         endif
         pre_tunbound(i) = tunbound
 15   CONTINUE
        print *,'nbods_new',nbods_new
      CLOSE(8)

C read the time of next file
      fstep=fstep+1 
      sstring(1:4)='SNAP'
      istring=1+MOD(fstep,1000)/100
      sstring(5:5)=nstring(istring:istring)
      istring=1+MOD(fstep,100)/10
      sstring(6:6)=nstring(istring:istring)
      istring=1+MOD(fstep,10)/1
      sstring(7:7)=nstring(istring:istring)
      fname_next=sstring(1:7)
      OPEN(UNIT=3,FILE=Nbody_dir//fname_next,STATUS='OLD')
      READ(3,*)tmp1,t_nbody_next
      t_nbody_next = t_nbody_next*tu*1000.
        print *,fname,t_nbody_next
      CLOSE(3)

C Read in subhalo positions, velocities, and properties
      IF(firstc)THEN
         firstc=.FALSE.
c         OPEN(UNIT=9,FILE='../frog_VLsubhalo.dat',STATUS='OLD')
         OPEN(UNIT=9,
     +        FILE='../close_encounter/frog_5Rs_encounter_e7.dat')
         DO 20 i=1,nsubhalos
            READ(9,*)x_subhalo(i),y_subhalo(i),z_subhalo(i),
     &            vx_subhalo(i),vy_subhalo(i),vz_subhalo(i),
     &            Mvir_subhalo(i),rs_subhalo(i),rvir_subhalo(i)
 20      CONTINUE
         CLOSE(9)
      ENDIF

      RETURN
      END

C***********************************************************************
      SUBROUTINE integrate
C***********************************************************************
C
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1,dt2

      t=zero+t_nbody
      dt1=dt * Myr2sec
      dt2=dt/2. * Myr2sec
      i=0
      CALL accel(i)
C          -----
      CALL out_snap
      CALL out_snap_subhalo
      CALL out
C          ---
      CALL stepvel(dt2)
C          -------
      DO 10 i=1,nsteps
C
         IF(MOD(i,1000).EQ.0)WRITE(6,*)i,' steps'
C   Read out more data from Nbody result
            print *,'t ',t
         if(t.gt.t_nbody_next) CALL get_data
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
C               ----
c            out_snap_dt = dt*i
c            IF(MOD(out_snap_dt,nout_snap).EQ.0) THEN
            out_snap_dt = dt*i
            IF(MOD(i,nout_snap).EQ.0) THEN
               CALL out_snap
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

      DO j=1,nbods_new
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

      OPEN(UNIT=997,FILE=out_dir//'part_peri',STATUS='UNKNOWN')
      DO i=1,nbods_new
         WRITE(997,99)t_peri(i),x_peri(i),y_peri(i),z_peri(i),
     &       Vx_peri(i),Vy_peri(i),Vz_peri(i),e_peri(i),
     &       J_peri(i)
      ENDDO
 99   FORMAT(9(1pe14.6))
      CLOSE(997)

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
      CHARACTER*30 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc

      IF(firstc)THEN
         firstc=.FALSE.
         DO 5 i=1,nbods_new
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

      fname='snapshot/'//sstring(1:9)
      OPEN(UNIT=10,FILE=out_dir//fname,STATUS='UNKNOWN')

      DO 10 i=1,nbods_new
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
      SUBROUTINE out_snap_subhalo
C output snapshot of orbits
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e_subhalo,de_subhalo
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
      istring=1+MOD(out_snap_dt,100000)/10000
      sstring(8:8)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,10000)/1000
      sstring(9:9)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,1000)/100
      sstring(10:10)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,100)/10
      sstring(11:11)=nstring(istring:istring)
      istring=1+MOD(out_snap_dt,10)/1
      sstring(12:12)=nstring(istring:istring)

      fname='snapshot_subhalo/'//sstring(1:12)
      OPEN(UNIT=10,FILE=out_dir//fname,STATUS='UNKNOWN')

      DO 10 i=1,nsubhalos
         e_subhalo=(vx_subhalo(i)**2.+vy_subhalo(i)**2.
     &             +vz_subhalo(i)**2.)/2.+pot_subhalo(i)
c         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2 +
c     &           (z(i)*vx(i)-x(i)*vz(i))**2 +
c     &           (x(i)*vy(i)-y(i)*vx(i))**2 )
C Fractional change
         de_subhalo=(e_subhalo-e0_subhalo(i))/e0_subhalo(i)
         IF(ABS(e_subhalo).GT.maxe_subhalo(i))maxe_subhalo(i)
     &                                       =ABS(de_subhalo)
         WRITE(10,99)t,x_subhalo(i),y_subhalo(i),z_subhalo(i),
     &               vx_subhalo(i),vy_subhalo(i),vz_subhalo(i),
     &               de_subhalo
 10   CONTINUE

 99   FORMAT(8(1pe14.6))
      CLOSE(10)
      RETURN
      END


C***********************************************************************
      SUBROUTINE out
C   print out the position & velocity of the first particle(Pal5 center)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e,J,dJ,J0
c,e0(0:nmax)
      CHARACTER*10 nstring,sstring
      CHARACTER*7 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc,J0
c,e0

      IF(firstc)THEN
         firstc=.FALSE.
C
         DO 5 i=1,1
            sstring(1:4)='part'
            sstring(5:5)=nstring(1+i/100:1+i/100)
            istring=1+MOD(i,100)/10
            sstring(6:6)=nstring(istring:istring)
            istring=1+MOD(i,10)
            sstring(7:7)=nstring(istring:istring)
C
            fname=sstring(1:7)
            OPEN(UNIT=10+i,FILE=out_dir//fname,STATUS='UNKNOWN')

C The particles should conserve their total energy
            maxe(i)=zero
            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
            J0=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2. +
     &               (z(i)*vx(i)-x(i)*vz(i))**2. +
     &               (x(i)*vy(i)-y(i)*vx(i))**2. )
 5       CONTINUE
      ENDIF

      DO 10 i=1,1
C The particles should conserve their total energy         
         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
         J=SQRT( (y(i)*vz(i)-z(i)*vy(i))**2. +
     &           (z(i)*vx(i)-x(i)*vz(i))**2. +
     &           (x(i)*vy(i)-y(i)*vx(i))**2. )
         dJ = (J-J0)/J0
C Fractional change
         de=(e-e0(i))/e0(i)
         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(de)
         WRITE(10+i,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),de,e,dJ
 10   CONTINUE

 99   FORMAT(10(1pe14.6))

      RETURN
      END


C***********************************************************************
      SUBROUTINE steppos(dt1)
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1

      DO 10 i=1,nbods_new
         x(i)=x(i)+dt1*vx(i)/kpc2km
         y(i)=y(i)+dt1*vy(i)/kpc2km
         z(i)=z(i)+dt1*vz(i)/kpc2km
 10   CONTINUE

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

      DO 10 i=1,nbods_new
         vx(i)=vx(i)+dt1*ax(i)
         vy(i)=vy(i)+dt1*ay(i)
         vz(i)=vz(i)+dt1*az(i)
 10   CONTINUE

      DO 20 i=1,nsubhalos
         vx_subhalo(i)=vx_subhalo(i)+dt1*ax_subhalo(i)
         vy_subhalo(i)=vy_subhalo(i)+dt1*ay_subhalo(i)
         vz_subhalo(i)=vz_subhalo(i)+dt1*az_subhalo(i)
 20   CONTINUE

      RETURN
      END
