      INTEGER nbods,nmax,nsteps,nout,nout_snap,nsubhalos
      PARAMETER(nmax=20000)
      REAL*8 x,y,z,vx,vy,vz,ax,ay,az,pot,t,dt,maxe,e0
      REAL*8 x_subhalo,y_subhalo,z_subhalo,vx_subhalo,vy_subhalo,
     &       vz_subhalo,ax_subhalo,ay_subhalo,az_subhalo,
     &       Mvir_subhalo,rs_subhalo,rvir_subhalo,pot_subhalo
     &       e0_subhalo,maxe_subhalo
      REAL*8 ax_by_subhalo,ay_by_subhalo,az_by_subhalo
      REAL*8 G,msun,zero,Mvir,rs,rvir,kpc2km,Myr2sec
      REAL*8 r_peri,t_peri,x_peri,y_peri,z_peri,Vx_peri,Vy_peri
      REAL*8 Vz_peri,pot_peri,e_peri,J_peri
      CHARACTER*33 out_dir
      DATA out_dir/'/scratch/jhyoon/Research/LG/inVL/'/
C     G = 4.3*10^-6 kpc M_solar^-1 (km/s)^2
      PARAMETER(G=4.3e-6,msun=1.0d0,zero=0.0d0) 
      PARAMETER(Mvir=1.77e+12,rs=24.6d0,rvir=389.0d0)
c      PARAMETER(Mvir=1.4e+12,rs=20.7406d0,rvir=287.7606d0)
      PARAMETER(kpc2km=3.0857e16,Myr2sec=3.1536e13)
      COMMON/ints/nbods,nsteps,nout,nout_snap,nsubhalos
      COMMON/time/t,dt,out_snap_dt
      COMMON/reals/x(0:nmax),y(0:nmax),z(0:nmax),
     &                vx(0:nmax),vy(0:nmax),vz(0:nmax),
     &                ax(0:nmax),ay(0:nmax),az(0:nmax),
     &                pot(0:nmax),maxe(0:nmax),e0(0:nmax)
      COMMON/subhalo/x_subhalo(0:nmax),y_subhalo(0:nmax),
     &       z_subhalo(0:nmax),vx_subhalo(0:nmax),vy_subhalo(0:nmax),
     &       vz_subhalo(0:nmax),ax_subhalo(0:nmax),ay_subhalo(0:nmax),
     &       az_subhalo(0:nmax),Mvir_subhalo(0:nmax),rs_subhalo(0:nmax),
     &       rvir_subhalo(0:nmax),pot_subhalo(0:nmax),
     &       e0_subhalo(0:nmax),maxe_subhalo(0:nmax),
     &       ax_by_subhalo(0:nmax),ay_by_subhalo(0:nmax),
     &       az_by_subhalo(0:nmax)
      COMMON/peri/r_peri(0:nmax),x_peri(0:nmax),y_peri(0:nmax),
     &                z_peri(0:nmax),e_peri(0:nmax),
     &                Vx_peri(0:nmax),Vy_peri(0:nmax),Vz_peri(0:nmax),
     &                t_peri(0:nmax),pot_peri(0:nmax),J_peri(0:nmax)


