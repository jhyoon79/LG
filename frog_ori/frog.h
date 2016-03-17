      INTEGER nbods,nmax,nsteps,nout
      PARAMETER(nmax=1000)
      REAL*8 x,y,z,vx,vy,vz,ax,ay,az,pot,t,dt,maxe
      REAL*8 G,msun,zero
      PARAMETER(G=1.0d0,msun=1.0d0,zero=0.0d0) 
      COMMON/ints/nbods,nsteps,nout
      COMMON/time/t,dt
      COMMON/reals/x(0:nmax),y(0:nmax),z(0:nmax),
     &                vx(0:nmax),vy(0:nmax),vz(0:nmax),
     &                ax(0:nmax),ay(0:nmax),az(0:nmax),
     &                pot(0:nmax),maxe(0:nmax)


