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
C
C
      SUBROUTINE accel(j)
C
C
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,j
      REAL*8 r,ar
      
      DO 10 i=1,nbods
C Acceleration due to Sun
         r=SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         ar=-G*msun/(r*r)
C
         ax(i)=ar*x(i)/r
         ay(i)=ar*y(i)/r
         az(i)=ar*z(i)/r
C
         pot(i)=-G*msun/r
C
 10   CONTINUE
      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE check
C
C
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
C
C
      SUBROUTINE get_data
C
C
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i

      OPEN(UNIT=8,FILE='frogin',STATUS='OLD')
      READ(8,*)nbods
      READ(8,*)nsteps,nout
      READ(8,*)dt
      CLOSE(8)

C Read in test particle positions and velocities
      OPEN(UNIT=8,FILE='frog.dat',STATUS='OLD')
      DO 10 i=1,nbods
         READ(8,*)x(i),y(i),z(i),vx(i),vy(i),vz(i)
 10   CONTINUE
      CLOSE(8)

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE integrate
C
C
C***********************************************************************
C
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1,dt2

      t=zero
      dt1=dt
      dt2=dt/2.

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
C
C
      SUBROUTINE out
C
C
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i,istring
      REAL*8 e,e0(0:nmax)
      CHARACTER*10 nstring,sstring
      CHARACTER*7 fname
      LOGICAL firstc
      DATA nstring/'0123456789'/
      DATA firstc/.TRUE./
      SAVE nstring,firstc,e0

      IF(firstc)THEN
         firstc=.FALSE.
C
         DO 5 i=1,nbods
            sstring(1:4)='part'
            sstring(5:5)=nstring(1+i/100:1+i/100)
            istring=1+MOD(i,100)/10
            sstring(6:6)=nstring(istring:istring)
            istring=1+MOD(i,10)
            sstring(7:7)=nstring(istring:istring)
C
            fname=sstring(1:7)
            OPEN(UNIT=10+i,FILE=fname,STATUS='UNKNOWN')

C The particles should conserve their total energy
            maxe(i)=zero
            e0(i)=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
 5       CONTINUE
      ENDIF

      DO 10 i=1,nbods
C The particles should conserve their total energy         
         e=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.+pot(i)
C Fractional change
         e=(e-e0(i))/e0(i)
         IF(ABS(e).GT.maxe(i))maxe(i)=ABS(e)
         WRITE(10+i,99)t,x(i),y(i),z(i),vx(i),vy(i),vz(i),e
 10   CONTINUE

 99   FORMAT(8(1pe12.4))

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE steppos(dt1)
C
C
C***********************************************************************
      INCLUDE 'frog.h'
      INTEGER i
      REAL*8 dt1

      DO 10 i=1,nbods
         x(i)=x(i)+dt1*vx(i)
         y(i)=y(i)+dt1*vy(i)
         z(i)=z(i)+dt1*vz(i)
 10   CONTINUE

      t=t+dt1

      RETURN
      END

C***********************************************************************
C
C
      SUBROUTINE stepvel(dt1)
C
C
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





