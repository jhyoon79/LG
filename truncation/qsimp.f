c234567
        SUBROUTINE qsimp(func,a,b,s)
        INTEGER JMAX,j
        REAL*8 a,b,func,s,EPS,os,ost,st
        EXTERNAL func
        PARAMETER (EPS=1.E-6,JMAX=20)
        ost=-1.e30
        os=-1.e30
        do j=1,JMAX
            call trapzd(func,a,b,st,j)
            s=(4.*st-ost)/3.
            if (j.gt.5) then
                if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.))
     +              return
            endif
            os=s
            ost=st
        enddo
        pause 'too many steps in qsimp'
        end
