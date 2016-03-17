c234567
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n,it,j
      REAL*8 a,b,s,func,del,summ,tnm,x
      EXTERNAL func
      if (n.eq.1) then
         s=0.5*(b-a)*(func(a)+func(b))
      else
         it=2**(n-2)
         tnm=it
         del=(b-a)/tnm
         x=a+0.5*del
         summ=0.
         do j=1,it
            summ=summ+func(x)
            x=x+del
         enddo
         s=0.5*(s+(b-a)*summ/tnm)
      endif
      return
      end
