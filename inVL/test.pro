pro test

N = 3000
x=randomu(1,N)*20
Bk = dblarr(N)
k = (findgen(N)+1)/N*(360d)

for m=0,N-1 do begin
  fact = k[m]*x*!dtor
  Bk[m] = abs(1./N*total(complex(cos(fact),sin(fact))))
endfor
k = 360/k
@plot_setting
device,file='test.ps'
plot,k,Bk,/xlog,xr=[1,40]
vline,max(x)/(findgen(20)+1)
device,/close
stop
END
