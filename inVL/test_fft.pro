pro test_fft

N = 100
x = findgen(N+1)/N*!dpi*2 - !dpi*1.
T = 2*!dpi
y1 = sin(1*x)
y2 = sin(3*x)
y = y1+y2
M = indgen(N)-(N/2-1d)
dx = x[2]-x[1]
F = M/(N*dx) / (1d/T) 
;y = replicate(1,101)
;y[45:55] = 1
;y[55] = -30
;y[45] = -30
;y = [y0,y1,y0]
c = complex(x,y)
a = fft(y)
a = shift(a,N/2-1)
rea = real_part(a)
ima = imaginary(a)
!p.multi=[0,2,2]
plot,x,y
oplot,x,y1,linestyle=1
oplot,x,y2,linestyle=2
plot,F,rea,xr=[-5,5]
plot,F,ima,xr=[-5,5]
plot,F,abs(a),xr=[-5,5]

stop

END
