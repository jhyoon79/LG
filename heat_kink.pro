pro heat_kink

sigmav=[0.,0.15,0.3,0.6,1.]
Pheat=[1.7888544,7.7653075,15.201974,38.288815,76.607441]
Pheat_radial=[2.9389316,19.303190,30.790077,48.891427,68.020794]
set_plot,'ps'
@plot_setting
loadct,13
device,file='heat_kink.ps',/color
plot,sigmav,Pheat,xr=[-0.02,1.02],yr=[0,100],xtitle=textoidl('\sigma_v [km/s]'),ytitle=textoidl('P_{heat}')
oplot,sigmav,Pheat,psym=1
oplot,sigmav,Pheat_radial,psym=1,color=255
oplot,sigmav,Pheat_radial,color=255


N = 3000
x = randomu(1,N)*10.
x1 = randomu(99,N/4)*10;randomn(99,N/4)+3.
x2 = randomn(98,N/4)/2.+5.
x3 = randomn(97,N/4)/3.+9.
x4 = x1;randomn(96,N/4)+7.
x = [x1,x2,x3,x4]
x2 = randomn(98,N/4)/3.+5.
x3 = randomn(97,N/4)/4.+9.
xa = [x1,x2,x3,x4]
x2 = randomn(98,N/4)/3.+5.
x3 = randomn(97,N/4)/4.1+9.
xb = [x1,x2,x3,x4]
x2 = randomn(98,N/4)/1.+5.
x3 = randomn(97,N/4)/2.+9.
xc = [x1,x2,x3,x4]
xc = randomu(iseed,N)*10

y = randomn(2,N)
z = randomn(3,N)
r = sqrt(x^2.+y^2.+z^2.)
sigma_v = 0.15
Vx = randomn(4,N)*sigma_v
Vy = randomn(5,N)*sigma_v
Vz = randomn(6,N)*sigma_v
Mvir=1.77d12
rvir = 389.d
rs = 24.6d
G = 4.3e-6
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
phi0 = G*Mhalo/rs
p = r/rs
pot = -phi0/p*alog(p+1.)
E = (Vx^2.+Vy^2.+Vz^2.)/2.+pot
sum_Ediff = 0.d
Ediff = dblarr(N)
min_d1 = dblarr(N)
min_d2 = dblarr(N)
min_d3 = dblarr(N)
min_d4 = dblarr(N)
min_d5 = dblarr(N)
min_d6 = dblarr(N)
min_d5a = dblarr(N)
min_d5b = dblarr(N)
min_d5c = dblarr(N)
for i=0,N-1 do begin
  d = sqrt((x[i]-x)^2.+(y[i]-y)^2.+(z[i]-z)^2.)
  da = sqrt((xa[i]-xa)^2.+(y[i]-y)^2.+(z[i]-z)^2.)
  db = sqrt((xb[i]-xb)^2.+(y[i]-y)^2.+(z[i]-z)^2.)
  dc = sqrt((xc[i]-xc)^2.+(y[i]-y)^2.+(z[i]-z)^2.)
  dsort = sort(d)
  dasort = sort(da)
  dbsort = sort(db)
  dcsort = sort(dc)
  Ediff[i] = E[i]-E[dsort[1]]
  sum_Ediff += Ediff[i]^2.
  min_d1[i] = d[dsort[1]]
  min_d2[i] = d[dsort[2]]
  min_d3[i] = d[dsort[3]]
  min_d4[i] = d[dsort[4]]
  min_d5[i] = d[dsort[5]]
  min_d6[i] = d[dsort[6]]
  min_d5a[i] = da[dasort[6]]
  min_d5b[i] = db[dbsort[6]]
  min_d5c[i] = dc[dcsort[6]]
endfor
P_heat = sqrt(sum_Ediff/N)
plot,x,y,psym=3,xr=[0,10],yr=[-10,10],xtitle='x',ytitle='y'
legend,[textoidl('P_{heat}=')+strtrim(P_heat,2)],box=0,/bottom
plothist,Ediff,bin=100
bsize = stddev(E)/6.
plothist,E,xh,yh,bin=bsize
erase
!p.multi=[0,3,2]
multiplot
plot,x,min_d1,psym=1,/ylog
multiplot
plot,x,min_d2,psym=1,/ylog
multiplot
plot,x,min_d3,psym=1,/ylog
multiplot
plot,x,min_d4,psym=1,ytitle='min d to the nearest particle',/ylog
multiplot
plot,x,min_d5,psym=1,/ylog
multiplot
plot,x,min_d6,psym=1,xtitle='x',/ylog
multiplot,/reset
erase

num=40
multiplot
plot,x,1./min_d1,psym=1,/ylog
eqnum_bin,x,1./min_d1,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(1./min_d1),color=200
hline,median(1./min_d1),color=100
multiplot
plot,x,1./min_d2,psym=1,/ylog
eqnum_bin,x,1./min_d2,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(1./min_d2),color=200
hline,median(1./min_d2),color=100
multiplot
plot,x,1./min_d3,psym=1,/ylog
eqnum_bin,x,1./min_d3,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(1./min_d3),color=200
hline,median(1./min_d3),color=100
multiplot
plot,x,1./min_d4,psym=1,ytitle='1/min d to the nearest particle',/ylog
eqnum_bin,x,1./min_d4,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(1./min_d4),color=200
hline,median(1./min_d4),color=100
multiplot
plot,x,1./min_d5,psym=1,xtitle='x',/ylog
eqnum_bin,x,1./min_d5,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(1./min_d5),color=200
hline,median(1./min_d5),color=100

multiplot
plot,x,1./min_d6,psym=1,/ylog
eqnum_bin,x,1./min_d6,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(1./min_d6),color=200
hline,median(1./min_d6),color=100
multiplot,/reset
erase
print,max(1./min_d1)/median(1./min_d1)
print,max(1./min_d2)/median(1./min_d2)
print,max(1./min_d3)/median(1./min_d3)
print,max(1./min_d4)/median(1./min_d4)
print,max(1./min_d5)/median(1./min_d5)
print,max(1./min_d6)/median(1./min_d6)


!p.multi=[0,2,2]
sort_min_d5 = sort(1./min_d5)
ten_percent = round(N_elements(x)*0.1)-1
bottom_d5 = mean(1./min_d5[sort_min_d5[0:ten_percent]])
top_d5 = mean(1./min_d5[(reverse(sort_min_d5))[0:ten_percent]])
clumpy = mean(top_d5)/mean(bottom_d5)
plot,x,1./min_d5,psym=1,xtitle='x',/ylog,title=strtrim(clumpy,2)
eqnum_bin,x,1./min_d5,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(top_d5),color=200,linestyle=2
hline,mean(bottom_d5),color=200,linestyle=2

sort_min_d5a = sort(1./min_d5a)
ten_percent = round(N_elements(x)*0.1)-1
bottom_d5a = mean(1./min_d5a[sort_min_d5a[0:ten_percent]])
top_d5a = mean(1./min_d5a[(reverse(sort_min_d5a))[0:ten_percent]])
clumpy = mean(top_d5a)/mean(bottom_d5a)
plot,x,1./min_d5a,psym=1,xtitle='x',/ylog,title=strtrim(clumpy,2)
eqnum_bin,x,1./min_d5a,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(top_d5a),color=200,linestyle=2
hline,mean(bottom_d5a),color=200,linestyle=2

sort_min_d5b = sort(1./min_d5b)
ten_percent = round(N_elements(x)*0.1)-1
bottom_d5b = mean(1./min_d5b[sort_min_d5b[0:ten_percent]])
top_d5b = mean(1./min_d5b[(reverse(sort_min_d5b))[0:ten_percent]])
clumpy = mean(top_d5b)/mean(bottom_d5b)
plot,x,1./min_d5b,psym=1,xtitle='x',/ylog,title=strtrim(clumpy,2)
eqnum_bin,x,1./min_d5b,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(top_d5b),color=200,linestyle=2

sort_min_d5c = sort(1./min_d5c)
ten_percent = round(N_elements(x)*0.1)-1
bottom_d5c = mean(1./min_d5c[sort_min_d5c[0:ten_percent]])
top_d5c = mean(1./min_d5c[(reverse(sort_min_d5c))[0:ten_percent]])
clumpy = mean(top_d5c)/mean(bottom_d5c)
plot,x,1./min_d5c,psym=1,xtitle='x',/ylog,title=strtrim(clumpy,2)
eqnum_bin,x,1./min_d5c,xm,ym,num=num
oplot,xm,ym,color=255
hline,mean(top_d5c),color=200,linestyle=2

device,/close
stop

device,file='heat_time.ps',/color

dir = '/scratch/jhyoon/Research/LG/snapshot_nosubhalo/'
spawn,'ls '+dir+' > temp'
chr_rdtbl,'temp',0,arr
fname = arr[0,*]
age = float(strmid(fname,4,7))/1000.	; in Gyr
spawn,'rm -f temp'

N_file = N_elements(fname)
P_heat = dblarr(N_file)

for k=0,N_file-1 do begin
  chr_rdtbl,dir+fname[k],0,arr
  arr = double(arr)
  t = arr[0,*]
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  Vx = arr[4,*]
  Vy = arr[5,*]
  Vz = arr[6,*]
  dE = arr[7,*]
  E_total = arr[8,*]
  J = arr[9,*]
  
  ;=== heating parameter ===
  sum_Ediff = 0.d
  N_particle = N_elements(x)
  Ediff = dblarr(N_particle)
  for i=1,N_particle-1 do begin
    d = sqrt((x[i]-x)^2.+(y[i]-y)^2.+(z[i]-z)^2.)
    dsort = sort(d)
    Ediff[i] = E_total[i]-E_total[dsort[1]]
    sum_Ediff += Ediff[i]^2.
  endfor
  P_heat[k] = sqrt(sum_Ediff/N_particle)
endfor
plot,age,P_heat,xtitle='t [Gyr]',ytitle=textoidl('P_{heat}')


device,/close

END
