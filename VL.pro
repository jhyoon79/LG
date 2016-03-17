pro VL

chr_rdtbl,'vltwosubs.txt',3,arr
arr = double(arr)
Mtidal = arr[5,*]
arr = arr[*,reverse(sort(Mtidal))]
id = arr[0,*]
GCdistance = arr[1,*]
peakVmax = arr[2,*]
Vmax = arr[3,*]
rVmax = arr[4,*]
rs = rVmax/2.
Mtidal = arr[5,*]
rtidal = arr[6,*]
x = arr[7,*]
y = arr[8,*]
z = arr[9,*]
Vx = arr[10,*]
Vy = arr[11,*]
Vz = arr[12,*]
Rgc = sqrt(x^2.+y^2.+z^2.)
print,'r_min=',min(Rgc),'    r_max=',max(Rgc)
print,'M_min=',min(Mtidal),'    M_max=',max(Mtidal)
;dir_out = '/scratch/jhyoon/Research/LG/snapshot_VLsubhalo/'
;chr_rdtbl,dir_out+'subhalo0000',0,arr
;arr = double(arr)
;x_subhalo0000 = arr[1,*]
;y_subhalo0000 = arr[2,*]
;z_subhalo0000 = arr[3,*]
;r_subhalo0000 = sqrt(x_subhalo0000^2.+y_subhalo0000^2.+z_subhalo0000^2.)
;chr_rdtbl,dir_out+'subhalo1500',0,arr
;arr = double(arr)
;x_subhalo1500 = arr[1,*]
;y_subhalo1500 = arr[2,*]
;z_subhalo1500 = arr[3,*]
;r_subhalo1500 = sqrt(x_subhalo1500^2.+y_subhalo1500^2.+z_subhalo1500^2.)
;chr_rdtbl,dir_out+'subhalo2950',0,arr
;arr = double(arr)
;x_subhalo2950 = arr[1,*]
;y_subhalo2950 = arr[2,*]
;z_subhalo2950 = arr[3,*]
;r_subhalo2950 = sqrt(x_subhalo2950^2.+y_subhalo2950^2.+z_subhalo2950^2.)

chr_rdtbl,'stepToTimeVL2.txt',3,arr
redshift = reform(double(arr[2,*]))

;=== readout past subhalo mass ===
chr_rdtbl,'progMtidal.txt',0,arr
progMtidal = double(arr)
chr_rdtbl,'progGCdistance.txt',0,arr
progRgc = double(arr)

z0 = redshift[0]
Mtidal_z0 = reform(progMtidal[0,*])
Rgc_z0 = reform(progRgc[0,*])
z1 = redshift[11]
Mtidal_z1 = reform(progMtidal[11,*])
Rgc_z1 = reform(progRgc[11,*])
z2 = redshift[15]
Mtidal_z2 = reform(progMtidal[15,*])
Rgc_z2 = reform(progRgc[15,*])
help,z0,z1,z2
sub0 = where(Mtidal_z0 gt 0)
sub1 = where(Mtidal_z1 gt 0)
sub2 = where(Mtidal_z2 gt 0)
sub0_R20 = where(Mtidal_z0 gt 0 and Rgc_z0 le 20)
sub1_R20 = where(Mtidal_z1 gt 0 and Rgc_z1 le 20)
sub2_R20 = where(Mtidal_z2 gt 0 and Rgc_z2 le 20)
sub0_R50 = where(Mtidal_z0 gt 0 and Rgc_z0 le 50)
sub1_R50 = where(Mtidal_z1 gt 0 and Rgc_z1 le 50)
sub2_R50 = where(Mtidal_z2 gt 0 and Rgc_z2 le 50)

G = 4.3e-6 ; (km/s)^2 M^-1 kpc
rvir = 389.
rvir_sub = where(Rgc le rvir)
rvir_sub1 = where(Rgc le 25)
rvir_sub2 = where(Rgc le 40)
rvir_sub3 = where(Rgc le 30)
rvir_sub4 = where(Rgc le 25)
rvir_sub5 = where(Rgc le 20)

V = sqrt(Vx^2.+Vy^2.+Vz^2.)
print,stddev(V[rvir_sub1])

@plot_setting
device,file='VL_f1.eps',/color,/enc,/cmyk,ysize=15;,yoffset=0.5
bsize = 0.5
!p.multi=[0,1,3]
!p.charsize=1.3
loadct,0
xt = textoidl('log M_{subhalo} [M')+sunsymbol()+']'
multiplot
plot,[0],[0],/nodata,xr=[2,11],yr=[0.6,2e4],ytitle='#',/ylog
plothist,alog10(Mtidal),bin=bsize,color=150,/overplot
plothist,alog10(Mtidal[rvir_sub1]),xh,yh,bin=bsize,color=0,/overplot
legend,['All','r < 25kpc'],color=[150,0],linestyle=0,box=0,charsize=1
print,xh,yh

multiplot
plot,[0],[0],/nodata,ytitle=textoidl('r_{tidal} [kpc]'),xr=[2,11],yr=[3e-2,2e2],/ylog
oplot,alog10(Mtidal),rtidal,psym=1,color=150
oplot,alog10(Mtidal[rvir_sub1]),rtidal[rvir_sub1],psym=1
legend,['All','r < 25kpc'],color=[150,0],psym=1,box=0,charsize=1
fit_sub = where(Rgc le 25 and alog10(Mtidal) gt 5)
fit= poly_fit(alog10(Mtidal[fit_sub]),alog10(rtidal[fit_sub]),1)
xfit = [2,11]
yfit = fit[0]+fit[1]*xfit
;oplot,xfit,10^yfit,color=100
print,'rtidal fit',fit

multiplot
plot,[0],[0],/nodata,xtitle=xt,ytitle=textoidl('r_{s} [kpc]'),xr=[2,11],yr=[1e-2,1e1],/ylog
oplot,alog10(Mtidal),rs,psym=1,color=150
oplot,alog10(Mtidal[rvir_sub1]),rs[rvir_sub1],psym=1
fit_sub = where(Rgc le 25 and alog10(Mtidal) gt 6)
fit= poly_fit(alog10(Mtidal[fit_sub]),alog10(rs[fit_sub]),1)
xfit = [2,11]
yfit = fit[0]+fit[1]*xfit
;oplot,xfit,10^yfit,color=100
print,'rs fit',fit
device,/close

spawn,'convert VL_f1.eps -resize 700x700 VL_f1.small.eps'

loadct,13
erase
!p.multi=0
device,file='VL.ps',/color,ysize=18,yoffset=1
plothist,alog10(Mtidal),xh0,yh0,bin=bsize,xr=[1.5,11.9],yr=[0.8,50000],ytitle='#',/ylog
plothist,alog10(Mtidal[rvir_sub]),xh1,yh1,bin=bsize,color=30,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub1]),xh1,yh1,bin=bsize,color=255,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub2]),xh2,yh2,bin=bsize,color=220,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub3]),xh3,yh3,bin=bsize,color=150,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub4]),xh4,yh4,bin=bsize,color=100,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub5]),xh5,yh5,bin=bsize,color=70,/overplot,/ylog
legend,['R<50kpc','R<20kpc'],box=0,color=[255,100],linestyle=0,/right,charsize=1
;plot,[0],[0],/nodata,xr=[1.5,11.9],yr=[0.8,50000],/ylog
;i1 = 8  &  i2 = 15
i1 = 6  &  i2 = 10
help,xh0,yh0
fit0 = poly_fit(xh0[i1:i2],alog10(yh0[i1:i2]),1,yfit=yfit)
oplot,xh0[i1:i2],10^yfit
i1 = 5  &  i2 = 8
fit1 = poly_fit(xh1[i1:i2],alog10(yh1[i1:i2]),1,yfit=yfit)
oplot,xh1[i1:i2],10^yfit,color=255
fit2 = poly_fit(xh2[i1:i2],alog10(yh2[i1:i2]),1,yfit=yfit)
oplot,xh2[i1:i2],10^yfit,color=220
fit3 = poly_fit(xh3[i1:i2],alog10(yh3[i1:i2]),1,yfit=yfit)
oplot,xh3[i1:i2],10^yfit,color=150
i1 = 4  &  i2 = 8
fit4 = poly_fit(xh4[i1:i2],alog10(yh4[i1:i2]),1,yfit=yfit4)
oplot,xh4[i1:i2],10^yfit4,color=100
i1 = 4  &  i2 = 6
fit5 = poly_fit(xh5[i1:i2],alog10(yh5[i1:i2]),1,yfit=yfit)
oplot,xh5[i1:i2],10^yfit,color=70
print,'slope',fit1[1],fit2[1],fit3[1],fit4[1],fit5[1]
xt1 = [5.5,6.5,7.5,8.5]
yt1 = [1824,950,187,28]
yt2 = [1824*48,950*3,187,28]
fit5 = poly_fit(xt1,alog10(yt1),1,yfit=yfit)
oplot,xt1,10^yfit,linestyle=2
fit5 = poly_fit(xt1,alog10(yt2),1,yfit=yfit)
oplot,xt1,10^yfit,linestyle=2
;oplot,xh4[1:4],10^yfit4*25,linestyle=1
vline,alog10(4e6),0.1,1e6,linestyle=1
erase

!p.multi=[0,2,2]
multiplot
plothist,alog10(Mtidal),xh0,yh0,bin=bsize,xr=[1.5,11.9],yr=[0.8,50000],ytitle='#',/ylog
plothist,alog10(Mtidal[rvir_sub1]),xh1,yh1,bin=bsize,color=255,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub2]),xh2,yh2,bin=bsize,color=220,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub3]),xh3,yh3,bin=bsize,color=150,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub4]),xh4,yh4,bin=bsize,color=100,/overplot,/ylog
plothist,alog10(Mtidal[rvir_sub5]),xh5,yh5,bin=bsize,color=70,/overplot,/ylog
legend,['R<50kpc','R<20kpc'],box=0,color=[255,100],linestyle=0,/right,charsize=1

multiplot
plothist,alog10(Mtidal_z0[sub0]),bin=1,xr=[1.5,11.9],yr=[0.8,50000],/ylog
plothist,alog10(Mtidal_z1[sub1]),bin=1,color=70,/overplot
plothist,alog10(Mtidal_z2[sub2]),bin=1,color=255,/overplot
legend,['z=0','z=1','z=2'],color=[0,70,255],linestyle=0,box=0,/right,charsize=1

multiplot
plothist,alog10(Mtidal_z0[sub0_R50]),bin=1,xr=[1.5,11.9],yr=[0.8,50000],$
  xtitle=textoidl('log M_{tidal}'),ytitle='#',/ylog
plothist,alog10(Mtidal_z1[sub1_R50]),bin=1,color=70,/overplot
plothist,alog10(Mtidal_z2[sub2_R50]),bin=1,color=255,/overplot
legend,['R<50kpc'],box=0,/right,charsize=1

multiplot
plothist,alog10(Mtidal_z0[sub0_R20]),bin=1,xr=[1.5,11.9],yr=[0.8,50000],$
  xtitle=textoidl('log M_{tidal}'),/ylog
plothist,alog10(Mtidal_z1[sub1_R20]),bin=1,color=70,/overplot
plothist,alog10(Mtidal_z2[sub2_R20]),bin=1,color=255,/overplot
legend,['R<20kpc'],box=0,/right,charsize=1
multiplot,/reset
erase


!p.multi=0
;plothist,alog10(Mtidal)/mean(alog10(Mtidal)),bin=0.1,xtitle=textoidl('log M_{tidal}'),ytitle='#'

plothist,alog10(rtidal),bin=0.2,xtitle=textoidl('log r_{tidal}'),ytitle='#'
plothist,alog10(rtidal)/mean(alog10(rtidal)),bin=0.2,xtitle=textoidl('log r_{tidal}'),ytitle='#'

plothist,Vmax,bin=1,xtitle=textoidl('V_{max}'),yr=[0.8,2000],/ylog
plothist,rVmax,bin=0.2,xtitle=textoidl('r_{Vmax}'),yr=[0.8,5000],/ylog

plothist,Rgc,bin=50,xtitle='Galactocentric radius [kpc]',ytitle='#'
vline,rvir,linestyle=2,color=90
plothist,Rgc,bin=10,xr=[0,400],xtitle='Galactocentric radius [kpc]',ytitle='#'
vline,rvir,linestyle=2,color=90

plot,Rgc,Mtidal,psym=1,xtitle='Galactocentric radius [kpc]', $
	ytitle=textoidl('M_{tidal}'),xr=[1,10000],yr=[1e4,1e12],/xlog,/ylog
vline,rvir,1e-4,1e20,linestyle=2,color=90
plot,Rgc,rtidal,psym=1,xtitle='Galactocentric radius [kpc]', $
	ytitle=textoidl('r_{tidal}'),xr=[1,10000],yr=[1e-2,1e3],/xlog,/ylog
vline,rvir,1e-4,1e20,linestyle=2,color=90

plot,Mtidal,Vmax,psym=1,xtitle=textoidl('M_{tidal}'),ytitle=textoidl('V_{max}'),/xlog,/ylog
plot,rtidal,Vmax,psym=1,xtitle=textoidl('r_{tidal}'),ytitle=textoidl('V_{max}'),/xlog,/ylog
plot,Mtidal,rVmax/2.,psym=1,xtitle=textoidl('M_{tidal}'),ytitle=textoidl('r_s (=r_{Vmax}/2)'),/xlog,/ylog
plot,rtidal,rVmax,psym=1,xtitle=textoidl('r_{tidal}'),ytitle=textoidl('r_{Vmax}'),/xlog,/ylog

plothist,rtidal/rs,xtitld='c (=rtidal/rs)',bin=1

plot,Rgc,Vmax,psym=1,xtitle='Galactocentric radius [kpc]',ytitle=textoidl('V_{max}'),/ylog
plot,Rgc,rVmax,psym=1,xtitle='Galactocentric radius [kpc]',ytitle=textoidl('r_{Vmax}'),/ylog



plot,Mtidal,rtidal,psym=1,ytitle=textoidl('r_{tidal} [kpc]'),xtitle=textoidl('M_{tidal}'), xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
fit = poly_fit(alog10(Mtidal),alog10(rtidal),1)
x4fit = alog10([1e1,1e15])
y4fit = fit[0] + fit[1]*x4fit
oplot,10^x4fit,10^y4fit,color=255
str_fit = strtrim(fit,2)
legend,['y='+str_fit[0]+'+'+str_fit[1]+'*x'],box=0
print, fit

r1 = where(Rgc ge 0 and Rgc lt 100)
r2 = where(Rgc ge 100 and Rgc lt 200)
r3 = where(Rgc ge 200 and Rgc lt 300)
r4 = where(Rgc ge 300 and Rgc lt 400)
r5 = where(Rgc ge 400 and Rgc lt 600)
r6 = where(Rgc ge 600 and Rgc lt 1000)
r7 = where(Rgc ge 1000 and Rgc lt 1400)


;=== estimate r_tidal of a subhalo with a certain mass from Via lactea ===
erase
!p.multi=[0,3,3]
multiplot
plot,[0],[0],/nodata,ytitle=textoidl('r_{tidal} [kpc]'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r1],rtidal[r1],psym=1,color=30
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['0<r<100'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r2],rtidal[r2],psym=1,color=70
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['100<r<200'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r3],rtidal[r3],psym=1,color=100
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['200<r<300'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,ytitle=textoidl('r_{tidal} [kpc]'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r4],rtidal[r4],psym=1,color=150
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['300<r<400'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r5],rtidal[r5],psym=1,color=200
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['400<r<600'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r6],rtidal[r6],psym=1,color=230
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['600<r<1000'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('M_{tidal}'), $
	ytitle=textoidl('r_{tidal} [kpc]'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r7],rtidal[r7],psym=1,color=255
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['1000<r<1400'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('M_{tidal}'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r1],rtidal[r1],psym=1,color=30
oplot,Mtidal[r2],rtidal[r2],psym=1,color=70
oplot,Mtidal[r3],rtidal[r3],psym=1,color=100
oplot,Mtidal[r4],rtidal[r4],psym=1,color=150
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['0<r<400'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('M_{tidal}'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r5],rtidal[r5],psym=1,color=200
oplot,Mtidal[r6],rtidal[r6],psym=1,color=230
oplot,Mtidal[r7],rtidal[r7],psym=1,color=255
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['400<r<1400'],box=0,charsize=1,/bottom,/right
multiplot,/reset
erase
!p.multi=0


;=== estimate r_s of a subhalo with a certain mass from Via lactea ===
erase
!p.multi=[0,3,3]
multiplot
plot,[0],[0],/nodata,ytitle=textoidl('r_s [kpc]'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r1],rs[r1],psym=1,color=30
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['0<r<100'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r2],rs[r2],psym=1,color=70
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['100<r<200'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r3],rs[r3],psym=1,color=100
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['200<r<300'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,ytitle=textoidl('r_s [kpc]'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r4],rs[r4],psym=1,color=150
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['300<r<400'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r5],rs[r5],psym=1,color=200
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['400<r<600'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r6],rs[r6],psym=1,color=230
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['600<r<1000'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('M_{tidal}'), $
	ytitle=textoidl('r_s [kpc]'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r7],rs[r7],psym=1,color=255
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['1000<r<1400'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('M_{tidal}'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r1],rs[r1],psym=1,color=30
oplot,Mtidal[r2],rs[r2],psym=1,color=70
oplot,Mtidal[r3],rs[r3],psym=1,color=100
oplot,Mtidal[r4],rs[r4],psym=1,color=150
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['0<r<400'],box=0,charsize=1,/bottom,/right
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('M_{tidal}'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,Mtidal[r5],rs[r5],psym=1,color=200
oplot,Mtidal[r6],rs[r6],psym=1,color=230
oplot,Mtidal[r7],rs[r7],psym=1,color=255
oplot,[1e4,1e12],[1e-2,1e3],linestyle=1
legend,['400<r<1400'],box=0,charsize=1,/bottom,/right
multiplot,/reset
erase
!p.multi=0


help,rvir_sub
plot,Mtidal[rvir_sub],rtidal[rvir_sub],psym=1,ytitle=textoidl('r_{tidal} [kpc]'),xtitle=textoidl('M_{tidal}'), $
	xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
fit = poly_fit(alog10(Mtidal[rvir_sub]),alog10(rtidal[rvir_sub]),1,yfit=yfit)
x4fit = alog10([1e1,1e15])
y4fit = fit[0] + fit[1]*x4fit
oplot,10^x4fit,10^y4fit,color=255
str_fit = strtrim(fit,2)
legend,['y='+str_fit[0]+'+'+str_fit[1]+'*x'],box=0
print, 'within rvir',fit

normalized_log_rtidal = alog10(rtidal[rvir_sub])-yfit
plot,Mtidal[rvir_sub],normalized_log_rtidal,psym=1,xr=[1e4,1e12],/xlog
tmp_sub = where(alog10(Mtidal[rvir_sub]) gt 5.8 and alog10(Mtidal[rvir_sub]) lt 6)
plothist,normalized_log_rtidal,bin=0.1
plothist,normalized_log_rtidal[tmp_sub],bin=0.1,color=255,/overplot


;=== estimate r_tidal of a certain mass from the Via Lactea Mtidal-rtidal relation ===
;chr_rdtbl,'BI5',1,arr
;arr = double(arr)
;rscale = 24.6/0.5
;x_BI5 = arr[1,*]
;y_BI5 = arr[2,*]
;z_BI5 = arr[3,*]
;r_BI5 = sqrt(x_BI5^2.+y_BI5^2.+z_BI5^2.)*rscale
m_VL = moment(alog10(Mtidal[rvir_sub]))
N_VL = N_elements(rvir_sub)

;device,file='tmp.ps',/color
;slope=0
;for k=0,100 do begin

art_mass = randomn(iseed,N_VL)*sqrt(m_VL[1]) + m_VL[0]
;art_r = r_BI5[round(randomu(iseed,N_VL)*(N_VL-1))]
art_r = Rgc[rvir_sub[round(randomu(iseed,N_VL)*(N_VL-1))]]
plothist,alog10(Mtidal[rvir_sub]),bin=0.5,xtitle=textoidl('log M_{tidal}'),ytitle='#'
plothist, art_mass,bin=0.5,/overplot,color=255
art_mass = 10^art_mass
art_rtidal = dblarr(N_VL)
mean_log_Mtidal = mean(alog10(Mtidal[rvir_sub]))
mean_r = mean(Rgc[rvir_sub])
for i=0,N_VL-1 do begin
  r_sub = where(abs(art_r[i]-Rgc[rvir_sub]) lt 90)
  tmp = min(abs(art_mass[i]-Mtidal[rvir_sub[r_sub]]),min_sub)
  art_rtidal[i] = rtidal[rvir_sub[r_sub[min_sub]]]
endfor
plot,Mtidal[rvir_sub],rtidal[rvir_sub],psym=1,ytitle=textoidl('r_{tidal} [kpc]'),xtitle=textoidl('M_{tidal}'),xr=[1e4,1e12],yr=[1e-2,1e3],/xlog,/ylog
oplot,10^x4fit,10^y4fit,color=255
oplot,art_mass,art_rtidal,psym=7,color=100
fit = poly_fit(alog10(art_mass),alog10(art_rtidal),1)
x4fit = alog10([1e1,1e15])
y4fit = fit[0] + fit[1]*x4fit
oplot,10^x4fit,10^y4fit,color=70
str_fit = strtrim(fit,2)
legend,['y='+str_fit[0]+'+'+str_fit[1]+'*x'],box=0
print, 'art',fit

;slope = [slope,fit[1]]
;endfor
;plothist,slope[1:N_elements(slope)-1],bin=0.002
;device,/close
;stop





rbin = 1.
dr = findgen(5000)*rbin
;Nrho_dr = dblarr(N_elements(dr))
rho_dr_all = dblarr(N_elements(dr))
;rho_dr0000 = dblarr(N_elements(dr))
;rho_dr1500 = dblarr(N_elements(dr))
;rho_dr2950 = dblarr(N_elements(dr))
for i=0,N_elements(dr)-1 do begin
  in_dr_all = where(Rgc ge dr[i] and Rgc lt dr[i]+rbin)
;  in_dr0000 = where(r_subhalo0000 ge dr[i] and r_subhalo0000 lt dr[i]+rbin)
;  in_dr1500 = where(r_subhalo1500 ge dr[i] and r_subhalo1500 lt dr[i]+rbin)
;  in_dr2950 = where(r_subhalo2950 ge dr[i] and r_subhalo2950 lt dr[i]+rbin)
  area_dr = 4./3.*!pi*((dr[i]+rbin)^3.-dr[i]^3.)
  if (in_dr_all[0] ne -1) then rho_dr_all[i] = total(Mtidal[in_dr_all])/area_dr
;  if (in_dr0000[0] ne -1) then rho_dr0000[i] = total(Mtidal[in_dr0000])/area_dr
;  if (in_dr1500[0] ne -1) then rho_dr1500[i] = total(Mtidal[in_dr1500])/area_dr
;  if (in_dr2950[0] ne -1) then rho_dr2950[i] = total(Mtidal[in_dr2950])/area_dr
;    Nrho_dr[i] = Ndr/area_dr
endfor

plot,[0],[0],/nodata,xtitle='r [kpc]',ytitle=textoidl('\rho / M'+sunsymbol()+'kpc^{-3}'),xr=[0.4,8000],yr=[1e-3,1e+10],/xlog,/ylog
oplot,dr+rbin/2.,rho_dr_all,psym=1
;oplot,dr+rbin/2.,rho_dr0000,psym=4
;oplot,dr+rbin/2.,rho_dr1500,psym=5,color=70
;oplot,dr+rbin/2.,rho_dr2950,psym=6,color=255

r = dr+rbin/2.
rs = 24.6
c = rvir/rs
delta_c = 200./3.*c^3./(alog(1.+c)-c/(1.+c))
H = 73. /1000. ; km/s/kpc
rho_cri = 3.*H^2./(8.*!pi*G)
rho = delta_c / ((r/rs)*(1.+r/rs)^2.) * rho_cri
oplot,r,rho,color=150
vline,rvir,1e-19,1e19,linestyle=2,color=90

rvir = 150.
c = rvir/rs
delta_c = 200./3.*c^3./(alog(1.+c)-c/(1.+c))
rho = delta_c / ((r/rs)*(1.+r/rs)^2.) * rho_cri
oplot,r,rho,color=100
legend,['NFW(rvir=389)','NFW(rvir=150)','0.00Gyr','1.50Gyr','2.95Gyr'], $
	psym=[0,0,1,5,6],color=[150,100,0,70,255],box=0,/right


;;=== orbits of Via Lactea subhalos ===
;spawn, 'ls snapshot_VLsubhalo > temp2'
;chr_rdtbl, 'temp2', 0, fname
;fname = 'snapshot_VLsubhalo/'+reform(fname)
;spawn, 'rm -f temp2'
;N_pts = N_elements(fname)-2
;N_subhalos = N_elements(Mtidal);6553
;t_subhalo = dblarr(N_pts,N_subhalos)
;x_subhalo = dblarr(N_pts,N_subhalos)
;y_subhalo = dblarr(N_pts,N_subhalos)
;z_subhalo = dblarr(N_pts,N_subhalos)
;vx_subhalo = dblarr(N_pts,N_subhalos)
;vy_subhalo = dblarr(N_pts,N_subhalos)
;vz_subhalo = dblarr(N_pts,N_subhalos)
;for i=0,N_pts-1 do begin
;  chr_rdtbl,fname[i],0,arr,/silent
;  arr = double(arr)
;  t_subhalo[i,*] = arr[0,*]
;  x_subhalo[i,*] = arr[1,*]
;  y_subhalo[i,*] = arr[2,*]
;  z_subhalo[i,*] = arr[3,*]
;  vx_subhalo[i,*] = arr[4,*]
;  vy_subhalo[i,*] = arr[5,*]
;  vz_subhalo[i,*] = arr[6,*]
;endfor
;
;r_subhalo = sqrt(x_subhalo^2.+y_subhalo^2.+z_subhalo^2.)
;rxy_subhalo = sqrt(x_subhalo^2.+y_subhalo^2.)
;  sub_2950 = where(t_subhalo eq 2950.)
;xmin = -1000
;xmax = 1000
;ymin = -1000
;ymax = 1000
;zmin = -1000
;zmax = 1000
;rmin = 0
;rmax = 1000
;;erase
;!p.multi=[0,5,5]
;;plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',/isotropic
;for i=0,N_subhalos-1 do begin
;;for i=0,1000-1 do begin
;  if r_subhalo[sub_2950[i]] le 500 then begin
;  plot,x_subhalo[*,i],y_subhalo[*,i],color=70,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',/isotropic
;;  oplot,x_subhalo[*,i],y_subhalo[*,i],color=70
;  oplot,[x_subhalo[sub_2950[i]]],[y_subhalo[sub_2950[i]]],psym=7,color=255
;  oplot,[x_subhalo[0,i]],[y_subhalo[0,i]],psym=5,color=215
;  endif
;endfor
;plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='y [kpc]',ytitle='z [kpc]',/isotropic
;for i=0,N_subhalos-1 do begin
;  oplot,y_subhalo[*,i],z_subhalo[*,i],color=70
;  oplot, [y_subhalo[sub_2950,i]],[z_subhalo[sub_2950,i]],psym=7,color=255
;  oplot, [y_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
;endfor
;plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='z [kpc]',/isotropic
;for i=0,N_subhalos-1 do begin
;  oplot,x_subhalo[*,i],z_subhalo[*,i],color=70
;  oplot, [x_subhalo[sub_2950,i]],[z_subhalo[sub_2950,i]],psym=7,color=255
;  oplot, [x_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
;endfor
;plot,[0],[0],/nodata,xr=[rmin,rmax],yr=[zmin,zmax],xtitle='r [kpc]', ytitle='z [kpc]', /isotropic
;for i=0,N_subhalos-1 do begin
;  oplot,rxy_subhalo[*,i],z_subhalo[*,i],color=70
;  oplot, [rxy_subhalo[sub_2950,i]],[z_subhalo[sub_2950,i]],psym=7,color=255
;  oplot, [rxy_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
;endfor
;erase

device,/close

END
