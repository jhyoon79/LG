pro stochastic

dir_inVL = '/scratch/jhyoon/Research/LG/inVL/'
spawn,'ls '+dir_inVL+'snap1* > tmp'
chr_rdtbl,'tmp',0,fname,/silent
spawn,'rm -f tmp'

set_plot,'ps'
@plot_setting
!p.charsize=1.3
!p.multi=[0,4,3]
device,file='stochastic.ps',/color,/landscape

for i=0,N_elements(fname)-1 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  vx = arr[4,*]
  vy = arr[5,*]
  vz = arr[6,*]
  E_total = arr[8,*]

;=== heating parameter ===
  sum_Ediff = 0.d
  Npts = N_elements(x)
  Ediff = dblarr(Npts)
  for j=1,Npts-1 do begin	; the first particle should be excluded.
    d = sqrt((x[j]-x)^2.+(y[j]-y)^2.+(z[j]-z)^2.)
    dsort = sort(d)
    Ediff[j] = E_total[j]-E_total[dsort[1]]
    sum_Ediff += Ediff[j]^2.
  endfor
  P_heat = sqrt(sum_Ediff/Npts)

  xmin = -14.99;x[0]-10
  xmax = 14.99;x[0]+10
  ymin = -14.99;y[0]-10
  ymax = 14.99;y[0]+10
  xt = ''  &  yt = ''
  if i ge 8 then xt = 'x'
  if i eq 0 or i eq 4 or i eq 8 then yt = 'y'
  multiplot
  plot,x,y,psym=3,xtitle=xt,ytitle=yt,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  legend,['t='+t+'Gyr'],box=0,charsize=1
;  legend,[textoidl('P_{heat}=')+strtrim(P_heat,2)],box=0,charsize=1,/bottom
endfor
plot,[0],[0],/nodata,xsty=4,ysty=4
multiplot,/reset
erase

;=== Observable signature ===
x_sun = -8.d
for i=0,N_elements(fname)-1 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  vx = arr[4,*]
  vy = arr[5,*]
  vz = arr[6,*]
  v = sqrt(vx^2.+vy^2.+vz^2.)
  cos_theta = ((x-x_sun)*vx+y*vy+z*vz) / (sqrt((x-x_sun)^2.+y^2.+z^2.)*abs(v))
  Vradial = v*cos_theta
  E_radial = Vradial^2.

  l = atan(y/(x+8.))*!radeg 
  b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
  lcosb = l*cos(b*!dtor)
  delta = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
  alpha = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(delta*!dtor) )*!radeg + 282.25

;=== heating parameter ===
  sum_Eradial_diff = 0.d
  Npts = N_elements(x)
  Eradial_diff = dblarr(Npts)
  for j=1,Npts-1 do begin	; the first particle should be excluded.
    d_2d = sqrt((alpha[j]-alpha)^2.+(delta[j]-delta)^2.)
    dsort_2d = sort(d_2d)
    Eradial_diff[j] = E_radial[j]-E_radial[dsort_2d[1]]
    sum_Eradial_diff += Eradial_diff[j]^2.
  endfor
  P_heat_radial = sqrt(sum_Eradial_diff/Npts)

  xmin = 359;x[0]-10
  xmax = 1;x[0]+10
  ymin = -89;y[0]-10
  ymax = 89;y[0]+10
  xt = 'RA'
  yt = 'Dec'
  plot,alpha,delta,psym=3,xtitle=xt,ytitle=yt,/isotropic
;  plot,alpha,delta,psym=3,xtitle=xt,ytitle=yt,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  legend,['t='+t+'Gyr'],box=0,charsize=1
;  legend,[textoidl('P_{heat}=')+strtrim(P_heat_radial,2)],box=0,charsize=1,/bottom
endfor
plot,[0],[0],/nodata,xsty=4,ysty=4
multiplot,/reset
erase

for i=0,N_elements(fname)-1 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  vx = arr[4,*]
  vy = arr[5,*]
  vz = arr[6,*]
  v = sqrt(vx^2.+vy^2.+vz^2.)
  cos_theta = ((x-x_sun)*vx+y*vy+z*vz) / (sqrt((x-x_sun)^2.+y^2.+z^2.)*abs(v))
  Vradial = v*cos_theta
  E_radial = Vradial^2.

  l = atan(y/(x+8.))*!radeg 
  b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
  lcosb = l*cos(b*!dtor)
  delta = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
  alpha = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(delta*!dtor) )*!radeg + 282.25
  xmin = 359;x[0]-10
  xmax = 1;x[0]+10
  ymin = -89;y[0]-10
  ymax = 89;y[0]+10
  xt = 'RA'
  yt = 'Dec'
  contour_bin = 10
  levels = indgen(300)*contour_bin + contour_bin/2
  pts_contour,alpha,delta,20,40,gv,x,y
  contour,gv,x,y,xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1
  legend,['t='+t+'Gyr'],box=0,charsize=1
;  legend,[textoidl('P_{heat}=')+strtrim(P_heat_radial,2)],box=0,charsize=1,/bottom
endfor
erase

!p.multi=[0,2,2]
sub = [0,5,8,9]
for ii=0,3 do begin
  i = sub[ii]
  chr_rdtbl,fname[i],0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  vx = arr[4,*]
  vy = arr[5,*]
  vz = arr[6,*]
  E_total = arr[8,*]

;=== heating parameter ===
  sum_Ediff = 0.d
  Npts = N_elements(x)
  Ediff = dblarr(Npts)
  for j=1,Npts-1 do begin	; the first particle should be excluded.
    d = sqrt((x[j]-x)^2.+(y[j]-y)^2.+(z[j]-z)^2.)
    dsort = sort(d)
    Ediff[j] = E_total[j]-E_total[dsort[1]]
    sum_Ediff += Ediff[j]^2.
  endfor
  P_heat = sqrt(sum_Ediff/Npts)

  xmin = -14.99;x[0]-10
  xmax = 14.99;x[0]+10
  ymin = -14.99;y[0]-10
  ymax = 14.99;y[0]+10
  xt = ''  &  yt = ''
  if i ge 8 then xt = 'x'
  if i eq 0 or i eq 4 or i eq 8 then yt = 'y'
  multiplot
  plot,x,y,psym=3,xtitle=xt,ytitle=yt,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  legend,['t='+t+'Gyr'],box=0,charsize=1
;  legend,[textoidl('P_{heat}=')+strtrim(P_heat,2)],box=0,charsize=1,/bottom
endfor
plot,[0],[0],/nodata,xsty=4,ysty=4
multiplot,/reset
erase




device,/close

END
