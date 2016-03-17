pro Pal5_plot

;=== print out starting time ===
readcol,'time_start',t1
t2 = systime(1)
dt = (t2 - t1)/60.	; min
openw,32,'tmp_time'
printf,32,dt
close,32

;=== read output of the first particle ===
chr_rdtbl, 'part001', 0, arr
arr = double(arr)
t_part1 = arr[0,*]
x_part1 = arr[1,*]
y_part1 = arr[2,*]
z_part1 = arr[3,*]
Vx_part1 = arr[4,*]
Vy_part1 = arr[5,*]
Vz_part1 = arr[6,*]
delE_part1 = arr[7,*]
Etot_part1 = arr[8,*]
delJ_part1 = arr[9,*]
r_part1 = sqrt(x_part1^2.+y_part1^2.)

;Vr viewed at the Sun
x_sun = -8.d
Vr_part1 = reform( (Vx_part1*(x_part1-x_sun)+Vy_part1*y_part1+Vz_part1*z_part1)/sqrt((x_part1-x_sun)^2.+y_part1^2.+z_part1^2.) )
print,'recession velocity of Pal5=',Vr_part1[where(t_part1 eq 3200)]

N_arr = N_elements(t_part1)
t_init = t_part1[0]
x_init = x_part1[0]  &  Vx_init = Vx_part1[0]
y_init = y_part1[0]  &  Vy_init = Vy_part1[0]
z_init = z_part1[0]  &  Vz_init = Vz_part1[0]
r_init = r_part1[0];  & ; Vr_init = Vr_part1[0]
t_1Gyr = where(abs(t_part1) eq 1000.)
t_2Gyr = where(abs(t_part1) eq 2000.)
t_3_20Gyr = where(abs(t_part1) eq 3200.)
;t_3_20Gyr = where(abs(t_part1) eq 3350.)
;t_3_20Gyr = where(abs(t_part1) eq 2950.)
t_500Myr = where(abs(t_part1) eq 500.)
x_1Gyr = x_part1[t_1Gyr]  &  Vx_1Gyr = Vx_part1[t_1Gyr]
y_1Gyr = y_part1[t_1Gyr]  &  Vy_1Gyr = Vy_part1[t_1Gyr]
z_1Gyr = z_part1[t_1Gyr]  &  Vz_1Gyr = Vz_part1[t_1Gyr]
r_1Gyr = r_part1[t_1Gyr];  &  Vr_1Gyr = Vr_part1[t_1Gyr]
x_2Gyr = x_part1[t_2Gyr]  &  Vx_2Gyr = Vx_part1[t_2Gyr]
y_2Gyr = y_part1[t_2Gyr]  &  Vy_2Gyr = Vy_part1[t_2Gyr]
z_2Gyr = z_part1[t_2Gyr]  &  Vz_2Gyr = Vz_part1[t_2Gyr]
r_2Gyr = r_part1[t_2Gyr];  &  Vr_2Gyr = Vr_part1[t_2Gyr]
x_500Myr = x_part1[t_500Myr]  &  Vx_500Myr = Vx_part1[t_500Myr]
y_500Myr = y_part1[t_500Myr]  &  Vy_500Myr = Vy_part1[t_500Myr]
z_500Myr = z_part1[t_500Myr]  &  Vz_500Myr = Vz_part1[t_500Myr]
r_500Myr = r_part1[t_500Myr];  &  Vr_500Myr = Vr_part1[t_500Myr]
x_3_20Gyr = x_part1[t_3_20Gyr]  &  Vx_3_20Gyr = Vx_part1[t_3_20Gyr]
y_3_20Gyr = y_part1[t_3_20Gyr]  &  Vy_3_20Gyr = Vy_part1[t_3_20Gyr]
z_3_20Gyr = z_part1[t_3_20Gyr]  &  Vz_3_20Gyr = Vz_part1[t_3_20Gyr]
r_3_20Gyr = r_part1[t_3_20Gyr];  &  Vr_3_20Gyr = Vr_part1[t_3_20Gyr]
t_final = t_part1[N_arr-1]
x_final = x_part1[N_arr-1]  &  Vx_final = Vx_part1[N_arr-1]
y_final = y_part1[N_arr-1]  &  Vy_final = Vy_part1[N_arr-1]
z_final = z_part1[N_arr-1]  &  Vz_final = Vz_part1[N_arr-1]
r_final = r_part1[N_arr-1]  &;  Vr_final = Vr_part1[N_arr-1]
t_mark_1Gyr = 't=1Gyr'
t_mark_2Gyr = 't=2Gyr'
t_mark_3_20Gyr = 't=3.20Gyr'
t_mark_500Myr = 't=500Myr'
t_mark_final = textoidl('t_{final}=')+strtrim(string(t_final/1000.,f='(f6.2)'),2)+'Gyr'
xmin = -20  &  xmax = 20
ymin = -20  &  ymax = 20
zmin = -20  &  zmax = 20
rmin = 0    &  rmax = 25
Vxmin = -250 &  Vxmax = 250
Vymin = -250 &  Vymax = 250
Vzmin = -250 &  Vzmax = 250
;Vrmin = 0  &  Vrmax = 250

openr, 99, 'frogin'
readf, 99, N_stars
readf, 99, tmp,tmp
readf, 99, tmp
readf, 99, N_subhalos
close, 99

;=== read out info of all test particles ===  
chr_rdtbl, 'frog_Pal5_xy.dat', 0, arr
x = strtrim(string(double(arr[0,0]),f='(f5.1)'),2)
y = strtrim(string(double(arr[1,0]),f='(f5.1)'),2)
z = strtrim(string(double(arr[2,0]),f='(f5.1)'),2)
Vx = strtrim(fix(arr[3,0]),2)
Vy = strtrim(fix(arr[4,0]),2)
Vz = strtrim(fix(arr[5,0]),2)
Info_posi = 'x='+x+' y='+y+' z='+z
Info_vel = textoidl('V_x='+Vx+' V_y='+Vy+' V_z='+Vz)
Info1 = 'Part0001'+textoidl(' N_{particle}=900')

arr = double(arr)
x_GC = arr[0,*]
y_GC = arr[1,*]
z_GC = arr[2,*]
Vx_GC = arr[3,*]
Vy_GC = arr[4,*]
Vz_GC = arr[5,*]
r_GC = sqrt( (x_GC-x_GC[0])^2.+(y_GC-y_GC[0])^2.+(z_GC-z_GC[0])^2. )
r_GC = r_GC/23.2*60.*!radeg
N_bin = 100
r_size = 1.
r_bin = findgen(N_bin)*r_size
surface_density = dblarr(N_bin)
surface_density_err = dblarr(N_bin)
for i=0,N_bin-1 do begin
  inR = where(r_GC ge r_bin[i] and r_GC lt r_bin[i]+r_size)
  dN_dr = (inR[0] eq -1)? 0:N_elements(inR)
  surface_density[i] = dN_dr / (!pi*(r_bin[i]+r_size/2.)^2.)
  surface_density_err[i] = sqrt(dN_dr) / (!pi*(r_bin[i]+r_size/2.)^2.)
endfor

@plot_setting
!p.charsize=1.2
device, file='Pal5_plot.ps', /color,/landscape

;=== initial position of particles: the shape of Pal5 at the initial position ===
!p.multi = [0,2,2]
plot, x_GC, y_GC, psym=3, xtitle='x [kpc]', ytitle='y [kpc]', /isotropic, title='initial position of particles'
plot, y_GC, z_GC, psym=3, xtitle='x [kpc]', ytitle='y [kpc]', /isotropic, title='initial position of particles'
plot, z_GC, x_GC, psym=3, xtitle='x [kpc]', ytitle='y [kpc]', /isotropic, title='initial position of particles'
multiplot,/reset
erase

;;=== reproduce of Odenkirchen+03 ===
;!p.multi = [0,2,2]
;tmp1 = (indgen(16)+1)*10		; *10 Myr
;tmp2 = (indgen(16)+1)*10 + 160	; *10 Myr
;t_Oden = where(abs(t_part1) ge 2950.-650. and abs(t_part1) le 2950.+110.)
;plot, r_part1[t_Oden], z_part1[t_Oden], xr=[rmin,rmax], yr=[zmin,zmax], xtitle='r [kpc]', ytitle='z [kpc]',/isotropic
;oplot, [r_3_20Gyr], [z_3_20Gyr], psym=7, color=255
;plot, y_part1[t_Oden], z_part1[t_Oden], xr=[ymax,ymin], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]',/isotropic
;oplot, [y_3_20Gyr], [z_3_20Gyr], psym=7, color=255
;plot, x_part1[t_Oden], z_part1[t_Oden], xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]',/isotropic
;oplot, [x_3_20Gyr], [z_3_20Gyr], psym=7, color=255
;plot, [0],[0],/nodata, xsty=4, ysty=4
;
;;=== the orbit of Pal5 on the meridional plane ===
r3D_part1 = sqrt(x_part1^2.+y_part1^2.+z_part1^2.)
peri = min(r3D_part1,sub_peri)
apo = max(r3D_part1,sub_apo)
print, 'peri=',peri, '   apo=',apo
;peri_apo = 'peri='+strtrim(string(peri,f='(f4.1)'),2)+'kpc apo='+strtrim(string(apo,f='(f4.1)'),2)+'kpc'
;help, sub_peri, sub_apo
;!p.multi = 0
;plot, r_part1, z_part1, psym=0, xr=[rmin,rmax], yr=[zmin,zmax], xtitle='r [kpc]', ytitle='z [kpc]', title=peri_apo, /isotropic
;xyouts, r_init, z_init, 't=0', charsize=1, color=255
;xyouts, r_final, z_final, t_mark_final, charsize=1, color=255
;oplot, [r_part1[tmp1[0]]], [z_part1[tmp1[0]]], psym=7, color=100
;oplot, r_part1[tmp1], z_part1[tmp1], psym=5, color=100
;oplot, r_part1[tmp2], z_part1[tmp2], psym=6, color=70
;oplot, [r_3_20Gyr], [z_3_20Gyr], psym=7, color=255
;erase
;
V_part1 = sqrt(Vx_part1^2.+Vy_part1^2.+Vz_part1^2.)
!p.multi=[0,1,2]
multiplot
plot, t_part1/1000., V_part1, xr=[0,3.5], yr=[80,250], ytitle=textoidl('V_{3D} [km/s]')
oplot, [t_part1[t_3_20Gyr]/1000.], [V_part1[t_3_20Gyr]], psym=7, color=255

multiplot
plot, t_part1/1000., r3D_part1, xr=[0,3.5],yr=[4.5,20.],xtitle='t [Gyr]', ytitle='r [kpc]'
oplot, [t_part1[t_3_20Gyr]/1000.], [r3D_part1[t_3_20Gyr]], psym=7, color=255
multiplot, /reset
erase

;=== the orbit of the central particle(part0001) for 2.95Gyr ===
!p.multi=[0,2,2]
plot, x_part1, y_part1, psym=0, xr=[xmin,xmax], yr=[ymin,ymax], xtitle='x [kpc]', ytitle='y [kpc]', title=Info1, /isotropic
oplot, [x_3_20Gyr], [y_3_20Gyr], psym=7, color=255
xyouts, x_init, y_init, 't=0', charsize=1, color=255
xyouts, x_1Gyr, y_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, x_2Gyr, y_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, x_500Myr, y_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, x_final, y_final, t_mark_final, charsize=1, color=255

plot, y_part1, z_part1, psym=0, xr=[ymin,ymax], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]', title=Info_posi, /isotropic
oplot, [y_3_20Gyr], [z_3_20Gyr], psym=7, color=255
xyouts, y_init, z_init, 't=0', charsize=1, color=255
xyouts, y_1Gyr, z_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, y_2Gyr, z_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, y_500Myr, z_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, y_final, z_final, t_mark_final, charsize=1, color=255

plot, x_part1, z_part1, psym=0, xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]', title=Info_vel, /isotropic
oplot, [x_3_20Gyr], [z_3_20Gyr], psym=7, color=255
xyouts, x_init, z_init, 't=0', charsize=1, color=255
xyouts, x_1Gyr, z_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, x_2Gyr, z_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, x_500Myr, z_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, x_final, z_final, t_mark_final, charsize=1, color=255

plot, r_part1, z_part1, psym=0, xr=[rmin,rmax], yr=[zmin,zmax], xtitle='r [kpc]', ytitle='z [kpc]', /isotropic
oplot, [r_3_20Gyr], [z_3_20Gyr], psym=7, color=255
xyouts, r_init, z_init, 't=0', charsize=1, color=255
xyouts, r_1Gyr, z_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, r_2Gyr, z_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, r_500Myr, z_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, r_final, z_final, t_mark_final, charsize=1, color=255


; read out the final position of all particles
;chr_rdtbl, 'snapshot/snap9990',0, arr
chr_rdtbl, 'snapshot/snap3200',0, arr
;chr_rdtbl, 'snapshot/snap3350',0, arr
;chr_rdtbl, 'snapshot/snap2950',0, arr
;chr_rdtbl, 'snapshot_dt_0.1/snap2950',0,arr
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
Vr = reform( (Vx*(x-x_sun)+Vy*y+Vz*z)/sqrt((x-x_sun)^2.+y^2.+z^2.) )
E_radial = Vr^2./2.

xmin =   6  &  xmax = 10
ymin =  -3  &  ymax =  3
zmin =  14  &  zmax = 17.5
Vxmin =  -60 &  Vxmax = -10
Vymin = -100 &  Vymax = -80
Vzmin =  -60 &  Vzmax =  20

;=== final position of particles: the shape of tidal tail ===
!p.multi=[0,2,2]
plot, x, y, psym=3, xr=[xmin,xmax], yr=[ymin,ymax], xtitle='x [kpc]', ytitle='y [kpc]', title='Pal 5', /isotropic
oplot, [x[0]], [y[0]], psym=1, color=255
plot, y, z, psym=3, xr=[ymin,ymax], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]', /isotropic
oplot, [y[0]], [z[0]], psym=1, color=255
plot, x, z, psym=3, xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]', /isotropic
oplot, [x[0]], [z[0]], psym=1, color=255

l = atan(y/(x+8.))*!radeg 
l_sub = where(l lt 0)
;if l_sub[0] ne -1 then l[l_sub] = l[l_sub] + 360.
b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
lcosb = l*cos(b*!dtor)


;=== estimate equatorial coord. ===
delta = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
alpha = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(delta*!dtor) )*!radeg + 282.25
alpha = reform(alpha)
delta = reform(delta)
alpha_sub = where(alpha lt 0)
if alpha_sub[0] ne -1 then alpha[alpha_sub] = alpha[alpha_sub] + 360.
alpha_sub2 = where(alpha gt 360)
if alpha_sub2[0] ne -1 then alpha[alpha_sub2] = alpha[alpha_sub2] - 360.

;delta=randomn(iseed,N_elements(alpha))/5+((alpha-228.5)/3)^5.-((alpha-230.5)/3)^3.

sub_t = where(t_part1 gt 3050 and t_part1 lt 3350)
;sub_t = where(t_part1 gt 3200 and t_part1 lt 3500)
;sub_t = where(t_part1 gt 2850 and t_part1 lt 3050)
l_part1 = atan(y_part1/(x_part1+8.))*!radeg 
l_part1_sub = where(l_part1 lt 0)
;if l_part1_sub[0] ne -1 then l_part1[l_part1_sub] = l_part1[l_part1_sub] + 360.
l_part1_sub2 = where(l_part1 gt 360)
;if l_part1_sub2[0] ne -1 then l_part1[l_part1_sub2] = l_part1[l_part1_sub2] - 360.
b_part1 = atan(z_part1/sqrt((x_part1+8.)^2.+y_part1^2.))*!radeg
lcosb_part1 = l_part1*cos(b_part1*!dtor)
;sub_b = where(lcosb_part1 gt -12 and lcosb_part1 lt 12 and t_part1 gt 2800 and t_part1 lt 3100)


delta_part1 = asin( cos(b_part1*!dtor)*sin((l_part1-33.)*!dtor)*sin(62.6*!dtor)+sin(b_part1*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_part1 = asin( (cos(b_part1*!dtor)*sin((l_part1-33.)*!dtor)*cos(62.6*!dtor)-sin(b_part1*!dtor)*sin(62.6*!dtor))/cos(delta_part1*!dtor) )*!radeg + 282.25
alpha_part1 = reform(alpha_part1)
delta_part1 = reform(delta_part1)
alpha_part1_sub = where(alpha_part1 lt 0)
if alpha_part1_sub[0] ne -1 then alpha_part1[alpha_part1_sub] = alpha_part1[alpha_part1_sub] + 360.
alpha_part1_sub2 = where(alpha_part1 gt 360)
if alpha_part1_sub2[0] ne -1 then alpha_part1[alpha_part1_sub2] = alpha_part1[alpha_part1_sub2] - 360.
sub_b = sub_t
print, 'max(lcosb)=',max(lcosb),'  min(lcosb)=',min(lcosb)
alpha_Pal5 = alpha[0]
delta_Pal5 = delta[0]
l_Pal5 = l[0]
lcosb_Pal5 = lcosb[0]
b_Pal5 = b[0]


;=== Pal5 in equatorial coord. ===
!p.multi=0
plot,[0],[0],/nodata,xr=[247,218],yr=[-10,10],xtitle=textoidl('\alpha [deg]'),ytitle=textoidl('\delta [deg]')
;plot,[0],[0],/nodata,xr=[247,218],yr=[-22,18],xtitle=textoidl('\alpha [deg]'),ytitle=textoidl('\delta [deg]'),title=sigma
oplot,alpha,delta,psym=8,symsize=0.2
;oplot,alpha_part1,delta_part1,color=100;,psym=1
oplot,alpha_part1[sub_t],delta_part1[sub_t],color=100;,psym=1
oplot,[alpha_Pal5],[delta_Pal5],psym=1,color=255,symsize=1.8,thick=8

;;=== polynomial fit of orbit ===
;north = where(alpha ge alpha_Pal5,complement=south)
;fit = poly_fit(alpha[north],delta[north],5)
;xfit = findgen(1000)/30.+220
;yfit = fit[0]+fit[1]*xfit+fit[2]*xfit^2.+fit[3]*xfit^3.+fit[4]*xfit^4.+fit[5]*xfit^5.
;oplot,xfit,yfit,color=200
;fit = poly_fit(alpha[south],delta[south],5)
;xfit = findgen(1000)/30.+220
;yfit = fit[0]+fit[1]*xfit+fit[2]*xfit^2.+fit[3]*xfit^3.+fit[4]*xfit^4.+fit[5]*xfit^5.
;oplot,xfit,yfit,color=220

N_particle = N_elements(alpha)
goto, skip
;=== heating parameter ===
sum_Ediff = 0.d
sum_Eradial_diff = 0.d
Ediff = dblarr(N_particle)
Eradial_diff = dblarr(N_particle)
for i=1,N_particle-1 do begin	; the first particle should be excluded.
  d = sqrt((x[i]-x)^2.+(y[i]-y)^2.+(z[i]-z)^2.)
  dsort = sort(d)
  Ediff[i] = E_total[i]-E_total[dsort[1]]
  sum_Ediff += Ediff[i]^2.

  d_2d = sqrt((alpha[i]-alpha)^2.+(delta[i]-delta)^2.)
  dsort_2d = sort(d_2d)
  Eradial_diff[i] = E_radial[i]-E_radial[dsort_2d[1]]
  sum_Eradial_diff += Eradial_diff[i]^2.
endfor
P_heat = sqrt(sum_Ediff/N_particle)
P_heat_radial = sqrt(sum_Eradial_diff/N_particle)


;=== kink parameter ===
sub = round(randomu(32,30)*(N_particle-1))
N_sub = N_elements(sub)
alpha_test = alpha[sub]
delta_test = delta[sub]
theta3 = 3.d
theta6 = 6.d
theta10 = 10.d
kink3 = dblarr(N_sub)
kinkv3 = dblarr(N_sub)
for i=0,N_sub-1 do begin
  d_2d = sqrt((alpha_test[i]-alpha)^2.+(delta_test[i]-delta)^2.)
  in_circle = where(d_2d le theta3)
  N_in = N_elements(in_circle)
  alpha_in = alpha[in_circle]
  delta_in = delta[in_circle]
  Vr_in = Vr[in_circle]
  fit = poly_fit(alpha_in,delta_in,1)
  fitv = poly_fit(alpha_in,Vr_in,1)
  d2line = abs(fit[1]*alpha_in-1.*delta_in+fit[0])/sqrt(fit[1]^2.+1.)
  d2linev = abs(fit[1]*alpha_in-1.*Vr_in+fit[0])/sqrt(fit[1]^2.+1.)
  sigma_y = stddev(d2line)
  sigma_v = stddev(d2linev)
  sum_del_y2 = 0.d
  sum_del_v2 = 0.d
  for j=0,N_in-1 do begin
    d2ortho_line = abs(1./fit[1]*alpha_in+1.*delta_in-alpha_in[j]/fit[1]-delta_in[j])/sqrt((1./fit[1])^2.+1.)
    del_y = d2line[j]-d2line[(sort(d2ortho_line))[1]]
    sum_del_y2 += del_y^2. 
    d2ortho_linev = abs(1./fit[1]*alpha_in+1.*Vr_in-alpha_in[j]/fit[1]-Vr_in[j])/sqrt((1./fit[1])^2.+1.)
    del_v = d2linev[j]-d2linev[(sort(d2ortho_linev))[1]]
    sum_del_v2 += del_v^2. 
  endfor
  mean_del_y2 = sum_del_y2/N_in
  mean_del_v2 = sum_del_v2/N_in
  kink3[i] = (2.*sigma_y^2./mean_del_y2-1.)
  kinkv3[i] = (2.*sigma_v^2./mean_del_v2-1.)
endfor

kink6 = dblarr(N_sub)
for i=0,N_sub-1 do begin
  d_2d = sqrt((alpha_test[i]-alpha)^2.+(delta_test[i]-delta)^2.)
  in_circle = where(d_2d le theta6)
  N_in = N_elements(in_circle)
  alpha_in = alpha[in_circle]
  delta_in = delta[in_circle]
  fit = poly_fit(alpha_in,delta_in,1)
  d2line = abs(fit[1]*alpha_in-1.*delta_in+fit[0])/sqrt(fit[1]^2.+1.)
  sigma_y = stddev(d2line)
  sum_del_y2 = 0.d
  for j=0,N_in-1 do begin
    d2ortho_line = abs(1./fit[1]*alpha_in+1.*delta_in-alpha_in[j]/fit[1]-delta_in[j])/sqrt((1./fit[1])^2.+1.)
    del_y = d2line[j]-d2line[(sort(d2ortho_line))[1]]
    sum_del_y2 += del_y^2. 
  endfor
  mean_del_y2 = sum_del_y2/N_in
  kink6[i] = (2.*sigma_y^2./mean_del_y2-1.)
endfor

kink10 = dblarr(N_sub)
for i=0,N_sub-1 do begin
  d_2d = sqrt((alpha_test[i]-alpha)^2.+(delta_test[i]-delta)^2.)
  in_circle = where(d_2d le theta10)
  N_in = N_elements(in_circle)
  alpha_in = alpha[in_circle]
  delta_in = delta[in_circle]
  fit = poly_fit(alpha_in,delta_in,1)
  d2line = abs(fit[1]*alpha_in-1.*delta_in+fit[0])/sqrt(fit[1]^2.+1.)
  sigma_y = stddev(d2line)
  sum_del_y2 = 0.d
  for j=0,N_in-1 do begin
    d2ortho_line = abs(1./fit[1]*alpha_in+1.*delta_in-alpha_in[j]/fit[1]-delta_in[j])/sqrt((1./fit[1])^2.+1.)
    del_y = d2line[j]-d2line[(sort(d2ortho_line))[1]]
    sum_del_y2 += del_y^2. 
  endfor
  mean_del_y2 = sum_del_y2/N_in
  kink10[i] = (2.*sigma_y^2./mean_del_y2-1.)
endfor

m_kink3 = mean(kink3)
m_kink6 = mean(kink6)
m_kink10 = mean(kink10)
s_kink3 = reverse(sort(kink3))
s_kink6 = reverse(sort(kink6))
s_kink10 = reverse(sort(kink10))
print, 'kink',kink3
max1_kink3 = kink3[s_kink3[0]]
max2_kink3 = kink3[s_kink3[1]]
max3_kink3 = kink3[s_kink3[2]]
max1_kink6 = kink6[s_kink6[0]]
max2_kink6 = kink6[s_kink6[1]]
max3_kink6 = kink6[s_kink6[2]]
max1_kink10 = kink10[s_kink10[0]]
max2_kink10 = kink10[s_kink10[1]]
max3_kink10 = kink10[s_kink10[2]]
legend,[textoidl('P_{heat}=')+strtrim(P_heat,2),$
        textoidl('P_{heat,K(r)}=')+strtrim(P_heat_radial,2),$
        textoidl('Kink_{\theta=3}=')+strtrim(m_kink3,2),$
        textoidl('Kink1_{\theta=3}=')+strtrim(max1_kink3,2),$
        textoidl('Kink2_{\theta=3}=')+strtrim(max2_kink3,2),$
        textoidl('Kink3_{\theta=3}=')+strtrim(max3_kink3,2),$
        textoidl('Kink_{\theta=6}=')+strtrim(m_kink6,2),$
        textoidl('Kink1_{\theta=6}=')+strtrim(max1_kink6,2),$
        textoidl('Kink2_{\theta=6}=')+strtrim(max2_kink6,2),$
        textoidl('Kink3_{\theta=6}=')+strtrim(max3_kink6,2),$
        textoidl('Kink_{\theta=10}=')+strtrim(m_kink10,2),$
        textoidl('Kink1_{\theta=10}=')+strtrim(max1_kink10,2),$
        textoidl('Kink2_{\theta=10}=')+strtrim(max2_kink10,2),$
        textoidl('Kink3_{\theta=10}=')+strtrim(max3_kink10,2)],box=0,/bottom

s_alpha = sort(alpha_test)
plot,[0],[0],/nodata,xr=[max(alpha_test)*1.01,min(alpha_test)*0.99],xtitle=textoidl('\alpha'),ystyle=4
axis,yaxis=0,yr=[min(kink3)-0.1,max(kink3)*1.1],ytitle=textoidl('Kink_y'),/save,color=70
oplot,alpha_test[s_alpha],kink3[s_alpha],psym=0,color=70
axis,yaxis=1,yr=[min(kinkv3)-1,max(kinkv3)*1.1],ytitle=textoidl('Kink_v'),/save,color=255
oplot,alpha_test[s_alpha],kinkv3[s_alpha],psym=0,color=255


goto, skip
;=== kink parameter test ===
!p.multi=[0,4,5]
m_kink = dblarr(20)
all_kink = dblarr(20,9999)
all_kink[*,*] = -9999
all_theta = dblarr(20)
for p=0,19 do begin
  delta=randomn(iseed,N_elements(alpha))/5+((alpha-228.5)/3)^5.-((alpha-230.5)/3)^3.
  sub20 = where(delta lt 20 and delta gt -20)  
  sub = sub20[round(randomu(32,30)*(N_elements(sub20)-1))]
  N_sub = N_elements(sub)
  alpha_test = alpha[sub]
  delta_test = delta[sub]
  theta = 1.d*p+1.
  kink = dblarr(N_sub)
  for i=0,N_sub-1 do begin
    d_2d = sqrt((alpha_test[i]-alpha)^2.+(delta_test[i]-delta)^2.)
    in_circle = where(d_2d le theta)
    N_in = N_elements(in_circle)
    alpha_in = alpha[in_circle]
    delta_in = delta[in_circle]
    fit = poly_fit(alpha_in,delta_in,1)
    d2line = abs(fit[1]*alpha_in-1.*delta_in+fit[0])/sqrt(fit[1]^2.+1.)
    sigma_y = stddev(d2line)
    sum_del_y2 = 0.d
    for j=0,N_in-1 do begin
      d2ortho_line = abs(1./fit[1]*alpha_in+1.*delta_in-alpha_in[j]/fit[1]-delta_in[j])/sqrt((1./fit[1])^2.+1.)
      del_y = d2line[j]-d2line[(sort(d2ortho_line))[1]]
      sum_del_y2 += del_y^2. 
    endfor
    mean_del_y2 = sum_del_y2/N_in
    kink[i] = (2.*sigma_y^2./mean_del_y2-1.)
  endfor
  plot,alpha_test,kink,psym=8,xr=[max(alpha_test)*1.01,min(alpha_test)*0.99],$
  ;    yr=[-0.3,0.3],xtitle=textoidl('\alpha'),ytitle='Kink'
      yr=[min(kink)-0.1,max(kink)+0.1],xtitle=textoidl('\alpha'),ytitle='Kink'
  hline,0,linestyle=2
  m_kink[p] = mean(kink)
  all_theta[p] = theta
  for q=0,N_sub-1 do all_kink[p,q] = kink[q]
endfor
!p.multi=0
plot,[0],[0],/nodata,xr=[0,21],yr=[-1,max(all_kink)*1.1],xtitle=textoidl('\theta'),ytitle='Kink'
for p=0,19 do begin
  sub = where(all_kink[p,*] ge -999)
  for q=0,N_elements(sub)-1 do oplot,[all_theta[p]],[all_kink[p,q]],psym=7
endfor
oplot,all_theta,all_kink,psym=8,color=255
skip:


;=== clumpy parameter ===
min_d2 = dblarr(N_particle)
for i=0,N_particle-1 do begin
  d2 = sqrt((alpha[i]-alpha)^2.+(delta[i]-delta)^2.)
  s_d2 = sort(d2)
  min_d2[i] = d2[s_d2[5]]
endfor
plot,alpha,1./min_d2,psym=1,/ylog
legend,['max/median='+strtrim(max(1./min_d2)/median(1./min_d2,/even),2)],box=0,/bottom


;=== Pal5 in lcosb & b plane ===
if (N_subhalos eq 0) then begin 
  plot, lcosb, b, psym=3, xr=[8,-4], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title=sigma
endif else begin
  plot, lcosb,b,xr=[max(lcosb),min(lcosb)],yr=[min(b),max(b)],psym=3,xtitle='l cos b [deg]',ytitle='b [deg]',title=sigma
endelse
oplot, [lcosb_Pal5], [b_Pal5], psym=1, color=255,symsize=1.8,thick=8

;=== sigma clipping fitting of test particles ===
;yerror0 = 0
;yerror=1
;lcosb_fit = lcosb
;b_fit = b
;while(yerror0 ne yerror) do begin
;  yerror0 = yerror
;  fit = poly_fit(lcosb_fit,b_fit,4,yerror=yerror,yfit=yfit)
;  sigma_clipping = where(b_fit ge (yfit-3.*yerror) and b_fit le (yfit+3.*yerror))
;  lcosb_fit = lcosb_fit[sigma_clipping]
;  b_fit = b_fit[sigma_clipping]
;endwhile

;yfit = fit[0]+fit[1]*lcosb+fit[2]*lcosb^2.+fit[3]*lcosb^3.+fit[4]*lcosb^4.
;lcosb_sub = sort(lcosb)
;oplot,lcosb[lcosb_sub],yfit[lcosb_sub],color=255
;rms = sqrt(total((b-yfit)^2.)/N_elements(b))

oplot,lcosb_part1[sub_t],b_part1[sub_t],color=70;,psym=1
b_interpol = interpol(b_part1[sub_t],lcosb_part1[sub_t],lcosb,/spline)	; interpolation of orbit
orbit_rms = sqrt(total((b_interpol-b)^2.)/N_elements(b))

;orbit_fit = poly_fit(lcosb_part1[sub_t],b_part1[sub_t],6,yfit=part1_orbit_yfit)
;orbit_yfit = orbit_fit[0]+orbit_fit[1]*lcosb[sub_t]+orbit_fit[2]*lcosb[sub_t]^2.+orbit_fit[3]*lcosb[sub_t]^3.+orbit_fit[4]*lcosb[sub_t]^4.+orbit_fit[5]*lcosb[sub_t]^5.+orbit_fit[6]*lcosb[sub_t]^6.
;orbit_rms = sqrt(total((b[sub_t]-orbit_yfit)^2.)/N_elements(sub_t))
;oplot,lcosb_part1[sub_t],part1_orbit_yfit,color=70
;legend, ['rms_orbit='+strtrim(string(orbit_rms,f='(f12.8)'),2),'rms='+strtrim(string(rms,f='(f12.8)'),2)],box=0
;legend, ['  orbit  ','orbit fit','   fit   '],color=[100,70,255],psym=[1,0,0],box=0,/right
legend, ['rms from orbit='+strtrim(string(orbit_rms,f='(f12.8)'),2)],box=0
legend, ['Pal5 center','  orbit   '],color=[255,70],psym=[1,0],box=0,/right

;=== contour plot of Pal5 ===
contour_bin = 30
levels = indgen(300)*contour_bin + contour_bin
pts_contour,lcosb,b,20,15,gv,x,y
if (N_subhalos eq 0) then begin 
  contour,gv,x,y,xr=[8,-4],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1
endif else begin
  contour,gv,x,y,xr=[max(lcosb),min(lcosb)],yr=[min(b),max(b)],xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1
endelse


;=== Pal5 in lcosb & b plane 2 ===
plot, lcosb, b, psym=3, xr=[8,-4], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title=sigma
oplot, [lcosb_Pal5], [b_Pal5], psym=1, color=255,symsize=1.8,thick=8

oplot,lcosb_part1[sub_t],b_part1[sub_t],color=70;,psym=1
legend, ['rms from orbit='+strtrim(string(orbit_rms,f='(f12.8)'),2)],box=0
legend, ['Pal5 center','  orbit   '],color=[255,70],psym=[1,0],box=0,/right
;=== contour plot of Pal5 ===
contour,gv,x,y,xr=[8,-4],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1


;=== density profile of Pal 5 ===
radius = 60.*sqrt((lcosb-lcosb_Pal5)^2.+(b-b_Pal5)^2.)
dr = 5.
annuli = findgen(1000)*dr
tmp_N = 0
N_annuli = N_elements(annuli)
cul_Ndr = intarr(N_annuli)
Ndr = intarr(N_annuli)
for i=0,N_annuli-1 do begin
  in_annuli = where(radius ge annuli[i] and radius lt annuli[i]+dr)
  Ndr[i] = (in_annuli[0] eq -1) ? 0:N_elements(in_annuli)
  tmp_N = tmp_N + Ndr[i]
  cul_Ndr[i] = tmp_N
endfor
surface_density = Ndr/(!pi*(annuli+dr)^2.)
surface_density_err = sqrt(Ndr)/(!pi*(annuli+dr)^2.)


;=== width profile ===
dw = 5./60. 	; 1arcmin
width = findgen(1000)*dw
N_w = N_elements(width)
Ndw = intarr(N_w)
for i=0,N_w-1 do begin
  in_dw = where(abs(b-b_interpol) ge width[i] and abs(b-b_interpol) lt width[i]+dw)
  Ndw[i] = (in_dw[0] eq -1) ? 0:N_elements(in_dw)
endfor

erase
!p.multi=[0,1,3]
;multiplot
ploterror,annuli+dr/2.,surface_density,surface_density_err, $
	xr=[0.6,max(annuli)*1.05],yr=[1e-4,2],ytitle=textoidl('\Sigma [arcmin^{-2}]'), /xlog, /ylog
;multiplot
ploterror, annuli+dr/2., cul_Ndr, sqrt(cul_Ndr),xr=[0.6,max(annuli)*1.05],xtitle='r [arcmin]',ytitle='N(<r)',/xlog
;multiplot
ploterror, width+dw/2.,Ndw,sqrt(Ndw),xr=[dw/2.*0.95,10],yr=[0.8,max(Ndw)*1.2],xtitle='width [deg]',ytitle='Ndw',/xlog,/ylog
;multiplot, /reset
erase

;Vtot = (Vx^2.+Vy^2.+Vz^2.)
;Vtot_part1 = (Vx_part1^2.+Vy_part1^2.+Vz_part1^2.)
;sub = where( Etot-Etot_part1[295] gt 0)
;oplot, lcosb[sub],b[sub],psym=1
;ub = where( Vtot-Vtot_part1[295] gt 0)
;oplot, lcosb[sub],b[sub],psym=1
;plot, x, Vx, psym=3, xr=[xmin,xmax], yr=[Vxmin,Vxmax], xtitle=textoidl('x [kpc]'), ytitle=textoidl('V_x [km/s]'), title='Pal 5'
;oplot, [x[0]], [Vx[0]], psym=1, color=255


;=== final position of particles: the distribution of particles ===
!p.multi=[0,2,2]
sig = stddev(lcosb)
bsize = sig/5.
plothist, lcosb, xh, yh, bin=bsize, xtitle='l cos b [deg]', ytitle='#'
oploterror, xh, yh, sqrt(yh), psym=3
oplot, [0,sig],[0,0],thick=12
legend, [textoidl('\sigma=')+strtrim(string(sig,f='(f12.8)'),2)],box=0

bsize = stddev(delE_part1)/5.
plothist, delE_part1, bin=bsize, xtitle=textoidl('\DeltaE / E (part1)'), ytitle='#'

bsize = stddev(delJ_part1)/5.
plothist, delJ_part1, bin=bsize, xtitle=textoidl('\DeltaJ / J (part1)'), ytitle='#'

plot,[0],[0],/nodata,xsty=4,ysty=4

;=== Energy & Angular momentum check ===
;Mvir = 1.4e+12
;rs = 20.7406d
;rvir = 287.7606d
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d

c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl, 'part_peri',0, arr
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
V_peri = sqrt(Vx_peri^2.+Vy_peri^2.+Vz_peri^2.)
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
r_tide = 0.108653
Msat = 20000.
epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)


;====================================
; check E&J when there is no subhalo 
;====================================
  ;=== Energy plot ===
  delE = E_total-E_total[0]
  q = (delE)/epsilon
  bsize = stddev(delE)/5.
  plothist, delE, bin=bsize,xtitle=textoidl('\DeltaE (E_{debris}-E_{sat})')
  sig = stddev(delE)
  in_delE = where(delE gt -1.*sig and delE lt sig,complement=out_delE)
  plothist, delE[in_delE], bin=bsize, /fill, /overplot, fcolor=100
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(round(sig),2),alignment=0.5
  dJ = J-J[0]
  bsize = stddev(dJ)/5.
  plothist, dJ, bin=bsize, xtitle=textoidl('\DeltaJ (J_{debris}-J_{sat})')
  sig = stddev(dJ)
  in_dJ = where(dJ gt -1.*sig and dJ lt sig,complement=out_dJ)
  plothist, dJ[in_dJ], bin=bsize, /fill, /overplot, fcolor=100
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f14.8)'),2),alignment=0.5
  
  scale_dJ = dJ / (r_tide/r_peri * J)
;  scale_dJ_sub1 = where(abs(scale_dJ) le 0.5,complement=scale_dJ_out)
;  scale_dJ_sub2 = where(abs(scale_dJ) gt 0.5 and abs(scale_dJ) le 1.)
;  scale_dJ_sub3 = where(abs(scale_dJ) gt 1. and abs(scale_dJ) le 1.5)
;  scale_dJ_sub4 = where(abs(scale_dJ) gt 1.5)
;  scale_dJ_sub1_p = where(scale_dJ gt 0 and scale_dJ le 0.5,complement=scale_dJ_out)
;  scale_dJ_sub2_p = where(scale_dJ gt 0.5 and scale_dJ le 1.)
;  scale_dJ_sub3_p = where(scale_dJ gt 1. and scale_dJ le 1.5)
;  scale_dJ_sub4_p = where(scale_dJ gt 1.5)
;  scale_dJ_sub1_m = where(scale_dJ gt -0.5 and scale_dJ le 0,complement=scale_dJ_out)
;  scale_dJ_sub2_m = where(scale_dJ le -0.5 and scale_dJ gt -1.)
;  scale_dJ_sub3_m = where(scale_dJ le -1. and scale_dJ gt -1.5)
;  scale_dJ_sub4_m = where(scale_dJ le -1.5)
  scale_dJ_sub1 = where(abs(scale_dJ) le 0.25,complement=scale_dJ_out)
  scale_dJ_sub2 = where(abs(scale_dJ) gt 0.25 and abs(scale_dJ) le 0.5)
  scale_dJ_sub3 = where(abs(scale_dJ) gt 0.5 and abs(scale_dJ) le 1.)
  scale_dJ_sub4 = where(abs(scale_dJ) gt 1.)
  scale_dJ_sub1_p = where(scale_dJ gt 0 and scale_dJ le 0.25,complement=scale_dJ_out)
  scale_dJ_sub2_p = where(scale_dJ gt 0.25 and scale_dJ le 0.5)
  scale_dJ_sub3_p = where(scale_dJ gt 0.5 and scale_dJ le 1.)
  scale_dJ_sub4_p = where(scale_dJ gt 1.)
  scale_dJ_sub1_m = where(scale_dJ gt -0.25 and scale_dJ le 0,complement=scale_dJ_out)
  scale_dJ_sub2_m = where(scale_dJ le -0.25 and scale_dJ gt -0.5)
  scale_dJ_sub3_m = where(scale_dJ le -0.5 and scale_dJ gt -1.)
  scale_dJ_sub4_m = where(scale_dJ le -1.)
  
  bsize = stddev(epsilon)/5.
  plothist, epsilon, bin=bsize,xtitle=textoidl('\epsilon')
  q_sub1 = where(abs(q) le 1,complement=q_out)
  q_sub2 = where(abs(q) gt 1 and abs(q) le 2)
  q_sub3 = where(abs(q) gt 2 and abs(q) le 3)
  q_sub4 = where(abs(q) gt 3)
;  q_sub1 = where(abs(q) le 0.5,complement=q_out)
;  q_sub2 = where(abs(q) gt 0.5 and abs(q) le 1)
;  q_sub3 = where(abs(q) gt 1 and abs(q) le 1.5)
;  q_sub4 = where(abs(q) gt 1.5)
  bsize = stddev(q)/5.
;  plot,[0],[0],/nodata,xsty=4,ysty=4

if N_subhalos eq 0 then begin
  plothist, q, bin=bsize, xtitle='q'
  plothist, q[q_sub1],bin=bsize,/fill,/overplot,fcolor=255
  plothist, q[q_sub2],bin=bsize,/fill,/overplot,fcolor=220
  plothist, q[q_sub3],bin=bsize,/fill,/overplot,fcolor=100
  plothist, q[q_sub4],bin=bsize,/fill,/overplot,fcolor=70
  sig = stddev(q)
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f12.8)'),2),alignment=0.5
;  plot,[0],[0],/nodata,xsty=4,ysty=4

  
;=== Angular momentum plot ===
  bsize = stddev(scale_dJ)/5.
  plothist, scale_dJ, bin=bsize, xtitle=textoidl('\DeltaJ / sJ (s=r_{tide}/R_{peri})')
  plothist, scale_dJ[scale_dJ_sub1],bin=bsize,/fill,/overplot,fcolor=255
  plothist, scale_dJ[scale_dJ_sub2],bin=bsize,/fill,/overplot,fcolor=220
  plothist, scale_dJ[scale_dJ_sub3],bin=bsize,/fill,/overplot,fcolor=100
  plothist, scale_dJ[scale_dJ_sub4],bin=bsize,/fill,/overplot,fcolor=70
  sig = stddev(scale_dJ)
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f14.8)'),2),alignment=0.5
  
  
;  ;=== z-Angular momentum plot ===
  cut_sub = 3.5
  lcosb_sub = where(lcosb gt cut_sub or lcosb lt -1.*cut_sub)
  bsize = stddev(delE)/5.
  cut_title = string(cut_sub,f='(f3.1)')
  plothist, delE[lcosb_sub], bin=bsize,xtitle=textoidl('\DeltaE (E_{debris}-E_{sat})'),title='lcosb<-'+cut_title+' or lcosb>'+cut_title
  plothist, delE[q_out], bin=bsize,xtitle=textoidl('\DeltaE (E_{debris}-E_{sat})'), title='|q|>1'
  bsize = stddev(q)/5.
  plothist, q[lcosb_sub], bin=bsize, xtitle='q',title='lcosb<-'+cut_title+' or lcosb>'+cut_title
  
  
  ;=== test particles with different energy ===
  plot, lcosb, b, psym=3, xr=[12,-10], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title=sigma
  oplot, [lcosb[0]], [b[0]], psym=1, color=255, thick=8
  fit = poly_fit(lcosb_part1[sub_b],b_part1[sub_b],2,yerror=yerror,yfit=yfit)
  y = fit[0] + fit[1]*lcosb + fit[2]*lcosb^2.
  oplot, lcosb_part1[sub_t], b_part1[sub_t], color=70
  rms = sqrt(total((b-y)^2.)/N_elements(b))
  legend, ['rms='+strtrim(string(rms,f='(f12.8)'),2)], box=0
  
  plot, lcosb[in_delE], b[in_delE], psym=3, xr=[12,-10], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title='in dE'
  oplot, [lcosb[0]], [b[0]], psym=1, color=255, thick=8
  
  plot, lcosb[out_delE], b[out_delE], psym=3, xr=[12,-10], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title='out dE'
  oplot, [lcosb[0]], [b[0]], psym=1, color=255, thick=8
  
  ;plot,[0],[0],/nodata,xr=[60,-5],yr=[-80,40],xtitle='l cos b [deg]',ytitle='b [deg]'
  plot, [0],[0],/nodata,xr=[12,-10],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]'
  oplot, lcosb[q_sub1],b[q_sub1],psym=8,symsize=0.2,color=255
  oplot, lcosb[q_sub2],b[q_sub2],psym=8,symsize=0.2,color=220
  oplot, lcosb[q_sub3],b[q_sub3],psym=8,symsize=0.2,color=100
  oplot, lcosb[q_sub4],b[q_sub4],psym=8,symsize=0.2,color=70
  oplot, [lcosb[0]], [b[0]], psym=1, thick=8
  legend, ['|q|<1','1<|q|<2','2<|q|<3','|q|>3'],psym=8,color=[255,220,100,70],box=0
  
  
  ;=== test particles with different angular momentum ===
  ;plot,[0],[0],/nodata,xr=[60,-5],yr=[-80,40],xtitle='l cos b [deg]',ytitle='b [deg]'
  plot,[0],[0],/nodata,xr=[12,-10],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]'
  oplot,lcosb[scale_dJ_sub1_p],b[scale_dJ_sub1_p],psym=1,symsize=0.5,color=255
  oplot,lcosb[scale_dJ_sub2_p],b[scale_dJ_sub2_p],psym=1,symsize=0.5,color=220
  oplot,lcosb[scale_dJ_sub3_p],b[scale_dJ_sub3_p],psym=1,symsize=0.5,color=100
  oplot,lcosb[scale_dJ_sub4_p],b[scale_dJ_sub4_p],psym=1,symsize=0.5,color=70
  oplot,[lcosb[0]],[b[0]],psym=1,thick=8
  
  ;plot,[0],[0],/nodata,xr=[60,-5],yr=[-80,40],xtitle='l cos b [deg]',ytitle='b [deg]'
  plot,[0],[0],/nodata,xr=[12,-10],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]'
  oplot,lcosb[scale_dJ_sub1_m],b[scale_dJ_sub1_m],psym=7,symsize=0.5,color=255
  oplot,lcosb[scale_dJ_sub2_m],b[scale_dJ_sub2_m],psym=7,symsize=0.5,color=220
  oplot,lcosb[scale_dJ_sub3_m],b[scale_dJ_sub3_m],psym=7,symsize=0.5,color=100
  oplot,lcosb[scale_dJ_sub4_m],b[scale_dJ_sub4_m],psym=7,symsize=0.5,color=70
  oplot,[lcosb[0]],[b[0]],psym=1,thick=8
  
  ;plot,[0],[0],/nodata,xr=[60,-5],yr=[-80,40],xtitle='l cos b [deg]',ytitle='b [deg]'
  plot,[0],[0],/nodata,xr=[12,-10],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]'
  oplot,lcosb[scale_dJ_sub1],b[scale_dJ_sub1],psym=8,symsize=0.2,color=255
  oplot,lcosb[scale_dJ_sub2],b[scale_dJ_sub2],psym=8,symsize=0.2,color=220
  oplot,lcosb[scale_dJ_sub3],b[scale_dJ_sub3],psym=8,symsize=0.2,color=100
  oplot,lcosb[scale_dJ_sub4],b[scale_dJ_sub4],psym=8,symsize=0.2,color=70
  oplot,[lcosb[0]],[b[0]],psym=1,thick=8
  legend, [textoidl('|\DeltaJ/sJ|<0.25'),textoidl('0.25<|\DeltaJ/sJ|<0.5'),textoidl('0.5<|\DeltaJ/sJ|<1'),textoidl('|\DeltaJ/sJ|>1')],psym=8,color=[255,220,100,70],box=0
;  legend, [textoidl('|\DeltaJ/sJ|<0.5'),textoidl('0.5<|\DeltaJ/sJ|<1'),textoidl('1<|\DeltaJ/sJ|<1.5'),textoidl('|\DeltaJ/sJ|>1.5')],psym=8,color=[255,220,100,70],box=0
  
  ;plot,[0],[0],/nodata,xsty=4,ysty=4
  plot,lcosb,b,xr=[max(lcosb),min(lcosb)],yr=[min(b),max(b)],xtitle='l cos b [deg]',ytitle='b [deg]',psym=3
  oplot,[lcosb[0]],[b[0]],psym=1,thick=8
  ;bottom_sub = where(b lt 45.5 and lcosb gt 0 and lcosb lt 1)
  ;oplot,lcosb[bottom_sub],b[bottom_sub],psym=7
  
endif
  ;=== comparison between energy & angular momentum ===
  plot,delE,dJ,psym=3,xtitle=textoidl('\DeltaE'),ytitle=textoidl('\DeltaJ')
  vline,0,linestyle=2
  hline,0,linestyle=2
  plot,q,scale_dJ,psym=3,xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  vline,0,linestyle=2
  hline,0,linestyle=2
  ;oplot,q[bottom_sub],scale_dJ[bottom_sub],psym=8,symsize=0.2,color=255
  
  ; region of q > 1
  plot,delE[q_out],dJ[q_out],psym=3,xtitle=textoidl('\DeltaE'),ytitle=textoidl('\DeltaJ')
  vline,0,linestyle=2
  hline,0,linestyle=2
  plot,q[q_out],scale_dJ[q_out],psym=3,xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  vline,0,linestyle=2
  hline,0,linestyle=2


;;=== radial mass distribution ===
;;!p.multi=0
r = findgen(10001)/10. + 1./10.
p = r/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
V = 3./4.*!pi*r^3.
;openr,6,'R_peri'
;readf,6,R_peri
;close,6
readcol,'R_peri',R_peri
print,R_peri
tmp = min(abs(r-R_peri[0]),min_sub)
print, 'M(R<R_peri)=',r[min_sub],mr[min_sub]
print, 'M(R<100kpc)=',mr[where(r eq 100)]
openw,11,'M_Rperi'
printf,11,mr[min_sub]
close,11


;=== tail growth rate -- Johnston(2001) ===
r_peri = 5.50114d
p_analytic = r_peri/rs
mr_analytic = Mhalo*(alog(p_analytic+1.)-p_analytic/(p_analytic+1.))
s = (20000. / mr_analytic )^(1./3.)
del_Psi = findgen(401)/10.
N_orb = del_Psi*!dtor / (2.*!pi) / s
Nmin = min(abs(N_orb-10.),min_sub)
print, 'Psi after 10 orbit=',del_Psi[min_sub],' degrees'
;
plot, del_Psi, N_orb, xtitle=textoidl('\Delta\Psi [deg]'), ytitle=textoidl('N_{orb}')
oplot, [del_Psi[min_sub]], [N_orb[min_sub]], psym=1, color=255

T_psi = 0.314925373 	; 1 orbit [Gyr]
t_orbit = findgen(35)/10. ;3.4
Psi = 4.*s*2.*!dpi*t_orbit/T_psi *!radeg
print, 'Psi after 3.4Gyr=',Psi
plot,t_orbit,Psi,xtitle='t [Gyr]', ytitle=textoidl('\Psi')

;=== subhalo orbit ===
spawn, 'ls snapshot_subhalo > temp2'
chr_rdtbl, 'temp2', 0, fname
fname = 'snapshot_subhalo/'+reform(fname)
spawn, 'rm -f temp2'
N_pts = N_elements(fname)-11
t_subhalo = dblarr(N_pts,N_subhalos)
x_subhalo = dblarr(N_pts,N_subhalos)
y_subhalo = dblarr(N_pts,N_subhalos)
z_subhalo = dblarr(N_pts,N_subhalos)
vx_subhalo = dblarr(N_pts,N_subhalos)
vy_subhalo = dblarr(N_pts,N_subhalos)
vz_subhalo = dblarr(N_pts,N_subhalos)
de_subhalo = dblarr(N_pts,N_subhalos)
for i=0,N_pts-1 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  arr = double(arr)
  t_subhalo[i,*] = arr[0,*]
  x_subhalo[i,*] = arr[1,*]
  y_subhalo[i,*] = arr[2,*]
  z_subhalo[i,*] = arr[3,*]
  vx_subhalo[i,*] = arr[4,*]
  vy_subhalo[i,*] = arr[5,*]
  vz_subhalo[i,*] = arr[6,*]
  de_subhalo[i,*] = arr[7,*]
endfor

r_subhalo = sqrt(x_subhalo^2.+y_subhalo^2.+z_subhalo^2.)
rxy_subhalo = sqrt(x_subhalo^2.+y_subhalo^2.)
sub_2950 = where(t_subhalo eq 3200.)
;sub_2950 = where(t_subhalo eq 3350.)
;sub_2950 = where(t_subhalo eq 2950.)
xmin = -100;min(x_subhalo)*1.1
xmax = 100;max(x_subhalo)*1.1
ymin = -100;min(y_subhalo)*1.1
ymax = 100;max(y_subhalo)*1.1
zmin = -100;min(z_subhalo)*1.1
zmax = 100;max(z_subhalo)*1.1
rmin = 0
rmax = 200;max(r_subhalo)*1.1

if N_subhalos le 100 then begin

plot, x_part1, y_part1, psym=0, xr=[xmin,xmax], yr=[ymin,ymax], xtitle='x [kpc]', ytitle='y [kpc]',/isotropic
oplot, [x_3_20Gyr], [y_3_20Gyr], psym=7, color=255
oplot, [x_init], [y_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,x_subhalo[*,i],y_subhalo[*,i],color=70
  oplot,[x_subhalo[sub_2950[i]]],[y_subhalo[sub_2950[i]]],psym=7,color=255
  oplot,[x_subhalo[0,i]],[y_subhalo[0,i]],psym=5,color=215
endfor

plot, y_part1, z_part1, psym=0, xr=[ymin,ymax], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]',/isotropic
oplot, [y_3_20Gyr], [z_3_20Gyr], psym=7, color=255
oplot, [y_init], [z_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,y_subhalo[*,i],z_subhalo[*,i],color=70
  oplot, [y_subhalo[sub_2950[i]]],[z_subhalo[sub_2950[i]]],psym=7,color=255
  oplot, [y_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
endfor

plot, x_part1, z_part1, psym=0, xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]',/isotropic
oplot, [x_3_20Gyr], [z_3_20Gyr], psym=7, color=255
oplot, [x_init], [z_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,x_subhalo[*,i],z_subhalo[*,i],color=70
  oplot, [x_subhalo[sub_2950[i]]],[z_subhalo[sub_2950[i]]],psym=7,color=255
  oplot, [x_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
endfor

plot, r_part1, z_part1, psym=0, xr=[rmin,rmax], yr=[zmin,zmax], xtitle='r [kpc]', ytitle='z [kpc]', /isotropic
oplot, [r_3_20Gyr], [z_3_20Gyr], psym=7, color=255
oplot, [r_init], [z_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,rxy_subhalo[*,i],z_subhalo[*,i],color=70
  oplot, [rxy_subhalo[sub_2950,i]],[z_subhalo[sub_2950,i]],psym=7,color=255
  oplot, [rxy_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
endfor

endif

;plot,t_subhalo/1000.,r_subhalo,xtitle='t [Gyr]',ytitle='r [kpc]'
;oplot, [t_subhalo[sub_2950]],[r_subhalo[sub_2950]],psym=7,color=255
;plothist,de_subhalo,bin=stddev(de_subhalo)/5.,xtitle='de/E0'

erase
!p.multi=0
d_sub_sat_min = dblarr(N_subhalos)
d_sub_sat_min[*] = 9e9
t_min = dblarr(N_subhalos)
v_sub_sat_min = dblarr(N_subhalos)

for i=0,N_subhalos-1 do begin
;  d_min = 9e9
  for k=0,N_elements(t_subhalo[*,i])-12 do begin
    t_same = where(t_part1 eq t_subhalo[k,i])
    d = sqrt((x_part1[t_same]-x_subhalo[k,i])^2.+(y_part1[t_same]-y_subhalo[k,i])^2.+(z_part1[t_same]-z_subhalo[k,i])^2.)
    if (d lt d_sub_sat_min[i]) then begin
       d_sub_sat_min[i] = d
       t_min[i] = t_part1[t_same]
       v_sub_sat_min[i] = sqrt((Vx_part1[t_same]-Vx_subhalo[k,i])^2.+(Vy_part1[t_same]-Vy_subhalo[k,i])^2.+(Vz_part1[t_same]-Vz_subhalo[k,i])^2.)
    endif
  endfor
;  d_subhalo_sat[i] = d_min
endfor

if N_subhalos ge 2 then begin
  plothist,d_sub_sat_min,bin=10,xtitle='minimum distance btw sub-sat [kpc]',xr=[0.1,max(d_sub_sat_min)*1.1],/xlog
;plot,t_subhalo/1000.,d_subhalo_sat,xtitle='t [Gyr]',ytitle='distance (subhalo-sat) [kpc]'
;oplot, [t_subhalo[sub_2950]],[d_subhalo_sat[sub_2950]],psym=7,color=255
endif


;=== check the energy change by subhalos ===
spawn, 'ls snapshot > temp2'
chr_rdtbl, 'temp2', 0, fname
fname = 'snapshot/'+reform(fname)
spawn, 'rm -f temp2'
N_pts = N_elements(fname)-11
t_part = dblarr(N_pts,N_stars)
x_part = dblarr(N_pts,N_stars)
y_part = dblarr(N_pts,N_stars)
z_part = dblarr(N_pts,N_stars)
vx_part = dblarr(N_pts,N_stars)
vy_part = dblarr(N_pts,N_stars)
vz_part = dblarr(N_pts,N_stars)
de_part = dblarr(N_pts,N_stars)

help,N_pts,N_subhalos,N_stars
d = dblarr(N_pts,N_subhalos,N_stars)
v = dblarr(N_pts,N_subhalos,N_stars)

for i=0,N_pts-1 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  arr = double(arr)
  t_part[i,*] = arr[0,*]
  x_part[i,*] = arr[1,*]
  y_part[i,*] = arr[2,*]
  z_part[i,*] = arr[3,*]
  vx_part[i,*] = arr[4,*]
  vy_part[i,*] = arr[5,*]
  vz_part[i,*] = arr[6,*]
  de_part[i,*] = arr[7,*]
  
  for k=0,N_subhalos-1 do begin
    d[i,k,*] = sqrt((x_part[i,*]-x_subhalo[i,k])^2.+(y_part[i,*]-y_subhalo[i,k])^2.+(z_part[i,*]-z_subhalo[i,k])^2.)
    v[i,k,*] = sqrt((vx_part[i,*]-vx_subhalo[i,k])^2.+(vy_part[i,*]-vy_subhalo[i,k])^2.+(vz_part[i,*]-vz_subhalo[i,k])^2.)
  endfor
endfor
min_d = min(d,dimension=1,sub_min)
sub_min_d = fix((double(sub_min/N_pts)-sub_min/N_pts)*N_pts)
help,min_d,sub_min_d

;G = 4.3e-6
;chr_rdtbl,'frog_subhalo.dat',0,arr
;M_subhalo = double(arr[6,*])
;
;delE_bysubhalo = dblarr(N_stars)
;for i=0,N_stars-1 do begin
;  delv_bysubhalo = 2.*G*M_subhalo / (min_d[*,i]*v[sub_min_d[*,i],*,i])
;  delE_bysubhalo[i] = total(1./2. * delv_bysubhalo^2.)
;endfor
;
;delv_bysubhalo = G*M_subhalo / (d_sub_sat_min*v_sub_sat_min)
;delE_bysubhalo = 1./2. * delv_bysubhalo^2.
;forprint, d_sub_sat_min,v_sub_sat_min,delv_bysubhalo,delE_bysubhalo,textout=2
;
;chr_rdtbl,'snapshot/snap2950',0,arr
;E_total = double(arr[8,*])
;chr_rdtbl,'snapshot_wo_subhalo/snap2950',0,arr
;E_total_wo_subhalo = double(arr[8,*])

;plot,E_total_wo_subhalo,E_total,xtitle='E w/o subhalos',ytitle='E w/ subhalos',psym=3;,/xlog,/ylog
;plot,E_total_wo_subhalo+delE_bysubhalo,E_total,xtitle='E w/o subhalos',ytitle='E w/ subhalos',psym=3;,/xlog,/ylog
;plothist,E_total_wo_subhalo,bin=1e2,xtitle='E_total w/o subhalo'


;=== radial distribution of subhalos ===
chr_rdtbl,'BI5',1,arr
arr = double(arr)
Msat_BI5 = arr[0,*] * 1.77e12
x_BI5 = arr[1,*] * 24.6/0.5
y_BI5 = arr[2,*] * 24.6/0.5
z_BI5 = arr[3,*] * 24.6/0.5
r_BI5 = sqrt(x_BI5^2.+y_BI5^2.+z_BI5^2.)
chr_rdtbl,'snapshot_subhalo/subhalo0000',0,arr
arr = double(arr)
x_subhalo0000 = arr[1,*]
y_subhalo0000 = arr[2,*]
z_subhalo0000 = arr[3,*]
r_subhalo0000 = sqrt(x_subhalo0000^2.+y_subhalo0000^2.+z_subhalo0000^2.)
chr_rdtbl,'snapshot_subhalo/subhalo1500',0,arr
arr = double(arr)
x_subhalo1500 = arr[1,*]
y_subhalo1500 = arr[2,*]
z_subhalo1500 = arr[3,*]
r_subhalo1500 = sqrt(x_subhalo1500^2.+y_subhalo1500^2.+z_subhalo1500^2.)
chr_rdtbl,'snapshot_subhalo/subhalo3200',0,arr
;chr_rdtbl,'snapshot_subhalo/subhalo3350',0,arr
;chr_rdtbl,'snapshot_subhalo/subhalo2950',0,arr
arr = double(arr)
x_subhalo2950 = arr[1,*]
y_subhalo2950 = arr[2,*]
z_subhalo2950 = arr[3,*]
r_subhalo2950 = sqrt(x_subhalo2950^2.+y_subhalo2950^2.+z_subhalo2950^2.)

;chr_rdtbl,'frog_subhalo.dat',0,arr
M_subhalo = double(arr[6,*])
;M_subhalo = replicate(2.349e-5*1.77e12,100);double(arr[6,*])
;tmp_r = sqrt(double(arr[0,*])^2.+double(arr[1,*])^2.+double(arr[2,*])^2.)
;forprint,r_subhalo0000,tmp_r,textout=2

rbin = 1.
dr = findgen(5000)*rbin
;Nrho_dr = dblarr(N_elements(dr))
rho_dr_all = dblarr(N_elements(dr))
rho_dr0000 = dblarr(N_elements(dr))
rho_dr1500 = dblarr(N_elements(dr))
rho_dr2950 = dblarr(N_elements(dr))
for i=0,N_elements(dr)-1 do begin
  in_dr_all = where(r_BI5 ge dr[i] and r_BI5 lt dr[i]+rbin)
  in_dr0000 = where(r_subhalo0000 ge dr[i] and r_subhalo0000 lt dr[i]+rbin)
  in_dr1500 = where(r_subhalo1500 ge dr[i] and r_subhalo1500 lt dr[i]+rbin)
  in_dr2950 = where(r_subhalo2950 ge dr[i] and r_subhalo2950 lt dr[i]+rbin)
  area_dr = 4./3.*!pi*((dr[i]+rbin)^3.-dr[i]^3.)
  if (in_dr_all[0] ne -1) then rho_dr_all[i] = total(Msat_BI5[in_dr_all])/area_dr
  if (in_dr0000[0] ne -1) then rho_dr0000[i] = total(M_subhalo[in_dr0000])/area_dr
  if (in_dr1500[0] ne -1) then rho_dr1500[i] = total(M_subhalo[in_dr1500])/area_dr
  if (in_dr2950[0] ne -1) then rho_dr2950[i] = total(M_subhalo[in_dr2950])/area_dr
endfor

!p.multi=0
plot,[0],[0],/nodata,ytitle=textoidl('\rho / M'+sunsymbol()+'kpc^{-3}'),xr=[0.4,8000],yr=[1e-8,1e+10],/xlog,/ylog
oplot,dr+rbin/2.,rho_dr_all,psym=1
oplot,dr+rbin/2.,rho_dr0000,psym=4,color=70
oplot,dr+rbin/2.,rho_dr1500,psym=5,color=120
oplot,dr+rbin/2.,rho_dr2950,psym=6,color=255
legend,['0.00Gyr','1.50Gyr','2.95Gyr'],psym=[4,5,6],color=[70,120,255],box=0,/right

;=== orbits of subhalos ===
!p.multi=[0,5,5]
for i=0,N_subhalos-1 do begin
  plot,r_subhalo[*,i],z_subhalo[*,i],xtitle='r [kpc]',ytitle='z [kpc]',/isotropic
  oplot,[r_subhalo[sub_2950,i]],[z_subhalo[sub_2950,i]],psym=7,color=255
  oplot,[r_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
endfor

device, /close

END
