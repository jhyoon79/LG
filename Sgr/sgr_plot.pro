pro sgr_plot

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
delJz_part1 = arr[10,*]
r_part1 = sqrt(x_part1^2.+y_part1^2.)
Vr_part1 = sqrt(Vx_part1^2.+Vy_part1^2.)
r3D_part1 = sqrt(x_part1^2.+y_part1^2.+z_part1^2.)
V_part1 = sqrt(Vx_part1^2.+Vy_part1^2.+Vz_part1^2.)
print,'Rmin=',min(r3D_part1),'    Rmax=',max(r3D_part1)

N_arr = N_elements(t_part1)
t_init = t_part1[0]
x_init = x_part1[0]  &  Vx_init = Vx_part1[0]
y_init = y_part1[0]  &  Vy_init = Vy_part1[0]
z_init = z_part1[0]  &  Vz_init = Vz_part1[0]
r_init = r_part1[0]  &  Vr_init = Vr_part1[0]
t_1Gyr = where(abs(t_part1) eq 1000.)
t_2Gyr = where(abs(t_part1) eq 2000.)
t_2_95Gyr = where(abs(t_part1) eq 3200.)
;t_2_95Gyr = where(abs(t_part1) eq 3350.)
;t_2_95Gyr = where(abs(t_part1) eq 2950.)
t_500Myr = where(abs(t_part1) eq 500.)
x_1Gyr = x_part1[t_1Gyr]  &  Vx_1Gyr = Vx_part1[t_1Gyr]
y_1Gyr = y_part1[t_1Gyr]  &  Vy_1Gyr = Vy_part1[t_1Gyr]
z_1Gyr = z_part1[t_1Gyr]  &  Vz_1Gyr = Vz_part1[t_1Gyr]
r_1Gyr = r_part1[t_1Gyr]  &  Vr_1Gyr = Vr_part1[t_1Gyr]
x_2Gyr = x_part1[t_2Gyr]  &  Vx_2Gyr = Vx_part1[t_2Gyr]
y_2Gyr = y_part1[t_2Gyr]  &  Vy_2Gyr = Vy_part1[t_2Gyr]
z_2Gyr = z_part1[t_2Gyr]  &  Vz_2Gyr = Vz_part1[t_2Gyr]
r_2Gyr = r_part1[t_2Gyr]  &  Vr_2Gyr = Vr_part1[t_2Gyr]
x_500Myr = x_part1[t_500Myr]  &  Vx_500Myr = Vx_part1[t_500Myr]
y_500Myr = y_part1[t_500Myr]  &  Vy_500Myr = Vy_part1[t_500Myr]
z_500Myr = z_part1[t_500Myr]  &  Vz_500Myr = Vz_part1[t_500Myr]
r_500Myr = r_part1[t_500Myr]  &  Vr_500Myr = Vr_part1[t_500Myr]
x_2_95Gyr = x_part1[t_2_95Gyr]  &  Vx_2_95Gyr = Vx_part1[t_2_95Gyr]
y_2_95Gyr = y_part1[t_2_95Gyr]  &  Vy_2_95Gyr = Vy_part1[t_2_95Gyr]
z_2_95Gyr = z_part1[t_2_95Gyr]  &  Vz_2_95Gyr = Vz_part1[t_2_95Gyr]
r_2_95Gyr = r_part1[t_2_95Gyr]  &  Vr_2_95Gyr = Vr_part1[t_2_95Gyr]
t_final = t_part1[N_arr-1]
x_final = x_part1[N_arr-1]  &  Vx_final = Vx_part1[N_arr-1]
y_final = y_part1[N_arr-1]  &  Vy_final = Vy_part1[N_arr-1]
z_final = z_part1[N_arr-1]  &  Vz_final = Vz_part1[N_arr-1]
r_final = r_part1[N_arr-1]  &  Vr_final = Vr_part1[N_arr-1]
t_mark_1Gyr = 't=1Gyr'
t_mark_2Gyr = 't=2Gyr'
t_mark_2_95Gyr = 't=2_95Gyr'
t_mark_500Myr = 't=500Myr'
t_mark_final = textoidl('t_{final}=')+strtrim(string(t_final/1000.,f='(f6.2)'),2)+'Gyr'
openr, 99, 'frogin'
readf, 99, N_stars
readf, 99, tmp,tmp
readf, 99, tmp
readf, 99, N_subhalos
close, 99

;=== read out info of all test particles ===  
chr_rdtbl, 'frog_sgr.dat', 0, arr
x = strtrim(string(double(arr[0,0]),f='(f5.1)'),2)
y = strtrim(string(double(arr[1,0]),f='(f5.1)'),2)
z = strtrim(string(double(arr[2,0]),f='(f5.1)'),2)
Vx = strtrim(fix(arr[3,0]),2)
Vy = strtrim(fix(arr[4,0]),2)
Vz = strtrim(fix(arr[5,0]),2)
Info_posi = 'x='+x+' y='+y+' z='+z
Info_vel = textoidl('V_x='+Vx+' V_y='+Vy+' V_z='+Vz)
Info1 = 'Part0001'+textoidl(' N_{particle}=900')

;chr_rdtbl, 'frog_sgr.dat', 0, arr
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

chr_rdtbl, 'sgr_sigma.dat', 0, arr
sigma_v = textoidl('\sigma_v=')+string(arr[0,0],f='(f4.2)')+'km/s'
sigma_r = textoidl('  \sigma_r=')+string(arr[1,0],f='(f5.3)')+'kpc'
sigma = sigma_v + sigma_r


;=== find apogee & perigee ===
t_part1 = t_part1/1000.
N_orbit = N_elements(t_part1)
dr1 = shift(r3D_part1,1)-r3D_part1
dr2 = r3D_part1 - shift(r3D_part1,-1)
dr1[0] = 0
dr2[N_orbit-1] = 0
sub = where(dr1*dr2 lt 0)
m_r3D_part1 = mean(r3D_part1[sub])
apogee = where(r3D_part1[sub] gt m_r3D_part1,complement=perigee)
r_apo = r3D_part1[sub[apogee]]
r_peri = r3D_part1[sub[perigee]]
t_apo = t_part1[sub[apogee]]
t_peri = t_part1[sub[perigee]]

tR_all = t_apo-shift(t_apo,1)
tR_all = tR_all[where(tR_all gt 0)]
tR = mean(tR_all)
openw,3,'t_apo_Sgr'
printf,3,'t_apo(Gyr)'
for i=0,N_elements(t_apo)-1 do printf,3,t_apo[i]
printf,3,'mean tR=',tR
close,3


@plot_setting
!p.charsize=1.2
device, file='sgr_plot.ps', /color,/landscape

;=== initial position of particles: the shape of sgr at the initial position ===
!p.multi = [0,2,2]
plot, x_GC, y_GC, psym=3,xtitle='x [kpc]', ytitle='y [kpc]', /isotropic, title='initial position of particles'
plot, y_GC, z_GC, psym=3,xtitle='x [kpc]', ytitle='y [kpc]', /isotropic, title='initial position of particles'
plot, z_GC, x_GC, psym=3,xtitle='x [kpc]', ytitle='y [kpc]', /isotropic, title='initial position of particles'
multiplot,/reset
erase

;=== the orbit of sgr on the meridional plane ===
!p.multi=[0,1,2]
multiplot
plot, t_part1, V_part1,yr=[min(V_part1)*0.9,max(V_part1)*1.1], ytitle=textoidl('V_{3D} [km/s]')
oplot, [t_part1[t_2_95Gyr]], [V_part1[t_2_95Gyr]], psym=7, color=255

multiplot
plot, t_part1, r3D_part1,yr=[min(r3D_part1)*0.9,max(r3D_part1)*1.1],xtitle='t [Gyr]', ytitle='r [kpc]'
oplot, [t_part1[t_2_95Gyr]], [r3D_part1[t_2_95Gyr]], psym=7, color=255
multiplot, /reset
erase

;=== the orbit of the central particle(part0001) for 2.95Gyr ===
!p.multi=[0,2,2]
xmin = -100  &  xmax = 100
ymin = -100  &  ymax = 100
zmin = -100  &  zmax = 100
rmin = 0    &  rmax = 100
plot, x_part1, y_part1, psym=0, xr=[xmin,xmax], yr=[ymin,ymax], xtitle='x [kpc]', ytitle='y [kpc]', title=Info1, /isotropic
oplot, [x_2_95Gyr], [y_2_95Gyr], psym=7, color=255
xyouts, x_init, y_init, 't=0', charsize=1, color=255
xyouts, x_1Gyr, y_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, x_2Gyr, y_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, x_500Myr, y_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, x_final, y_final, t_mark_final, charsize=1, color=255

plot, y_part1, z_part1, psym=0, xr=[ymin,ymax], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]', title=Info_posi, /isotropic
oplot, [y_2_95Gyr], [z_2_95Gyr], psym=7, color=255
xyouts, y_init, z_init, 't=0', charsize=1, color=255
xyouts, y_1Gyr, z_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, y_2Gyr, z_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, y_500Myr, z_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, y_final, z_final, t_mark_final, charsize=1, color=255

plot, x_part1, z_part1, psym=0, xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]', title=Info_vel, /isotropic
oplot, [x_2_95Gyr], [z_2_95Gyr], psym=7, color=255
xyouts, x_init, z_init, 't=0', charsize=1, color=255
xyouts, x_1Gyr, z_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, x_2Gyr, z_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, x_500Myr, z_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, x_final, z_final, t_mark_final, charsize=1, color=255

plot, r_part1, z_part1, psym=0, xr=[rmin,rmax], yr=[zmin,zmax], xtitle='r [kpc]', ytitle='z [kpc]', /isotropic
oplot, [r_2_95Gyr], [z_2_95Gyr], psym=7, color=255
xyouts, r_init, z_init, 't=0', charsize=1, color=255
xyouts, r_1Gyr, z_1Gyr, t_mark_1Gyr, charsize=1, color=255
xyouts, r_2Gyr, z_2Gyr, t_mark_2Gyr, charsize=1, color=255
xyouts, r_500Myr, z_500Myr, t_mark_500Myr, charsize=1, color=255
xyouts, r_final, z_final, t_mark_final, charsize=1, color=255

xr = abs(min(x_part1)) > abs(max(x_part1))
yr = abs(min(y_part1)) > abs(max(y_part1))
zr = abs(min(z_part1)) > abs(max(z_part1))
xmin = -xr  &  xmax = xr
ymin = -yr  &  ymax = yr
zmin = -zr  &  zmax = zr
rmin = 0    &  rmax = max(r_part1)
Vxmin = -250 &  Vxmax = 250
Vymin = -250 &  Vymax = 250
Vzmin = -250 &  Vzmax = 250
Vrmin = 0  &  Vrmax = 250



; read out the final position of all particles
chr_rdtbl, 'snapshot/snap8000',0, arr
;chr_rdtbl, 'snapshot/snap3200',0, arr
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
Jz = arr[10,*]

;xmin =   6  &  xmax = 10
;ymin =  -3  &  ymax =  3
;zmin =  14  &  zmax = 17.5
;Vxmin =  -60 &  Vxmax = -10
;Vymin = -100 &  Vymax = -80
;Vzmin =  -60 &  Vzmax =  20
xr = abs(min(x)) > abs(max(x))
yr = abs(min(y)) > abs(max(y))
zr = abs(min(z)) > abs(max(z))
xmin = -xr  &  xmax = xr
ymin = -yr  &  ymax = yr
zmin = -zr  &  zmax = zr


;=== final position of particles: the shape of tidal tail ===
!p.multi=[0,2,2]
plot, x, y, psym=3, xr=[xmin,xmax], yr=[ymin,ymax], xtitle='x [kpc]', ytitle='y [kpc]', title='Pal 5', /isotropic
oplot, [x[0]], [y[0]], psym=1, color=255
plot, y, z, psym=3, xr=[ymin,ymax], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]', /isotropic
oplot, [y[0]], [z[0]], psym=1, color=255
plot, x, z, psym=3, xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]', /isotropic
oplot, [x[0]], [z[0]], psym=1, color=255


;=== estimate Galactic coord. ===
l = atan(y/(x+8.))*!radeg 
l_sub = where(l lt 0)
if l_sub[0] ne -1 then l[l_sub] = l[l_sub] + 360.
l_sub2 = where(l gt 360)
if l_sub2[0] ne -1 then l[l_sub2] = l[l_sub2] - 360.
b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
lcosb = l*cos(b*!dtor)


;=== estimate equatorial coord. ===
delta = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
alpha = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(delta*!dtor) )*!radeg + 282.25

alpha_sub = where(alpha lt 0)
if alpha_sub[0] ne -1 then alpha[alpha_sub] = alpha[alpha_sub] + 360.
alpha_sub2 = where(alpha gt 360)
if alpha_sub2[0] ne -1 then alpha[alpha_sub2] = alpha[alpha_sub2] - 360.


sub_t = where(t_part1 gt 3.050 and t_part1 lt 3.350)
;sub_t = where(t_part1 gt 3200 and t_part1 lt 3500)
;sub_t = where(t_part1 gt 2850 and t_part1 lt 3050)
l_part1 = atan(y_part1/(x_part1+8.))*!radeg 
l_part1_sub = where(l_part1 lt 0)
if l_part1_sub[0] ne -1 then l_part1[l_part1_sub] = l_part1[l_part1_sub] + 360.
l_part1_sub2 = where(l_part1 gt 360)
if l_part1_sub2[0] ne -1 then l_part1[l_part1_sub2] = l_part1[l_part1_sub2] - 360.
b_part1 = atan(z_part1/sqrt((x_part1+8.)^2.+y_part1^2.))*!radeg
lcosb_part1 = l_part1*cos(b_part1*!dtor)


delta_part1 = asin( cos(b_part1*!dtor)*sin((l_part1-33.)*!dtor)*sin(62.6*!dtor)+sin(b_part1*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_part1 = asin( (cos(b_part1*!dtor)*sin((l_part1-33.)*!dtor)*cos(62.6*!dtor)-sin(b_part1*!dtor)*sin(62.6*!dtor))/cos(delta_part1*!dtor) )*!radeg + 282.25
alpha_part1_sub = where(alpha_part1 lt 0)
if alpha_part1_sub[0] ne -1 then alpha_part1[alpha_part1_sub] = alpha_part1[alpha_part1_sub] + 360.
alpha_part1_sub2 = where(alpha_part1 gt 360)
if alpha_part1_sub2[0] ne -1 then alpha_part1[alpha_part1_sub2] = alpha_part1[alpha_part1_sub2] - 360.
;sub_b = where(lcosb_part1 gt -12 and lcosb_part1 lt 12 and t_part1 gt 2800 and t_part1 lt 3100)
sub_b = sub_t
alpha_sgr = alpha[0]
delta_sgr = delta[0]
l_sgr = l[0]
lcosb_sgr = lcosb[0]
b_sgr = b[0]


;=== Sgr in equatorial coord. ===
!p.multi=0
plot,[0],[0],/nodata,xr=[360,0],yr=[-90,90],xtitle=textoidl('\alpha [deg]'),ytitle=textoidl('\delta [deg]'),title=sigma
oplot,alpha_part1,delta_part1,color=100;,psym=1
;oplot,alpha_part1[sub_t],delta_part1[sub_t],color=100;,psym=1
oplot,alpha,delta,psym=8,symsize=0.2
oplot,[alpha_sgr],[delta_sgr], psym=1, color=255,symsize=1.8,thick=8
print,'RA=',alpha_sgr,'   Dec=',delta_sgr

;=== Sgr in l & b plane ===
plot, l,b,xr=[360,0],yr=[-90,90],psym=3,xtitle='l [deg]',ytitle='b [deg]',title=sigma
oplot,l_part1[sub_t],b_part1[sub_t],color=70;,psym=1
oplot,[l_sgr],[b_sgr], psym=1, color=255,symsize=1.8,thick=8


;=== sgr in lcosb & b plane ===
;if (N_subhalos eq 0) then begin 
;  plot, lcosb, b, psym=3, xr=[8,-4], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title=sigma
;endif else begin
  plot, lcosb,b,xr=[max(lcosb),min(lcosb)],yr=[min(b),max(b)],psym=3,xtitle='l cos b [deg]',ytitle='b [deg]',title=sigma
;endelse
oplot, [lcosb_sgr], [b_sgr], psym=1, color=255,symsize=1.8,thick=8

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
legend, ['sgr center','  orbit   '],color=[255,70],psym=[1,0],box=0,/right

;=== contour plot of sgr ===
contour_bin = 30
levels = indgen(300)*contour_bin + contour_bin
pts_contour,lcosb,b,20,15,gv,x,y
;if (N_subhalos eq 0) then begin 
;  contour,gv,x,y,xr=[8,-4],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1
;endif else begin
  contour,gv,x,y,xr=[max(lcosb),min(lcosb)],yr=[min(b),max(b)],xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1
;endelse


;=== sgr in lcosb & b plane 2 ===
plot, lcosb, b, psym=3, xr=[8,-4], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title=sigma
oplot, [lcosb_sgr], [b_sgr], psym=1, color=255,symsize=1.8,thick=8

oplot,lcosb_part1[sub_t],b_part1[sub_t],color=70;,psym=1
legend, ['rms from orbit='+strtrim(string(orbit_rms,f='(f12.8)'),2)],box=0
legend, ['sgr center','  orbit   '],color=[255,70],psym=[1,0],box=0,/right
;=== contour plot of sgr ===
contour,gv,x,y,xr=[8,-4],yr=[42,49],xtitle='l cos b [deg]',ytitle='b [deg]',levels=levels,ticklen=1, xgridstyle=1, ygridstyle=1


;=== density profile of Pal 5 ===
radius = 60.*sqrt((lcosb-lcosb_sgr)^2.+(b-b_sgr)^2.)
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

;bsize = stddev(delJz_part1)/5.
;plothist, delJz_part1, bin=bsize, xtitle=textoidl('\DeltaJ_z / J_z (part1)'), ytitle='#'

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
Msat = 3.e8
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
  
  bsize = stddev(epsilon)/5.
  plothist, epsilon, bin=bsize,xtitle=textoidl('\epsilon')
  bsize = stddev(q)/5.
;  plot,[0],[0],/nodata,xsty=4,ysty=4

if N_subhalos eq 0 then begin
  plothist, q, bin=bsize, xtitle='q'
  sig = stddev(q)
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f12.8)'),2),alignment=0.5
;  plot,[0],[0],/nodata,xsty=4,ysty=4
  
  cut_sub = 3.5
  lcosb_sub = where(lcosb gt cut_sub or lcosb lt -1.*cut_sub)
  bsize = stddev(delE)/5.
  cut_title = string(cut_sub,f='(f3.1)')
  plothist, delE[lcosb_sub], bin=bsize,xtitle=textoidl('\DeltaE (E_{debris}-E_{sat})'),title='lcosb<-'+cut_title+' or lcosb>'+cut_title
  bsize = stddev(q)/5.
  plothist, q[lcosb_sub], bin=bsize, xtitle='q',title='lcosb<-'+cut_title+' or lcosb>'+cut_title
  plot,[0],[0],/nodata,xsty=4,ysty=4
  
  
  ;=== test particles with different energy ===
  plot, lcosb, b, psym=3, xr=[12,-10], yr=[42,49], xtitle='l cos b [deg]', ytitle='b [deg]', title=sigma
  oplot, [lcosb[0]], [b[0]], psym=1, color=255, thick=8
  fit = poly_fit(lcosb_part1[sub_b],b_part1[sub_b],2,yerror=yerror,yfit=yfit)
  y = fit[0] + fit[1]*lcosb + fit[2]*lcosb^2.
  oplot, lcosb_part1[sub_t], b_part1[sub_t], color=70
  rms = sqrt(total((b-y)^2.)/N_elements(b))
  legend, ['rms='+strtrim(string(rms,f='(f12.8)'),2)], box=0
  
  
  ;=== Angular momentum plot ===
  bsize = stddev(scale_dJ)/5.
  plothist, scale_dJ, bin=bsize, xtitle=textoidl('\DeltaJ / sJ (s=r_{tide}/R_{peri})')
  sig = stddev(scale_dJ)
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f14.8)'),2),alignment=0.5
  
  
  ;=== z-Angular momentum plot ===
  dJz = Jz-Jz[0]
  bsize = stddev(dJz)/5.
  plothist, dJz, bin=bsize, xtitle=textoidl('\DeltaJ_z (J_{z,debris}-J_{z,sat})')
  sig = stddev(dJz)
  in_dJz = where(dJz gt -1.*sig and dJz lt sig,complement=out_dJz)
  plothist, dJz[in_dJz], bin=bsize, /fill, /overplot, fcolor=100
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f8.3)'),2),alignment=0.5
  
  ;r_tide = 0.061
  r_tide = 0.108653
  scale_dJz = dJz / (r_tide/r_peri * Jz)
  bsize = stddev(scale_dJz)/5.
  plothist, scale_dJz, bin=bsize, xtitle=textoidl('\DeltaJ_z / sJ_z (s=r_{tide}/R_{peri})')
  sig = stddev(scale_dJz)
  oplot, [0,sig],[0,0],thick=12
  xyouts, 0,15, textoidl('\sigma=')+strtrim(string(sig,f='(f8.3)'),2),alignment=0.5
  
  
  ;=== test particles with different angular momentum ===
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
t_subhalo = t_subhalo/1000.

if N_subhalos le 100 then begin

plot, x_part1, y_part1, psym=0, xr=[xmin,xmax], yr=[ymin,ymax], xtitle='x [kpc]', ytitle='y [kpc]',/isotropic
oplot, [x_2_95Gyr], [y_2_95Gyr], psym=7, color=255
oplot, [x_init], [y_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,x_subhalo[*,i],y_subhalo[*,i],color=70
  oplot,[x_subhalo[sub_2950[i]]],[y_subhalo[sub_2950[i]]],psym=7,color=255
  oplot,[x_subhalo[0,i]],[y_subhalo[0,i]],psym=5,color=215
endfor

plot, y_part1, z_part1, psym=0, xr=[ymin,ymax], yr=[zmin,zmax], xtitle='y [kpc]', ytitle='z [kpc]',/isotropic
oplot, [y_2_95Gyr], [z_2_95Gyr], psym=7, color=255
oplot, [y_init], [z_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,y_subhalo[*,i],z_subhalo[*,i],color=70
  oplot, [y_subhalo[sub_2950[i]]],[z_subhalo[sub_2950[i]]],psym=7,color=255
  oplot, [y_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
endfor

plot, x_part1, z_part1, psym=0, xr=[xmin,xmax], yr=[zmin,zmax], xtitle='x [kpc]', ytitle='z [kpc]',/isotropic
oplot, [x_2_95Gyr], [z_2_95Gyr], psym=7, color=255
oplot, [x_init], [z_init], psym=5, color=215
for i=0,N_subhalos-1 do begin
  oplot,x_subhalo[*,i],z_subhalo[*,i],color=70
  oplot, [x_subhalo[sub_2950[i]]],[z_subhalo[sub_2950[i]]],psym=7,color=255
  oplot, [x_subhalo[0,i]],[z_subhalo[0,i]],psym=5,color=215
endfor

plot, r_part1, z_part1, psym=0, xr=[rmin,rmax], yr=[zmin,zmax], xtitle='r [kpc]', ytitle='z [kpc]', /isotropic
oplot, [r_2_95Gyr], [z_2_95Gyr], psym=7, color=255
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
