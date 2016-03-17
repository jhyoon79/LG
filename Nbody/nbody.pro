pro nbody


;=== scale ===
kmpkpc = 3.240756d-17
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
  ;Msun = 1.989d33
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
mu = 1.e4
ru = 0.0075
tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec

dir_data = '/media/SEADISK/LG/Nbody/pal5_final/'
dir_out = '/media/SEADISK/LG/Nbody/'

chr_rdtbl,dir_data+'SCFCEN',1,arr
arr = double(arr)
orbit_t = arr[0,*] * tu
orbit_dt = arr[1,*] * tu
orbit_x = arr[2,*] * ru
orbit_y = arr[3,*] * ru
orbit_z = arr[4,*] * ru
orbit_vx = arr[5,*] * vu
orbit_vy = arr[6,*] * vu
orbit_vz = arr[7,*] * vu
orbit_r = sqrt(orbit_x^2.+orbit_y^2.+orbit_z^2.)

chr_rdtbl,dir_data+'SCFLOG',1,arr
arr = double(arr)
log_t = arr[0,*] * tu
log_Mtot = arr[1,*] * mu
log_Mbnd = arr[2,*] * mu

spawn,'ls '+dir_data+'SNAP* > tmp5'
chr_rdtbl,'tmp5',0,arr
snap_fname = reform(arr[0,*])
spawn,'rm -f tmp5'


;=== plot t & r and mass loss ===
@plot_setting
device,file=dir_out+'nbody.ps',/color
!p.multi=[0,1,2]
multiplot
plot,orbit_t,orbit_r,ytitle='r [kpc]'
multiplot
plot,log_t,log_Mbnd,xtitle='t [Gyr]',ytitle=textoidl('M_{bnd} [M'+odot+']')
multiplot,/reset

;=== find apogee & perigee ===
N_orbit = N_elements(orbit_t)
dr1 = shift(orbit_r,1)-orbit_r
dr2 = orbit_r - shift(orbit_r,-1)
dr1[0] = 0
dr2[N_orbit-1] = 0
sub = where(dr1*dr2 lt 0)
m_orbit_r = mean(orbit_r[sub])
apogee = where(orbit_r[sub] gt m_orbit_r,complement=perigee)
Rapo = orbit_r[sub[apogee]]
Rperi = orbit_r[sub[perigee]]
t_apo = orbit_t[sub[apogee]]
t_peri = orbit_t[sub[perigee]]
if mu eq 5e4 then begin
  tR = mean([(t_peri-shift(t_peri,1))[1:8],(t_apo-shift(t_apo,1))[1:7]])
endif else begin
  tR = mean([(t_peri-shift(t_peri,1))[1:8],(t_apo-shift(t_apo,1))[1:8]])
endelse
print,'t_apo',t_apo
help,t_apo
print,Rperi,Rapo


Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
Msat = 10000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
Rperi = mean(Rperi)
p = Rperi/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s = (mu/mr)^(1./3.)
r_tide = s*Rperi	;0.108653
epsilon_nbody = s*(4.3e-6*mr/Rperi)

!p.multi=[0,2,2]
plot,orbit_x,orbit_y,xr=[-20,20],yr=[-20,20],/isotropic
plot,orbit_y,orbit_z,xr=[-20,20],yr=[-20,20],/isotropic
plot,orbit_x,orbit_z,xr=[-20,20],yr=[-20,20],/isotropic
plot,sqrt(orbit_x^2.+orbit_y^2.),orbit_z,xr=[0,25],yr=[-20,20],/isotropic
erase

!p.multi=0
ocolor = color_code(orbit_t)
plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='y',/isotropic
for i=0,N_orbit-1 do oplot,[orbit_x[i]],[orbit_y[i]],psym=8,symsize=0.3,color=ocolor[i]
nn = findgen(5) + 1
t_snap = [0,N_orbit*nn/max(nn)-1]
for i=0,N_elements(t_snap)-1 do begin
  tt = string(orbit_t[t_snap[i]],f='(f6.4)')+'Gyr'
  xyouts,orbit_x[t_snap[i]],orbit_y[t_snap[i]],tt,charsize=0.6,color=ocolor[t_snap[i]]
endfor

device,/close

N_file = N_elements(snap_fname);27
snap_time = dblarr(N_file)
Npart = 10000
b_nbody = dblarr(N_file,Npart)
lcosb_nbody = dblarr(N_file,Npart)
RA_nbody = dblarr(N_file,Npart)
Dec_nbody = dblarr(N_file,Npart)
; coord. projected in the orbital plane
Etot_nbody_cons = dblarr(N_file,Npart)
dE_cons = dblarr(N_file,Npart)
Npts = intarr(N_file,Npart)

!p.multi=[0,1,3]
device,file=dir_out+'nbody_movie.ps',/color,ysize=20,yoffset=0.5
for i=0,N_file-1 do begin
  openr,1,snap_fname[i]
  readf,1,Npart,time
  close,1
  time0 = round(time*tu*1000.)/1000.*1000

  chr_rdtbl,snap_fname[i],1,arr,/silent
  arr = double(arr)
  m_nbody = arr[0,*] * mu
  x_nbody = reform(arr[1,*] * ru)
  y_nbody = reform(arr[2,*] * ru)
  z_nbody = reform(arr[3,*] * ru)
  vx_nbody = reform(arr[4,*] * vu)
  vy_nbody = reform(arr[5,*] * vu)
  vz_nbody = reform(arr[6,*] * vu)
  pot_int_nbody = reform(arr[7,*] * vu^2.)
  pot_ext_nbody = reform(arr[8,*] * vu^2.)
  if (i eq 0) then begin
    tunbound_nbody = reform(arr[9,*])
    bounded = where(tunbound_nbody eq 0)
    escaped = where(tunbound_nbody ne 0)
    count = 0
  endif else begin
    pre_tunbound_nbody = tunbound_nbody 
    tunbound_nbody = reform(arr[9,*])
    bounded = where(tunbound_nbody eq 0)
    escaped = where(tunbound_nbody ne 0)
    recently_escaped = where(tunbound_nbody ne 0 and pre_tunbound_nbody eq 0,count)
  endelse

  Etot_nbody = (vx_nbody^2.+vy_nbody^2.+vz_nbody^2.)/2.+pot_ext_nbody+pot_int_nbody
  J_nbody = sqrt( (y_nbody*vz_nbody-z_nbody*vy_nbody)^2. + $
                  (z_nbody*vx_nbody-x_nbody*vz_nbody)^2. + $
                  (x_nbody*vy_nbody-y_nbody*vx_nbody)^2.)
  Emin = min(Etot_Nbody)
  Emax = max(Etot_Nbody)
  Ecolor = round((Etot_Nbody-Emin)/(Emax-Emin)*245.)+10

  l_nbody = atan(y_nbody/(x_nbody+8.))*!radeg 
  snap_time[i] = time0
  b_nbody[i,*] = atan(z_nbody/sqrt((x_nbody+8.)^2.+y_nbody^2.))*!radeg
  lcosb_nbody[i,*] = l_nbody*cos(b_nbody[i,*]*!dtor)
  Dec_nbody[i,*] = asin( cos(b_nbody[i,*]*!dtor)*sin((l_nbody-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody[i,*]*!dtor)*cos(62.6*!dtor) )*!radeg
  RA_nbody[i,*] = asin( (cos(b_nbody[i,*]*!dtor)*sin((l_nbody-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody[i,*]*!dtor)*sin(62.6*!dtor))/cos(Dec_nbody[i,*]*!dtor) )*!radeg + 282.25

  xmin = -20;min([min(x_nosub),min(x_nbody)])
  xmax = 20;max([max(x_nosub),max(x_nbody)])
  ymin = -20;min([min(y_nosub),min(y_nbody)])
  ymax = 20;max([max(y_nosub),max(y_nbody)])
  zmin = -20;min([min(z_nosub),min(z_nbody)])
  zmax = 20;max([max(z_nosub),max(z_nbody)])
  if (i eq 0) then begin
    m_Etot_nbody = mean(Etot_nbody)
    m_J_nbody = mean(J_nbody)
  endif
  dE_nbody = Etot_nbody-m_Etot_nbody
  q_nbody = dE_nbody/epsilon_nbody
  sJ_nbody = s*J_nbody
  dJ_nbody = J_nbody-m_J_nbody
  dJ_sJ_nbody = dJ_nbody/sJ_nbody

  plot,x_nbody,y_nbody,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',title='Nbody',/isotropic
  legend,[string(time0/1000,f='(f7.4)')+'Gyr'],box=0,charsize=1

  plot,y_nbody,z_nbody,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z [kpc]',/isotropic

  qmin = -5  &  qmax = 5
  Jmin = -4.5  &  Jmax = 4.5
;  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),title='Msat='+strtrim(mu,2)+'M'+odot
  plot,q_nbody,dJ_sJ_nbody,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),title='Msat='+strtrim(mu,2)+'M'+odot
  oplot,q_nbody[escaped],dJ_sJ_nbody[escaped],psym=3,color=255
  if (i ge 1 and count gt 0) then oplot,q_nbody[recently_escaped],dJ_sJ_nbody[recently_escaped],psym=3,color=70
  hline,0,linestyle=2
  vline,0,linestyle=2

endfor
out:device,/close
erase


device,file=dir_out+'nbody_snap.ps',/color,ysize=16
!p.multi=[0,4,7]
!p.charsize=1
for i=0,N_file-1 do begin
  x_m = median(ox_nbody[i,*])
  y_m = median(oy_nbody[i,*])
  sub_orbit = where(abs(orbit_t-snap_time[i]/1000) lt 0.05)
  new_xo = orbit_x3[sub_orbit]-x_m
  new_yo = orbit_y3[sub_orbit]-y_m
  fit = poly_fit(new_xo,new_yo,1,yfit=yfit)
  theta = -atan(fit[1])
  if mu eq 1e4 then begin
    tmp = [2,3,6,7,10,11,15,16,19,20,23]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
    if i eq 9 then theta = theta - 5*!dtor 
    if i eq 12 then theta = theta + 10*!dtor 
    if i eq 14 then theta = theta - 10*!dtor 
    if i eq 16 then theta = theta - 10*!dtor 
    if i eq 17 then theta = theta - 5*!dtor 
    if i eq 18 then theta = theta - 10*!dtor 
    if i eq 20 then theta = theta - 10*!dtor 
    if i eq 22 then theta = theta - 10*!dtor 
    if i eq 23 then theta = theta - 10*!dtor 
    if i eq 24 then theta = theta - 180*!dtor 
  endif else if mu eq 5e4 then begin
    tmp = [2,3,6,7,11,12,16,17,20,21,23]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
    if i eq 2 then theta = theta + 5*!dtor 
    if i eq 7 then theta = theta - 5*!dtor 
    if i eq 10 then theta = theta - 10*!dtor 
    if i eq 13 then theta = theta + 10*!dtor 
    if i eq 15 then theta = theta - 10*!dtor 
    if i eq 17 then theta = theta - 5*!dtor 
    if i eq 18 then theta = theta - 5*!dtor 
    if i eq 21 then theta = theta - 5*!dtor 
    if i eq 22 then theta = theta + 5*!dtor 
    if i eq 23 then theta = theta - 5*!dtor 
  endif else if mu eq 1e5 then begin
    tmp = [2,3,6,7,10,11,15,16,18,19,20,23]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
    if i eq 9 then theta = theta - 5*!dtor 
    if i eq 12 then theta = theta + 10*!dtor 
    if i eq 14 then theta = theta - 20*!dtor 
    if i eq 16 then theta = theta - 20*!dtor 
    if i eq 18 then theta = theta + 20*!dtor 
    if i eq 20 then theta = theta - 90*!dtor 
    if i eq 22 then theta = theta - 10*!dtor 
    if i eq 23 then theta = theta - 10*!dtor 
  endif
  rot_xo = cos(theta)*new_xo-sin(theta)*new_yo
  rot_yo = sin(theta)*new_xo+cos(theta)*new_yo
  rot_ox_nbody = cos(theta)*(ox_nbody[i,*]-x_m)-sin(theta)*(oy_nbody[i,*]-y_m)
  rot_oy_nbody = sin(theta)*(ox_nbody[i,*]-x_m)+cos(theta)*(oy_nbody[i,*]-y_m)
  rot_ox_nbody_ex_1peri = cos(theta)*(ox_nbody_ex_1peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_1peri[i,*]-y_m)
  rot_oy_nbody_ex_1peri = sin(theta)*(ox_nbody_ex_1peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_1peri[i,*]-y_m)
  rot_ox_nbody_ex_2peri = cos(theta)*(ox_nbody_ex_2peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_2peri[i,*]-y_m)
  rot_oy_nbody_ex_2peri = sin(theta)*(ox_nbody_ex_2peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_2peri[i,*]-y_m)
  rot_ox_nbody_ex_3peri = cos(theta)*(ox_nbody_ex_3peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_3peri[i,*]-y_m)
  rot_oy_nbody_ex_3peri = sin(theta)*(ox_nbody_ex_3peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_3peri[i,*]-y_m)
  rot_ox_nbody_ex_4peri = cos(theta)*(ox_nbody_ex_4peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_4peri[i,*]-y_m)
  rot_oy_nbody_ex_4peri = sin(theta)*(ox_nbody_ex_4peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_4peri[i,*]-y_m)
  multiplot
  tname = (i eq 0) ? textoidl('M_{sat}=')+strtrim(mu,2)+'M'+odot:''
  plot,[0],[0],/nodata,xr=[-3.2,3.2],yr=[-1.6,0.9],title=tname,/isotropic
  oplot,rot_xo,rot_yo,color=100,thick=1
  oplot,rot_ox_nbody,rot_oy_nbody,psym=3
  oplot,rot_ox_nbody_ex_4peri,rot_oy_nbody_ex_4peri,psym=3,color=220
  oplot,rot_ox_nbody_ex_3peri,rot_oy_nbody_ex_3peri,psym=3,color=150
  oplot,rot_ox_nbody_ex_2peri,rot_oy_nbody_ex_2peri,psym=3,color=70
  oplot,rot_ox_nbody_ex_1peri,rot_oy_nbody_ex_1peri,psym=3,color=30
;  oplot,new_xo,new_yo,color=220
;  oplot,new_xo,yfit,color=255
  legend,[string(snap_time[i]/1000,f='(f7.4)')+'Gyr'],box=0,charsize=0.6,position=[0,-0.5]
endfor
multiplot,/reset
device,/close
erase

device,file=dir_out+'nbody_current_snap.ps',/color,ysize=25
!p.multi=[0,1,6]
!p.charsize=1
;multiplot
plot,orbit_t,orbit_r,ytitle='r [kpc]'
vline,snap_time[N_file-10]/1000,linestyle=1
vline,snap_time[N_file-7]/1000,linestyle=1
vline,snap_time[N_file-5]/1000,linestyle=1
vline,snap_time[N_file-2]/1000,linestyle=1
;vline,snap_time[N_file-18]/1000,linestyle=1
;vline,snap_time[N_file-15]/1000,linestyle=1
;vline,snap_time[N_file-12]/1000,linestyle=1
;vline,snap_time[N_file-10]/1000,linestyle=1
;vline,snap_time[N_file-26]/1000,linestyle=1
;vline,snap_time[N_file-23]/1000,linestyle=1
;vline,snap_time[N_file-20]/1000,linestyle=1
;vline,snap_time[N_file-18]/1000,linestyle=1
;multiplot
plot,log_t,log_Mbnd,xtitle='t [Gyr]',ytitle=textoidl('M_{bnd} [M'+odot+']')
legend,[textoidl('t_R=')+string(tR,f='(f6.4)')+'Gyr',textoidl('t_{peri,recent}=')+string(max(t_peri),f='(f6.4)')+'Gyr'],box=0,/right

for j=0,3 do begin
  if j eq 0 then i = N_file-10
  if j eq 1 then i = N_file-7
  if j eq 2 then i = N_file-5
  if j eq 3 then i = N_file-2
;  if j eq 0 then i = N_file-79
;  if j eq 1 then i = N_file-77
;  if j eq 2 then i = N_file-74
;  if j eq 3 then i = N_file-2
;  if j eq 0 then i = N_file-26
;  if j eq 1 then i = N_file-23
;  if j eq 2 then i = N_file-20
;  if j eq 3 then i = N_file-18
  x_m = median(ox_nbody[i,*])
  y_m = median(oy_nbody[i,*])
  sub_orbit = where(abs(orbit_t-snap_time[i]/1000) lt 0.05)
  new_xo = orbit_x3[sub_orbit]-x_m
  new_yo = orbit_y3[sub_orbit]-y_m
  fit = poly_fit(new_xo,new_yo,1,yfit=yfit)
  theta = -atan(fit[1])
  if mu eq 1e4 then begin
    tmp = [2,3,6,7,10,11,15,16,19,20,23]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
    if j eq 0 then theta = theta-!dpi; - 10*!dtor 
    if j eq 1 then theta = theta + 5*!dtor 
    if j eq 2 then theta = theta + 5*!dtor 
  endif else if mu eq 5e4 then begin
    tmp = [2,3,6,7,11,12,16,17,20,21,23]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
    if j eq 1 then theta = theta + 5*!dtor 
    if j eq 2 then theta = theta - 5*!dtor 
  endif else if mu eq 1e5 then begin
    tmp = [2,3,6,7,10,11,15,16,18,19,20,23]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
    if j eq 1 then theta = theta - 10*!dtor 
    if j eq 2 then theta = theta - 10*!dtor 
  endif
  rot_xo = cos(theta)*new_xo-sin(theta)*new_yo
  rot_yo = sin(theta)*new_xo+cos(theta)*new_yo
  rot_ox_nbody = cos(theta)*(ox_nbody[i,*]-x_m)-sin(theta)*(oy_nbody[i,*]-y_m)
  rot_oy_nbody = sin(theta)*(ox_nbody[i,*]-x_m)+cos(theta)*(oy_nbody[i,*]-y_m)
  rot_ox_nbody_ex_1peri = cos(theta)*(ox_nbody_ex_1peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_1peri[i,*]-y_m)
  rot_oy_nbody_ex_1peri = sin(theta)*(ox_nbody_ex_1peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_1peri[i,*]-y_m)
  rot_ox_nbody_ex_2peri = cos(theta)*(ox_nbody_ex_2peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_2peri[i,*]-y_m)
  rot_oy_nbody_ex_2peri = sin(theta)*(ox_nbody_ex_2peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_2peri[i,*]-y_m)
  rot_ox_nbody_ex_3peri = cos(theta)*(ox_nbody_ex_3peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_3peri[i,*]-y_m)
  rot_oy_nbody_ex_3peri = sin(theta)*(ox_nbody_ex_3peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_3peri[i,*]-y_m)
  rot_ox_nbody_ex_4peri = cos(theta)*(ox_nbody_ex_4peri[i,*]-x_m)-sin(theta)*(oy_nbody_ex_4peri[i,*]-y_m)
  rot_oy_nbody_ex_4peri = sin(theta)*(ox_nbody_ex_4peri[i,*]-x_m)+cos(theta)*(oy_nbody_ex_4peri[i,*]-y_m)

;  multiplot
  tname = textoidl('M_{sat}=')+strtrim(mu,2)+'M'+odot
  plot,[0],[0],/nodata,xr=[-8.9,10.9],yr=[-5.4,2.3],/isotropic
;  plot,[0],[0],/nodata,xr=[-2.1,2.1],yr=[-0.7,0.5],/isotropic
  oplot,rot_xo,rot_yo,color=100,thick=1
  oplot,rot_ox_nbody,rot_oy_nbody,psym=3
  oplot,rot_ox_nbody_ex_4peri,rot_oy_nbody_ex_4peri,psym=3,color=230
  oplot,rot_ox_nbody_ex_3peri,rot_oy_nbody_ex_3peri,psym=3,color=210
  oplot,rot_ox_nbody_ex_2peri,rot_oy_nbody_ex_2peri,psym=3,color=150
  oplot,rot_ox_nbody_ex_1peri,rot_oy_nbody_ex_1peri,psym=3,color=70
;  oplot,new_xo,new_yo,color=220
;  oplot,new_xo,yfit,color=255
  legend,[string(snap_time[i]/1000,f='(f7.4)')+'Gyr',tname],box=0
endfor
multiplot,/reset
device,/close
erase

print,'begin sky'
device,file=dir_out+'nbody_sky.ps',/color,ysize=25
!p.multi=[0,1,3]
for i=0,N_file-1 do begin
;  m_lcosb_nbody = median(lcosb_nbody[i,*])
;  if (i eq N_file-2) then begin
;    plot,lcosb_nbody[i,*]-m_lcosb_nbody,b_nbody[i,*],psym=3,xtitle=textoidl('\Delta(lcosb) [deg]'),ytitle='b [deg]',xr=[7,-7],yr=[43.2,46.47],/isotropic
;  endif else begin
;    plot,lcosb_nbody[i,*]-m_lcosb_nbody,b_nbody[i,*],psym=3,xtitle=textoidl('\Delta(lcosb) [deg]'),ytitle='b [deg]',/isotropic
;  endelse
;  legend,[string(snap_time[i]/1000,f='(f7.4)')+'Gyr',tname],box=0
  if (snap_time[i]/1000 gt 12.600) then goto, out2
  plot,RA_nbody[i,*],Dec_nbody[i,*],xr=[RA_max,RA_min],yr=[Dec_min,Dec_max],psym=3,xtitle='RA',ytitle='Dec'
  legend,[string(snap_time[i]/1000,f='(f7.4)')+'Gyr',tname],box=0
endfor
out2:
device,/close
erase


END
