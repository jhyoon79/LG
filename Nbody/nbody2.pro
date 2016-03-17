pro nbody2

select = 2

;=== scale ===
kmpkpc = 3.240756d-17
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
  ;Msun = 1.989d33
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
mu = [1.e7,1.e6,1.e5,5.e4,1.e4]
ru = [0.075,0.0348,0.0162,0.0123,0.0075]
mu = mu[select]
ru = ru[select]

tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec

mass = strtrim(fix(alog10(mu)),2)
dir_data = '/scratch/jhyoon/Research/LG/Nbody/outer/ten'+mass+'/'
if mu eq 5.e4 then begin
  mass = '5e4'
  dir_data = '/scratch/jhyoon/Research/LG/Nbody/outer/5ten4/'
endif
dir_nosub = '/scratch/jhyoon/Research/LG/'

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

ox1 = orbit_x[1]  &  ox2 = orbit_x[2]  &  ox3 = orbit_x[3]
oy1 = orbit_y[1]  &  oy2 = orbit_y[2]  &  oy3 = orbit_y[3]
oz1 = orbit_z[1]  &  oz2 = orbit_z[2]  &  oz3 = orbit_z[3]
; normal vector of orbital plane (Ax+By+Cz+dd=0)
Ax = (oy2-oy1)*(oz3-oz1)-(oz2-oz1)*(oy3-oy1)
By = (oz2-oz1)*(ox3-ox1)-(ox2-ox1)*(oz3-oz1)
Cz = (ox2-ox1)*(oy3-oy1)-(oy2-oy1)*(ox3-ox1)
dd = -1.*(ox1*Ax+oy1*By+oz1*Cz)

; x-axis:theta    y-axis:-phi
costheta = Cz/sqrt(By^2.+Cz^2.)
sintheta = -By/sqrt(By^2.+Cz^2.)
cosphi = (-sintheta*By + costheta*Cz)/sqrt(Ax^2.+By^2.+Cz^2.)
sinphi = sqrt(1-cosphi^2.)
xx = orbit_x
yy = costheta*orbit_y + sintheta*orbit_z
zz = -sintheta*orbit_y + costheta*orbit_z
; orbit projected in orbital plane.
orbit_x3 = cosphi*xx - sinphi*zz
orbit_y3 = yy
orbit_z3 = sinphi*xx + cosphi*zz

; check whether the normal vector of the orbital plane corresponds to z-axis
Ax2 = Ax
By2 = costheta*By + sintheta*Cz
Cz2 = -sintheta*By + costheta*Cz
cosphi = (-sintheta*By + costheta*Cz)/sqrt(Ax^2.+By^2.+Cz^2.)
sinphi = sqrt(1-cosphi^2.);-Ax/sqrt(Ax^2.+By^2.+Cz^2.)
Ax3 = cosphi*Ax2 - sinphi*Cz2
By3 = By2
Cz3 = sinphi*Ax2 + cosphi*Cz2
print,Ax,By,Cz
print,Ax2,By2,Cz2
print,Ax3,By3,Cz3	; Ax3 & By3 should be zero


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
set_plot,'ps'
@plot_setting
device,file='nbody2.ps',/color
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
r_apo = orbit_r[sub[apogee]]
r_peri = orbit_r[sub[perigee]]
t_apo = orbit_t[sub[apogee]]
t_peri = orbit_t[sub[perigee]]
if mu eq 5e4 then begin
  tR = mean([(t_peri-shift(t_peri,1))[1:8],(t_apo-shift(t_apo,1))[1:7]])
endif else begin
  tR = mean([(t_peri-shift(t_peri,1))[1:4],(t_apo-shift(t_apo,1))[1:3]])
endelse

Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
r_tide = 0.108653
Msat = 20000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
chr_rdtbl, dir_nosub+'part_peri_nosubhalo',0, arr,/silent
arr = double(arr)
t_peri_mine = reform(arr[0,0])
x_peri_mine = reform(arr[1,0])
y_peri_mine = reform(arr[2,0])
z_peri_mine = reform(arr[3,0])
r_peri_mine = (sqrt(x_peri_mine^2.+y_peri_mine^2.+z_peri_mine^2.))[0]
p = r_peri_mine/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
;epsilon_mine = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri_mine)
epsilon_mine = (mu/mr)^(1./3.)*(4.3e-6*mr/r_peri_mine)
epsilon_nbody = (mu/mr)^(1./3.)*(4.3e-6*mr/r_peri_mine)
s = (mu/mr)^(1./3.)
r_tide = s*r_peri_mine

!p.multi=[0,2,2]
plot,orbit_x,orbit_y,xr=[-50,50],yr=[-50,50],/isotropic
plot,orbit_y,orbit_z,xr=[-50,50],yr=[-50,50],/isotropic
plot,orbit_x,orbit_z,xr=[-50,50],yr=[-50,50],/isotropic
plot,sqrt(orbit_x^2.+orbit_y^2.),orbit_z,xr=[0,55],yr=[-50,50],/isotropic

!p.multi=0
erase
ocolor = color_code(orbit_t)
plot,[0],[0],/nodata,xr=[-50,50],yr=[-50,50],xtitle='x',ytitle='y',/isotropic
for i=0,N_orbit-1 do oplot,[orbit_x[i]],[orbit_y[i]],psym=8,symsize=0.3,color=ocolor[i]
nn = findgen(5) + 1
t_snap = [0,N_orbit*nn/max(nn)-1]
for i=0,N_elements(t_snap)-1 do begin
  tt = string(orbit_t[t_snap[i]],f='(f5.3)')+'Gyr'
  xyouts,orbit_x[t_snap[i]],orbit_y[t_snap[i]],tt,charsize=0.6,color=ocolor[t_snap[i]]
endfor

;plot,orbit_x3,orbit_y3,xr=[-20,20],yr=[-20,20],/isotropic
;plot,orbit_y3,orbit_z3,xr=[-20,20],yr=[-20,20],/isotropic
;plot,orbit_x3,orbit_z3,xr=[-20,20],yr=[-20,20],/isotropic
;plot,[0],[0],/nodata
device,/close

N_file = N_elements(snap_fname)
dir_nosub = '/scratch/jhyoon/Research/LG/snapshot_nosubhalo_1Myr/'
snap_time = dblarr(N_file)
Npart = 3000
b_nbody = dblarr(N_file,Npart)
lcosb_nbody = dblarr(N_file,Npart)
; coord. projected in the orbital plane
ox_nbody = dblarr(N_file,Npart)
oy_nbody = dblarr(N_file,Npart)
oz_nbody = dblarr(N_file,Npart)
ox_nbody2 = dblarr(N_file,Npart)
oy_nbody2 = dblarr(N_file,Npart)
oz_nbody2 = dblarr(N_file,Npart)
ox_nbody_ex = dblarr(N_file,Npart)
oy_nbody_ex = dblarr(N_file,Npart)
oz_nbody_ex = dblarr(N_file,Npart)
ox_nbody_ex_1peri = dblarr(N_file,Npart)
oy_nbody_ex_1peri = dblarr(N_file,Npart)
oz_nbody_ex_1peri = dblarr(N_file,Npart)
ox_nbody_ex_2peri = dblarr(N_file,Npart)
oy_nbody_ex_2peri = dblarr(N_file,Npart)
oz_nbody_ex_2peri = dblarr(N_file,Npart)
ox_nbody_ex_3peri = dblarr(N_file,Npart)
oy_nbody_ex_3peri = dblarr(N_file,Npart)
oz_nbody_ex_3peri = dblarr(N_file,Npart)
Etot_nbody_cons = dblarr(N_file,Npart)
dE_cons = dblarr(N_file,Npart)

device,file='nbody2_movie.ps',/color,ysize=20,yoffset=0.5
for i=0,N_file-1 do begin
  openr,1,snap_fname[i]
  readf,1,Npart,time
  close,1
  time = round(time*tu*1000.)/1000.*1000
  fname_mine = 'snap'+string(time,f='(i4.4)')
  chr_rdtbl,dir_nosub+fname_mine,0,arr,/silent
  arr = double(arr)
  t_mine = arr[0,*]
  x_mine = arr[1,*]; / r_tide
  y_mine = arr[2,*]; / r_tide
  z_mine = arr[3,*]; / r_tide
  Vx_mine = arr[4,*]
  Vy_mine = arr[5,*]
  Vz_mine = arr[6,*]
  dE_mine = arr[7,*]
  Etot_mine = arr[8,*]
  J_mine = arr[9,*]

  chr_rdtbl,snap_fname[i],1,arr,/silent
  arr = double(arr)
  if (Npart ge 3000) then begin
    sub_sample = round(randomu(11,3000)*(N_elements(arr[0,*])-1))
    arr = arr[*,sub_sample]
  endif
  m_nbody = arr[0,*] * mu
  x_nbody = arr[1,*] * ru; / r_tide
  y_nbody = arr[2,*] * ru; / r_tide
  z_nbody = arr[3,*] * ru; / r_tide
  vx_nbody = arr[4,*] * vu
  vy_nbody = arr[5,*] * vu 
  vz_nbody = arr[6,*] * vu
  pot_int_nbody = arr[7,*] * vu^2.
  pot_ext_nbody = arr[8,*] * vu^2.
  tunbound_nbody = arr[9,*]
  bounded = where(tunbound_nbody eq 0)
  escaped = where(tunbound_nbody ne 0)
  Etot_nbody = (vx_nbody^2.+vy_nbody^2.+vz_nbody^2.)/2.+pot_ext_nbody+pot_int_nbody
  J_nbody = sqrt( (y_nbody*vz_nbody-z_nbody*vy_nbody)^2. + $
                  (z_nbody*vx_nbody-x_nbody*vz_nbody)^2. + $
                  (x_nbody*vy_nbody-y_nbody*vx_nbody)^2.)
  Emin = min(Etot_Nbody)
  Emax = max(Etot_Nbody)
  Ecolor = round((Etot_Nbody-Emin)/(Emax-Emin)*245.)+10

  l_nbody = atan(y_nbody/(x_nbody+8.))*!radeg 
  l_nbody_sub = where(l_nbody lt 0)
  snap_time[i] = time
  b_nbody[i,*] = atan(z_nbody/sqrt((x_nbody+8.)^2.+y_nbody^2.))*!radeg
  lcosb_nbody[i,*] = l_nbody*cos(b_nbody[i,*]*!dtor)
  delta_nbody = asin( cos(b_nbody[i,*]*!dtor)*sin((l_nbody-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody[i,*]*!dtor)*cos(62.6*!dtor) )*!radeg
  alpha_nbody = asin( (cos(b_nbody[i,*]*!dtor)*sin((l_nbody-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody[i,*]*!dtor)*sin(62.6*!dtor))/cos(delta_nbody*!dtor) )*!radeg + 282.25

;=== particle coord. in the orbital plane ===
  x2 = x_nbody
  y2 = costheta*y_nbody + sintheta*z_nbody
  z2 = -sintheta*y_nbody + costheta*z_nbody
  ox_nbody[i,*] = cosphi*x2 - sinphi*z2
  oy_nbody[i,*] = y2
  oz_nbody[i,*] = sinphi*x2 + cosphi*z2

  xmin = -50;min([min(x_mine),min(x_nbody)])
  xmax = 50;max([max(x_mine),max(x_nbody)])
  ymin = -50;min([min(y_mine),min(y_nbody)])
  ymax = 50;max([max(y_mine),max(y_nbody)])
  zmin = -50;min([min(z_mine),min(z_nbody)])
  zmax = 50;max([max(z_mine),max(z_nbody)])
;  xmin = -500;min([min(x_mine),min(x_nbody)])
;  xmax = 500;max([max(x_mine),max(x_nbody)])
;  ymin = -500;min([min(y_mine),min(y_nbody)])
;  ymax = 500;max([max(y_mine),max(y_nbody)])
;  zmin = -500;min([min(z_mine),min(z_nbody)])
;  zmax = 500;max([max(z_mine),max(z_nbody)])
  if (i eq 0) then begin
    m_Etot_nbody = mean(Etot_nbody)
    m_J_nbody = mean(J_nbody)
  endif

 
;=== find particle escaped at the 3rd latest perigee ===
t_apo_cri1 = t_apo[0]
t_apo_cri2 = t_apo[1]
  if (time/1000. lt t_apo_cri1) then pre_tunbound_nbody_ex_3peri = tunbound_nbody
  if (time/1000. gt t_apo_cri1 and time/1000. lt t_apo_cri2) then $
    sub_3peri = where(tunbound_nbody ne 0 and pre_tunbound_nbody_ex_3peri eq 0)
  if (time/1000. gt t_apo_cri1) then begin
    x_nbody_ex_3peri = x_nbody[sub_3peri]
    y_nbody_ex_3peri = y_nbody[sub_3peri]
    z_nbody_ex_3peri = z_nbody[sub_3peri]
    vx_nbody_ex_3peri = vx_nbody[sub_3peri]
    vy_nbody_ex_3peri = vy_nbody[sub_3peri]
    vz_nbody_ex_3peri = vz_nbody[sub_3peri]
    pot_int_nbody_ex_3peri = pot_int_nbody[sub_3peri]
    pot_ext_nbody_ex_3peri = pot_ext_nbody[sub_3peri]
    tunbound_nbody_ex_3peri = tunbound_nbody[sub_3peri]
    Etot_nbody_ex_3peri = (vx_nbody_ex_3peri^2.+vy_nbody_ex_3peri^2.+vz_nbody_ex_3peri^2.)/2.+pot_ext_nbody_ex_3peri+pot_int_nbody_ex_3peri
    J_nbody_ex_3peri = sqrt((y_nbody_ex_3peri*vz_nbody_ex_3peri-z_nbody_ex_3peri*vy_nbody_ex_3peri)^2. + $
                    (z_nbody_ex_3peri*vx_nbody_ex_3peri-x_nbody_ex_3peri*vz_nbody_ex_3peri)^2. + $
                    (x_nbody_ex_3peri*vy_nbody_ex_3peri-y_nbody_ex_3peri*vx_nbody_ex_3peri)^2.)
    dE_nbody_ex_3peri = Etot_nbody_ex_3peri-m_Etot_nbody
    dJ_nbody_ex_3peri = J_nbody_ex_3peri-m_J_nbody
    q_nbody_ex_3peri = dE_nbody_ex_3peri/epsilon_nbody
    sJ_nbody_ex_3peri = s*J_nbody_ex_3peri
;    sJ_nbody_ex_3peri = r_tide/r_peri_mine*J_nbody_ex_3peri
    dJ_sJ_nbody_ex_3peri = dJ_nbody_ex_3peri/sJ_nbody_ex_3peri

; particle coord. in the orbital plane ===
    x2 = x_nbody_ex_3peri
    y2 = costheta*y_nbody_ex_3peri + sintheta*z_nbody_ex_3peri
    z2 = -sintheta*y_nbody_ex_3peri + costheta*z_nbody_ex_3peri
    Nex = N_elements(x2)
    ox_nbody_ex_3peri[i,0:Nex-1] = cosphi*x2 - sinphi*z2
    oy_nbody_ex_3peri[i,0:Nex-1] = y2
    oz_nbody_ex_3peri[i,0:Nex-1] = sinphi*x2 + cosphi*z2
  endif


;=== find particle escaped at the 2nd latest perigee ===
t_apo_cri1 = t_apo[1]
t_apo_cri2 = t_apo[2]
  if (time/1000. lt t_apo_cri1) then pre_tunbound_nbody_ex_2peri = tunbound_nbody
  if (time/1000. gt t_apo_cri1 and time/1000. lt t_apo_cri2) then $
    sub_2peri = where(tunbound_nbody ne 0 and pre_tunbound_nbody_ex_2peri eq 0)
  if (time/1000. gt t_apo_cri1) then begin
    x_nbody_ex_2peri = x_nbody[sub_2peri]
    y_nbody_ex_2peri = y_nbody[sub_2peri]
    z_nbody_ex_2peri = z_nbody[sub_2peri]
    vx_nbody_ex_2peri = vx_nbody[sub_2peri]
    vy_nbody_ex_2peri = vy_nbody[sub_2peri]
    vz_nbody_ex_2peri = vz_nbody[sub_2peri]
    pot_int_nbody_ex_2peri = pot_int_nbody[sub_2peri]
    pot_ext_nbody_ex_2peri = pot_ext_nbody[sub_2peri]
    tunbound_nbody_ex_2peri = tunbound_nbody[sub_2peri]
    Etot_nbody_ex_2peri = (vx_nbody_ex_2peri^2.+vy_nbody_ex_2peri^2.+vz_nbody_ex_2peri^2.)/2.+pot_ext_nbody_ex_2peri+pot_int_nbody_ex_2peri
    J_nbody_ex_2peri = sqrt((y_nbody_ex_2peri*vz_nbody_ex_2peri-z_nbody_ex_2peri*vy_nbody_ex_2peri)^2. + $
                    (z_nbody_ex_2peri*vx_nbody_ex_2peri-x_nbody_ex_2peri*vz_nbody_ex_2peri)^2. + $
                    (x_nbody_ex_2peri*vy_nbody_ex_2peri-y_nbody_ex_2peri*vx_nbody_ex_2peri)^2.)
    dE_nbody_ex_2peri = Etot_nbody_ex_2peri-m_Etot_nbody
    dJ_nbody_ex_2peri = J_nbody_ex_2peri-m_J_nbody
    q_nbody_ex_2peri = dE_nbody_ex_2peri/epsilon_nbody
    sJ_nbody_ex_2peri = s*J_nbody_ex_2peri
;    sJ_nbody_ex_2peri = r_tide/r_peri_mine*J_nbody_ex_2peri
    dJ_sJ_nbody_ex_2peri = dJ_nbody_ex_2peri/sJ_nbody_ex_2peri

; particle coord. in the orbital plane ===
    x2 = x_nbody_ex_2peri
    y2 = costheta*y_nbody_ex_2peri + sintheta*z_nbody_ex_2peri
    z2 = -sintheta*y_nbody_ex_2peri + costheta*z_nbody_ex_2peri
    Nex = N_elements(x2)
    ox_nbody_ex_2peri[i,0:Nex-1] = cosphi*x2 - sinphi*z2
    oy_nbody_ex_2peri[i,0:Nex-1] = y2
    oz_nbody_ex_2peri[i,0:Nex-1] = sinphi*x2 + cosphi*z2
  endif


;=== find particle escaped at the latest perigee ===
t_apo_cri1 = t_apo[2]
t_apo_cri2 = t_apo[3]
  if (time/1000. lt t_apo_cri1) then pre_tunbound_nbody_ex_1peri = tunbound_nbody
  if (time/1000. gt t_apo_cri1 and time/1000. lt t_apo_cri2) then $
    sub_1peri = where(tunbound_nbody ne 0 and pre_tunbound_nbody_ex_1peri eq 0)
  if (time/1000. gt t_apo_cri1) then begin
    x_nbody_ex_1peri = x_nbody[sub_1peri]
    y_nbody_ex_1peri = y_nbody[sub_1peri]
    z_nbody_ex_1peri = z_nbody[sub_1peri]
    vx_nbody_ex_1peri = vx_nbody[sub_1peri]
    vy_nbody_ex_1peri = vy_nbody[sub_1peri]
    vz_nbody_ex_1peri = vz_nbody[sub_1peri]
    pot_int_nbody_ex_1peri = pot_int_nbody[sub_1peri]
    pot_ext_nbody_ex_1peri = pot_ext_nbody[sub_1peri]
    tunbound_nbody_ex_1peri = tunbound_nbody[sub_1peri]
    Etot_nbody_ex_1peri = (vx_nbody_ex_1peri^2.+vy_nbody_ex_1peri^2.+vz_nbody_ex_1peri^2.)/2.+pot_ext_nbody_ex_1peri+pot_int_nbody_ex_1peri
    J_nbody_ex_1peri = sqrt((y_nbody_ex_1peri*vz_nbody_ex_1peri-z_nbody_ex_1peri*vy_nbody_ex_1peri)^2. + $
                    (z_nbody_ex_1peri*vx_nbody_ex_1peri-x_nbody_ex_1peri*vz_nbody_ex_1peri)^2. + $
                    (x_nbody_ex_1peri*vy_nbody_ex_1peri-y_nbody_ex_1peri*vx_nbody_ex_1peri)^2.)
    dE_nbody_ex_1peri = Etot_nbody_ex_1peri-m_Etot_nbody
    dJ_nbody_ex_1peri = J_nbody_ex_1peri-m_J_nbody
    q_nbody_ex_1peri = dE_nbody_ex_1peri/epsilon_nbody
    sJ_nbody_ex_1peri = s*J_nbody_ex_1peri
;    sJ_nbody_ex_1peri = r_tide/r_peri_mine*J_nbody_ex_1peri
    dJ_sJ_nbody_ex_1peri = dJ_nbody_ex_1peri/sJ_nbody_ex_1peri

; particle coord. in the orbital plane ===
    x2 = x_nbody_ex_1peri
    y2 = costheta*y_nbody_ex_1peri + sintheta*z_nbody_ex_1peri
    z2 = -sintheta*y_nbody_ex_1peri + costheta*z_nbody_ex_1peri
    Nex = N_elements(x2)
    ox_nbody_ex_1peri[i,0:Nex-1] = cosphi*x2 - sinphi*z2
    oy_nbody_ex_1peri[i,0:Nex-1] = y2
    oz_nbody_ex_1peri[i,0:Nex-1] = sinphi*x2 + cosphi*z2
  endif


  dE_nbody = Etot_nbody-m_Etot_nbody
  q_nbody = dE_nbody/epsilon_nbody
  dJ_nbody = J_nbody-m_J_nbody
  dE_mine = Etot_mine-Etot_mine[0]
  q_mine = dE_mine/epsilon_mine

;  sJ_mine = r_tide/r_peri_mine*J_mine
;  sJ_nbody = r_tide/r_peri_mine*J_nbody
  sJ_mine = s*J_mine
  sJ_nbody = s*J_nbody
  dJ_mine = J_mine-J_mine[0]
  dJ_sJ_mine = dJ_mine/sJ_mine
  dJ_sJ_nbody = dJ_nbody/sJ_nbody

  !p.multi=[0,2,3]
  plot,x_nbody,y_nbody,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y = [kpc]',/isotropic
;  if (time/1000. gt t_apo[9]) then oplot,x_nbody2,y_nbody2,psym=8,symsize=0.3,color=100
;  if (time/1000. gt t_peri[0]) then oplot,x_nbody_ex_1peri,y_nbody_ex_1peri,psym=3,color=150
  legend,[string(time/1000,f='(f5.3)')+'Gyr'],box=0

  plot,x_mine,y_mine,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y = [kpc]',/isotropic

;  plot,x_nbody,y_nbody,psym=3,xtitle='x [kpc]',ytitle='y = [kpc]',/isotropic
;  if (time/1000. gt t_cri) then oplot,x_nbody2,y_nbody2,psym=3,color=100

  plot,y_nbody,z_nbody,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z = [kpc]',/isotropic
;  if (time/1000. gt t_cri) then oplot,y_nbody2,z_nbody2,psym=8,symsize=0.3,color=100
;  if (time/1000. gt t_peri[0]) then oplot,y_nbody_ex_1peri,z_nbody_ex_1peri,psym=3,color=150
;  plot,y_nbody,z_nbody,psym=3,xtitle='y [kpc]',ytitle='z = [kpc]',/isotropic
;  if (time/1000. gt t_cri) then oplot,y_nbody2,z_nbody2,psym=3,color=100
  plot,y_mine,z_mine,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z = [kpc]',/isotropic

  qm = max(abs([min(q_nbody),min(q_mine),max(q_nbody),max(q_mine)]))
  Jm = max(abs([min(dJ_sJ_nbody),min(dJ_sJ_mine),max(dJ_sJ_nbody),max(dJ_sJ_mine)]))
  qmin = -1.*qm  &  qmax = qm
  Jmin = -1.*Jm  &  Jmax = Jm

  plot,q_nbody,dJ_sJ_nbody,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),title='Msat='+strtrim(mu,2)+'M'+odot
  oplot,q_nbody[escaped],dJ_sJ_nbody[escaped],psym=3,color=255
;  if (time/1000. gt t_peri[0]) then oplot,q_nbody_ex_1peri,dJ_sJ_nbody_ex_1peri,psym=3,color=150
;  if (time/1000. gt t_cri) then oplot,q_nbody2,dJ_sJ_nbody2,psym=8,symsize=0.3,color=100

;  plot,q_mine,dJ_sJ_mine,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')

  levels = indgen(3000)*10 + 5
  pts_contour,q_nbody[escaped],dJ_sJ_nbody[escaped],15,15,gv,xx,yy
  contour,gv,xx,yy,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),levels=levels,xgridstyle=1,ygridstyle=1
  



;;=== E conservation check ===
;  if i eq 0 then begin
;    Etot_nbody0 = Etot_nbody
;    J_nbody0 = J_nbody
;  endif
;  Etot_nbody_cons[i,*] = Etot_nbody
;  dE_cons[i,*] = (Etot_nbody-Etot_nbody0)/Etot_nbody0
;plot,Etot_nbody0/mu,dE_cons[i,*],psym=3
;;  dJ = (J_nbody-J_nbody0)/J_nbody0*100.
;;  if i ne 0 then plothist,dJ,bin=0.05,xtitle='(J-J0)/J0 * 100'

endfor
device,/close
erase

device,file='nbody2_snap.ps',/color,ysize=16
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
    tmp = [3,4,5,6,10,11,12,13,17,18,19,20]
    if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
  ;  if i eq 2 then theta = theta - 180*!dtor 
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
  multiplot
  tname = (i eq 0) ? textoidl('M_{sat}=')+strtrim(mu,2)+'M'+odot:''
  plot,[0],[0],/nodata,xr=[-3.2,3.2],yr=[-1.6,0.9],title=tname,/isotropic
  oplot,rot_xo,rot_yo,color=100,thick=1
  oplot,rot_ox_nbody,rot_oy_nbody,psym=3
  oplot,rot_ox_nbody_ex_3peri,rot_oy_nbody_ex_3peri,psym=3,color=150
  oplot,rot_ox_nbody_ex_2peri,rot_oy_nbody_ex_2peri,psym=3,color=70
  oplot,rot_ox_nbody_ex_1peri,rot_oy_nbody_ex_1peri,psym=3,color=30
;  oplot,new_xo,new_yo,color=220
;  oplot,new_xo,yfit,color=255
  legend,[string(snap_time[i]/1000-3.2,f='(f6.3)')+'Gyr'],box=0,charsize=0.6,position=[0,-0.5]
endfor
multiplot,/reset
device,/close
erase

device,file='nbody2_current_snap.ps',/color,ysize=25
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
legend,[textoidl('t_R=')+string(tR,f='(f5.3)')+'Gyr',textoidl('t_{peri,recent}=')+string(max(t_peri),f='(f5.3)')+'Gyr'],box=0,/right

for j=0,3 do begin
  if j eq 0 then i = N_file-10
  if j eq 1 then i = N_file-7
  if j eq 2 then i = N_file-5
  if j eq 3 then i = N_file-2
;  if j eq 0 then i = N_file-18
;  if j eq 1 then i = N_file-15
;  if j eq 2 then i = N_file-12
;  if j eq 3 then i = N_file-10
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
    if j eq 0 then theta = theta - 10*!dtor 
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

;  multiplot
  tname = textoidl('M_{sat}=')+strtrim(mu,2)+'M'+odot
  plot,[0],[0],/nodata,xr=[-0.9,0.9],yr=[-0.4,0.3],/isotropic
;  plot,[0],[0],/nodata,xr=[-2.1,2.1],yr=[-0.7,0.5],/isotropic
  oplot,rot_xo,rot_yo,color=100,thick=1
  oplot,rot_ox_nbody,rot_oy_nbody,psym=3
  oplot,rot_ox_nbody_ex_3peri,rot_oy_nbody_ex_3peri,psym=3,color=210
  oplot,rot_ox_nbody_ex_2peri,rot_oy_nbody_ex_2peri,psym=3,color=150
  oplot,rot_ox_nbody_ex_1peri,rot_oy_nbody_ex_1peri,psym=3,color=70
;  oplot,new_xo,new_yo,color=220
;  oplot,new_xo,yfit,color=255
  legend,[string(snap_time[i]/1000,f='(f6.3)')+'Gyr',tname],box=0
endfor
multiplot,/reset
device,/close
erase

device,file='nbody2_sky.ps',/color,ysize=15
!p.multi=0
for i=0,N_file-1 do begin
  m_lcosb_nbody = median(lcosb_nbody[i,*])
  if (i eq N_file-2) then begin
    plot,lcosb_nbody[i,*]-m_lcosb_nbody,b_nbody[i,*],psym=3,xtitle=textoidl('\Delta(lcosb) [deg]'),ytitle='b [deg]',xr=[7,-7],yr=[43.2,46.47],/isotropic
  endif else begin
    plot,lcosb_nbody[i,*]-m_lcosb_nbody,b_nbody[i,*],psym=3,xtitle=textoidl('\Delta(lcosb) [deg]'),ytitle='b [deg]',/isotropic
  endelse
  legend,[string(snap_time[i]/1000,f='(f6.3)')+'Gyr',tname],box=0
endfor
device,/close
erase

;device,file='nbody_compare.ps',/color
;for i=0,N_file2-1 do begin
;  x_m = median(ox_nbody[i,*])
;  y_m = median(oy_nbody[i,*])
;  sub_orbit = where(abs(orbit_t-snap_time[i]/1000) lt 0.05)
;  new_xo = orbit_x3[sub_orbit]-x_m
;  new_yo = orbit_y3[sub_orbit]-y_m
;  fit = poly_fit(new_xo,new_yo,1,yfit=yfit)
;  theta = -atan(fit[1])
;  tmp = [2,3,6,7,10,11,15,16,19,20,23]
;  if ((where(tmp eq i))[0] ne -1) then theta = theta+!dpi
;  rot_xo = cos(theta)*new_xo-sin(theta)*new_yo
;  rot_yo = sin(theta)*new_xo+cos(theta)*new_yo
;  rot_ox_nbody = cos(theta)*(ox_nbody[i,*]-x_m)-sin(theta)*(oy_nbody[i,*]-y_m)
;  rot_oy_nbody = sin(theta)*(ox_nbody[i,*]-x_m)+cos(theta)*(oy_nbody[i,*]-y_m)
;  rot_ox_nbody2 = cos(theta)*(ox_nbody2[i,*]-x_m)-sin(theta)*(oy_nbody2[i,*]-y_m)
;  rot_oy_nbody2 = sin(theta)*(ox_nbody2[i,*]-x_m)+cos(theta)*(oy_nbody2[i,*]-y_m)
;  multiplot
;  tname = (i eq 0) ? 'Msat='+strtrim(mu,2)+'M'+odot:''
;
;  !p.multi=[0,1,2]
;  multiplot
;  plot,[0],[0],/nodata,xr=[-3.2,3.2],yr=[-1.6,0.9],title=tname
;  oplot,rot_xo,rot_yo,color=100,thick=1
;  oplot,rot_ox_nbody,rot_oy_nbody,psym=3
;  legend,[string(snap_time[i]/1000-3.2,f='(f6.3)')+'Gyr'],box=0,/right
;  legend,['old'],box=0
;  multiplot
;  plot,[0],[0],/nodata,xr=[-3.2,3.2],yr=[-1.6,0.9]
;  oplot,rot_xo,rot_yo,color=100,thick=1
;  oplot,rot_ox_nbody2,rot_oy_nbody2,psym=3
;  legend,['new'],box=0
;  erase
;endfor
;multiplot,/reset
;device,/close

END
