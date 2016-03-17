pro nbody

select = 4

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
;dir_data = '/scratch/jhyoon/Research/LG/Nbody_new/pal5/ten4_new/'
;dir_out = '/scratch/jhyoon/Research/LG/Nbody_new/'
dir_data = '/media/SEADISK/LG/Nbody/orb3/'
dir_out = '/media/SEADISK/LG/Nbody/'


if mu eq 5.e4 then begin
  mass = '5e4'
  dir_data = '/scratch/jhyoon/Research/LG/Nbody/pal5/5ten4/'
endif
dir_nbodynew = '/scratch/jhyoon/Research/LG/Nbody_new/'
dir_nbodynew = '/scratch/jhyoon/Research/LG/Nbody_new/'

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
;print,Ax,By,Cz
;print,Ax2,By2,Cz2
;print,Ax3,By3,Cz3	; Ax3 & By3 should be zero


chr_rdtbl,dir_data+'SCFLOG',1,arr
arr = double(arr)
log_t = arr[0,*] * tu
log_Mtot = arr[1,*] * mu
log_Mbnd = arr[2,*] * mu

spawn,'ls '+dir_data+'SNAP* > tmp5'
chr_rdtbl,'tmp5',0,arr
snap_fname = reform(arr[0,*])
spawn,'rm -f tmp5'

;chr_rdtbl,'../close_encounter/frog_5Rs_encounter_e7.dat',0,arr
;arr = double(arr)
;Mtidal_subhalo = reform(arr[6,*])
;Rs_subhalo = reform(arr[7,*])
;Rtidal_subhalo = reform(arr[8,*])

;=== plot t & r and mass loss ===
set_plot,'ps'
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
r_apo = orbit_r[sub[apogee]]
r_peri = orbit_r[sub[perigee]]
t_apo = orbit_t[sub[apogee]]
t_peri = orbit_t[sub[perigee]]
if mu eq 5e4 then begin
  tR = mean([(t_peri-shift(t_peri,1))[1:8],(t_apo-shift(t_apo,1))[1:7]])
endif else begin
  tR = mean([(t_peri-shift(t_peri,1))[1:8],(t_apo-shift(t_apo,1))[1:8]])
endelse
print,'t_apo',t_apo
help,t_apo

Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
r_tide = 0.108653
Msat = 20000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
;chr_rdtbl, dir_nbodynew+'part_peri_nosub',0, arr,/silent
chr_rdtbl, dir_nbodynew+'part_peri',0, arr,/silent
arr = double(arr)
t_peri_nosub = reform(arr[0,0])
x_peri_nosub = reform(arr[1,0])
y_peri_nosub = reform(arr[2,0])
z_peri_nosub = reform(arr[3,0])
r_peri_nosub = (sqrt(x_peri_nosub^2.+y_peri_nosub^2.+z_peri_nosub^2.))[0]
p = r_peri_nosub/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
;epsilon_nosub = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri_nosub)
epsilon_nosub = (mu/mr)^(1./3.)*(4.3e-6*mr/r_peri_nosub)
epsilon_35sub = (mu/mr)^(1./3.)*(4.3e-6*mr/r_peri_nosub)
epsilon_nbody = (mu/mr)^(1./3.)*(4.3e-6*mr/r_peri_nosub)
s = (mu/mr)^(1./3.)
r_tide = s*r_peri_nosub


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
b_nosub = dblarr(N_file,Npart)
lcosb_nosub = dblarr(N_file,Npart)
b_35sub = dblarr(N_file,Npart)
lcosb_35sub = dblarr(N_file,Npart)
alpha_nbody = dblarr(N_file,Npart)
delta_nbody = dblarr(N_file,Npart)
alpha_nosub = dblarr(N_file,Npart)
delta_nosub = dblarr(N_file,Npart)
alpha_35sub = dblarr(N_file,Npart)
delta_35sub = dblarr(N_file,Npart)
;alpha_nbody[*,*] = -999.
;delta_nbody[*,*] = -999.
alpha_nosub[*,*] = -999.
delta_nosub[*,*] = -999.
alpha_35sub[*,*] = -999.
delta_35sub[*,*] = -999.
; coord. projected in the orbital plane
ox_nbody = dblarr(N_file,Npart)
oy_nbody = dblarr(N_file,Npart)
oz_nbody = dblarr(N_file,Npart)
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
ox_nbody_ex_4peri = dblarr(N_file,Npart)
oy_nbody_ex_4peri = dblarr(N_file,Npart)
oz_nbody_ex_4peri = dblarr(N_file,Npart)
Etot_nbody_cons = dblarr(N_file,Npart)
dE_cons = dblarr(N_file,Npart)
Npts = intarr(N_file,Npart)

device,file=dir_out+'nbody_movie.ps',/color,ysize=20,yoffset=0.5
for i=0,N_file-1 do begin
  openr,1,snap_fname[i]
  readf,1,Npart,time
  close,1
  time0 = round(time*tu*1000.)/1000.*1000
  time = round(time*tu*1000.-00158)/1000.*1000
  fname_nosub = 'snap'+string(time,f='(i5.5)')
  fname_nosub_sub = 'subhalo'+string(time0,f='(i5.5)')
print,'name  ',fname_nosub
  if (time0 gt 12758) then goto, out
  chr_rdtbl,dir_nbodynew+'snapshot_nosub/'+fname_nosub,0,arr,/silent
  arr = double(arr)
  t_nosub = reform(arr[0,*])
  x_nosub = reform(arr[1,*])
  y_nosub = reform(arr[2,*])
  z_nosub = reform(arr[3,*])
  Vx_nosub = reform(arr[4,*])
  Vy_nosub = reform(arr[5,*])
  Vz_nosub = reform(arr[6,*])
  dE_nosub = reform(arr[7,*])
  Etot_nosub = reform(arr[8,*])
  J_nosub = reform(arr[9,*])
  Npts[i] = N_elements(x_nosub)
  Npts_nosub = Npts[i]
  if(i eq 0) then begin
    Npts_pre = 0
  endif else if(Npts[i-1] ne Npts[i]) then begin
    Npts_pre = Npts[i-1]
  endif 

;  chr_rdtbl,dir_nbodynew+'snapshot_35sub/'+fname_nosub,0,arr,/silent
;  arr = double(arr)
;  t_35sub = reform(arr[0,*])
 ; x_35sub = reform(arr[1,*])
;  y_35sub = reform(arr[2,*])
;  z_35sub = reform(arr[3,*])
;  Vx_35sub = reform(arr[4,*])
;  Vy_35sub = reform(arr[5,*])
 ; Vz_35sub = reform(arr[6,*])
;  dE_35sub = reform(arr[7,*])
;  Etot_35sub = reform(arr[8,*])
;  J_35sub = reform(arr[9,*])
;  Npts_35sub = N_elements(x_35sub)

  chr_rdtbl,dir_nbodynew+'snapshot_allsub/'+fname_nosub,0,arr,/silent
  arr = double(arr)
  t_allsub = reform(arr[0,*])
  x_allsub = reform(arr[1,*])
  y_allsub = reform(arr[2,*])
  z_allsub = reform(arr[3,*])
  Vx_allsub = reform(arr[4,*])
  Vy_allsub = reform(arr[5,*])
  Vz_allsub = reform(arr[6,*])
  dE_allsub = reform(arr[7,*])
  Etot_allsub = reform(arr[8,*])
  J_allsub = reform(arr[9,*])
  Npts_allsub = N_elements(x_allsub)


;  chr_rdtbl,dir_nbodynew+'snapshot_subhalo_35sub/'+fname_nosub_sub,0,arr,/silent
;  arr = double(arr)
;  t_subhalo = reform(arr[0,*])
;  x_subhalo = reform(arr[1,*])
;  y_subhalo = reform(arr[2,*])
;  z_subhalo = reform(arr[3,*])
;  Vx_subhalo = reform(arr[4,*])
;  Vy_subhalo = reform(arr[5,*])
;  Vz_subhalo = reform(arr[6,*])


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
  endif else begin
    pre_tunbound_nbody = tunbound_nbody 
    tunbound_nbody = reform(arr[9,*])
    bounded = where(tunbound_nbody eq 0)
    escaped = where(tunbound_nbody ne 0)
    recently_escaped = where(tunbound_nbody ne 0 and pre_tunbound_nbody eq 0)
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
  delta_nbody[i,*] = asin( cos(b_nbody[i,*]*!dtor)*sin((l_nbody-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody[i,*]*!dtor)*cos(62.6*!dtor) )*!radeg
  alpha_nbody[i,*] = asin( (cos(b_nbody[i,*]*!dtor)*sin((l_nbody-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody[i,*]*!dtor)*sin(62.6*!dtor))/cos(delta_nbody[i,*]*!dtor) )*!radeg + 282.25

  l_nosub = atan(y_nosub/(x_nosub+8.))*!radeg 
  b_nosub[i,0:Npts_nosub-1] = atan(z_nosub/sqrt((x_nosub+8.)^2.+y_nosub^2.))*!radeg
  lcosb_nosub[i,0:Npts_nosub-1] = l_nosub*cos(b_nosub[i,0:Npts_nosub-1]*!dtor)
  delta_nosub[i,0:Npts_nosub-1] = asin( cos(b_nosub[i,0:Npts_nosub-1]*!dtor)*sin((l_nosub-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nosub[i,0:Npts_nosub-1]*!dtor)*cos(62.6*!dtor) )*!radeg
  alpha_nosub[i,0:Npts_nosub-1] = asin( (cos(b_nosub[i,0:Npts_nosub-1]*!dtor)*sin((l_nosub-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nosub[i,0:Npts_nosub-1]*!dtor)*sin(62.6*!dtor))/cos(delta_nosub[i,0:Npts_nosub-1]*!dtor) )*!radeg + 282.25

;  l_35sub = atan(y_35sub/(x_35sub+8.))*!radeg 
;  b_35sub[i,0:Npts_35sub-1] = atan(z_35sub/sqrt((x_35sub+8.)^2.+y_35sub^2.))*!radeg
;  lcosb_35sub[i,0:Npts_35sub-1] = l_35sub*cos(b_35sub[i,0:Npts_35sub-1]*!dtor)
;  delta_35sub[i,0:Npts_35sub-1] = asin( cos(b_35sub[i,0:Npts_35sub-1]*!dtor)*sin((l_35sub-33.)*!dtor)*sin(62.6*!dtor)+sin(b_35sub[i,0:Npts_35sub-1]*!dtor)*cos(62.6*!dtor) )*!radeg
;  alpha_35sub[i,0:Npts_35sub-1] = asin( (cos(b_35sub[i,0:Npts_35sub-1]*!dtor)*sin((l_35sub-33.)*!dtor)*cos(62.6*!dtor)-sin(b_35sub[i,0:Npts_35sub-1]*!dtor)*sin(62.6*!dtor))/cos(delta_35sub[i,0:Npts_35sub-1]*!dtor) )*!radeg + 282.25


;=== particle coord. in the orbital plane ===
  x2 = x_nbody
  y2 = costheta*y_nbody + sintheta*z_nbody
  z2 = -sintheta*y_nbody + costheta*z_nbody
  ox_nbody[i,*] = cosphi*x2 - sinphi*z2
  oy_nbody[i,*] = y2
  oz_nbody[i,*] = sinphi*x2 + cosphi*z2

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
  dE_nosub = Etot_nosub-m_Etot_nbody
  q_nosub = dE_nosub/epsilon_nosub
;  dE_35sub = Etot_35sub-m_Etot_nbody
;  q_35sub = dE_35sub/epsilon_35sub
  dE_allsub = Etot_allsub-m_Etot_nbody
  q_allsub = dE_allsub/epsilon_nosub
;  sJ_nosub = r_tide/r_peri_nosub*J_nosub
;  sJ_nbody = r_tide/r_peri_nosub*J_nbody
  sJ_nbody = s*J_nbody
  dJ_nbody = J_nbody-m_J_nbody
  dJ_sJ_nbody = dJ_nbody/sJ_nbody
  sJ_nosub = s*J_nosub
  dJ_nosub = J_nosub-m_J_nbody
  dJ_sJ_nosub = dJ_nosub/sJ_nosub
;  sJ_35sub = s*J_35sub
;  dJ_35sub = J_35sub-m_J_nbody
;  dJ_sJ_35sub = dJ_35sub/sJ_35sub
  sJ_allsub = s*J_allsub
  dJ_allsub = J_allsub-m_J_nbody
  dJ_sJ_allsub = dJ_allsub/sJ_allsub


  !p.multi=[0,4,3]
  plot,x_nbody,y_nbody,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',title='Nbody',/isotropic
  legend,[string(time0/1000,f='(f7.4)')+'Gyr'],box=0,charsize=1
  plot,x_nosub,y_nosub,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',title='nosub',/isotropic
  legend,['n='+strtrim(Npts_nosub,2)],box=0,charsize=1
;  plot,x_35sub,y_35sub,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',title='35sub',/isotropic
;  for l=0L,N_elements(x_subhalo)-1 do begin
;    Rs_xy_circle = circle(x_subhalo[l],y_subhalo[l],Rs_subhalo[l])
;    oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
;  endfor
  plot,x_allsub,y_allsub,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x [kpc]',ytitle='y [kpc]',title='allsub',/isotropic


  plot,y_nbody,z_nbody,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z [kpc]',/isotropic
  plot,y_nosub,z_nosub,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z [kpc]',/isotropic
;  plot,y_35sub,z_35sub,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z [kpc]',/isotropic
;  for l=0L,N_elements(x_subhalo)-1 do begin
;    Rs_yz_circle = circle(y_subhalo[l],z_subhalo[l],Rs_subhalo[l])
;    oplot,Rs_yz_circle[0,*],Rs_yz_circle[1,*]
;  endfor
  plot,y_allsub,z_allsub,psym=3,xr=[ymin,ymax],yr=[zmin,zmax],xtitle='y [kpc]',ytitle='z [kpc]',/isotropic


;  qm = max(abs([min(q_nbody),min(q_nosub),max(q_nbody),max(q_nosub)]))
;  Jm = max(abs([min(dJ_sJ_nbody),min(dJ_sJ_nosub),max(dJ_sJ_nbody),max(dJ_sJ_nosub)]))
;  qmin = -1.*qm  &  qmax = qm
;  Jmin = -1.*Jm  &  Jmax = Jm
  qmin = -5  &  qmax = 5
  Jmin = -4.5  &  Jmax = 4.5

;  plot,q_nbody,dJ_sJ_nbody,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),title='Msat='+strtrim(mu,2)+'M'+odot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),title='Msat='+strtrim(mu,2)+'M'+odot
  oplot,q_nbody[escaped],dJ_sJ_nbody[escaped],psym=3,color=255
  if (i ge 1) then oplot,q_nbody[recently_escaped],dJ_sJ_nbody[recently_escaped],psym=3,color=70
  hline,0,linestyle=2
  vline,0,linestyle=2

  plot,q_nosub,dJ_sJ_nosub,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  hline,0,linestyle=2
  vline,0,linestyle=2
  oplot,q_nosub[Npts_pre:*],dJ_sJ_nosub[Npts_pre:*],psym=3,color=255
print,q_nosub

;  plot,q_35sub,dJ_sJ_35sub,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  hline,0,linestyle=2
  vline,0,linestyle=2

  qmin = -15  &  qmax = 15
  Jmin = -14.5  &  Jmax = 14.5
  plot,q_allsub,dJ_sJ_allsub,psym=3,xr=[qmin,qmax],yr=[Jmin,Jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  hline,0,linestyle=2
  vline,0,linestyle=2
endfor
out:device,/close
erase


device,file=dir_out+'nbody_snap.ps',/color,ysize=16
!p.multi=[0,4,7]
!p.charsize=1
;for i=0,N_file-1 do begin
for i=N_file-28,N_file-1,1 do begin
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
  sub_nosub = where(alpha_nosub[i,*] gt -999.)
;  sub_35sub = where(alpha_35sub[i,*] gt -999.)
;  alpha_max = max([reform(alpha_nbody[i,*]),reform(alpha_nosub[i,sub_nosub]),reform(alpha_35sub[i,sub_35sub])])
;  alpha_min = max([reform(alpha_nbody[i,*]),reform(alpha_nosub[i,sub_nosub]),reform(alpha_35sub[i,sub_35sub])])
;  delta_max = max([reform(delta_nbody[i,*]),reform(delta_nosub[i,sub_nosub]),reform(delta_35sub[i,sub_35sub])])
;  delta_min = max([reform(delta_nbody[i,*]),reform(delta_nosub[i,sub_nosub]),reform(delta_35sub[i,sub_35sub])])
  plot,alpha_nbody[i,*],delta_nbody[i,*],xr=[alpha_max,alpha_min],yr=[delta_min,delta_max],psym=3,xtitle='RA',ytitle='Dec'
  plot,alpha_nosub[i,sub_nosub],delta_nosub[i,sub_nosub],xr=[alpha_max,alpha_min],yr=[delta_min,delta_max],psym=3,xtitle='RA',ytitle='Dec'
;  plot,alpha_35sub[i,sub_35sub],delta_35sub[i,sub_35sub],xr=[alpha_max,alpha_min],yr=[delta_min,delta_max],psym=3,xtitle='RA',ytitle='Dec'
  legend,[string(snap_time[i]/1000,f='(f7.4)')+'Gyr',tname],box=0
endfor
out2:
device,/close
erase


END
