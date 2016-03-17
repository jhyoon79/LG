pro figure1

dir_Pal5 = '/media/SEADISK/LG/realistic/'
;dir_out = '~/Research/LG/Report/'
dir_out = './'
;dir_pal5_nbody = '/media/SEADISK/LG/Nbody/pal5/ten4/neworb/'
dir_pal5_nbody = '/media/SEADISK/LG/Nbody/pal5_final/'

openr, 99, 'frogin'
readf, 99, N_stars
close, 99
openr, 99, 'frogin_initial_Pal5'
readf, 99, tmp
readf, 99, t_final
close, 99
t_now = t_final/1000.

openr,99,'Pal5_sigma.dat'
readf,99,sigma_v,stddev_v
close,99

;=== read output of the first particle ===
chr_rdtbl,dir_Pal5+'part001',0,arr,/silent
arr = double(arr)
t_part1 = arr[0,*]/1000.
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

;=== the orbit of the central particle ===
;chr_rdtbl,dir_Pal5+'snapshot_nosub/snap12758',0,arr,/silent
fname = 'snapshot/snap'+string(t_final,f='(i5.5)')
chr_rdtbl,dir_Pal5+fname,0,arr,/silent
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

l = atan(y/(x+8.))*!radeg 
l_sub = where(l lt 0)
;if l_sub[0] ne -1 then l[l_sub] = l[l_sub] + 360.
b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
lcosb = l*cos(b*!dtor)
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,dir_Pal5+'part_peri',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
V_peri = sqrt(Vx_peri^2.+Vy_peri^2.+Vz_peri^2.)
print,'Rperi',r_peri
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
Msat = 10000.
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
s = (Msat/mr)^(1./3.)
epsilon = s*(G*mr/r_peri)

;=== N-body ===
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
mu = 1.e4
ru = 0.0075
tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec
epsilon_nbody = (mu/mr)^(1./3.)*(4.3e-6*mr/r_peri)

chr_rdtbl,dir_pal5_nbody+'SCFCEN',1,arr,/silent
arr = double(arr)
orbit_pal5_t = arr[0,*] * tu
orbit_pal5_dt = arr[1,*] * tu
orbit_pal5_x = arr[2,*] * ru
orbit_pal5_y = arr[3,*] * ru
orbit_pal5_z = arr[4,*] * ru
orbit_pal5_vx = arr[5,*] * vu
orbit_pal5_vy = arr[6,*] * vu
orbit_pal5_vz = arr[7,*] * vu
orbit_pal5_r = sqrt(orbit_pal5_x^2.+orbit_pal5_y^2.+orbit_pal5_z^2.)
orbit_pal5_rxy = sqrt(orbit_pal5_x^2.+orbit_pal5_y^2.)
sub_t_nbody_pal5 = where(orbit_pal5_t gt t_now-.4 and orbit_pal5_t lt t_now+.4)

N_orbit_pal5 = N_elements(orbit_pal5_t)
dr1 = shift(orbit_pal5_r,1)-orbit_pal5_r
dr2 = orbit_pal5_r - shift(orbit_pal5_r,-1)
dr1[0] = 0
dr2[N_orbit_pal5-1] = 0
sub = where(dr1*dr2 lt 0)
m_orbit_pal5_r = mean(orbit_pal5_r[sub])
apogee = where(orbit_pal5_r[sub] gt m_orbit_pal5_r,complement=perigee)
r_apo_pal5 = orbit_pal5_r[sub[apogee]]
r_peri_pal5 = orbit_pal5_r[sub[perigee]]
t_apo_pal5 = orbit_pal5_t[sub[apogee]]
t_peri_pal5 = orbit_pal5_t[sub[perigee]]

  openr,1,dir_pal5_nbody+'SNAP034'
  readf,1,Npart,time
  close,1
  time = time*tu
print,'t=',time

chr_rdtbl,dir_pal5_nbody+'SNAP034',1,arr,/silent
arr = double(arr)
x_nbody = reform(arr[1,*])*ru
y_nbody = reform(arr[2,*])*ru
z_nbody = reform(arr[3,*])*ru
vx_nbody = reform(arr[4,*])*vu
vy_nbody = reform(arr[5,*])*vu
vz_nbody = reform(arr[6,*])*vu
pot_int_nbody = reform(arr[7,*])*vu^2.
pot_ext_nbody = reform(arr[8,*])*vu^2.
v_nbody = sqrt(Vx_nbody^2.+Vy_nbody^2.+Vz_nbody^2.)
tunbound_nbody = reform(arr[9,*])
Etot_nbody = v_nbody^2./2.+pot_ext_nbody+pot_int_nbody
;bounded = where(tunbound_nbody eq 0)
;escaped = where(tunbound_nbody ne 0)
J_nbody = sqrt( (y_nbody*vz_nbody-z_nbody*vy_nbody)^2. + $
                (z_nbody*vx_nbody-x_nbody*vz_nbody)^2. + $
                (x_nbody*vy_nbody-y_nbody*vx_nbody)^2.)
Emin = min(Etot_nbody)
Emax = max(Etot_nbody)

l_nbody = atan(y_nbody/(x_nbody+8.))*!radeg 
b_nbody = atan(z_nbody/sqrt((x_nbody+8.)^2.+y_nbody^2.))*!radeg
lcosb_nbody = l_nbody*cos(b_nbody*!dtor)
delta_nbody = asin( cos(b_nbody*!dtor)*sin((l_nbody-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_nbody = asin( (cos(b_nbody*!dtor)*sin((l_nbody-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody*!dtor)*sin(62.6*!dtor))/cos(delta_nbody*!dtor) )*!radeg + 282.25

l_nbody_orbit_pal5 = atan(orbit_pal5_y/(orbit_pal5_x+8.))*!radeg 
b_nbody_orbit_pal5 = atan(orbit_pal5_z/sqrt((orbit_pal5_x+8.)^2.+orbit_pal5_y^2.))*!radeg
lcosb_nbody_orbit_pal5 = l_nbody_orbit_pal5*cos(b_nbody_orbit_pal5*!dtor)
delta_nbody_orbit_pal5 = asin( cos(b_nbody_orbit_pal5*!dtor)*sin((l_nbody_orbit_pal5-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody_orbit_pal5*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_nbody_orbit_pal5 = asin( (cos(b_nbody_orbit_pal5*!dtor)*sin((l_nbody_orbit_pal5-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody_orbit_pal5*!dtor)*sin(62.6*!dtor))/cos(delta_nbody_orbit_pal5*!dtor) )*!radeg + 282.25

s = (mu/mr)^(1./3.)
m_Etot_nbody = mean(Etot_nbody)
dE_nbody = Etot_nbody-m_Etot_nbody
q_nbody = dE_nbody/epsilon_nbody
sJ_nbody = s*J_nbody
m_J_nbody = mean(J_nbody)
dJ_nbody = J_nbody-m_J_nbody
dJ_sJ_nbody = dJ_nbody/sJ_nbody


;=== estimate equatorial coord. ===
delta = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
alpha = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(delta*!dtor) )*!radeg + 282.25
alpha = reform(alpha)
delta = reform(delta)
alpha_sub = where(alpha lt 0)
if alpha_sub[0] ne -1 then alpha[alpha_sub] = alpha[alpha_sub] + 360.
alpha_sub2 = where(alpha gt 360)
if alpha_sub2[0] ne -1 then alpha[alpha_sub2] = alpha[alpha_sub2] - 360.

;sub_t = where(t_part1 gt 12.658 and t_part1 lt 12.958)
sub_t = where(t_part1 gt t_now-.4 and t_part1 lt t_now+.4)
l_part1 = atan(y_part1/(x_part1+8.))*!radeg 
l_part1_sub = where(l_part1 lt 0)
l_part1_sub2 = where(l_part1 gt 360)
b_part1 = atan(z_part1/sqrt((x_part1+8.)^2.+y_part1^2.))*!radeg
lcosb_part1 = l_part1*cos(b_part1*!dtor)


delta_part1 = asin( cos(b_part1*!dtor)*sin((l_part1-33.)*!dtor)*sin(62.6*!dtor)+sin(b_part1*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_part1 = asin( (cos(b_part1*!dtor)*sin((l_part1-33.)*!dtor)*cos(62.6*!dtor)-sin(b_part1*!dtor)*sin(62.6*!dtor))/cos(delta_part1*!dtor) )*!radeg + 282.25
alpha_part1 = reform(alpha_part1)
delta_part1 = reform(delta_part1)
alpha_part1_sub = where(alpha_part1 lt 0)
if alpha_part1_sub[0] ne -1 then alpha_part1[alpha_part1_sub] = alpha_part1[alpha_part1_sub] + 360.
alpha_part1_sub2 = where(alpha_part1 gt 360)
if alpha_part1_sub2[0] ne -1 then alpha_part1[alpha_part1_sub2] = alpha_part1[alpha_part1_sub2] - 360.
sub_b = sub_t
alpha_Pal5 = alpha[0]
delta_Pal5 = delta[0]
l_Pal5 = l[0]
lcosb_Pal5 = lcosb[0]
b_Pal5 = b[0]

set_plot,'ps'
@plot_setting
;!p.charsize=1.2
loadct,13
device, file=dir_out+'figure1.eps',/color,/landscape,/enc,/cmyk,yoffset=24
!p.multi=[0,2,2]

symsize=0.1
gray = 150
;=== Energy plot ===
delE = E_total-E_total[0]
q = (delE)/epsilon
dJ = J-J[0]
scale_dJ = dJ / (s * J)

;=== comparison between energy & angular momentum ===
qmin = 4  &  qmax = -4
dJmin = -4  &  dJmax = 4
plot,q_nbody,dJ_sJ_nbody,psym=8,symsize=symsize,xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),xr=[qmin,qmax],yr=[dJmin,dJmax]
;oplot,q,scale_dJ,psym=3,color=255
vline,0,linestyle=2
hline,0,linestyle=2

plot,[0],[0],/nodata,xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),xr=[qmin,qmax],yr=[dJmin,dJmax]
oplot,q_nbody,dJ_sJ_nbody,psym=8,symsize=symsize,color=gray
oplot,q,scale_dJ,psym=8,symsize=symsize
tmp1 = where(q gt 2 and q lt 2.5)
tmp2 = where(q gt -2.5 and q lt -2)
tmp3 = where(q gt 2 and q lt 2.5 and scale_dJ gt 2.3)
oplot,q[tmp1],scale_dJ[tmp1],psym=8,symsize=.3,color=100
oplot,q[tmp2],scale_dJ[tmp2],psym=8,symsize=.3,color=255
oplot,q[tmp3],scale_dJ[tmp3],psym=8,symsize=.3,color=70
vline,0,linestyle=2
hline,0,linestyle=2

;=== Pal5 in equatorial coord. ===
xmin = 247  &  xmax = 216
ymin = -17  &  ymax = 10
plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='RA [deg]',ytitle='Dec [deg]'
oplot,alpha_part1[sub_t],delta_part1[sub_t],linestyle=2,color=150
oplot,alpha_nbody,delta_nbody,psym=8,symsize=symsize
oplot,[alpha_Pal5],[delta_Pal5],psym=1,symsize=2,thick=6
arrow,227,2,223,-2,/data
;oplot,alpha_nbody_orbit_pal5[sub_t_nbody_pal5],delta_nbody_orbit_pal5[sub_t_nbody_pal5],color=100


plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='RA [deg]',ytitle='Dec [deg]'
oplot,alpha_part1[sub_t],delta_part1[sub_t],linestyle=2,color=150
oplot,alpha_nbody,delta_nbody,psym=8,color=gray,symsize=symsize
oplot,alpha,delta,psym=8,symsize=symsize
oplot,[alpha_Pal5],[delta_Pal5],psym=1,symsize=2,thick=6
oplot,alpha[tmp1],delta[tmp1],psym=8,symsize=0.3,color=100
oplot,alpha[tmp2],delta[tmp2],psym=8,symsize=0.3,color=255
oplot,alpha[tmp3],delta[tmp3],psym=8,symsize=0.3,color=70
;legend,[textoidl('\sigma_v=')+string(sigma_v,f='(f5.3)'),$
;        textoidl('\sigma_{v2}=')+string(stddev_v,f='(f5.3)')],box=0,/right

device, /close

END
