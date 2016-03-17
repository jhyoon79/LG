pro Etest

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
set_plot,'ps'
@plot_setting
loadct,13
device, file='Etest.ps',/color,ysize=22
!p.multi=[0,1,2]

for i=1,34 do begin
t_final = 8437-(34-i)*20
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

ff = 'SNAP'+string(i,f='(i3.3)')
  openr,1,dir_pal5_nbody+ff
  readf,1,Npart,time
  close,1
  time = time*tu
print,'t=',time

chr_rdtbl,dir_pal5_nbody+ff,1,arr,/silent
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

s = (mu/mr)^(1./3.)
m_Etot_nbody = mean(Etot_nbody)
dE_nbody = Etot_nbody-m_Etot_nbody
q_nbody = dE_nbody/epsilon_nbody
sJ_nbody = s*J_nbody
m_J_nbody = mean(J_nbody)
dJ_nbody = J_nbody-m_J_nbody
dJ_sJ_nbody = dJ_nbody/sJ_nbody


symsize=0.1
gray = 150
;=== Energy plot ===
delE = E_total-E_total[0]
q = (delE)/epsilon
dJ = J-J[0]
scale_dJ = dJ / (s * J)

;=== comparison between energy & angular momentum ===
tmp1 = where(q gt 2 and q lt 2.5)
tmp2 = where(q gt -2.5 and q lt -2)
tmp3 = where(q gt 2 and q lt 2.5 and scale_dJ gt 2.0)
tmp4 = where(q gt 2 and q lt 2.5 and scale_dJ lt 2.0 and scale_dJ gt 1.5)
plot,x,y,psym=3
oplot,x[tmp1],y[tmp1],psym=8,symsize=0.3,color=100
oplot,x[tmp2],y[tmp2],psym=8,symsize=0.3,color=255
oplot,x[tmp3],y[tmp3],psym=8,symsize=0.3,color=70
oplot,x[tmp4],y[tmp4],psym=8,symsize=0.3,color=210

if i eq 2 then begin
;tm1 = where(q_nbody gt 2 and q_nbody lt 2.2 and dJ_sJ_nbody gt 1)
;tm2 = where(q_nbody gt -2.2 and q_nbody lt -2 and dJ_sJ_nbody lt -1)
tm1 = where(dJ_sJ_nbody gt 1.2 and dJ_sJ_nbody lt 1.5)
tm2 = where(dJ_sJ_nbody gt 1.5 and dJ_sJ_nbody lt 1.8)
plot,q_nbody,dJ_sJ_nbody,psym=3
endif

plot,x_nbody,y_nbody,psym=3
if i ge 2 then begin
oplot,x_nbody[tm1],y_nbody[tm1],psym=8,symsize=0.3,color=100
oplot,x_nbody[tm2],y_nbody[tm2],psym=8,symsize=0.3,color=255
endif




endfor

device, /close

END
