pro one_encounter_snap

dir_out = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_nosub = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_nosub/'
dir_part_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy/'
dir_subhalo_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy/'
;dir_part_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v200/'
;dir_subhalo_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v200/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
Msat = 10000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,dir_out+'part_peri_xy_nosub',0, arr
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
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s = (Msat/mr)^(1./3.)
epsilon = s*(4.3e-6*mr/r_peri)

;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part_v200+'snap00000',0, arr
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (s * J_initial)
qmin = min(q_initial)*0.9
qmax = max(q_initial)*0.9
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)
qmin = -4
qmax = 4
dJmin = -4
dJmax = 4

;=== read out simulation info ===
chr_rdtbl,'orbit_subhalo_e7.dat',0,arr
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)
N_time = N_elements(fname_part_v200)


;=== v200 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_v200+'snap00480',0,arr,/silent
arr = double(arr)
t_part_v200 = reform(arr[0,*])
x_part_v200 = reform(arr[1,*])
y_part_v200 = reform(arr[2,*])
z_part_v200 = reform(arr[3,*])
Vx_part_v200 = reform(arr[4,*])
Vy_part_v200 = reform(arr[5,*])
Vz_part_v200 = reform(arr[6,*])
dE_part_v200 = reform(arr[7,*])
Etot_part_v200 = reform(arr[8,*])
J_part_v200 = reform(arr[9,*])
delE = Etot_part_v200-Etot_part_v200[0]
q_part_v200 = delE/epsilon
dJ_part_v200 = J_part_v200-J_part_v200[0]
scale_dJ_part_v200 = dJ_part_v200 / (s * J_part_v200)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v200+'subhalo00480',0,arr,/silent
arr = double(arr)
x_subhalo_v200 = reform(arr[1,*])
y_subhalo_v200 = reform(arr[2,*])
z_subhalo_v200 = reform(arr[3,*])
Vx_subhalo_v200 = reform(arr[4,*])
Vy_subhalo_v200 = reform(arr[5,*])
Vz_subhalo_v200 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v200 = circle(x_subhalo_v200[0],y_subhalo_v200[0],double(Rs[0]))
Rtidal_xy_circle_v200 = circle(x_subhalo_v200[0],y_subhalo_v200[0],double(Rtidal[0]))
Rs_yz_circle_v200 = circle(y_subhalo_v200[0],z_subhalo_v200[0],double(Rs[0]))
Rtidal_yz_circle_v200 = circle(y_subhalo_v200[0],z_subhalo_v200[0],double(Rtidal[0]))
Rs_xz_circle_v200 = circle(x_subhalo_v200[0],z_subhalo_v200[0],double(Rs[0]))
Rtidal_xz_circle_v200 = circle(x_subhalo_v200[0],z_subhalo_v200[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_v200+'snap00500',0,arr,/silent
arr = double(arr)
t_part_v200_2 = reform(arr[0,*])
x_part_v200_2 = reform(arr[1,*])
y_part_v200_2 = reform(arr[2,*])
z_part_v200_2 = reform(arr[3,*])
Vx_part_v200_2 = reform(arr[4,*])
Vy_part_v200_2 = reform(arr[5,*])
Vz_part_v200_2 = reform(arr[6,*])
dE_part_v200_2 = reform(arr[7,*])
Etot_part_v200_2 = reform(arr[8,*])
J_part_v200_2 = reform(arr[9,*])
delE = Etot_part_v200_2-Etot_part_v200_2[0]
q_part_v200_2 = delE/epsilon
dJ_part_v200_2 = J_part_v200_2-J_part_v200_2[0]
scale_dJ_part_v200_2 = dJ_part_v200_2 / (s * J_part_v200_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v200+'subhalo00500',0,arr,/silent
arr = double(arr)
x_subhalo_v200_2 = reform(arr[1,*])
y_subhalo_v200_2 = reform(arr[2,*])
z_subhalo_v200_2 = reform(arr[3,*])
Vx_subhalo_v200_2 = reform(arr[4,*])
Vy_subhalo_v200_2 = reform(arr[5,*])
Vz_subhalo_v200_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v200_2 = circle(x_subhalo_v200_2[0],y_subhalo_v200_2[0],double(Rs[0]))
Rtidal_xy_circle_v200_2 = circle(x_subhalo_v200_2[0],y_subhalo_v200_2[0],double(Rtidal[0]))
Rs_yz_circle_v200_2 = circle(y_subhalo_v200_2[0],z_subhalo_v200_2[0],double(Rs[0]))
Rtidal_yz_circle_v200_2 = circle(y_subhalo_v200_2[0],z_subhalo_v200_2[0],double(Rtidal[0]))
Rs_xz_circle_v200_2 = circle(x_subhalo_v200_2[0],z_subhalo_v200_2[0],double(Rs[0]))
Rtidal_xz_circle_v200_2 = circle(x_subhalo_v200_2[0],z_subhalo_v200_2[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_v200+'snap00520',0,arr,/silent
arr = double(arr)
t_part_v200_3 = reform(arr[0,*])
x_part_v200_3 = reform(arr[1,*])
y_part_v200_3 = reform(arr[2,*])
z_part_v200_3 = reform(arr[3,*])
Vx_part_v200_3 = reform(arr[4,*])
Vy_part_v200_3 = reform(arr[5,*])
Vz_part_v200_3 = reform(arr[6,*])
dE_part_v200_3 = reform(arr[7,*])
Etot_part_v200_3 = reform(arr[8,*])
J_part_v200_3 = reform(arr[9,*])
delE = Etot_part_v200_3-Etot_part_v200_3[0]
q_part_v200_3 = delE/epsilon
dJ_part_v200_3 = J_part_v200_3-J_part_v200_3[0]
scale_dJ_part_v200_3 = dJ_part_v200_3 / (s * J_part_v200_3)

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_v200+'subhalo00520',0,arr,/silent
arr = double(arr)
x_subhalo_v200_3 = reform(arr[1,*])
y_subhalo_v200_3 = reform(arr[2,*])
z_subhalo_v200_3 = reform(arr[3,*])
Vx_subhalo_v200_3 = reform(arr[4,*])
Vy_subhalo_v200_3 = reform(arr[5,*])
Vz_subhalo_v200_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v200_3 = circle(x_subhalo_v200_3[0],y_subhalo_v200_3[0],double(Rs[0]))
Rtidal_xy_circle_v200_3 = circle(x_subhalo_v200_3[0],y_subhalo_v200_3[0],double(Rtidal[0]))
Rs_yz_circle_v200_3 = circle(y_subhalo_v200_3[0],z_subhalo_v200_3[0],double(Rs[0]))
Rtidal_yz_circle_v200_3 = circle(y_subhalo_v200_3[0],z_subhalo_v200_3[0],double(Rtidal[0]))
Rs_xz_circle_v200_3 = circle(x_subhalo_v200_3[0],z_subhalo_v200_3[0],double(Rs[0]))
Rtidal_xz_circle_v200_3 = circle(x_subhalo_v200_3[0],z_subhalo_v200_3[0],double(Rtidal[0]))


;============
;=== plot ===
;============
set_plot,'ps'
@plot_setting
!p.charsize=1.8
!p.multi=[0,3,3]
file_out = dir_out+'one_encounter_xy_snap1.ps'
device,file=file_out,/color,/landscape;,ysize=25;/landscape

;=== 20Myr earlier than the closest approach ===
plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_v200[k]-x_part_v200[0]],[z_part_v200[k]-z_part_v200[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_v200[0,*]-x_part_v200[0],Rs_xz_circle_v200[1,*]-z_part_v200[0]
oplot,Rtidal_xz_circle_v200[0,*]-x_part_v200[0],Rtidal_xz_circle_v200[1,*]-z_part_v200[0],linestyle=2
legend,['dt=-20Myr'],box=0,charsize=0.8
arrow,-6.6,3,-6.6,-3,/data
xyouts,-6.9,-5,textoidl('V_z=-200km/s'),charsize=0.8

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_v200[k]],[scale_dJ_part_v200[k]],psym=3,color=qcolor[k]

; q histogram
plothist,q_part_v200,xr=[qmin,qmax],bin=0.1,yr=[0,150],xtitle='q',ytitle='N'

;=== at the moment of closest approach, 10Gyr ===
plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_v200_2[k]-x_part_v200_2[0]],[z_part_v200_2[k]-z_part_v200_2[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_v200_2[0,*]-x_part_v200_2[0],Rs_xz_circle_v200_2[1,*]-z_part_v200_2[0]
oplot,Rtidal_xz_circle_v200_2[0,*]-x_part_v200_2[0],Rtidal_xz_circle_v200_2[1,*]-z_part_v200_2[0],linestyle=2
legend,['dt=0Myr'],box=0,charsize=0.8

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_v200_2[k]],[scale_dJ_part_v200_2[k]],psym=3,color=qcolor[k]

; q histogram
plothist,q_part_v200_2,xr=[qmin,qmax],bin=0.1,yr=[0,150],xtitle='q',ytitle='N'

;=== 20Myr later than the closest approach ===
plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_v200_3[k]-x_part_v200_3[0]],[z_part_v200_3[k]-z_part_v200_3[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_v200_3[0,*]-x_part_v200_3[0],Rs_xz_circle_v200_3[1,*]-z_part_v200_3[0]
oplot,Rtidal_xz_circle_v200_3[0,*]-x_part_v200_3[0],Rtidal_xz_circle_v200_3[1,*]-z_part_v200_3[0],linestyle=2
legend,['dt=20Myr'],box=0,charsize=0.8

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_v200_3[k]],[scale_dJ_part_v200_3[k]],psym=3,color=qcolor[k]

; q histogram
plothist,q_part_v200_3,xr=[qmin,qmax],bin=0.1,yr=[0,150],xtitle='q',ytitle='N'




;=====================================
;=== encounter pass along the tail ===
;=====================================
dir_part_VxVy = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_VxVy/'
dir_subhalo_VxVy = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_VxVy/'


;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part_VxVy+'snap00000',0, arr
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (s * J_initial)
qmin = min(q_initial)*0.9
qmax = max(q_initial)*0.9
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)
qmin = -4
qmax = 4
dJmin = -4
dJmax = 4


;=== VxVy ===
; read out data of the tail particles 
chr_rdtbl,dir_part_VxVy+'snap00960',0,arr,/silent
arr = double(arr)
t_part_VxVy = reform(arr[0,*])
x_part_VxVy = reform(arr[1,*])
y_part_VxVy = reform(arr[2,*])
z_part_VxVy = reform(arr[3,*])
Vx_part_VxVy = reform(arr[4,*])
Vy_part_VxVy = reform(arr[5,*])
Vz_part_VxVy = reform(arr[6,*])
dE_part_VxVy = reform(arr[7,*])
Etot_part_VxVy = reform(arr[8,*])
J_part_VxVy = reform(arr[9,*])
delE = Etot_part_VxVy-Etot_part_VxVy[0]
q_part_VxVy = delE/epsilon
dJ_part_VxVy = J_part_VxVy-J_part_VxVy[0]
scale_dJ_part_VxVy = dJ_part_VxVy / (s * J_part_VxVy)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_VxVy+'subhalo00980',0,arr,/silent
arr = double(arr)
x_subhalo_VxVy = reform(arr[1,*])
y_subhalo_VxVy = reform(arr[2,*])
z_subhalo_VxVy = reform(arr[3,*])
Vx_subhalo_VxVy = reform(arr[4,*])
Vy_subhalo_VxVy = reform(arr[5,*])
Vz_subhalo_VxVy = reform(arr[6,*])
print,sqrt(Vx_subhalo_VxVy^2.+Vy_subhalo_VxVy^2.+Vz_subhalo_VxVy^2.)

; make a circle of a subhalo 
Rs_xy_circle_VxVy = circle(x_subhalo_VxVy[0],y_subhalo_VxVy[0],double(Rs[0]))
Rtidal_xy_circle_VxVy = circle(x_subhalo_VxVy[0],y_subhalo_VxVy[0],double(Rtidal[0]))
Rs_yz_circle_VxVy = circle(y_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rs[0]))
Rtidal_yz_circle_VxVy = circle(y_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rtidal[0]))
Rs_xz_circle_VxVy = circle(x_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rs[0]))
Rtidal_xz_circle_VxVy = circle(x_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_VxVy+'snap01000',0,arr,/silent
arr = double(arr)
t_part_VxVy_2 = reform(arr[0,*])
x_part_VxVy_2 = reform(arr[1,*])
y_part_VxVy_2 = reform(arr[2,*])
z_part_VxVy_2 = reform(arr[3,*])
Vx_part_VxVy_2 = reform(arr[4,*])
Vy_part_VxVy_2 = reform(arr[5,*])
Vz_part_VxVy_2 = reform(arr[6,*])
dE_part_VxVy_2 = reform(arr[7,*])
Etot_part_VxVy_2 = reform(arr[8,*])
J_part_VxVy_2 = reform(arr[9,*])
delE = Etot_part_VxVy_2-Etot_part_VxVy_2[0]
q_part_VxVy_2 = delE/epsilon
dJ_part_VxVy_2 = J_part_VxVy_2-J_part_VxVy_2[0]
scale_dJ_part_VxVy_2 = dJ_part_VxVy_2 / (s * J_part_VxVy_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_VxVy+'subhalo00500',0,arr,/silent
arr = double(arr)
x_subhalo_VxVy_2 = reform(arr[1,*])
y_subhalo_VxVy_2 = reform(arr[2,*])
z_subhalo_VxVy_2 = reform(arr[3,*])
Vx_subhalo_VxVy_2 = reform(arr[4,*])
Vy_subhalo_VxVy_2 = reform(arr[5,*])
Vz_subhalo_VxVy_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_VxVy_2 = circle(x_subhalo_VxVy_2[0],y_subhalo_VxVy_2[0],double(Rs[0]))
Rtidal_xy_circle_VxVy_2 = circle(x_subhalo_VxVy_2[0],y_subhalo_VxVy_2[0],double(Rtidal[0]))
Rs_yz_circle_VxVy_2 = circle(y_subhalo_VxVy_2[0],z_subhalo_VxVy_2[0],double(Rs[0]))
Rtidal_yz_circle_VxVy_2 = circle(y_subhalo_VxVy_2[0],z_subhalo_VxVy_2[0],double(Rtidal[0]))
Rs_xz_circle_VxVy_2 = circle(x_subhalo_VxVy_2[0],z_subhalo_VxVy_2[0],double(Rs[0]))
Rtidal_xz_circle_VxVy_2 = circle(x_subhalo_VxVy_2[0],z_subhalo_VxVy_2[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_VxVy+'snap01040',0,arr,/silent
arr = double(arr)
t_part_VxVy_3 = reform(arr[0,*])
x_part_VxVy_3 = reform(arr[1,*])
y_part_VxVy_3 = reform(arr[2,*])
z_part_VxVy_3 = reform(arr[3,*])
Vx_part_VxVy_3 = reform(arr[4,*])
Vy_part_VxVy_3 = reform(arr[5,*])
Vz_part_VxVy_3 = reform(arr[6,*])
dE_part_VxVy_3 = reform(arr[7,*])
Etot_part_VxVy_3 = reform(arr[8,*])
J_part_VxVy_3 = reform(arr[9,*])
delE = Etot_part_VxVy_3-Etot_part_VxVy_3[0]
q_part_VxVy_3 = delE/epsilon
dJ_part_VxVy_3 = J_part_VxVy_3-J_part_VxVy_3[0]
scale_dJ_part_VxVy_3 = dJ_part_VxVy_3 / (s * J_part_VxVy_3)

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_VxVy+'subhalo01040',0,arr,/silent
arr = double(arr)
x_subhalo_VxVy_3 = reform(arr[1,*])
y_subhalo_VxVy_3 = reform(arr[2,*])
z_subhalo_VxVy_3 = reform(arr[3,*])
Vx_subhalo_VxVy_3 = reform(arr[4,*])
Vy_subhalo_VxVy_3 = reform(arr[5,*])
Vz_subhalo_VxVy_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_VxVy_3 = circle(x_subhalo_VxVy_3[0],y_subhalo_VxVy_3[0],double(Rs[0]))
Rtidal_xy_circle_VxVy_3 = circle(x_subhalo_VxVy_3[0],y_subhalo_VxVy_3[0],double(Rtidal[0]))
Rs_yz_circle_VxVy_3 = circle(y_subhalo_VxVy_3[0],z_subhalo_VxVy_3[0],double(Rs[0]))
Rtidal_yz_circle_VxVy_3 = circle(y_subhalo_VxVy_3[0],z_subhalo_VxVy_3[0],double(Rtidal[0]))
Rs_xz_circle_VxVy_3 = circle(x_subhalo_VxVy_3[0],z_subhalo_VxVy_3[0],double(Rs[0]))
Rtidal_xz_circle_VxVy_3 = circle(x_subhalo_VxVy_3[0],z_subhalo_VxVy_3[0],double(Rtidal[0]))


;============
;=== plot ===
;============
!p.multi=[0,3,3]
file_out = dir_out+'one_encounter_xy_snap2.ps'
device,file=file_out,/color;,ysize=25;/landscape
;=== 40Myr earlier than the closest approach ===
plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_VxVy[k]-x_part_VxVy[0]],[y_part_VxVy[k]-y_part_VxVy[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_VxVy[0,*]-x_part_VxVy[0],Rs_xy_circle_VxVy[1,*]-y_part_VxVy[0]
oplot,Rtidal_xy_circle_VxVy[0,*]-x_part_VxVy[0],Rtidal_xy_circle_VxVy[1,*]-y_part_VxVy[0],linestyle=2
legend,['dt=-40Myr'],box=0,charsize=0.8
arrow,-10.,4,-5,0,/data
xyouts,-12,-4,textoidl('V=200km/s'),charsize=0.8

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_VxVy[k]],[scale_dJ_part_VxVy[k]],psym=3,color=qcolor[k]

; q histogram
plothist,q_part_VxVy,xr=[qmin,qmax],bin=0.1,yr=[0,290],xtitle='q',ytitle='N'

;=== at the moment of closest approach, 10Gyr ===
plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_VxVy_2[k]-x_part_VxVy_2[0]],[y_part_VxVy_2[k]-y_part_VxVy_2[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_VxVy_2[0,*]-x_part_VxVy_2[0],Rs_xy_circle_VxVy_2[1,*]-y_part_VxVy_2[0]
oplot,Rtidal_xy_circle_VxVy_2[0,*]-x_part_VxVy_2[0],Rtidal_xy_circle_VxVy_2[1,*]-y_part_VxVy_2[0],linestyle=2
legend,['dt=0Myr'],box=0,charsize=0.8

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_VxVy_2[k]],[scale_dJ_part_VxVy_2[k]],psym=3,color=qcolor[k]

; q histogram
plothist,q_part_VxVy_2,xr=[qmin,qmax],bin=0.1,yr=[0,290],xtitle='q',ytitle='N'

;=== 40Myr later than the closest approach ===
plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_VxVy_3[k]-x_part_VxVy_3[0]],[y_part_VxVy_3[k]-y_part_VxVy_3[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_VxVy_3[0,*]-x_part_VxVy_3[0],Rs_xy_circle_VxVy_3[1,*]-y_part_VxVy_3[0]
oplot,Rtidal_xy_circle_VxVy_3[0,*]-x_part_VxVy_3[0],Rtidal_xy_circle_VxVy_3[1,*]-y_part_VxVy_3[0],linestyle=2
legend,['dt=40Myr'],box=0,charsize=0.8

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_VxVy_3[k]],[scale_dJ_part_VxVy_3[k]],psym=3,color=qcolor[k]

; q histogram
plothist,q_part_VxVy_3,xr=[qmin,qmax],bin=0.1,yr=[0,290],xtitle='q',ytitle='N'

device,/close

END
