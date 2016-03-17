pro bcrit_test

dir_out = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_nosub = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_nosub/'

dir_part_e6_1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e6_1Rs/'
dir_subhalo_e6_1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e6_1Rs/'
dir_part_e6_2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e6_2Rs/'
dir_subhalo_e6_2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e6_2Rs/'
dir_part_e6_4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e6_4Rs/'
dir_subhalo_e6_4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e6_4Rs/'

dir_part_e7_1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e7_1Rs/'
dir_subhalo_e7_1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e7_1Rs/'
dir_part_e7_2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e7_2Rs/'
dir_subhalo_e7_2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e7_2Rs/'
dir_part_e7_4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e7_4Rs/'
dir_subhalo_e7_4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e7_4Rs/'

dir_part_e8_1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e8_1Rs/'
dir_subhalo_e8_1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e8_1Rs/'
dir_part_e8_2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e8_2Rs/'
dir_subhalo_e8_2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e8_2Rs/'
dir_part_e8_4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_close_e8_4Rs/'
dir_subhalo_e8_4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_close_e8_4Rs/'


G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
r_tide = 0.108653
Msat = 20000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

;==============
;=== v test ===
;==============

;=== e6_1Rs pericenter ===
;chr_rdtbl,dir_out+'part_peri_xy_e6_1Rs',0, arr
chr_rdtbl,dir_out+'part_peri_nosub',0, arr
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
epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)


;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part_e6_1Rs+'snap00000',0, arr
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (r_tide/r_peri * J_initial)
qmin = min(q_initial)
qmax = max(q_initial)
dJmin = min(scale_dJ_initial)
dJmax = max(scale_dJ_initial)
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)
qmin = qmin-0.5
qmax = qmax+0.5
dJmin = dJmin-0.1
dJmax = dJmax+0.1


;=== read out simulation info ===
chr_rdtbl,'orbit_subhalo_e7.dat',0,arr
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)
N_time = N_elements(fname_part_e6_1Rs)

Rs_e6 = 0.11058515
Rtidal_e6 = 0.87319397
Rs_e8 = 0.78954220
Rtidal_e8 = 16.358786


;=== e6_1Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e6_1Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e6_1Rs = reform(arr[0,*])
x_part_e6_1Rs = reform(arr[1,*])
y_part_e6_1Rs = reform(arr[2,*])
z_part_e6_1Rs = reform(arr[3,*])
Vx_part_e6_1Rs = reform(arr[4,*])
Vy_part_e6_1Rs = reform(arr[5,*])
Vz_part_e6_1Rs = reform(arr[6,*])
dE_part_e6_1Rs = reform(arr[7,*])
Etot_part_e6_1Rs = reform(arr[8,*])
J_part_e6_1Rs = reform(arr[9,*])
delE = Etot_part_e6_1Rs-Etot_part_e6_1Rs[0]
q_part_e6_1Rs = delE/epsilon
dJ_part_e6_1Rs = J_part_e6_1Rs-J_part_e6_1Rs[0]
scale_dJ_part_e6_1Rs = dJ_part_e6_1Rs / (r_tide/r_peri * J_part_e6_1Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6_1Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e6_1Rs = reform(arr[1,*])
y_subhalo_e6_1Rs = reform(arr[2,*])
z_subhalo_e6_1Rs = reform(arr[3,*])
Vx_subhalo_e6_1Rs = reform(arr[4,*])
Vy_subhalo_e6_1Rs = reform(arr[5,*])
Vz_subhalo_e6_1Rs = reform(arr[6,*])

t_step = string(t_part_e6_1Rs[0]/1000.,f='(f7.4)')

; make a circle of a subhalo 
Rs_xy_circle_e6_1Rs = circle(x_subhalo_e6_1Rs[0],y_subhalo_e6_1Rs[0],double(Rs[0]))
Rtidal_xy_circle_e6_1Rs = circle(x_subhalo_e6_1Rs[0],y_subhalo_e6_1Rs[0],double(Rtidal[0]))
Rs_yz_circle_e6_1Rs = circle(y_subhalo_e6_1Rs[0],z_subhalo_e6_1Rs[0],double(Rs[0]))
Rtidal_yz_circle_e6_1Rs = circle(y_subhalo_e6_1Rs[0],z_subhalo_e6_1Rs[0],double(Rtidal[0]))
Rs_xz_circle_e6_1Rs = circle(x_subhalo_e6_1Rs[0],z_subhalo_e6_1Rs[0],double(Rs[0]))
Rtidal_xz_circle_e6_1Rs = circle(x_subhalo_e6_1Rs[0],z_subhalo_e6_1Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e6_1Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e6_1Rs_2 = reform(arr[0,*])
x_part_e6_1Rs_2 = reform(arr[1,*])
y_part_e6_1Rs_2 = reform(arr[2,*])
z_part_e6_1Rs_2 = reform(arr[3,*])
Vx_part_e6_1Rs_2 = reform(arr[4,*])
Vy_part_e6_1Rs_2 = reform(arr[5,*])
Vz_part_e6_1Rs_2 = reform(arr[6,*])
dE_part_e6_1Rs_2 = reform(arr[7,*])
Etot_part_e6_1Rs_2 = reform(arr[8,*])
J_part_e6_1Rs_2 = reform(arr[9,*])
delE = Etot_part_e6_1Rs_2-Etot_part_e6_1Rs_2[0]
q_part_e6_1Rs_2 = delE/epsilon
dJ_part_e6_1Rs_2 = J_part_e6_1Rs_2-J_part_e6_1Rs_2[0]
scale_dJ_part_e6_1Rs_2 = dJ_part_e6_1Rs_2 / (r_tide/r_peri * J_part_e6_1Rs_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6_1Rs+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_e6_1Rs_2 = reform(arr[1,*])
y_subhalo_e6_1Rs_2 = reform(arr[2,*])
z_subhalo_e6_1Rs_2 = reform(arr[3,*])
Vx_subhalo_e6_1Rs_2 = reform(arr[4,*])
Vy_subhalo_e6_1Rs_2 = reform(arr[5,*])
Vz_subhalo_e6_1Rs_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_1Rs_2 = circle(x_subhalo_e6_1Rs_2[0],y_subhalo_e6_1Rs_2[0],double(Rs[0]))
Rtidal_xy_circle_e6_1Rs_2 = circle(x_subhalo_e6_1Rs_2[0],y_subhalo_e6_1Rs_2[0],double(Rtidal[0]))
Rs_yz_circle_e6_1Rs_2 = circle(y_subhalo_e6_1Rs_2[0],z_subhalo_e6_1Rs_2[0],double(Rs[0]))
Rtidal_yz_circle_e6_1Rs_2 = circle(y_subhalo_e6_1Rs_2[0],z_subhalo_e6_1Rs_2[0],double(Rtidal[0]))
Rs_xz_circle_e6_1Rs_2 = circle(x_subhalo_e6_1Rs_2[0],z_subhalo_e6_1Rs_2[0],double(Rs[0]))
Rtidal_xz_circle_e6_1Rs_2 = circle(x_subhalo_e6_1Rs_2[0],z_subhalo_e6_1Rs_2[0],double(Rtidal[0]))

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_e6_1Rs+'subhalo00980',0,arr,/silent
arr = double(arr)
x_subhalo_e6_1Rs_3 = reform(arr[1,*])
y_subhalo_e6_1Rs_3 = reform(arr[2,*])
z_subhalo_e6_1Rs_3 = reform(arr[3,*])
Vx_subhalo_e6_1Rs_3 = reform(arr[4,*])
Vy_subhalo_e6_1Rs_3 = reform(arr[5,*])
Vz_subhalo_e6_1Rs_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_1Rs_3 = circle(x_subhalo_e6_1Rs_3[0],y_subhalo_e6_1Rs_3[0],double(Rs[0]))
Rtidal_xy_circle_e6_1Rs_3 = circle(x_subhalo_e6_1Rs_3[0],y_subhalo_e6_1Rs_3[0],double(Rtidal[0]))
Rs_yz_circle_e6_1Rs_3 = circle(y_subhalo_e6_1Rs_3[0],z_subhalo_e6_1Rs_3[0],double(Rs[0]))
Rtidal_yz_circle_e6_1Rs_3 = circle(y_subhalo_e6_1Rs_3[0],z_subhalo_e6_1Rs_3[0],double(Rtidal[0]))
Rs_xz_circle_e6_1Rs_3 = circle(x_subhalo_e6_1Rs_3[0],z_subhalo_e6_1Rs_3[0],double(Rs[0]))
Rtidal_xz_circle_e6_1Rs_3 = circle(x_subhalo_e6_1Rs_3[0],z_subhalo_e6_1Rs_3[0],double(Rtidal[0]))


;=== e6_2Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e6_2Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e6_2Rs = reform(arr[0,*])
x_part_e6_2Rs = reform(arr[1,*])
y_part_e6_2Rs = reform(arr[2,*])
z_part_e6_2Rs = reform(arr[3,*])
Vx_part_e6_2Rs = reform(arr[4,*])
Vy_part_e6_2Rs = reform(arr[5,*])
Vz_part_e6_2Rs = reform(arr[6,*])
dE_part_e6_2Rs = reform(arr[7,*])
Etot_part_e6_2Rs = reform(arr[8,*])
J_part_e6_2Rs = reform(arr[9,*])
delE = Etot_part_e6_2Rs-Etot_part_e6_2Rs[0]
q_part_e6_2Rs = delE/epsilon
dJ_part_e6_2Rs = J_part_e6_2Rs-J_part_e6_2Rs[0]
scale_dJ_part_e6_2Rs = dJ_part_e6_2Rs / (r_tide/r_peri * J_part_e6_2Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6_2Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e6_2Rs = reform(arr[1,*])
y_subhalo_e6_2Rs = reform(arr[2,*])
z_subhalo_e6_2Rs = reform(arr[3,*])
Vx_subhalo_e6_2Rs = reform(arr[4,*])
Vy_subhalo_e6_2Rs = reform(arr[5,*])
Vz_subhalo_e6_2Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_2Rs = circle(x_subhalo_e6_2Rs[0],y_subhalo_e6_2Rs[0],double(Rs[0]))
Rtidal_xy_circle_e6_2Rs = circle(x_subhalo_e6_2Rs[0],y_subhalo_e6_2Rs[0],double(Rtidal[0]))
Rs_yz_circle_e6_2Rs = circle(y_subhalo_e6_2Rs[0],z_subhalo_e6_2Rs[0],double(Rs[0]))
Rtidal_yz_circle_e6_2Rs = circle(y_subhalo_e6_2Rs[0],z_subhalo_e6_2Rs[0],double(Rtidal[0]))
Rs_xz_circle_e6_2Rs = circle(x_subhalo_e6_2Rs[0],z_subhalo_e6_2Rs[0],double(Rs[0]))
Rtidal_xz_circle_e6_2Rs = circle(x_subhalo_e6_2Rs[0],z_subhalo_e6_2Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e6_2Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e6_2Rs_2 = reform(arr[0,*])
x_part_e6_2Rs_2 = reform(arr[1,*])
y_part_e6_2Rs_2 = reform(arr[2,*])
z_part_e6_2Rs_2 = reform(arr[3,*])
Vx_part_e6_2Rs_2 = reform(arr[4,*])
Vy_part_e6_2Rs_2 = reform(arr[5,*])
Vz_part_e6_2Rs_2 = reform(arr[6,*])
dE_part_e6_2Rs_2 = reform(arr[7,*])
Etot_part_e6_2Rs_2 = reform(arr[8,*])
J_part_e6_2Rs_2 = reform(arr[9,*])
delE = Etot_part_e6_2Rs_2-Etot_part_e6_2Rs_2[0]
q_part_e6_2Rs_2 = delE/epsilon
dJ_part_e6_2Rs_2 = J_part_e6_2Rs_2-J_part_e6_2Rs_2[0]
scale_dJ_part_e6_2Rs_2 = dJ_part_e6_2Rs_2 / (r_tide/r_peri * J_part_e6_2Rs_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6_2Rs+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_e6_2Rs_2 = reform(arr[1,*])
y_subhalo_e6_2Rs_2 = reform(arr[2,*])
z_subhalo_e6_2Rs_2 = reform(arr[3,*])
Vx_subhalo_e6_2Rs_2 = reform(arr[4,*])
Vy_subhalo_e6_2Rs_2 = reform(arr[5,*])
Vz_subhalo_e6_2Rs_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_2Rs_2 = circle(x_subhalo_e6_2Rs_2[0],y_subhalo_e6_2Rs_2[0],double(Rs[0]))
Rtidal_xy_circle_e6_2Rs_2 = circle(x_subhalo_e6_2Rs_2[0],y_subhalo_e6_2Rs_2[0],double(Rtidal[0]))
Rs_yz_circle_e6_2Rs_2 = circle(y_subhalo_e6_2Rs_2[0],z_subhalo_e6_2Rs_2[0],double(Rs[0]))
Rtidal_yz_circle_e6_2Rs_2 = circle(y_subhalo_e6_2Rs_2[0],z_subhalo_e6_2Rs_2[0],double(Rtidal[0]))
Rs_xz_circle_e6_2Rs_2 = circle(x_subhalo_e6_2Rs_2[0],z_subhalo_e6_2Rs_2[0],double(Rs[0]))
Rtidal_xz_circle_e6_2Rs_2 = circle(x_subhalo_e6_2Rs_2[0],z_subhalo_e6_2Rs_2[0],double(Rtidal[0]))

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_e6_2Rs+'subhalo00980',0,arr,/silent
arr = double(arr)
x_subhalo_e6_2Rs_3 = reform(arr[1,*])
y_subhalo_e6_2Rs_3 = reform(arr[2,*])
z_subhalo_e6_2Rs_3 = reform(arr[3,*])
Vx_subhalo_e6_2Rs_3 = reform(arr[4,*])
Vy_subhalo_e6_2Rs_3 = reform(arr[5,*])
Vz_subhalo_e6_2Rs_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_2Rs_3 = circle(x_subhalo_e6_2Rs_3[0],y_subhalo_e6_2Rs_3[0],double(Rs[0]))
Rtidal_xy_circle_e6_2Rs_3 = circle(x_subhalo_e6_2Rs_3[0],y_subhalo_e6_2Rs_3[0],double(Rtidal[0]))
Rs_yz_circle_e6_2Rs_3 = circle(y_subhalo_e6_2Rs_3[0],z_subhalo_e6_2Rs_3[0],double(Rs[0]))
Rtidal_yz_circle_e6_2Rs_3 = circle(y_subhalo_e6_2Rs_3[0],z_subhalo_e6_2Rs_3[0],double(Rtidal[0]))
Rs_xz_circle_e6_2Rs_3 = circle(x_subhalo_e6_2Rs_3[0],z_subhalo_e6_2Rs_3[0],double(Rs[0]))
Rtidal_xz_circle_e6_2Rs_3 = circle(x_subhalo_e6_2Rs_3[0],z_subhalo_e6_2Rs_3[0],double(Rtidal[0]))


;=== e6_4Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e6_4Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e6_4Rs = reform(arr[0,*])
x_part_e6_4Rs = reform(arr[1,*])
y_part_e6_4Rs = reform(arr[2,*])
z_part_e6_4Rs = reform(arr[3,*])
Vx_part_e6_4Rs = reform(arr[4,*])
Vy_part_e6_4Rs = reform(arr[5,*])
Vz_part_e6_4Rs = reform(arr[6,*])
dE_part_e6_4Rs = reform(arr[7,*])
Etot_part_e6_4Rs = reform(arr[8,*])
J_part_e6_4Rs = reform(arr[9,*])
delE = Etot_part_e6_4Rs-Etot_part_e6_4Rs[0]
q_part_e6_4Rs = delE/epsilon
dJ_part_e6_4Rs = J_part_e6_4Rs-J_part_e6_4Rs[0]
scale_dJ_part_e6_4Rs = dJ_part_e6_4Rs / (r_tide/r_peri * J_part_e6_4Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6_4Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e6_4Rs = reform(arr[1,*])
y_subhalo_e6_4Rs = reform(arr[2,*])
z_subhalo_e6_4Rs = reform(arr[3,*])
Vx_subhalo_e6_4Rs = reform(arr[4,*])
Vy_subhalo_e6_4Rs = reform(arr[5,*])
Vz_subhalo_e6_4Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_4Rs = circle(x_subhalo_e6_4Rs[0],y_subhalo_e6_4Rs[0],double(Rs[0]))
Rtidal_xy_circle_e6_4Rs = circle(x_subhalo_e6_4Rs[0],y_subhalo_e6_4Rs[0],double(Rtidal[0]))
Rs_yz_circle_e6_4Rs = circle(y_subhalo_e6_4Rs[0],z_subhalo_e6_4Rs[0],double(Rs[0]))
Rtidal_yz_circle_e6_4Rs = circle(y_subhalo_e6_4Rs[0],z_subhalo_e6_4Rs[0],double(Rtidal[0]))
Rs_xz_circle_e6_4Rs = circle(x_subhalo_e6_4Rs[0],z_subhalo_e6_4Rs[0],double(Rs[0]))
Rtidal_xz_circle_e6_4Rs = circle(x_subhalo_e6_4Rs[0],z_subhalo_e6_4Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e6_4Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e6_4Rs_2 = reform(arr[0,*])
x_part_e6_4Rs_2 = reform(arr[1,*])
y_part_e6_4Rs_2 = reform(arr[2,*])
z_part_e6_4Rs_2 = reform(arr[3,*])
Vx_part_e6_4Rs_2 = reform(arr[4,*])
Vy_part_e6_4Rs_2 = reform(arr[5,*])
Vz_part_e6_4Rs_2 = reform(arr[6,*])
dE_part_e6_4Rs_2 = reform(arr[7,*])
Etot_part_e6_4Rs_2 = reform(arr[8,*])
J_part_e6_4Rs_2 = reform(arr[9,*])
delE = Etot_part_e6_4Rs_2-Etot_part_e6_4Rs_2[0]
q_part_e6_4Rs_2 = delE/epsilon
dJ_part_e6_4Rs_2 = J_part_e6_4Rs_2-J_part_e6_4Rs_2[0]
scale_dJ_part_e6_4Rs_2 = dJ_part_e6_4Rs_2 / (r_tide/r_peri * J_part_e6_4Rs_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6_4Rs+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_e6_4Rs_2 = reform(arr[1,*])
y_subhalo_e6_4Rs_2 = reform(arr[2,*])
z_subhalo_e6_4Rs_2 = reform(arr[3,*])
Vx_subhalo_e6_4Rs_2 = reform(arr[4,*])
Vy_subhalo_e6_4Rs_2 = reform(arr[5,*])
Vz_subhalo_e6_4Rs_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_4Rs_2 = circle(x_subhalo_e6_4Rs_2[0],y_subhalo_e6_4Rs_2[0],double(Rs[0]))
Rtidal_xy_circle_e6_4Rs_2 = circle(x_subhalo_e6_4Rs_2[0],y_subhalo_e6_4Rs_2[0],double(Rtidal[0]))
Rs_yz_circle_e6_4Rs_2 = circle(y_subhalo_e6_4Rs_2[0],z_subhalo_e6_4Rs_2[0],double(Rs[0]))
Rtidal_yz_circle_e6_4Rs_2 = circle(y_subhalo_e6_4Rs_2[0],z_subhalo_e6_4Rs_2[0],double(Rtidal[0]))
Rs_xz_circle_e6_4Rs_2 = circle(x_subhalo_e6_4Rs_2[0],z_subhalo_e6_4Rs_2[0],double(Rs[0]))
Rtidal_xz_circle_e6_4Rs_2 = circle(x_subhalo_e6_4Rs_2[0],z_subhalo_e6_4Rs_2[0],double(Rtidal[0]))

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_e6_4Rs+'subhalo00980',0,arr,/silent
arr = double(arr)
x_subhalo_e6_4Rs_3 = reform(arr[1,*])
y_subhalo_e6_4Rs_3 = reform(arr[2,*])
z_subhalo_e6_4Rs_3 = reform(arr[3,*])
Vx_subhalo_e6_4Rs_3 = reform(arr[4,*])
Vy_subhalo_e6_4Rs_3 = reform(arr[5,*])
Vz_subhalo_e6_4Rs_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6_4Rs_3 = circle(x_subhalo_e6_4Rs_3[0],y_subhalo_e6_4Rs_3[0],double(Rs[0]))
Rtidal_xy_circle_e6_4Rs_3 = circle(x_subhalo_e6_4Rs_3[0],y_subhalo_e6_4Rs_3[0],double(Rtidal[0]))
Rs_yz_circle_e6_4Rs_3 = circle(y_subhalo_e6_4Rs_3[0],z_subhalo_e6_4Rs_3[0],double(Rs[0]))
Rtidal_yz_circle_e6_4Rs_3 = circle(y_subhalo_e6_4Rs_3[0],z_subhalo_e6_4Rs_3[0],double(Rtidal[0]))
Rs_xz_circle_e6_4Rs_3 = circle(x_subhalo_e6_4Rs_3[0],z_subhalo_e6_4Rs_3[0],double(Rs[0]))
Rtidal_xz_circle_e6_4Rs_3 = circle(x_subhalo_e6_4Rs_3[0],z_subhalo_e6_4Rs_3[0],double(Rtidal[0]))


;=====================
;=== b_impact test ===
;=====================
;=== b=0 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e7_1Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e7_1Rs = reform(arr[0,*])
x_part_e7_1Rs = reform(arr[1,*])
y_part_e7_1Rs = reform(arr[2,*])
z_part_e7_1Rs = reform(arr[3,*])
Vx_part_e7_1Rs = reform(arr[4,*])
Vy_part_e7_1Rs = reform(arr[5,*])
Vz_part_e7_1Rs = reform(arr[6,*])
dE_part_e7_1Rs = reform(arr[7,*])
Etot_part_e7_1Rs = reform(arr[8,*])
J_part_e7_1Rs = reform(arr[9,*])
delE = Etot_part_e7_1Rs-Etot_part_e7_1Rs[0]
q_part_e7_1Rs = delE/epsilon
dJ_part_e7_1Rs = J_part_e7_1Rs-J_part_e7_1Rs[0]
scale_dJ_part_e7_1Rs = dJ_part_e7_1Rs / (r_tide/r_peri * J_part_e7_1Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e7_1Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e7_1Rs = reform(arr[1,*])
y_subhalo_e7_1Rs = reform(arr[2,*])
z_subhalo_e7_1Rs = reform(arr[3,*])
Vx_subhalo_e7_1Rs = reform(arr[4,*])
Vy_subhalo_e7_1Rs = reform(arr[5,*])
Vz_subhalo_e7_1Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e7_1Rs = circle(x_subhalo_e7_1Rs[0],y_subhalo_e7_1Rs[0],double(Rs[0]))
Rtidal_xy_circle_e7_1Rs = circle(x_subhalo_e7_1Rs[0],y_subhalo_e7_1Rs[0],double(Rtidal[0]))
Rs_yz_circle_e7_1Rs = circle(y_subhalo_e7_1Rs[0],z_subhalo_e7_1Rs[0],double(Rs[0]))
Rtidal_yz_circle_e7_1Rs = circle(y_subhalo_e7_1Rs[0],z_subhalo_e7_1Rs[0],double(Rtidal[0]))
Rs_xz_circle_e7_1Rs = circle(x_subhalo_e7_1Rs[0],z_subhalo_e7_1Rs[0],double(Rs[0]))
Rtidal_xz_circle_e7_1Rs = circle(x_subhalo_e7_1Rs[0],z_subhalo_e7_1Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e7_1Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e7_1Rs_2 = reform(arr[0,*])
x_part_e7_1Rs_2 = reform(arr[1,*])
y_part_e7_1Rs_2 = reform(arr[2,*])
z_part_e7_1Rs_2 = reform(arr[3,*])
Vx_part_e7_1Rs_2 = reform(arr[4,*])
Vy_part_e7_1Rs_2 = reform(arr[5,*])
Vz_part_e7_1Rs_2 = reform(arr[6,*])
dE_part_e7_1Rs_2 = reform(arr[7,*])
Etot_part_e7_1Rs_2 = reform(arr[8,*])
J_part_e7_1Rs_2 = reform(arr[9,*])
delE = Etot_part_e7_1Rs_2-Etot_part_e7_1Rs_2[0]
q_part_e7_1Rs_2 = delE/epsilon
dJ_part_e7_1Rs_2 = J_part_e7_1Rs_2-J_part_e7_1Rs_2[0]
scale_dJ_part_e7_1Rs_2 = dJ_part_e7_1Rs_2 / (r_tide/r_peri * J_part_e7_1Rs_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e7_1Rs+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_e7_1Rs_2 = reform(arr[1,*])
y_subhalo_e7_1Rs_2 = reform(arr[2,*])
z_subhalo_e7_1Rs_2 = reform(arr[3,*])
Vx_subhalo_e7_1Rs_2 = reform(arr[4,*])
Vy_subhalo_e7_1Rs_2 = reform(arr[5,*])
Vz_subhalo_e7_1Rs_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e7_1Rs_2 = circle(x_subhalo_e7_1Rs_2[0],y_subhalo_e7_1Rs_2[0],double(Rs[0]))
Rtidal_xy_circle_e7_1Rs_2 = circle(x_subhalo_e7_1Rs_2[0],y_subhalo_e7_1Rs_2[0],double(Rtidal[0]))
Rs_yz_circle_e7_1Rs_2 = circle(y_subhalo_e7_1Rs_2[0],z_subhalo_e7_1Rs_2[0],double(Rs[0]))
Rtidal_yz_circle_e7_1Rs_2 = circle(y_subhalo_e7_1Rs_2[0],z_subhalo_e7_1Rs_2[0],double(Rtidal[0]))
Rs_xz_circle_e7_1Rs_2 = circle(x_subhalo_e7_1Rs_2[0],z_subhalo_e7_1Rs_2[0],double(Rs[0]))
Rtidal_xz_circle_e7_1Rs_2 = circle(x_subhalo_e7_1Rs_2[0],z_subhalo_e7_1Rs_2[0],double(Rtidal[0]))


;=== b=1Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e7_2Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e7_2Rs = reform(arr[0,*])
x_part_e7_2Rs = reform(arr[1,*])
y_part_e7_2Rs = reform(arr[2,*])
z_part_e7_2Rs = reform(arr[3,*])
Vx_part_e7_2Rs = reform(arr[4,*])
Vy_part_e7_2Rs = reform(arr[5,*])
Vz_part_e7_2Rs = reform(arr[6,*])
dE_part_e7_2Rs = reform(arr[7,*])
Etot_part_e7_2Rs = reform(arr[8,*])
J_part_e7_2Rs = reform(arr[9,*])
delE = Etot_part_e7_2Rs-Etot_part_e7_2Rs[0]
q_part_e7_2Rs = delE/epsilon
dJ_part_e7_2Rs = J_part_e7_2Rs-J_part_e7_2Rs[0]
scale_dJ_part_e7_2Rs = dJ_part_e7_2Rs / (r_tide/r_peri * J_part_e7_2Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e7_2Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e7_2Rs = reform(arr[1,*])
y_subhalo_e7_2Rs = reform(arr[2,*])
z_subhalo_e7_2Rs = reform(arr[3,*])
Vx_subhalo_e7_2Rs = reform(arr[4,*])
Vy_subhalo_e7_2Rs = reform(arr[5,*])
Vz_subhalo_e7_2Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e7_2Rs = circle(x_subhalo_e7_2Rs[0],y_subhalo_e7_2Rs[0],double(Rs[0]))
Rtidal_xy_circle_e7_2Rs = circle(x_subhalo_e7_2Rs[0],y_subhalo_e7_2Rs[0],double(Rtidal[0]))
Rs_yz_circle_e7_2Rs = circle(y_subhalo_e7_2Rs[0],z_subhalo_e7_2Rs[0],double(Rs[0]))
Rtidal_yz_circle_e7_2Rs = circle(y_subhalo_e7_2Rs[0],z_subhalo_e7_2Rs[0],double(Rtidal[0]))
Rs_xz_circle_e7_2Rs = circle(x_subhalo_e7_2Rs[0],z_subhalo_e7_2Rs[0],double(Rs[0]))
Rtidal_xz_circle_e7_2Rs = circle(x_subhalo_e7_2Rs[0],z_subhalo_e7_2Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e7_2Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e7_2Rs_2 = reform(arr[0,*])
x_part_e7_2Rs_2 = reform(arr[1,*])
y_part_e7_2Rs_2 = reform(arr[2,*])
z_part_e7_2Rs_2 = reform(arr[3,*])
Vx_part_e7_2Rs_2 = reform(arr[4,*])
Vy_part_e7_2Rs_2 = reform(arr[5,*])
Vz_part_e7_2Rs_2 = reform(arr[6,*])
dE_part_e7_2Rs_2 = reform(arr[7,*])
Etot_part_e7_2Rs_2 = reform(arr[8,*])
J_part_e7_2Rs_2 = reform(arr[9,*])
delE = Etot_part_e7_2Rs_2-Etot_part_e7_2Rs_2[0]
q_part_e7_2Rs_2 = delE/epsilon
dJ_part_e7_2Rs_2 = J_part_e7_2Rs_2-J_part_e7_2Rs_2[0]
scale_dJ_part_e7_2Rs_2 = dJ_part_e7_2Rs_2 / (r_tide/r_peri * J_part_e7_2Rs_2)


;=== b=2Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e7_4Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e7_4Rs = reform(arr[0,*])
x_part_e7_4Rs = reform(arr[1,*])
y_part_e7_4Rs = reform(arr[2,*])
z_part_e7_4Rs = reform(arr[3,*])
Vx_part_e7_4Rs = reform(arr[4,*])
Vy_part_e7_4Rs = reform(arr[5,*])
Vz_part_e7_4Rs = reform(arr[6,*])
dE_part_e7_4Rs = reform(arr[7,*])
Etot_part_e7_4Rs = reform(arr[8,*])
J_part_e7_4Rs = reform(arr[9,*])
delE = Etot_part_e7_4Rs-Etot_part_e7_4Rs[0]
q_part_e7_4Rs = delE/epsilon
dJ_part_e7_4Rs = J_part_e7_4Rs-J_part_e7_4Rs[0]
scale_dJ_part_e7_4Rs = dJ_part_e7_4Rs / (r_tide/r_peri * J_part_e7_4Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e7_4Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e7_4Rs = reform(arr[1,*])
y_subhalo_e7_4Rs = reform(arr[2,*])
z_subhalo_e7_4Rs = reform(arr[3,*])
Vx_subhalo_e7_4Rs = reform(arr[4,*])
Vy_subhalo_e7_4Rs = reform(arr[5,*])
Vz_subhalo_e7_4Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e7_4Rs = circle(x_subhalo_e7_4Rs[0],y_subhalo_e7_4Rs[0],double(Rs[0]))
Rtidal_xy_circle_e7_4Rs = circle(x_subhalo_e7_4Rs[0],y_subhalo_e7_4Rs[0],double(Rtidal[0]))
Rs_yz_circle_e7_4Rs = circle(y_subhalo_e7_4Rs[0],z_subhalo_e7_4Rs[0],double(Rs[0]))
Rtidal_yz_circle_e7_4Rs = circle(y_subhalo_e7_4Rs[0],z_subhalo_e7_4Rs[0],double(Rtidal[0]))
Rs_xz_circle_e7_4Rs = circle(x_subhalo_e7_4Rs[0],z_subhalo_e7_4Rs[0],double(Rs[0]))
Rtidal_xz_circle_e7_4Rs = circle(x_subhalo_e7_4Rs[0],z_subhalo_e7_4Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e7_4Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e7_4Rs_2 = reform(arr[0,*])
x_part_e7_4Rs_2 = reform(arr[1,*])
y_part_e7_4Rs_2 = reform(arr[2,*])
z_part_e7_4Rs_2 = reform(arr[3,*])
Vx_part_e7_4Rs_2 = reform(arr[4,*])
Vy_part_e7_4Rs_2 = reform(arr[5,*])
Vz_part_e7_4Rs_2 = reform(arr[6,*])
dE_part_e7_4Rs_2 = reform(arr[7,*])
Etot_part_e7_4Rs_2 = reform(arr[8,*])
J_part_e7_4Rs_2 = reform(arr[9,*])
delE = Etot_part_e7_4Rs_2-Etot_part_e7_4Rs_2[0]
q_part_e7_4Rs_2 = delE/epsilon
dJ_part_e7_4Rs_2 = J_part_e7_4Rs_2-J_part_e7_4Rs_2[0]
scale_dJ_part_e7_4Rs_2 = dJ_part_e7_4Rs_2 / (r_tide/r_peri * J_part_e7_4Rs_2)


;=====================
;=== Msubhalo test ===
;=====================
;=== Msubhalo=e8_1Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e8_1Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e8_1Rs = reform(arr[0,*])
x_part_e8_1Rs = reform(arr[1,*])
y_part_e8_1Rs = reform(arr[2,*])
z_part_e8_1Rs = reform(arr[3,*])
Vx_part_e8_1Rs = reform(arr[4,*])
Vy_part_e8_1Rs = reform(arr[5,*])
Vz_part_e8_1Rs = reform(arr[6,*])
dE_part_e8_1Rs = reform(arr[7,*])
Etot_part_e8_1Rs = reform(arr[8,*])
J_part_e8_1Rs = reform(arr[9,*])
delE = Etot_part_e8_1Rs-Etot_part_e8_1Rs[0]
q_part_e8_1Rs = delE/epsilon
dJ_part_e8_1Rs = J_part_e8_1Rs-J_part_e8_1Rs[0]
scale_dJ_part_e8_1Rs = dJ_part_e8_1Rs / (r_tide/r_peri * J_part_e8_1Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e8_1Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e8_1Rs = reform(arr[1,*])
y_subhalo_e8_1Rs = reform(arr[2,*])
z_subhalo_e8_1Rs = reform(arr[3,*])
Vx_subhalo_e8_1Rs = reform(arr[4,*])
Vy_subhalo_e8_1Rs = reform(arr[5,*])
Vz_subhalo_e8_1Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e8_1Rs = circle(x_subhalo_e8_1Rs[0],y_subhalo_e8_1Rs[0],double(Rs[0]))
Rtidal_xy_circle_e8_1Rs = circle(x_subhalo_e8_1Rs[0],y_subhalo_e8_1Rs[0],double(Rtidal[0]))
Rs_yz_circle_e8_1Rs = circle(y_subhalo_e8_1Rs[0],z_subhalo_e8_1Rs[0],double(Rs[0]))
Rtidal_yz_circle_e8_1Rs = circle(y_subhalo_e8_1Rs[0],z_subhalo_e8_1Rs[0],double(Rtidal[0]))
Rs_xz_circle_e8_1Rs = circle(x_subhalo_e8_1Rs[0],z_subhalo_e8_1Rs[0],double(Rs[0]))
Rtidal_xz_circle_e8_1Rs = circle(x_subhalo_e8_1Rs[0],z_subhalo_e8_1Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e8_1Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e8_1Rs_2 = reform(arr[0,*])
x_part_e8_1Rs_2 = reform(arr[1,*])
y_part_e8_1Rs_2 = reform(arr[2,*])
z_part_e8_1Rs_2 = reform(arr[3,*])
Vx_part_e8_1Rs_2 = reform(arr[4,*])
Vy_part_e8_1Rs_2 = reform(arr[5,*])
Vz_part_e8_1Rs_2 = reform(arr[6,*])
dE_part_e8_1Rs_2 = reform(arr[7,*])
Etot_part_e8_1Rs_2 = reform(arr[8,*])
J_part_e8_1Rs_2 = reform(arr[9,*])
delE = Etot_part_e8_1Rs_2-Etot_part_e8_1Rs_2[0]
q_part_e8_1Rs_2 = delE/epsilon
dJ_part_e8_1Rs_2 = J_part_e8_1Rs_2-J_part_e8_1Rs_2[0]
scale_dJ_part_e8_1Rs_2 = dJ_part_e8_1Rs_2 / (r_tide/r_peri * J_part_e8_1Rs_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e8_1Rs+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_e8_1Rs_2 = reform(arr[1,*])
y_subhalo_e8_1Rs_2 = reform(arr[2,*])
z_subhalo_e8_1Rs_2 = reform(arr[3,*])
Vx_subhalo_e8_1Rs_2 = reform(arr[4,*])
Vy_subhalo_e8_1Rs_2 = reform(arr[5,*])
Vz_subhalo_e8_1Rs_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e8_1Rs_2 = circle(x_subhalo_e8_1Rs_2[0],y_subhalo_e8_1Rs_2[0],double(Rs[0]))
Rtidal_xy_circle_e8_1Rs_2 = circle(x_subhalo_e8_1Rs_2[0],y_subhalo_e8_1Rs_2[0],double(Rtidal[0]))
Rs_yz_circle_e8_1Rs_2 = circle(y_subhalo_e8_1Rs_2[0],z_subhalo_e8_1Rs_2[0],double(Rs[0]))
Rtidal_yz_circle_e8_1Rs_2 = circle(y_subhalo_e8_1Rs_2[0],z_subhalo_e8_1Rs_2[0],double(Rtidal[0]))
Rs_xz_circle_e8_1Rs_2 = circle(x_subhalo_e8_1Rs_2[0],z_subhalo_e8_1Rs_2[0],double(Rs[0]))
Rtidal_xz_circle_e8_1Rs_2 = circle(x_subhalo_e8_1Rs_2[0],z_subhalo_e8_1Rs_2[0],double(Rtidal[0]))


;=== Msubhalo=e8_2Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e8_2Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e8_2Rs = reform(arr[0,*])
x_part_e8_2Rs = reform(arr[1,*])
y_part_e8_2Rs = reform(arr[2,*])
z_part_e8_2Rs = reform(arr[3,*])
Vx_part_e8_2Rs = reform(arr[4,*])
Vy_part_e8_2Rs = reform(arr[5,*])
Vz_part_e8_2Rs = reform(arr[6,*])
dE_part_e8_2Rs = reform(arr[7,*])
Etot_part_e8_2Rs = reform(arr[8,*])
J_part_e8_2Rs = reform(arr[9,*])
delE = Etot_part_e8_2Rs-Etot_part_e8_2Rs[0]
q_part_e8_2Rs = delE/epsilon
dJ_part_e8_2Rs = J_part_e8_2Rs-J_part_e8_2Rs[0]
scale_dJ_part_e8_2Rs = dJ_part_e8_2Rs / (r_tide/r_peri * J_part_e8_2Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e8_2Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e8_2Rs = reform(arr[1,*])
y_subhalo_e8_2Rs = reform(arr[2,*])
z_subhalo_e8_2Rs = reform(arr[3,*])
Vx_subhalo_e8_2Rs = reform(arr[4,*])
Vy_subhalo_e8_2Rs = reform(arr[5,*])
Vz_subhalo_e8_2Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e8_2Rs = circle(x_subhalo_e8_2Rs[0],y_subhalo_e8_2Rs[0],double(Rs_e8[0]))
Rtidal_xy_circle_e8_2Rs = circle(x_subhalo_e8_2Rs[0],y_subhalo_e8_2Rs[0],double(Rtidal_e8[0]))
Rs_yz_circle_e8_2Rs = circle(y_subhalo_e8_2Rs[0],z_subhalo_e8_2Rs[0],double(Rs_e8[0]))
Rtidal_yz_circle_e8_2Rs = circle(y_subhalo_e8_2Rs[0],z_subhalo_e8_2Rs[0],double(Rtidal_e8[0]))
Rs_xz_circle_e8_2Rs = circle(x_subhalo_e8_2Rs[0],z_subhalo_e8_2Rs[0],double(Rs_e8[0]))
Rtidal_xz_circle_e8_2Rs = circle(x_subhalo_e8_2Rs[0],z_subhalo_e8_2Rs[0],double(Rtidal_e8[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e8_2Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e8_2Rs_2 = reform(arr[0,*])
x_part_e8_2Rs_2 = reform(arr[1,*])
y_part_e8_2Rs_2 = reform(arr[2,*])
z_part_e8_2Rs_2 = reform(arr[3,*])
Vx_part_e8_2Rs_2 = reform(arr[4,*])
Vy_part_e8_2Rs_2 = reform(arr[5,*])
Vz_part_e8_2Rs_2 = reform(arr[6,*])
dE_part_e8_2Rs_2 = reform(arr[7,*])
Etot_part_e8_2Rs_2 = reform(arr[8,*])
J_part_e8_2Rs_2 = reform(arr[9,*])
delE = Etot_part_e8_2Rs_2-Etot_part_e8_2Rs_2[0]
q_part_e8_2Rs_2 = delE/epsilon
dJ_part_e8_2Rs_2 = J_part_e8_2Rs_2-J_part_e8_2Rs_2[0]
scale_dJ_part_e8_2Rs_2 = dJ_part_e8_2Rs_2 / (r_tide/r_peri * J_part_e8_2Rs_2)


;=== Msubhalo=e8_4Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e8_4Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e8_4Rs = reform(arr[0,*])
x_part_e8_4Rs = reform(arr[1,*])
y_part_e8_4Rs = reform(arr[2,*])
z_part_e8_4Rs = reform(arr[3,*])
Vx_part_e8_4Rs = reform(arr[4,*])
Vy_part_e8_4Rs = reform(arr[5,*])
Vz_part_e8_4Rs = reform(arr[6,*])
dE_part_e8_4Rs = reform(arr[7,*])
Etot_part_e8_4Rs = reform(arr[8,*])
J_part_e8_4Rs = reform(arr[9,*])
delE = Etot_part_e8_4Rs-Etot_part_e8_4Rs[0]
q_part_e8_4Rs = delE/epsilon
dJ_part_e8_4Rs = J_part_e8_4Rs-J_part_e8_4Rs[0]
scale_dJ_part_e8_4Rs = dJ_part_e8_4Rs / (r_tide/r_peri * J_part_e8_4Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e8_4Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e8_4Rs = reform(arr[1,*])
y_subhalo_e8_4Rs = reform(arr[2,*])
z_subhalo_e8_4Rs = reform(arr[3,*])
Vx_subhalo_e8_4Rs = reform(arr[4,*])
Vy_subhalo_e8_4Rs = reform(arr[5,*])
Vz_subhalo_e8_4Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e8_4Rs = circle(x_subhalo_e8_4Rs[0],y_subhalo_e8_4Rs[0],double(Rs_e8[0]))
Rtidal_xy_circle_e8_4Rs = circle(x_subhalo_e8_4Rs[0],y_subhalo_e8_4Rs[0],double(Rtidal_e8[0]))
Rs_yz_circle_e8_4Rs = circle(y_subhalo_e8_4Rs[0],z_subhalo_e8_4Rs[0],double(Rs_e8[0]))
Rtidal_yz_circle_e8_4Rs = circle(y_subhalo_e8_4Rs[0],z_subhalo_e8_4Rs[0],double(Rtidal_e8[0]))
Rs_xz_circle_e8_4Rs = circle(x_subhalo_e8_4Rs[0],z_subhalo_e8_4Rs[0],double(Rs_e8[0]))
Rtidal_xz_circle_e8_4Rs = circle(x_subhalo_e8_4Rs[0],z_subhalo_e8_4Rs[0],double(Rtidal_e8[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e8_4Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e8_4Rs_2 = reform(arr[0,*])
x_part_e8_4Rs_2 = reform(arr[1,*])
y_part_e8_4Rs_2 = reform(arr[2,*])
z_part_e8_4Rs_2 = reform(arr[3,*])
Vx_part_e8_4Rs_2 = reform(arr[4,*])
Vy_part_e8_4Rs_2 = reform(arr[5,*])
Vz_part_e8_4Rs_2 = reform(arr[6,*])
dE_part_e8_4Rs_2 = reform(arr[7,*])
Etot_part_e8_4Rs_2 = reform(arr[8,*])
J_part_e8_4Rs_2 = reform(arr[9,*])
delE = Etot_part_e8_4Rs_2-Etot_part_e8_4Rs_2[0]
q_part_e8_4Rs_2 = delE/epsilon
dJ_part_e8_4Rs_2 = J_part_e8_4Rs_2-J_part_e8_4Rs_2[0]
scale_dJ_part_e8_4Rs_2 = dJ_part_e8_4Rs_2 / (r_tide/r_peri * J_part_e8_4Rs_2)


;============
;=== plot ===
;============
set_plot,'ps'
@plot_setting
;!p.charsize=1
!p.multi=[0,4,3]
file_out = dir_out+'bcrit_test.ps'
device,file=file_out,/color,/landscape;,xoffset=0.5,yoffset=0.5


;=== plot v test ===
; at the moment of closest approach, 10Gyr
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e6_1Rs[k]-x_part_e6_1Rs[0]],[y_part_e6_1Rs[k]-y_part_e6_1Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_e6_1Rs[0,*]-x_part_e6_1Rs[0],Rs_xy_circle_e6_1Rs[1,*]-y_part_e6_1Rs[0]
oplot,Rtidal_xy_circle_e6_1Rs[0,*]-x_part_e6_1Rs[0],Rtidal_xy_circle_e6_1Rs[1,*]-y_part_e6_1Rs[0],linestyle=2

plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e6_1Rs[k]-x_part_e6_1Rs[0]],[z_part_e6_1Rs[k]-z_part_e6_1Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_e6_1Rs[0,*]-x_part_e6_1Rs[0],Rs_xz_circle_e6_1Rs[1,*]-z_part_e6_1Rs[0]
oplot,Rtidal_xz_circle_e6_1Rs[0,*]-x_part_e6_1Rs[0],Rtidal_xz_circle_e6_1Rs[1,*]-z_part_e6_1Rs[0],linestyle=2
arrow,-8.3,-3,-8.3,-8,/data

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin-0.6,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_e6_1Rs_2[k]],[scale_dJ_part_e6_1Rs_2[k]],psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e6_2Rs_2[k]],[scale_dJ_part_e6_2Rs_2[k]]-0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e6_4Rs_2[k]],[scale_dJ_part_e6_4Rs_2[k]]-0.6,psym=3,color=qcolor[k]
xyouts,0,0,textoidl('1Rs'),alignment=0.5,charsize=0.6
xyouts,0,0-0.3,textoidl('1Rs'),alignment=0.5,charsize=0.6
xyouts,0,0-0.6,textoidl('1Rs'),alignment=0.5,charsize=0.6

; q histogram
plothist,q_part_e6_1Rs_2,xr=[qmin,qmax],bin=0.1,xtitle='q',yr=[0,510]
plothist,q_part_e6_2Rs_2,xr=[qmin,qmax],bin=0.1,color=70,/overplot;,/peak;yr=[0,210]
plothist,q_part_e6_4Rs_2,xr=[qmin,qmax],bin=0.1,color=255,/overplot;,/peak;yr=[0,210]
legend,['e6_1Rs=0.2065','e6_2Rs=0.1003','e6_4Rs=0.0494'],box=0,color=[0,70,255],psym=8,/right,charsize=0.6


;=== plot b test ===
; xy plane : at the moment of closest approach, 10Gyr
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e7_1Rs[k]-x_part_e7_1Rs[0]],[y_part_e7_1Rs[k]-y_part_e7_1Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_e7_1Rs[0,*]-x_part_e7_1Rs[0],Rs_xy_circle_e7_1Rs[1,*]-y_part_e7_1Rs[0]
oplot,Rtidal_xy_circle_e7_1Rs[0,*]-x_part_e7_1Rs[0],Rtidal_xy_circle_e7_1Rs[1,*]-y_part_e7_1Rs[0],linestyle=2
oplot,Rs_xy_circle_e7_2Rs[0,*]-x_part_e7_2Rs[0],Rs_xy_circle_e7_2Rs[1,*]-y_part_e7_2Rs[0],color=70
oplot,Rtidal_xy_circle_e7_2Rs[0,*]-x_part_e7_2Rs[0],Rtidal_xy_circle_e7_2Rs[1,*]-y_part_e7_2Rs[0],linestyle=2,color=70
oplot,Rs_xy_circle_e7_4Rs[0,*]-x_part_e7_4Rs[0],Rs_xy_circle_e7_4Rs[1,*]-y_part_e7_4Rs[0],color=150
oplot,Rtidal_xy_circle_e7_4Rs[0,*]-x_part_e7_4Rs[0],Rtidal_xy_circle_e7_4Rs[1,*]-y_part_e7_4Rs[0],linestyle=2,color=150

; xz plane
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e7_1Rs[k]-x_part_e7_1Rs[0]],[z_part_e7_1Rs[k]-z_part_e7_1Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_e7_1Rs[0,*]-x_part_e7_1Rs[0],Rs_xz_circle_e7_1Rs[1,*]-z_part_e7_1Rs[0]
oplot,Rtidal_xz_circle_e7_1Rs[0,*]-x_part_e7_1Rs[0],Rtidal_xz_circle_e7_1Rs[1,*]-z_part_e7_1Rs[0],linestyle=2
arrow,-8.3,-3,-8.3,-8,/data

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin-0.9,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_e7_1Rs_2[k]],[scale_dJ_part_e7_1Rs_2[k]],psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e7_2Rs_2[k]],[scale_dJ_part_e7_2Rs_2[k]]-0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e7_4Rs_2[k]],[scale_dJ_part_e7_4Rs_2[k]]-0.6,psym=3,color=qcolor[k]
xyouts,0,0,textoidl('b=1R_s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.3,textoidl('b=2R_s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.6,textoidl('b=4R_s'),alignment=0.5,charsize=0.6

; q histogram
plothist,q_part_e7_1Rs_2,xr=[qmin,qmax],bin=0.1,xtitle='q',yr=[0,240]
plothist,q_part_e7_2Rs_2,xr=[qmin,qmax],bin=0.1,color=70,/overplot;,/peak;yr=[0,210]
plothist,q_part_e7_4Rs_2,xr=[qmin,qmax],bin=0.1,color=150,/overplot;,/peak;yr=[0,210]
legend,['e7_1Rs=0.7452','e7_2Rs=0.3685','e7_4Rs=0.1860'],box=0,color=[0,70,255],psym=8,/right,charsize=0.6


;=== plot Msubhalo test ===
; xy plane : at the moment of closest approach, 10Gyr
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e8_1Rs[k]-x_part_e8_1Rs[0]],[y_part_e8_1Rs[k]-y_part_e8_1Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_e8_1Rs[0,*]-x_part_e8_1Rs[0],Rs_xy_circle_e8_1Rs[1,*]-y_part_e8_1Rs[0],color=70
oplot,Rtidal_xy_circle_e8_1Rs[0,*]-x_part_e8_1Rs[0],Rtidal_xy_circle_e8_1Rs[1,*]-y_part_e8_1Rs[0],linestyle=2,color=70
oplot,Rs_xy_circle_e8_2Rs[0,*]-x_part_e8_2Rs[0],Rs_xy_circle_e8_2Rs[1,*]-y_part_e8_2Rs[0],color=0
oplot,Rtidal_xy_circle_e8_2Rs[0,*]-x_part_e8_2Rs[0],Rtidal_xy_circle_e8_2Rs[1,*]-y_part_e8_2Rs[0],linestyle=2,color=0
oplot,Rs_xy_circle_e8_4Rs[0,*]-x_part_e8_4Rs[0],Rs_xy_circle_e8_4Rs[1,*]-y_part_e8_4Rs[0],color=150
oplot,Rtidal_xy_circle_e8_4Rs[0,*]-x_part_e8_4Rs[0],Rtidal_xy_circle_e8_4Rs[1,*]-y_part_e8_4Rs[0],linestyle=2,color=150

; xz plane
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e8_1Rs[k]-x_part_e8_1Rs[0]],[z_part_e8_1Rs[k]-z_part_e8_1Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_e8_1Rs[0,*]-x_part_e8_1Rs[0],Rs_xz_circle_e8_1Rs[1,*]-z_part_e8_1Rs[0],color=70
oplot,Rtidal_xz_circle_e8_1Rs[0,*]-x_part_e8_1Rs[0],Rtidal_xz_circle_e8_1Rs[1,*]-z_part_e8_1Rs[0],linestyle=2,color=70
arrow,-8.3,-3,-8.3,-8,/data
oplot,Rs_xz_circle_e8_2Rs[0,*]-x_part_e8_2Rs[0],Rs_xz_circle_e8_2Rs[1,*]-z_part_e8_2Rs[0],color=0
oplot,Rtidal_xz_circle_e8_2Rs[0,*]-x_part_e8_2Rs[0],Rtidal_xz_circle_e8_2Rs[1,*]-z_part_e8_2Rs[0],linestyle=2,color=0
oplot,Rs_xz_circle_e8_4Rs[0,*]-x_part_e8_4Rs[0],Rs_xz_circle_e8_4Rs[1,*]-z_part_e8_4Rs[0],color=150
oplot,Rtidal_xz_circle_e8_4Rs[0,*]-x_part_e8_4Rs[0],Rtidal_xz_circle_e8_4Rs[1,*]-z_part_e8_4Rs[0],linestyle=2,color=150

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin-2.0,qmax],yr=[dJmin-0.8,dJmax+0.3],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_e8_1Rs_2[k]],[scale_dJ_part_e8_1Rs_2[k]],psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e8_2Rs_2[k]],[scale_dJ_part_e8_2Rs_2[k]]-0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e8_4Rs_2[k]],[scale_dJ_part_e8_4Rs_2[k]]-0.3,psym=3,color=qcolor[k]
xyouts,0,0,textoidl('1Rs'),alignment=0.5,charsize=0.6
xyouts,0,0-0.3,textoidl('2Rs'),alignment=0.5,charsize=0.6
xyouts,0,0-0.6,textoidl('4Rs'),alignment=0.5,charsize=0.6

; q histogram
plothist,q_part_e8_1Rs_2,xr=[qmin-2.0,qmax],bin=0.1,xtitle='q',yr=[0,780]
plothist,q_part_e8_2Rs_2,bin=0.1,color=70,/overplot
plothist,q_part_e8_4Rs_2,bin=0.1,color=150,/overplot
legend,['e8_1Rs=2.7718','e8_2Rs=1.3898','e8_4Rs=0.6897'],box=0,color=[0,70,255],psym=8,/right,charsize=0.6

device,/close


END
