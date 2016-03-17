pro one_encounter

dir_out = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_nosub = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_nosub/'
dir_part_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v200/'
dir_subhalo_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v200/'
dir_part_v400 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v400/'
dir_subhalo_v400 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v400/'
dir_part_v800 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v800/'
dir_subhalo_v800 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v800/'


dir_part_b0Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v400/'
dir_subhalo_b0Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v400/'
dir_part_b1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_b_1Rs/'
dir_subhalo_b1Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_b_1Rs/'
dir_part_b2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_b_2Rs/'
dir_subhalo_b2Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_b_2Rs/'
dir_part_b4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_b_4Rs/'
dir_subhalo_b4Rs = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_b_4Rs/'


dir_part_e7 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v400/'
dir_subhalo_e7 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v400/'
dir_part_e6 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_e6/'
dir_subhalo_e6 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_e6/'
dir_part_e8 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_e8/'
dir_subhalo_e8 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_e8/'
dir_part_e9 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_e9/'
dir_subhalo_e9 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_e9/'


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

;=== v200 pericenter ===
;chr_rdtbl,dir_out+'part_peri_xy_v200',0, arr
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
chr_rdtbl, dir_part_v200+'snap00000',0, arr
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
N_time = N_elements(fname_part_v200)
Rs_e6 = 0.11058515
Rtidal_e6 = 0.87319397
Rs_e8 = 0.78954220
Rtidal_e8 = 16.358786
Rs_e9 = 2.1096704
Rtidal_e9 = 70.806194


;=== v200 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_v200+'snap01000',0,arr,/silent
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
scale_dJ_part_v200 = dJ_part_v200 / (r_tide/r_peri * J_part_v200)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v200+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_v200 = reform(arr[1,*])
y_subhalo_v200 = reform(arr[2,*])
z_subhalo_v200 = reform(arr[3,*])
Vx_subhalo_v200 = reform(arr[4,*])
Vy_subhalo_v200 = reform(arr[5,*])
Vz_subhalo_v200 = reform(arr[6,*])

t_step = string(t_part_v200[0]/1000.,f='(f7.4)')

; make a circle of a subhalo 
Rs_xy_circle_v200 = circle(x_subhalo_v200[0],y_subhalo_v200[0],double(Rs[0]))
Rtidal_xy_circle_v200 = circle(x_subhalo_v200[0],y_subhalo_v200[0],double(Rtidal[0]))
Rs_yz_circle_v200 = circle(y_subhalo_v200[0],z_subhalo_v200[0],double(Rs[0]))
Rtidal_yz_circle_v200 = circle(y_subhalo_v200[0],z_subhalo_v200[0],double(Rtidal[0]))
Rs_xz_circle_v200 = circle(x_subhalo_v200[0],z_subhalo_v200[0],double(Rs[0]))
Rtidal_xz_circle_v200 = circle(x_subhalo_v200[0],z_subhalo_v200[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_v200+'snap01020',0,arr,/silent
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
scale_dJ_part_v200_2 = dJ_part_v200_2 / (r_tide/r_peri * J_part_v200_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v200+'subhalo01020',0,arr,/silent
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

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_v200+'subhalo00980',0,arr,/silent
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


;=== v400 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_v400+'snap01000',0,arr,/silent
arr = double(arr)
t_part_v400 = reform(arr[0,*])
x_part_v400 = reform(arr[1,*])
y_part_v400 = reform(arr[2,*])
z_part_v400 = reform(arr[3,*])
Vx_part_v400 = reform(arr[4,*])
Vy_part_v400 = reform(arr[5,*])
Vz_part_v400 = reform(arr[6,*])
dE_part_v400 = reform(arr[7,*])
Etot_part_v400 = reform(arr[8,*])
J_part_v400 = reform(arr[9,*])
delE = Etot_part_v400-Etot_part_v400[0]
q_part_v400 = delE/epsilon
dJ_part_v400 = J_part_v400-J_part_v400[0]
scale_dJ_part_v400 = dJ_part_v400 / (r_tide/r_peri * J_part_v400)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v400+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_v400 = reform(arr[1,*])
y_subhalo_v400 = reform(arr[2,*])
z_subhalo_v400 = reform(arr[3,*])
Vx_subhalo_v400 = reform(arr[4,*])
Vy_subhalo_v400 = reform(arr[5,*])
Vz_subhalo_v400 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v400 = circle(x_subhalo_v400[0],y_subhalo_v400[0],double(Rs[0]))
Rtidal_xy_circle_v400 = circle(x_subhalo_v400[0],y_subhalo_v400[0],double(Rtidal[0]))
Rs_yz_circle_v400 = circle(y_subhalo_v400[0],z_subhalo_v400[0],double(Rs[0]))
Rtidal_yz_circle_v400 = circle(y_subhalo_v400[0],z_subhalo_v400[0],double(Rtidal[0]))
Rs_xz_circle_v400 = circle(x_subhalo_v400[0],z_subhalo_v400[0],double(Rs[0]))
Rtidal_xz_circle_v400 = circle(x_subhalo_v400[0],z_subhalo_v400[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_v400+'snap01020',0,arr,/silent
arr = double(arr)
t_part_v400_2 = reform(arr[0,*])
x_part_v400_2 = reform(arr[1,*])
y_part_v400_2 = reform(arr[2,*])
z_part_v400_2 = reform(arr[3,*])
Vx_part_v400_2 = reform(arr[4,*])
Vy_part_v400_2 = reform(arr[5,*])
Vz_part_v400_2 = reform(arr[6,*])
dE_part_v400_2 = reform(arr[7,*])
Etot_part_v400_2 = reform(arr[8,*])
J_part_v400_2 = reform(arr[9,*])
delE = Etot_part_v400_2-Etot_part_v400_2[0]
q_part_v400_2 = delE/epsilon
dJ_part_v400_2 = J_part_v400_2-J_part_v400_2[0]
scale_dJ_part_v400_2 = dJ_part_v400_2 / (r_tide/r_peri * J_part_v400_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v400+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_v400_2 = reform(arr[1,*])
y_subhalo_v400_2 = reform(arr[2,*])
z_subhalo_v400_2 = reform(arr[3,*])
Vx_subhalo_v400_2 = reform(arr[4,*])
Vy_subhalo_v400_2 = reform(arr[5,*])
Vz_subhalo_v400_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v400_2 = circle(x_subhalo_v400_2[0],y_subhalo_v400_2[0],double(Rs[0]))
Rtidal_xy_circle_v400_2 = circle(x_subhalo_v400_2[0],y_subhalo_v400_2[0],double(Rtidal[0]))
Rs_yz_circle_v400_2 = circle(y_subhalo_v400_2[0],z_subhalo_v400_2[0],double(Rs[0]))
Rtidal_yz_circle_v400_2 = circle(y_subhalo_v400_2[0],z_subhalo_v400_2[0],double(Rtidal[0]))
Rs_xz_circle_v400_2 = circle(x_subhalo_v400_2[0],z_subhalo_v400_2[0],double(Rs[0]))
Rtidal_xz_circle_v400_2 = circle(x_subhalo_v400_2[0],z_subhalo_v400_2[0],double(Rtidal[0]))

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_v400+'subhalo00980',0,arr,/silent
arr = double(arr)
x_subhalo_v400_3 = reform(arr[1,*])
y_subhalo_v400_3 = reform(arr[2,*])
z_subhalo_v400_3 = reform(arr[3,*])
Vx_subhalo_v400_3 = reform(arr[4,*])
Vy_subhalo_v400_3 = reform(arr[5,*])
Vz_subhalo_v400_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v400_3 = circle(x_subhalo_v400_3[0],y_subhalo_v400_3[0],double(Rs[0]))
Rtidal_xy_circle_v400_3 = circle(x_subhalo_v400_3[0],y_subhalo_v400_3[0],double(Rtidal[0]))
Rs_yz_circle_v400_3 = circle(y_subhalo_v400_3[0],z_subhalo_v400_3[0],double(Rs[0]))
Rtidal_yz_circle_v400_3 = circle(y_subhalo_v400_3[0],z_subhalo_v400_3[0],double(Rtidal[0]))
Rs_xz_circle_v400_3 = circle(x_subhalo_v400_3[0],z_subhalo_v400_3[0],double(Rs[0]))
Rtidal_xz_circle_v400_3 = circle(x_subhalo_v400_3[0],z_subhalo_v400_3[0],double(Rtidal[0]))


;=== v800 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_v800+'snap01000',0,arr,/silent
arr = double(arr)
t_part_v800 = reform(arr[0,*])
x_part_v800 = reform(arr[1,*])
y_part_v800 = reform(arr[2,*])
z_part_v800 = reform(arr[3,*])
Vx_part_v800 = reform(arr[4,*])
Vy_part_v800 = reform(arr[5,*])
Vz_part_v800 = reform(arr[6,*])
dE_part_v800 = reform(arr[7,*])
Etot_part_v800 = reform(arr[8,*])
J_part_v800 = reform(arr[9,*])
delE = Etot_part_v800-Etot_part_v800[0]
q_part_v800 = delE/epsilon
dJ_part_v800 = J_part_v800-J_part_v800[0]
scale_dJ_part_v800 = dJ_part_v800 / (r_tide/r_peri * J_part_v800)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v800+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_v800 = reform(arr[1,*])
y_subhalo_v800 = reform(arr[2,*])
z_subhalo_v800 = reform(arr[3,*])
Vx_subhalo_v800 = reform(arr[4,*])
Vy_subhalo_v800 = reform(arr[5,*])
Vz_subhalo_v800 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v800 = circle(x_subhalo_v800[0],y_subhalo_v800[0],double(Rs[0]))
Rtidal_xy_circle_v800 = circle(x_subhalo_v800[0],y_subhalo_v800[0],double(Rtidal[0]))
Rs_yz_circle_v800 = circle(y_subhalo_v800[0],z_subhalo_v800[0],double(Rs[0]))
Rtidal_yz_circle_v800 = circle(y_subhalo_v800[0],z_subhalo_v800[0],double(Rtidal[0]))
Rs_xz_circle_v800 = circle(x_subhalo_v800[0],z_subhalo_v800[0],double(Rs[0]))
Rtidal_xz_circle_v800 = circle(x_subhalo_v800[0],z_subhalo_v800[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_v800+'snap01020',0,arr,/silent
arr = double(arr)
t_part_v800_2 = reform(arr[0,*])
x_part_v800_2 = reform(arr[1,*])
y_part_v800_2 = reform(arr[2,*])
z_part_v800_2 = reform(arr[3,*])
Vx_part_v800_2 = reform(arr[4,*])
Vy_part_v800_2 = reform(arr[5,*])
Vz_part_v800_2 = reform(arr[6,*])
dE_part_v800_2 = reform(arr[7,*])
Etot_part_v800_2 = reform(arr[8,*])
J_part_v800_2 = reform(arr[9,*])
delE = Etot_part_v800_2-Etot_part_v800_2[0]
q_part_v800_2 = delE/epsilon
dJ_part_v800_2 = J_part_v800_2-J_part_v800_2[0]
scale_dJ_part_v800_2 = dJ_part_v800_2 / (r_tide/r_peri * J_part_v800_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_v800+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_v800_2 = reform(arr[1,*])
y_subhalo_v800_2 = reform(arr[2,*])
z_subhalo_v800_2 = reform(arr[3,*])
Vx_subhalo_v800_2 = reform(arr[4,*])
Vy_subhalo_v800_2 = reform(arr[5,*])
Vz_subhalo_v800_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v800_2 = circle(x_subhalo_v800_2[0],y_subhalo_v800_2[0],double(Rs[0]))
Rtidal_xy_circle_v800_2 = circle(x_subhalo_v800_2[0],y_subhalo_v800_2[0],double(Rtidal[0]))
Rs_yz_circle_v800_2 = circle(y_subhalo_v800_2[0],z_subhalo_v800_2[0],double(Rs[0]))
Rtidal_yz_circle_v800_2 = circle(y_subhalo_v800_2[0],z_subhalo_v800_2[0],double(Rtidal[0]))
Rs_xz_circle_v800_2 = circle(x_subhalo_v800_2[0],z_subhalo_v800_2[0],double(Rs[0]))
Rtidal_xz_circle_v800_2 = circle(x_subhalo_v800_2[0],z_subhalo_v800_2[0],double(Rtidal[0]))

; read out data of 20Myr before
chr_rdtbl,dir_subhalo_v800+'subhalo00980',0,arr,/silent
arr = double(arr)
x_subhalo_v800_3 = reform(arr[1,*])
y_subhalo_v800_3 = reform(arr[2,*])
z_subhalo_v800_3 = reform(arr[3,*])
Vx_subhalo_v800_3 = reform(arr[4,*])
Vy_subhalo_v800_3 = reform(arr[5,*])
Vz_subhalo_v800_3 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_v800_3 = circle(x_subhalo_v800_3[0],y_subhalo_v800_3[0],double(Rs[0]))
Rtidal_xy_circle_v800_3 = circle(x_subhalo_v800_3[0],y_subhalo_v800_3[0],double(Rtidal[0]))
Rs_yz_circle_v800_3 = circle(y_subhalo_v800_3[0],z_subhalo_v800_3[0],double(Rs[0]))
Rtidal_yz_circle_v800_3 = circle(y_subhalo_v800_3[0],z_subhalo_v800_3[0],double(Rtidal[0]))
Rs_xz_circle_v800_3 = circle(x_subhalo_v800_3[0],z_subhalo_v800_3[0],double(Rs[0]))
Rtidal_xz_circle_v800_3 = circle(x_subhalo_v800_3[0],z_subhalo_v800_3[0],double(Rtidal[0]))


;=====================
;=== b_impact test ===
;=====================
;=== b=0 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_b0Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_b0Rs = reform(arr[0,*])
x_part_b0Rs = reform(arr[1,*])
y_part_b0Rs = reform(arr[2,*])
z_part_b0Rs = reform(arr[3,*])
Vx_part_b0Rs = reform(arr[4,*])
Vy_part_b0Rs = reform(arr[5,*])
Vz_part_b0Rs = reform(arr[6,*])
dE_part_b0Rs = reform(arr[7,*])
Etot_part_b0Rs = reform(arr[8,*])
J_part_b0Rs = reform(arr[9,*])
delE = Etot_part_b0Rs-Etot_part_b0Rs[0]
q_part_b0Rs = delE/epsilon
dJ_part_b0Rs = J_part_b0Rs-J_part_b0Rs[0]
scale_dJ_part_b0Rs = dJ_part_b0Rs / (r_tide/r_peri * J_part_b0Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_b0Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_b0Rs = reform(arr[1,*])
y_subhalo_b0Rs = reform(arr[2,*])
z_subhalo_b0Rs = reform(arr[3,*])
Vx_subhalo_b0Rs = reform(arr[4,*])
Vy_subhalo_b0Rs = reform(arr[5,*])
Vz_subhalo_b0Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_b0Rs = circle(x_subhalo_b0Rs[0],y_subhalo_b0Rs[0],double(Rs[0]))
Rtidal_xy_circle_b0Rs = circle(x_subhalo_b0Rs[0],y_subhalo_b0Rs[0],double(Rtidal[0]))
Rs_yz_circle_b0Rs = circle(y_subhalo_b0Rs[0],z_subhalo_b0Rs[0],double(Rs[0]))
Rtidal_yz_circle_b0Rs = circle(y_subhalo_b0Rs[0],z_subhalo_b0Rs[0],double(Rtidal[0]))
Rs_xz_circle_b0Rs = circle(x_subhalo_b0Rs[0],z_subhalo_b0Rs[0],double(Rs[0]))
Rtidal_xz_circle_b0Rs = circle(x_subhalo_b0Rs[0],z_subhalo_b0Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_b0Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_b0Rs_2 = reform(arr[0,*])
x_part_b0Rs_2 = reform(arr[1,*])
y_part_b0Rs_2 = reform(arr[2,*])
z_part_b0Rs_2 = reform(arr[3,*])
Vx_part_b0Rs_2 = reform(arr[4,*])
Vy_part_b0Rs_2 = reform(arr[5,*])
Vz_part_b0Rs_2 = reform(arr[6,*])
dE_part_b0Rs_2 = reform(arr[7,*])
Etot_part_b0Rs_2 = reform(arr[8,*])
J_part_b0Rs_2 = reform(arr[9,*])
delE = Etot_part_b0Rs_2-Etot_part_b0Rs_2[0]
q_part_b0Rs_2 = delE/epsilon
dJ_part_b0Rs_2 = J_part_b0Rs_2-J_part_b0Rs_2[0]
scale_dJ_part_b0Rs_2 = dJ_part_b0Rs_2 / (r_tide/r_peri * J_part_b0Rs_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_b0Rs+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_b0Rs_2 = reform(arr[1,*])
y_subhalo_b0Rs_2 = reform(arr[2,*])
z_subhalo_b0Rs_2 = reform(arr[3,*])
Vx_subhalo_b0Rs_2 = reform(arr[4,*])
Vy_subhalo_b0Rs_2 = reform(arr[5,*])
Vz_subhalo_b0Rs_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_b0Rs_2 = circle(x_subhalo_b0Rs_2[0],y_subhalo_b0Rs_2[0],double(Rs[0]))
Rtidal_xy_circle_b0Rs_2 = circle(x_subhalo_b0Rs_2[0],y_subhalo_b0Rs_2[0],double(Rtidal[0]))
Rs_yz_circle_b0Rs_2 = circle(y_subhalo_b0Rs_2[0],z_subhalo_b0Rs_2[0],double(Rs[0]))
Rtidal_yz_circle_b0Rs_2 = circle(y_subhalo_b0Rs_2[0],z_subhalo_b0Rs_2[0],double(Rtidal[0]))
Rs_xz_circle_b0Rs_2 = circle(x_subhalo_b0Rs_2[0],z_subhalo_b0Rs_2[0],double(Rs[0]))
Rtidal_xz_circle_b0Rs_2 = circle(x_subhalo_b0Rs_2[0],z_subhalo_b0Rs_2[0],double(Rtidal[0]))


;=== b=1Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_b1Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_b1Rs = reform(arr[0,*])
x_part_b1Rs = reform(arr[1,*])
y_part_b1Rs = reform(arr[2,*])
z_part_b1Rs = reform(arr[3,*])
Vx_part_b1Rs = reform(arr[4,*])
Vy_part_b1Rs = reform(arr[5,*])
Vz_part_b1Rs = reform(arr[6,*])
dE_part_b1Rs = reform(arr[7,*])
Etot_part_b1Rs = reform(arr[8,*])
J_part_b1Rs = reform(arr[9,*])
delE = Etot_part_b1Rs-Etot_part_b1Rs[0]
q_part_b1Rs = delE/epsilon
dJ_part_b1Rs = J_part_b1Rs-J_part_b1Rs[0]
scale_dJ_part_b1Rs = dJ_part_b1Rs / (r_tide/r_peri * J_part_b1Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_b1Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_b1Rs = reform(arr[1,*])
y_subhalo_b1Rs = reform(arr[2,*])
z_subhalo_b1Rs = reform(arr[3,*])
Vx_subhalo_b1Rs = reform(arr[4,*])
Vy_subhalo_b1Rs = reform(arr[5,*])
Vz_subhalo_b1Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_b1Rs = circle(x_subhalo_b1Rs[0],y_subhalo_b1Rs[0],double(Rs[0]))
Rtidal_xy_circle_b1Rs = circle(x_subhalo_b1Rs[0],y_subhalo_b1Rs[0],double(Rtidal[0]))
Rs_yz_circle_b1Rs = circle(y_subhalo_b1Rs[0],z_subhalo_b1Rs[0],double(Rs[0]))
Rtidal_yz_circle_b1Rs = circle(y_subhalo_b1Rs[0],z_subhalo_b1Rs[0],double(Rtidal[0]))
Rs_xz_circle_b1Rs = circle(x_subhalo_b1Rs[0],z_subhalo_b1Rs[0],double(Rs[0]))
Rtidal_xz_circle_b1Rs = circle(x_subhalo_b1Rs[0],z_subhalo_b1Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_b1Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_b1Rs_2 = reform(arr[0,*])
x_part_b1Rs_2 = reform(arr[1,*])
y_part_b1Rs_2 = reform(arr[2,*])
z_part_b1Rs_2 = reform(arr[3,*])
Vx_part_b1Rs_2 = reform(arr[4,*])
Vy_part_b1Rs_2 = reform(arr[5,*])
Vz_part_b1Rs_2 = reform(arr[6,*])
dE_part_b1Rs_2 = reform(arr[7,*])
Etot_part_b1Rs_2 = reform(arr[8,*])
J_part_b1Rs_2 = reform(arr[9,*])
delE = Etot_part_b1Rs_2-Etot_part_b1Rs_2[0]
q_part_b1Rs_2 = delE/epsilon
dJ_part_b1Rs_2 = J_part_b1Rs_2-J_part_b1Rs_2[0]
scale_dJ_part_b1Rs_2 = dJ_part_b1Rs_2 / (r_tide/r_peri * J_part_b1Rs_2)


;=== b=2Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_b2Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_b2Rs = reform(arr[0,*])
x_part_b2Rs = reform(arr[1,*])
y_part_b2Rs = reform(arr[2,*])
z_part_b2Rs = reform(arr[3,*])
Vx_part_b2Rs = reform(arr[4,*])
Vy_part_b2Rs = reform(arr[5,*])
Vz_part_b2Rs = reform(arr[6,*])
dE_part_b2Rs = reform(arr[7,*])
Etot_part_b2Rs = reform(arr[8,*])
J_part_b2Rs = reform(arr[9,*])
delE = Etot_part_b2Rs-Etot_part_b2Rs[0]
q_part_b2Rs = delE/epsilon
dJ_part_b2Rs = J_part_b2Rs-J_part_b2Rs[0]
scale_dJ_part_b2Rs = dJ_part_b2Rs / (r_tide/r_peri * J_part_b2Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_b2Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_b2Rs = reform(arr[1,*])
y_subhalo_b2Rs = reform(arr[2,*])
z_subhalo_b2Rs = reform(arr[3,*])
Vx_subhalo_b2Rs = reform(arr[4,*])
Vy_subhalo_b2Rs = reform(arr[5,*])
Vz_subhalo_b2Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_b2Rs = circle(x_subhalo_b2Rs[0],y_subhalo_b2Rs[0],double(Rs[0]))
Rtidal_xy_circle_b2Rs = circle(x_subhalo_b2Rs[0],y_subhalo_b2Rs[0],double(Rtidal[0]))
Rs_yz_circle_b2Rs = circle(y_subhalo_b2Rs[0],z_subhalo_b2Rs[0],double(Rs[0]))
Rtidal_yz_circle_b2Rs = circle(y_subhalo_b2Rs[0],z_subhalo_b2Rs[0],double(Rtidal[0]))
Rs_xz_circle_b2Rs = circle(x_subhalo_b2Rs[0],z_subhalo_b2Rs[0],double(Rs[0]))
Rtidal_xz_circle_b2Rs = circle(x_subhalo_b2Rs[0],z_subhalo_b2Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_b2Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_b2Rs_2 = reform(arr[0,*])
x_part_b2Rs_2 = reform(arr[1,*])
y_part_b2Rs_2 = reform(arr[2,*])
z_part_b2Rs_2 = reform(arr[3,*])
Vx_part_b2Rs_2 = reform(arr[4,*])
Vy_part_b2Rs_2 = reform(arr[5,*])
Vz_part_b2Rs_2 = reform(arr[6,*])
dE_part_b2Rs_2 = reform(arr[7,*])
Etot_part_b2Rs_2 = reform(arr[8,*])
J_part_b2Rs_2 = reform(arr[9,*])
delE = Etot_part_b2Rs_2-Etot_part_b2Rs_2[0]
q_part_b2Rs_2 = delE/epsilon
dJ_part_b2Rs_2 = J_part_b2Rs_2-J_part_b2Rs_2[0]
scale_dJ_part_b2Rs_2 = dJ_part_b2Rs_2 / (r_tide/r_peri * J_part_b2Rs_2)


;=== b=4Rs ===
; read out data of the tail particles 
chr_rdtbl,dir_part_b4Rs+'snap01000',0,arr,/silent
arr = double(arr)
t_part_b4Rs = reform(arr[0,*])
x_part_b4Rs = reform(arr[1,*])
y_part_b4Rs = reform(arr[2,*])
z_part_b4Rs = reform(arr[3,*])
Vx_part_b4Rs = reform(arr[4,*])
Vy_part_b4Rs = reform(arr[5,*])
Vz_part_b4Rs = reform(arr[6,*])
dE_part_b4Rs = reform(arr[7,*])
Etot_part_b4Rs = reform(arr[8,*])
J_part_b4Rs = reform(arr[9,*])
delE = Etot_part_b4Rs-Etot_part_b4Rs[0]
q_part_b4Rs = delE/epsilon
dJ_part_b4Rs = J_part_b4Rs-J_part_b4Rs[0]
scale_dJ_part_b4Rs = dJ_part_b4Rs / (r_tide/r_peri * J_part_b4Rs)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_b4Rs+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_b4Rs = reform(arr[1,*])
y_subhalo_b4Rs = reform(arr[2,*])
z_subhalo_b4Rs = reform(arr[3,*])
Vx_subhalo_b4Rs = reform(arr[4,*])
Vy_subhalo_b4Rs = reform(arr[5,*])
Vz_subhalo_b4Rs = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_b4Rs = circle(x_subhalo_b4Rs[0],y_subhalo_b4Rs[0],double(Rs[0]))
Rtidal_xy_circle_b4Rs = circle(x_subhalo_b4Rs[0],y_subhalo_b4Rs[0],double(Rtidal[0]))
Rs_yz_circle_b4Rs = circle(y_subhalo_b4Rs[0],z_subhalo_b4Rs[0],double(Rs[0]))
Rtidal_yz_circle_b4Rs = circle(y_subhalo_b4Rs[0],z_subhalo_b4Rs[0],double(Rtidal[0]))
Rs_xz_circle_b4Rs = circle(x_subhalo_b4Rs[0],z_subhalo_b4Rs[0],double(Rs[0]))
Rtidal_xz_circle_b4Rs = circle(x_subhalo_b4Rs[0],z_subhalo_b4Rs[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_b4Rs+'snap01020',0,arr,/silent
arr = double(arr)
t_part_b4Rs_2 = reform(arr[0,*])
x_part_b4Rs_2 = reform(arr[1,*])
y_part_b4Rs_2 = reform(arr[2,*])
z_part_b4Rs_2 = reform(arr[3,*])
Vx_part_b4Rs_2 = reform(arr[4,*])
Vy_part_b4Rs_2 = reform(arr[5,*])
Vz_part_b4Rs_2 = reform(arr[6,*])
dE_part_b4Rs_2 = reform(arr[7,*])
Etot_part_b4Rs_2 = reform(arr[8,*])
J_part_b4Rs_2 = reform(arr[9,*])
delE = Etot_part_b4Rs_2-Etot_part_b4Rs_2[0]
q_part_b4Rs_2 = delE/epsilon
dJ_part_b4Rs_2 = J_part_b4Rs_2-J_part_b4Rs_2[0]
scale_dJ_part_b4Rs_2 = dJ_part_b4Rs_2 / (r_tide/r_peri * J_part_b4Rs_2)


;=====================
;=== Msubhalo test ===
;=====================
;=== Msubhalo=e7 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e7+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e7 = reform(arr[0,*])
x_part_e7 = reform(arr[1,*])
y_part_e7 = reform(arr[2,*])
z_part_e7 = reform(arr[3,*])
Vx_part_e7 = reform(arr[4,*])
Vy_part_e7 = reform(arr[5,*])
Vz_part_e7 = reform(arr[6,*])
dE_part_e7 = reform(arr[7,*])
Etot_part_e7 = reform(arr[8,*])
J_part_e7 = reform(arr[9,*])
delE = Etot_part_e7-Etot_part_e7[0]
q_part_e7 = delE/epsilon
dJ_part_e7 = J_part_e7-J_part_e7[0]
scale_dJ_part_e7 = dJ_part_e7 / (r_tide/r_peri * J_part_e7)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e7+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e7 = reform(arr[1,*])
y_subhalo_e7 = reform(arr[2,*])
z_subhalo_e7 = reform(arr[3,*])
Vx_subhalo_e7 = reform(arr[4,*])
Vy_subhalo_e7 = reform(arr[5,*])
Vz_subhalo_e7 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e7 = circle(x_subhalo_e7[0],y_subhalo_e7[0],double(Rs[0]))
Rtidal_xy_circle_e7 = circle(x_subhalo_e7[0],y_subhalo_e7[0],double(Rtidal[0]))
Rs_yz_circle_e7 = circle(y_subhalo_e7[0],z_subhalo_e7[0],double(Rs[0]))
Rtidal_yz_circle_e7 = circle(y_subhalo_e7[0],z_subhalo_e7[0],double(Rtidal[0]))
Rs_xz_circle_e7 = circle(x_subhalo_e7[0],z_subhalo_e7[0],double(Rs[0]))
Rtidal_xz_circle_e7 = circle(x_subhalo_e7[0],z_subhalo_e7[0],double(Rtidal[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e7+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e7_2 = reform(arr[0,*])
x_part_e7_2 = reform(arr[1,*])
y_part_e7_2 = reform(arr[2,*])
z_part_e7_2 = reform(arr[3,*])
Vx_part_e7_2 = reform(arr[4,*])
Vy_part_e7_2 = reform(arr[5,*])
Vz_part_e7_2 = reform(arr[6,*])
dE_part_e7_2 = reform(arr[7,*])
Etot_part_e7_2 = reform(arr[8,*])
J_part_e7_2 = reform(arr[9,*])
delE = Etot_part_e7_2-Etot_part_e7_2[0]
q_part_e7_2 = delE/epsilon
dJ_part_e7_2 = J_part_e7_2-J_part_e7_2[0]
scale_dJ_part_e7_2 = dJ_part_e7_2 / (r_tide/r_peri * J_part_e7_2)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e7+'subhalo01020',0,arr,/silent
arr = double(arr)
x_subhalo_e7_2 = reform(arr[1,*])
y_subhalo_e7_2 = reform(arr[2,*])
z_subhalo_e7_2 = reform(arr[3,*])
Vx_subhalo_e7_2 = reform(arr[4,*])
Vy_subhalo_e7_2 = reform(arr[5,*])
Vz_subhalo_e7_2 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e7_2 = circle(x_subhalo_e7_2[0],y_subhalo_e7_2[0],double(Rs[0]))
Rtidal_xy_circle_e7_2 = circle(x_subhalo_e7_2[0],y_subhalo_e7_2[0],double(Rtidal[0]))
Rs_yz_circle_e7_2 = circle(y_subhalo_e7_2[0],z_subhalo_e7_2[0],double(Rs[0]))
Rtidal_yz_circle_e7_2 = circle(y_subhalo_e7_2[0],z_subhalo_e7_2[0],double(Rtidal[0]))
Rs_xz_circle_e7_2 = circle(x_subhalo_e7_2[0],z_subhalo_e7_2[0],double(Rs[0]))
Rtidal_xz_circle_e7_2 = circle(x_subhalo_e7_2[0],z_subhalo_e7_2[0],double(Rtidal[0]))


;=== Msubhalo=e6 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e6+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e6 = reform(arr[0,*])
x_part_e6 = reform(arr[1,*])
y_part_e6 = reform(arr[2,*])
z_part_e6 = reform(arr[3,*])
Vx_part_e6 = reform(arr[4,*])
Vy_part_e6 = reform(arr[5,*])
Vz_part_e6 = reform(arr[6,*])
dE_part_e6 = reform(arr[7,*])
Etot_part_e6 = reform(arr[8,*])
J_part_e6 = reform(arr[9,*])
delE = Etot_part_e6-Etot_part_e6[0]
q_part_e6 = delE/epsilon
dJ_part_e6 = J_part_e6-J_part_e6[0]
scale_dJ_part_e6 = dJ_part_e6 / (r_tide/r_peri * J_part_e6)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e6+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e6 = reform(arr[1,*])
y_subhalo_e6 = reform(arr[2,*])
z_subhalo_e6 = reform(arr[3,*])
Vx_subhalo_e6 = reform(arr[4,*])
Vy_subhalo_e6 = reform(arr[5,*])
Vz_subhalo_e6 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e6 = circle(x_subhalo_e6[0],y_subhalo_e6[0],double(Rs_e6[0]))
Rtidal_xy_circle_e6 = circle(x_subhalo_e6[0],y_subhalo_e6[0],double(Rtidal_e6[0]))
Rs_yz_circle_e6 = circle(y_subhalo_e6[0],z_subhalo_e6[0],double(Rs_e6[0]))
Rtidal_yz_circle_e6 = circle(y_subhalo_e6[0],z_subhalo_e6[0],double(Rtidal_e6[0]))
Rs_xz_circle_e6 = circle(x_subhalo_e6[0],z_subhalo_e6[0],double(Rs_e6[0]))
Rtidal_xz_circle_e6 = circle(x_subhalo_e6[0],z_subhalo_e6[0],double(Rtidal_e6[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e6+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e6_2 = reform(arr[0,*])
x_part_e6_2 = reform(arr[1,*])
y_part_e6_2 = reform(arr[2,*])
z_part_e6_2 = reform(arr[3,*])
Vx_part_e6_2 = reform(arr[4,*])
Vy_part_e6_2 = reform(arr[5,*])
Vz_part_e6_2 = reform(arr[6,*])
dE_part_e6_2 = reform(arr[7,*])
Etot_part_e6_2 = reform(arr[8,*])
J_part_e6_2 = reform(arr[9,*])
delE = Etot_part_e6_2-Etot_part_e6_2[0]
q_part_e6_2 = delE/epsilon
dJ_part_e6_2 = J_part_e6_2-J_part_e6_2[0]
scale_dJ_part_e6_2 = dJ_part_e6_2 / (r_tide/r_peri * J_part_e6_2)


;=== Msubhalo=e8 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e8+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e8 = reform(arr[0,*])
x_part_e8 = reform(arr[1,*])
y_part_e8 = reform(arr[2,*])
z_part_e8 = reform(arr[3,*])
Vx_part_e8 = reform(arr[4,*])
Vy_part_e8 = reform(arr[5,*])
Vz_part_e8 = reform(arr[6,*])
dE_part_e8 = reform(arr[7,*])
Etot_part_e8 = reform(arr[8,*])
J_part_e8 = reform(arr[9,*])
delE = Etot_part_e8-Etot_part_e8[0]
q_part_e8 = delE/epsilon
dJ_part_e8 = J_part_e8-J_part_e8[0]
scale_dJ_part_e8 = dJ_part_e8 / (r_tide/r_peri * J_part_e8)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e8+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e8 = reform(arr[1,*])
y_subhalo_e8 = reform(arr[2,*])
z_subhalo_e8 = reform(arr[3,*])
Vx_subhalo_e8 = reform(arr[4,*])
Vy_subhalo_e8 = reform(arr[5,*])
Vz_subhalo_e8 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e8 = circle(x_subhalo_e8[0],y_subhalo_e8[0],double(Rs_e8[0]))
Rtidal_xy_circle_e8 = circle(x_subhalo_e8[0],y_subhalo_e8[0],double(Rtidal_e8[0]))
Rs_yz_circle_e8 = circle(y_subhalo_e8[0],z_subhalo_e8[0],double(Rs_e8[0]))
Rtidal_yz_circle_e8 = circle(y_subhalo_e8[0],z_subhalo_e8[0],double(Rtidal_e8[0]))
Rs_xz_circle_e8 = circle(x_subhalo_e8[0],z_subhalo_e8[0],double(Rs_e8[0]))
Rtidal_xz_circle_e8 = circle(x_subhalo_e8[0],z_subhalo_e8[0],double(Rtidal_e8[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e8+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e8_2 = reform(arr[0,*])
x_part_e8_2 = reform(arr[1,*])
y_part_e8_2 = reform(arr[2,*])
z_part_e8_2 = reform(arr[3,*])
Vx_part_e8_2 = reform(arr[4,*])
Vy_part_e8_2 = reform(arr[5,*])
Vz_part_e8_2 = reform(arr[6,*])
dE_part_e8_2 = reform(arr[7,*])
Etot_part_e8_2 = reform(arr[8,*])
J_part_e8_2 = reform(arr[9,*])
delE = Etot_part_e8_2-Etot_part_e8_2[0]
q_part_e8_2 = delE/epsilon
dJ_part_e8_2 = J_part_e8_2-J_part_e8_2[0]
scale_dJ_part_e8_2 = dJ_part_e8_2 / (r_tide/r_peri * J_part_e8_2)


;=== Msubhalo=e9 ===
; read out data of the tail particles 
chr_rdtbl,dir_part_e9+'snap01000',0,arr,/silent
arr = double(arr)
t_part_e9 = reform(arr[0,*])
x_part_e9 = reform(arr[1,*])
y_part_e9 = reform(arr[2,*])
z_part_e9 = reform(arr[3,*])
Vx_part_e9 = reform(arr[4,*])
Vy_part_e9 = reform(arr[5,*])
Vz_part_e9 = reform(arr[6,*])
dE_part_e9 = reform(arr[7,*])
Etot_part_e9 = reform(arr[8,*])
J_part_e9 = reform(arr[9,*])
delE = Etot_part_e9-Etot_part_e9[0]
q_part_e9 = delE/epsilon
dJ_part_e9 = J_part_e9-J_part_e9[0]
scale_dJ_part_e9 = dJ_part_e9 / (r_tide/r_peri * J_part_e9)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_e9+'subhalo01000',0,arr,/silent
arr = double(arr)
x_subhalo_e9 = reform(arr[1,*])
y_subhalo_e9 = reform(arr[2,*])
z_subhalo_e9 = reform(arr[3,*])
Vx_subhalo_e9 = reform(arr[4,*])
Vy_subhalo_e9 = reform(arr[5,*])
Vz_subhalo_e9 = reform(arr[6,*])

; make a circle of a subhalo 
Rs_xy_circle_e9 = circle(x_subhalo_e9[0],y_subhalo_e9[0],double(Rs_e9[0]))
Rtidal_xy_circle_e9 = circle(x_subhalo_e9[0],y_subhalo_e9[0],double(Rtidal_e9[0]))
Rs_yz_circle_e9 = circle(y_subhalo_e9[0],z_subhalo_e9[0],double(Rs_e9[0]))
Rtidal_yz_circle_e9 = circle(y_subhalo_e9[0],z_subhalo_e9[0],double(Rtidal_e9[0]))
Rs_xz_circle_e9 = circle(x_subhalo_e9[0],z_subhalo_e9[0],double(Rs_e9[0]))
Rtidal_xz_circle_e9 = circle(x_subhalo_e9[0],z_subhalo_e9[0],double(Rtidal_e9[0]))

; read data of 20Myr later 
chr_rdtbl,dir_part_e9+'snap01020',0,arr,/silent
arr = double(arr)
t_part_e9_2 = reform(arr[0,*])
x_part_e9_2 = reform(arr[1,*])
y_part_e9_2 = reform(arr[2,*])
z_part_e9_2 = reform(arr[3,*])
Vx_part_e9_2 = reform(arr[4,*])
Vy_part_e9_2 = reform(arr[5,*])
Vz_part_e9_2 = reform(arr[6,*])
dE_part_e9_2 = reform(arr[7,*])
Etot_part_e9_2 = reform(arr[8,*])
J_part_e9_2 = reform(arr[9,*])
delE = Etot_part_e9_2-Etot_part_e9_2[0]
q_part_e9_2 = delE/epsilon
dJ_part_e9_2 = J_part_e9_2-J_part_e9_2[0]
scale_dJ_part_e9_2 = dJ_part_e9_2 / (r_tide/r_peri * J_part_e9_2)


;============
;=== plot ===
;============
set_plot,'ps'
@plot_setting
;!p.charsize=1
!p.multi=[0,4,3]
file_out = dir_out+'one_encounter_xy.ps'
device,file=file_out,/color,/landscape;,xoffset=0.5,yoffset=0.5


;=== plot v test ===
; at the moment of closest approach, 10Gyr
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_v200[k]-x_part_v200[0]],[y_part_v200[k]-y_part_v200[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_v200[0,*]-x_part_v200[0],Rs_xy_circle_v200[1,*]-y_part_v200[0]
oplot,Rtidal_xy_circle_v200[0,*]-x_part_v200[0],Rtidal_xy_circle_v200[1,*]-y_part_v200[0],linestyle=2

plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_v200[k]-x_part_v200[0]],[z_part_v200[k]-z_part_v200[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_v200[0,*]-x_part_v200[0],Rs_xz_circle_v200[1,*]-z_part_v200[0]
oplot,Rtidal_xz_circle_v200[0,*]-x_part_v200[0],Rtidal_xz_circle_v200[1,*]-z_part_v200[0],linestyle=2
arrow,-8.3,-3,-8.3,-8,/data

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin-0.6,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_v200_2[k]],[scale_dJ_part_v200_2[k]],psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_v400_2[k]],[scale_dJ_part_v400_2[k]]-0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_v800_2[k]],[scale_dJ_part_v800_2[k]]-0.6,psym=3,color=qcolor[k]
xyouts,0,0,textoidl('v_z=200km/s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.3,textoidl('v_z=400km/s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.6,textoidl('v_z=800km/s'),alignment=0.5,charsize=0.6

; q histogram
plothist,q_part_v200_2,xr=[qmin,qmax],bin=0.1,xtitle='q',yr=[0,510]
plothist,q_part_v400_2,xr=[qmin,qmax],bin=0.1,color=70,/overplot;,/peak;yr=[0,210]
plothist,q_part_v800_2,xr=[qmin,qmax],bin=0.1,color=255,/overplot;,/peak;yr=[0,210]
legend,[textoidl('V_z=200km/s'),textoidl('V_z=400km/s'),textoidl('V_z=800km/s')],box=0,color=[0,70,255],psym=8,/right,charsize=0.6


;=== plot b test ===
; xy plane : at the moment of closest approach, 10Gyr
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_b0Rs[k]-x_part_b0Rs[0]],[y_part_b0Rs[k]-y_part_b0Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_b0Rs[0,*]-x_part_b0Rs[0],Rs_xy_circle_b0Rs[1,*]-y_part_b0Rs[0]
oplot,Rtidal_xy_circle_b0Rs[0,*]-x_part_b0Rs[0],Rtidal_xy_circle_b0Rs[1,*]-y_part_b0Rs[0],linestyle=2
oplot,Rs_xy_circle_b1Rs[0,*]-x_part_b1Rs[0],Rs_xy_circle_b1Rs[1,*]-y_part_b1Rs[0],color=70
oplot,Rtidal_xy_circle_b1Rs[0,*]-x_part_b1Rs[0],Rtidal_xy_circle_b1Rs[1,*]-y_part_b1Rs[0],linestyle=2,color=70
oplot,Rs_xy_circle_b2Rs[0,*]-x_part_b2Rs[0],Rs_xy_circle_b2Rs[1,*]-y_part_b2Rs[0],color=150
oplot,Rtidal_xy_circle_b2Rs[0,*]-x_part_b2Rs[0],Rtidal_xy_circle_b2Rs[1,*]-y_part_b2Rs[0],linestyle=2,color=150
oplot,Rs_xy_circle_b4Rs[0,*]-x_part_b4Rs[0],Rs_xy_circle_b4Rs[1,*]-y_part_b4Rs[0],color=255
oplot,Rtidal_xy_circle_b4Rs[0,*]-x_part_b4Rs[0],Rtidal_xy_circle_b4Rs[1,*]-y_part_b4Rs[0],linestyle=2,color=255
legend,[textoidl('b=0R_s'),textoidl('b=1R_s'),textoidl('b=2R_s'),textoidl('b=4R_s')],box=0,color=[0,70,150,255],psym=8,/right,charsize=0.6

; xz plane
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_b0Rs[k]-x_part_b0Rs[0]],[z_part_b0Rs[k]-z_part_b0Rs[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_b0Rs[0,*]-x_part_b0Rs[0],Rs_xz_circle_b0Rs[1,*]-z_part_b0Rs[0]
oplot,Rtidal_xz_circle_b0Rs[0,*]-x_part_b0Rs[0],Rtidal_xz_circle_b0Rs[1,*]-z_part_b0Rs[0],linestyle=2
arrow,-8.3,-3,-8.3,-8,/data

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin-0.9,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_b0Rs_2[k]],[scale_dJ_part_b0Rs_2[k]],psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_b1Rs_2[k]],[scale_dJ_part_b1Rs_2[k]]-0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_b2Rs_2[k]],[scale_dJ_part_b2Rs_2[k]]-0.6,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_b4Rs_2[k]],[scale_dJ_part_b4Rs_2[k]]-0.9,psym=3,color=qcolor[k]
xyouts,0,0,textoidl('b=0R_s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.3,textoidl('b=1R_s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.6,textoidl('b=2R_s'),alignment=0.5,charsize=0.6
xyouts,0,0-0.9,textoidl('b=4R_s'),alignment=0.5,charsize=0.6

; q histogram
plothist,q_part_b0Rs_2,xr=[qmin,qmax],bin=0.1,xtitle='q',yr=[0,240]
plothist,q_part_b1Rs_2,xr=[qmin,qmax],bin=0.1,color=70,/overplot;,/peak;yr=[0,210]
plothist,q_part_b2Rs_2,xr=[qmin,qmax],bin=0.1,color=150,/overplot;,/peak;yr=[0,210]
plothist,q_part_b4Rs_2,xr=[qmin,qmax],bin=0.1,color=255,/overplot;,/peak;yr=[0,210]
legend,[textoidl('b=0R_s'),textoidl('b=1R_s'),textoidl('b=2R_s'),textoidl('b=4R_s')],box=0,color=[0,70,150,255],psym=8,/right,charsize=0.6


;=== plot Msubhalo test ===
; xy plane : at the moment of closest approach, 10Gyr
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e7[k]-x_part_e7[0]],[y_part_e7[k]-y_part_e7[0]],psym=3,color=qcolor[k]
oplot,Rs_xy_circle_e7[0,*]-x_part_e7[0],Rs_xy_circle_e7[1,*]-y_part_e7[0],color=70
oplot,Rtidal_xy_circle_e7[0,*]-x_part_e7[0],Rtidal_xy_circle_e7[1,*]-y_part_e7[0],linestyle=2,color=70
oplot,Rs_xy_circle_e6[0,*]-x_part_e6[0],Rs_xy_circle_e6[1,*]-y_part_e6[0],color=0
oplot,Rtidal_xy_circle_e6[0,*]-x_part_e6[0],Rtidal_xy_circle_e6[1,*]-y_part_e6[0],linestyle=2,color=0
oplot,Rs_xy_circle_e8[0,*]-x_part_e8[0],Rs_xy_circle_e8[1,*]-y_part_e8[0],color=150
oplot,Rtidal_xy_circle_e8[0,*]-x_part_e8[0],Rtidal_xy_circle_e8[1,*]-y_part_e8[0],linestyle=2,color=150
oplot,Rs_xy_circle_e9[0,*]-x_part_e9[0],Rs_xy_circle_e9[1,*]-y_part_e9[0],color=255
oplot,Rtidal_xy_circle_e9[0,*]-x_part_e9[0],Rtidal_xy_circle_e9[1,*]-y_part_e9[0],linestyle=2,color=255
;legend,[textoidl('e6')textoidl('e7'),,textoidl('e8')],box=0,color=[0,70,255],psym=8,/right,charsize=0.6
legend,[textoidl('e6'),textoidl('e7'),textoidl('e8'),textoidl('e9')],box=0,color=[0,70,150,255],psym=8,/right,charsize=0.6

; xz plane
plot,[0],[0],/nodata,xr=[-14,14],yr=[-14,14],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
for k=0,N_particle-1 do oplot,[x_part_e7[k]-x_part_e7[0]],[z_part_e7[k]-z_part_e7[0]],psym=3,color=qcolor[k]
oplot,Rs_xz_circle_e7[0,*]-x_part_e7[0],Rs_xz_circle_e7[1,*]-z_part_e7[0],color=70
oplot,Rtidal_xz_circle_e7[0,*]-x_part_e7[0],Rtidal_xz_circle_e7[1,*]-z_part_e7[0],linestyle=2,color=70
arrow,-8.3,-3,-8.3,-8,/data
oplot,Rs_xz_circle_e6[0,*]-x_part_e6[0],Rs_xz_circle_e6[1,*]-z_part_e6[0],color=0
oplot,Rtidal_xz_circle_e6[0,*]-x_part_e6[0],Rtidal_xz_circle_e6[1,*]-z_part_e6[0],linestyle=2,color=0
oplot,Rs_xz_circle_e8[0,*]-x_part_e8[0],Rs_xz_circle_e8[1,*]-z_part_e8[0],color=150
oplot,Rtidal_xz_circle_e8[0,*]-x_part_e8[0],Rtidal_xz_circle_e8[1,*]-z_part_e8[0],linestyle=2,color=150
oplot,Rs_xz_circle_e9[0,*]-x_part_e9[0],Rs_xz_circle_e9[1,*]-z_part_e9[0],color=255
oplot,Rtidal_xz_circle_e9[0,*]-x_part_e9[0],Rtidal_xz_circle_e9[1,*]-z_part_e9[0],linestyle=2,color=255

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin-2.0,qmax],yr=[dJmin-0.8,dJmax+0.3],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,N_particle-1 do oplot,[q_part_e7_2[k]],[scale_dJ_part_e7_2[k]],psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e6_2[k]],[scale_dJ_part_e6_2[k]]+0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e8_2[k]],[scale_dJ_part_e8_2[k]]-0.3,psym=3,color=qcolor[k]
for k=0,N_particle-1 do oplot,[q_part_e9_2[k]],[scale_dJ_part_e9_2[k]]-0.6,psym=3,color=qcolor[k]
xyouts,0,0,textoidl('M_{sub}=e7'),alignment=0.5,charsize=0.6
xyouts,0,0+0.3,textoidl('M_{sub}=e6'),alignment=0.5,charsize=0.6
xyouts,0,0-0.3,textoidl('M_{sub}=e8'),alignment=0.5,charsize=0.6
xyouts,0,0-0.6,textoidl('M_{sub}=e9'),alignment=0.5,charsize=0.6

; q histogram
plothist,q_part_e7_2,xr=[qmin-2.0,qmax],bin=0.1,xtitle='q',yr=[0,780]
plothist,q_part_e6_2,bin=0.1,color=70,/overplot
plothist,q_part_e8_2,bin=0.1,color=150,/overplot
plothist,q_part_e9_2,bin=0.1,color=255,/overplot
legend,[textoidl('e6'),textoidl('e7'),textoidl('e8'),textoidl('e9')],box=0,color=[0,70,150,255],psym=8,/right,charsize=0.6

device,/close


END
