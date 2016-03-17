pro one_encounter_snap_movie

dir_out = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_nosub = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_nosub/'
dir_part_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_v200/'
dir_subhalo_v200 = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_v200/'
dir_part_VxVy = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_VxVy/'
dir_subhalo_VxVy = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy_VxVy/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
r_tide = 0.108653
Msat = 20000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,dir_out+'part_peri_nosub',0,arr,/silent
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
chr_rdtbl, dir_part_v200+'snap00000',0,arr,/silent
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
qmin = qmin-0.2
qmax = qmax+0.2
dJmin = dJmin-0.05
dJmax = dJmax+0.05

;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part_VxVy+'snap00000',0,arr,/silent
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial_vxvy = delE/epsilon
dJ_initial_vxvy = J_initial-J_initial[0]
scale_dJ_initial_vxvy = dJ_initial / (r_tide/r_peri * J_initial)
qmin_vxvy = min(q_initial)
qmax_vxvy = max(q_initial)
dJmin_vxvy = min(scale_dJ_initial)
dJmax_vxvy = max(scale_dJ_initial)
qcolor_vxvy = round((q_initial_vxvy-qmin_vxvy)/(qmax_vxvy-qmin_vxvy)*245.)+10
N_particle = N_elements(qcolor_vxvy)
qmin_vxvy = qmin-0.2
qmax_vxvy = qmax+0.2
dJmin_vxvy = dJmin-0.05
dJmax_vxvy = dJmax+0.05


;=== read out simulation info ===
chr_rdtbl,'orbit_subhalo_e7.dat',0,arr,/silent
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)
N_time = N_elements(fname_part_v200)


spawn,'ls '+dir_part_v200+' > temp'
chr_rdtbl,'temp',0,fname_part,/silent
spawn,'ls '+dir_subhalo_v200+' > temp'
chr_rdtbl,'temp',0,fname_subhalo,/silent
spawn,'ls '+dir_part_VxVy+' > temp'
chr_rdtbl,'temp',0,fname_part_vxvy,/silent
spawn,'ls '+dir_subhalo_VxVy+' > temp'
chr_rdtbl,'temp',0,fname_subhalo_vxvy,/silent
spawn,'rm -f temp'


set_plot,'ps'
@plot_setting
;p.charsize=2.0
!p.multi=[0,3,2]
file_out = dir_out+'one_encounter_xy_snap_movie.ps'
device,file=file_out,/color,xoffset=0.2,yoffset=0.2

;for i=0,N_elements(fname_part)-1 do begin
for i=950,1050 do begin
  ;=== v200 ===
  chr_rdtbl,dir_part_v200+fname_part[i],0,arr,/silent
  arr = double(arr)
  t_part_v200 = strtrim(string(reform(arr[0,0])/1000.-10.,f='(f8.3)'),2)+'Gyr'
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
  chr_rdtbl,dir_subhalo_v200+fname_subhalo[i],0,arr,/silent
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
  
  ;============
  ;=== plot ===
  ;============
  
  ;=== 20Myr earlier than the closest approach ===
  plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltaz'),/isotropic
  for k=0,N_particle-1 do oplot,[x_part_v200[k]-x_part_v200[0]],[z_part_v200[k]-z_part_v200[0]],psym=3,color=qcolor[k]
  oplot,Rs_xz_circle_v200[0,*]-x_part_v200[0],Rs_xz_circle_v200[1,*]-z_part_v200[0]
  oplot,Rtidal_xz_circle_v200[0,*]-x_part_v200[0],Rtidal_xz_circle_v200[1,*]-z_part_v200[0],linestyle=2
  legend,['dt='+t_part_v200],box=0,charsize=0.8,/bottom,/right
  
  ; energy & angular momentum 
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part_v200[k]],[scale_dJ_part_v200[k]],psym=3,color=qcolor[k]
  
  ; q histogram
  plothist,q_part_v200,xr=[qmin,qmax],bin=0.1,yr=[0,510],xtitle='q'


;=====================================
;=== encounter pass along the tail ===
;=====================================
  ;=== VxVy ===
  ; read out data of the tail particles 
  chr_rdtbl,dir_part_VxVy+fname_part_VxVy[i],0,arr,/silent
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
  scale_dJ_part_VxVy = dJ_part_VxVy / (r_tide/r_peri * J_part_VxVy)
  
  ; read out data of the subhalo 
  chr_rdtbl,dir_subhalo_VxVy+fname_subhalo_VxVy[i],0,arr,/silent
  arr = double(arr)
  x_subhalo_VxVy = reform(arr[1,*])
  y_subhalo_VxVy = reform(arr[2,*])
  z_subhalo_VxVy = reform(arr[3,*])
  Vx_subhalo_VxVy = reform(arr[4,*])
  Vy_subhalo_VxVy = reform(arr[5,*])
  Vz_subhalo_VxVy = reform(arr[6,*])
  
  ; make a circle of a subhalo 
  Rs_xy_circle_VxVy = circle(x_subhalo_VxVy[0],y_subhalo_VxVy[0],double(Rs[0]))
  Rtidal_xy_circle_VxVy = circle(x_subhalo_VxVy[0],y_subhalo_VxVy[0],double(Rtidal[0]))
  Rs_yz_circle_VxVy = circle(y_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rs[0]))
  Rtidal_yz_circle_VxVy = circle(y_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rtidal[0]))
  Rs_xz_circle_VxVy = circle(x_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rs[0]))
  Rtidal_xz_circle_VxVy = circle(x_subhalo_VxVy[0],z_subhalo_VxVy[0],double(Rtidal[0]))
  
  ;============
  ;=== plot ===
  ;============
  ;=== 40Myr earlier than the closest approach ===
  plot,[0],[0],/nodata,xr=[-15,15],yr=[-15,15],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_part_VxVy[k]-x_part_VxVy[0]],[y_part_VxVy[k]-y_part_VxVy[0]],psym=3,color=qcolor_vxvy[k]
  oplot,Rs_xy_circle_VxVy[0,*]-x_part_VxVy[0],Rs_xy_circle_VxVy[1,*]-y_part_VxVy[0]
  oplot,Rtidal_xy_circle_VxVy[0,*]-x_part_VxVy[0],Rtidal_xy_circle_VxVy[1,*]-y_part_VxVy[0],linestyle=2
  
  ; energy & angular momentum 
  plot,[0],[0],/nodata,xr=[qmin_vxvy,qmax_vxvy],yr=[dJmin_vxvy,dJmax_vxvy],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part_VxVy[k]],[scale_dJ_part_VxVy[k]],psym=3,color=qcolor_vxvy[k]
  
  ; q histogram
  plothist,q_part_VxVy,xr=[qmin_vxvy,qmax_vxvy],bin=0.1,yr=[0,410],xtitle='q'
endfor
device,/close

spawn,'convert '+dir_out+'one_encounter_xy_snap_movie.ps '+dir_out+'one_encounter_xy_snap_movie.gif'

END
